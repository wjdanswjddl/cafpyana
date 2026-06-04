#!/usr/bin/env python
"""Process a SINGLE input file through the numuCC 1p0pi event selection.

HDF5 ``.df`` files store splits ``evt_0``, ``evt_1``, ... Default (**splits** mode)
loads one split at a time and accumulates histograms. **concat** mode uses ``load_dfs``.

Use ``--trace`` when jobs vanish without a Python traceback: it tees a log file under
``--out_dir`` and prints every pipeline stage so you can see the last step executed.

Usage
-----
    python event_selection_chunk.py --df_file PATH.df \\
                                     --sample {mc,data,intime,offbeam,dirt} \\
                                     --out_dir OUT_DIR \\
                                     --trace

GENIE MC uses unit ``pot_weight`` per event by default; aggregation scales MC to data POT.
For GiBUU (or other generators that need ``mc.genweight`` per event), add ``--use-mc-genweight``.
"""
from __future__ import annotations

import argparse
import os

# Non-interactive backend if anything downstream touches matplotlib.
os.environ.setdefault("MPLBACKEND", "Agg")

import faulthandler
import gc
import platform
import resource
import signal
import sys
import traceback
from os import path

import numpy as np
import pandas as pd

import warnings

warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

sys.path.append(path.dirname(path.dirname(path.dirname(path.dirname(path.abspath(__file__))))))


TRACE_EPILOG = """
Diagnostics (when the process dies with no Python traceback):
-------------------------------------------------------------
Exit code (run: echo $? right after the job):
  137 (128+SIGKILL)  Often Linux OOM killer or cgroup memory limit — check:
                       dmesg -T | tail -50
                       Or Slurm: sacct -j JOBID -o State,MaxRSS,Elapsed,ExitCode
  139 (128+SIGSEGV)  Native fault (HDF5/PyTables/BLAS). Try the same file with
                       gdb --args python …/event_selection_chunk.py …
  143 (128+SIGTERM)  Scheduler preemption, walltime, or manual kill.

Tracing:
  --trace writes chunk_trace__<sample>__<stem>.log under --out_dir (tee stdout/stderr).
  The LAST printed "[pipeline] …" line is the last completed step; death happened
  before the following line appeared.

If the job hangs, find PID and run:  kill -USR1 <pid>
  (with --trace, registers a handler that dumps thread stacks to stderr/log.)

Keep the trace log + exact command line when reporting the issue.
"""


class _TeeTextIO:
    """Duplicate writes to several text streams (unbuffered flush per write)."""

    __slots__ = ("streams",)

    def __init__(self, *streams):
        self.streams = streams

    def write(self, data: str):
        for s in self.streams:
            s.write(data)
            try:
                s.flush()
            except Exception:
                pass
        return len(data)

    def flush(self):
        for s in self.streams:
            try:
                s.flush()
            except Exception:
                pass

    def isatty(self):
        return False


def _proc_status_mb() -> tuple[float, float, float]:
    rss_kb = hwm_kb = peak_kb = 0
    try:
        with open("/proc/self/status") as f:
            for line in f:
                if line.startswith("VmRSS:"):
                    rss_kb = int(line.split()[1])
                elif line.startswith("VmHWM:"):
                    hwm_kb = int(line.split()[1])
                elif line.startswith("VmPeak:"):
                    peak_kb = int(line.split()[1])
    except OSError:
        pass
    fac = 1.0 / 1024.0
    return rss_kb * fac, hwm_kb * fac, peak_kb * fac


def _mem_line(tag: str, scale_mb: float) -> str:
    rss, hwm, vpeak = _proc_status_mb()
    width = 28
    if scale_mb > 0:
        frac = max(0.0, min(1.0, rss / scale_mb))
        filled = int(round(frac * width))
        bar = "#" * filled + "-" * (width - filled)
    else:
        bar = "?" * width
    return (
        f"[chunk][mem] {tag:36s} |{bar}| "
        f"RSS={rss:8.1f} MB  VmHWM={hwm:8.1f}  VmPeak={vpeak:8.1f}"
    )


def _print_diag_banner(args):
    print(
        f"[chunk][diag] pid={os.getpid()} ppid={os.getppid()} host={platform.node()}",
        flush=True,
    )
    print(f"[chunk][diag] python={sys.version.split()[0]} exe={sys.executable}", flush=True)
    print(f"[chunk][diag] cwd={os.getcwd()}", flush=True)
    for attr in ("RLIMIT_AS", "RLIMIT_STACK", "RLIMIT_CPU", "RLIMIT_DATA"):
        if hasattr(resource, attr):
            soft, hard = resource.getrlimit(getattr(resource, attr))
            print(f"[chunk][diag] rlimit {attr} soft={soft} hard={hard}", flush=True)


def _install_signal_loggers():
    def _log(sig, frame):
        try:
            name = signal.Signals(sig).name
        except Exception:
            name = str(sig)
        sys.__stderr__.write(f"[chunk][signal] caught {name} ({sig})\n")
        sys.__stderr__.flush()

    for sig in (signal.SIGTERM, signal.SIGINT, signal.SIGQUIT):
        try:
            signal.signal(sig, _log)
        except Exception:
            pass


from analysis_village.numucc_1p0pi.event_selection_pipeline_def import build_runner
from analysis_village.numucc_1p0pi.evt_derived_kinematics import (
    ensure_derived_trk_kinematics_cols,
    ensure_mc_level_phi_mcnu,
)
from analysis_village.numucc_1p0pi.selection_framework import multicol_resolve_column_key
from pyanalib.pandas_helpers import pad_column_name
from pyanalib.split_df_helpers import get_n_split, load_dfs


def _prefix_mcnu_columns(mc_nu_df: pd.DataFrame) -> None:
    """Ensure ``mcnu`` columns have a leading ``mc`` level (same as GENIE / legacy chunk-map)."""
    if isinstance(mc_nu_df.columns, pd.MultiIndex):
        try:
            first_level = mc_nu_df.columns.get_level_values(0)
            need_prefix = not np.all(first_level == "mc")
        except Exception:
            need_prefix = True
        if need_prefix:
            mc_nu_df.columns = pd.MultiIndex.from_tuples(
                [tuple(["mc"] + list(c)) for c in mc_nu_df.columns]
            )


def _ensure_trk_phi_col(trk_df: pd.DataFrame | None) -> None:
    """Add ``pfp.trk.phi`` (degrees) from ``pfp.trk.dir.{x,y}`` when missing (same as pipeline)."""
    if trk_df is None or len(trk_df) == 0:
        return
    if not isinstance(trk_df.columns, pd.MultiIndex):
        return
    if multicol_resolve_column_key(trk_df, ("pfp", "trk", "phi", "", "", "")) is not None:
        return
    kx = multicol_resolve_column_key(trk_df, ("pfp", "trk", "dir", "x", ""))
    ky = multicol_resolve_column_key(trk_df, ("pfp", "trk", "dir", "y", ""))
    if kx is None or ky is None:
        return
    phi_col = pad_column_name(("pfp", "trk", "phi", "", "", ""), trk_df)
    trk_df.loc[:, phi_col] = np.degrees(
        np.arctan2(
            np.asarray(trk_df.loc[:, kx], dtype=float),
            np.asarray(trk_df.loc[:, ky], dtype=float),
        )
    )


def _ensure_phi_and_kinematics_cols(
    evt_df: pd.DataFrame,
    trk_df: pd.DataFrame | None,
    mcnu_df: pd.DataFrame | None,
) -> tuple[pd.DataFrame, pd.DataFrame | None]:
    """Add reco/truth/track ``phi`` and related columns so ``VariableConfig`` keys resolve in plots."""
    evt_df = ensure_derived_trk_kinematics_cols(evt_df)
    _ensure_trk_phi_col(trk_df)
    if mcnu_df is not None and len(mcnu_df) > 0:
        _prefix_mcnu_columns(mcnu_df)
        mcnu_df = ensure_mc_level_phi_mcnu(mcnu_df)
    return evt_df, mcnu_df


def _hdf_has_mcnu(df_file: str) -> bool:
    """Return True if this split HDF file stores ``mcnu_<i>`` datasets."""
    try:
        with pd.HDFStore(df_file, mode="r") as store:
            keys = store.keys()
        return any(str(k).startswith("/mcnu_") for k in keys)
    except Exception:
        return False


def _hdr_chunk_pot(hdr_df: pd.DataFrame | None) -> float:
    if hdr_df is None or "pot" not in hdr_df.columns:
        return 0.0
    return float(hdr_df["pot"].sum())


def _hdr_data_gates_bnb(hdr_df: pd.DataFrame | None) -> float:
    if hdr_df is None or "nbnbinfo" not in hdr_df.columns:
        return 0.0
    return float(hdr_df["nbnbinfo"].sum())


def _hdr_cosmic_gates_intime(hdr_df: pd.DataFrame | None) -> float:
    if hdr_df is None or "ngenevt" not in hdr_df.columns:
        return 0.0
    return float(hdr_df.loc[hdr_df["first_in_subrun"] == 1, "ngenevt"].sum())


def _hdr_cosmic_gates_offbeam(hdr_df: pd.DataFrame | None) -> float:
    if hdr_df is None or "noffbeambnb" not in hdr_df.columns:
        return 0.0
    return float(hdr_df.loc[hdr_df["first_in_subrun"] == 1, "noffbeambnb"].sum())


def _intrinsic_weight_series(
    df: pd.DataFrame | None,
    sample: str,
    use_mc_genweight: bool = False,
) -> np.ndarray:
    """Per-event weight before global POT / gates scaling in the aggregator.

    Default (GENIE-style MC): MC and dirt use **1.0** per event; compare to data via
    ``apply_global_exposure_scales``. Optional ``use_mc_genweight`` multiplies MC/dirt
    by ``mc.genweight`` (e.g. GiBUU).
    """
    if df is None or len(df) == 0:
        return np.ones(0, dtype=float)
    if sample == "data":
        return np.ones(len(df), dtype=float)
    if sample in ("intime", "offbeam"):
        return np.ones(len(df), dtype=float)
    if sample in ("mc", "dirt"):
        if not use_mc_genweight:
            return np.ones(len(df), dtype=float)
        try:
            gw = df["mc"]["genweight"]
            w = np.asarray(gw, dtype=float).reshape(-1)
            # Any NaN weight poisons every bin of np.histogram(..., weights=w).
            return np.nan_to_num(w, nan=0.0, posinf=0.0, neginf=0.0)
        except Exception:
            return np.ones(len(df), dtype=float)
    raise ValueError(sample)


def attach_intrinsic_weights(
    evt_df: pd.DataFrame | None,
    trk_df: pd.DataFrame | None,
    sample: str,
    use_mc_genweight: bool = False,
):
    if evt_df is not None and len(evt_df) > 0:
        evt_df["pot_weight"] = _intrinsic_weight_series(evt_df, sample, use_mc_genweight)
    if trk_df is not None and len(trk_df) > 0:
        trk_df["pot_weight"] = _intrinsic_weight_series(trk_df, sample, use_mc_genweight)


def _accumulate_hdr_meta(sample: str, hdr_df: pd.DataFrame | None,
                         chunk_pot: list[float], chunk_gates_bnb: list[float],
                         chunk_cosmic_gates_intime: list[float],
                         chunk_cosmic_gates_offbeam: list[float]):
    chunk_pot[0] += _hdr_chunk_pot(hdr_df)
    if sample == "data":
        chunk_gates_bnb[0] += _hdr_data_gates_bnb(hdr_df)
    elif sample == "intime":
        chunk_cosmic_gates_intime[0] += _hdr_cosmic_gates_intime(hdr_df)
    elif sample == "offbeam":
        chunk_cosmic_gates_offbeam[0] += _hdr_cosmic_gates_offbeam(hdr_df)


def parse_args():
    p = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=TRACE_EPILOG,
    )
    p.add_argument("--df_file", required=True, help="Path to a single .df file to process")
    p.add_argument("--sample", required=True,
                   choices=["mc", "data", "intime", "offbeam", "dirt"],
                   help="Which sample type this chunk represents")
    p.add_argument("--out_dir", required=True, help="Output directory for the pickle")
    p.add_argument("--out_tag", default="",
                   help="Optional tag appended to the output pickle name")
    p.add_argument("--max_splits", type=int, default=0,
                   help="Cap HDF5 splits merged in concat mode, or processed in splits mode; "
                        "0 means all splits in file.")
    p.add_argument(
        "--load_mode",
        choices=("splits", "concat"),
        default="splits",
        help="splits: sequential slices (lower peak RAM). concat: single pd.concat via load_dfs.",
    )
    p.add_argument("--no_mem_diag", action="store_true",
                   help="Disable /proc/self/status RSS lines")
    p.add_argument("--mem_bar_scale_mb", type=float, default=0.0,
                   help="RSS value that fills the ASCII bar (0 -> 32768 MB)")
    p.add_argument("--trace", action="store_true",
                   help="Pipeline stage logging + diagnostic banner + trace log file under out_dir")
    p.add_argument("--trace_log", default=None,
                   help="With --trace, explicit log path (default: out_dir/chunk_trace__...log)")
    p.add_argument(
        "--use-mc-genweight",
        action="store_true",
        dest="use_mc_genweight",
        help="MC/dirt: multiply pot_weight by mc.genweight (e.g. GiBUU). "
             "Default is unit weight per event for GENIE MC; scale to data POT in aggregation.",
    )
    p.add_argument(
        "--mc-univ-syst",
        default="",
        help="MC only: comma-separated multi-universe syst folder names under mc.* whose "
             "univ_i weights are accumulated for chunked covariance (e.g. Flux,G4,GENIE). "
             "Empty disables (smaller pickles). Aggregate uses these by default for overlay syst bands.",
    )
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    log_f = None
    old_stdout, old_stderr = sys.stdout, sys.stderr
    trace_log_path = None
    if args.trace:
        trace_log_path = args.trace_log
        if not trace_log_path:
            stem = path.splitext(path.basename(args.df_file))[0]
            trace_log_path = path.join(args.out_dir, f"chunk_trace__{args.sample}__{stem}.log")
        log_f = open(trace_log_path, "w", buffering=1)
        sys.stdout = _TeeTextIO(sys.__stdout__, log_f)
        sys.stderr = _TeeTextIO(sys.__stderr__, log_f)
        print(f"[chunk][diag] trace log: {trace_log_path}", flush=True)

    try:
        _main_body(args)
    except SystemExit as e:
        if args.trace:
            print(f"[chunk][diag] SystemExit code={e.code}", flush=True)
        raise
    except Exception:
        print(
            "[chunk][diag] --- uncaught Python exception (see stderr below) ---",
            file=sys.stderr,
            flush=True,
        )
        traceback.print_exc(file=sys.stderr)
        sys.stderr.flush()
        sys.stdout.flush()
        sys.exit(1)
    finally:
        if args.trace:
            print("[chunk][diag] leaving main() (tee teardown next)", flush=True)
        if log_f is not None:
            sys.stdout, sys.stderr = old_stdout, old_stderr
            log_f.close()
            if trace_log_path:
                print(f"[chunk] trace log saved: {trace_log_path}", file=sys.__stdout__, flush=True)


def _main_body(args) -> int:
    faulthandler.enable()
    if args.trace:
        try:
            faulthandler.register(signal.SIGUSR1)
            print("[chunk][diag] SIGUSR1 registered -> dump traceback on demand", flush=True)
        except Exception as e:
            print(f"[chunk][diag] could not register SIGUSR1 handler: {e}", flush=True)
        _install_signal_loggers()

    mem_diag = not args.no_mem_diag
    bar_scale = args.mem_bar_scale_mb if args.mem_bar_scale_mb > 0 else 32768.0

    pipeline_trace = None
    if args.trace:
        _print_diag_banner(args)

        def pipeline_trace(msg: str) -> None:
            print(msg, flush=True)

    print(f"[chunk] sample={args.sample}  file={args.df_file}  load_mode={args.load_mode}", flush=True)
    if args.use_mc_genweight and args.sample in ("mc", "dirt"):
        print(
            "[chunk] mc/dirt pot_weight includes mc.genweight (--use-mc-genweight)",
            flush=True,
        )
    if mem_diag:
        print(_mem_line("after interpreter + imports", bar_scale), flush=True)

    keys = ["evt", "trk", "hdr"]
    load_mcnu = args.sample == "mc" and _hdf_has_mcnu(args.df_file)
    if args.sample == "mc" and not load_mcnu:
        print(
            "[chunk] Skipping MC efficiency calculation: this HDF file has no mcnu table "
            "(efficiency vs generated neutrinos requires mcnu; evt-only tables are reconstructed events).",
            flush=True,
        )
    keys_load = keys + (["mcnu"] if load_mcnu else [])

    n_keys = int(get_n_split(args.df_file))
    if args.max_splits and args.max_splits > 0:
        n_use = min(args.max_splits, n_keys)
    else:
        n_use = n_keys
    if n_use <= 0:
        raise SystemExit("[chunk] ERROR: no HDF5 splits to process")

    print(f"[chunk] HDF5 n_split={n_keys}  using {n_use} split(s) ({args.load_mode})", flush=True)
    if load_mcnu:
        print("[chunk] mcnu present → MC efficiency accumulators enabled", flush=True)

    mc_univ_tags = tuple(
        x.strip() for x in (args.mc_univ_syst or "").split(",") if x.strip()
    )
    if mc_univ_tags and args.sample != "mc":
        print(
            f"[chunk] WARN: --mc-univ-syst ignored for sample={args.sample!r} (MC only)",
            flush=True,
        )
        mc_univ_tags = ()
    runner = build_runner(
        args.sample,
        mc_univ_syst_tags=mc_univ_tags if args.sample == "mc" else None,
    )
    if mc_univ_tags:
        print(f"[chunk] mc_univ_syst_tags={mc_univ_tags}", flush=True)
    if mem_diag:
        print(_mem_line("after build_runner", bar_scale), flush=True)

    chunk_pot = [0.0]
    chunk_gates_bnb = [0.0]
    chunk_cosmic_gates_intime = [0.0]
    chunk_cosmic_gates_offbeam = [0.0]
    n_evt_total = 0

    if args.load_mode == "concat":
        cap = args.max_splits if args.max_splits > 0 else 999
        print("[chunk] load_dfs (concat) starting …", flush=True)
        dfs = load_dfs(args.df_file, keys2load=keys_load, n_max_concat=min(cap, n_keys))
        evt_df = dfs["evt"]
        trk_df = dfs["trk"]
        hdr_df = dfs["hdr"]
        mcnu_df = dfs["mcnu"] if load_mcnu else None
        print("[chunk] load_dfs (concat) finished", flush=True)
        if mem_diag:
            print(_mem_line("after load_dfs (concat)", bar_scale), flush=True)

        _accumulate_hdr_meta(args.sample, hdr_df, chunk_pot, chunk_gates_bnb,
                             chunk_cosmic_gates_intime, chunk_cosmic_gates_offbeam)
        attach_intrinsic_weights(
            evt_df, trk_df, args.sample, use_mc_genweight=args.use_mc_genweight
        )
        evt_df, mcnu_df = _ensure_phi_and_kinematics_cols(evt_df, trk_df, mcnu_df)
        n_evt_total = int(len(evt_df))
        if mem_diag:
            print(_mem_line("after attach_intrinsic_weights", bar_scale), flush=True)

        print("[chunk] runner.run (concat) …", flush=True)
        runner.run(
            {"evt": evt_df, "trk": trk_df, "hdr": hdr_df, "mcnu": mcnu_df},
            pipeline_trace=pipeline_trace,
        )
        print("[chunk] runner.run (concat) finished", flush=True)
        if mem_diag:
            print(_mem_line("after runner.run (concat)", bar_scale), flush=True)

        meta_splits_used = min(cap, n_keys)
    else:
        for i in range(n_use):
            print(f"[chunk] read_hdf split {i + 1}/{n_use} …", flush=True)
            dfs = {k: pd.read_hdf(args.df_file, key=f"{k}_{i}") for k in keys}
            if load_mcnu:
                dfs["mcnu"] = pd.read_hdf(args.df_file, key=f"mcnu_{i}")
            evt_df = dfs["evt"]
            trk_df = dfs["trk"]
            hdr_df = dfs["hdr"]
            mcnu_df = dfs["mcnu"] if load_mcnu else None
            print(f"[chunk] read_hdf split {i + 1}/{n_use} done", flush=True)

            print(
                f"[chunk] split {i + 1}/{n_use}  evt={len(evt_df)}  trk={len(trk_df)}  hdr={len(hdr_df)}",
                flush=True,
            )
            if mem_diag:
                print(_mem_line(f"split {i + 1}: after read_hdf", bar_scale), flush=True)

            _accumulate_hdr_meta(
                args.sample, hdr_df,
                chunk_pot, chunk_gates_bnb,
                chunk_cosmic_gates_intime, chunk_cosmic_gates_offbeam,
            )
            attach_intrinsic_weights(
                evt_df, trk_df, args.sample, use_mc_genweight=args.use_mc_genweight
            )
            evt_df, mcnu_df = _ensure_phi_and_kinematics_cols(evt_df, trk_df, mcnu_df)
            n_evt_total += int(len(evt_df))
            if mem_diag:
                print(_mem_line(f"split {i + 1}: after weights", bar_scale), flush=True)

            if args.trace:
                print(f"[chunk][trace] ===== split {i + 1}/{n_use} pipeline begin =====", flush=True)
            runner.run(
                {"evt": evt_df, "trk": trk_df, "hdr": hdr_df, "mcnu": mcnu_df},
                pipeline_trace=pipeline_trace,
            )
            if args.trace:
                print(f"[chunk][trace] ===== split {i + 1}/{n_use} pipeline end =====", flush=True)
            if mem_diag:
                print(_mem_line(f"split {i + 1}: after runner.run", bar_scale), flush=True)

            del dfs, evt_df, trk_df, hdr_df
            gc.collect()

        meta_splits_used = n_use

    base = path.splitext(path.basename(args.df_file))[0]
    tag = ("_" + args.out_tag) if args.out_tag else ""
    out_path = path.join(args.out_dir, f"{args.sample}__{base}{tag}.pkl")

    meta = {
        "weight_scheme": "intrinsic",
        "use_mc_genweight": bool(args.use_mc_genweight),
        "mc_univ_syst_tags": list(mc_univ_tags) if mc_univ_tags else [],
        "df_file": args.df_file,
        "sample": args.sample,
        "chunk_pot": chunk_pot[0],
        "chunk_gates_bnb": chunk_gates_bnb[0],
        "chunk_cosmic_gates_intime": chunk_cosmic_gates_intime[0],
        "chunk_cosmic_gates_offbeam": chunk_cosmic_gates_offbeam[0],
        "n_evt": n_evt_total,
        "n_splits_processed": meta_splits_used,
        "n_splits_in_file": n_keys,
        "load_mode": args.load_mode,
        "mc_efficiency_enabled": load_mcnu,
    }
    print(f"[chunk] pickle.dump → {out_path} …", flush=True)
    runner.save(out_path, extra_meta=meta)
    print("[chunk] pickle.dump finished", flush=True)
    if mem_diag:
        print(_mem_line("after pickle save", bar_scale), flush=True)

    print(
        f"[chunk] wrote {out_path}  total_evt={n_evt_total}  meta POT/gates: pot={chunk_pot[0]:.3e} "
        f"bnb_gates={chunk_gates_bnb[0]:.3e} "
        f"cosmic_i={chunk_cosmic_gates_intime[0]:.3e} cosmic_ob={chunk_cosmic_gates_offbeam[0]:.3e}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    main()
