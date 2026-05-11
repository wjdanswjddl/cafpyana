#!/usr/bin/env python3
"""
Smoke-test numucc_1p0pi dataframe configs that build MC multi/unisim weights.

  # Show what will run (no imports / CAFs)
  python analysis_village/numucc_1p0pi/scripts/test_wgt_df_configs.py --list-cases

  # Fast: load every config and check DFS / ARGS / NAMES (needs same Python deps as batch jobs)
  python analysis_village/numucc_1p0pi/scripts/test_wgt_df_configs.py --dry-run

  # Integration: one flatcaf file, subset of configs (-nfile 1, pool mode)
  python analysis_village/numucc_1p0pi/scripts/test_wgt_df_configs.py \\
      --caf /path/to/sample.flat.caf.root

  # Include expensive \"all knob groups\" GENIE / flux configs when testing locally
  python analysis_village/numucc_1p0pi/scripts/test_wgt_df_configs.py --caf sample.root --full-multisim

  # Site-specific preprocess sample (optional)
  python analysis_village/numucc_1p0pi/scripts/test_wgt_df_configs.py --dry-run --with-add-ar23p

Run from anywhere; discovers cafpyana root from this script path (must contain run_df_maker.py).
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
import textwrap
from pathlib import Path

_ENV_KEYS_ROTATING = ("GENIE_KNOB_GROUP", "FLUX_GROUP")


def _cafpyana_root(cli_root: Path | None) -> Path:
    if cli_root is not None:
        r = cli_root.resolve()
    else:
        r = Path(__file__).resolve().parents[3]
    if not (r / "run_df_maker.py").is_file():
        raise SystemExit("cafpyana root does not contain run_df_maker.py: %s" % r)
    return r


def _cases():
    """(config_relative_path, env_extra dict, label, df_quick_ok).

    df_quick_ok: included in --caf runs unless --full-multisim is set (flux-all/genie-all are heavy).
    """
    flux_groups = ["beam", "hadron", "xsec"]
    genie_one_offs = ["CCQE", "MEC"]

    out = []

    out.append(("configs/numucc_1p0pi/sel_mup-g4wgts.py", {}, "mup-g4wgts", True))
    out.append(("configs/numucc_1p0pi/sel_mup-mcstatwgts.py", {}, "mup-mcstatwgts", True))
    out.append(("configs/numucc_1p0pi/sel_2prong-wgts-mc.py", {}, "2prong-wgts-mc", True))

    out.append(("configs/numucc_1p0pi/sel_mup-fluxwgts-knobgroups.py", {}, "flux-all-groups", False))
    for g in flux_groups:
        out.append(
            (
                "configs/numucc_1p0pi/sel_mup-fluxwgts-knobgroups.py",
                {"FLUX_GROUP": g},
                "flux-group-%s" % g.lower(),
                True,
            )
        )

    out.append(("configs/numucc_1p0pi/sel_mup-geniewgts-knobgroups.py", {}, "genie-all-groups", False))
    for g in genie_one_offs:
        out.append(
            (
                "configs/numucc_1p0pi/sel_mup-geniewgts-knobgroups.py",
                {"GENIE_KNOB_GROUP": g},
                "genie-group-%s" % g.lower(),
                True,
            )
        )

    # Optional preprocess + GENIE Ar23p weights (paths inside config must exist on your machine)
    out.append(("configs/numucc_1p0pi/add_ar23p.py", {}, "add-ar23p-preprocess", False))

    return out


def _env_for_subprocess(extra: dict) -> dict:
    e = os.environ.copy()
    for k in _ENV_KEYS_ROTATING:
        e.pop(k, None)
    e.update(extra)
    return e


def _validate_once(root: Path, cfg_rel: str, env_extra: dict, label: str) -> None:
    cfg_abs = root / cfg_rel
    if not cfg_abs.is_file():
        raise FileNotFoundError(cfg_abs)

    env_json = json.dumps(env_extra)
    root_s = str(root)
    cfg_s = str(cfg_abs)

    snippet = textwrap.dedent(
        """
        import json, os, sys
        from pathlib import Path

        root = Path(%r)
        os.chdir(root)
        if str(root) not in sys.path:
            sys.path.insert(0, str(root))

        for k in %r:
            os.environ.pop(k, None)
        os.environ.update(json.loads(%r))

        g = {"__name__": "__main__"}
        exec(open(%r, encoding="utf-8").read(), g)
        assert "DFS" in g and "ARGS" in g and "NAMES" in g
        assert len(g["DFS"]) == len(g["ARGS"]) == len(g["NAMES"])
        print("validated:", %r, "n_tables=", len(g["NAMES"]))
        """
        % (root_s, list(_ENV_KEYS_ROTATING), env_json, cfg_s, label)
    )

    subprocess.run(
        [sys.executable, "-c", snippet],
        cwd=root_s,
        check=True,
        timeout=600,
    )


def _run_df_maker(
    root: Path,
    cfg_rel: str,
    env_extra: dict,
    label: str,
    caf_path: Path,
    work_dir: Path,
) -> None:
    out_prefix = work_dir / ("smoke_%s" % label.replace("/", "_"))
    cmd = [
        sys.executable,
        str(root / "run_df_maker.py"),
        "-c",
        str(cfg_rel),
        "-o",
        str(out_prefix),
        "-i",
        str(caf_path),
        "-nfile",
        "1",
        "-ncpu",
        "1",
    ]
    subprocess.run(
        cmd,
        cwd=str(root),
        env=_env_for_subprocess(env_extra),
        check=True,
        timeout=7200,
    )


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cafpyana-root", type=Path, default=None, help="cafpyana checkout (default: infer from script location)")
    ap.add_argument("--list-cases", action="store_true", help="Print the config/env matrix and exit")
    ap.add_argument("--dry-run", action="store_true", help="Only exec configs and assert DFS/ARGS/NAMES (recommended before batch)")
    ap.add_argument("--caf", type=Path, default=None, help="One flatcaf ROOT file for pool-mode smoke tests (-nfile 1)")
    ap.add_argument(
        "--full-multisim",
        action="store_true",
        help="With --caf, also run flux-all-groups and genie-all-groups (slow).",
    )
    ap.add_argument(
        "--with-add-ar23p",
        action="store_true",
        help="Include configs/numucc_1p0pi/add_ar23p.py (site-specific preprocess paths; dry-run and/or --caf).",
    )
    args = ap.parse_args()

    root = _cafpyana_root(args.cafpyana_root)
    cases = _cases()

    if args.list_cases:
        print("config_rel | env | label | df_quick_default")
        for cfg_rel, env_extra, label, df_quick in cases:
            print(" ", cfg_rel, "|", env_extra or "{}", "|", label, "|", df_quick)
        return 0

    failed = []

    if args.dry_run:
        n_dry = sum(
            1
            for _cr, _env, lab, _dq in cases
            if lab != "add-ar23p-preprocess" or args.with_add_ar23p
        )
        print("Dry-run validation (%d cases), cwd=%s" % (n_dry, root))
        for cfg_rel, env_extra, label, _dfq in cases:
            if label == "add-ar23p-preprocess" and not args.with_add_ar23p:
                continue
            try:
                print("  [%s] %s env=%s" % (label, cfg_rel, env_extra or "{}"))
                _validate_once(root, cfg_rel, env_extra, label)
            except Exception as ex:
                print("  FAIL:", ex)
                failed.append(("dry-run", label, str(ex)))

    if args.caf is not None:
        caf_path = args.caf.expanduser().resolve()
        if not caf_path.is_file():
            print("CAF file not found:", caf_path, file=sys.stderr)
            return 1

        df_cases = []
        for cfg_rel, env_extra, label, df_quick in cases:
            if label == "add-ar23p-preprocess":
                if args.with_add_ar23p:
                    df_cases.append((cfg_rel, env_extra, label))
                continue
            if df_quick or args.full_multisim:
                df_cases.append((cfg_rel, env_extra, label))

        print("Pool-mode smoke tests (%d cases), caf=%s" % (len(df_cases), caf_path))
        tmp = Path(tempfile.mkdtemp(prefix="cafpyana_wgt_smoke_", dir=str(Path.cwd())))
        try:
            for cfg_rel, env_extra, label in df_cases:
                try:
                    print("  [%s] %s env=%s" % (label, cfg_rel, env_extra or "{}"))
                    _run_df_maker(root, cfg_rel, env_extra, label, caf_path, tmp)
                    df_out = tmp / ("smoke_%s.df" % label.replace("/", "_"))
                    if not df_out.is_file():
                        raise RuntimeError("expected output missing: %s" % df_out)
                except Exception as ex:
                    print("  FAIL:", ex)
                    failed.append(("df-run", label, str(ex)))
        finally:
            shutil.rmtree(tmp, ignore_errors=True)

    if not args.dry_run and args.caf is None:
        ap.print_help()
        print("\nProvide --dry-run and/or --caf <file.flat.caf.root>", file=sys.stderr)
        return 1

    if failed:
        print("\nFailures:")
        for mode, label, msg in failed:
            print(" ", mode, label, "->", msg)
        return 1

    print("\nAll requested checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
