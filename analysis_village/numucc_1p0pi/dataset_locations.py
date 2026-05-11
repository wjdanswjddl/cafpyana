"""Central directory and glob registry for numu CC 1p0pi workflows.

Edit paths **here only** so drivers stay thin:

- ``run_event_selection_chunked.sh`` — map shards: MC/data/intime/offbeam/dirt ``.df`` files for
  ``event_selection_chunk.py`` / ``event_selection_aggregate.py``.
- ``run_syst_multisim_chunked.sh`` — per systematic + CAF shard (each syst may use its own glob).
- ``syst_multisim_aggregate.py`` — merge map outputs into a ``syst_disk_layout`` tree.
- ``syst_detvar_chunk.py`` / ``syst_detvar_aggregate.py`` — WireMod + calo variants;
  input globs are listed in ``DETVAR_WIREMOD_GLOBS`` / :func:`iter_detvar_chunk_jobs`.
- ``get_systematics_genie.py`` / ``run_syst_genie_chunked.sh`` — GENIE reweights: one glob per
  knob **group** (``GENIE_GROUP_GLOBS``) / :func:`iter_genie_chunk_map_tasks`.
- ``run_syst_cosmics_chunked.sh`` — ``syst_cosmics_chunk.py`` / ``syst_cosmics_aggregate.py``;
  globs reuse ``EVENT_SELECTION_GLOBS`` ``offbeam`` / ``intime``.

**Naming:** *HDF splits*, *map shards* (one ``.df`` file), and *exposure batches* (time-ordered
data slices for staged access) are different concepts — see ``exposure_access``.

Relative globs are resolved from ``SPRING_GEN1_ROOT``. Override any constant by
setting environment variables of the same name before importing (advanced).

After producers fill the syst disk tree (``syst_disk_layout``), point ``utils.get_syst_unc`` /
``event_selection_aggregate.py`` at the **root** via ``NUMUCC_SYST_DISK_ROOT`` or
``--syst-disk-root`` (subfolders ``MCstat``, ``Flux``, ``G4``, ``GENIE``, ``Cosmics``,
``Detector``). Loaders require every expected file; there are no alternate search paths.
"""

from __future__ import annotations

import glob
import os
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

import sys
sys.path.append(os.path.join(os.path.dirname(__file__), "..", ".."))
from makedf.geniesyst import *

# -----------------------------------------------------------------------------
# Base release (Spring Gen 1 — edit for other campaigns)
# -----------------------------------------------------------------------------
SPRING_GEN1_ROOT = Path(
    os.environ.get(
        "NUMUCC_SPRING_GEN1_ROOT",
        "/pnfs/sbnd/scratch/users/munjung/cafpyana_out/dfs",
    )
)

PLOTS_BASE = Path(
    os.environ.get(
        "NUMUCC_PLOTS_BASE",
        "/exp/sbnd/data/users/munjung/plots/numucc1p0pi",
    )
)

# -----------------------------------------------------------------------------
# Event selection chunked driver (same logical samples as the old bash arrays)
# Keys: mc, data, intime, offbeam, dirt
# -----------------------------------------------------------------------------
EVENT_SELECTION_GLOBS: Dict[str, str] = {
    "mc": str(SPRING_GEN1_ROOT / "MC/BNB_cosmics/*-sel_all-wgts.df"),
    "data": str(SPRING_GEN1_ROOT / "data/BNB/_Fixed_all.df"),
    "intime": str(SPRING_GEN1_ROOT / "MC/intime/*_all.df"),
    "offbeam": str(SPRING_GEN1_ROOT / "data/OffBeam/*_all.df"),
    "dirt": str(SPRING_GEN1_ROOT / "MC/lowE/*_all.df"),
}

# -----------------------------------------------------------------------------
# Selected-events (final selection) bundles — used by selected_events*.py
# Keys: mc, data, intime, dirt
#
# These are produced by the "updated workflow" which writes timestamped directories
# like: 2026_05_01_102310__sel_mup-data-BNB_cosmics/sel_mup-data-BNB_cosmics_*.df
# -----------------------------------------------------------------------------
#
# By default we target the ``sel_mup`` campaign; override via:
#   export NUMUCC_SELECTED_EVENTS_TAG=sel_2prong
# (or any other selection tag that matches the directory naming convention).
SELECTED_EVENTS_TAG = os.environ.get("NUMUCC_SELECTED_EVENTS_TAG", "sel_mup")

SELECTED_EVENTS_GLOBS: Dict[str, str] = {
    "mc": str(SPRING_GEN1_ROOT / f"*__{SELECTED_EVENTS_TAG}-mc-BNB_cosmics/*.df"),
    "data": str(SPRING_GEN1_ROOT / f"*__{SELECTED_EVENTS_TAG}-data-BNB_cosmics/*.df"),
    # "intime" here means the off-beam light/offbeam sample used as the cosmics estimate
    "intime": str(SPRING_GEN1_ROOT / f"*__{SELECTED_EVENTS_TAG}-data-OffBeamLight/*.df"),
    # optional (some workflows don't produce it for selected-events plots)
    "dirt": str(SPRING_GEN1_ROOT / f"*__{SELECTED_EVENTS_TAG}-mc-dirt/*.df"),
}

# -----------------------------------------------------------------------------
# Multisim syst chunk inputs — **one glob per systematic** (often different dirs).
# Keys must match ``syst_multisim_common.NEUTRINO_SYST_ORDER``: MCstat, Flux, G4.
# Defaults repeat the same patterns as the old single-glob workflow; edit per-syst paths here.
# ``final``: tight-selection-style bundles; ``sel_all``: loose + wgts.
# -----------------------------------------------------------------------------
MULTISIM_SYST_GLOBS_FINAL: Dict[str, str] = {
    "MCstat": str(SPRING_GEN1_ROOT / "MC/BNB_cosmics/*.df"),
    "Flux": str(SPRING_GEN1_ROOT / "MC/BNB_cosmics/*.df"),
    "G4": str(SPRING_GEN1_ROOT / "MC/BNB_cosmics/*.df"),
}
MULTISIM_SYST_GLOBS_SEL_ALL: Dict[str, str] = {
    "MCstat": str(SPRING_GEN1_ROOT / "MC/BNB_cosmics/*-sel_all-wgts.df"),
    "Flux": str(SPRING_GEN1_ROOT / "MC/BNB_cosmics/*-sel_all-wgts.df"),
    "G4": str(SPRING_GEN1_ROOT / "MC/BNB_cosmics/*-sel_all-wgts.df"),
}

# For ``mc_df_stage == "final"``: drop paths whose basename contains this substring.
MULTISIM_FINAL_EXCLUDE_SUBSTRING = "-sel_all-wgts"

# Legacy single-glob names (union of per-syst globs); kept for scripts that need one pattern string.
MULTISIM_MC_GLOB_FINAL = MULTISIM_SYST_GLOBS_FINAL["Flux"]
MULTISIM_MC_GLOB_SEL_ALL = MULTISIM_SYST_GLOBS_SEL_ALL["Flux"]

# -----------------------------------------------------------------------------
# Detvar (WireMod + calo unisim) — typical chunk input dirs / globs
# -----------------------------------------------------------------------------
DETVAR_DF_GLOB_EXAMPLE_WIREMOD_YZ = str(
    SPRING_GEN1_ROOT / "2026_05_09_223419__sel_2prong-mc-BNB_cosmics-WireModYZ/*.df"
)
DETVAR_DF_GLOB_EXAMPLE_WIREMOD_XTXW = str(
    SPRING_GEN1_ROOT / "2026_05_09_223419__sel_2prong-mc-BNB_cosmics-WireModXTXW/*.df"
)

# -----------------------------------------------------------------------------
# DetVar chunked driver — ``(WireMod tag, glob)`` pairs for ``syst_detvar_chunk.py``
# -----------------------------------------------------------------------------
DETVAR_WIREMOD_GLOBS: List[Tuple[str, str]] = [
    ("wiremod_yz", DETVAR_DF_GLOB_EXAMPLE_WIREMOD_YZ),
    ("wiremod_xtxw", DETVAR_DF_GLOB_EXAMPLE_WIREMOD_XTXW),
]

# -----------------------------------------------------------------------------
# GENIE reweight samples — one glob per **group** (same layout as monolithic drivers).
# CCQE lives under ``genie_wgts-CCQE`` with ``*_geniewgts_CCQE.df`` filenames; other
# groups use ``genie_wgts-<Tag>/*.df``.
# -----------------------------------------------------------------------------
GENIE_GROUP_ORDER: Tuple[str, ...] = ("CCQE", "MEC", "RES", "DIS", "Other")

GENIE_GROUP_GLOBS: Dict[str, str] = {
    "CCQE": str(SPRING_GEN1_ROOT / "2026_05_10_235559__sel_mup-wgts_genie_CCQE/*.df"),
    "MEC": str(SPRING_GEN1_ROOT / "2026_05_10_235711__sel_mup-wgts_genie_MEC/*.df"),
    "RES": str(SPRING_GEN1_ROOT / "2026_05_10_235751__sel_mup-wgts_genie_RES/*.df"),
    # "nonRES": str(SPRING_GEN1_ROOT / "2026_05_10_235831__sel_mup-wgts_genie_nonRES/*.df"),
    "DIS": str(SPRING_GEN1_ROOT / "2026_05_10_235943__sel_mup-wgts_genie_DIS/*.df"),
    "Other": str(SPRING_GEN1_ROOT / "2026_05_11_000024__sel_mup-wgts_genie_Other/*.df"),
    # "Ar23p": str(SPRING_GEN1_ROOT / "MC/BNB_cosmics/genie_wgts-Ar23p/*.df"),
}


GENIE_GROUP_KNOBS: Dict[str, List[str]] = dict(
    zip(
        GENIE_GROUP_ORDER,
        [
            list(qe_genie_systematics),
            list(mec_genie_systematics),
            list(res_genie_systematics),
            list(dis_genie_systematics),
            list(other_genie_systematics),
        ],
    )
)


def iter_detvar_chunk_jobs(
    wiremod_globs: Optional[Sequence[Tuple[str, str]]] = None,
) -> Iterator[Tuple[str, str]]:
    """Yield ``(wiremod_tag, df_path)`` for ``syst_detvar_chunk.py`` (paths sorted per glob)."""
    pairs: Sequence[Tuple[str, str]] = wiremod_globs if wiremod_globs is not None else DETVAR_WIREMOD_GLOBS
    for tag, pattern in pairs:
        for p in sorted_glob(pattern):
            yield tag, p


def iter_genie_group_df_paths(
    genie_group: str,
    group_globs: Optional[Dict[str, str]] = None,
) -> Iterator[str]:
    """Yield sorted ``.df`` paths for one GENIE knob group (see ``GENIE_GROUP_GLOBS``)."""
    gmap = group_globs if group_globs is not None else GENIE_GROUP_GLOBS
    if genie_group not in gmap:
        raise KeyError(
            "unknown genie_group %r; expected one of %s" % (genie_group, tuple(gmap.keys()))
        )
    yield from sorted_glob(gmap[genie_group])


def iter_genie_chunk_map_tasks(
    group_globs: Optional[Dict[str, str]] = None,
) -> Iterator[Tuple[str, str]]:
    """Yield ``(genie_group_tag, df_path)`` for ``get_systematics_genie.py chunk-map``."""
    gmap = group_globs if group_globs is not None else GENIE_GROUP_GLOBS
    for tag in GENIE_GROUP_ORDER:
        if tag not in gmap:
            continue
        for p in sorted_glob(gmap[tag]):
            yield tag, p


def iter_cosmics_chunk_df_paths(sample: str) -> Iterator[str]:
    """Yield ``.df`` paths for cosmics chunk map (``sample`` is ``offbeam`` or ``intime``)."""
    if sample not in ("offbeam", "intime"):
        raise ValueError("sample must be 'offbeam' or 'intime', got %r" % sample)
    yield from iter_event_selection_df_paths(sample)


# -----------------------------------------------------------------------------
# Default work directories (override with env or shell ``WORK_BASE``)
# -----------------------------------------------------------------------------
def default_event_selection_work_root(tag: str | None = None) -> Path:
    from datetime import datetime

    t = tag or datetime.now().strftime("%Y%m%d")
    base = os.environ.get("NUMUCC_EVENT_SELECTION_WORK_BASE")
    if base:
        return Path(base)
    return Path(
        f"/exp/sbnd/data/users/{os.environ.get('USER', 'user')}/xsec/numucc_1p0pi/"
        f"event_selection-chunked-{t}"
    )


def default_genie_syst_work_root(tag: str | None = None) -> Path:
    from datetime import datetime

    t = tag or datetime.now().strftime("%Y%m%d")
    base = os.environ.get("NUMUCC_GENIE_SYST_WORK_BASE")
    if base:
        return Path(base)
    return Path(
        f"/exp/sbnd/data/users/{os.environ.get('USER', 'user')}/xsec/numucc_1p0pi/"
        f"genie_syst-chunked-{t}"
    )


def default_cosmics_syst_work_root(tag: str | None = None) -> Path:
    from datetime import datetime

    t = tag or datetime.now().strftime("%Y%m%d")
    base = os.environ.get("NUMUCC_COSMICS_SYST_WORK_BASE")
    if base:
        return Path(base)
    return Path(
        f"/exp/sbnd/data/users/{os.environ.get('USER', 'user')}/xsec/numucc_1p0pi/"
        f"cosmics_syst-chunked-{t}"
    )


def default_multisim_syst_work_root(tag: str | None = None) -> Path:
    from datetime import datetime

    t = tag or datetime.now().strftime("%Y%m%d")
    base = os.environ.get("NUMUCC_MULTISIM_SYST_WORK_BASE")
    if base:
        return Path(base)
    return Path(
        f"/exp/sbnd/data/users/{os.environ.get('USER', 'user')}/xsec/numucc_1p0pi/"
        f"multisim_syst-chunked-{t}"
    )


# -----------------------------------------------------------------------------
# Glob helpers
# -----------------------------------------------------------------------------
def sorted_glob(pattern: str) -> List[str]:
    paths = sorted(glob.glob(pattern))
    return [p for p in paths if Path(p).is_file()]


def iter_event_selection_df_paths(sample: str) -> Iterator[str]:
    if sample not in EVENT_SELECTION_GLOBS:
        raise KeyError("unknown sample %r; choose from %s" % (sample, tuple(EVENT_SELECTION_GLOBS)))
    for p in sorted_glob(EVENT_SELECTION_GLOBS[sample]):
        yield p


def _multisim_glob_map(mc_df_stage: str) -> Dict[str, str]:
    if mc_df_stage == "final":
        return MULTISIM_SYST_GLOBS_FINAL
    if mc_df_stage == "sel_all":
        return MULTISIM_SYST_GLOBS_SEL_ALL
    raise ValueError("mc_df_stage must be 'final' or 'sel_all', got %r" % mc_df_stage)


def iter_multisim_syst_df_paths(mc_df_stage: str, syst_name: str) -> Iterator[str]:
    """Yield ``.df`` paths for one neutrino systematic (its dedicated glob / directory)."""
    from analysis_village.numucc_1p0pi.syst_multisim_common import NEUTRINO_SYST_ORDER

    if syst_name not in NEUTRINO_SYST_ORDER:
        raise KeyError(
            "unknown syst_name %r; expected one of %s" % (syst_name, NEUTRINO_SYST_ORDER)
        )
    glob_map = _multisim_glob_map(mc_df_stage)
    pattern = glob_map[syst_name]
    if mc_df_stage == "final":
        excl = MULTISIM_FINAL_EXCLUDE_SUBSTRING
        for p in sorted_glob(pattern):
            if excl and excl in Path(p).name:
                continue
            yield p
        return
    yield from sorted_glob(pattern)


def iter_multisim_chunk_tasks(mc_df_stage: str) -> Iterator[Tuple[str, str]]:
    """Yield ``(syst_name, df_path)`` for the multisim **map** phase (one CAF shard per path).

    Legacy name retained for drivers. Prefer :func:`iter_multisim_map_tasks`.

    When every systematic shares the **same** glob string for this stage, yields
    ``("COMBINED", path)`` once per file so the driver runs a single HDF pass with
    all weights (same as legacy). Different globs emit one row per (syst, path).
    """
    from analysis_village.numucc_1p0pi.syst_multisim_common import NEUTRINO_SYST_ORDER

    glob_map = _multisim_glob_map(mc_df_stage)
    patterns = tuple(glob_map[sn] for sn in NEUTRINO_SYST_ORDER)
    if len(set(patterns)) == 1:
        sn0 = NEUTRINO_SYST_ORDER[0]
        for p in iter_multisim_syst_df_paths(mc_df_stage, sn0):
            yield "COMBINED", p
        return
    for sn in NEUTRINO_SYST_ORDER:
        for p in iter_multisim_syst_df_paths(mc_df_stage, sn):
            yield sn, p


def iter_multisim_mc_df_paths(mc_df_stage: str) -> Iterator[str]:
    """Yield unique ``.df`` paths across all systematics (legacy / convenience)."""
    seen = set()
    for _, p in iter_multisim_map_tasks(mc_df_stage):
        if p not in seen:
            seen.add(p)
            yield p


def iter_multisim_map_tasks(mc_df_stage: str) -> Iterator[Tuple[str, str]]:
    """Alias for :func:`iter_multisim_chunk_tasks` — clearer name (map shard, not exposure batch)."""
    yield from iter_multisim_chunk_tasks(mc_df_stage)


def summary_lines() -> Iterable[str]:
    """Human-readable listing for logs."""
    yield "# dataset_locations (numucc_1p0pi)"
    yield "SPRING_GEN1_ROOT=%s" % SPRING_GEN1_ROOT
    yield "PLOTS_BASE=%s" % PLOTS_BASE
    yield ""
    yield "## event_selection"
    for k, pat in EVENT_SELECTION_GLOBS.items():
        n = len(sorted_glob(pat))
        yield "  %s: %d file(s)  glob=%s" % (k, n, pat)
    yield ""
    yield "## multisim_mc (per-syst globs)"
    for stage, dmap in (("final", MULTISIM_SYST_GLOBS_FINAL), ("sel_all", MULTISIM_SYST_GLOBS_SEL_ALL)):
        yield "  [%s]" % stage
        for sn, pat in dmap.items():
            n = len(list(iter_multisim_syst_df_paths(stage, sn)))
            yield "    %s: %d file(s)  glob=%s" % (sn, n, pat)
    yield ""
    yield "## detvar (WireMod tags)"
    for tag, pat in DETVAR_WIREMOD_GLOBS:
        n = sum(1 for _ in sorted_glob(pat))
        yield "  %s: %d file(s)  glob=%s" % (tag, n, pat)
    yield ""
    yield "## genie (knob groups)"
    for tag, pat in GENIE_GROUP_GLOBS.items():
        n = len(list(iter_genie_group_df_paths(tag)))
        yield "  %s: %d file(s)  glob=%s" % (tag, n, pat)
    yield ""
    yield "## cosmics chunk inputs (offbeam / intime)"
    for sample in ("offbeam", "intime"):
        pat = EVENT_SELECTION_GLOBS[sample]
        n = len(sorted_glob(pat))
        yield "  %s: %d file(s)  glob=%s" % (sample, n, pat)


def print_summary() -> None:
    for line in summary_lines():
        print(line)


if __name__ == "__main__":
    print_summary()
