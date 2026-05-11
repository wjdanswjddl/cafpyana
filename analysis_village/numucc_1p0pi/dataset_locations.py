"""Central directory and glob registry for numu CC 1p0pi workflows.

Edit paths **here only** so drivers stay thin:

- ``run_event_selection_chunked.sh`` — map shards: MC/data/intime/offbeam/dirt ``.df`` files for
  ``event_selection_chunk.py`` / ``event_selection_aggregate.py``.
- ``run_syst_multisim_chunked.sh`` — per systematic + CAF shard (each syst may use its own glob).
- ``syst_multisim_aggregate.py`` — merge map outputs into a ``syst_disk_layout`` tree.
- ``syst_detvar_chunk.py`` / ``syst_detvar_aggregate.py`` — WireMod + calo variants;
  input globs are listed in ``DETVAR_WIREMOD_GLOBS`` / :func:`iter_detvar_chunk_jobs`.

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

# -----------------------------------------------------------------------------
# Base release (Spring Gen 1 — edit for other campaigns)
# -----------------------------------------------------------------------------
SPRING_GEN1_ROOT = Path(
    os.environ.get(
        "NUMUCC_SPRING_GEN1_ROOT",
        "/exp/sbnd/data/users/munjung/xsec/2025Spring_v10_06_00_09",
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
    SPRING_GEN1_ROOT / "MC/BNB_cosmics/wiremod_yz/*.df"
)
DETVAR_DF_GLOB_EXAMPLE_WIREMOD_XTXW = str(
    SPRING_GEN1_ROOT / "MC/BNB_cosmics/wiremod_xtxw/*.df"
)

# -----------------------------------------------------------------------------
# DetVar chunked driver — ``(WireMod tag, glob)`` pairs for ``syst_detvar_chunk.py``
# -----------------------------------------------------------------------------
DETVAR_WIREMOD_GLOBS: List[Tuple[str, str]] = [
    ("wiremod_yz", DETVAR_DF_GLOB_EXAMPLE_WIREMOD_YZ),
    ("wiremod_xtxw", DETVAR_DF_GLOB_EXAMPLE_WIREMOD_XTXW),
]


def iter_detvar_chunk_jobs(
    wiremod_globs: Optional[Sequence[Tuple[str, str]]] = None,
) -> Iterator[Tuple[str, str]]:
    """Yield ``(wiremod_tag, df_path)`` for ``syst_detvar_chunk.py`` (paths sorted per glob)."""
    pairs: Sequence[Tuple[str, str]] = wiremod_globs if wiremod_globs is not None else DETVAR_WIREMOD_GLOBS
    for tag, pattern in pairs:
        for p in sorted_glob(pattern):
            yield tag, p


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


def print_summary() -> None:
    for line in summary_lines():
        print(line)


if __name__ == "__main__":
    print_summary()
