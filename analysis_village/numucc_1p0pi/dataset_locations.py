"""Central directory and glob registry for numu CC 1p0pi workflows.

Edit paths **here only** so drivers stay thin:

- ``run_event_selection_chunked.sh`` — map shards: MC/data/intime/offbeam/dirt ``.df`` files for
  ``event_selection_chunk.py`` / ``event_selection_aggregate.py``.
- ``run_syst_multisim_chunked.sh`` — per systematic + CAF shard (default syst subset: Flux+G4;
  MCstat opt-in via ``MULTISIM_SYST_TYPES`` / ``--syst-types``); **Combined** chunks live under
  ``multisim_syst-chunked-*`` / ``chunks/Combined/``; **MCstat**, **Flux**, and **G4** map shards use
  parallel ``mcstat_syst-chunked-*``, ``flux_syst-chunked-*``, and ``g4_syst-chunked-*`` (see
  :func:`default_mcstat_syst_work_root`, :func:`default_flux_syst_work_root`,
  :func:`default_g4_syst_work_root`).
- ``syst_multisim_aggregate.py`` — merge neutrino multisim map outputs into ``MCstat/``, ``Flux/``, ``G4/``.
- ``syst_detvar_chunk.py`` / ``syst_detvar_aggregate.py`` — WireMod + calo variants;
  input globs are listed in ``DETVAR_WIREMOD_GLOBS`` / :func:`iter_detvar_chunk_jobs`.
- ``get_systematics_genie.py`` / ``run_syst_genie_chunked.sh`` — GENIE reweights: one glob per
  knob **group** (``GENIE_GROUP_GLOBS``) / :func:`iter_genie_chunk_map_tasks`; ``syst_genie_aggregate.py``
  publishes ``GENIE/cov_mat_dict.pkl`` on the syst disk (phase 3 of the shell driver).
- ``run_syst_cosmics_chunked.sh`` — ``syst_cosmics_chunk.py`` / ``syst_cosmics_aggregate.py``;
  globs reuse ``EVENT_SELECTION_GLOBS`` ``offbeam`` / ``intime``.
- ``default_syst_disk_root()`` — unified ``syst_disk_layout`` root (``Cosmics/``, ``MCstat/``,
  ``Detector/``, …) used by the ``run_syst_*`` drivers unless ``NUMUCC_SYST_DISK_ROOT`` is set
  or a script passes an explicit override.

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

SPRING_GEN1_ROOT_EAF = Path(
    os.environ.get(
        "NUMUCC_SPRING_GEN1_ROOT",
        "/scratch/7DayLifetime/munjung/xsec",
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
    # "mc": str(SPRING_GEN1_ROOT / "2026_05_11_041007__sel_all-mc-BNB_cosmics/*df"),
    #"mc": str(SPRING_GEN1_ROOT / "2026_05_11_183347__sel_all-mc-BNB_cosmics-EField_R00/*df"),
    "mc": str(SPRING_GEN1_ROOT / "2026_05_11_041007__sel_all-mc-BNB_cosmics/*.df"),
    "data": str(SPRING_GEN1_ROOT / "2026_05_16_230859__sel_all-data-1e20/*.df"),
    "intime": str(SPRING_GEN1_ROOT / "2026_05_11_040132__sel_all-mc-Intime/*.df"),
    "offbeam": str(SPRING_GEN1_ROOT / "2026_05_11_035756__sel_all-data-OffBeamLight/*.df"),
    "dirt": str(SPRING_GEN1_ROOT / "2026_05_11_040638__sel_all-mc-dirt/*df"),
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
    "data": str(SPRING_GEN1_ROOT / f"*__{SELECTED_EVENTS_TAG}-data-Gen1/*.df"),
    "intime": str(SPRING_GEN1_ROOT / f"*__{SELECTED_EVENTS_TAG}-mc-Intime/*.df"),
    "offbeam": str(SPRING_GEN1_ROOT / f"*__{SELECTED_EVENTS_TAG}-data-OffBeamLight/*.df"),
    "dirt": str(SPRING_GEN1_ROOT / f"*__{SELECTED_EVENTS_TAG}-mc-dirt/*.df"),
}

# -----------------------------------------------------------------------------
# Multisim syst chunk inputs — **one glob per systematic** (often different dirs).
# Keys must match ``syst_multisim_common.NEUTRINO_SYST_ORDER``: MCstat, Flux, G4.
# Defaults repeat the same patterns as the old single-glob workflow; edit per-syst paths here.
# ``final``: tight-selection-style bundles; ``sel_all``: loose + wgts.
# -----------------------------------------------------------------------------
MULTISIM_SYST_GLOBS_FINAL: Dict[str, str] = {
    "MCstat": str(SPRING_GEN1_ROOT / "2026_05_18_145611__sel_mup-wgts_mcstat/merged_perTPC/*.df"),
    "Flux": str(SPRING_GEN1_ROOT / "2026_05_11_155745__sel_mup-wgts_flux/merged_perTPC/*.df"),
    # "Flux": str(SPRING_GEN1_ROOT_EAF / "2026_05_11_155745__sel_mup-wgts_flux/*.df"),
    "G4": str(SPRING_GEN1_ROOT / "2026_05_11_031351__sel_mup-wgts_g4/merged_perTPC/*.df"),
    # "G4": str(SPRING_GEN1_ROOT_EAF / "2026_05_11_031351__sel_mup-wgts_g4/*.df"),
}
MULTISIM_SYST_GLOBS_SEL_ALL: Dict[str, str] = {
    # "MCstat": str(SPRING_GEN1_ROOT / "2026_05_11_084007__sel_mup-wgts_mcstat/*.df"),
    # "Flux": str(SPRING_GEN1_ROOT / "2026_05_11_031846__sel_mup-wgts_flux/*.df"),
    # "G4": str(SPRING_GEN1_ROOT / "2026_05_11_031351__sel_mup-wgts_g4/*.df"),
}

# Legacy single-glob names (union of per-syst globs); kept for scripts that need one pattern string.
# MULTISIM_MC_GLOB_FINAL = MULTISIM_SYST_GLOBS_FINAL["Flux"]
# MULTISIM_MC_GLOB_SEL_ALL = MULTISIM_SYST_GLOBS_SEL_ALL["Flux"]

# -----------------------------------------------------------------------------
# Detvar (WireMod + calo unisim) — typical chunk input dirs / globs
# -----------------------------------------------------------------------------
DETVAR_DF_GLOB_EXAMPLE_WIREMOD_YZ = str(
    SPRING_GEN1_ROOT / "2026_05_09_223419__sel_2prong-mc-BNB_cosmics-WireModYZ/*.df"
)
DETVAR_DF_GLOB_EXAMPLE_WIREMOD_XTXW = str(
    SPRING_GEN1_ROOT / "2026_05_11_103733__sel_2prong-mc-BNB_cosmics-WireModXTXW/*df"
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
GENIE_GROUP_ORDER: Tuple[str, ...] = ("CCQE", "MEC", "RES", "nonRES", "DIS", "Other", "Ar23p")

# Final-selection-style GENIE weight bundles (same convention as ``MULTISIM_SYST_GLOBS_FINAL``).
GENIE_GROUP_GLOBS: Dict[str, str] = {
    "CCQE": str(SPRING_GEN1_ROOT / "2026_05_11_024530__sel_mup-wgts_genie_CCQE/merged_perTPC/*.df"),
    "MEC": str(SPRING_GEN1_ROOT / "2026_05_11_030314__sel_mup-wgts_genie_MEC/merged_perTPC/*.df"),
    "RES": str(SPRING_GEN1_ROOT / "2026_05_11_030547__sel_mup-wgts_genie_RES/merged_perTPC/*.df"),
    "nonRES": str(SPRING_GEN1_ROOT / "2026_05_11_030906__sel_mup-wgts_genie_nonRES/merged_perTPC/*.df"),
    "DIS": str(SPRING_GEN1_ROOT / "2026_05_11_031206__sel_mup-wgts_genie_DIS/merged_perTPC/*.df"),
    "Other": str(SPRING_GEN1_ROOT / "2026_05_11_031520__sel_mup-wgts_genie_Other/merged_perTPC/*.df"),
    "Ar23p": str(SPRING_GEN1_ROOT / "2026_05_12_010953__sel_mup-wgts_genie_Ar23p/merged_perTPC/*.df"),
}

# Loose ``sel_all``-style MC + GENIE weights (evt / trk / hdr / mcnu). Fill when running
# ``get_systematics_genie.py chunk-map --input-stage sel_all``; empty groups are skipped.
GENIE_GROUP_GLOBS_SEL_ALL: Dict[str, str] = {}


GENIE_GROUP_KNOBS: Dict[str, List[str]] = dict(
    zip(
        GENIE_GROUP_ORDER,
        [
            list(qe_genie_systematics),
            list(mec_genie_systematics),
            list(res_genie_systematics),
            list(nonres_genie_systematics),
            list(dis_genie_systematics),
            list(other_genie_systematics),
            list(ar23p_genie_systematics),
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


def _genie_glob_map(mc_df_stage: str) -> Dict[str, str]:
    if mc_df_stage == "final":
        return GENIE_GROUP_GLOBS
    if mc_df_stage == "sel_all":
        return GENIE_GROUP_GLOBS_SEL_ALL
    raise ValueError("mc_df_stage must be 'final' or 'sel_all', got %r" % mc_df_stage)


def iter_genie_group_df_paths(
    genie_group: str,
    group_globs: Optional[Dict[str, str]] = None,
    *,
    mc_df_stage: str = "final",
) -> Iterator[str]:
    """Yield sorted ``.df`` paths for one GENIE knob group.

    ``mc_df_stage`` selects :data:`GENIE_GROUP_GLOBS` vs :data:`GENIE_GROUP_GLOBS_SEL_ALL`
    when ``group_globs`` is omitted.
    """
    gmap = group_globs if group_globs is not None else _genie_glob_map(mc_df_stage)
    if genie_group not in gmap:
        raise KeyError(
            "unknown genie_group %r; expected one of %s" % (genie_group, tuple(gmap.keys()))
        )
    yield from sorted_glob(gmap[genie_group])


def iter_genie_chunk_map_tasks(
    group_globs: Optional[Dict[str, str]] = None,
    mc_df_stage: str = "final",
) -> Iterator[Tuple[str, str]]:
    """Yield ``(genie_group_tag, df_path)`` for ``get_systematics_genie.py chunk-map``.

    ``mc_df_stage`` (``final`` | ``sel_all``) picks the default glob map; pass
    ``group_globs`` explicitly to override.

    Tags are every key present in ``gmap``, ordered by :data:`GENIE_GROUP_ORDER` first,
    then any remaining keys (sorted) so ad-hoc entries in ``GENIE_GROUP_GLOBS`` are included.
    """
    gmap = group_globs if group_globs is not None else _genie_glob_map(mc_df_stage)
    seen: set[str] = set()
    for tag in GENIE_GROUP_ORDER:
        if tag not in gmap:
            continue
        seen.add(tag)
        for p in sorted_glob(gmap[tag]):
            yield tag, p
    for tag in sorted(k for k in gmap if k not in seen):
        for p in sorted_glob(gmap[tag]):
            yield tag, p


def iter_cosmics_chunk_df_paths(sample: str, input_stage: str = "sel_all") -> Iterator[str]:
    """Yield ``.df`` paths for cosmics chunk map.

    ``sample`` is ``offbeam`` or ``intime``. ``input_stage`` selects the input glob:

    * ``sel_all`` (default): :data:`EVENT_SELECTION_GLOBS` — raw evt/trk/hdr dfs.
    * ``final``: :data:`SELECTED_EVENTS_GLOBS` — already-final-selected dfs.
    """
    if sample not in ("offbeam", "intime"):
        raise ValueError("sample must be 'offbeam' or 'intime', got %r" % sample)
    if input_stage == "sel_all":
        yield from iter_event_selection_df_paths(sample)
        return
    if input_stage == "final":
        if sample not in SELECTED_EVENTS_GLOBS:
            raise KeyError("sample %r not in SELECTED_EVENTS_GLOBS" % sample)
        for p in sorted_glob(SELECTED_EVENTS_GLOBS[sample]):
            yield p
        return
    raise ValueError(
        "input_stage must be 'sel_all' or 'final', got %r" % input_stage
    )


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


def default_g4_syst_work_root(tag: str | None = None) -> Path:
    """Default map-shard root for G4-only neutrino multisim chunks (parallel to multisim work dir)."""
    from datetime import datetime

    t = tag or datetime.now().strftime("%Y%m%d")
    base = os.environ.get("NUMUCC_G4_SYST_WORK_BASE")
    if base:
        return Path(base)
    return Path(
        f"/exp/sbnd/data/users/{os.environ.get('USER', 'user')}/xsec/numucc_1p0pi/"
        f"g4_syst-chunked-{t}"
    )


def default_flux_syst_work_root(tag: str | None = None) -> Path:
    """Default map-shard root for Flux-only neutrino multisim chunks (parallel to multisim work dir)."""
    from datetime import datetime

    t = tag or datetime.now().strftime("%Y%m%d")
    base = os.environ.get("NUMUCC_FLUX_SYST_WORK_BASE")
    if base:
        return Path(base)
    return Path(
        f"/exp/sbnd/data/users/{os.environ.get('USER', 'user')}/xsec/numucc_1p0pi/"
        f"flux_syst-chunked-{t}"
    )


def default_mcstat_syst_work_root(tag: str | None = None) -> Path:
    """Default map-shard root for MCstat-only neutrino multisim chunks (parallel to multisim work dir)."""
    from datetime import datetime

    t = tag or datetime.now().strftime("%Y%m%d")
    base = os.environ.get("NUMUCC_MCSTAT_SYST_WORK_BASE")
    if base:
        return Path(base)
    return Path(
        f"/exp/sbnd/data/users/{os.environ.get('USER', 'user')}/xsec/numucc_1p0pi/"
        f"mcstat_syst-chunked-{t}"
    )


def default_detvar_syst_work_root(tag: str | None = None) -> Path:
    """Scratch/output root for chunked detvar map pickles (``chunks/`` under here)."""
    from datetime import datetime

    t = tag or datetime.now().strftime("%Y%m%d")
    base = os.environ.get("NUMUCC_DETVAR_SYST_WORK_BASE")
    if base:
        return Path(base)
    return Path(
        f"/exp/sbnd/data/users/{os.environ.get('USER', 'user')}/xsec/numucc_1p0pi/"
        f"detvar_systematics-{t}"
    )


def default_syst_disk_root() -> Path:
    """Default root for the unified ``syst_disk_layout`` tree (``Cosmics/``, ``MCstat/``, …).

    Same logical tree that ``utils.get_syst_unc`` reads when ``NUMUCC_SYST_DISK_ROOT`` is set.
    If that environment variable is set, this function returns that path (expanded). If not,
    returns a stable per-user default so ``run_syst_*`` scripts can aggregate without extra args.
    """
    env = os.environ.get("NUMUCC_SYST_DISK_ROOT")
    if env:
        return Path(env).expanduser()
    return Path(
        f"/exp/sbnd/data/users/{os.environ.get('USER', 'user')}/xsec/numucc_1p0pi/syst_disk"
    )


def default_syst_disk_cc_root() -> Path:
    """Default root for **joint** (cross-variable) syst outputs (``syst_disk_CC`` tree).

    Set ``NUMUCC_SYST_DISK_CC_ROOT`` to override. Otherwise uses ``syst_disk_CC`` as a sibling
    directory next to :func:`default_syst_disk_root` when that path ends with ``syst_disk``,
    else ``<parent>/syst_disk_CC`` alongside the same parent as ``default_syst_disk_root``.
    """
    env = os.environ.get("NUMUCC_SYST_DISK_CC_ROOT")
    if env:
        return Path(env).expanduser()
    base = default_syst_disk_root()
    name = base.name
    if name == "syst_disk":
        return base.parent / "syst_disk_CC"
    return base.parent / "syst_disk_CC"


def default_joint_genie_cc_work_root(tag: str | None = None) -> Path:
    """Default map-shard root for joint (cross-variable) GENIE CC chunks."""
    from datetime import datetime

    t = tag or datetime.now().strftime("%Y%m%d")
    base = os.environ.get("NUMUCC_JOINT_GENIE_CC_WORK_BASE")
    if base:
        return Path(base)
    return Path(
        f"/exp/sbnd/data/users/{os.environ.get('USER', 'user')}/xsec/numucc_1p0pi/"
        f"joint_genie_cc-chunked-{t}"
    )


def default_joint_multisim_cc_work_root(tag: str | None = None) -> Path:
    """Default map-shard root for joint (cross-variable) multisim CC chunks."""
    from datetime import datetime

    t = tag or datetime.now().strftime("%Y%m%d")
    base = os.environ.get("NUMUCC_JOINT_MULTISIM_CC_WORK_BASE")
    if base:
        return Path(base)
    return Path(
        f"/exp/sbnd/data/users/{os.environ.get('USER', 'user')}/xsec/numucc_1p0pi/"
        f"joint_multisim_cc-chunked-{t}"
    )


# -----------------------------------------------------------------------------
# Glob helpers
# -----------------------------------------------------------------------------
def sorted_glob(pattern: str) -> List[str]:
    """Sorted ``glob.glob`` hits that exist and are not directories (after ``realpath``).

    Avoids ``Path.is_file()``, which can be false on some PNFS/dCache shards.
    """
    out: List[str] = []
    for p in sorted(glob.glob(pattern)):
        if not os.path.exists(p):
            continue
        try:
            resolved = os.path.realpath(p)
        except OSError:
            continue
        if os.path.isdir(resolved):
            continue
        out.append(p)
    return out


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
    yield "default_syst_disk_root=%s" % default_syst_disk_root()
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
