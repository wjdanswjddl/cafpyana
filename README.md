# cafpyana

A set of scripts for analyzing SBN CAF files with the python.
Goal is to provide an easy starting point for SBN analysis to everyone.
For more details and instructions, please check [wiki](https://github.com/sungbinoh/cafpyana/wiki).

## Releases
Please check tags for releases!

## Version compatibility
Current `main` branch is based on python v3.9.15.
It is for running the repository without any issue at gpvm servers with the `spack`.
There is no need to open an SL7 image.
If python version is updated in gpvm servers, compatibility issue should be revisited.

---

## νμ CC 1p0π analysis workflow (`numucc_1p0pi`)

This section describes the **structural** pipeline around `analysis_village/numucc_1p0pi/`.
Input directories are centralized in `analysis_village/numucc_1p0pi/dataset_locations.py`
(override `NUMUCC_SPRING_GEN1_ROOT` or edit `SPRING_GEN1_ROOT` as needed). Naming:
**map shard** = one production `.df` file; **exposure batch** = time-ordered slice of data for staged access (`exposure_access.py`).

### Phase 1 — Systematic covariance outputs

Outputs are lightweight **compressed NPZs** (`*_syst_dict.npz`) plus diagnostic plots
from `syst_multisim_aggregate.py`, and a machine-readable **`covariance_manifest.json`**
(variable list and map-shard metadata).

| Source | Role | Entry script |
| --- | --- | --- |
| MC statistics | Poisson-style universe weights (`makedf/mcstat.py`) | Folded into Flux/G4 map via `syst_multisim_chunk.py` or dedicated runners |
| Flux & G4 | Multisim weights from MC CAFs | `scripts/run_syst_multisim_chunked.sh` → `syst_multisim_aggregate.py` |
| GENIE | Multisim / unisim knobs | `scripts/get_systematics_genie.py` (`chunk-map` / `chunk-merge` isolates rate vs **xsec** paths in code) |
| Cosmics | Offbeam CV vs intime unisim (`get_systematics_cosmics.py`) | Optional inside `syst_multisim_aggregate.py`, or standalone script |
| Detector | DetVar CAFs (WireMod + calo) | `scripts/syst_detvar_chunk.py` → `syst_detvar_aggregate.py` |

**Example — Flux/G4/MCstat map + aggregate:**

```bash
export PYTHONPATH="/path/to/cafpyana${PYTHONPATH:+:$PYTHONPATH}"
cd /path/to/cafpyana/analysis_village/numucc_1p0pi/scripts
export NUMUCC_SPRING_GEN1_ROOT="/exp/sbnd/data/users/<you>/xsec/2025Spring_v10_06_00_09"   # optional
./run_syst_multisim_chunked.sh
# Outputs under WORK_BASE (default dated tree): syst_disk_layout folders MCstat/, Flux/, G4/,
# optional Cosmics/, plus covariance_manifest.json at WORK_BASE.
```

**Example — GENIE chunk map/merge (after editing knob lists / paths as needed):**

```bash
export PYTHONPATH="/path/to/cafpyana${PYTHONPATH:+:$PYTHONPATH}"
python get_systematics_genie.py chunk-map --df-file /path/to/MC/file.df --out-dir ./genie_chunks
python get_systematics_genie.py chunk-merge --chunks-dir ./genie_chunks --out-dir ./genie_cov --npz ./genie_syst_dict.npz
```

**Example — cosmics-only NPZ:**

```bash
export NUMUCC_SYST_DISK_ROOT=/path/to/syst_disk
python get_systematics_cosmics.py
```

Point downstream plotting at precomputed covariances with **`NUMUCC_SYST_DISK_ROOT`**
(one directory containing `MCstat/`, `Flux/`, `G4/`, `GENIE/`, `Cosmics/`, `Detector/` —
see `analysis_village/numucc_1p0pi/syst_disk_layout.py`) or pass
`event_selection_aggregate.py --syst-disk-root <dir>`.

### Phase 2 — Event selection, overlays, and χ²

1. **Map/reduce histograms:** `run_event_selection_chunked.sh` runs `event_selection_chunk.py`
   per sample shard, then `event_selection_aggregate.py` merges and plots.
2. **Asymmetric data error bars on plots** use `pyanalib.stat_helpers.return_data_stat_err`
   (gamma-based central intervals for count-like bins), wired through
   `utils.overlay_hists_from_histdata` / `overlay_histdata`.
3. **χ²** uses `pyanalib.stat_helpers.get_chi2` with a **combined** covariance (systematic +
   diagonal data term derived from those asymmetric intervals).

```bash
export PYTHONPATH="/path/to/cafpyana${PYTHONPATH:+:$PYTHONPATH}"
cd /path/to/cafpyana/analysis_village/numucc_1p0pi/scripts
export WORK_BASE="/exp/sbnd/data/users/$(whoami)/xsec/numucc_1p0pi/event_selection-$(date +%Y%m%d)"
./run_event_selection_chunked.sh
```

Optional: pass `--syst-disk-root` to `event_selection_aggregate.py` if MC-universe bands
are not embedded in the pickles (``utils.get_syst_unc`` loads **all** category files from that
tree and **fails loudly** if any are missing).

Producer scripts (`syst_multisim_aggregate.py`, `syst_detvar_aggregate.py`, cosmics/GENIE drivers)
each write into their subdirectory under the same root.

**Integrated smoke test** (multisim + DetVar + event selection, capped file counts):

```bash
python analysis_village/numucc_1p0pi/scripts/run_workflow_test.py -o /path/to/workflow_out --max-files 2
```

### Phase 3 — Staged data access (exposure batches)

Policy stages `DataAccessStage` in `exposure_access.py`:

| Stage | Meaning | Typical driver |
| --- | --- | --- |
| 1 | Fixed Dev sample only | Point `EVENT_SELECTION_GLOBS["data"]` / files config at dev sample |
| 2 | Gen 1 independent exposure batches | `selected_events.py --n_time_splits N --exposure-batch-index K` |
| 3 | Gen 1 cumulative batches | `selected_events_cumulative.py --n_time_splits N [--exposure-batch-indices ...]` |
| 4 | Full Gen 1 | Single batch spanning all data (`n_time_splits=1`) or full concat policy |

**Examples:**

```bash
# Stage 2 — one batch (legacy flag still works)
python selected_events.py --n_time_splits 15 --exposure-batch-index 3

# Stage 3 — cumulative indices 0,1,2 only
python selected_events_cumulative.py --n_time_splits 15 --exposure-batch-indices 0 1 2
```

### Artifact reference

- `syst_multisim_aggregate.py` output: under `--syst-disk-root`, writes `MCstat/`, `Flux/`, `G4/`,
  optional `Cosmics/`, plus `covariance_manifest.json` at the root of that tree.
- Chunked map pickles: `nu__*.pkl` (multisim), `mc__*.pkl` / `data__*.pkl` (event selection).