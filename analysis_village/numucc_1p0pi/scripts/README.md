# `numucc_1p0pi/scripts`

Analysis drivers for systematic uncertainties, event selection, and related workflows. Below focuses on **multisim** (MCstat, flux, Geant4 reinteraction weights) and how the Python entry points relate to each other.

---

## Multisim (MCstat / Flux / G4)

### `get_systematics_multisim.py` — primary monolithic driver

Use this for **new work** when the full MC sample fits in memory.

- Builds universe histograms for **total selected event rate** and **background-subtracted signal rate** via `analysis_village.numucc_1p0pi.utils.get_univ_rates`.
- Writes covariance / fractional covariance / correlation matrices and optional plots.
- Outputs NPZs consumed downstream (including layouts compatible with `utils.get_syst_unc` when keyed appropriately):
  - `mcstat_syst_dict.npz`
  - `flux_syst_dict.npz`
  - `g4_syst_dict.npz`
  - `multisim_syst_dict_all.npz` (combined convenience bundle)

Supports bundled `mc.Flux` / `mc.G4` blocks or **per-knob** flux/G4 modes via CLI (`--flux-mode`, `--g4-mode`, etc.). See the script docstring for examples.

**Supersedes** the removed standalone scripts `get_systematics-flux.py` and `get_systematics-G4.py`, which duplicated this logic and used obsolete call signatures.

---

## Chunked multisim (low memory / grid jobs)

Use this path when CAFs are split across many `.df` HDF files or HDF splits inside files.

| File | Role |
|------|------|
| **`syst_multisim_chunk.py`** | **Map phase:** reads one `.df` sequentially (`evt_0`, `evt_1`, …), accumulates summed `univ_events` / `cv_events` per systematic and variable, writes **`nu__*.pkl`** (or `nu__<syst>__<stem>.pkl` when `--syst-names` is a subset). |
| **`syst_multisim_aggregate.py`** | **Reduce phase:** merges all `nu__*.pkl` under a directory, builds covariances, writes `MCstat/`, `Flux/`, `G4/` under `--syst-disk-root` (neutrino multisim only; cosmics use `run_syst_cosmics_chunked.sh`). |

### `--input-stage` (`final` vs `sel_all`)

Both `syst_multisim_chunk.py` and `syst_cosmics_chunk.py` accept
`--input-stage {final,sel_all}`:

- `final` *(default)*: the `.df` is at the final-selection level (e.g. `SELECTED_EVENTS_GLOBS` / `MULTISIM_SYST_GLOBS_FINAL`). The chunk reads only the `evt_{i}` table and histograms the final-selected variables. Multisim uses the legacy `get_univ_rates` (signal + background-subtracted rate covariance).
- `sel_all`: the `.df` is a raw sel_all bundle (`EVENT_SELECTION_GLOBS` / `MULTISIM_SYST_GLOBS_SEL_ALL`) carrying `evt+trk+hdr`. The chunk re-runs the full numuCC 1p0pi event-selection pipeline (`event_selection_pipeline_def.build_pipeline`) on every split and accumulates histograms at every **cut stage** (`nu_score`, `n_trks`, `track_score`, `vtx_dist`, `trk_len`, `mcs_range_diff`, `chi2_mu`, `chi2_p`) *and* at the **final stage** for the final-selected variables.

The aggregate scripts auto-detect the layout from the chunk pickles. Cut-stage and final-stage `var_save_name`s are disjoint, so the output NPZs (`cosmics_syst_dict.npz`, `mcstat_syst_dict.npz`, `flux_syst_dict.npz`, `g4_syst_dict.npz`) keep their existing flat layout while gaining the cut-variable covariances.

Shell convenience: **`run_syst_multisim_chunked.sh`** — loops chunk tasks from `dataset_locations` (`MC_DF_STAGE={final,sel_all}` drives both the input glob *and* `--input-stage` on the chunk); **`run_syst_cosmics_chunked.sh`** — same, controlled by `INPUT_STAGE={final,sel_all}` (default `sel_all`).

Shared helpers live in **`../syst_multisim_common.py`** (`build_var_configs`, `drop_bad_g4_weights`, `save_neutrino_multisim_npzs`, `save_cosmics_legacy_npz`, …) and **`../syst_disk_layout.py`**.

---

## Legacy orchestrator

### `get_systematics_mcstat_flux_g4.py`

Older **single-metric** workflow (`bkgd_subtract=True` only for neutrino multisim), optional **cosmics** folded into the same `syst_dict`, and a shortcut **`--chunks-dir`** that invokes `syst_multisim_aggregate.py` instead of loading all MC in RAM.

Prefer **`get_systematics_multisim.py`** for neutrino-only multisim with two metrics; prefer **`get_systematics_cosmics.py`** for cosmic unisim alone; keep this script when you need exact parity with older notebooks or the combined legacy pipeline.

---

## Cosmics

### `get_systematics_cosmics.py`

Cosmic background unisim (e.g. off-beam vs intime MC), separate from MCstat/Flux/G4 multisim. Requires **`--syst-disk-root`** or **`NUMUCC_SYST_DISK_ROOT`**; writes **`Cosmics/`** under that root.

---

## GENIE multisim

### `get_systematics_genie.py`

GENIE knob uncertainties require **`evt`** and **`mcnu`** dataframes and separate **rate** vs **cross-section** covariance logic (see script docstring). Supports monolithic use or **`chunk-map`** / **`chunk-merge`** for HDF-split processing.

---

## Quick reference: which script when?

| Goal | Script(s) |
|------|-----------|
| Monolithic MCstat + Flux + G4, two rate metrics, modern NPZs | `get_systematics_multisim.py` |
| Chunked neutrino multisim → legacy NPZs + plots | `syst_multisim_chunk.py` → `syst_multisim_aggregate.py` (or `run_syst_multisim_chunked.sh`) |
| Legacy single-metric + cosmics or `--chunks-dir` wrapper | `get_systematics_mcstat_flux_g4.py` |
| Cosmics only | `get_systematics_cosmics.py` |
| GENIE rate + xsec (with `mcnu`) | `get_systematics_genie.py` |

---

## Related code

- **`utils.get_syst_unc`** — reads the **full** `syst_disk_layout` tree (`NUMUCC_SYST_DISK_ROOT`) and aborts if any expected file is missing.
- **`event_selection_aggregate.py`** — chunked selection reduce; passes `--syst-disk-root` into `get_syst_unc` when chunk pickles omit universe histograms.
