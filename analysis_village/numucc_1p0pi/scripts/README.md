# `numucc_1p0pi/scripts`

Analysis drivers for systematic uncertainties, event selection, and related workflows.
All production paths are **map/reduce** ("chunked") so they run within memory limits on
large samples; monolithic legacy drivers have been retired.

---

## Multisim (MCstat / Flux / G4)

Chunked map/reduce over `.df` HDF files (or HDF splits inside files):

| File | Role |
|------|------|
| **`syst_multisim_chunk.py`** | **Map phase:** reads one `.df` sequentially (`evt_0`, `evt_1`, …), accumulates summed `univ_events` / `cv_events` per systematic and variable, writes **`nu__*.pkl`** (or `nu__<syst>__<stem>.pkl` when `--syst-names` is a subset). |
| **`syst_multisim_parallel.py`** | Runs many chunk tasks in parallel workers. |
| **`syst_multisim_aggregate.py`** | **Reduce phase:** merges all `nu__*.pkl` under a directory, builds covariances, writes `MCstat/`, `Flux/`, `G4/` under `--syst-disk-root` (neutrino multisim only; cosmics use `run_syst_cosmics_chunked.sh`). |

### `--input-stage` (`final` vs `sel_all`)

Both `syst_multisim_chunk.py` and `syst_cosmics_chunk.py` accept
`--input-stage {final,sel_all}`:

- `final` *(default)*: the `.df` is at the final-selection level (e.g. `SELECTED_EVENTS_GLOBS` / `MULTISIM_SYST_GLOBS_FINAL`). The chunk reads only the `evt_{i}` table and histograms the final-selected variables.
- `sel_all`: the `.df` is a raw sel_all bundle (`EVENT_SELECTION_GLOBS` / `MULTISIM_SYST_GLOBS_SEL_ALL`) carrying `evt+trk+hdr`. The chunk re-runs the full numuCC 1p0pi event-selection pipeline (`event_selection_pipeline_def.build_pipeline`) on every split and accumulates histograms at every **cut stage** (`nu_score`, `n_trks`, `track_score`, `vtx_dist`, `trk_len`, `mcs_range_diff`, `chi2_mu`, `chi2_p`) *and* at the **final stage** for the final-selected variables.

The aggregate scripts auto-detect the layout from the chunk pickles. Cut-stage and final-stage `var_save_name`s are disjoint, so the output NPZs (`cosmics_syst_dict.npz`, `mcstat_syst_dict.npz`, `flux_syst_dict.npz`, `g4_syst_dict.npz`) keep their existing flat layout while gaining the cut-variable covariances.

Shell convenience: **`run_syst_multisim_chunked.sh`** — loops chunk tasks from `dataset_locations` (`MC_DF_STAGE={final,sel_all}` drives both the input glob *and* `--input-stage` on the chunk); **`run_syst_cosmics_chunked.sh`** — same, controlled by `INPUT_STAGE={final,sel_all}` (default `sel_all`).

Shared helpers live in **`../syst_multisim_common.py`** (`build_var_configs`, `drop_bad_g4_weights`, `save_neutrino_multisim_npzs`, `save_cosmics_legacy_npz`, …) and **`../syst_disk_layout.py`**.

---

## Cosmics

- **`syst_cosmics_chunk.py`** → **`syst_cosmics_aggregate.py`** (driver: `run_syst_cosmics_chunked.sh`): chunked off-beam CV vs intime unisim, writes `Cosmics/` under the syst disk root.
- **`get_systematics_cosmics.py`**: cosmics-only NPZ in one pass. Requires **`--syst-disk-root`** or **`NUMUCC_SYST_DISK_ROOT`**.

---

## GENIE multisim

- **`get_systematics_genie.py`**: GENIE knob uncertainties. Requires **`evt`** and **`mcnu`** dataframes and separate **rate** vs **cross-section** covariance logic (see script docstring). Supports monolithic use or **`chunk-map`** / **`chunk-merge`** for HDF-split processing.
- **`syst_genie_parallel.py`** / **`syst_genie_aggregate.py`** (drivers: `run_syst_genie_chunked.sh`, `run_genie_mp.sh`): parallel chunked orchestration.
- **`merge_integrated_genie_into_final.py`**: folds integrated GENIE covariances into the final covariance tree.

---

## Detector variations

- **`syst_detvar_chunk.py`** → **`syst_detvar_aggregate.py`** (driver: `run_syst_detvar_chunked.sh`): WireMod / calorimetry / E-field unisim covariances from DetVar CAF `.df` files.
- **`wiremod_match_common_events.py`** (driver: `run_match_detvars.sh`): match common events between CV and WireMod samples.
- **`sce_match_common_events.py`** (driver: `run_match_sce.sh`): same for 0x/2x SCE samples.

---

## Joint (cross-variable) covariances for the conditional constraint

- **`syst_cc_joint_multisim_{chunk,parallel,aggregate}.py`** (driver: `run_syst_cc_joint_multisim_chunked.sh`) and
  **`syst_cc_joint_genie_{chunk,parallel,aggregate}.py`** (driver: `run_syst_cc_joint_genie_chunked.sh`, combined: `run_cc_systs.sh`):
  build joint covariances across variable pairs, consumed by `../cc_joint_cov.py`.
- **`conditional_constraint_validation.py`**: conditional Gaussian constraint (muon → proton) validation plots and JSON diagnostics.

---

## Event selection

- **`event_selection_batch_survey.py`** → **`event_selection_batch_map.py`** → **`event_selection_aggregate.py`**
  (driver: **`run_event_selection_batched.sh`**): groups input `.df` files into ≤1 GiB jobs,
  runs the selection pipeline per job, aggregates histograms, and renders final plots.
- **`selected_events.py`** / **`selected_events_cumulative.py`** (drivers: `run_event_rate_comp*.sh`):
  final-selection data/MC rate comparisons, per exposure batch or cumulative.
- **`unfolding_data.py`**: scripted Wiener-SVD unfolding of beam data (twin of `notebooks/unfolding-data.ipynb`).

---

## Utilities / tests

- **`merge_grid_job_dfs.py`** (driver: `run_merge_job_outputs.sh`): merge per-grid-job `.df` outputs.
- **`test_wgt_df_configs.py`** (driver: repo-root `test_wgt_jobs.sh`): smoke test of the weight df configs on a single CAF.
- **`run_workflow_test.py`**: integrated multisim + DetVar + selection smoke test with capped file counts.

---

## Quick reference: which script when?

| Goal | Script(s) |
|------|-----------|
| Chunked neutrino multisim (MCstat/Flux/G4) → NPZs + plots | `syst_multisim_chunk.py` → `syst_multisim_aggregate.py` (or `run_syst_multisim_chunked.sh`) |
| Cosmics only | `run_syst_cosmics_chunked.sh` or `get_systematics_cosmics.py` |
| GENIE rate + xsec (with `mcnu`) | `get_systematics_genie.py` (or `run_syst_genie_chunked.sh`) |
| Detector variations | `run_syst_detvar_chunked.sh` |
| Joint covariances for the constraint | `run_cc_systs.sh` |
| Event selection map/reduce + plots | `run_event_selection_batched.sh` |
| Unfolded cross section from data | `unfolding_data.py` |

---

## Related code

- **`utils.get_syst_unc`** — reads the **full** `syst_disk_layout` tree (`NUMUCC_SYST_DISK_ROOT`) and aborts if any expected file is missing.
- **`event_selection_aggregate.py`** — chunked selection reduce; passes `--syst-disk-root` into `get_syst_unc` when chunk pickles omit universe histograms.
