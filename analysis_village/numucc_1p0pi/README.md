# νμ CC 1p0π cross-section analysis (`numucc_1p0pi`)

Code release for the SBND νμ CC 1p0π differential cross-section analysis.
Workflow: CAF ROOT files → pandas dataframes (`.df`, via `run_df_maker.py` +
`configs/numucc_1p0pi/`) → event selection → systematics → unfolding.
See the repo-root `README.md` for the phase-by-phase run instructions and
`scripts/README.md` for the systematics/selection drivers.

## Event selection (read this first)

Selection is defined in **two** places only. Everything else (CAF makers,
batched notebook, syst walkers) should call into them — do not copy thresholds
or cut order into notebooks.

| Edit this… | When you want to… |
|---|---|
| `makedf/selections.py` | Change a **threshold** or cut formula (`NU_SCORE_TH`, PID, kinematics, FV) |
| `event_selection_pipeline_def.py` (`build_pipeline`) | Change **cut order**, stage names, or which plots/syst vars appear at each stage |

CAF products (written by `makedf/makedf.py` → `make_pandora_evtdf`):

| Product | `sel_level` | Use for |
|---|---|---|
| `sel_all` | `all` | Cut-stage studies / systematics (loose slices + tracks) |
| `sel_2prong` | `2prong` | Mid-selection checks |
| `sel_mup` | `mup` | Final rates, xsec, unfolding |

Stage names in the pipeline use hyphens (`2prong-mup`); CAF `sel_level` uses
underscores / short names (`mup`, `muX`). Mapping:
`CAF_SEL_LEVEL_TO_STAGE` in `event_selection_pipeline_def.py`.

The frozen predecessor under `analysis_village/numucc1p0pi-old/` uses different
thresholds and FV — do not mix it with Gen-1 / per-TPC DFs.

Notebooks:

| Notebook | Use for |
|---|---|
| `notebooks/event_selection.ipynb` | Develop / tune cuts on a few files (tweak thresholds in-notebook) |
| `notebooks/event_selection_batched.ipynb` | Live batched overlays on full samples + summary / efficiency |

## Systematics products

Weight DF production (`makedf/{getsyst,geniesyst,bnbsyst,g4syst,mcstat}.py`) is
centralized — **do not fork it**. Everything below is about *using* those weights.

| Product | Variables | Inputs | GENIE xsec? | Default path |
|---|---|---|---|---|
| **A — Selection-stage** | Cut-stage vars (nu_score, vtxdist, χ², …) + optional topology breakdown | Large tables | **No** (rate only) | **A2** histcounts |
| **B — Measurement** | Integrated, μ/p kinematics, TKIs, … | `sel_mup` (+ `mcnu` for GENIE xsec) | **Yes** | Chunked / GENIE scripts with `MC_DF_STAGE=final` |

### Product B — Measurement (`sel_mup` + `mcnu` for GENIE xsec)

Entry points (all force the final / measurement stage):

| Script | What it runs |
|---|---|
| `scripts/run_syst_measurement_multisim.sh` | Flux / G4 / MCstat (`MC_DF_STAGE=final`) |
| `scripts/run_syst_measurement_genie.sh` | GENIE rate **+ xsec** (`MC_DF_STAGE=final`) |
| `scripts/run_syst_measurement_cosmics.sh` | Cosmics (`INPUT_STAGE=final`) |

For cosmics, set `COSMICS_SELECTED_MC_DF=/path/to/sel_mup.df` (HDF `evt`) so aggregate
attaches contamination-scaled `SelectedRate` (needed by `systematics-summary.ipynb`).

Interactive GENIE inspection: `notebooks/systematics-genie-inspect.ipynb`
(retired group/mixed-summary notebooks live under `notebooks/archive_syst/`).

Interactive Product **B** multisim (notebook production):
`systematics-mcstat.ipynb`, `systematics-flux.ipynb`, `systematics-g4.ipynb`
(helpers in `syst_multisim_inspect.py`).

### Shared covariance math

- `syst_genie_cov.py` — GENIE univ alias (`/cv`) + xsec accumulate/finalize
- `syst_multisim_common.combine_indep_knob_cov_packs` — independent-knob total
  (**sum of `cov_frac`**); histcounts / notebooks call this (no local twins)

### Product A2 (default) — CAF / DF histcounts

- Fill: `configs/numucc_1p0pi/syst_histcounts.py` → `make_syst_histcounts*`
- HDF keys are **per variable**: `syst_hists__<var_save_name>_<split>`
  (legacy monolith `syst_hists_<split>` still loads). Selective read::

      load_syst_hists_from_df_file(path, vars=["nu_score"])

- Sum / cov: `scripts/syst_histcounts_stream_sum.py`, `notebooks/systematics-histcounts.ipynb`

### Product A1 (legacy) — existing sel_all weight DFs

If you already built `sel_all` weight tables and need to walk the selection
pipeline offline, use:

`analysis_village/numucc_1p0pi/legacy_sel_all_syst/`

(wrappers set `MC_DF_STAGE=sel_all`; see that README).

## Package modules (`analysis_village/numucc_1p0pi/`)

- `categories.py` — truth topology, fiducial-volume, and signal-definition masks.
- `constants.py` — analysis constants (POT, number of targets, bin edges, etc.).
- `variable_configs.py` — `VariableConfig` registry: binning, labels, and save names for all histogrammed variables.
- `final_selected_evt_vars.py` — registries of final-selection and intermediate-cut variables.
- `evt_derived_kinematics.py` — derived μ/p kinematics columns (momenta, angles, TKI inputs).
- `dataset_locations.py` — central input globs, work roots, and syst-disk roots (`default_syst_disk_root` → PRL `productB_sel_mup`; `prl_syst_disk_root("A"|"B")`); override with `NUMUCC_SPRING_GEN1_ROOT` / `NUMUCC_SYST_DISK_ROOT`.
- `files_config.py` — monolithic `.df` loaders (`get_ana_dfs`) used by scripts and most notebooks.
- `files_config_new.py` — split-file loaders used by a subset of notebooks (with `pyanalib.split_df_helpers_new`).
- `exposure_access.py` — staged data-access policy (`DataAccessStage`) and exposure-batch definitions.
- `utils.py` — plotting/overlay helpers, event-rate builders, `get_syst_unc`, χ² wiring; the hub imported by nearly everything.
- `selection_framework.py` — chunked selection engine (`ChunkRunner`, histogram accumulation/merging).
- `event_selection_pipeline_def.py` — **cut order + stage plots** (`build_pipeline`); CAF / notebook / syst all consume this.
- `event_selection_batched.py` — batched (≤1 GiB per job) selection orchestrator.
- `event_selection_batch_core.py` — per-batch selection runner used by the map jobs.
- `legacy_samples.py` — in-memory `SampleBundle` for small tests (`run_pipeline()` wraps `build_pipeline`).
- `syst_disk_layout.py` — on-disk layout of the systematics NPZ/pickle tree (`MCstat/`, `Flux/`, `G4/`, `GENIE/`, `Cosmics/`, `Detector/`); canonical consumer roots under `PRL/systematics/product{A,B}_*`.
- `syst_disk_cc_layout.py` — same for the joint (cross-variable) covariance tree.
- `syst_multisim_common.py` — shared MCstat/Flux/G4 multisim helpers; **canonical** `combine_indep_knob_cov_packs`.
- `syst_genie_cov.py` — shared GENIE univ alias + xsec accumulate/finalize (used by `get_systematics_genie` and `syst_histcounts`).
- `syst_genie_inspect.py` — load/plot helpers for GENIE `cov_mat_dict` inspection notebook.
- `syst_cosmics_common.py` — cosmics variable registry, flat unc, **SelectedRate** contamination helpers.
- `syst_detvar_common.py` — WireMod/DENT matching wrappers, batched hist fill, envelope/unisim packs, Detector combine/plots.
- `syst_summary_inspect.py` — load multi-root syst disks + total-by-source plots for `systematics-summary.ipynb`.
- `syst_cc_joint_multisim_common.py` — joint-pair layout and filename helpers.
- `syst_pipeline_walker.py` — walks `build_pipeline` on sel_all dfs (Product A1 / shared fill).
- `syst_histcounts.py` — Product **A2** histcounts pack/unpack; per-variable HDF keys.
- `legacy_sel_all_syst/` — Product **A1** runners + notebooks for existing sel_all weight DFs.
- `syst_category_summary.py` — pack/load per-category systematic summary NPZ.
- `cc_joint_cov.py` — builds the joint covariance for the conditional (muon → proton) constraint.
- `genie_flat_helpers.py` — flat-GENIE / generator-comparison cross-section helpers.
- `fake_data_test_configs.py` — reweighting configurations for fake-data unfolding tests.

## Dataframe makers (`makedf/` in this directory)

- `makedf/makedf.py` — CAF→dataframe makers; selection via `apply_selection_pipeline` (`sel_all` → `sel_mup`).
- `makedf/selections.py` — **thresholds + cut primitives** shared by makers and the selection pipeline.

## Dataframe configs (`configs/numucc_1p0pi/`)

- `sel_all-mc.py` / `sel_all-data.py` — loose preselection dataframes (MC / data), inputs to cut-stage studies.
- `sel_2prong-mc.py` / `sel_2prong-data.py` / `sel_2prong-wgts-mc.py` — 2-prong stage dataframes (MC, data, MC with systematic weights).
- `sel_mup.py` / `sel_mup-data.py` — final μ+p selection dataframes (MC / data).
- `sel_mup-fluxwgts-knobgroups.py` — BNB flux multisim weights.
- `sel_mup-g4wgts.py` — Geant4 reinteraction multisim weights.
- `sel_mup-geniewgts-knobgroups.py` / `sel_mup-genieslimwgts.py` — GENIE knob-group / slimmed GENIE weights (`GENIE_KNOB_GROUP` env selects the group).
- `sel_mup-mcstatwgts.py` — Poisson MC-statistics universe weights.
- `sel_all-updatecalo.py` — detector variations at `sel_all`: CV + ±1σ calo universes (`evt_*` / `trk_*`) plus E-field redo (`evt_efield` / `trk_efield`).
- `add_ar23p.py` — preprocess hook adding AR23_20i_00_000 reweight knobs to CAFs.

## Scripts (`scripts/`) — see `scripts/README.md` for details

- `run_event_selection_batched.sh` + `event_selection_batch_{survey,map}.py` + `event_selection_aggregate.py` — map/reduce event selection and plotting.
- `run_syst_measurement_{multisim,genie,cosmics}.sh` — Product **B** entry points (`MC_DF_STAGE`/`INPUT_STAGE=final`).
- `run_syst_multisim_chunked.sh` + `syst_multisim_{chunk,parallel,aggregate}.py` — MCstat/Flux/G4 multisim covariances.
- `run_syst_cosmics_chunked.sh` + `syst_cosmics_{chunk,aggregate}.py`, `get_systematics_cosmics.py` — cosmic-background unisim covariances.
- `run_syst_genie_chunked.sh`, `run_genie_mp.sh` + `get_systematics_genie.py`, `syst_genie_{parallel,aggregate}.py`, `merge_integrated_genie_into_final.py` — GENIE covariances (rate and cross-section).
- `run_syst_detvar_chunked.sh` + `syst_detvar_{chunk,aggregate}.py` — detector-variation covariances (legacy chunk path).
- `run_match_detvars.sh` / `run_match_dent.sh` / `run_match_sce.sh` + `wiremod_match_common_events.py` / `dent_match_common_events.py` / `sce_match_common_events.py` — backends for event matching; prefer `notebooks/systematics-detector-match.ipynb`.
- `run_cc_systs.sh`, `run_syst_cc_joint_{multisim,genie}_chunked.sh` + `syst_cc_joint_*` — joint covariances for the conditional constraint.
- `conditional_constraint_validation.py` — constraint validation plots and diagnostics.
- `selected_events.py` / `selected_events_cumulative.py` + `run_event_rate_comp*.sh` — data/MC rate comparisons per exposure batch.
- `unfolding_data.py` — scripted Wiener-SVD unfolding of beam data.
- `merge_grid_job_dfs.py` + `run_merge_job_outputs.sh` — merge per-grid-job `.df` outputs.
- `test_wgt_df_configs.py`, `run_workflow_test.py` — smoke tests for weight configs and the integrated workflow.

## Notebooks (`notebooks/`)

Exposure and flux:

- `exposure_accounting.ipynb` — POT / trigger accounting and absolute normalization.
- `beam_quality.ipynb` — spill-level beam-quality metrics and cuts.
- `flux_closure.ipynb` — flux closure: dk2nu / gsimple vs reported SBND flux.
- `genie_vs_production_xsec.ipynb` — GENIE flat-tree vs production rate closure.

Event selection and PID:

- `event_selection.ipynb` — develop / tune selection in-memory (few files, tweakable thresholds).
- `event_selection_batched.ipynb` — live batched overlays on full (or partial) samples + summary/efficiency.

Data/MC comparison and validation:

- `data_mc_comparison.ipynb` — final-selection data vs MC overlays with uncertainties.
- `data_mc_comparison-chi2_summary.ipynb` — χ²/ndof summary tables across variables.
- `data_mc_comparison_gibuu.ipynb` — same overlays with GiBUU as the MC model.
- `homongeneity.ipynb` — spatial (octant) and temporal homogeneity of the selection.
- `data_driven_validation.ipynb` — conditional constraint (muon → proton) validation.

Systematics:

- `systematics-histcounts.ipynb` — Product **A2** (default): load/sum per-variable `syst_hists`, build covs.
- `systematics-genie-inspect.ipynb` — **inspect** GENIE syst-disk outputs (per-knob / per-mode / top-10 plots; rate+xsec).
- `systematics-mcstat.ipynb` / `systematics-flux.ipynb` / `systematics-g4.ipynb` — Product **B** multisim covariances (helpers in `syst_multisim_inspect.py`). Flux notebook includes integrated asymmetry (former `systematics-flux_asymmetry`).
- `systematics-multisim-live.ipynb` — quick live walk of final-selected weight dfs (debug / spot-check).
- `systematics-cosmic.ipynb` — cosmic-background systematics (`SelectedRate` via `syst_cosmics_common`).
- `systematics-detector-match.ipynb` — match detector-variation events at **sel_all** by `(E, run, subrun, evt)` (WireMod + DENT).
- `wiremod.ipynb` / `dent.ipynb` — Product **A** (cut-stage) + **B** (measurement) from matched sel_all pipeline walks; WireMod calo envelope / DENT unisim (`syst_detvar_common.py`).
- `systematics-detector.ipynb` — combine WireMod YZ/XTXW + DENT → `Detector/detector_syst_dict.npz` + overlay plots.
- `systematics-summary.ipynb` — combine source disks → CategorySummary + total-by-source plots; **Detector** is one source from `systematics-detector.ipynb`.
- `prl-genie-syst-summary.ipynb` / `total_uncertainty_del_Tp.ipynb` — specialized GENIE / δp_T summaries.
- `multisigma_to_multisim.ipynb` — multisigma → multisim conversion for GENIE knobs.
- `notebooks/archive_syst/` — retired near-duplicates (`systematics-genie.FULL`, group/mixed-summary, …).
- `legacy_sel_all_syst/` — Product **A1** (existing sel_all weight DFs); see that README.
- `sce.ipynb` — SCE matched-event study (legacy; Detector total now uses WireMod + DENT).
- `detector_Efield_doubleanode.ipynb` — validation of the in-repo double-anode E-field map.

Unfolding and generators:

- `unfolding-prepare.ipynb` — **PRL Product B step 1:** load Sep-1 `sel_mup` DFs (beam-quality data + MC `evt`/`mcnu`), recompute data–MC overlays and assert vs `PRL/data_mc_overlays/productB_sel_mup/counts_report.npz`, build efficiency/response matrices → `PRL/response_matrices/`.
- `unfolding.ipynb` — **PRL Product B step 2:** load response pack + Product B CategorySummary `total_xsec`, MC closure test, data Wiener-SVD (`C_type=2`), save flux + unfolded products under `PRL/unfolded/`.
- `unfolding-legacy-gen1.ipynb` — May Gen1 recovered-cov rebuild (`CovRotation` recovery; χ²≈34.5/12 for `tki-del_Tp`). Do not use for the Product B data release.
- `unfolding-genie-comparison.ipynb` — currently wired to Gen1 ingredients; Old vs New GENIE `total_xsec` extracted xsecs.
- `generator_comparison.ipynb` — unfolded data vs generator predictions.
- `notebooks/archive_unfolding/` — retired unfold notebooks (`unfolding-data`, fake-data tests, …).

Style:

- `presentation.mplstyle` — shared matplotlib style (loaded by `utils.py`).

## Flux estimation (`analysis_village/flux/`)

- `raytrace_volume_defs.py` — detector / fiducial volume definitions for flux ray tracing (imported by the unfolding code).
- `build_voxel_flux.py` — voxelized flux map builder from gsimple files.
- `gsimple_batch.py` / `gsimple_raytrace_batch.py` + `run_gsimple.sh` / `run_raytrace.sh` — batch gsimple readers and FV ray-trace flux integration.
- `dk2nu.ipynb` / `gsimple.ipynb` / `gsimple_raytrace.ipynb` / `compare_gsimple_raytrace_batch.ipynb` — dk2nu and gsimple flux inspection, ray-traced flux estimation, and cross-checks.
- `raytrace_figures/` — example ray-trace illustration figures.

## Unfolding library (`analysis_village/unfolding/`)

- `wienersvd.py` — updated Wiener-SVD unfolding implementation used by `utils.py` and the unfolding notebooks/scripts.

## Common-code updates included in this release

- `run_df_maker.py` — `-ncpu` option, per-file preprocess/args hooks, POT/genevt histogram tables for empty recTrees, `GENIE_KNOB_GROUP` forwarding to grid workers, grid-submission fixes.
- `setup.sh` — spack/venv environment setup (py3.10, SAMWeb, hdf5, xrootd).
- `bin/grid_executable.sh`, `bin/init_grid.sh` — grid worker bootstrap for `.df` production jobs.
- `makedf/makedf.py` — calorimetry/E-field variation support, weight wiring (GENIE/flux/G4/MCstat), POT/genevt histograms.
- `makedf/chi2pid.py`, `makedf/calo.py` — χ² PID recalculation with `CALO_VARIATIONS` parameter grid; Wion fix.
- `makedf/branches.py` — extended branch lists (PID, shower, true-hit info).
- `makedf/bnbsyst.py`, `makedf/geniesyst.py`, `makedf/g4syst.py`, `makedf/getsyst.py` — slimmed/grouped multisim weight readers (`GENIE_KNOB_GROUPS`, flux knob groups, neutron knob).
- `makedf/mcstat.py` — Poisson MC-statistics universe weights.
- `makedf/getenv.py` — CAF `env` tree reader (file provenance).
- `makedf/constants.py`, `makedf/util.py` — detector constants, SBND fiducial-volume variants, `avg_chi2`, track-matching helpers.
- `pyanalib/covariance.py` — covariance / fractional-covariance / correlation utilities.
- `pyanalib/stat_helpers.py` — χ² with combined covariance, asymmetric (gamma-based) data errors.
- `pyanalib/split_df_helpers.py`, `pyanalib/split_df_helpers_new.py` — loaders for split/chunked `.df` files.
- `pyanalib/variable_calculator.py` — CC 1p0π TKI variable calculators.
- `pyanalib/ntuple_glob.py`, `pyanalib/pandas_helpers.py` — per-file args/preprocess support, faster branch loading.
- `data/efield/sbnd_sce_doubleanode_2d_v10c.root` — double-anode SCE / E-field map shipped for `updateefield` dataframe jobs.
- `preprocess/` — CAF preprocessing hook + AR23_20i_00_000 reweight script and fcl.
- `flux_closure/` — GENIE GHEP → dataframe utility for flux/xsec closure (requires the external `pyGENIE` package).

## Root-level driver scripts (dataframe production bookkeeping)

- `submit_mc_jobs.sh` — BNB+cosmics MC (and GiBUU) dataframe grid jobs.
- `submit_data_jobs.sh` — beam-on data dataframes.
- `submit_offbeam_jobs.sh` / `submit_intime_jobs.sh` / `submit_dirt_jobs.sh` — off-beam data, in-time cosmics, and dirt MC dataframes.
- `submit_mc_jobs_{flux,g4,GENIE,GENIEslim,mcstat}.sh` — systematic-weight dataframe jobs.
- `submit_mc_jobs_detvar.sh` — detector-variation dataframe jobs (`sel_all-updatecalo` includes efield).
- `makedf-lowE.sh` — low-energy dirt sample dataframes.
- `merge_selection_chunks.sh` — template for aggregating selection chunk outputs.
- `test_wgt_jobs.sh` — weight-config smoke test on a single CAF.
