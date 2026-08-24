# νμ CC 1p0π cross-section analysis (`numucc_1p0pi`)

Code release for the SBND νμ CC 1p0π differential cross-section analysis.
Workflow: CAF ROOT files → pandas dataframes (`.df`, via `run_df_maker.py` +
`configs/numucc_1p0pi/`) → event selection → systematics → unfolding.
See the repo-root `README.md` for the phase-by-phase run instructions and
`scripts/README.md` for the systematics/selection drivers.

## Package modules (`analysis_village/numucc_1p0pi/`)

- `categories.py` — truth topology, fiducial-volume, and signal-definition masks.
- `constants.py` — analysis constants (POT, number of targets, bin edges, etc.).
- `variable_configs.py` — `VariableConfig` registry: binning, labels, and save names for all histogrammed variables.
- `final_selected_evt_vars.py` — registries of final-selection and intermediate-cut variables.
- `evt_derived_kinematics.py` — derived μ/p kinematics columns (momenta, angles, TKI inputs).
- `dataset_locations.py` — central input globs, work roots, and syst-disk roots; override with `NUMUCC_SPRING_GEN1_ROOT`.
- `files_config.py` — monolithic `.df` loaders (`get_ana_dfs`) used by scripts and most notebooks.
- `files_config_new.py` — split-file loaders used by a subset of notebooks (with `pyanalib.split_df_helpers_new`).
- `exposure_access.py` — staged data-access policy (`DataAccessStage`) and exposure-batch definitions.
- `utils.py` — plotting/overlay helpers, event-rate builders, `get_syst_unc`, χ² wiring; the hub imported by nearly everything.
- `selection_framework.py` — chunked selection engine (`ChunkRunner`, histogram accumulation/merging).
- `event_selection_pipeline_def.py` — the single definition of the cut/plot pipeline (edit cuts here).
- `event_selection_batched.py` — batched (≤1 GiB per job) selection orchestrator.
- `event_selection_batch_core.py` — per-batch selection runner used by the map jobs.
- `syst_disk_layout.py` — on-disk layout of the systematics NPZ/pickle tree (`MCstat/`, `Flux/`, `G4/`, `GENIE/`, `Cosmics/`, `Detector/`).
- `syst_disk_cc_layout.py` — same for the joint (cross-variable) covariance tree.
- `syst_multisim_common.py` — shared MCstat/Flux/G4 multisim helpers.
- `syst_cosmics_common.py` — cosmics variable registry and NPZ helpers.
- `syst_cc_joint_multisim_common.py` — joint-pair layout and filename helpers.
- `syst_pipeline_walker.py` — walks the selection pipeline stage-by-stage on sel_all dataframes (for cut-stage systematics).
- `syst_category_summary.py` — pack/load per-category systematic summary NPZ.
- `cc_joint_cov.py` — builds the joint covariance for the conditional (muon → proton) constraint.
- `genie_flat_helpers.py` — flat-GENIE / generator-comparison cross-section helpers.
- `fake_data_test_configs.py` — reweighting configurations for fake-data unfolding tests.

## Dataframe makers (`makedf/` in this directory)

- `makedf/makedf.py` — CAF→dataframe makers at each selection stage (`all` → `2prong` → `mup`), plus weight, calorimetry-variation, and E-field-variation variants.
- `makedf/selections.py` — cut primitives and μ/p candidate identification shared by the makers and the selection pipeline.

## Dataframe configs (`configs/numucc_1p0pi/`)

- `sel_all-mc.py` / `sel_all-data.py` — loose preselection dataframes (MC / data), inputs to cut-stage studies.
- `sel_2prong-mc.py` / `sel_2prong-data.py` / `sel_2prong-wgts-mc.py` — 2-prong stage dataframes (MC, data, MC with systematic weights).
- `sel_mup.py` / `sel_mup-data.py` — final μ+p selection dataframes (MC / data).
- `sel_mup-fluxwgts-knobgroups.py` — BNB flux multisim weights.
- `sel_mup-g4wgts.py` — Geant4 reinteraction multisim weights.
- `sel_mup-geniewgts-knobgroups.py` / `sel_mup-genieslimwgts.py` — GENIE knob-group / slimmed GENIE weights (`GENIE_KNOB_GROUP` env selects the group).
- `sel_mup-mcstatwgts.py` — Poisson MC-statistics universe weights.
- `sel_mup-updatecalo.py` / `sel_2prong-updatecalo.py` — recomputed calorimetry (±1σ parameter variations).
- `sel_2prong-updateefield.py` — recomputed PID with the double-anode E-field map.
- `add_ar23p.py` — preprocess hook adding AR23_20i_00_000 reweight knobs to CAFs.

## Scripts (`scripts/`) — see `scripts/README.md` for details

- `run_event_selection_batched.sh` + `event_selection_batch_{survey,map}.py` + `event_selection_aggregate.py` — map/reduce event selection and plotting.
- `run_syst_multisim_chunked.sh` + `syst_multisim_{chunk,parallel,aggregate}.py` — MCstat/Flux/G4 multisim covariances.
- `run_syst_cosmics_chunked.sh` + `syst_cosmics_{chunk,aggregate}.py`, `get_systematics_cosmics.py` — cosmic-background unisim covariances.
- `run_syst_genie_chunked.sh`, `run_genie_mp.sh` + `get_systematics_genie.py`, `syst_genie_{parallel,aggregate}.py`, `merge_integrated_genie_into_final.py` — GENIE covariances (rate and cross-section).
- `run_syst_detvar_chunked.sh` + `syst_detvar_{chunk,aggregate}.py` — detector-variation covariances.
- `run_match_detvars.sh` / `run_match_sce.sh` + `wiremod_match_common_events.py` / `sce_match_common_events.py` — event matching between CV and detector-variation samples.
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

- `event_selection.ipynb` — production event selection (batched workflow).

Data/MC comparison and validation:

- `data_mc_comparison.ipynb` — final-selection data vs MC overlays with uncertainties.
- `data_mc_comparison-chi2_summary.ipynb` — χ²/ndof summary tables across variables.
- `data_mc_comparison_gibuu.ipynb` — same overlays with GiBUU as the MC model.
- `homongeneity.ipynb` — spatial (octant) and temporal homogeneity of the selection.
- `data_driven_validation.ipynb` — conditional constraint (muon → proton) validation.

Systematics:

- `systematics.ipynb` — MCstat / Flux / G4 multisim covariance production.
- `systematics-genie.ipynb` — GENIE covariance production.
- `systematics-cosmic.ipynb` — cosmic-background systematics.
- `systematics-detector.ipynb` — detector-variation systematics orchestration.
- `systematics-summary.ipynb` — per-category uncertainty breakdown export.
- `systematics-flux_asymmetry.ipynb` — flux-universe asymmetry inspection.
- `multisigma_to_multisim.ipynb` — multisigma → multisim conversion for GENIE knobs.
- `wiremod.ipynb` / `sce.ipynb` / `dent.ipynb` — WireMod, SCE, and DENT matched-event detector studies.
- `detector_Efield_doubleanode.ipynb` — validation of the in-repo double-anode E-field map.
- `total_uncertainty_del_Tp.ipynb` — total uncertainty example for δp_T.

Unfolding and generators:

- `unfolding-data.ipynb` — Wiener-SVD unfolding of beam data.
- `unfolding-fake_data_tests.ipynb` — Asimov closure and fake-data tests.
- `unfolding.ipynb` — unfolded data vs generator truth comparisons.
- `generator_comparison.ipynb` — cross-section model comparisons on the measurement.

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
- `submit_mc_jobs_detvar.sh` / `run_efieldvar.sh` — detector-variation dataframe jobs.
- `makedf-lowE.sh` — low-energy dirt sample dataframes.
- `merge_selection_chunks.sh` — template for aggregating selection chunk outputs.
- `test_wgt_jobs.sh` — weight-config smoke test on a single CAF.
