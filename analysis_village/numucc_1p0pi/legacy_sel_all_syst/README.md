# Legacy Product A1 — sel_all chunked walk

Use this directory when you already have **weight-bearing `sel_all` `.df` files**
and want selection-stage systematics **without** regenerating CAF histcounts
(Product A2).

## What this is

| | Product **A1** (here) | Product **A2** (default going forward) |
|---|---|---|
| Input | Existing `sel_all` weight `.df` (`evt`+`trk`+`hdr`[+`mcnu`]) | CAF → `syst_hists__<var>_*` via `configs/numucc_1p0pi/syst_histcounts.py` |
| Engine | Re-walk `build_pipeline` per file → pickle chunks | Fill long histcounts on the grid (per-variable HDF keys) |
| Output | `nu__*.pkl` / `genie__*.pkl` → aggregate → syst disk | Stream-sum `syst_hists` → covs |
| GENIE xsec on cut-stage vars | No (rate only at cut stages) | No for selection-stage; xsec is Product B |

Weight DF production (`makedf/*syst*`) is unchanged and shared.

## How to run

From the **repo root** (scripts call into `analysis_village/numucc_1p0pi/scripts/`):

```bash
# Flux / G4 (and optional MCstat) on sel_all
./analysis_village/numucc_1p0pi/legacy_sel_all_syst/scripts/run_syst_multisim_sel_all.sh

# GENIE groups on sel_all (cut-stage = rate only; final vars still get xsec)
./analysis_village/numucc_1p0pi/legacy_sel_all_syst/scripts/run_syst_genie_sel_all.sh

# Cosmics on sel_all
./analysis_village/numucc_1p0pi/legacy_sel_all_syst/scripts/run_syst_cosmics_sel_all.sh
```

Or set `MC_DF_STAGE=sel_all` yourself when calling the main chunked drivers under
`../scripts/`.

## Notebooks

- `notebooks/systematics-selection-genie-Ar23p.ipynb` — interactive GENIE Ar23p on
  sel_all / selected tables (rate + optional xsec for measurement vars).

For **new** selection-stage campaigns, prefer Product A2:

- Config: `configs/numucc_1p0pi/syst_histcounts.py`
- Sum/cov: `scripts/syst_histcounts_stream_sum.py`, `notebooks/systematics-histcounts.ipynb`
- Selective load: `load_syst_hists_from_df_file(path, vars=["nu_score"])`

See the package README → **Systematics products**.
