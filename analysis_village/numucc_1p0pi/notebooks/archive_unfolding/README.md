# Archived unfolding notebooks

Prefer the live notebooks in `notebooks/`:

| Archived | Use instead |
|---|---|
| `unfolding.ipynb` (generator overlays) | `../generator_comparison.ipynb` + unfold products from `../unfolding.ipynb` |
| `unfolding-data.ipynb` | `../unfolding-prepare.ipynb` + `../unfolding.ipynb` (PRL Product B) |
| `unfolding-fake_data_tests.ipynb` | keep for archaeology; live closure is in `../unfolding.ipynb` |
| `unfolding-sanity-overlay.ipynb` | `../unfolding.ipynb` + `../unfolding-genie-comparison.ipynb` |

**Canonical PRL Product B workflow**

1. `../unfolding-prepare.ipynb` — Sep-1 DFs, data–MC assert vs `counts_report.npz`, response matrices → `PRL/response_matrices/`
2. `../unfolding.ipynb` — Product B `total_xsec` cov, MC closure, data Wiener-SVD → `PRL/unfolded/`
3. `../generator_comparison.ipynb` — unfolded data vs generator predictions

**Legacy Gen1 (May recovered cov)**

- `../unfolding-legacy-gen1.ipynb` — Gen1 cache + `CovRotation` recovery (χ²≈34.5/12 for `tki-del_Tp`)
- `../unfolding-genie-comparison.ipynb` — still points at Gen1 ingredients until updated
