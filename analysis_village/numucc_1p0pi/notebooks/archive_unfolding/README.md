# Archived unfolding notebooks

Prefer the live notebooks in `notebooks/`:

| Archived | Use instead |
|---|---|
| `unfolding.ipynb` (generator overlays) | `../generator_comparison.ipynb` + unfold products from `../unfolding.ipynb` |
| `unfolding-data.ipynb` | `../unfolding.ipynb` (Gen1 cache Wiener-SVD) |
| `unfolding-fake_data_tests.ipynb` | keep for archaeology / closure tests only |
| `unfolding-sanity-overlay.ipynb` | `../unfolding.ipynb` + `../unfolding-genie-comparison.ipynb` |

**Canonical workflow**

1. `../unfolding.ipynb` — rebuild Gen1 Wiener-SVD unfold; write results pickle/npz  
2. `../unfolding-genie-comparison.ipynb` — same ingredients, old vs new GENIE `total_xsec`  
3. `../generator_comparison.ipynb` — unfolded data vs generator predictions  
