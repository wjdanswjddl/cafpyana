# ARCHIVED — `numucc1p0pi-old`

This tree is a **frozen predecessor** of the Gen-1 / per-TPC analysis in
`analysis_village/numucc_1p0pi/`.

Do **not** use it for new work or mix its cut values with Gen-1 DFs.

| Item | This tree (old) | Current (`numucc_1p0pi`) |
|---|---|---|
| `nu_score` | 0.5 | 0.45 |
| μ χ²_μ | 30 | 25 |
| MCS/range qual | 0.3 | 0.2 |
| FV / containment | nominal SBND (often ends-only) | Gen-1 + per-TPC start+end |
| Cut definitions | Hardcoded in maker + notebook | `makedf/selections.py` + `build_pipeline()` |

Canonical selection documentation:
`analysis_village/numucc_1p0pi/README.md` → **Event selection**.
