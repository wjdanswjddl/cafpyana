"""Diff the Product A / Product B variable lists against a GENIE cov_mat_dict.pkl.

Product A = CUT_STAGE_RATE_ONLY_SLUGS (rate only by construction).
Product B = final_stage_var_configs() (rate + xsec).
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, "/exp/sbnd/app/users/munjung/xsec/freeze/cafpyana")

from analysis_village.numucc_1p0pi.scripts._genie_pkl_compat import (  # noqa: E402
    load_cov_mat_dict_compat as load_cov_mat_dict,
)
from analysis_village.numucc_1p0pi.syst_genie_inspect import (  # noqa: E402
    iter_knob_cov_fracs,
)
from analysis_village.numucc_1p0pi.syst_pipeline_walker import (  # noqa: E402
    CUT_STAGE_RATE_ONLY_SLUGS,
    final_stage_var_configs,
)


def _counts(row) -> tuple[int, int]:
    return (
        sum(1 for _ in iter_knob_cov_fracs(row, "rate")),
        sum(1 for _ in iter_knob_cov_fracs(row, "xsec")),
    )


def report(path: str) -> None:
    cov = load_cov_mat_dict(Path(path))
    prod_a = sorted(CUT_STAGE_RATE_ONLY_SLUGS)
    prod_b = sorted({vc.var_save_name for vc in final_stage_var_configs() if vc.var_save_name})

    print(f"### {path}")
    print(f"    vars in file={len(cov)}  productA expects={len(prod_a)}  productB expects={len(prod_b)}")

    for label, want, need_xsec in (("A", prod_a, False), ("B", prod_b, True)):
        ok, no_rate, no_xsec, missing = [], [], [], []
        for slug in want:
            row = cov.get(slug)
            if row is None:
                missing.append(slug)
                continue
            nr, nx = _counts(row)
            if nr == 0:
                no_rate.append(slug)
            if need_xsec and nx == 0:
                no_xsec.append(slug)
            if nr and (nx or not need_xsec):
                ok.append(slug)
        print(f"  product {label}: complete={len(ok)}/{len(want)}")
        for name, lst in (("absent", missing), ("no rate", no_rate), ("no xsec", no_xsec)):
            if lst:
                print(f"    {name} ({len(lst)}): {', '.join(lst)}")

    extra = sorted(set(cov) - set(prod_a) - set(prod_b))
    if extra:
        print(f"  in file but in neither list ({len(extra)}): {', '.join(extra)}")


if __name__ == "__main__":
    for p in sys.argv[1:]:
        report(p)
        print()
