"""Report variable x (rate|xsec) coverage for one GENIE cov_mat_dict.pkl.

One file per process so the arrays are freed between products (shared machine).
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

TOTAL_RATE = "genie_rate"
TOTAL_XSEC = "genie"


def main(path: str) -> None:
    cov = load_cov_mat_dict(Path(path))
    print(f"path: {path}")
    print(f"variables: {len(cov)}")
    for var in sorted(cov):
        row = cov[var]
        n_rate = sum(1 for _ in iter_knob_cov_fracs(row, "rate"))
        n_xsec = sum(1 for _ in iter_knob_cov_fracs(row, "xsec"))
        has_tr = TOTAL_RATE in row
        has_tx = TOTAL_XSEC in row
        nbin = None
        for _k, m in iter_knob_cov_fracs(row, "rate" if n_rate else "xsec"):
            nbin = m.shape[0]
            break
        print(
            f"  {var:18s} rate_knobs={n_rate:3d} tot_rate={'Y' if has_tr else 'n'}"
            f"   xsec_knobs={n_xsec:3d} tot_xsec={'Y' if has_tx else 'n'}"
            f"   nbins={nbin}"
        )


if __name__ == "__main__":
    main(sys.argv[1])
