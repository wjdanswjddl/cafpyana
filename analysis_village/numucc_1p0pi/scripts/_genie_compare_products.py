"""Compare total GENIE fractional uncertainty between two cov_mat_dict.pkl files."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/exp/sbnd/app/users/munjung/xsec/freeze/cafpyana")

from analysis_village.numucc_1p0pi.scripts._genie_pkl_compat import (  # noqa: E402
    load_cov_mat_dict_compat,
)


def _tot(row, key) -> np.ndarray | None:
    m = row.get(key)
    return None if m is None else np.sqrt(np.clip(np.diag(np.asarray(m)), 0, None)) * 100.0


def main(a: str, b: str) -> None:
    ca, cb = load_cov_mat_dict_compat(a), load_cov_mat_dict_compat(b)
    print(f"A = {a}")
    print(f"B = {b}")
    print(f"{'variable':18s} {'kind':5s} {'maxA%':>8s} {'maxB%':>8s} {'max|A-B|%':>10s}")
    for var in sorted(set(ca) & set(cb)):
        for kind, key in (("rate", "genie_rate"), ("xsec", "genie")):
            fa, fb = _tot(ca[var], key), _tot(cb[var], key)
            if fa is None or fb is None or fa.shape != fb.shape:
                continue
            if not fa.any() and not fb.any():
                continue
            print(
                f"{var:18s} {kind:5s} {fa.max():8.3f} {fb.max():8.3f}"
                f" {np.abs(fa - fb).max():10.4f}"
            )
    # Knob sets should match too, otherwise the totals agreeing means little.
    ka = {k for k in ca["integrated"] if not k.startswith("genie")}
    kb = {k for k in cb["integrated"] if not k.startswith("genie")}
    print(f"\nknob keys: A={len(ka)} B={len(kb)} only-in-A={sorted(ka - kb)} only-in-B={sorted(kb - ka)}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
