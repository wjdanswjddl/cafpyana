import numpy as np

# Floor for nominal counts when building **fractional** covariances (avoids inf/NaN in empty bins).
_FRAC_CV_EPS = 1e-12


def cov_from_fraccov(cov_frac, cv_vals):
    """Absolute covariance from fractional matrix and per-bin CV (counts).

    ``0 * inf`` can appear when ``cov_frac`` came from upstream code that divided by zero
    nominal CV; those entries are cleared to zero so empty bins do not poison linear algebra.
    """
    cf = np.asarray(cov_frac, dtype=float)
    v = np.asarray(cv_vals, dtype=float)
    cov = cf * np.outer(v, v)
    return np.nan_to_num(cov, nan=0.0, posinf=0.0, neginf=0.0)


def fraccov_from_cov(cov, cv_vals):
    cov = np.asarray(cov, dtype=float)
    v = np.asarray(cv_vals, dtype=float).reshape(-1)
    if cov.shape[0] != cov.shape[1] or cov.shape[0] != v.shape[0]:
        raise ValueError("cov shape %s incompatible with cv_vals length %d" % (cov.shape, v.shape[0]))
    denom = np.maximum(np.abs(np.outer(v, v)), _FRAC_CV_EPS)
    return cov / denom

def corr_from_fraccov(cov_frac):
    corr = np.zeros_like(cov_frac)
    for i in range(cov_frac.shape[0]):
        for j in range(cov_frac.shape[1]):
            di = max(float(cov_frac[i, i]), 0.0)
            dj = max(float(cov_frac[j, j]), 0.0)
            denom = np.sqrt(di * dj)
            corr[i, j] = (cov_frac[i, j] / denom) if denom > _FRAC_CV_EPS else 0.0
    return corr

def get_covariance_matrix(univ_events, 
                          cv_events):
    n_univ, n_bins = univ_events.shape

    cov_frac = np.zeros((n_bins, n_bins))
    cov = np.zeros((n_bins, n_bins))

    # looping & calculating with the CV value for clarity, 
    # but techincally np.cov should also be fine under the assumption of gaussian universes that we're using
    for uidx in range(n_univ):
        for i in range(univ_events.shape[1]):
            for j in range(univ_events.shape[1]):
                nom_i = cv_events[i] 
                nom_j = cv_events[j] 

                univ_i = univ_events[uidx, i] 
                univ_j = univ_events[uidx, j] 

                cov_entry = (univ_i - nom_i) * (univ_j - nom_j)
                den_i = max(float(nom_i), _FRAC_CV_EPS)
                den_j = max(float(nom_j), _FRAC_CV_EPS)
                frac_cov_entry = ((univ_i - nom_i) / den_i) * ((univ_j - nom_j) / den_j)

                # TODO: uboone code has clipping that I'm not sure why.. investigate later
                # if cov_entry > 0:
                #     this_cov = max( cov_entry, eps * scale_factor)
                # else:
                #     this_cov = min( cov_entry, eps * scale_factor)

                # if frac_cov_entry > 0:
                #     this_frac_cov = max( frac_cov_entry, eps * scale_factor)
                # else:
                #     this_frac_cov = min( frac_cov_entry, eps * scale_factor)

                cov[i, j] += cov_entry
                cov_frac[i, j] += frac_cov_entry

    cov = cov / n_univ
    cov_frac = cov_frac / n_univ
    corr = np.zeros_like(cov)
    for i in range(len(cv_events)):
        for j in range(len(cv_events)):
            di = max(float(cov[i, i]), 0.0)
            dj = max(float(cov[j, j]), 0.0)
            denom = np.sqrt(di * dj)
            corr[i, j] = (cov[i, j] / denom) if denom > _FRAC_CV_EPS else 0.0

    return {"cov_frac": cov_frac, 
            "cov": cov,
            "corr": corr,
            }
