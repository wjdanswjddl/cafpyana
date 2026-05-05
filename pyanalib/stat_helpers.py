import numpy as np
from scipy.stats import gamma
from scipy.stats import chi2

alpha = 0.3173  # for 68.27% CL, alpha = 1 - 0.6827

def return_data_stat_err(data_array):
    L = np.where(data_array == 0, 0, gamma.ppf(alpha / 2, data_array))
    U = np.where(data_array == 0,
             gamma.isf(alpha, data_array + 1),
             gamma.isf(alpha / 2, data_array + 1))

    yerr_low = np.where(data_array == 0, 0.1, data_array - L)
    yerr_high = np.where(data_array == 0, 1.8, U - data_array)

    return yerr_low, yerr_high


def get_chi2(data, model, cov):
    chi2_value = (data - model) @ np.linalg.inv(cov) @ (data - model)
    ndof = len(data)
    p_value = 1 - chi2.cdf(chi2_value, df=ndof)
    return chi2_value, p_value


def get_chi2_shape(data, model, cov):
    """Shape chi2 using full covariance matrix: model normalized to data integral."""
    scale = np.sum(data) / np.sum(model)
    model_shape = model * scale
    chi2_value = (data - model_shape) @ np.linalg.inv(cov) @ (data - model_shape)
    ndof = len(data)
    p_value = 1 - chi2.cdf(chi2_value, df=ndof)
    return chi2_value, p_value


def get_chi2_avg(data, model, cov):
    """Chi2 averaged over bins, treating bins as independent (diagonal cov only)."""
    var = np.diag(cov)
    valid_idx = var > 0
    data = data[valid_idx]
    model = model[valid_idx]
    var = var[valid_idx]
    ndof = int(np.sum(valid_idx))
    chi2_per_bin = (data - model) ** 2 / var
    chi2_total = np.sum(chi2_per_bin)
    p_value = 1 - chi2.cdf(chi2_total, df=ndof)
    return chi2_total / ndof, p_value


def get_chi2_shape_avg(data, model, cov):
    """Shape chi2 averaged over bins: model normalized to data integral, bins treated independently."""
    scale = np.sum(data) / np.sum(model)
    model_shape = model * scale
    var = np.diag(cov)
    valid_idx = var > 0
    data = data[valid_idx]
    model_shape = model_shape[valid_idx]
    var = var[valid_idx]
    ndof = int(np.sum(valid_idx))
    chi2_per_bin = (data - model_shape) ** 2 / var
    chi2_total = np.sum(chi2_per_bin)
    p_value = 1 - chi2.cdf(chi2_total, df=ndof)
    return chi2_total / ndof, p_value