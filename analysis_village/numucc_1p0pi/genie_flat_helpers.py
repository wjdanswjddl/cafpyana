"""GENIE flat ROOT loading, integrated closure, and comparison plotting for numucc_1p0pi."""

from __future__ import annotations

import warnings
from os import path

import awkward as ak
import numpy as np
import pandas as pd
import uproot

from makedf.util import InFV
from pyanalib.pandas_helpers import pad_column_name
from analysis_village.unfolding.wienersvd import Matrix_Decomp
from pyanalib.variable_calculator import get_cc1p0pi_tki

from analysis_village.numucc_1p0pi.categories import IsNuInFV_NumuCC_1p0pi
from analysis_village.numucc_1p0pi.utils import (
    _covariance_per_bin_width,
    add_approval_text,
    add_chi2_text,
    add_genie_version_text,
    get_chi2,
    get_clipped_evts,
    get_textloc_x,
)

BRANCHES_NU = [
    "PDGnu", "cc", "Enu_true", "tgt", "ELep", "fScaleFactor", "RWWeight", "Mode",
    "Q2", "q0", "q3", "x", "y", "W_nuc_rest", "W", "W_genie", "flagCC0pi", "Weight",
    "CosThetaAdler", "PhiAdler", "dalphat", "dpt", "dphit", "pnreco_C",
]
BRANCHES_TRK = ["pdg", "px", "py", "pz", "E"]
BRANCHES_VTX = ["px_vert", "py_vert", "pz_vert"]

GIBUU_EXTRA_SCALE = 1000.0

GENIE_VAR_MAP = {
    "muon-p": ("mu_p", 1.0),
    "muon-dir_z": ("mu_dirz", 1.0),
    "proton-p": ("proton_p", 1.0),
    "proton-dir_z": ("proton_dirz", 1.0),
    "tki-del_Tp": ("del_Tp", 1.0),
    "tki-del_Tp_x": ("del_Tp_x", 1.0),
    "tki-del_Tp_y": ("del_Tp_y", 1.0),
    "tki-del_p": ("del_p", 1.0),
    "tki-del_alpha": ("del_alpha", 1.0),
    "tki-del_phi": ("del_phi", 1.0),
}


def genie_trk_arrays_to_df(trk_arr):
    counts = ak.num(trk_arr["pdg"])
    n_ev = len(counts)
    event_ids = np.repeat(np.arange(n_ev), counts)
    subentry = ak.flatten(ak.unflatten(np.arange(ak.sum(counts)), counts))
    flat = ak.zip(
        {
            "event": event_ids,
            "subentry": subentry,
            **{c: ak.flatten(trk_arr[c]) for c in BRANCHES_TRK},
        }
    )
    return ak.to_dataframe(flat).set_index(["event", "subentry"])


def prepare_genie_trk_df(trk):
    if isinstance(trk, pd.DataFrame):
        if isinstance(trk.index, pd.MultiIndex) and pd.api.types.is_numeric_dtype(trk["pdg"]):
            return trk
        trk_arr = ak.Array(
            {c: trk[c].array if hasattr(trk[c], "array") else trk[c].tolist() for c in BRANCHES_TRK}
        )
    else:
        trk_arr = trk
    return genie_trk_arrays_to_df(trk_arr)


def attach_genie_vtx_from_flat_branches(nu_df, vert_arr):
    nu_df = nu_df.copy()
    nu_df["vtx_x"] = ak.fill_none(ak.firsts(vert_arr["px_vert"]), np.nan)
    nu_df["vtx_y"] = ak.fill_none(ak.firsts(vert_arr["py_vert"]), np.nan)
    nu_df["vtx_z"] = ak.fill_none(ak.firsts(vert_arr["pz_vert"]), np.nan)
    return nu_df


def get_trk_info(nudf, trkdf, det="SBND"):
    """Multiplicity counts (|pdg|, p > threshold) + leading μ/p by highest momentum."""
    nudf = nudf.copy()
    trkdf = trkdf.copy()
    trkdf["_p"] = np.sqrt(trkdf.px**2 + trkdf.py**2 + trkdf.pz**2)
    ntrks = trkdf.pdg.groupby(level=[0]).count()
    nudf["ntrks"] = ntrks
    if det == "SBND":
        pid_pth = zip([13, 2212, 211, 111], ["mu", "proton", "pi", "pi0"], [0.22, 0.3, 0.07, 0])
    else:
        pid_pth = zip([13, 2212, 211, 111], ["mu", "proton", "pi", "pi0"], [0.22, 0.3, 0.07, 0])
    for pid, pname, pth in pid_pth:
        ntrks_pid = trkdf[(np.abs(trkdf.pdg) == pid) & (trkdf._p > pth)].pdg.groupby(level=[0]).count()
        nudf[f"n{pname}s"] = ntrks_pid.fillna(0)
        species = trkdf[trkdf.pdg == pid] if pname == "proton" else trkdf[np.abs(trkdf.pdg) == pid]
        leading = species.sort_values("_p", ascending=False).groupby(level=[0]).head(1)
        leading_p = leading["_p"]
        leading_p.name = f"{pname}_p"
        nudf = nudf.join(leading_p.reset_index(level=[1])[f"{pname}_p"])
        for comp, trk_col in zip(("x", "y", "z"), ("px", "py", "pz")):
            leading_dir = leading[trk_col] / leading_p
            leading_dir.name = f"{pname}_dir{comp}"
            nudf = nudf.join(leading_dir.reset_index(level=[1])[f"{pname}_dir{comp}"])
    return nudf


def get_genie_topology_kinematics_mask(df, det="SBND"):
    # nmus==1 / nprotons==1 encode the MeV thresholds; mu_p/proton_p are highest-momentum leading
    mu_kin = df.mu_p < 1.0
    p_kin = df.proton_p < 1.0
    topology = (
        (df.nmus == 1)
        & (df.nprotons == 1)  # np_300MeVc == 1; protons below 0.3 GeV allowed
        & (np.nan_to_num(df.npis, nan=0) == 0)
        & (np.nan_to_num(df.npi0s, nan=0) == 0)
    )
    return topology & mu_kin & p_kin


def get_genie_signal_mask(df):
    """νμ CC + 1μ1p0π topology (no vertex FV on GENIE flat files)."""
    return (df.PDGnu == 14) & (df.cc == 1) & get_genie_topology_kinematics_mask(df)


def get_genie_numu_cc_mask(df):
    return (df.PDGnu == 14) & (df.cc == 1)


def get_genie_vertex_fv_mask(df, genie_flat_fv_det, flat_aligned_x_cm=10.0, apply_x_strip=True):
    if not all(c in df.columns for c in ("vtx_x", "vtx_y", "vtx_z")):
        return None
    vtx_cm = pd.DataFrame(
        {
            "x": df["vtx_x"].to_numpy(dtype=float) * 100.0,
            "y": df["vtx_y"].to_numpy(dtype=float) * 100.0,
            "z": df["vtx_z"].to_numpy(dtype=float) * 100.0,
        }
    )
    mask = InFV(vtx_cm, det=genie_flat_fv_det)
    if apply_x_strip:
        mask = mask & (np.abs(vtx_cm.x) > flat_aligned_x_cm)
    return mask


def get_genie_closure_1p0pi_mask(df, genie_flat_fv_det, flat_aligned_x_cm=10.0):
    """GENIE flat 1p0π: vertex FV + νμ CC + 1μ1p0π topology (genie_vs_production_xsec).

    Topology matches production ``np_300MeVc == 1`` (one proton above 0.3 GeV; lower-p
    protons allowed) and ``mc.p.genp < 1`` on highest-momentum proton.
    """
    numu_cc = (df.PDGnu == 14) & (df.cc == 1)
    topo = get_genie_topology_kinematics_mask(df)
    fv = get_genie_vertex_fv_mask(df, genie_flat_fv_det, flat_aligned_x_cm)
    if fv is None:
        return None
    return numu_cc & topo & fv


def get_mcnu_numu_cc_mask(nudf):
    return (nudf.mc.pdg == 14) & (nudf.mc.iscc == 1)


def get_mcnu_vertex_fv_mask(nudf, genie_flat_fv_det, flat_aligned_x_cm=10.0, apply_x_strip=True):
    in_fv = InFV(nudf.mc.position, det=genie_flat_fv_det)
    if apply_x_strip:
        in_fv = in_fv & (np.abs(nudf.mc.position.x) > flat_aligned_x_cm)
    return in_fv


def get_production_signal_mask(nudf, signal_truth_fv="per_tpc"):
    return IsNuInFV_NumuCC_1p0pi(nudf, signal_truth_fv=signal_truth_fv)


def infer_genie_ref_pot(nudf, ref_pot):
    fsf = np.asarray(nudf["fScaleFactor"].unique(), dtype=float)
    if fsf.size != 1:
        warnings.warn(f"{fsf.size} unique fScaleFactor values; using GENIE_REF_POT.", stacklevel=2)
    return float(ref_pot)


def _genie_mu_p_frames(nu_df):
    mudf = pd.DataFrame(index=nu_df.index)
    pdf = pd.DataFrame(index=nu_df.index)
    for target, prefix in (("mudf", "mu"), ("pdf", "proton")):
        frame = mudf if target == "mudf" else pdf
        frame[("totp", "", "", "", "", "", "")] = nu_df[f"{prefix}_p"].to_numpy(dtype=float)
        frame[("dir", "x", "", "", "", "", "")] = nu_df[f"{prefix}_dirx"].to_numpy(dtype=float)
        frame[("dir", "y", "", "", "", "", "")] = nu_df[f"{prefix}_diry"].to_numpy(dtype=float)
        frame[("dir", "z", "", "", "", "", "")] = nu_df[f"{prefix}_dirz"].to_numpy(dtype=float)
        frame.columns = pd.MultiIndex.from_tuples(frame.columns)
    return mudf, pdf


def add_genie_tki_columns(nu_df):
    nu_df = nu_df.copy()
    need = [
        "mu_p", "mu_dirx", "mu_diry", "mu_dirz",
        "proton_p", "proton_dirx", "proton_diry", "proton_dirz",
    ]
    if not all(c in nu_df.columns for c in need):
        return nu_df
    mudf, pdf = _genie_mu_p_frames(nu_df)
    tki = get_cc1p0pi_tki(
        mudf, pdf,
        pad_column_name(("totp",), mudf),
        pad_column_name(("totp",), pdf),
    )
    for name, values in tki.items():
        nu_df[name] = np.asarray(values, dtype=float)
    return nu_df


def load_genie_flat_file(
    flat_path,
    data_pot,
    genie_ref_pot,
    attach_vertex=True,
    genie_flat_fv_det="SBND_nohighyz",
    flat_aligned_x_cm=10.0,
):
    events = uproot.open(flat_path + ":FlatTree_VARS")
    nu_df = events.arrays(BRANCHES_NU, library="pd")
    trk_df = prepare_genie_trk_df(events.arrays(BRANCHES_TRK, library="ak"))
    nu_df = get_trk_info(nu_df, trk_df)
    nu_df = add_genie_tki_columns(nu_df)
    if attach_vertex:
        vert_arr = events.arrays(BRANCHES_VTX, library="ak")
        nu_df = attach_genie_vtx_from_flat_branches(nu_df, vert_arr)
    ref_pot = infer_genie_ref_pot(nu_df, genie_ref_pot)
    pot_scale = data_pot / ref_pot
    if "GiBUU" in flat_path or "gibuu" in flat_path.lower():
        pot_scale = pot_scale / GIBUU_EXTRA_SCALE
    numu_cc_mask = get_genie_numu_cc_mask(nu_df)
    sig_mask = get_genie_signal_mask(nu_df)
    closure_1p0pi_mask = None
    if attach_vertex:
        closure_1p0pi_mask = get_genie_closure_1p0pi_mask(
            nu_df, genie_flat_fv_det, flat_aligned_x_cm
        )
    return {
        "nu_df": nu_df,
        "sig_mask": sig_mask,
        "numu_cc_mask": numu_cc_mask,
        "closure_1p0pi_mask": closure_1p0pi_mask,
        "ref_pot": ref_pot,
        "pot_scale": pot_scale,
        "flat_path": flat_path,
    }


def build_genie_flat_cache(files, data_pot, genie_ref_pot, **load_kw):
    if isinstance(files, (list, tuple)):
        files = {path.basename(p).replace(".flat.root", ""): p for p in files}
    cache = {}
    for label, flat_path in files.items():
        print(f"Loading {label}: {flat_path}")
        pack = load_genie_flat_file(flat_path, data_pot, genie_ref_pot, **load_kw)
        n_sig = int(pack["sig_mask"].sum())
        n_cc = int(pack["numu_cc_mask"].sum())
        n_fv = (
            int(pack["closure_1p0pi_mask"].sum())
            if pack.get("closure_1p0pi_mask") is not None
            else None
        )
        print(
            f"  entries={len(pack['nu_df']):,}  νμCC={n_cc:,}  signal(topo)={n_sig:,}"
            + (f"  1p0π(FV+topo)={n_fv:,}" if n_fv is not None else "")
            + f"  GENIE_REF_POT={pack['ref_pot']:.3e}  pot_scale={pack['pot_scale']:.3e}"
        )
        cache[label] = pack
    return cache


def genie_flat_integrated_sigma(nu_df, mask, pot_scale):
    m = np.asarray(mask, dtype=bool)
    w = (
        40.0
        * nu_df.loc[m, "fScaleFactor"].to_numpy(dtype=float)
        * nu_df.loc[m, "Weight"].to_numpy(dtype=float)
    )
    return float(w.sum() * float(pot_scale))


def mcnu_integrated_sigma(mc_nu_df, mask, xsec_unit):
    m = np.asarray(mask, dtype=bool)
    return float(mc_nu_df.loc[m, "pot_weight"].sum() * xsec_unit)


def print_integrated_closure(label, sigma_mc, sigma_genie, n_mc=None, n_genie=None):
    ratio = sigma_genie / sigma_mc if sigma_mc else np.nan
    print(f"  [{label}]")
    if n_mc is not None and n_genie is not None and n_genie:
        print(f"    N_mc={n_mc:,}  N_genie={n_genie:,}  N_mc/N_genie={n_mc/n_genie:.4f}")
    print(f"    ∫σ_mc   = {sigma_mc:.6e}")
    print(f"    ∫σ_flat = {sigma_genie:.6e}")
    print(f"    flat/MC = {ratio:.4f}")


def genie_flat_differential_xsec(nu_df, sig_mask, var_config, pot_scale):
    branch, scale = GENIE_VAR_MAP[var_config.var_save_name]
    if branch not in nu_df.columns:
        raise KeyError(f"GENIE branch '{branch}' missing for {var_config.var_save_name}")
    sig = nu_df[sig_mask]
    var = np.clip(sig[branch].to_numpy(dtype=float) * scale, var_config.bins[0], var_config.bins[-1])
    wgt = 40.0 * sig["fScaleFactor"].to_numpy(dtype=float) * sig["Weight"].to_numpy(dtype=float)
    wgt = wgt * float(pot_scale)
    bins = var_config.bins
    sigma_bin, _ = np.histogram(var, bins=bins, weights=wgt)
    bin_widths = np.diff(bins)
    if len(bins) == 2:
        bin_widths = np.array([1.0])
    return sigma_bin, sigma_bin / bin_widths


def apply_ac_differential_xsec(sigma_bin, unfold, bin_widths):
    smeared = np.asarray(unfold["AddSmear"], dtype=float) @ np.asarray(sigma_bin, dtype=float)
    return smeared / bin_widths


def mc_truth_differential_xsec(ret_signal_hists, var_config, xsec_unit):
    model = np.asarray(ret_signal_hists["nevts_allmc"], dtype=float) * xsec_unit
    bin_widths = np.diff(var_config.bins)
    if len(var_config.bins) == 2:
        bin_widths = np.array([1.0])
    return model, model / bin_widths


def unfolded_differential_xsec(unfold, var_config):
    unfolded = np.asarray(unfold["unfold"], dtype=float)
    bin_widths = np.diff(var_config.bins)
    if len(var_config.bins) == 2:
        bin_widths = np.array([1.0])
    return unfolded, unfolded / bin_widths


def integrated_xsec_from_differential(sigma_bin):
    return float(np.asarray(sigma_bin, dtype=float).sum())


def run_integrated_closure_checks(
    mc_nu_df,
    genie_flat_cache,
    xsec_unit,
    data_tot_pot,
    genie_ref_pot,
    mc_pot_scale,
    signal_truth_fv="per_tpc",
    genie_flat_fv_det="SBND_nohighyz",
    flat_aligned_x_cm=10.0,
    unfold_cache=None,
):
    """Sanity checks mirroring genie_vs_production_xsec integrated σ recipe."""
    mc_fv = get_mcnu_vertex_fv_mask(mc_nu_df, genie_flat_fv_det, flat_aligned_x_cm)
    mc_numu_cc_fv = get_mcnu_numu_cc_mask(mc_nu_df) & mc_fv
    mc_numu_cc_all = get_mcnu_numu_cc_mask(mc_nu_df)
    mc_sig_unfold = get_production_signal_mask(mc_nu_df, signal_truth_fv)

    ar23 = genie_flat_cache.get("GENIE AR23 (CC)")
    if ar23 is None:
        ar23 = next(iter(genie_flat_cache.values()))
    nu_df = ar23["nu_df"]
    ps = ar23["pot_scale"]

    print("=== Integrated σ closure (AR23 flat vs production mcnu) ===")
    print(f"data_tot_pot={data_tot_pot:.3e}  GENIE_REF_POT={genie_ref_pot:.3e}  mc_pot_scale={mc_pot_scale:.3e}")
    print(f"xsec_unit={xsec_unit:.6e} cm²/nucleon\n")

    sigma_mc_cc_fv = mcnu_integrated_sigma(mc_nu_df, mc_numu_cc_fv, xsec_unit)
    sigma_g_cc = genie_flat_integrated_sigma(nu_df, ar23["numu_cc_mask"], ps)
    print_integrated_closure(
        "νμ CC (MC vtx FV + |x|>10 cm; GENIE all flat CC)",
        sigma_mc_cc_fv,
        sigma_g_cc,
        n_mc=int(mc_numu_cc_fv.sum()),
        n_genie=int(ar23["numu_cc_mask"].sum()),
    )

    sigma_mc_cc_all = mcnu_integrated_sigma(mc_nu_df, mc_numu_cc_all, xsec_unit)
    print_integrated_closure(
        "νμ CC (MC all mcnu rows, no vtx FV cut)",
        sigma_mc_cc_all,
        sigma_g_cc,
        n_mc=int(mc_numu_cc_all.sum()),
        n_genie=int(ar23["numu_cc_mask"].sum()),
    )

    sigma_mc_sig = mcnu_integrated_sigma(mc_nu_df, mc_sig_unfold, xsec_unit)
    sigma_g_topo = genie_flat_integrated_sigma(nu_df, ar23["sig_mask"], ps)
    print_integrated_closure(
        f"1p0π signal (MC IsNuInFV_NumuCC_1p0pi {signal_truth_fv}; GENIE topo, no vtx FV)",
        sigma_mc_sig,
        sigma_g_topo,
        n_mc=int(mc_sig_unfold.sum()),
        n_genie=int(ar23["sig_mask"].sum()),
    )

    if ar23.get("closure_1p0pi_mask") is not None:
        sigma_g_fv_topo = genie_flat_integrated_sigma(nu_df, ar23["closure_1p0pi_mask"], ps)
        print_integrated_closure(
            "1p0π [optional] GENIE vtx FV + topo (not used for plots)",
            sigma_mc_sig,
            sigma_g_fv_topo,
            n_mc=int(mc_sig_unfold.sum()),
            n_genie=int(ar23["closure_1p0pi_mask"].sum()),
        )
    else:
        print("  [GENIE vtx+topo diagnostic] skipped — vtx branches missing on flat file")

    if unfold_cache:
        from analysis_village.numucc_1p0pi.variable_configs import VariableConfig

        vc = VariableConfig.muon_momentum()
        if vc.var_save_name in unfold_cache:
            cache = unfold_cache[vc.var_save_name]
            sigma_mc_bins, _ = mc_truth_differential_xsec(
                cache["ret_signal_hists"], vc, xsec_unit
            )
            sigma_g_bins, _ = genie_flat_differential_xsec(nu_df, ar23["sig_mask"], vc, ps)
            print("\n=== Binned sum check (μ momentum, 1p0π topology on flat, no vtx FV) ===")
            print(f"  ∫σ_mc (sum bins, unfold signal_hists) = {integrated_xsec_from_differential(sigma_mc_bins):.6e}")
            print(f"  ∫σ_flat (sum bins)                   = {integrated_xsec_from_differential(sigma_g_bins):.6e}")
            s_mc = integrated_xsec_from_differential(sigma_mc_bins)
            print(f"  flat/MC = {integrated_xsec_from_differential(sigma_g_bins)/s_mc:.4f}")

    w_cc_1e20 = (
        40.0
        * nu_df.loc[ar23["numu_cc_mask"], "fScaleFactor"]
        * nu_df.loc[ar23["numu_cc_mask"], "Weight"]
    ).sum()
    print(
        f"\nGENIE sum(40×fSF×W) νμ CC @ GENIE_REF_POT = {w_cc_1e20:.6e} "
        "(genie_vs_production_xsec)"
    )


def plot_generator_comparison(
    var_config,
    cache,
    genie_flat_cache,
    xsec_unit,
    save_fig=False,
    save_fig_dir=".",
    genie_cache=None,
    show_ratio=False,
    show_mc_truth=True,
    apply_ac=True,
    genie_scale=1.0,
    show_chi2=True,
    show_title=True,
    plot_labels=None,
    textloc=(0.05, 0.55),
    approval="internal",
    fig_suffix="",
):
    import matplotlib.pyplot as plt

    genie_cache = genie_flat_cache if genie_cache is None else genie_cache
    unfold = cache["unfold"]
    measured = cache["measured"]
    model = cache["model"]

    bins = var_config.bins
    centers = var_config.bin_centers
    bin_widths = np.diff(bins)
    if len(bins) == 2:
        bin_widths = np.array([1.0])

    Unfolded = np.asarray(unfold["unfold"], dtype=float)
    Unfolded_perwidth = Unfolded / bin_widths
    cov_unfold_perwidth = _covariance_per_bin_width(unfold["UnfoldCov"], bin_widths)

    UnfoldCov_stat = unfold["StatUnfoldCov"]
    Unfold_uncert_stat = np.diag(UnfoldCov_stat)
    UnfoldCov_syst = unfold["SystUnfoldCov"]
    SystUnfoldCov_norm, _, SystUnfoldCov_shape = Matrix_Decomp(model, UnfoldCov_syst)
    Unfold_uncert_norm = np.sqrt(np.abs(np.diag(SystUnfoldCov_norm)))
    Unfold_uncert_shape = np.sqrt(np.abs(np.diag(SystUnfoldCov_shape)))
    Unfold_uncert_stat_shape_perwidth = Unfold_uncert_shape / bin_widths
    Unfold_uncert_norm_perwidth = Unfold_uncert_norm / bin_widths

    Data_frac_unc = 1.0 / np.sqrt(np.maximum(measured / xsec_unit, 1e-300))
    if len(bins) == 2:
        tot_err = np.sqrt(
            (Unfolded_perwidth * Data_frac_unc) ** 2 + Unfold_uncert_norm_perwidth**2
        )
    else:
        tot_err = np.sqrt(
            (Unfolded_perwidth * Data_frac_unc) ** 2 + Unfold_uncert_stat_shape_perwidth**2
        )

    models = {}
    if show_mc_truth:
        models["GENIE"] = np.asarray(model, dtype=float)
    for glabel, gpack in genie_cache.items():
        try:
            sigma_bin, _ = genie_flat_differential_xsec(
                gpack["nu_df"], gpack["sig_mask"], var_config, gpack["pot_scale"]
            )
            models[glabel] = np.asarray(sigma_bin, dtype=float)
        except KeyError as ex:
            print(f"  skip {glabel}: {ex}")

    nrows = 2 if show_ratio else 1
    fig, axes = plt.subplots(nrows, 1, figsize=(8.5, 7 if nrows == 1 else 9), sharex=True)
    if nrows == 1:
        axes = [axes]
    ax = axes[0]

    bar_handle = ax.errorbar(
        centers, Unfolded_perwidth,
        yerr=Unfold_uncert_stat_shape_perwidth,
        fmt="o", color="black", capsize=3,
    )
    norm_handle = None
    if len(bins) != 2:
        norm_handle = ax.bar(
            centers, Unfold_uncert_norm_perwidth, width=bin_widths,
            label="Syst. error (norm)", alpha=0.5, color="gray",
        )
    ax.errorbar(centers, Unfolded_perwidth, yerr=tot_err, fmt="o", color="black", capsize=3)

    add_smear = unfold["AddSmear"]
    model_handles = []
    model_labels = []
    pw_for_ylim = [Unfolded_perwidth]

    for mkey, mvec in models.items():
        mvec = np.asarray(mvec, dtype=float)
        m_smeared = add_smear @ mvec if apply_ac else mvec
        scale = 1.0 if mkey == "GENIE" else float(genie_scale)
        m_perwidth = m_smeared / bin_widths * scale
        pw_for_ylim.append(m_perwidth)

        label = mkey
        if show_chi2:
            mask = (Unfolded_perwidth > 0) & (m_perwidth > 0)
            if np.any(mask):
                chi2_val, _ = get_chi2(
                    Unfolded_perwidth[mask], m_perwidth[mask],
                    cov_unfold_perwidth[np.ix_(mask, mask)],
                )
                label = f"{mkey} ($\\chi^2$/ndof = {chi2_val:.1f}/{int(np.sum(mask))})"
            else:
                label = f"{mkey} ($\\chi^2$/ndof = n/a)"

        mh, = ax.step(
            bins, np.append(m_perwidth, m_perwidth[-1]), where="post", lw=2,
        )
        model_handles.append(mh)
        model_labels.append(label)

    if len(bins) == 2:
        handles = [bar_handle] + model_handles
        labels = ["Data (Syst. Unc. + Stat. Unc.)"] + model_labels
    else:
        handles = [bar_handle, norm_handle] + model_handles
        labels = ["Data (Shape Syst. Unc. + Stat. Unc.)", "Norm. Syst. Unc."] + model_labels

    ax.legend(handles, labels, loc="upper left", fontsize=12, frameon=False,
              ncol=1, bbox_to_anchor=(0.02, 0.98))
    ax.set_xlim(bins[0], bins[-1])
    ymax = np.nanmax(np.concatenate([np.ravel(x) for x in pw_for_ylim if x is not None]))
    ax.set_ylim(0.0, ymax * 1.7)
    ax.set_xlabel(var_config.var_labels[0], fontsize=20)
    ax.set_ylabel(var_config.xsec_label, fontsize=20)
    if show_title:
        title = (plot_labels or ["", "", ""])[2] or var_config.var_plot_name
        ax.set_title(title, fontsize=20)

    textloc_x, textloc_ha = get_textloc_x(Unfolded_perwidth, bins, list(textloc))
    textloc_y = textloc[1]
    add_approval_text(approval, textloc_x, textloc_y, textloc_ha)
    add_genie_version_text(textloc_x, textloc_y - 0.1, textloc_ha)

    if show_ratio and show_mc_truth and "GENIE" in models:
        genie_pw = (add_smear @ models["GENIE"] if apply_ac else models["GENIE"]) / bin_widths
        axr = axes[1]
        axr.axhline(1.0, color="gray", ls=":", lw=1)
        ratio = Unfolded_perwidth / np.maximum(genie_pw, 1e-40)
        axr.step(bins, np.append(ratio, ratio[-1]), where="post", lw=2, label="Data / GENIE")
        axr.set_ylabel("Ratio to GENIE")
        axr.set_xlabel(var_config.var_labels[0], fontsize=20)
        axr.legend(loc="best", fontsize=10, frameon=False)
        axr.set_ylim(0.5, 1.5)

    fig.tight_layout()
    if save_fig:
        suffix = fig_suffix or ("-ratio" if show_ratio else "")
        out = path.join(save_fig_dir, f"{var_config.var_save_name}-genie_generator_comparison{suffix}.pdf")
        fig.savefig(out, bbox_inches="tight", dpi=300)
        print("saved", out)
    plt.show()
