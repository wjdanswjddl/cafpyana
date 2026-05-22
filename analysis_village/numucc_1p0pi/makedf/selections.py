import pandas as pd
import numpy as np
from pyanalib.variable_calculator import *
from pyanalib.pandas_helpers import *
from makedf.constants import *
from makedf.util import *
from analysis_village.numucc_1p0pi.categories import PER_TPC_INCATHODE_CM


# ==== events selection cuts ====
# slice cuts
NU_SCORE_TH   = 0.45
SAVE_NTRKS    = 2
# track quality cuts
TRACKSCORE_TH = 0.5
VTXDIST_TH    = 1.2
# pid cuts
MU_CHI2MU_TH  = 25
MU_CHI2P_TH   = 100
MU_LEN_TH     = 50
QUAL_TH       = 0.2
P_CHI2P_TH    = 90
P_LEN_TH      = 0
# kinematic cuts
MU_PLO_TH     = 0.22
MU_PHI_TH     = 1
P_PLO_TH      = 0.3
P_PHI_TH      = 1


def cut_clear_cosmic(df):
    return df[df.slc.is_clear_cosmic == 0]

def cut_vertex_in_fv(df, det="SBND"):
    return df[InFV(df.slc.vertex, det=det)]

def cut_nu_score(df, th=0.5):
    return df[df.slc.nu_score > th]

def get_valid_trks(df):
    return df[df.pfp.trk.producer != 4294967295]

def cut_good_trks(trkdf):
    mask = (trkdf.pfp.trk.len > 0) &\
         (trkdf.pfp.pfochar.vtxdist < 100) #&\
    return trkdf[mask]


def _dedupe_event_level_index(trk_df: pd.DataFrame, n_event_levels: int) -> pd.DataFrame:
    """Drop duplicate slice keys before merging tracks onto ``evt``.

    ``groupby(...).nth(i)`` can return multiple rows per slice when track ordering
    ties; ``multicol_merge(..., validate='one_to_one')`` then raises and ``trk1`` /
    ``trk2`` never get attached.
    """
    if trk_df is None or len(trk_df) == 0:
        return trk_df
    if trk_df.index.duplicated().any():
        return trk_df[~trk_df.index.duplicated(keep="first")]
    return trk_df


def _is_merged_trk_block_name(name) -> bool:
    """True for per-slice track blocks added by :func:`get_trk_info` (incl. merge suffixes)."""
    s = str(name)
    if s in ("mu", "p"):
        return True
    if s.startswith("nocut_trk") or (s.startswith("trk") and len(s) > 3 and s[3].isdigit()):
        return True
    return False


def evt_has_trk1_trk2(evtdf: pd.DataFrame) -> bool:
    """True when per-slice ``trk1`` / ``trk2`` column blocks are present on ``evt``."""
    if evtdf is None or len(evtdf) == 0:
        return False
    try:
        top = evtdf.columns.get_level_values(0).unique()
    except Exception:
        return False
    return ("trk1" in top) and ("trk2" in top)


def _drop_merged_trk_blocks(evtdf: pd.DataFrame) -> pd.DataFrame:
    """Remove prior ``trk*`` / ``nocut_trk*`` / ``mu`` / ``p`` blocks before re-merging tracks.

    Re-running :func:`get_trk_info` without this leaves duplicate top-level names; pandas
    then suffixes columns (``trk1_x``) and ``evt.trk1`` attribute access breaks.
    """
    if evtdf is None or len(evtdf.columns) == 0:
        return evtdf
    if not isinstance(evtdf.columns, pd.MultiIndex):
        return evtdf
    lev0 = evtdf.columns.get_level_values(0)
    keep = ~pd.Index(lev0).map(_is_merged_trk_block_name)
    if keep.all():
        return evtdf
    return evtdf.loc[:, keep]


def get_trk_info(evtdf, trkdf, save_ntrks=3):
    evtdf = _drop_merged_trk_blocks(evtdf)
    nlevels = len(trkdf.index.names)
    ntrks = trkdf.pfp.id.groupby(level=list(range(nlevels-1))).count()
    ntrks.reindex(evtdf.index, fill_value=0)
    evtdf.loc[:, "n_trks"] = ntrks.copy()

    good_trks = cut_good_trks(trkdf).copy()
    ntrks = good_trks.pfp.id.groupby(level=list(range(nlevels-1))).count()
    ntrks.reindex(evtdf.index, fill_value=0)
    evtdf.loc[:, "n_good_trks"] = ntrks.copy()

    trks_sorted = trkdf.sort_values(by=('pfp','trk','len'), ascending=False)
    good_trks_sorted = good_trks.sort_values(by=('pfp','trk','len'), ascending=False)
    # get 'ntrks' longest tracks
    evt_levels = nlevels - 1
    for i in range(save_ntrks):
        trk_i = good_trks_sorted.groupby(level=list(range(evt_levels))).nth(i)
        trk_i.columns = pd.MultiIndex.from_tuples([tuple(["nocut_trk" + str(i+1)] + list(c)) for c in trk_i.columns])
        trk_i = _dedupe_event_level_index(trk_i.droplevel(-1), evt_levels)
        evtdf = multicol_merge(evtdf, trk_i, left_index=True, right_index=True, how="left", validate="one_to_one")

        good_trk_i = good_trks_sorted.groupby(level=list(range(evt_levels))).nth(i)
        good_trk_i.columns = pd.MultiIndex.from_tuples([tuple(["trk" + str(i+1)] + list(c)) for c in good_trk_i.columns])
        good_trk_i = _dedupe_event_level_index(good_trk_i.droplevel(-1), evt_levels)
        evtdf = multicol_merge(evtdf, good_trk_i, left_index=True, right_index=True, how="left", validate="one_to_one")

    return evtdf


def cut_2prong(df):
    # return df[(df.n_good_trks == 2) & (df.n_trks <= 3)]
    return df[(df.n_good_trks == 2)]

def cut_2prong_contained(df, det="SBND"):
    if det == "SBND_Gen1":
        in_TPC1_cut = InFV(df.slc.vertex, det="SBND_TPC1", incathode=PER_TPC_INCATHODE_CM) \
                & InFV(df.trk1.pfp.trk.end, det="SBND_TPC1", incathode=PER_TPC_INCATHODE_CM) \
                & InFV(df.trk2.pfp.trk.end, det="SBND_TPC1", incathode=PER_TPC_INCATHODE_CM)
        in_TPC2_cut = InFV(df.slc.vertex, det="SBND_TPC2", incathode=PER_TPC_INCATHODE_CM) \
                & InFV(df.trk1.pfp.trk.end, det="SBND_TPC2", incathode=PER_TPC_INCATHODE_CM) \
                & InFV(df.trk2.pfp.trk.end, det="SBND_TPC2", incathode=PER_TPC_INCATHODE_CM)
        perTPC_cut = in_TPC1_cut | in_TPC2_cut
        return df[perTPC_cut]

    else:
        return df[InFV(df.trk1.pfp.trk.start, det=det) & InFV(df.trk1.pfp.trk.end, det=det) \
            & InFV(df.trk2.pfp.trk.start, det=det) & InFV(df.trk2.pfp.trk.end, det=det)]

def cut_2prong_trackscore(df, trackscore_th=0.5):
    return df[(df.trk1.pfp.trackScore > trackscore_th) & (df.trk2.pfp.trackScore > trackscore_th)]

def cut_2prong_vtxdist(df, vtxdist_th=1.5):
    return df[(df.trk1.pfp.pfochar.vtxdist < vtxdist_th) & (df.trk2.pfp.pfochar.vtxdist < vtxdist_th)]

def get_mu_p_candidate(df, 
                       mu_chi2mu_th=30, mu_chi2p_th=100, mu_len_th=50, qual_th=0.25,
                       p_chi2mu_th=30, p_chi2p_th=90, p_len_th=0, score_tag=""):

    nlevels = len(df.index.names)

    if not evt_has_trk1_trk2(df):
        raise KeyError(
            "evt is missing trk1/trk2 columns required for mu/p PID — "
            "call get_trk_info(evt, trk) after matching tracks to the current slice table"
        )
    trks = pd.concat([df.trk1, df.trk2])

    chimu_avg = avg_chi2(trks, f"chi2_muon{score_tag}")
    chip_avg = avg_chi2(trks, f"chi2_proton{score_tag}")

    mcs_range_diff = np.abs((trks.pfp.trk.rangeP.p_muon - trks.pfp.trk.mcsP.fwdP_muon) / trks.pfp.trk.rangeP.p_muon)

    mu_cut = (chimu_avg > 0) & (chimu_avg < mu_chi2mu_th) & \
            (chip_avg > mu_chi2p_th) & \
            (trks.pfp.trk.len > mu_len_th) & \
            (mcs_range_diff < qual_th)

    mu_candidate = trks[mu_cut]
    mu_candidate = mu_candidate.groupby(level=list(range(nlevels))).nth(0)

    mu_candidate.columns = pd.MultiIndex.from_tuples([tuple(["mu"] + list(c)) for c in mu_candidate.columns])
    df = multicol_merge(df, mu_candidate, left_index=True, right_index=True, how="left", validate="one_to_one")

    # TODO: keep & use original trk index?
    not_mu_candidate = pd.concat([trks[~mu_cut], trks[mu_cut].groupby(level=list(range(nlevels))).nth(1)])
    chip_avg = avg_chi2(not_mu_candidate, f"chi2_proton{score_tag}")
    p_candidate = not_mu_candidate[(chip_avg > 0) & (chip_avg < p_chi2p_th) & (not_mu_candidate.pfp.trk.len > p_len_th)]
    p_candidate = p_candidate.groupby(level=list(range(nlevels))).nth(0)

    p_candidate.columns = pd.MultiIndex.from_tuples([tuple(["p"] + list(c)) for c in p_candidate.columns])
    df = multicol_merge(df, p_candidate, left_index=True, right_index=True, how="left", validate="one_to_one")

    return df

def cut_has_mu(df):
    return df[~np.isnan(df.mu.pfp.trk.producer)]

def cut_has_p(df):
    return df[~np.isnan(df.p.pfp.trk.producer)]

def cut_mu_kinematics(df, mu_Plo_th=0.22, mu_Phi_th= 1):
    return df[(df.mu.pfp.trk.rangeP.p_muon > mu_Plo_th) & (df.mu.pfp.trk.rangeP.p_muon < mu_Phi_th)]

def cut_p_kinematics(df, p_Plo_th=0.3, p_Phi_th= 1):
    return df[(df.p.pfp.trk.rangeP.p_proton > p_Plo_th) & (df.p.pfp.trk.rangeP.p_proton < p_Phi_th)]