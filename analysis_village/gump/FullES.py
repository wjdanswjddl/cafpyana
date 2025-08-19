import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from plot_tools import *
# Add the head direcoty to sys.path
workspace_root = os.getcwd()  
sys.path.insert(0, workspace_root + "/../../")

# import this repo's classes
from makedf.util import *

import kinematics

def load_data(file, start=0, stop=1000):
    """Load event, header, and mcnu data from HDF file."""
    df_evt = pd.read_hdf(file, "evt", start=start, stop=stop)
    df_hdr = pd.read_hdf(file, "hdr", start=start, stop=stop)
    df_mcnu = pd.read_hdf(file, "mcnu", start=start, stop=stop)
    return df_evt, df_hdr, df_mcnu

def scale_pot(df, df_hdr, desired_pot):
    """Scale DataFrame by desired POT."""
    pot = sum(df_hdr.pot.tolist())
    print(f"POT: {pot}\nScaling to: {desired_pot}")
    scale = desired_pot / pot
    df['glob_scale'] = scale
    return pot, scale


def fv_cut(df_nd, df_fd):
    plot_cut(df_nd.slc.vertex.x, [SBNDFVCuts['x']['min'], SBNDFVCuts['x']['max']],
             df_fd.slc.vertex.x, [ICARUSFVCuts['C0']['x']['min'], ICARUSFVCuts['C0']['x']['max'],
                                 ICARUSFVCuts['C1']['x']['min'], ICARUSFVCuts['C1']['x']['max']],
             "Vertex X", "FVCutVarX.png", [-300, 300])
    plot_cut(df_nd.slc.vertex.y, [SBNDFVCuts['y']['min'], SBNDFVCuts['y']['max']],
             df_fd.slc.vertex.y, [ICARUSFVCuts['C0']['y']['min'], ICARUSFVCuts['C0']['y']['max'],
                                 ICARUSFVCuts['C1']['y']['min'], ICARUSFVCuts['C1']['y']['max']],
             "Vertex Y", "FVCutVarY.png", [-250, 250])
    plot_cut(df_nd.slc.vertex.z, [SBNDFVCuts['z']['min'], SBNDFVCuts['z']['max']],
             df_fd.slc.vertex.z, [ICARUSFVCuts['C0']['z']['min'], ICARUSFVCuts['C0']['z']['max'],
                                 ICARUSFVCuts['C1']['z']['min'], ICARUSFVCuts['C1']['z']['max']],
             "Vertex Z", "FVCutVarZ.png", [-1000, 1000])
    df_nd = df_nd[SelFV(df_nd.slc.vertex, "SBND")]
    df_fd = df_fd[SelFV(df_fd.slc.vertex, "ICARUS")]
    return df_nd, df_fd

def cosmic_cut(df_nd, df_fd):
    plot_nuscore_cut([
        df_nd.slc.nu_score[is_cosmic(df_nd)], df_nd.slc.nu_score[np.invert(is_cosmic(df_nd))]
    ], [0.5], [
        df_fd.slc.nu_score[is_cosmic(df_fd)], df_fd.slc.nu_score[np.invert(is_cosmic(df_fd))]
    ], [0.5], "NuScore", "NuScore.png", [0, 1],
        ['SBND Cosmic', 'SBND Neutrino', 'ICARUS Cosmic', 'ICARUS Neutrino'], 'right', 'select as neutrino')
    df_nd = df_nd[df_nd.slc.nu_score > 0.5]
    df_fd = df_fd[df_fd.slc.nu_score > 0.5]
    return df_nd, df_fd

def twoprong_cut(df_nd, df_fd):
    df_nd = df_nd[np.isnan(df_nd.other_shw_length) & np.isnan(df_nd.other_trk_length)]
    df_fd = df_fd[np.isnan(df_fd.other_shw_length) & np.isnan(df_fd.other_trk_length)]
    return df_nd, df_fd

def pid_cut(df_nd, df_fd):

    plot_PID_cut([[df_nd.p.pfp.trk.chi2pid.I2.chi2_muon, df_nd.mu.pfp.trk.chi2pid.I2.chi2_muon], 
                  [df_fd.p.pfp.trk.chi2pid.I2.chi2_muon, df_fd.mu.pfp.trk.chi2pid.I2.chi2_muon]], 
                  [25, 25], "MuScore", "MuScore.png", [0, 65],
                  ['SBND, Protons', 'SBND, Muons', 'ICARUS, Protons', 'ICARUS, Muons'], 
                  'left', 'select as muon')

    plot_PID_cut([[df_nd.p.pfp.trk.chi2pid.I2.chi2_proton, df_nd.mu.pfp.trk.chi2pid.I2.chi2_proton], 
                  [df_fd.p.pfp.trk.chi2pid.I2.chi2_proton, df_fd.mu.pfp.trk.chi2pid.I2.chi2_proton]],
                  [100, 100], "ProtonScore", "ProtonScore.png", [0, 200],
                  ['SBND, Protons', 'SBND, Muons', 'ICARUS, Protons', 'ICARUS, Muons'], 
                  'right', 'select as muon')

    plot_cut(df_nd.mu.pfp.trk.len, [50], df_fd.mu.pfp.trk.len, [50], "MuLen", "MuonMuLen.png", [0, 500])
    MUSEL_MUSCORE_TH, MUSEL_PSCORE_TH, MUSEL_LEN_TH = 25, 100, 50
    mu_cut_nd = (df_nd.mu.pfp.trk.chi2pid.I2.chi2_muon < MUSEL_MUSCORE_TH) & \
                (df_nd.mu.pfp.trk.chi2pid.I2.chi2_proton > MUSEL_PSCORE_TH) & \
                (df_nd.mu.pfp.trk.len > MUSEL_LEN_TH)
    mu_cut_fd = (df_fd.mu.pfp.trk.chi2pid.I2.chi2_muon < MUSEL_MUSCORE_TH) & \
                (df_fd.mu.pfp.trk.chi2pid.I2.chi2_proton > MUSEL_PSCORE_TH) & \
                (df_fd.mu.pfp.trk.len > MUSEL_LEN_TH)
    PSEL_MUSCORE_TH, PSEL_PSCORE_TH = 0, 90
    p_cut_nd = (df_nd.p.pfp.trk.chi2pid.I2.chi2_muon > PSEL_MUSCORE_TH) & \
               (df_nd.p.pfp.trk.chi2pid.I2.chi2_proton < PSEL_PSCORE_TH)
    p_cut_fd = (df_fd.p.pfp.trk.chi2pid.I2.chi2_muon > PSEL_MUSCORE_TH) & \
               (df_fd.p.pfp.trk.chi2pid.I2.chi2_proton < PSEL_PSCORE_TH)
    df_nd = df_nd.loc[mu_cut_nd & p_cut_nd]
    df_fd = df_fd.loc[mu_cut_fd & p_cut_fd]
    return df_nd, df_fd

def stub_cut(df_nd, df_fd):
    l_cuts = [5.5e5, 4e5, 3e5, 0]
    plot_stub_2d(df_nd, l_cuts, "SBND.png", "SBND Stub Cut")
    plot_stub_2d(df_nd[np.invert(is_signal(df_nd, "SBND"))], l_cuts, "FPSBND.png", "False Positive SBND Stub Cut")
    plot_stub_2d(df_fd, l_cuts, "ICARUS.png", "ICARUS Stub Cut")
    plot_stub_2d(df_fd[np.invert(is_signal(df_fd, "ICARUS"))], l_cuts, "FPICARUS.png", "False Positive ICARUS Stub Cut")
    for i, l in enumerate(["0_5", "1", "2", "3"]):
        nd_cut = np.invert((df_nd.stub[f"l{l}cm"].Q / df_nd.stub[f"l{l}cm"].length) > l_cuts[i])
        fd_cut = np.invert((df_fd.stub[f"l{l}cm"].Q / df_fd.stub[f"l{l}cm"].length) > l_cuts[i])
        df_nd = df_nd.loc[nd_cut]
        df_fd = df_fd.loc[fd_cut]
    return df_nd, df_fd


def save_output(df_nd, df_fd, df_nd_mcnu, df_fd_mcnu, nd_POT, fd_POT, dffile_nd, dffile_fd):
    for df, mcnu, pot, file, label in [
        (df_nd, df_nd_mcnu, nd_POT, dffile_nd, "SBND"),
        (df_fd, df_fd_mcnu, fd_POT, dffile_fd, "ICARUS")
    ]:
        ret = recodf(df)
        savefile = file.split(".")[0] + ".gump.df"
        ret = tmatchdf(ret, mcnu, label)
        weightcol = [c for c in ret.columns if is_weightcol(c[0])]
        varcol = [c for c in ret.columns if not is_weightcol(c[0])]
        ret[varcol].droplevel(1, 1).reset_index(drop=True).to_hdf(savefile, key="var")
        ret[weightcol].reset_index(drop=True).to_hdf(savefile, key="wgt")
        pd.DataFrame([pot]).to_hdf(savefile, key="pot")

def apply_cuts(df_nd, df_fd, df_nd_mcnu, df_fd_mcnu, nd_POT, fd_POT, dffile_nd, dffile_fd):
    comp = []
    # FV Cut
    df_nd, df_fd = fv_cut(df_nd, df_fd)
    comp.append(plot_int(df_nd, df_fd, 'del_p', "0. FV Cut", "int0delp_fv.png", r"$\delta p$ [GeV/c]", eff_bool=True))
    plot_int(df_nd, df_fd, 'del_Tp', "0. FV Cut", "int0delpT_fv.png", r"$\delta p_{T}$ [GeV/c]")
    plot_fs(df_nd, df_fd, 'del_p', "0. FV Cut", "fs0_fv.png")

    # Cosmic cut
    df_nd, df_fd = cosmic_cut(df_nd, df_fd)
    comp.append(plot_int(df_nd, df_fd, 'del_p', "1. Cosmic Rejection", "int1delp_cosmic_rej.png", r"$\delta p$ [GeV/c]", eff_bool=True))
    plot_int(df_nd, df_fd, 'del_Tp', "1. Cosmic Rejection", "int1delpT_cosmic_rej.png", r"$\delta p_{T}$ [GeV/c]")
    plot_fs(df_nd, df_fd, 'del_p', "1. Cosmic Rejection", "fs1_cosmic_rej.png")

    # Two prong cut
    df_nd, df_fd = twoprong_cut(df_nd, df_fd)
    comp.append(plot_int(df_nd, df_fd, 'del_p', "2. Two Prong", "int2delp_two_prong.png", r"$\delta p$ [GeV/c]", eff_bool=True))
    plot_int(df_nd, df_fd, 'del_Tp', "2. Two Prong", "int2delpT_two_prong.png", r"$\delta p_{T}$ [GeV/c]")
    plot_fs(df_nd, df_fd, 'del_p', "2. Two Prong", "fs2_two_prong.png")

    # PID cuts
    df_nd, df_fd = pid_cut(df_nd, df_fd)
    comp.append(plot_int(df_nd, df_fd, 'del_p', "3. 1mu+1p", "int3delp_1mu1p.png", r"$\delta p$ [GeV/c]", eff_bool=True))
    plot_int(df_nd, df_fd, 'del_Tp', "3. 1mu+1p", "int3delpT_1mu1p.png", r"$\delta p_{T}$ [GeV/c]")
    plot_fs(df_nd, df_fd, 'del_p', "3. 1mu+1p", "fs3_1mu1p.png")

    # Stub cut
    df_nd, df_fd = stub_cut(df_nd, df_fd)
    comp.append(plot_int(df_nd, df_fd, 'del_p', "4. Stub", "int4delp_stub.png", r"$\delta p$ [GeV/c]", eff_bool=True))
    plot_int(df_nd, df_fd, 'del_Tp', "4. Stub", "int4delpT_stub.png", r"$\delta p_{T}$ [GeV/c]")
    plot_fs(df_nd, df_fd, 'del_p', "4. Stub", "fs4_stub.png")

    # Composition plots
    comp = np.array(comp)
    cut_labels = ["FV Cut", "Cosmic Cut", "Two Prong", "PID", "Stub"]
    mode_labels = ['QE', 'MEC', 'RES', 'SIS/DIS', 'COH', "other"]
    plot_composition(comp[:, 0, :], cut_labels, mode_labels, 'SBND', 'SBNDComp.png')
    plot_composition(comp[:, 1, :], cut_labels, mode_labels, 'ICARUS', 'ICARUSComp.png')

    return df_nd, df_fd

def main():
    """Main analysis pipeline."""
    dffile_nd = "/exp/sbnd/data/users/nrowe/old_sbnd.df"
    dffile_fd = "/exp/sbnd/app/users/nrowe/gen_ana/sbnana/sbnana/SBNAna/osc-village/gump/makedf_outputs/new_icarus.df"
    df_nd, df_nd_hdr, df_nd_mcnu = load_data(dffile_nd)
    df_fd, df_fd_hdr, df_fd_mcnu = load_data(dffile_fd)
    df_nd['og_sig_ct'] = len(df_nd[is_signal(df_nd, "SBND")])
    df_fd['og_sig_ct'] = len(df_fd[is_signal(df_fd, "ICARUS")])
    des_nd_POT = 1e20
    des_fd_POT = 5e20
    nd_POT, _ = scale_pot(df_nd, df_nd_hdr, des_nd_POT)
    fd_POT, _ = scale_pot(df_fd, df_fd_hdr, des_fd_POT)
    df_nd, df_fd = apply_cuts(df_nd, df_fd, df_nd_mcnu, df_fd_mcnu, nd_POT, fd_POT, dffile_nd, dffile_fd)
    # outputs from here can be fed into PyROOTConverter.py
    save_output(df_nd, df_fd, df_nd_mcnu, df_fd_mcnu, nd_POT, fd_POT, dffile_nd, dffile_fd)

if __name__ == "__main__":
    main()
