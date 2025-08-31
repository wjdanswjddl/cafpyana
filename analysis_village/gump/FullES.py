import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from memory_profiler import profile

from plot_tools import *
# Add the head direcoty to sys.path
workspace_root = os.getcwd()  
sys.path.insert(0, workspace_root + "/../../")

# import this repo's classes
import pyanalib.pandas_helpers as ph
from makedf.util import *

import kinematics

def load_data(file, nfiles=1):
    """Load event, header, and mcnu data from HDF file."""
    for s in range(nfiles):
        print("df index:"+str(s))
        df_evt = pd.read_hdf(file, "evt_"+str(s))
        df_hdr = pd.read_hdf(file, "hdr_"+str(s))
        df_mcnu = pd.read_hdf(file, "mcnu_"+str(s))
        df_stub = pd.read_hdf(file, "stub_"+str(s))

        matchdf = df_evt.copy()
        matchdf.columns = pd.MultiIndex.from_tuples([(col, '') for col in matchdf.columns])
        df_evt = ph.multicol_merge(matchdf.reset_index(), df_mcnu.reset_index(),
                                    left_on=[("__ntuple", ""), ("entry", ""), ("tmatch_idx", "")],
                                    right_on=[("__ntuple", ""), ("entry", ""), ("rec.mc.nu..index", "")],
                                    how="left") ## -- save all sllices

        cols_to_drop = []
        for c in df_evt.columns:
            if 'GENIE' in c[0] or 'Flux' in c[0]:
                cols_to_drop.append(c)

        df_evt.drop(columns=cols_to_drop, inplace=True)
        del df_mcnu

        if s == 0:
            res_df_evt = df_evt
            res_df_hdr = df_hdr
            res_df_stub = df_stub
        else:
            res_df_evt = pd.concat([res_df_evt, df_evt])
            res_df_hdr = pd.concat([res_df_hdr, df_hdr])
            res_df_stub = pd.concat([res_df_stub, df_stub])

        del df_evt
        del df_hdr
        del df_stub

    return res_df_evt, res_df_hdr, res_df_stub

def scale_pot(df, df_hdr, desired_pot):
    """Scale DataFrame by desired POT."""
    pot = sum(df_hdr.pot.tolist())
    print(f"POT: {pot}\nScaling to: {desired_pot}")
    scale = desired_pot / pot
    df['glob_scale'] = scale
    return pot, scale

def apply_cuts(df_nd, df_nd_stub, df_fd, df_fd_stub, nd_POT, fd_POT, dffile_nd, dffile_fd):
    comp_sbnd = []
    comp_icarus = []
    cut_labels = ["FV Cut", "Cosmic Cut", "Two Prong", "PID", "Stub"]
    mode_labels = ['QE', 'MEC', 'RES', 'SIS/DIS', 'COH', "other"]

    # FV Cut
    nd_vtx = pd.DataFrame({'x': df_nd.slc_vtx_x,
                           'y': df_nd.slc_vtx_y,
                           'z': df_nd.slc_vtx_z})

    fd_vtx = pd.DataFrame({'x': df_fd.slc_vtx_x,
                           'y': df_fd.slc_vtx_y,
                           'z': df_fd.slc_vtx_z})

    df_nd = df_nd[fv_cut(nd_vtx, "SBND")]
    df_fd = df_fd[fv_cut(fd_vtx, "ICARUS")]

    comp_sbnd.append(plot_int(df_nd, 'del_p', "FV Cut", "int0delp_fv_sbnd.png", r"$\delta p$ [GeV/c]", mode_labels, 'SBND', eff_bool=True))
    comp_icarus.append(plot_int(df_fd, 'del_p', "FV Cut", "int0delp_fv_icarus.png", r"$\delta p$ [GeV/c]", mode_labels, 'ICARUS', eff_bool=True))
    plot_int(df_nd, 'del_Tp', "FV Cut", "int0delpT_fv_sbnd.png", r"$\delta p_{T}$ [GeV/c]", mode_labels, 'SBND')
    plot_int(df_fd, 'del_Tp', "FV Cut", "int0delpT_fv_icarus.png", r"$\delta p_{T}$ [GeV/c]", mode_labels, 'ICARUS')
    # plot_fs(df_nd, 'del_p', "FV Cut", "fs0_fv_sbnd.png", 'SBND')
    # plot_fs(df_fd, 'del_p', "FV Cut", "fs0_fv_icarus.png", 'ICARUS')

    plot_nuscore_cut([
        df_nd.nu_score[cosmic_cut(df_nd)], df_nd.nu_score[np.invert(cosmic_cut(df_nd))]
    ], [0.5], [
        df_fd.nu_score[cosmic_cut(df_fd)], df_fd.nu_score[np.invert(cosmic_cut(df_fd))]
    ], [0.5], "NuScore", "NuScore.png", [0, 1],
        ['SBND Cosmic', 'SBND Neutrino', 'ICARUS Cosmic', 'ICARUS Neutrino'], 'right', 'select as neutrino')

    # Cosmic cut
    df_nd = df_nd[cosmic_cut(df_nd)]
    df_fd = df_fd[cosmic_cut(df_fd)]
    comp_sbnd.append(plot_int(df_nd, 'del_p', "Cosmic Rejection Cut", "int1delp_cosmic_rej_sbnd.png", r"$\delta p$ [GeV/c]", mode_labels, 'SBND', eff_bool=True))
    comp_icarus.append(plot_int(df_fd, 'del_p', "Cosmic Rejection Cut", "int1delp_cosmic_rej_icarus.png", r"$\delta p$ [GeV/c]", mode_labels, 'ICARUS', eff_bool=True))
    plot_int(df_nd, 'del_Tp', "Cosmic Rejection Cut", "int1delpT_cosmic_rej_sbnd.png", r"$\delta p_{T}$ [GeV/c]", mode_labels, 'SBND')
    plot_int(df_fd, 'del_Tp', "Cosmic Rejection Cut", "int1delpT_cosmic_rej_icarus.png", r"$\delta p_{T}$ [GeV/c]", mode_labels, 'ICARUS')
    # plot_fs(df_nd, 'del_p', "Cosmic Rejection Cut", "fs1_cosmic_rej_sbnd.png", 'SBND')
    # plot_fs(df_fd, 'del_p', "Cosmic Rejection Cut", "fs1_cosmic_rej_icarus.png", 'ICARUS')

    # Two prong cut
    df_nd = df_nd[twoprong_cut(df_nd)]
    df_fd = df_fd[twoprong_cut(df_fd)]
    comp_sbnd.append(plot_int(df_nd, 'del_p', "Two Prong Cut", "int2delp_two_prong_sbnd.png", r"$\delta p$ [GeV/c]", mode_labels, 'SBND', eff_bool=True))
    comp_icarus.append(plot_int(df_fd, 'del_p', "Two Prong Cut", "int2delp_two_prong_icarus.png", r"$\delta p$ [GeV/c]", mode_labels, 'ICARUS', eff_bool=True))
    plot_int(df_nd, 'del_Tp', "Two Prong Cut", "int2delpT_two_prong_sbnd.png", r"$\delta p_{T}$ [GeV/c]", mode_labels, 'SBND')
    plot_int(df_fd, 'del_Tp', "Two Prong Cut", "int2delpT_two_prong_icarus.png", r"$\delta p_{T}$ [GeV/c]", mode_labels, 'ICARUS')
    # plot_fs(df_nd, 'del_p', "Two Prong Cut", "fs2_two_prong_sbnd.png", 'SBND')
    # plot_fs(df_fd, 'del_p', "Two Prong Cut", "fs2_two_prong_icarus.png", 'ICARUS')

    plot_PID_cut([[df_nd.mu_chi2_of_prot_cand, df_nd.mu_chi2_of_mu_cand],
                  [df_fd.mu_chi2_of_prot_cand, df_fd.mu_chi2_of_mu_cand]],
                  [25, 25], "MuScore", "MuScore.png", [0, 65],
                  ['SBND, Protons', 'SBND, Muons', 'ICARUS, Protons', 'ICARUS, Muons'],
                  'left', 'select as muon')

    plot_PID_cut([[df_nd.prot_chi2_of_prot_cand, df_nd.prot_chi2_of_mu_cand],
                  [df_fd.prot_chi2_of_prot_cand, df_fd.prot_chi2_of_mu_cand]],
                  [100, 100], "ProtonScore", "ProtonScore.png", [0, 200],
                  ['SBND, Protons', 'SBND, Muons', 'ICARUS, Protons', 'ICARUS, Muons'],
                  'right', 'select as muon')

    # PID cuts
    bmu = np.linspace(0, 70, 30)
    plt.hist(df_nd.mu_chi2_of_mu_cand, bins=bmu, histtype='step', label='mu cand')
    plt.hist(df_nd.mu_chi2_of_prot_cand, bins=bmu, histtype='step', label='prot cand')
    plt.legend()
    plt.savefig("mu_score.png")
    plt.clf()
    bprot = np.linspace(0, 200, 30)
    plt.hist(df_nd.prot_chi2_of_mu_cand, bins=bprot, histtype='step', label='mu cand')
    plt.hist(df_nd.prot_chi2_of_prot_cand, bins=bprot, histtype='step', label='prot cand')
    plt.legend()
    plt.savefig("prot_score.png")
    plt.clf()

    df_nd = df_nd[pid_cut(df_nd.mu_chi2_of_mu_cand, df_nd.mu_chi2_of_prot_cand, 
                          df_nd.prot_chi2_of_mu_cand, df_nd.prot_chi2_of_prot_cand, df_nd.mu_len)]
    df_fd = df_fd[pid_cut(df_fd.mu_chi2_of_mu_cand, df_fd.mu_chi2_of_prot_cand, 
                          df_fd.prot_chi2_of_mu_cand, df_fd.prot_chi2_of_prot_cand, df_fd.mu_len)]
    comp_sbnd.append(plot_int(df_nd, 'del_p', "1mu+1p Cut", "int3delp_1mu1p_sbnd.png", r"$\delta p$ [GeV/c]", mode_labels, 'SBND', eff_bool=True))
    comp_icarus.append(plot_int(df_fd, 'del_p', "1mu+1p Cut", "int3delp_1mu1p_icarus.png", r"$\delta p$ [GeV/c]", mode_labels, 'ICARUS', eff_bool=True))
    plot_int(df_nd, 'del_Tp', "1mu+1p Cut", "int3delpT_1mu1p_sbnd.png", r"$\delta p_{T}$ [GeV/c]", mode_labels, 'SBND')
    plot_int(df_fd, 'del_Tp', "1mu+1p Cut", "int3delpT_1mu1p_icarus.png", r"$\delta p_{T}$ [GeV/c]", mode_labels, 'ICARUS')
    # plot_fs(df_nd, 'del_p', "1mu+1p Cut", "fs3_1mu1p_sbnd.png", 'SBND')
    # plot_fs(df_fd, 'del_p', "1mu+1p Cut", "fs3_1mu1p_icarus.png", 'ICARUS')


    # Stub cut
    df_nd_stub = df_nd_stub.reset_index('rec.slc.reco.stub..index', drop=True)
    sel = df_nd.set_index(['__ntuple', 'entry', 'rec.slc..index'], drop=True).index
    sel = sel.intersection(df_nd_stub.index)
    df_nd_stub = df_nd_stub.loc[sel]
    plot_stub_2d(df_nd_stub['length'], df_nd_stub['dqdx'], "SBND.png", "SBND Stub Cut")
    # plot_stub_2d(df_nd[(df_nd.is_sig != True)], "FPSBND.png", "False Positive SBND Stub Cut")
    plot_stub_2d(df_fd_stub['length'], df_fd_stub['dqdx'], "ICARUS.png", "ICARUS Stub Cut")
    # plot_stub_2d(df_fd[(df_fd.is_sig != True)], "FPICARUS.png", "False Positive ICARUS Stub Cut")

    df_nd = df_nd[stub_cut(df_nd)]
    df_fd = df_fd[stub_cut(df_fd)]
    comp_sbnd.append(plot_int(df_nd, 'del_p', "Stub Cut", "int4delp_stub_sbnd.png", r"$\delta p$ [GeV/c]", mode_labels, 'SBND', eff_bool=True))
    comp_icarus.append(plot_int(df_fd, 'del_p', "Stub Cut", "int4delp_stub_icarus.png", r"$\delta p$ [GeV/c]", mode_labels, 'ICARUS', eff_bool=True))
    plot_int(df_nd, 'del_Tp', "Stub Cut", "int4delpT_stub_sbnd.png", r"$\delta p_{T}$ [GeV/c]", mode_labels, 'SBND')
    plot_int(df_fd, 'del_Tp', "Stub Cut", "int4delpT_stub_icarus.png", r"$\delta p_{T}$ [GeV/c]", mode_labels, 'ICARUS')
    # plot_fs(df_nd, 'del_p', "Stub Cut", "fs4_stub_sbnd.png", 'SBND')
    # plot_fs(df_fd, 'del_p', "Stub Cut", "fs4_stub_icarus.png", 'ICARUS')

    # Composition plots
    comp_sbnd = np.array(comp_sbnd)
    comp_icarus = np.array(comp_icarus)
    plot_composition(comp_sbnd, cut_labels, mode_labels, 'SBND', 'SBNDComp.png')
    plot_composition(comp_icarus, cut_labels, mode_labels, 'ICARUS', 'ICARUSComp.png')

    return df_nd, df_fd

def main():
    """Main analysis pipeline."""
    dffile_nd = "/home/nathanielerowe/SBN/cafpyana_gump/sbnd_no_cuts.df"
    dffile_fd = "/home/nathanielerowe/SBN/cafpyana_gump/icarus_no_cuts.df"
    df_nd, df_nd_hdr, df_nd_stub = load_data(dffile_nd)
    df_fd, df_fd_hdr, df_fd_stub = load_data(dffile_fd)

    print('data loaded!')
    df_nd['og_sig_ct'] = len(df_nd.is_sig[df_nd.is_sig == True])
    df_fd['og_sig_ct'] = len(df_fd.is_sig[df_fd.is_sig == True])

    des_nd_POT = 1e20
    des_fd_POT = 5e20
    nd_POT = scale_pot(df_nd, df_nd_hdr, des_nd_POT)
    fd_POT = scale_pot(df_fd, df_fd_hdr, des_fd_POT)
    print("Starting to apply cuts!")
    df_nd, df_fd = apply_cuts(df_nd, df_nd_stub, df_fd, df_fd_stub, nd_POT, fd_POT, dffile_nd, dffile_fd)

if __name__ == "__main__":
    main()
