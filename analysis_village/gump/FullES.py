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

def load_data(file, start=0, stop=10000):
    """Load event, header, and mcnu data from HDF file."""
    df_evt = pd.read_hdf(file, "evt_0")
    df_hdr = pd.read_hdf(file, "hdr_0")
    return df_evt, df_hdr

def scale_pot(df, df_hdr, desired_pot):
    """Scale DataFrame by desired POT."""
    pot = sum(df_hdr.pot.tolist())
    print(f"POT: {pot}\nScaling to: {desired_pot}")
    scale = desired_pot / pot
    df['glob_scale'] = scale
    return pot, scale

def apply_cuts(df_nd, df_fd, nd_POT, fd_POT, dffile_nd, dffile_fd):
    print("Applying cuts...")
    comp_sbnd = []
    comp_icarus = []
    cut_labels = ["FV Cut", "Cosmic Cut", "Two Prong", "PID", "Stub"]
    mode_labels = ['QE', 'MEC', 'RES', 'SIS/DIS', 'COH', "other"]

    # FV Cut
    print("FV Cut")
    df_nd, df_fd = fv_cut(df_nd, df_fd)
    comp_sbnd.append(plot_int(df_nd, 'del_p', "FV Cut", "int0delp_fv_sbnd.png", r"$\delta p$ [GeV/c]", mode_labels, 'SBND', eff_bool=True))
    comp_icarus.append(plot_int(df_fd, 'del_p', "FV Cut", "int0delp_fv_icarus.png", r"$\delta p$ [GeV/c]", mode_labels, 'ICARUS', eff_bool=True))
    plot_int(df_nd, 'del_Tp', "FV Cut", "int0delpT_fv_sbnd.png", r"$\delta p_{T}$ [GeV/c]", mode_labels, 'SBND')
    plot_int(df_fd, 'del_Tp', "FV Cut", "int0delpT_fv_icarus.png", r"$\delta p_{T}$ [GeV/c]", mode_labels, 'ICARUS')
    plot_fs(df_nd, 'del_p', "FV Cut", "fs0_fv_sbnd.png", 'SBND')
    plot_fs(df_fd, 'del_p', "FV Cut", "fs0_fv_icarus.png", 'ICARUS')

    # Cosmic cut
    print("Cosmic Cut")
    df_nd, df_fd = cosmic_cut(df_nd, df_fd)
    comp_sbnd.append(plot_int(df_nd, 'del_p', "Cosmic Rejection Cut", "int1delp_cosmic_rej_sbnd.png", r"$\delta p$ [GeV/c]", mode_labels, 'SBND', eff_bool=True))
    comp_icarus.append(plot_int(df_fd, 'del_p', "Cosmic Rejection Cut", "int1delp_cosmic_rej_icarus.png", r"$\delta p$ [GeV/c]", mode_labels, 'ICARUS', eff_bool=True))
    plot_int(df_nd, 'del_Tp', "Cosmic Rejection Cut", "int1delpT_cosmic_rej_sbnd.png", r"$\delta p_{T}$ [GeV/c]", mode_labels, 'SBND')
    plot_int(df_fd, 'del_Tp', "Cosmic Rejection Cut", "int1delpT_cosmic_rej_icarus.png", r"$\delta p_{T}$ [GeV/c]", mode_labels, 'ICARUS')
    plot_fs(df_nd, 'del_p', "Cosmic Rejection Cut", "fs1_cosmic_rej_sbnd.png", 'SBND')
    plot_fs(df_fd, 'del_p', "Cosmic Rejection Cut", "fs1_cosmic_rej_icarus.png", 'ICARUS')

    # Two prong cut
    print("Two prong Cut")
    df_nd, df_fd = twoprong_cut(df_nd, df_fd)
    comp_sbnd.append(plot_int(df_nd, 'del_p', "Two Prong Cut", "int2delp_two_prong_sbnd.png", r"$\delta p$ [GeV/c]", mode_labels, 'SBND', eff_bool=True))
    comp_icarus.append(plot_int(df_fd, 'del_p', "Two Prong Cut", "int2delp_two_prong_icarus.png", r"$\delta p$ [GeV/c]", mode_labels, 'ICARUS', eff_bool=True))
    plot_int(df_nd, 'del_Tp', "Two Prong Cut", "int2delpT_two_prong_sbnd.png", r"$\delta p_{T}$ [GeV/c]", mode_labels, 'SBND')
    plot_int(df_fd, 'del_Tp', "Two Prong Cut", "int2delpT_two_prong_icarus.png", r"$\delta p_{T}$ [GeV/c]", mode_labels, 'ICARUS')
    plot_fs(df_nd, 'del_p', "Two Prong Cut", "fs2_two_prong_sbnd.png", 'SBND')
    plot_fs(df_fd, 'del_p', "Two Prong Cut", "fs2_two_prong_icarus.png", 'ICARUS')

    # PID cuts
    print("PID Cut")
    df_nd, df_fd = pid_cut(df_nd, df_fd)
    comp_sbnd.append(plot_int(df_nd, 'del_p', "1mu+1p Cut", "int3delp_1mu1p_sbnd.png", r"$\delta p$ [GeV/c]", mode_labels, 'SBND', eff_bool=True))
    comp_icarus.append(plot_int(df_fd, 'del_p', "1mu+1p Cut", "int3delp_1mu1p_icarus.png", r"$\delta p$ [GeV/c]", mode_labels, 'ICARUS', eff_bool=True))
    plot_int(df_nd, 'del_Tp', "1mu+1p Cut", "int3delpT_1mu1p_sbnd.png", r"$\delta p_{T}$ [GeV/c]", mode_labels, 'SBND')
    plot_int(df_fd, 'del_Tp', "1mu+1p Cut", "int3delpT_1mu1p_icarus.png", r"$\delta p_{T}$ [GeV/c]", mode_labels, 'ICARUS')
    plot_fs(df_nd, 'del_p', "1mu+1p Cut", "fs3_1mu1p_sbnd.png", 'SBND')
    plot_fs(df_fd, 'del_p', "1mu+1p Cut", "fs3_1mu1p_icarus.png", 'ICARUS')

    # Stub cut
    print("Stub Cut")
    df_nd, df_fd = stub_cut(df_nd, df_fd)
    comp_sbnd.append(plot_int(df_nd, 'del_p', "Stub Cut", "int4delp_stub_sbnd.png", r"$\delta p$ [GeV/c]", mode_labels, 'SBND', eff_bool=True))
    comp_icarus.append(plot_int(df_fd, 'del_p', "Stub Cut", "int4delp_stub_icarus.png", r"$\delta p$ [GeV/c]", mode_labels, 'ICARUS', eff_bool=True))
    plot_int(df_nd, 'del_Tp', "Stub Cut", "int4delpT_stub_sbnd.png", r"$\delta p_{T}$ [GeV/c]", mode_labels, 'SBND')
    plot_int(df_fd, 'del_Tp', "Stub Cut", "int4delpT_stub_icarus.png", r"$\delta p_{T}$ [GeV/c]", mode_labels, 'ICARUS')
    plot_fs(df_nd, 'del_p', "Stub Cut", "fs4_stub_sbnd.png", 'SBND')
    plot_fs(df_fd, 'del_p', "Stub Cut", "fs4_stub_icarus.png", 'ICARUS')

    # Composition plots
    print("Final Composition")
    comp_sbnd = np.array(comp_sbnd)
    comp_icarus = np.array(comp_icarus)
    plot_composition(comp_sbnd, cut_labels, mode_labels, 'SBND', 'SBNDComp.png')
    plot_composition(comp_icarus, cut_labels, mode_labels, 'ICARUS', 'ICARUSComp.png')

    return df_nd, df_fd

def main():
    """Main analysis pipeline."""
    dffile_nd = "no_cuts.df"
    dffile_fd = "no_cuts.df"
    df_nd, df_nd_hdr = load_data(dffile_nd)
    df_fd, df_fd_hdr = load_data(dffile_fd)
    df_nd['og_sig_ct'] = 1 # len(df_nd[is_signal(df_nd, "SBND")])
    df_fd['og_sig_ct'] = 1 # len(df_fd[is_signal(df_fd, "ICARUS")])
    des_nd_POT = 1e20
    des_fd_POT = 5e20
    nd_POT = 1 #, _ = scale_pot(df_nd, df_nd_hdr, des_nd_POT)
    fd_POT = 1 #, _ = scale_pot(df_fd, df_fd_hdr, des_fd_POT)
    df_nd, df_fd = apply_cuts(df_nd, df_fd, nd_POT, fd_POT, dffile_nd, dffile_fd)

if __name__ == "__main__":
    main()
