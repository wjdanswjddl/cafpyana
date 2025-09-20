import numpy as np

# tki bins: https://arxiv.org/abs/2301.03700
# mu, p kinematics bins: 
class VariableConfig:
    """
    A configurable class for setting up unfolding variable configurations.
    Choose a configuration using one of the provided class methods,
    or instantiate directly with custom parameters.
    """
    def __init__(self, var_save_name, var_plot_name, var_labels, bins, var_evt_reco_col, var_evt_truth_col, var_nu_col, xsec_label):
        self.var_save_name = var_save_name
        self.var_plot_name = var_plot_name
        self.var_labels = var_labels
        self.bins = bins
        self.bin_centers = (bins[:-1] + bins[1:]) / 2.
        self.var_evt_reco_col = var_evt_reco_col
        self.var_evt_truth_col = var_evt_truth_col
        self.var_nu_col = var_nu_col
        self.xsec_label = xsec_label


    @classmethod
    def muon_momentum(cls):
        return cls(
            var_save_name="muon-p",
            var_plot_name="P_\mu",
            var_labels=[r"$\mathrm{P_\mu}$ [GeV/c]", 
            r"$\mathrm{P_\mu^{reco.}}$ [GeV/c]", 
            r"$\mathrm{P_\mu^{true}}$ [GeV/c]"],
            bins=np.array([0.22, 0.27, 0.32, 0.37, 0.42, 0.47, 0.52, 0.57, 0.62, 0.67, 0.72, 0.77, 0.82, 0.91, 1.0]),
            var_evt_reco_col=('mu', 'pfp', 'trk', 'P', 'p_muon', '', '', ''),
            var_evt_truth_col=('mu', 'pfp', 'trk', 'truth', 'p', 'totp', '', ''),
            var_nu_col=('mu', 'totp', ''),
            xsec_label=r"$\frac{d\sigma}{dP_\mu}$ [$cm^2 / GeV/c$]"
        )

    @classmethod
    def muon_direction(cls):
        return cls(
            var_save_name="muon-dir_z",
            var_plot_name="cos(\\theta_\mu)",
            var_labels=[r"$\mathrm{cos(\theta_\mu)}$", 
            r"$\mathrm{cos(\theta_\mu^{reco.})}$", 
            r"$\mathrm{cos(\theta_\mu^{true})}$"],
            bins=np.array([-1, -0.75, -0.6, -0.45, -0.3, -0.15, 0, 0.15, 0.3, 0.45, 0.6, 0.7, 0.8, 0.9, 1]),
            var_evt_reco_col=('mu', 'pfp', 'trk', 'dir', 'z', '', '', ''),
            var_evt_truth_col=('mu', 'pfp', 'trk', 'truth', 'p', 'dir', 'z', ''),
            var_nu_col=('mu', 'dir', 'z'),
            xsec_label=r"$\frac{d\sigma}{dcos(\theta_\mu)}$ [$cm^2$]"
        )

    @classmethod
    def proton_momentum(cls):
        return cls(
            var_save_name="proton-p",
            var_plot_name="P_p",
            var_labels=[r"$\mathrm{P_p}$ [GeV/c]", 
            r"$\mathrm{P_p^{reco.}}$ [GeV/c]", 
            r"$\mathrm{P_p^{true}}$ [GeV/c]"],
            bins=np.array([0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.]),
            var_evt_reco_col=('p', 'pfp', 'trk', 'P', 'p_proton', '', '', ''),
            var_evt_truth_col=('p', 'pfp', 'trk', 'truth', 'p', 'totp', '', ''),
            var_nu_col=('p', 'totp', ''),
            xsec_label=r"$\frac{d\sigma}{dP_p}$ [$cm^2 / GeV/c$]"
        )

    @classmethod
    def proton_direction(cls):
        return cls(
            var_save_name="proton-dir_z",
            var_plot_name="cos(\\theta_p)",
            var_labels=[r"$\mathrm{cos(\theta_p)}$", 
            r"$\mathrm{cos(\theta_p^{reco.})}$", 
            r"$\mathrm{cos(\theta_p^{true})}$"],
            bins=np.array([-1, -0.75, -0.6, -0.45, -0.3, -0.15, 0, 0.15, 0.3, 0.45, 0.6, 0.7, 0.8, 0.9, 1]),
            var_evt_reco_col=('p', 'pfp', 'trk', 'dir', 'z', '', '', ''),
            var_evt_truth_col=('p', 'pfp', 'trk', 'truth', 'p', 'dir', 'z', ''),
            var_nu_col=('p', 'dir', 'z'),
            xsec_label=r"$\frac{d\sigma}{dcos(\theta_p)}$ [$cm^2$]"
        )

    @classmethod
    def tki_del_p(cls):
        return cls(
            var_save_name="tki-del_p",
            var_plot_name="$\\delta p$",
            var_labels=[r"$\mathrm{\delta p}$ [GeV/c]", 
            r"$\mathrm{\delta p^{reco.}}$ [GeV/c]", 
            r"$\mathrm{\delta p^{true}}$ [GeV/c]"],
            bins=np.array([0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.40, 0.45, 0.5, 0.6, 0.7, 0.8, 1.]),
            var_evt_reco_col=('del_p', '', '', '', '', '', '', ''),
            var_evt_truth_col=('mc_del_p', '', '', '', '', '', '', ''),
            var_nu_col=('del_p', '', ''),
            xsec_label=r"$\frac{d\sigma}{d\delta p}$ [$cm^2 / GeV/c$]"
        )
    
    @classmethod
    def tki_del_Tp(cls):
        return cls(
            var_save_name="tki-del_Tp",
            var_plot_name="$\\delta p_T$",
            var_labels=[r"$\mathrm{\delta p_T}$ [GeV/c]", 
            r"$\mathrm{\delta p_T^{reco.}}$ [GeV/c]", 
            r"$\mathrm{\delta p_T^{true}}$ [GeV/c]"],
            bins=np.array([0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.47, 0.55, 0.65, 0.75, 0.9]),
            var_evt_reco_col=('del_Tp', '', '', '', '', '', '', ''),
            var_evt_truth_col=('mc_del_Tp', '', '', '', '', '', '', ''),
            var_nu_col=('del_Tp', '', ''),
            xsec_label=r"$\frac{d\sigma}{d\delta p_T}$ [$cm^2 / GeV/c$]"
        )

    @classmethod
    def tki_del_Tp_x(cls):
        return cls(
            var_save_name="tki-del_Tp_x",
            var_plot_name="$\\delta p_T^x$",
            var_labels=[r"$\mathrm{\delta p_{T, x}}$ [GeV/c]", 
            r"$\mathrm{\delta p_{T, x}^{reco.}}$ [GeV/c]", 
            r"$\mathrm{\delta p_{T, x}^{true}}$ [GeV/c]"],
            bins=np.array([-0.55, -0.45, -0.35, -0.25, -0.15, -0.05, 0.05, 0.15, 0.25, 0.35, 0.45, 0.55]),
            var_evt_reco_col=('del_Tp_x', '', '', '', '', '', '', ''),
            var_evt_truth_col=('mc_del_Tp_x', '', '', '', '', '', '', ''),
            var_nu_col=('del_Tp_x', '', ''),
            xsec_label=r"$\frac{d\sigma}{d\delta p_T^x}$ [$cm^2 / GeV/c$]"
        )

    @classmethod
    def tki_del_Tp_y(cls):
        return cls(
            var_save_name="tki-del_Tp_y",
            var_plot_name="$\\delta p_T^y$",
            var_labels=[r"$\mathrm{\delta p_{T, y}}$ [GeV/c]", 
            r"$\mathrm{\delta p_{T, y}^{reco.}}$ [GeV/c]", 
            r"$\mathrm{\delta p_{T, y}^{true}}$ [GeV/c]"],
            bins=np.array([-0.55, -0.45, -0.35, -0.25, -0.15, -0.05, 0.05, 0.15, 0.25, 0.35, 0.45, 0.55]),
            var_evt_reco_col=('del_Tp_y', '', '', '', '', '', '', ''),
            var_evt_truth_col=('mc_del_Tp_y', '', '', '', '', '', '', ''),
            var_nu_col=('del_Tp_y', '', ''),
            xsec_label=r"$\frac{d\sigma}{d\delta p_T^y}$ [$cm^2 / GeV/c$]"
        )
    
    @classmethod
    def tki_del_alpha(cls):
        return cls(
            var_save_name="tki-del_alpha",
            var_plot_name="$\\delta \\alpha_T$",
            var_labels=[r"$\mathrm{\delta \alpha_T}$ [deg]", 
            r"$\mathrm{\delta \alpha_T^{reco.}}$ [deg]", 
            r"$\mathrm{\delta \alpha_T^{true}}$ [deg]"],
            bins=np.array([0,22,44,66,88,110,145,180]),
            var_evt_reco_col=('del_alpha', '', '', '', '', '', '', ''),
            var_evt_truth_col=('mc_del_alpha', '', '', '', '', '', '', ''),
            var_nu_col=('del_alpha', '', ''),
            xsec_label=r"$\frac{d\sigma}{d\delta \alpha_T}$ [$cm^2/deg$]"
        )
    
    @classmethod
    def tki_del_phi(cls):
        return cls(
            var_save_name="tki-del_phi",
            var_plot_name="$\\delta \\phi_T$",
            var_labels=[r"$\mathrm{\delta \phi_T}$ [deg]", 
            r"$\mathrm{\delta \phi_T^{reco.}}$ [deg]", 
            r"$\mathrm{\delta \phi_T^{true}}$ [deg]"],
            bins=np.array([0,12.5,25,37.5,50,60,75,90,106,126,145,162,180]),
            var_evt_reco_col=('del_phi', '', '', '', '', '', '', ''),
            var_evt_truth_col=('mc_del_phi', '', '', '', '', '', '', ''),
            var_nu_col=('del_phi', '', ''),
            xsec_label=r"$\frac{d\sigma}{d\delta \phi_T}$ [$cm^2/deg$]"
        )
