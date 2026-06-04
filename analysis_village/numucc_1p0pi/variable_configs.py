import numpy as np
import inspect

# ===== References
# MicroBooNE tki bins: https://arxiv.org/abs/2301.03700

# Constant axis for the single-bin ``integrated`` measurement (all events in one bin).
# Matches wiremod/sce cov production (``np.full(..., 500.)``), not the ``iscc`` placeholder columns.
INTEGRATED_HIST_DUMMY = 500.0

INTEGRATED_VAR_SAVE_NAME = "integrated"




def is_integrated_var_config(var_config) -> bool:
    return getattr(var_config, "var_save_name", None) == INTEGRATED_VAR_SAVE_NAME


class VariableConfig:
    """
    A configurable class for setting up unfolding variable configurations.
    Choose a configuration using one of the provided class methods,
    or instantiate directly with custom parameters.
    """
    def __init__(
        self,
        var_save_name,
        var_plot_name,
        var_labels,
        bins,
        var_evt_reco_col,
        var_evt_truth_col,
        var_nu_col,
        xsec_label,
        category_syst_var_save_name=None,
    ):
        self.var_save_name = var_save_name
        self.var_plot_name = var_plot_name
        self.var_labels = var_labels
        self.bins = bins
        self.bin_centers = (bins[:-1] + bins[1:]) / 2.
        self.var_evt_reco_col = var_evt_reco_col
        self.var_evt_truth_col = var_evt_truth_col
        self.var_nu_col = var_nu_col
        self.xsec_label = xsec_label
        # Optional: load category-summary syst for a different variable (same bin count).
        self.category_syst_var_save_name = category_syst_var_save_name


    # ==== variables for xsec measurement ====
    # for a integrated single-bin measurement of all events
    @classmethod
    def all_events(cls):
        return cls(
            var_save_name="integrated",
            var_plot_name="All Events",
            var_labels=[r"All Events", 
            r"All Events", ""],
            # use slc.producer as the dummy variable
            bins=np.linspace(0., 1000., 2),
            var_evt_reco_col=('mc', 'iscc', '', '', '', '', ''),
            var_evt_truth_col=('mc', 'iscc', '', '', '', '', ''),
            var_nu_col=('mc', 'iscc', ''),
            xsec_label=r"$\sigma$ $\left[\frac{\mathrm{cm}^2}{\mathrm{Ar}}\right)$"
        )

    @classmethod
    def muon_momentum(cls):
        return cls(
            var_save_name="muon-p",
            var_plot_name="P_\mu",
            var_labels=[r"$\mathrm{P_\mu}$ [GeV/c]", 
            r"$\mathrm{P_\mu^{reco.}}$ [GeV/c]", 
            r"$\mathrm{P_\mu^{true}}$ [GeV/c]"],
            bins=np.array([0.22, 0.27, 0.32, 0.37, 0.42, 0.47, 0.52, 0.57, 0.62, 0.7, 0.8, 0.9, 1.0]),
            var_evt_reco_col=('mu', 'pfp', 'trk', 'P', 'p_muon', '', ''),
            var_evt_truth_col=('mu', 'pfp', 'trk', 'truth', 'p', 'totp', ''),
            var_nu_col=('mc', 'mu', 'totp'),
            xsec_label=r"$\frac{d\sigma}{dP_\mu}$ $\left[\frac{\mathrm{cm}^2}{(\mathrm{GeV}/c)\ \mathrm{Ar}}\right]$"
        )

    @classmethod
    def muon_momentum_mcs(cls):
        return cls(
            var_save_name="muon-p",
            var_plot_name="P_\mu",
            var_labels=[r"$\mathrm{P_\mu}$ [GeV/c]", 
            r"$\mathrm{P_\mu^{reco.}}$ [GeV/c]", 
            r"$\mathrm{P_\mu^{true}}$ [GeV/c]"],
            bins=np.array([0.22, 0.27, 0.32, 0.37, 0.42, 0.47, 0.52, 0.57, 0.62, 0.7, 0.8, 0.9, 1.0]),
            var_evt_reco_col=('mu', 'pfp', 'trk', 'mcsP', 'fwdP_muon', '', ''),
            var_evt_truth_col=('mu', 'pfp', 'trk', 'truth', 'p', 'totp', ''),
            var_nu_col=('mc', 'mu', 'totp'),
            xsec_label=r"$\frac{d\sigma}{dP_\mu}$ $\left[\frac{\mathrm{cm}^2}{(\mathrm{GeV}/c)\ \mathrm{Ar}}\right]$"
        )

    @classmethod
    def muon_momentum_hybrid(cls):
        return cls(
            var_save_name="muon-p",
            var_plot_name="P_\mu",
            var_labels=[r"$\mathrm{P_\mu}$ [GeV/c]", 
            r"$\mathrm{P_\mu^{reco.}}$ [GeV/c]", 
            r"$\mathrm{P_\mu^{true}}$ [GeV/c]"],
            bins=np.array([0.22, 0.27, 0.32, 0.37, 0.42, 0.47, 0.52, 0.57, 0.62, 0.7, 0.8, 0.9, 1.0]),
            var_evt_reco_col=('mu', 'pfp', 'trk', 'hybridP', 'p_muon', '', ''),
            var_evt_truth_col=('mu', 'pfp', 'trk', 'truth', 'p', 'totp', ''),
            var_nu_col=('mc', 'mu', 'totp'),
            xsec_label=r"$\frac{d\sigma}{dP_\mu}$ $\left[\frac{\mathrm{cm}^2}{(\mathrm{GeV}/c)\ \mathrm{Ar}}\right]$"
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
            var_evt_reco_col=('mu', 'pfp', 'trk', 'dir', 'z', '', ''),
            var_evt_truth_col=('mu', 'pfp', 'trk', 'truth', 'p', 'dir', 'z'),
            var_nu_col=('mc', 'mu', 'dir', 'z'),
            xsec_label=r"$\frac{d\sigma}{dcos(\theta_\mu)}$ $\left[\frac{\mathrm{cm}^2}{\mathrm{Ar}}\right]$"
        )

    @classmethod
    def proton_momentum(cls):
        return cls(
            var_save_name="proton-p",
            var_plot_name="P_p",
            var_labels=[r"$\mathrm{P_p}$ [GeV/c]", 
            r"$\mathrm{P_p^{reco.}}$ [GeV/c]", 
            r"$\mathrm{P_p^{true}}$ [GeV/c]"],
            bins=np.array([0.3, 0.37, 0.44, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.]),
            var_evt_reco_col=('p', 'pfp', 'trk', 'P', 'p_proton', '', ''),
            var_evt_truth_col=('p', 'pfp', 'trk', 'truth', 'p', 'totp', ''),
            var_nu_col=('mc', 'p', 'totp'),
            xsec_label=r"$\frac{d\sigma}{dP_p}$ $\left[\frac{\mathrm{cm}^2}{(\mathrm{GeV}/c)\ \mathrm{Ar}}\right]$"
        )

    @classmethod
    def proton_direction(cls):
        return cls(
            var_save_name="proton-dir_z",
            var_plot_name="cos(\\theta_p)",
            var_labels=[r"$\mathrm{cos(\theta_p)}$", 
            r"$\mathrm{cos(\theta_p^{reco.})}$", 
            r"$\mathrm{cos(\theta_p^{true})}$"],
            bins=np.array([-1,-0.7, -0.4, -0.1, 0.1, 0.3, 0.5, 0.7, 0.8, 0.9, 1]),
            var_evt_reco_col=('p', 'pfp', 'trk', 'dir', 'z', '', ''),
            var_evt_truth_col=('p', 'pfp', 'trk', 'truth', 'p', 'dir', 'z'),
            var_nu_col=('mc', 'p', 'dir', 'z'),
            xsec_label=r"$\frac{d\sigma}{dcos(\theta_p)}$ $\left[\frac{\mathrm{cm}^2}{\mathrm{Ar}}\right]$"
        )

    @classmethod
    def tki_del_p(cls):
        return cls(
            var_save_name="tki-del_p",
            var_plot_name="$\\delta p$",
            var_labels=[r"$\mathrm{\delta p}$ [GeV/c]", 
            r"$\mathrm{\delta p^{reco.}}$ [GeV/c]", 
            r"$\mathrm{\delta p^{true}}$ [GeV/c]"],
            bins=np.array([0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.45, 0.55, 0.65, 0.8, 1.0]),
            var_evt_reco_col=('del_p', '', '', '', '', '', ''),
            # var_evt_truth_col=('mc_del_p', '', '', '', '', '', ''),
            var_evt_truth_col=('mc', 'del_p', '', '', '', '', ''),
            var_nu_col=('del_p', '', ''),
            xsec_label=r"$\frac{d\sigma}{d\delta p}$ $\left[\frac{\mathrm{cm}^2}{(\mathrm{GeV}/c)\ \mathrm{Ar}}\right]$"
        )
    
    @classmethod
    def tki_del_Tp(cls):
        return cls(
            var_save_name="tki-del_Tp",
            var_plot_name="$\\delta p_T$",
            var_labels=[r"$\mathrm{\delta p_T}$ [GeV/c]", 
            r"$\mathrm{\delta p_T^{reco.}}$ [GeV/c]", 
            r"$\mathrm{\delta p_T^{true}}$ [GeV/c]"],
            bins=np.array([0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.45, 0.55, 0.65, 0.8, 1.0]),
            var_evt_reco_col=('del_Tp', '', '', '', '', '', ''),
            # var_evt_truth_col=('mc_del_Tp', '', '', '', '', '', ''),
            var_evt_truth_col=('mc', 'del_Tp', '', '', '', '', ''),
            var_nu_col=('del_Tp', '', ''),
            xsec_label=r"$\frac{d\sigma}{d\delta p_T}$ $\left[\frac{\mathrm{cm}^2}{(\mathrm{GeV}/c)\ \mathrm{Ar}}\right]$"
        )

    @classmethod
    def tki_del_Tp_x(cls):
        return cls(
            var_save_name="tki-del_Tp_x",
            var_plot_name="$\\delta p_T^x$",
            var_labels=[r"$\mathrm{\delta p_{T, x}}$ [GeV/c]", 
            r"$\mathrm{\delta p_{T, x}^{reco.}}$ [GeV/c]", 
            r"$\mathrm{\delta p_{T, x}^{true}}$ [GeV/c]"],
            bins=np.array([-0.7, -0.55, -0.45, -0.35, -0.25, -0.15, -0.05, 0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.7]),
            var_evt_reco_col=('del_Tp_x', '', '', '', '', '', ''),
            # var_evt_truth_col=('mc_del_Tp_x', '', '', '', '', '', ''),
            var_evt_truth_col=('mc', 'del_Tp_x', '', '', '', '', ''),
            var_nu_col=('del_Tp_x', '', ''),
            xsec_label=r"$\frac{d\sigma}{d\delta p_T^x}$ $\left[\frac{\mathrm{cm}^2}{(\mathrm{GeV}/c)\ \mathrm{Ar}}\right]$"
        )

    @classmethod
    def tki_del_Tp_y(cls):
        return cls(
            var_save_name="tki-del_Tp_y",
            var_plot_name="$\\delta p_T^y$",
            var_labels=[r"$\mathrm{\delta p_{T, y}}$ [GeV/c]", 
            r"$\mathrm{\delta p_{T, y}^{reco.}}$ [GeV/c]", 
            r"$\mathrm{\delta p_{T, y}^{true}}$ [GeV/c]"],
            bins=np.array([-0.9, -0.7, -0.55, -0.45, -0.35, -0.25, -0.15, -0.05, 0.05, 0.15, 0.25, 0.35, 0.45, 0.55]),
            var_evt_reco_col=('del_Tp_y', '', '', '', '', '', ''),
            # var_evt_truth_col=('mc_del_Tp_y', '', '', '', '', '', ''),
            var_evt_truth_col=('mc', 'del_Tp_y', '', '', '', '', ''),
            var_nu_col=('del_Tp_y', '', ''),
            xsec_label=r"$\frac{d\sigma}{d\delta p_T^y}$ $\left[\frac{\mathrm{cm}^2}{(\mathrm{GeV}/c)\ \mathrm{Ar}}\right]$"
        )
    
    @classmethod
    def tki_del_alpha(cls):
        return cls(
            var_save_name="tki-del_alpha",
            var_plot_name="$\\delta \\alpha_T$",
            var_labels=[r"$\mathrm{\delta \alpha_T}$ [deg]", 
            r"$\mathrm{\delta \alpha_T^{reco.}}$ [deg]", 
            r"$\mathrm{\delta \alpha_T^{true}}$ [deg]"],
            bins=np.array([0,25,50,75,100,120,140,160,180]),
            var_evt_reco_col=('del_alpha', '', '', '', '', '', ''),
            # var_evt_truth_col=('mc_del_alpha', '', '', '', '', '', ''),
            var_evt_truth_col=('mc', 'del_alpha', '', '', '', '', ''),
            var_nu_col=('del_alpha', '', ''),
            xsec_label=r"$\frac{d\sigma}{d\delta \alpha_T}$ $\left[\frac{\mathrm{cm}^2}{\mathrm{deg} \ \mathrm{Ar}}\right]$"
        )
    
    @classmethod
    def tki_del_phi(cls):
        return cls(
            var_save_name="tki-del_phi",
            var_plot_name="$\\delta \\phi_T$",
            var_labels=[r"$\mathrm{\delta \phi_T}$ [deg]", 
            r"$\mathrm{\delta \phi_T^{reco.}}$ [deg]", 
            r"$\mathrm{\delta \phi_T^{true}}$ [deg]"],
            bins=np.array([0,10,20,30,40,55,70,90,110,130,150,180]),
            var_evt_reco_col=('del_phi', '', '', '', '', '', ''),
            # var_evt_truth_col=('mc_del_phi', '', '', '', '', '', ''),
            var_evt_truth_col=('mc', 'del_phi', '', '', '', '', ''),
            var_nu_col=('del_phi', '', ''),
            xsec_label=r"$\frac{d\sigma}{d\delta \phi_T}$ $\left[\frac{\mathrm{cm}^2}{\mathrm{deg} \ \mathrm{Ar}}\right]$"
        )

    @classmethod
    def opening_angle(cls):
        return cls(
            var_save_name="opening_angle",
            var_plot_name="$\\theta_{\\mu, p}$",
            var_labels=[r"$\mathrm{\theta_{\mu, p}}$ [deg]", 
            r"$\mathrm{\theta_{\mu, p}^{reco.}}$ [deg]", 
            r"$\mathrm{\theta_{\mu, p}^{true}}$ [deg]"],
            bins=np.linspace(0., np.pi, 16),
            var_evt_reco_col=('theta_mu_p', '', '', '', '', '', ''),
            var_evt_truth_col=('mc_theta_mu_p', '', '', '', '', '', ''),
            var_nu_col=('theta_mu_p', ''),
            xsec_label=r"$\frac{d\sigma}{d\theta_{\\mu p}}$ $\left[\frac{\mathrm{cm}^2}{\mathrm{deg}}\right]$"
        )

    # ==== additional variables for efficiency inspection ====


    # only use var_nu_col
    @classmethod
    def neutrino_energy(cls):
        return cls(
            var_save_name="E_nu",
            var_plot_name="$E_{\\nu}$",
            var_labels=["$\mathrm{E_{\\nu}}$ [GeV]", 
            "", 
            ""],
            bins=np.linspace(0.2, 1.7, 16),
            var_evt_reco_col=('mc', 'E', '', '', '', '', ''),
            var_evt_truth_col=('mc', 'E', '', '', '', '', ''),
            var_nu_col=('mc', 'E', ''),
            xsec_label=r"$\frac{d\sigma}{dE_{\nu}}$ $\left(\frac{\mathrm{cm}^2}{\mathrm{GeV}}\right)$"
        )

    @classmethod
    def vertex_x(cls):
        return cls(
            var_save_name="vertex_x",
            var_plot_name="Neutrino Vertex X [cm]",
            var_labels=["Neutrino Vertex X [cm]", 
            "Slice Vertex X [cm]", 
            ""],
            bins=np.linspace(-190, 190, 26),
            # var_evt_reco_col=('slc', 'vertex', 'x', '', '', '', ''),
            # var_evt_truth_col=('mc', 'position', 'x', '', '', '', ''),
            var_evt_reco_col=('slc', 'vertex', 'x', '', ''),
            var_evt_truth_col=('mc', 'position', 'x', '', ''),
            var_nu_col=('mc', 'position', 'x'),
            xsec_label=r"$\frac{d\sigma}{d\mathrm{Vertex X}}$ $\left(\frac{\mathrm{cm}^2}{\mathrm{cm}}\right)$"
        )

    @classmethod
    def vertex_y(cls):
        return cls(
            var_save_name="vertex_y",
            var_plot_name="Neutrino Vertex Y [cm]",
            var_labels=["Neutrino Vertex Y [cm]", 
            "Slice Vertex Y [cm]", 
            ""],
            bins=np.linspace(-190, 190, 51),
            var_evt_reco_col=('slc', 'vertex', 'y', '', ''),
            var_evt_truth_col=('mc', 'position', 'y', '', ''),
            var_nu_col=('mc', 'position', 'y'),
            xsec_label=r"$\frac{d\sigma}{d\mathrm{Vertex X}}$ ($\mathrm{cm}^2$ / cm)"
        )

    @classmethod
    def vertex_z(cls):
        return cls(
            var_save_name="vertex_z",
            var_plot_name="Neutrino Vertex Z [cm]",
            var_labels=["Neutrino Vertex Z [cm]", 
            "Slice Vertex Z [cm]", 
            ""],
            bins=np.linspace(0, 450, 51),
            var_evt_reco_col=('slc', 'vertex', 'z', '', ''),
            var_evt_truth_col=('mc', 'position', 'z', '', ''),
            var_nu_col=('mc', 'position', 'z'),
            xsec_label=r"$\frac{d\sigma}{d\mathrm{Vertex X}}$ ($\mathrm{cm}^2$ / cm)"
        )

    @classmethod
    def muon_direction_x(cls):
        return cls(
            var_save_name="muon-dir_x",
            var_plot_name="cos($\\theta_\mu^x$)",
            var_labels=[r"$\mathrm{cos(\theta_\mu^x)}$", 
            r"$\mathrm{cos(\theta_\mu^{x, reco.})}$", 
            r"$\mathrm{cos(\theta_\mu^{x, true})}$"],
            bins=np.linspace(-1, 1, 21),
            var_evt_reco_col=('mu', 'pfp', 'trk', 'dir', 'x', '', ''),
            var_evt_truth_col=('mu', 'pfp', 'trk', 'truth', 'p', 'dir', 'x'),
            var_nu_col=('mc', 'mu', 'dir', 'x'),
            xsec_label=r"$\frac{d\sigma}{dcos(\theta_{\\mu}^x)}$ $\left(\frac{\mathrm{cm}^2}{\mathrm{Ar}}\right)$"
        )

    @classmethod
    def muon_direction_y(cls):
        return cls(
            var_save_name="muon-dir_y",
            var_plot_name="cos($\\theta_\mu^y$)",
            var_labels=[r"$\mathrm{cos(\theta_\mu^y)}$", 
            r"$\mathrm{cos(\theta_\mu^{y, reco.})}$", 
            r"$\mathrm{cos(\theta_\mu^{y, true})}$"],
            bins=np.linspace(-1, 1, 21),
            var_evt_reco_col=('mu', 'pfp', 'trk', 'dir', 'y', '', ''),
            var_evt_truth_col=('mu', 'pfp', 'trk', 'truth', 'p', 'dir', 'y'),
            var_nu_col=('mc', 'mu', 'dir', 'y'),
            xsec_label=r"$\frac{d\sigma}{dcos(\theta_{\\mu}^y)}$ $\left(\frac{\mathrm{cm}^2}{\mathrm{Ar}}\right)$"
        )

    @classmethod
    def muon_direction_phi(cls):
        return cls(
            var_save_name="muon-dir_phi",
            var_plot_name="$\\phi_\mu$",
            var_labels=[r"$\mathrm{\phi_\mu}$ [deg]", 
            r"$\mathrm{\phi_\mu^{reco.}}$ [deg]", 
            r"$\mathrm{\phi_\mu^{true}}$ [deg]"],
            bins=np.linspace(-180, 180, 41),
            var_evt_reco_col=('mu', 'pfp', 'trk', 'phi', '', '', ''),
            var_evt_truth_col=('mu', 'pfp', 'trk', 'truth', 'p', 'phi', ''),
            var_nu_col=('mc', 'mu', 'phi', '', '', '', ''),
            xsec_label=r"$\frac{d\sigma}{d\phi_{\\mu}}$ $\left(\frac{\mathrm{cm}^2}{\mathrm{deg}}\right)$",
            category_syst_var_save_name="muon-dir_x",
        )

    @classmethod
    def muon_end_x(cls):
        return cls(
            var_save_name="muon-end_x",
            var_plot_name="x_\mu",
            var_labels=[r"Muon End X [cm]", 
            r"Muon End X [cm]", 
            r"Muon End X [cm]"],
            bins=np.linspace(-200, 200, 21),
            var_evt_reco_col=('mu', 'pfp', 'trk', 'end', 'x', '', ''),
            var_evt_truth_col=('mu', 'pfp', 'trk', 'truth', 'p', 'end', 'x'),
            var_nu_col=('mc', 'mu', 'end', 'x'),
            xsec_label=r"$\frac{d\sigma}{dx_\mu}$ $\left(\frac{\mathrm{cm}^2}{\mathrm{Ar}}\right)$"
        )

    @classmethod
    def muon_end_y(cls):
        return cls(
            var_save_name="muon-end_y",
            var_plot_name="y_\mu",
            var_labels=[r"Muon End Y [cm]", 
            r"Muon End Y [cm]", 
            r"Muon End Y [cm]"],
            bins=np.linspace(-200, 200, 21),
            var_evt_reco_col=('mu', 'pfp', 'trk', 'end', 'y', '', ''),
            var_evt_truth_col=('mu', 'pfp', 'trk', 'truth', 'p', 'end', 'y'),
            var_nu_col=('mc', 'mu', 'end', 'y'),
            xsec_label=r"$\frac{d\sigma}{dy_\mu}$ $\left(\frac{\mathrm{cm}^2}{\mathrm{Ar}}\right)$"
        )

    @classmethod
    def muon_end_z(cls):
        return cls(
            var_save_name="muon-end_z",
            var_plot_name="z_\mu",
            var_labels=[r"Muon End Z [cm]", 
            r"Muon End Z [cm]", 
            r"Muon End Z [cm]"],
            bins=np.linspace(0, 500, 21),
            var_evt_reco_col=('mu', 'pfp', 'trk', 'end', 'z', '', ''),
            var_evt_truth_col=('mu', 'pfp', 'trk', 'truth', 'p', 'end', 'z'),
            var_nu_col=('mc', 'mu', 'end', 'z'),
            xsec_label=r"$\frac{d\sigma}{dz_\mu}$ $\left(\frac{\mathrm{cm}^2}{\mathrm{Ar}}\right)$"
        )

    @classmethod
    def trk1_end_x(cls):
        return cls(
            var_save_name="trk1-end_x",
            var_plot_name="x_\mu",
            var_labels=[r"trk1 End X [cm]", 
            r"trk1 End X [cm]", 
            r"trk1 End X [cm]"],
            bins=np.linspace(-200, 200, 51),
            var_evt_reco_col=('trk1', 'pfp', 'trk', 'end', 'x', '', ''),
            var_evt_truth_col=('trk1', 'pfp', 'trk', 'truth', 'p', 'end', 'x'),
            var_nu_col=('mc', 'trk1', 'end', 'x'),
            xsec_label=r"$\frac{d\sigma}{dx_trk1}$ $\left(\frac{\mathrm{cm}^2}{\mathrm{Ar}}\right)$"
        )

    @classmethod
    def trk1_end_y(cls):
        return cls(
            var_save_name="trk1-end_y",
            var_plot_name="y_\trk1",
            var_labels=[r"trk1 End Y [cm]", 
            r"trk1 End Y [cm]", 
            r"trk1 End Y [cm]"],
            bins=np.linspace(-200, 200, 51),
            var_evt_reco_col=('trk1', 'pfp', 'trk', 'end', 'y', '', ''),
            var_evt_truth_col=('trk1', 'pfp', 'trk', 'truth', 'p', 'end', 'y'),
            var_nu_col=('mc', 'trk1', 'end', 'y'),
            xsec_label=r"$\frac{d\sigma}{dy_trk1}$ $\left(\frac{\mathrm{cm}^2}{\mathrm{Ar}}\right)$"
        )

    @classmethod
    def trk1_end_z(cls):
        return cls(
            var_save_name="trk1-end_z",
            var_plot_name="z_\trk1",
            var_labels=[r"trk1 End Z [cm]", 
            r"trk1 End Z [cm]", 
            r"trk1 End Z [cm]"],
            bins=np.linspace(0, 500, 51),
            var_evt_reco_col=('trk1', 'pfp', 'trk', 'end', 'z', '', ''),
            var_evt_truth_col=('trk1', 'pfp', 'trk', 'truth', 'p', 'end', 'z'),
            var_nu_col=('mc', 'trk1', 'end', 'z'),
            xsec_label=r"$\frac{d\sigma}{dz_\mu}$ $\left(\frac{\mathrm{cm}^2}{\mathrm{Ar}}\right)$"
        )

    @classmethod
    def muon_momentum_mcs(cls):
        return cls(
            var_save_name="muon-p",
            var_plot_name="P_\mu",
            var_labels=[r"$\mathrm{P_\mu}$ [GeV/c]", 
            r"$\mathrm{P_\mu^{reco.}}$ [GeV/c]", 
            r"$\mathrm{P_\mu^{true}}$ [GeV/c]"],
            bins=np.array([0.22, 0.27, 0.32, 0.37, 0.42, 0.47, 0.52, 0.57, 0.62, 0.7, 0.8, 0.9, 1.0]),
            var_evt_reco_col=('mu', 'pfp', 'trk', 'mcsP', 'fwdP_muon', '', ''),
            var_evt_truth_col=('mu', 'pfp', 'trk', 'truth', 'p', 'totp', ''),
            var_nu_col=('mc', 'mu', 'totp'),
            xsec_label=r"$\frac{d\sigma}{dP_\mu}$ $\left(\frac{\mathrm{cm}^2}{(\mathrm{GeV}/c)\ \mathrm{Ar}}\right)$"
        )

    @classmethod
    def proton_direction_x(cls):
        return cls(
            var_save_name="proton-dir_x",
            var_plot_name="cos($\\theta_p^x$)",
            var_labels=[r"$\mathrm{cos(\theta_p^x)}$", 
            r"$\mathrm{cos(\theta_p^{x, reco.})}$", 
            r"$\mathrm{cos(\theta_p^{x, true})}$"],
            bins=np.linspace(-1, 1, 21),
            var_evt_reco_col=('p', 'pfp', 'trk', 'dir', 'x', '', ''),
            var_evt_truth_col=('p', 'pfp', 'trk', 'truth', 'p', 'dir', 'x'),
            var_nu_col=('mc', 'p', 'dir', 'x'),
            xsec_label=r"$\frac{d\sigma}{dcos(\theta_{\\p}^x)}$ $\left(\frac{\mathrm{cm}^2}{\mathrm{Ar}}\right)$"
        )

    @classmethod
    def proton_direction_y(cls):
        return cls(
            var_save_name="proton-dir_y",
            var_plot_name="cos($\\theta_p^y$)",
            var_labels=[r"$\mathrm{cos(\theta_p^y)}$", 
            r"$\mathrm{cos(\theta_p^{y, reco.})}$", 
            r"$\mathrm{cos(\theta_p^{y, true})}$"],
            bins=np.linspace(-1, 1, 21),
            var_evt_reco_col=('p', 'pfp', 'trk', 'dir', 'y', '', ''),
            var_evt_truth_col=('p', 'pfp', 'trk', 'truth', 'p', 'dir', 'y'),
            var_nu_col=('mc', 'p', 'dir', 'y'),
            xsec_label=r"$\frac{d\sigma}{dcos(\theta_{\\p}^y)}$ $\left(\frac{\mathrm{cm}^2}{\mathrm{Ar}}\right)$"
        )


    @classmethod
    def proton_direction_phi(cls):
        return cls(
            var_save_name="proton-dir_phi",
            var_plot_name="$\\phi_p$",
            var_labels=[r"$\mathrm{\phi_p}$ [deg]", 
            r"$\mathrm{\phi_p^{reco.}}$ [deg]", 
            r"$\mathrm{\phi_p^{true}}$ [deg]"],
            bins=np.linspace(-180, 180, 21),
            var_evt_reco_col=('p', 'pfp', 'trk', 'phi', '', '', ''),
            var_evt_truth_col=('p', 'pfp', 'trk', 'truth', 'p', 'phi', ''),
            var_nu_col=('mc', 'p', 'phi', ''),
            xsec_label=r"$\frac{d\sigma}{d\phi_{\\p}}$ $\left(\frac{\mathrm{cm}^2}{\mathrm{deg}}\right)$"
        )

    @classmethod
    def trk1_len(cls):
        return cls(
            var_save_name="trk1_len",
            var_plot_name="Length [cm]",
            var_labels=[r"trk1 Length [cm]", 
            r"trk1 Length [cm]", 
            r"trk1 Length [cm]"],
            bins=np.linspace(0, 500, 51),
            var_evt_reco_col=('trk1', 'pfp', 'trk', 'len', '', '', ''),
            var_evt_truth_col=('trk1', 'pfp', 'trk', 'truth', 'p', 'len', ''),
            var_nu_col=('mc', 'trk1', 'len', ''),
            xsec_label=r"$\frac{d\sigma}{d\mathrm{Length_{\\trk1}}}$ $\left(\frac{\mathrm{cm}^2}{\mathrm{Ar}}\right)$"
        )

    @classmethod
    def trk1_direction_phi(cls):
        return cls(
            var_save_name="trk1-dir_phi",
            var_plot_name="$\\phi$",
            var_labels=[r"$\mathrm{\phi}$ [deg]", 
            r"$\mathrm{\phi^{reco.}}$ [deg]", 
            r"$\mathrm{\phi^{true}}$ [deg]"],
            bins=np.linspace(-180, 180, 41),
            var_evt_reco_col=('trk1', 'pfp', 'trk', 'phi', '', '', ''),
            var_evt_truth_col=('trk1', 'pfp', 'trk', 'truth', 'p', 'phi', ''),
            var_nu_col=('mc', 'trk1', 'phi', '', '', '', ''),
            xsec_label=r"$\frac{d\sigma}{d\phi_{\\mu}}$ $\left(\frac{\mathrm{cm}^2}{\mathrm{deg}}\right)$"
        )

    @classmethod
    def trk1_direction_x(cls):
        return cls(
            var_save_name="trk1-dir_x",
            var_plot_name="x",
            var_labels=[r"x [cm]", 
            r"x [cm]", 
            r"x [cm]"],
            bins=np.linspace(-1, 1, 21),
            var_evt_reco_col=('trk1', 'pfp', 'trk', 'dir', 'x', '', ''),
            var_evt_truth_col=('trk1', 'pfp', 'trk', 'truth', 'p', 'dir', 'x'),
            var_nu_col=('mc', 'trk1', 'dir', 'x'),
            xsec_label=r"$\frac{d\sigma}{dx_{\\trk1}}$ $\left(\frac{\mathrm{cm}^2}{\mathrm{Ar}}\right)$"
        )

    @classmethod
    def trk1_direction_y(cls):
        return cls(
            var_save_name="trk1-dir_y",
            var_plot_name="y",
            var_labels=[r"y [cm]", 
            r"y [cm]", 
            r"y [cm]"],
            bins=np.linspace(-1, 1, 21),
            var_evt_reco_col=('trk1', 'pfp', 'trk', 'dir', 'y', '', ''),
            var_evt_truth_col=('trk1', 'pfp', 'trk', 'truth', 'p', 'dir', 'y'),
            var_nu_col=('mc', 'trk1', 'dir', 'y'),
            xsec_label=r"$\frac{d\sigma}{dy_{\\trk1}}$ $\left(\frac{\mathrm{cm}^2}{\mathrm{Ar}}\right)$"
        )

    @classmethod
    def costh_trk1_trk2(cls):
        return cls(
            var_save_name="costh_trk1_trk2",
            var_plot_name=r"cos($\theta_{trk1,trk2}$)",
            var_labels=[r"cos($\theta_{trk1,trk2}$)", 
            r"cos($\theta_{trk1,trk2}^{reco.}$)", 
            r"cos($\theta_{trk1,trk2}^{true}$)"],
            bins=np.linspace(-1, 1, 41),
            var_evt_reco_col=('tracks_cos_theta', '', '', '', '', '', ''),
            var_evt_truth_col=('tracks_cos_theta', '', '', '', '', '', ''),
            var_nu_col=('tracks_cos_theta', '', '', '', '', '', ''),
            xsec_label=r"$\frac{d\sigma}{dcos(\theta_{trk1,trk2})}$ $\left(\frac{\mathrm{cm}^2}{\mathrm{Ar}}\right)$"
        )

    @classmethod
    def trk2_direction_phi(cls):
        return cls(
            var_save_name="trk2-dir_phi",
            var_plot_name="$\\phi$",
            var_labels=[r"$\mathrm{\phi}$ [deg]", 
            r"$\mathrm{\phi^{reco.}}$ [deg]", 
            r"$\mathrm{\phi^{true}}$ [deg]"],
            bins=np.linspace(-180, 180, 31),
            var_evt_reco_col=('trk2', 'pfp', 'trk', 'phi', '', '', ''),
            var_evt_truth_col=('trk2', 'pfp', 'trk', 'truth', 'p', 'phi', ''),
            var_nu_col=('mc', 'trk2', 'phi', '', '', '', ''),
            xsec_label=r"$\frac{d\sigma}{d\phi_{\\mu}}$ $\left(\frac{\mathrm{cm}^2}{\mathrm{deg}}\right)$"
        )

    @classmethod
    def trk1_end_x(cls):
        return cls(
            var_save_name="trk1-end_x",
            var_plot_name="x_\mu",
            var_labels=[r"Track 1 End X [cm]", 
            r"Track 1 End X [cm]", 
            r"Track 1 End X [cm]"],
            bins=np.linspace(-200, 200, 21),
            var_evt_reco_col=('trk1', 'pfp', 'trk', 'end', 'x', '', ''),
            var_evt_truth_col=('trk1', 'pfp', 'trk', 'truth', 'p', 'end', 'x'),
            var_nu_col=('mc', 'trk1', 'end', 'x'),
            xsec_label=r"$\frac{d\sigma}{dx_{\\trk1}}$ $\left(\frac{\mathrm{cm}^2}{\mathrm{Ar}}\right)$"
        )

    @classmethod
    def trk2_end_x(cls):
        return cls(
            var_save_name="trk2-end_x",
            var_plot_name="x_\mu",
            var_labels=[r"Track 2 End X [cm]", 
            r"Track 2 End X [cm]", 
            r"Track 2 End X [cm]"],
            bins=np.linspace(-200, 200, 21),
            var_evt_reco_col=('trk2', 'pfp', 'trk', 'end', 'x', '', ''),
            var_evt_truth_col=('trk2', 'pfp', 'trk', 'truth', 'p', 'end', 'x'),
            var_nu_col=('mc', 'trk2', 'end', 'x'),
            xsec_label=r"$\frac{d\sigma}{dx_{\\trk2}}$ $\left(\frac{\mathrm{cm}^2}{\mathrm{Ar}}\right)$"
        )

    # ==== additional variables for event selection ====
    @classmethod
    def nu_score(cls):
        return cls(
            var_save_name="nu_score",
            var_plot_name="$\\nu_{\\mathrm{score}}$",
            var_labels=[r"Neutrino-like Score", 
            "", 
            ""],
            bins=np.linspace(0, 1, 101),
            var_evt_reco_col=('slc', 'nu_score', '', '', ''),
            var_evt_truth_col=('slc', 'nu_score', '', '', ''),
            var_nu_col=('slc', 'nu_score', ''),
            xsec_label=r""
        )

    @classmethod
    def n_trks(cls):
        return cls(
            var_save_name="n_trks",
            var_plot_name="# of Tracks",
            var_labels=[r"", 
            "", 
            ""],
            bins=np.linspace(1, 7, 7),
            var_evt_reco_col=('n_trks', '', '', '', '', '', ''),
            var_evt_truth_col=('n_trks', '', '', '', '', '', ''),
            var_nu_col=('n_trks', '', ''),
            xsec_label=r""
        )

    @classmethod
    def track_score(cls):
        return cls(
            var_save_name="track_score",
            var_plot_name="PFP Track-like Score",
            var_labels=[r"PFP Track-like Score", 
            "", 
            ""],
            bins=np.linspace(0.2, 0.9, 71),
            var_evt_reco_col=('trk1', 'pfp', 'trackScore', '', '', '', ''),
            var_evt_truth_col=('trk1', 'pfp', 'trackScore', '', '', '', ''),
            var_nu_col=('', '', ''),
            xsec_label=r""
        )

    @classmethod
    def track_end_x(cls):
        return cls(
            var_save_name="track_end_x",
            var_plot_name="Track End X",
            var_labels=[r"Track End X [cm]", 
            "", 
            ""],
            bins=np.linspace(-10, 10, 41),
            var_evt_reco_col=('pfp', 'trk', 'end', 'x', '', ''),
            var_evt_truth_col=('pfp', 'trk', 'truth', 'end', 'x', ''),
            var_nu_col=('trk', 'end', 'x'),
            xsec_label=r"$\frac{d\sigma}{d\mathrm{Track End X}}$ $\left(\frac{\mathrm{cm}^2}{\mathrm{cm}}\right)$"
        )

    @classmethod
    def vtx_dist(cls):
        return cls(
            var_save_name="vtx_dist",
            var_plot_name="Slice Vertex - Track Start Distance [cm]",
            var_labels=[r"Slice Vertex - Track Start Distance [cm]", 
            "", 
            ""],
            bins=np.linspace(0, 6, 41),
            var_evt_reco_col=('pfp', 'pfochar', 'vtxdist', '', '', ''),
            var_evt_truth_col=('pfp', 'pfochar', 'vtxdist', '', '', ''),
            var_nu_col=('trk', 'vtxdist', ''),
            xsec_label=r""
        )

    @classmethod
    def trk_len(cls):
        return cls(
            var_save_name="trk_len",
            var_plot_name="Track Length [cm]",
            var_labels=[r"$\mathrm{Track \, \, Length}$ (cm)", 
            r"$\mathrm{Track Length^{reco.}}$ [cm]", 
            r"$\mathrm{Track Length^{true}}$ [cm]"],
            bins=np.linspace(0, 100, 26),
            var_evt_reco_col=('pfp', 'trk', 'len', '', '', ''),
            var_evt_truth_col=('pfp', 'trk', 'truth', 'len', '', ''),
            var_nu_col=('trk', 'len', ''),
            xsec_label=r""
        )

    @classmethod
    def mcs_range_diff(cls):
        return cls(
            var_save_name="mcs_range_diff",
            var_plot_name="MCS Range Difference",
            var_labels=[r"$\mathrm{(MCS - Range) \, / \, Range}$", 
            r"$\mathrm{(MCS - Range)^{reco.}}$", 
            r"$\mathrm{(MCS - Range)^{true}}$"],
            bins=np.linspace(-0.6, 0.4, 51),
            var_evt_reco_col=('pfp', 'trk', 'mcs_range_diff', '', '', ''),
            var_evt_truth_col=('pfp', 'trk', 'truth', 'mcs_range_diff', '', '', ''),
            var_nu_col=('trk', 'mcs_range_diff', '', ''),
            xsec_label=r""
        )

    @classmethod
    def chi2_mu(cls):
        return cls(
            var_save_name="chi2_mu",
            var_plot_name="$\\chi^2_{\\mu}$",
            var_labels=[r"$\mathrm{\chi^{2}_{\mu}}$",
            r"$\mathrm{\chi^{2}_{\mu,\mathrm{reco.}}}$",
            r"$\mathrm{\chi^{2}_{\mu,\mathrm{true}}}$"],
            bins=np.linspace(0, 55, 31),
            var_evt_reco_col=('pfp', 'trk', 'chi2pid', 'avg', 'chi2_muon',  ''),
            var_evt_truth_col=('', '', '', '', '', ''),
            var_nu_col=('', '', ''),
            xsec_label=r"$\frac{d\sigma}{d\chi^2_{\\mu}}$ ($\mathrm{cm}^2$)"
        )

    @classmethod
    def chi2_proton(cls):
        return cls(
            var_save_name="chi2_p",
            var_plot_name="$\\chi^2_{p}$",
            var_labels=[r"$\mathrm{\chi^{2}_{p}}$",
            r"$\mathrm{\chi^{2}_{p,\mathrm{reco.}}}$",
            r"$\mathrm{\chi^{2}_{p,\mathrm{true}}}$"],
            bins=np.linspace(0, 300, 41),
            var_evt_reco_col=('mu', 'pfp', 'trk', 'chi2pid', 'I2', 'chi2_proton',  ''),
            var_evt_truth_col=('', '', '', '', '', ''),
            var_nu_col=('', '', ''),
            xsec_label=r"$\frac{d\sigma}{d\chi^2_{\\p}}$ ($\mathrm{cm}^2$)"
        )

    @classmethod
    def trk1_chi2_mu(cls):
        return cls(
            var_save_name="trk1_chi2_mu",
            var_plot_name="$\\chi^2_{\\mu}$",
            var_labels=[r"$\mathrm{\chi^{2}_{\mu}}$",
            r"$\mathrm{\chi^{2}_{\mu,\mathrm{reco.}}}$"],
            bins=np.linspace(0, 60, 61),
            var_evt_reco_col=('trk1', 'pfp', 'trk', 'chi2pid', 'I2', 'chi2_muon',  ''),
            var_evt_truth_col=('trk1', 'pfp', 'trk', 'truth', 'chi2pid', 'I2', 'chi2_muon'),
            var_nu_col=('trk1', 'pfp', 'trk', 'truth', 'chi2pid', 'I2', 'chi2_muon'),
            xsec_label=r"$\frac{d\sigma}{d\chi^2_{\\mu}}$ ($\mathrm{cm}^2$)"
        )

    @classmethod
    def trk1_chi2_proton(cls):
        return cls(
            var_save_name="trk1_chi2_p",
            var_plot_name="$\\chi^2_{p}$",
            var_labels=[r"$\mathrm{\chi^{2}_{p}}$",
            r"$\mathrm{\chi^{2}_{p,\mathrm{reco.}}}$",
            r"$\mathrm{\chi^{2}_{p,\mathrm{true}}}$"],
            bins=np.linspace(0, 300, 61),
            var_evt_reco_col=('trk1', 'pfp', 'trk', 'chi2pid', 'I2', 'chi2_proton',  ''),
            var_evt_truth_col=('trk1', 'pfp', 'trk', 'truth', 'chi2pid', 'I2', 'chi2_proton'),
            var_nu_col=('trk1', 'pfp', 'trk', 'truth', 'chi2pid', 'I2', 'chi2_proton'),
            xsec_label=r"$\frac{d\sigma}{d\chi^2_{\\p}}$ ($\mathrm{cm}^2$)"
        )

    @classmethod
    def trk2_chi2_mu(cls):
        return cls(
            var_save_name="trk2_chi2_mu",
            var_plot_name="$\\chi^2_{\\mu}$",
            var_labels=[r"$\mathrm{\chi^{2}_{\mu}}$",
            r"$\mathrm{\chi^{2}_{\mu,\mathrm{reco.}}}$"],
            bins=np.linspace(0, 60, 61),
            var_evt_reco_col=('trk2', 'pfp', 'trk', 'chi2pid', 'I2', 'chi2_muon',  ''),
            var_evt_truth_col=('trk2', 'pfp', 'trk', 'truth', 'chi2pid', 'I2', 'chi2_muon'),
            var_nu_col=('trk2', 'pfp', 'trk', 'truth', 'chi2pid', 'I2', 'chi2_muon'),
            xsec_label=r"$\frac{d\sigma}{d\chi^2_{\\mu}}$ ($\mathrm{cm}^2$)"
        )

    @classmethod
    def trk2_chi2_proton(cls):
        return cls(
            var_save_name="trk2_chi2_p",
            var_plot_name="$\\chi^2_{p}$",
            var_labels=[r"$\mathrm{\chi^{2}_{p}}$",
            r"$\mathrm{\chi^{2}_{p,\mathrm{reco.}}}$",
            r"$\mathrm{\chi^{2}_{p,\mathrm{true}}}$"],
            bins=np.linspace(0, 300, 61),
            var_evt_reco_col=('trk2', 'pfp', 'trk', 'chi2pid', 'I2', 'chi2_proton',  ''),
            var_evt_truth_col=('trk2', 'pfp', 'trk', 'truth', 'chi2pid', 'I2', 'chi2_proton'),
            var_nu_col=('trk2', 'pfp', 'trk', 'truth', 'chi2pid', 'I2', 'chi2_proton'),
            xsec_label=r"$\frac{d\sigma}{d\chi^2_{\\p}}$ ($\mathrm{cm}^2$)"
        )

    @classmethod
    def prim_direction_phi(cls):
        return cls(
            var_save_name="trk_dir_phi",
            var_plot_name="$\\phi_{\\mathrm{trk}}$",
            var_labels=[r"$\mathrm{\phi}$ [deg]", 
            r"$\mathrm{\phi^{reco.}}$ [deg]", 
            r"$\mathrm{\phi^{true}}$ [deg]"],
            bins=np.linspace(-180, 180, 21),
            # var_evt_reco_col=('prim', 'pfp', 'trk', 'phi', '', ''),
            var_evt_reco_col=('prim_trk_phi', '', '', '', ''),
            var_evt_truth_col=('prim', 'pfp', 'trk', 'truth', 'phi', ''),
            var_nu_col=('prim', 'trk', 'phi', '', '', '', ''),
            xsec_label=r"$\frac{d\sigma}{d\phi_{\\mathrm{trk}}}$ $\left(\frac{\mathrm{cm}^2}{\mathrm{deg}}\right)$"
        )

    @classmethod
    def prim_direction_y(cls):
        return cls(
            var_save_name="trk_dir_y",
            var_plot_name="costh y",
            var_labels=[r"costh y", 
            r"costh y reco.", 
            r"costh y true"],
            bins=np.linspace(-1, 1, 21),
            var_evt_reco_col=('prim_trk_dir_y', '', '', '', ''),
            var_evt_truth_col=('prim_trk_dir_y', '', '', '', ''),
            var_nu_col=('prim_trk_dir_y', '', ''),
            xsec_label=r"$\frac{d\sigma}{dy_{\\mathrm{trk}}}$ $\left(\frac{\mathrm{cm}^2}{\mathrm{cm}}\right)$"
        )

    @classmethod
    def prim_start_x(cls):
        return cls(
            var_save_name="trk_start_x",
            var_plot_name="Start X (cm)",
            var_labels=[r"Track Start X [cm]", 
            "", 
            ""],
            bins=np.linspace(-200, 200, 101),
            var_evt_reco_col=('prim_trk_start_x', '', '', '', ''),
            var_evt_truth_col=('prim_trk_start_x', '', '', '', ''),
            var_nu_col=('prim_trk_start_x', '', ''),
            xsec_label=r"$\frac{d\sigma}{d\mathrm{Track Start X}}$ $\left(\frac{\mathrm{cm}^2}{\mathrm{cm}}\right)$"
        )

    @classmethod
    def prim_end_x(cls):
        return cls(
            var_save_name="trk_end_x",
            var_plot_name="End X (cm)",
            var_labels=[r"Track End X [cm]", 
            "", 
            ""],
            bins=np.linspace(-50, 50, 101),
            var_evt_reco_col=('prim_trk_end_x', '', '', '', ''),
            var_evt_truth_col=('prim_trk_end_x', '', '', '', ''),
            var_nu_col=('prim_trk_end_x', '', ''),
            xsec_label=r"$\frac{d\sigma}{d\mathrm{Track End X}}$ $\left(\frac{\mathrm{cm}^2}{\mathrm{cm}}\right)$"
        )

    @classmethod
    def prim_P_frac_diff(cls):
        return cls(
            var_save_name="prim_P_frac_diff",
            var_plot_name="(mcs P - range P) / range P",
            var_labels=[r"(mcs P - range P) / range P", 
            "", 
            ""],
            bins=np.linspace(-1, 0.5, 31),
            var_evt_reco_col=('prim_trk_P_frac_diff', '', '', '', ''),
            var_evt_truth_col=('prim_trk_P_frac_diff', '', '', '', ''),
            var_nu_col=('prim_trk_P_frac_diff', '', ''),
            xsec_label=r"$\frac{d\sigma}{d\mathrm{Track End X}}$ $\left(\frac{\mathrm{cm}^2}{\mathrm{cm}}\right)$"
        )

    @classmethod
    def prim_chi2pid_I0_muon(cls):
        return cls(
            var_save_name="prim_chi2pid_I0_muon",
            var_plot_name="$\\chi^2_{\\mu, I0}$",
            var_labels=[r"$\mathrm{\chi^2_{\mu, I0}}$", 
            "", 
            ""],
            bins=np.linspace(0, 60, 61),
            var_evt_reco_col=('prim_trk_chi2pid_I0_muon', '', '', '', ''),
            var_evt_truth_col=('prim_trk_chi2pid_I0_muon', '', '', '', ''),
            var_nu_col=('prim_trk_chi2pid_I0_muon', '', ''),
            xsec_label=r"$\frac{d\sigma}{d\chi^2_{\\mu, I0}}$ ($\mathrm{cm}^2$)"
        )

    @classmethod
    def prim_chi2pid_I0_proton(cls):
        return cls(
            var_save_name="prim_chi2pid_I0_proton",
            var_plot_name="$\\chi^2_{\\p, I0}$",
            var_labels=[r"$\mathrm{\chi^2_{p, I0}}$", 
            "", 
            ""],
            bins=np.linspace(0, 300, 61),
            var_evt_reco_col=('prim_trk_chi2pid_I0_proton', '', '', '', ''),
            var_evt_truth_col=('prim_trk_chi2pid_I0_proton', '', '', '', ''),
            var_nu_col=('prim_trk_chi2pid_I0_proton', '', ''),
            xsec_label=r"$\frac{d\sigma}{d\chi^2_{\\p, I0}}$ ($\mathrm{cm}^2$)"
        )

    @classmethod
    def prim_chi2pid_I1_muon(cls):
        return cls(
            var_save_name="prim_chi2pid_I1_muon",
            var_plot_name="$\\chi^2_{\\mu, I1}$",
            var_labels=[r"$\mathrm{\chi^2_{\mu, I1}}$", 
            "", 
            ""],
            bins=np.linspace(0, 60, 61),
            var_evt_reco_col=('prim_trk_chi2pid_I1_muon', '', '', '', ''),
            var_evt_truth_col=('prim_trk_chi2pid_I1_muon', '', '', '', ''),
            var_nu_col=('prim_trk_chi2pid_I1_muon', '', ''),
            xsec_label=r"$\frac{d\sigma}{d\chi^2_{\\mu, I1}}$ ($\mathrm{cm}^2$)"
        )

    @classmethod
    def prim_chi2pid_I1_proton(cls):
        return cls(
            var_save_name="prim_chi2pid_I1_proton",
            var_plot_name="$\\chi^2_{\\p, I1}$",
            var_labels=[r"$\mathrm{\chi^2_{p, I1}}$", 
            "", 
            ""],
            bins=np.linspace(0, 300, 61),
            var_evt_reco_col=('prim_trk_chi2pid_I1_proton', '', '', '', ''),
            var_evt_truth_col=('prim_trk_chi2pid_I1_proton', '', '', '', ''),
            var_nu_col=('prim_trk_chi2pid_I1_proton', '', ''),
            xsec_label=r"$\frac{d\sigma}{d\chi^2_{\\p, I1}}$ ($\mathrm{cm}^2$)" 
        )

    @classmethod
    def prim_chi2pid_I2_muon(cls):
        return cls(
            var_save_name="prim_chi2pid_I2_muon",
            var_plot_name="$\\chi^2_{\\mu, I2}$",
            var_labels=[r"$\mathrm{\chi^2_{\mu, I2}}$", 
            "", 
            ""],
            bins=np.linspace(0, 60, 61),
            var_evt_reco_col=('prim_trk_chi2pid_I2_muon', '', '', '', ''),
            var_evt_truth_col=('prim_trk_chi2pid_I2_muon', '', '', '', ''),
            var_nu_col=('prim_trk_chi2pid_I2_muon', '', ''),
            xsec_label=r"$\frac{d\sigma}{d\chi^2_{\\mu, I2}}$ ($\mathrm{cm}^2$)"
        )

    @classmethod
    def prim_chi2pid_I2_proton(cls):
        return cls(
            var_save_name="prim_chi2pid_I2_proton",
            var_plot_name="$\\chi^2_{\\p, I2}$",
            var_labels=[r"$\mathrm{\chi^2_{p, I2}}$", 
            "", 
            ""],
            bins=np.linspace(0, 300, 61),
            var_evt_reco_col=('prim_trk_chi2pid_I2_proton', '', '', '', ''),
            var_evt_truth_col=('prim_trk_chi2pid_I2_proton', '', '', '', ''),
            var_nu_col=('prim_trk_chi2pid_I2_proton', '', ''),
            xsec_label=r"$\frac{d\sigma}{d\chi^2_{\\p, I2}}$ ($\mathrm{cm}^2$)"
        )

    @classmethod
    def trk_direction_phi(cls):
        return cls(
            var_save_name="trk_dir_phi",
            var_plot_name="$\\phi_{\\mathrm{trk}}$",
            var_labels=[r"$\mathrm{\phi}$ [deg]", 
            r"$\mathrm{\phi^{reco.}}$ [deg]", 
            r"$\mathrm{\phi^{true}}$ [deg]"],
            bins=np.linspace(-180, 180, 21),
            var_evt_reco_col=('pfp', 'trk', 'phi', '', '', ''),
            var_evt_truth_col=('pfp', 'trk', 'truth', 'phi', '', ''),
            var_nu_col=('trk', 'phi', '', '', '', ''),
            xsec_label=r"$\frac{d\sigma}{d\phi_{\\mathrm{trk}}}$ $\left(\frac{\mathrm{cm}^2}{\mathrm{deg}}\right)$"
        )

    @classmethod
    def list_all_configs(cls, print_summary=True):
        members = inspect.getmembers(cls, predicate=inspect.isroutine)
        config_methods = [name for name, func in members
                          if getattr(func, "__self__", None) == cls
                          and not name.startswith("_")
                          and name not in ("list_all_configs")]
        print("Available VariableConfig options:")
        for name in config_methods:
            config = getattr(cls, name)()
            if print_summary:
                print(f"  {name:20}")
                print(f"    var_save_name   : {config.var_save_name}")
                print(f"    var_plot_name   : {config.var_plot_name}")
                print(f"    var_labels      : {config.var_labels}")
                print(f"    bins            : {np.array2string(config.bins, separator=', ')}")
                print(f"    bin_centers     : {np.array2string(config.bin_centers, separator=', ')}")
                print(f"    var_evt_reco_col: {config.var_evt_reco_col}")
                print(f"    var_evt_truth_col: {config.var_evt_truth_col}")
                print(f"    var_nu_col      : {config.var_nu_col}")
                print(f"    xsec_label      : {config.xsec_label}")
                print()
            else:
                print(f"  {name:20}")
        return config_methods


var_configs_measurement = [
                VariableConfig.all_events(),
                VariableConfig.muon_momentum(),
                VariableConfig.muon_direction(),
                VariableConfig.proton_momentum(),
                VariableConfig.proton_direction(),
                VariableConfig.tki_del_alpha(),
                VariableConfig.tki_del_phi(),
                VariableConfig.tki_del_Tp(),
                VariableConfig.tki_del_p(),
                VariableConfig.tki_del_Tp_x(),
                VariableConfig.tki_del_Tp_y(),
                ]

var_configs_extra_finalstate = [
                VariableConfig.muon_direction_x(),
                VariableConfig.muon_direction_y(),
                VariableConfig.muon_direction_phi(),
                VariableConfig.proton_direction_x(),
                VariableConfig.proton_direction_y(),
                VariableConfig.proton_direction_phi(),
                VariableConfig.opening_angle(),
                ]

var_configs_extra_slc = [
                VariableConfig.vertex_x(),
                VariableConfig.vertex_y(),
                VariableConfig.vertex_z(),
                ]