import numpy as np


def bump_center_position(var_config, bump_pos=None):
    """Truth-axis position of a center-bin Gaussian bump (matches ``FakeDataWeights``)."""
    bins_arr = np.asarray(var_config.bins, dtype=float)
    centers = np.asarray(var_config.bin_centers, dtype=float)
    axis_mid = 0.5 * (bins_arr[0] + bins_arr[-1])
    idx_center = int(np.argmin(np.abs(centers - axis_mid)))
    lo, hi = bins_arr[idx_center], bins_arr[idx_center + 1]
    if bump_pos is None:
        return float(0.5 * (float(lo) + float(hi)))
    return float(bump_pos)


def format_bump_test_label(var_config, bump_area_bin_fraction, bump_pos=None):
    pos = bump_center_position(var_config, bump_pos=bump_pos)
    pct = float(bump_area_bin_fraction) * 100.0
    return f"Injected Bump at {pos:.1f}, {pct:.0f}% of bin area"


class FakeDataWeights:
    def __init__(self, mc_evt_df, mc_nu_df, var_config):
        self.mc_evt_df = mc_evt_df
        self.mc_nu_df = mc_nu_df
        self.var_config = var_config

    def get_weights(self, test_name, **kwargs):
        weights_fake_data = np.ones(len(self.mc_evt_df))
        weight_fakedata_signal_truth = np.ones(len(self.mc_nu_df[self.mc_nu_df.topo_categ == 1]))
        
        # MEC normalization
        if test_name == "mec_test":
            scale_factor = kwargs.get("scale_factor", 0.5)
            weights_fake_data[self.mc_evt_df.mc.genie_mode == 10] *= scale_factor
            weight_fakedata_signal_truth[self.mc_nu_df[self.mc_nu_df.topo_categ == 1].mc.genie_mode == 10] *= scale_factor

        # QE normalization
        elif test_name == "qe_test":
            scale_factor = kwargs.get("scale_factor", 1.2)
            weights_fake_data[self.mc_evt_df.mc.genie_mode == 0] *= scale_factor
            weight_fakedata_signal_truth[self.mc_nu_df[self.mc_nu_df.topo_categ == 1].mc.genie_mode == 0] *= scale_factor

        # np normalization
        elif test_name == "np_test":
            scale_factor = kwargs.get("scale_factor", 2)
            weights_fake_data[self.mc_evt_df.topo_categ == 2] *= scale_factor
            weight_fakedata_signal_truth[self.mc_nu_df[self.mc_nu_df.topo_categ == 1].topo_categ == 2] *= scale_factor

        # sig normalization
        elif test_name == "sig_test":
            scale_factor = kwargs.get("scale_factor", 1.2)
            weights_fake_data[self.mc_evt_df.topo_categ == 1] *= scale_factor
            weight_fakedata_signal_truth[self.mc_nu_df[self.mc_nu_df.topo_categ == 1].topo_categ == 1] *= scale_factor

        # Q2 tilt
        elif test_name.startswith("q2_test_alpha_"):
            scale_factor = kwargs.get("scale_factor", 0.3)
            Q2 = 2 * self.mc_evt_df.mc.E * self.mc_evt_df.mu.pfp.trk.truth.p.startE * (1 - self.mc_evt_df.mu.pfp.trk.truth.p.dir.z)
            weights_fake_data *= np.ones(len(self.mc_evt_df)) + scale_factor * (Q2 - Q2.mean())/Q2.mean()
            weights_fake_data[np.isnan(weights_fake_data)] = 1
            assert np.isnan(weights_fake_data).sum() == 0
            Q2_nu = self.mc_nu_df[self.mc_nu_df.topo_categ == 1].mc.Q2
            weight_fakedata_signal_truth *= np.ones(len(Q2_nu)) + scale_factor * (Q2_nu - Q2_nu.mean())/Q2_nu.mean()
        
        # cos(theta) scale
        elif test_name.startswith("costh_weight_scale_"):
            scale_factor = kwargs.get("scale_factor", 0.7)
            weights_fake_data[self.mc_evt_df.mu.pfp.trk.truth.p.dir.z > 0.9] *= scale_factor
            weight_fakedata_signal_truth[self.mc_nu_df[self.mc_nu_df.topo_categ == 1].mc.mu.dir.z > 0.9] *= scale_factor

        # Proton P tilt
        elif test_name.startswith("proton_P_tilt_alpha_"):
            scale_factor = kwargs.get("scale_factor", 0.3)
            P_p_evt = self.mc_evt_df.p.pfp.trk.truth.p.totp
            weights_fake_data = np.ones(len(self.mc_evt_df)) + scale_factor * (P_p_evt - P_p_evt.mean())/P_p_evt.mean()
            P_p_nu = self.mc_nu_df[self.mc_nu_df.topo_categ == 1].mc.p.totp
            weight_fakedata_signal_truth = np.ones(len(P_p_nu)) + scale_factor * (P_p_nu - P_p_nu.mean())/P_p_nu.mean()

        # Bump: optional ``bump_area_bin_fraction`` sets total weighted excess
        # ``sum_i w_i (W_i - 1)`` to ``fraction * N_center``, where ``N_center`` is the
        # nominal weighted event count in the **center bin** (bin whose center is closest
        # to the axis midpoint). Bump is centered on that bin by default (``bump_pos``).
        elif test_name.startswith("bump_"):
            bins_arr = np.asarray(self.var_config.bins, dtype=float)
            centers = np.asarray(self.var_config.bin_centers, dtype=float)
            nb = len(centers)

            bump_var_evt = self.mc_evt_df[self.var_config.var_evt_truth_col]
            x_evt = np.asarray(bump_var_evt, dtype=float).ravel()

            if "pot_weight" in self.mc_evt_df.columns:
                w_evt = np.asarray(self.mc_evt_df["pot_weight"], dtype=float).ravel()
            else:
                w_evt = np.ones(len(x_evt), dtype=float)

            bump_area_frac = kwargs.get("bump_area_bin_fraction", None)
            if bump_area_frac is not None:
                axis_mid = 0.5 * (bins_arr[0] + bins_arr[-1])
                idx_center = int(np.argmin(np.abs(centers - axis_mid)))
                lo, hi = bins_arr[idx_center], bins_arr[idx_center + 1]
                bump_pos = kwargs.get("bump_pos", bump_center_position(self.var_config))
                bin_width = float(hi - lo)
                bump_width = kwargs.get("bump_width", None)
                if bump_width is None:
                    bump_width = max(0.15 * bin_width, 1e-12)
                else:
                    bump_width = float(bump_width)

                valid = np.isfinite(x_evt)
                x_safe = np.where(valid, x_evt, bins_arr[0])
                did = np.digitize(x_safe, bins_arr, right=False) - 1
                did = np.clip(did, 0, nb - 1)
                did = np.where(valid & (x_evt >= bins_arr[-1]), nb - 1, did)
                in_center = (did == idx_center) & valid
                N_c = float(np.sum(w_evt * in_center))

                G_evt = np.exp(-0.5 * ((x_evt - bump_pos) / bump_width) ** 2)
                G_evt = np.where(valid, G_evt, 0.0)
                sum_wG = float(np.sum(w_evt * G_evt))

                if sum_wG <= 0.0 or N_c <= 0.0:
                    bump_height_user = float(kwargs.get("bump_height", 0.001))
                else:
                    target = float(bump_area_frac) * N_c
                    bump_height_user = target * nb / (sum_wG * float(len(self.mc_evt_df)))
            else:
                bump_pos = float(kwargs.get("bump_pos", 0.6))
                bump_width = float(kwargs.get("bump_width", 0.0015))
                bump_height_user = float(kwargs.get("bump_height", 0.001))

            bump_pos = float(bump_pos)
            bump_amp = bump_height_user * float(len(self.mc_evt_df)) / float(nb)

            G_evt = np.exp(-0.5 * ((x_evt - bump_pos) / bump_width) ** 2)
            weights_fake_data = np.ones(len(self.mc_evt_df), dtype=float) + bump_amp * G_evt
            weights_fake_data[np.isnan(weights_fake_data)] = 1.0

            sig = self.mc_nu_df.topo_categ == 1
            bump_var_nu = self.mc_nu_df[sig][self.var_config.var_nu_col]
            x_nu = np.asarray(bump_var_nu, dtype=float).ravel()
            G_nu = np.exp(-0.5 * ((x_nu - bump_pos) / bump_width) ** 2)
            G_nu = np.where(np.isfinite(x_nu), G_nu, 0.0)
            weight_fakedata_signal_truth = np.ones(len(x_nu), dtype=float) + bump_amp * G_nu
        
        # TODO
        # # Enhance tail
        # elif test_name.startswith("enhance_tail_alpha_"):
        #     scale_factor = kwargs.get("scale_factor", 0.7)
        #     weights_fake_data = np.ones(len(self.mc_evt_df))
        #     weights_fake_data[self.mc_evt_df.del_p > 0.25] = scale_factor
        #     weights_fake_data[np.isnan(weights_fake_data)] = 1.
        #     weight_fakedata_signal_truth = np.ones(len(self.mc_nu_df[self.mc_nu_df.topo_categ == 1]))
        #     weight_fakedata_signal_truth[self.mc_nu_df[self.mc_nu_df.topo_categ == 1].mc.del_p > 0.25] = scale_factor
        #     weight_fakedata_signal_truth[np.isnan(weight_fakedata_signal_truth)] = 1.

        else:
            raise ValueError(f"Unknown test_name '{test_name}' provided.")

        return weights_fake_data, weight_fakedata_signal_truth