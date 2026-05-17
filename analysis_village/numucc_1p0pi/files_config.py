"""Monolithic HDF loaders for notebooks and cosmics blocks in multisim aggregate.

Chunked drivers take input paths from :mod:`analysis_village.numucc_1p0pi.dataset_locations`.
"""

from pyanalib.split_df_helpers import *
from analysis_village.numucc_1p0pi.utils import *


# ==== save configs ====
save_fig_base_dir = "/exp/sbnd/data/users/munjung/plots/numucc1p0pi"


# ==== Spring Gen 1 samples, besides detvar samples ====
file_dir = "/exp/sbnd/data/users/munjung/xsec/2025Spring_v10_06_00_09"

n_max_concat = 999


def get_ana_dfs(option="", syst_tag="", systs_mc_df_tag="", systs_chunk_tags=None):
    """Load analysis HDF bundles.

    Parameters
    ----------
    systs_mc_df_tag : str
        When ``option == "systs"``: suffix on each chunk file (e.g. ``""`` for
        nominal mup-weighted dfs, ``"-sel_all-wgts"`` for loose selection + weights).
    systs_chunk_tags : sequence of str or None
        When ``option == "systs"`` and not None: chunk tags to concatenate;
        default ``generate_tags("bl")[-5:]`` (final-selection-style chunks).
        Use e.g. ``generate_tags("ah")[1:]`` with ``systs_mc_df_tag="-sel_all-wgts"``
        to match :func:`get_ana_dfs` ``event_selection`` MC paths.
    """

    if option == "systs": # flux, g4, mcstat systs on evt dfs (mup-final or sel_all, etc.)
        chunk_tags = systs_chunk_tags if systs_chunk_tags is not None else generate_tags("bl")[-5:]
        ret_dfs = load_and_concat_mc_dfs(
            file_dir=file_dir,
            sub_dir="MC",
            sample_dir="BNB_cosmics",
            df_tag=systs_mc_df_tag,
            chunk_tags=chunk_tags,
            keys2load=['hdr', 'evt'],
            n_max_concat=n_max_concat
        )
        return ret_dfs

    elif option == "genie_systs": # need mcnu
        ret_dfs = load_and_concat_mc_dfs(
            file_dir=file_dir,
            sub_dir="MC",
            sample_dir="BNB_cosmics/genie_wgts-zexp",
            df_tag="",
            # df_tag="_geniewgts_CCQE",
            chunk_tags=generate_tags("ad"),
            # chunk_tags=["ac"],
            keys2load=['hdr', 'mcnu', 'evt'],
            n_max_concat=n_max_concat
        )
        return ret_dfs

    elif option == "cosmics_systs": # intime MC and offbeam data
        data_dfs = load_and_concat_mc_dfs(
            file_dir=file_dir,
            chunk_tags=generate_tags("ad"),
            df_tag="",
            keys2load=['evt', 'hdr'],
            n_max_concat=999,
            sub_dir="data",
            sample_dir="OffBeam"
        )
        data_evt_df = data_dfs['evt']
        data_hdr_df = data_dfs['hdr']
        data_gates = data_hdr_df[data_hdr_df['first_in_subrun'] == 1]['noffbeambnb'].sum()
        print("intime cosmics data gates: {:.2e}".format(data_gates))

        mc_dfs = load_and_concat_mc_dfs(
            file_dir=file_dir,
            chunk_tags=generate_tags("au"),
            df_tag="",
            keys2load=['evt', 'hdr'],
            n_max_concat=999,
            sub_dir="MC",
            sample_dir="intime"
        )
        mc_evt_df = mc_dfs['evt']
        mc_hdr_df = mc_dfs['hdr']
        mc_gates = mc_hdr_df[mc_hdr_df['first_in_subrun'] == 1]['ngenevt'].sum()
        print("intime cosmics MC gates: {:.2e}".format(mc_gates))

        scale = data_gates/mc_gates
        print("offbeam scale: {:.2f}".format(scale))

        mc_evt_df["pot_scale"] = scale
        data_evt_df["pot_scale"] = 1.

        return {"mc": mc_evt_df, "data": data_evt_df}

    elif option == "fake_data_test":
        concat_dfs = load_and_concat_mc_dfs(
            # file_dir=file_dir,
            file_dir="/pnfs/sbnd/scratch/users/munjung/xsec/2025Spring_v10_06_00_09",
            chunk_tags=generate_tags("ad"),
            # df_tag="_sel_mup-geniewgts",
            df_tag="",
            keys2load=['hdr', 'mcnu', 'evt'],
            n_max_concat=n_max_concat,
            sub_dir="MC",
            sample_dir="BNB_cosmics/genie_wgts-RES"
        )
        mc_hdr_df = concat_dfs['hdr']
        mc_nu_df = concat_dfs['mcnu']
        mc_evt_df = concat_dfs['evt']        

        mc_tot_pot = mc_hdr_df['pot'].sum()
        print("mc_tot_pot: %.3e" %(mc_tot_pot))
        target_pot = mc_tot_pot
        mc_pot_scale = target_pot / mc_tot_pot
        print("mc_pot_scale: %.3e" %(mc_pot_scale))

        mc_evt_df["pot_weight"] = mc_pot_scale * np.ones(len(mc_evt_df))
        mc_nu_df["pot_weight"] = mc_pot_scale * np.ones(len(mc_nu_df))

        pot_str = get_pot_str(mc_tot_pot)
        return {"evt": mc_evt_df, "mcnu": mc_nu_df, "hdr": mc_hdr_df, "pot_str": pot_str}

    elif option == "selected_events":
        mc_dfs = load_and_concat_mc_dfs(
            # file_dir="/pnfs/sbnd/scratch/users/munjung/xsec/2025spring_v10_06_00_10",
            file_dir="/exp/sbnd/data/users/munjung/xsec/2025Spring_v10_06_00_09",
            chunk_tags=generate_tags("ba")[1:],
            df_tag="",
            keys2load=['hdr', 'evt'],
            n_max_concat=n_max_concat,
            sub_dir="MC",
            sample_dir="BNB_cosmics"
        )
        mc_hdr_df = mc_dfs['hdr']
        mc_evt_df = mc_dfs['evt']

        ## -- low E MC
        dirt_file = path.join(file_dir, "MC", "lowE", "lowE_dirt.df")
        dirt_dfs = load_dfs(dirt_file, 
                            ['hdr', 'evt'], 
                            n_max_concat=1)
        dirt_evt_df = dirt_dfs['evt']
        dirt_hdr_df = dirt_dfs['hdr']

        ## -- Data
        # data_file = path.join(file_dir, "data", "BNB", "_Fixed_mup.df")
        # data_file = "/exp/sbnd/data/users/munjung/xsec/2025Spring_v10_06_00_09/data/BNB/Gen1_mup.df"
        data_file = "/exp/sbnd/data/users/munjung/xsec/2025Spring_v10_06_00_09/data/BNB/sel_mup-data-1e20.df"
        # data_file = "/exp/sbnd/data/users/munjung/xsec/2025Spring_v10_06_00_09/data/BNB/Gen1_2prong_wcandidates.df"
        data_dfs = load_dfs(data_file, 
                            ['evt', 'hdr'], #, 'bnbpot'], 
                            n_max_concat=n_max_concat)
        data_evt_df = data_dfs['evt']
        data_hdr_df = data_dfs['hdr']

        ## -- Intime Data
        intime_dfs = load_and_concat_mc_dfs(
            file_dir=file_dir,
            chunk_tags=generate_tags("ab"),
            df_tag="",
            keys2load=['hdr', 'evt'],
            n_max_concat=n_max_concat,
            sub_dir="data",
            sample_dir="OffBeam"
        )
        intime_hdr_df = intime_dfs['hdr']
        intime_evt_df = intime_dfs['evt']

        # TODO
        # Data
        # data_bnbpot_df = data_dfs['bnbpot']
        data_tot_pot = data_hdr_df['pot'].sum()
        data_evt_df["pot_weight"] = np.ones(len(data_evt_df))
        print("data_tot_pot: %.3e" %(data_tot_pot))
        pot_str = get_pot_str(data_tot_pot)
        pot_label = f"Events / Bin (POT={pot_str})"

        data_gates = data_hdr_df.nbnbinfo.sum()
        print("data tot gates : %.3e" %(data_gates))

        # BNB MC
        mc_tot_pot = mc_hdr_df['pot'].sum()
        mc_pot_scale = data_tot_pot / mc_tot_pot
        print("mc_pot_scale: %.3e" %(mc_pot_scale))
        mc_evt_df["pot_weight"] = mc_pot_scale * np.ones(len(mc_evt_df))

        dirt_tot_pot = dirt_hdr_df['pot'].sum()
        print("dirt_tot_pot: %.3e" %(dirt_tot_pot))
        dirt_pot_scale = data_tot_pot / dirt_tot_pot
        print("dirt_pot_scale: %.3e" %(dirt_pot_scale))
        dirt_evt_df["pot_weight"] = dirt_pot_scale * np.ones(len(dirt_evt_df))

        # Intime Data
        # intime_gates = intime_hdr.noffbeambnb.sum()
        intime_gates = intime_hdr_df[intime_hdr_df['first_in_subrun'] == 1]['noffbeambnb'].sum()
        f = 0.0753
        scale_intime_to_lightdata = (1-f)*data_gates/intime_gates
        print("intime data scale: {:.2f}".format(scale_intime_to_lightdata))
        intime_evt_df["gates_weight"] = scale_intime_to_lightdata * np.ones(len(intime_evt_df))
        intime_evt_df["pot_weight"] = scale_intime_to_lightdata * np.ones(len(intime_evt_df))

        return {"mc": mc_evt_df, "data": data_evt_df, "intime": intime_evt_df, 
                "mc_hdr": mc_hdr_df, "data_hdr": data_hdr_df, "intime_hdr": intime_hdr_df,
                "dirt_hdr": dirt_hdr_df, "dirt": dirt_evt_df,
                "pot_label": pot_label}

    elif option == "data_unfolding":
        concat_dfs = load_and_concat_mc_dfs(
            file_dir=file_dir,
            # file_dir="/pnfs/sbnd/scratch/users/munjung/xsec/2025Spring_v10_06_00_09",
            # chunk_tags=generate_tags("bl"),
            # df_tag="_sel_mup-geniewgts",
            chunk_tags=generate_tags("ad")[1:],
            df_tag="-sel_all-wgts",
            # df_tag="",
            keys2load=['hdr', 'mcnu', 'evt'],
            n_max_concat=n_max_concat,
            sub_dir="MC",
            # sample_dir="BNB_cosmics/genie_wgts-MEC"
            sample_dir="BNB_cosmics"
        )
        mc_hdr_df = concat_dfs['hdr']
        mc_nu_df = concat_dfs['mcnu']
        mc_evt_df = concat_dfs['evt']        

        ## -- Data
        # data_file = path.join(file_dir, "data", "BNB", "_Fixed_mup.df")
        # data_file = path.join(file_dir, "data", "BNB", "_Rolling_mup.df")
        data_file = "/exp/sbnd/data/users/munjung/xsec/2025Spring_v10_06_00_09/data/BNB/Gen1_mup.df"
        data_dfs = load_dfs(data_file, 
                            ['evt', 'hdr'], #, 'bnbpot'], 
                            n_max_concat=n_max_concat)
        data_evt_df = data_dfs['evt']
        data_hdr_df = data_dfs['hdr']

        # TODO
        # Data
        # data_bnbpot_df = data_dfs['bnbpot']
        data_tot_pot = data_hdr_df['pot'].sum()
        data_evt_df["pot_weight"] = np.ones(len(data_evt_df))
        print("data_tot_pot: %.3e" %(data_tot_pot))
        pot_str = get_pot_str(data_tot_pot)
        pot_label = f"Events / Bin (POT={pot_str})"

        data_gates = data_hdr_df.nbnbinfo.sum()
        print("data tot gates : %.3e" %(data_gates))

        mc_tot_pot = mc_hdr_df['pot'].sum()
        print("mc_tot_pot: %.3e" %(mc_tot_pot))
        mc_pot_scale = data_tot_pot / mc_tot_pot
        print("mc_pot_scale: %.3e" %(mc_pot_scale))
        mc_evt_df["pot_weight"] = mc_pot_scale * np.ones(len(mc_evt_df))
        mc_nu_df["pot_weight"] = mc_pot_scale * np.ones(len(mc_nu_df))

        pot_str = get_pot_str(mc_tot_pot)
        return {"evt": mc_evt_df, "mcnu": mc_nu_df, "hdr": mc_hdr_df, "pot_str": pot_str,
                "data": data_evt_df, "data_hdr": data_hdr_df, "pot_label": pot_label}

    elif option == "data_unfolding_full":
        concat_dfs = load_and_concat_mc_dfs(
            # file_dir=file_dir,
            file_dir="/pnfs/sbnd/scratch/users/munjung/xsec/2025Spring_v10_06_00_09",
            chunk_tags=generate_tags("bl"),
            # df_tag="_sel_mup-geniewgts",
            df_tag="",
            keys2load=['hdr', 'mcnu', 'evt'],
            n_max_concat=n_max_concat,
            sub_dir="MC",
            sample_dir="BNB_cosmics/genie_wgts-MEC"
        )
        mc_hdr_df = concat_dfs['hdr']
        mc_nu_df = concat_dfs['mcnu']
        mc_evt_df = concat_dfs['evt']        

        ## -- Data
        data_file = path.join(file_dir, "data", "BNB", "_Gen1_mup.df")
        data_dfs = load_dfs(data_file, 
                            ['evt', 'hdr', 'bnbpot'], 
                            n_max_concat=n_max_concat)
        data_evt_df = data_dfs['evt']
        data_hdr_df = data_dfs['hdr']

        # TODO
        # Data
        # data_bnbpot_df = data_dfs['bnbpot']
        data_tot_pot = data_hdr_df['pot'].sum()
        data_evt_df["pot_weight"] = np.ones(len(data_evt_df))
        print("data_tot_pot: %.3e" %(data_tot_pot))
        pot_str = get_pot_str(data_tot_pot)
        pot_label = f"Events / Bin (POT={pot_str})"

        data_gates = data_hdr_df.nbnbinfo.sum()
        print("data tot gates : %.3e" %(data_gates))

        mc_tot_pot = mc_hdr_df['pot'].sum()
        print("mc_tot_pot: %.3e" %(mc_tot_pot))
        mc_pot_scale = data_tot_pot / mc_tot_pot
        print("mc_pot_scale: %.3e" %(mc_pot_scale))
        mc_evt_df["pot_weight"] = mc_pot_scale * np.ones(len(mc_evt_df))
        mc_nu_df["pot_weight"] = mc_pot_scale * np.ones(len(mc_nu_df))

        pot_str = get_pot_str(mc_tot_pot)
        return {"evt": mc_evt_df, "mcnu": mc_nu_df, "hdr": mc_hdr_df, "pot_str": pot_str,
                "data": data_evt_df, "data_hdr": data_hdr_df, "pot_label": pot_label}


    elif option == "event_selection":
        ## -- MC 
        mc_keys2load = ['hdr', 'evt', 'trk', 'mcnu']
        if syst_tag == "":
            mc_dfs = load_and_concat_mc_dfs(
                file_dir=file_dir,
                # file_dir="/pnfs/sbnd/scratch/users/munjung/xsec/2025Spring_v10_06_00_10/MC/BNB_cosmics/all-wgts",
                chunk_tags=generate_tags("ac")[1:],
                df_tag="-sel_all-wgts",
                # df_tag="",
                # chunk_tags=[""],
                # df_tag="evt_sel-test",
                keys2load=mc_keys2load,
                n_max_concat=n_max_concat,
                sub_dir="MC",
                sample_dir="BNB_cosmics"
                # sub_dir="",
                # sample_dir=""
            )

        elif syst_tag == "GiBUU":
            mc_dfs = load_and_concat_mc_dfs(
                file_dir=file_dir,
                chunk_tags=[""],
                df_tag="GiBUU-sel_all",
                keys2load=mc_keys2load,
                n_max_concat=n_max_concat,
                sub_dir="MC",
                sample_dir="BNB_cosmics"
            )

        else:
            mc_dfs = load_and_concat_mc_dfs(
                file_dir="/exp/sbnd/data/users/munjung/xsec/2025Spring_v10_06_00_10",
                chunk_tags=[f"SystVar_{syst_tag}_wtrks"],
                df_tag="",
                keys2load=['hdr', 'evt', 'trk'],
                n_max_concat=5,
                sub_dir="",
                sample_dir=""
            )

        mc_hdr_df = mc_dfs['hdr']
        mc_evt_df = mc_dfs['evt']
        mc_trk_df = mc_dfs['trk']
        # mc_mcnu_df = mc_dfs['mcnu']

        ## -- low E MC
        dirt_dfs = load_and_concat_mc_dfs(
            file_dir=file_dir,
            chunk_tags=[t for t in generate_tags("ad") if t != "ah"],
            df_tag="_all",
            keys2load=mc_keys2load,
            n_max_concat=n_max_concat,
            sub_dir="MC",
            sample_dir="lowE"
        )
        dirt_hdr_df = dirt_dfs['hdr']
        dirt_evt_df = dirt_dfs['evt']
        dirt_trk_df = dirt_dfs['trk']
        dirt_mcnu_df = dirt_dfs['mcnu']

        ## -- Data
        data_file = path.join(file_dir, "data", "BNB", "old", "_Fixed_all.df")
        data_dfs = load_dfs(data_file, 
                            ['evt', 'trk', 'hdr', 'bnbpot'], 
                            n_max_concat=n_max_concat)
        data_evt_df = data_dfs['evt']
        data_trk_df = data_dfs['trk']
        data_hdr_df = data_dfs['hdr']
        data_bnbpot_df = data_dfs['bnbpot']

        # -- Offbeam Data
        offbeam_keys2load = ['hdr', 'evt', 'trk']
        offbeam_dfs = load_and_concat_mc_dfs(
            file_dir=file_dir,
            chunk_tags=["test"],
            df_tag="_all",
            keys2load=offbeam_keys2load,
            n_max_concat=n_max_concat,
            sub_dir="data",
            sample_dir="OffBeam"
        )
        offbeam_hdr_df = offbeam_dfs['hdr']
        offbeam_evt_df = offbeam_dfs['evt']
        offbeam_trk_df = offbeam_dfs['trk']

        ## -- Intime MC
        intime_keys2load = ['hdr', 'evt', 'trk']
        intime_dfs = load_and_concat_mc_dfs(
            file_dir=file_dir,
            chunk_tags=generate_tags("ag"),
            df_tag="_all",
            keys2load=intime_keys2load,
            n_max_concat=n_max_concat,
            sub_dir="MC",
            sample_dir="intime"
        )
        intime_hdr_df = intime_dfs['hdr']
        intime_evt_df = intime_dfs['evt']
        intime_trk_df = intime_dfs['trk']

        ## exposure accounting

        # BNB data
        # data_tot_pot = data_hdr_df['TOR875'].sum()
        data_tot_pot = data_hdr_df['pot'].sum()
        print("data_tot_pot: %.3e" %(data_tot_pot))
        pot_str = get_pot_str(data_tot_pot)
        data_evt_df["pot_weight"] = np.ones(len(data_evt_df))
        data_trk_df["pot_weight"] = np.ones(len(data_trk_df))
        data_gates = data_hdr_df.nbnbinfo.sum()
        print("data tot gates : %.3e" %(data_gates))

        # BNB MC
        mc_tot_pot = mc_hdr_df['pot'].sum()
        print("mc_tot_pot: %.3e" %(mc_tot_pot))
        mc_pot_scale = data_tot_pot / mc_tot_pot
        print("mc_pot_scale: %.3e" %(mc_pot_scale))
        mc_evt_df["pot_weight"] = mc_pot_scale * np.ones(len(mc_evt_df))
        mc_trk_df["pot_weight"] = mc_pot_scale * np.ones(len(mc_trk_df))
        if syst_tag == "GiBUU":
            mc_evt_df["pot_weight"] *= mc_evt_df.mc.genweight.copy()

        dirt_tot_pot = dirt_hdr_df['pot'].sum()
        print("dirt_tot_pot: %.3e" %(dirt_tot_pot))
        dirt_pot_scale = data_tot_pot / dirt_tot_pot
        print("dirt_pot_scale: %.3e" %(dirt_pot_scale))
        dirt_evt_df["pot_weight"] = dirt_pot_scale * np.ones(len(dirt_evt_df))
        dirt_trk_df["pot_weight"] = dirt_pot_scale * np.ones(len(dirt_trk_df))

        # TODO
        f = 0.08
        offbeam_gates = offbeam_hdr_df.noffbeambnb.sum()
        offbeam_gates = offbeam_hdr_df[offbeam_hdr_df['first_in_subrun'] == 1]['noffbeambnb'].sum()
        print("offbeam cosmics data gates: {:.2e}".format(offbeam_gates))
        scale_offbeam_to_lightdata = (1-f)*data_gates/offbeam_gates
        print("goal scale: {:.2f}".format(scale_offbeam_to_lightdata))
        offbeam_evt_df["gates_weight"] = scale_offbeam_to_lightdata * np.ones(len(offbeam_evt_df))
        offbeam_evt_df["pot_weight"] = scale_offbeam_to_lightdata * np.ones(len(offbeam_evt_df))
        offbeam_trk_df["pot_weight"] = scale_offbeam_to_lightdata * np.ones(len(offbeam_trk_df))

        intime_gates = intime_hdr_df[intime_hdr_df['first_in_subrun'] == 1]['ngenevt'].sum()
        print("intime cosmics data gates: {:.2e}".format(intime_gates))
        scale_intime_to_lightdata = (1-f)*data_gates/intime_gates
        print("goal scale: {:.2f}".format(scale_intime_to_lightdata))
        intime_evt_df["gates_weight"] = scale_intime_to_lightdata * np.ones(len(intime_evt_df))
        intime_evt_df["pot_weight"] = scale_intime_to_lightdata * np.ones(len(intime_evt_df))
        intime_trk_df["pot_weight"] = scale_intime_to_lightdata * np.ones(len(intime_trk_df))

        pot_label = f"Events / Bin (POT={pot_str})"
        return {
            "mc": mc_evt_df,
            "data": data_evt_df,
            "offbeam": offbeam_evt_df,
            "intime": intime_evt_df,
            "dirt": dirt_evt_df,
            "mc_trk": mc_trk_df,
            "data_trk": data_trk_df,
            "offbeam_trk": offbeam_trk_df,
            "intime_trk": intime_trk_df,
            "dirt_trk": dirt_trk_df,
            "mc_hdr": mc_hdr_df,
            "data_hdr": data_hdr_df,
            "intime_hdr": intime_hdr_df,
            "pot_str": pot_str,
            "pot_label": pot_label,
        }

    else:
        raise ValueError("Invalid option: {}".format(option))
    


# # ==== Spring Gen 1 detvar samples ====
# file_dir = "/exp/sbnd/data/users/munjung/xsec/2025Spring_v10_06_00_10"

# mc_keys2load = ['hdr', 'evt'] 

# concat_dfs = load_and_concat_mc_dfs(
#     file_dir=file_dir,
#     chunk_tags=generate_tags("ad"),
#     df_tag="",
#     keys2load=mc_keys2load,
#     n_max_concat=n_max_concat,
#     sub_dir="MC",
#     sample_dir="BNB_cosmics"
# )

# mc_hdr_df = concat_dfs['hdr']
# mc_evt_df = concat_dfs['evt']
