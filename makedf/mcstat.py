import numpy as np
import pandas as pd
from tqdm import tqdm


def mcstatsyst(hdr_df, nuind, multisim_nuniv=100, slim=True, poisson_mean=1.0):
    """
    MC statistical uncertainty as Poisson throws per neutrino interaction row.

    Uses the same pseudo-random convention as :func:`get_MCstat_unc` (per-row
    metadata seed XOR universe seed, then ``numpy.random.poisson``).  Rows are
    keyed by ``run``, ``subrun``, ``evt`` from ``hdr_df`` and neutrino index
    ``inu`` from ``nuind``, analogous to slice id in :func:`get_MCstat_unc`.

    Aligns with ``makedf.getsyst`` slim output: index ``nuind.index``
    ``(entry, inu)``, columns ``("MCstat", "univ_<i>")``.  Only MCstat
    universes exist, so column layout is unchanged when ``slim`` is False.

    Parameters
    ----------
    hdr_df : pandas.DataFrame
        Headers indexed by ``entry``, with ``run``, ``subrun``, ``evt``.
    nuind : pandas.Series
        Neutrino index series (e.g. ``mcdf.ind`` from :func:`make_mcnudf`).
    multisim_nuniv : int
        Number of Poisson universes.
    slim : bool
        Accepted for API parity with ``g4syst`` / ``bnbsyst``; unused (layout
        is always the slim MCstat multisim shape).
    poisson_mean : float
        Mean of the Poisson distribution (typically 1 for nominal-weight throws).
    """
    _ = slim  # API parity with g4syst/bnbsyst; MCstat layout matches slim multisim.

    index = nuind.index
    entry = index.get_level_values(0)
    hdr_aligned = hdr_df.reindex(entry)

    runs = np.asarray(hdr_aligned.run.values, dtype=object)
    subruns = np.asarray(hdr_aligned.subrun.values, dtype=object)
    evts = np.asarray(hdr_aligned.evt.values, dtype=object)
    inus = np.asarray(nuind.values)

    n_rows = len(index)
    meta_seeds = np.empty(n_rows, dtype=np.int64)
    for j in range(n_rows):
        meta_seeds[j] = hash(
            f"run_{runs[j]}_subrun_{subruns[j]}_evt_{evts[j]}_inu_{inus[j]}"
        ) % (2**32)

    cols = pd.MultiIndex.from_product(
        [["MCstat"], [f"univ_{i}" for i in range(multisim_nuniv)]],
    )
    data = {}
    for uidx in range(multisim_nuniv):
        universe_seed = hash(f"universe_{uidx}") % (2**32)
        col = np.empty(n_rows, dtype=np.float64)
        for j in range(n_rows):
            combined_seed = int((universe_seed + int(meta_seeds[j])) % (2**32))
            np.random.seed(combined_seed)
            col[j] = np.random.poisson(poisson_mean)
        data[("MCstat", f"univ_{uidx}")] = col

    return pd.DataFrame(data, index=index, columns=cols)


def get_MCstat_unc(evt_df, hdr_df, n_universes=100):
    """
    Create a unique seed based on event metadata
    Using a hash function that's deterministic
    """

    meta_seeds = []
    for i in tqdm(range(len(evt_df))):
        this_hdr_df = hdr_df.loc[evt_df.reset_index(level=[2]).index[i]]
        runno = this_hdr_df.run
        subrunno = this_hdr_df.subrun
        evtno = this_hdr_df.evt
        slcid = evt_df.loc[evt_df.index[i]].slc.self
        unique_seed = hash(f"run_{runno}_subrun_{subrunno}_evt_{evtno}_slcid_{slcid}") % (2**32)  # Ensure it's a 32-bit integer
        if unique_seed in meta_seeds:
            print("duplicate seed found", unique_seed)
            break
        meta_seeds.append(unique_seed)
    # make sure the seeds are unique!
    assert len(meta_seeds) == len(set(meta_seeds))

    # generate universes
    MCstat_univ_events = np.zeros((n_universes, len(evt_df)))
    poisson_mean = 1.0
    # get Poisson weights and save to "MCstat.univ_"
    # TODO: generalize column level padding
    mcstat_univ_cols = pd.MultiIndex.from_product(
        [["MCstat"], [f"univ_{i}" for i in range(n_universes)], [""], [""], [""], [""], [""]],
    )
    mcstat_univ_wgt = pd.DataFrame(
        1.0,
        index=evt_df.index,
        columns=mcstat_univ_cols,
    )
    for uidx in range(n_universes):
        universe_seed = hash(f"universe_{uidx}") % (2**32)
        
        poisson_weights = []
        for meta_seed in meta_seeds:
            # Combine universe seed with event seed for unique randomness -- per event, per universe
            combined_seed = (universe_seed + meta_seed) % (2**32)
            np.random.seed(combined_seed)
            
            poisson_val = np.random.poisson(poisson_mean)
            poisson_weights.append(poisson_val)
        
        mcstat_univ_wgt[("MCstat", "univ_{}".format(uidx), "", "", "", "", "")] = np.array(poisson_weights)
        MCstat_univ_events[uidx, :] = np.array(poisson_weights)

    evt_df = evt_df.join(mcstat_univ_wgt)

    return evt_df, MCstat_univ_events
