import os
import sys
from pathlib import Path

import awkward as ak
import matplotlib.pyplot as plt
import numpy as np
import uproot
import pickle

from multiprocessing import Pool
from tqdm import tqdm


def load_dk2nu(filename):
    dkf = uproot.open(filename)
    dk_tree = dkf["dk2nuTree"]
    pots = float(dkf["dkmetaTree"]["dkmeta/pots"].array()[0])
    # print("dk2nu simulated POT (dkmeta): %.6g" % pots)

    branches = [
        "dk2nu/decay/decay.ntype",
        "dk2nu/decay/decay.nimpwt",
        "dk2nu/nuray/nuray.E",
        "dk2nu/nuray/nuray.wgt",
        'dk2nu/nuray/nuray.px',
        'dk2nu/nuray/nuray.py',
        'dk2nu/nuray/nuray.pz',
        'dk2nu/decay/decay.vx',
        'dk2nu/decay/decay.vy',
        'dk2nu/decay/decay.vz',
    ]
    arrays = dk_tree.arrays(branches, library="ak")
    mask_numu = arrays["dk2nu/decay/decay.ntype"] == 14
    nimpwt = arrays["dk2nu/decay/decay.nimpwt"][mask_numu]
    vx = arrays["dk2nu/decay/decay.vx"][mask_numu]
    vy = arrays["dk2nu/decay/decay.vy"][mask_numu]
    vz = arrays["dk2nu/decay/decay.vz"][mask_numu]
    E_flat = arrays["dk2nu/nuray/nuray.E"][mask_numu] #[:, i]
    w_flat = arrays["dk2nu/nuray/nuray.wgt"][mask_numu] #[:, i]
    px = arrays["dk2nu/nuray/nuray.px"][mask_numu]
    py = arrays["dk2nu/nuray/nuray.py"][mask_numu]
    pz = arrays["dk2nu/nuray/nuray.pz"][mask_numu]

    dkf.close()

    return {"E_flat": np.array(E_flat), 
            "px": np.array(px),
            "py": np.array(py),
            "pz": np.array(pz),
            "vx": np.array(vx),
            "vy": np.array(vy),
            "vz": np.array(vz),
            "w_flat": np.array(w_flat), "nimpwt": np.array(nimpwt)}


if __name__ == "__main__":

    DK2NU_DIR =  "/cvmfs/sbnd.osgstorage.org/pnfs/fnal.gov/usr/sbnd/persistent/stash/fluxFiles/bnb/G4BNB/v1.1.1/fhc/a"
    dk2nu_files = [f for f in os.listdir(DK2NU_DIR) if f.endswith(".dk2nu.root")]
    dk2nu_files_fullpath = [os.path.join(DK2NU_DIR, f) for f in dk2nu_files][:100]

    def load_dk2nu_helper(args):
        idx, fname = args
        ret = load_dk2nu(fname)
        return idx, ret

    with Pool(10) as pool:
        results = list(tqdm(pool.imap(load_dk2nu_helper, enumerate(dk2nu_files_fullpath)), total=len(dk2nu_files_fullpath)))

    # Sort results by file index for stable concatenation order
    results.sort(key=lambda x: x[0])
    rets = [r[1] for r in results]

    E_flat_arr = np.concatenate([r["E_flat"] for r in rets], axis=0)
    nimpwt_arr = np.concatenate([r["nimpwt"] for r in rets], axis=0)
    w_flat_arr = np.concatenate([r["w_flat"] for r in rets], axis=0)
    px_arr = np.concatenate([r["px"] for r in rets], axis=0)
    py_arr = np.concatenate([r["py"] for r in rets], axis=0)
    pz_arr = np.concatenate([r["pz"] for r in rets], axis=0)
    vx_arr = np.concatenate([r["vx"] for r in rets], axis=0)
    vy_arr = np.concatenate([r["vy"] for r in rets], axis=0)
    vz_arr = np.concatenate([r["vz"] for r in rets], axis=0)

    save_filename = "/exp/sbnd/data/users/munjung/xsec/flux_closure/sbnd_dk2nu_arrays.pkl"
    with open(save_filename, "wb") as f:
        pickle.dump({
            "E_flat_arr": E_flat_arr,
            "nimpwt_arr": nimpwt_arr,
            "w_flat_arr": w_flat_arr,
            "px_arr": px_arr,
            "py_arr": py_arr,
            "pz_arr": pz_arr,
            "vx_arr": vx_arr,
            "vy_arr": vy_arr,
            "vz_arr": vz_arr,
        }, f)

    print(f"Saved flux arrays to {save_filename}")
