import numpy as np
import pandas as pd
import argparse
import os
import sys
import glob
import warnings
import matplotlib.pyplot as plt

# Try to import pyGENIE; ensure your environment is sourced first
try:
    import pyGENIE as pg
except ImportError:
    print("Error: pyGENIE not found. Please setup your LArSoft/GENIE environment.")
    sys.exit(1)

# --- Global Configurations ---
pd.options.mode.chained_assignment = None # Silence SettingWithCopy
warnings.filterwarnings("ignore", category=RuntimeWarning) # Silence NaN in sqrt

def get_args():
    parser = argparse.ArgumentParser(description='GENIE GHEP to Event-Level Pandas DataFrame')
    parser.add_argument('-i', '--input', type=str, required=True,
                        help="Input: .ghep.root file, .txt list of files, or a directory path.")
    parser.add_argument('-n', '--max_files', type=int, default=-1,
                        help="Limit the number of files to process (Default: -1 for all).")
    parser.add_argument('-e', '--max_events', type=int, default=-1,
                        help="Limit events per file for testing (Default: -1 for all).")
    parser.add_argument('-o', '--output', type=str, default="event_summary.pkl",
                        help="Output filename (e.g., summary.pkl or summary.csv).")
    return parser.parse_args()

def main():
    args = get_args()
    print(f"Starting Analysis: {args.input}")

    # 1. Resolve Input Files
    input_files = []
    if os.path.isdir(args.input):
        input_files = sorted(glob.glob(os.path.join(args.input, "*.ghep.root")))
    elif args.input.endswith('.txt'):
        with open(args.input, 'r') as f:
            input_files = [line.strip() for line in f if line.strip() and not line.startswith('#')]
    else:
        input_files = [args.input]

    if args.max_files > 0:
        input_files = input_files[:args.max_files]

    if not input_files:
        print("No valid GHEP files found. Check your input path.")
        return

    all_events = []
    f_idx = 0

    # 2. Process Files
    for ghep in input_files:
        if not os.path.exists(ghep):
            print(f"Skipping missing file: {ghep}")
            continue

        with pg.open_file(ghep, "gtree") as (fin, gtree):
            nEntries = gtree.GetEntries()
            loop_limit = nEntries if args.max_events < 0 else min(nEntries, args.max_events)
            
            print(f"File [{f_idx+1}/{len(input_files)}]: {os.path.basename(ghep)} ({loop_limit} entries)")

            with pg.silence_cpp():
                # Attach the record holders (GENIE + Flux)
                gtree, holder, fholder = pg.attach_record_with_flux(gtree, "gmcrec", "dk2nu")

            # 3. Event Loop
            for iev in range(loop_limit):
                with pg.silence_cpp():
                    gtree.GetEntry(iev)
                    rec, dk2nu = holder.rec, fholder.dk2nu
                    # Convert GHEP record to Particle-Level DataFrame
                    rec_df = pg.ghep_to_pandas(rec, event_id=iev, file_id=f_idx).copy()

                # --- Kinematics Calculation ---
                # Filter for Initial (0) and Final State (1) only
                prim_df = rec_df[rec_df.status.isin([0, 1])].copy()
                
                # Momentum p = sqrt(px^2 + py^2 + pz^2)
                prim_df['p'] = np.sqrt(prim_df.px**2 + prim_df.py**2 + prim_df.pz**2)
                prim_df['cos_th'] = prim_df.pz / prim_df.p
                
                # Mass and Kinetic Energy (KE = E - m)
                m2 = prim_df.E**2 - prim_df.px**2 - prim_df.py**2 - prim_df.pz**2
                prim_df['mass'] = np.sqrt(m2.clip(lower=0)) 
                prim_df['ke'] = prim_df.E - prim_df.mass

                # --- Event-Level Metadata ---
                nu_pdg = prim_df.pdg.iloc[0]
                target_pdg = prim_df.pdg.iloc[1]
                vtx = {'x': prim_df.vx.iloc[0], 'y': prim_df.vy.iloc[0], 'z': prim_df.vz.iloc[0]}

                # --- Final State Topology (Status 1) ---
                final_state = prim_df[prim_df.status == 1]
                
                # Muon Selection
                muon_data = final_state[np.abs(final_state.pdg) == 13]
                muon_p = muon_data.p.iloc[0] if not muon_data.empty else np.nan
                muon_cth = muon_data.cos_th.iloc[0] if not muon_data.empty else np.nan

                # Multiplicities
                n_mu = ( (final_state.pdg == 13) & (final_state.p > 0.22) & (final_state.p < 1.0)).sum()
                n_p = ( (np.abs(final_state.pdg) == 2212) & (final_state.p > 0.3) & (final_state.p < 1.0)).sum()
                n_p_add = ( (np.abs(final_state.pdg) == 2212) & (final_state.p > 0.3)).sum()
                n_cpi = ( (np.abs(final_state.pdg) == 211) & (final_state.p > 0.07) ).sum()
                n_pi0 = ( (np.abs(final_state.pdg) == 111) & (final_state.p > 0.0) ).sum()

                # --- Aggregate to Event Data ---
                event_data = {
                    'file_id': f_idx,
                    'event_id': iev,
                    'nu_pdg': nu_pdg,
                    'target_pdg': target_pdg,
                    'muon_p': muon_p,
                    'muon_cos_th': muon_cth,
                    'n_mu_220MeVc': n_mu,
                    'n_p_300MeVc': n_p,
                    'n_p_300MeVc_add': n_p_add,
                    'n_cpi_70MeVc': n_cpi,
                    'n_pi0_0MeV': n_pi0,
                    'vtx_x': vtx['x'],
                    'vtx_y': vtx['y'],
                    'vtx_z': vtx['z'],
                    # 'pot': dk2nu.potnum
                }
                all_events.append(event_data)

        f_idx += 1

    # 4. Final Output
    event_df = pd.DataFrame(all_events)
    print("\n--- Final Event-Level DataFrame ---")
    print(event_df.head())
    
    # Save the dataframe
    if args.output.endswith('.csv'):
        event_df.to_csv(args.output, index=False)
    else:
        event_df.to_pickle(args.output)
    
    print(f"\nSaved {len(event_df)} events to {args.output}")
    print("Goodbye!")

if __name__ == "__main__":
    main()
