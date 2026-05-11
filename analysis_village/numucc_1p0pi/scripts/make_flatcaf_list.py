import uproot
import os
from concurrent.futures import ProcessPoolExecutor

# Configuration
#input_dir = "/pnfs/sbnd/scratch/users/jaz8600/TestEfieldSim/cathodeSim_DriftVelSim_2026"
#output_list = "/exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/DENT/cathodeSim_DriftVelSim_2026.txt"
#input_dir = "/pnfs/sbnd/scratch/users/jaz8600/EField_R00/"
#output_list = "/exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/DENT/EField_R00.txt"
input_dir = "/pnfs/sbnd/scratch/users/jaz8600/EField_R30_Short/"
output_list = "/exp/sbnd/app/users/munjung/misc/filelists/MC/SBND/DENT/EField_R30_Short.txt"
MIN_BRANCHES = 2000  # Threshold to distinguish data from skeletons
NUM_WORKERS = 30  # Uses all available CPU cores

def check_file(filepath):
    """Function to check a single file; runs in parallel processes."""
    try:
        # We use 'with' to ensure the file handle is closed immediately
        with uproot.open(filepath) as file:
            if "recTree" in file:
                tree = file["recTree"]
                # len(tree.keys()) is metadata-only and fast
                if len(tree.keys()) > MIN_BRANCHES:
                    return filepath
    except Exception:
        pass
    return None

def main():
    all_root_files = []

    print(f"Gathering file list from: {input_dir}")
    for root, _, files in os.walk(input_dir):
        for filename in files:
            if filename.endswith(".root"):
                all_root_files.append(os.path.join(root, filename))

    print(f"Checking {len(all_root_files)} files using {NUM_WORKERS} cores...")

    # Start the parallel pool
    flat_caf_files = []
    with ProcessPoolExecutor(max_workers=NUM_WORKERS) as executor:
        # map handles the distribution of work automatically
        results = list(executor.map(check_file, all_root_files))

        # Filter out the None results
        flat_caf_files = [path for path in results if path is not None]

    # Save absolute paths
    with open(output_list, "w") as f:
        for path in sorted(flat_caf_files):
            f.write(os.path.abspath(path) + "\n")

    print(f"\nDone! Found {len(flat_caf_files)} validated flat CAFs.")

if __name__ == "__main__":
    main()
