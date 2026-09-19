#!/usr/bin/env python3 
import os
import sys
import time
import shlex
import datetime
import pathlib
#from TimeTools import *
import argparse
import tables
from pyanalib.ntuple_glob import NTupleGlob
import pandas as pd
import warnings

warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)
warnings.filterwarnings("ignore", category=tables.exceptions.NaturalNameWarning)
pd.set_option('future.no_silent_downcasting', True)

## Arguments
parser = argparse.ArgumentParser(
    description="Data frame maker command: process input flatcaf files and generate output dataframes.",
    epilog="""\
Examples:

  -- Use Pool
  $ python run_df_maker.py -c ./configs/cohpi_slcdf.py -o test_cohpi_slcdf -i input_0.root,input_1.root,...

  -- Use Grid (adding -ngrid to an integer > 0 will automatically submit grid jobs)
  $ python run_df_maker.py -ngrid 2 -c ./configs/cohpi_slcdf.py -o test_cohpi_slcdf -i input_0.root,input_1.root,...

  -- Note!!
  Output df files are sent to /pnfs/<exp>/scratch/users/<User>/cafpyana_out in Grid mode

  -- Knob-group configs (GENIE)
  When using configs that read GENIE_KNOB_GROUP from the environment, set it in the shell before
  invoking run_df_maker; it is forwarded into each grid worker so single-group jobs (HDF keys
  evt, mcnu, hdr) work under -ngrid. Flux multisim uses a single evt table (sel_mup-fluxwgts-knobgroups.py);
  FLUX_GROUP is not used.

  -- Cut / campaign nesting
  Pass -gsubdir <tag> (or export CAFPYANA_GRID_SUBDIR=<tag>) so all sample outputs for one cut
  land under $CAFPYANA_GRID_OUT_DIR/dfs/<tag>/<timestamp>__<output>/ (logs mirror under logs/<tag>/).
""",
    formatter_class=argparse.RawTextHelpFormatter  # Ensures line breaks are preserved
)
parser.add_argument('-c', dest='config', default="", help="Path to the data frame maker configuration file in ./configs, i.e.) -c ./configs/mcnu.py.")
parser.add_argument('-o', dest='output', default="", help="output data frame name prefix")
parser.add_argument('-i', dest='inputfiles', default="", help="input root file path, you can submit multiple files using comma, i.e.) -i input_0.root,input_1.root")
parser.add_argument('-l', dest='inputfilelist', default="", help="a file of list for input root files")
parser.add_argument('-ncpu', dest='NCPU', default=-1, type=int, help="Number of CPUs to run on. Default is set to number on server.")
parser.add_argument('-ngrid', dest='NGridJobs', default=0, type=int, help="Number of grid jobs. Default = 0, no grid submission.")
parser.add_argument('-nfile', dest='NFiles', default=0, type=int, help="Number of files to run. Default = 0, run all input files.")
parser.add_argument('-split', dest='SplitSize', default=1.0, type=float, help="Split size in GB before writing to HDF5. Default = 1.0 GB.")
parser.add_argument(
    '-gsubdir',
    dest='grid_subdir',
    default="",
    help=(
        "Optional subdirectory under dfs/ and logs/ for this campaign "
        "(e.g. a cut tag like fvfix). Overrides CAFPYANA_GRID_SUBDIR when set.\n"
        "Layout: $CAFPYANA_GRID_OUT_DIR/dfs/<subdir>/<timestamp>__<output>/"
    ),
)

args = parser.parse_args()


def _resolve_grid_subdir():
    """Return a safe relative path segment for nesting grid dfs/logs, or \"\"."""
    raw = (args.grid_subdir or os.environ.get("CAFPYANA_GRID_SUBDIR", "") or "").strip()
    if not raw:
        return ""
    # Allow nested tags (a/b) but reject absolute paths and ``..``.
    cleaned = raw.replace("\\", "/").strip("/")
    parts = []
    for part in cleaned.split("/"):
        part = part.strip()
        if not part or part == ".":
            continue
        if part == ".." or part.startswith("."):
            raise ValueError(
                "Invalid grid subdir %r (no '..' or hidden segments)" % raw
            )
        parts.append(part)
    if not parts:
        return ""
    return "/".join(parts)


def _maybe_write_syst_hist_var_config_snapshot(dest_dir):
    """For syst_histcounts jobs: freeze VariableConfig next to outputs at submit/run time."""
    cfg = getattr(args, "config", "") or ""
    if "syst_histcounts" not in cfg.replace("\\", "/"):
        return None
    try:
        from analysis_village.numucc_1p0pi.syst_histcounts import (
            VAR_CONFIG_SNAPSHOT_JSON_NAME,
            write_var_config_snapshot_json,
        )
    except Exception as ex:
        print("[run_df_maker] skip VariableConfig snapshot import (%s)" % ex)
        return None
    os.makedirs(dest_dir, exist_ok=True)
    out = os.path.join(dest_dir, VAR_CONFIG_SNAPSHOT_JSON_NAME)
    try:
        write_var_config_snapshot_json(out)
        print("[run_df_maker] wrote VariableConfig snapshot: %s" % out)
        return out
    except Exception as ex:
        print("[run_df_maker] failed writing VariableConfig snapshot (%s)" % ex)
        return None


def run_pool(output, inputs, nproc):
    os.nice(10)
    # Freeze binning next to the pool output (same snapshot also lands in each HDF).
    _maybe_write_syst_hist_var_config_snapshot(str(pathlib.Path(output).resolve().parent))
    ntuples = NTupleGlob(inputs, None)

    # if PREPROCESS doesn't exist, set it to None
    global PREPROCESS
    try:
        PREPROCESS
    except:
        PREPROCESS = []

    dfss = ntuples.dataframes(nproc=nproc, args=ARGS, fs=DFS, preprocess=PREPROCESS)
    output = pathlib.Path(output).with_suffix('.df')
    k_idx = 0
    split_margin = args.SplitSize
    with pd.HDFStore(output) as hdf_pd:
        NAMES.append("histpotdf")
        NAMES.append("histgenevtdf")
        size_counters = {k: 0 for k in NAMES}
        df_buffers = {k: [] for k in NAMES}

        for dfs in dfss:
            if not dfs:
                continue

            # Align table names with returned frames. Never zip unequal lengths —
            # that mislabeled hdr as var_configs when makers failed silently.
            if len(dfs) == 2:
                # no / empty recTree: only histpot + histgenevt from ntuple_glob
                name_df_pairs = list(zip(["histpotdf", "histgenevtdf"], dfs))
            elif len(dfs) == len(NAMES):
                name_df_pairs = list(zip(NAMES, dfs))
            else:
                print(
                    "[run_df_maker] ERROR: len(dfs)=%d != len(NAMES)=%d (NAMES=%s); "
                    "skipping this CAF result to avoid mislabeled HDF keys"
                    % (len(dfs), len(NAMES), NAMES)
                )
                continue

            for k, df in name_df_pairs:
                if df is None:
                    continue
                size_bytes = df.memory_usage(deep=True).sum()
                size_gb = size_bytes / (1024**3)
                size_counters[k] += size_gb
                df_buffers[k].append(df)
                del df

            if any(val > split_margin for val in size_counters.values()):
                # Concatenate and save accumulated DataFrames
                for k, buffer in df_buffers.items():
                    if buffer:  # only if buffer has data
                        concat_df = pd.concat(buffer, ignore_index=False)
                        this_key = k + "_" + str(k_idx)
                        try:
                            if k == "syst_hists":
                                from analysis_village.numucc_1p0pi.syst_histcounts import (
                                    put_syst_hists_by_var,
                                )
                                written = put_syst_hists_by_var(
                                    hdf_pd, concat_df, split_idx=k_idx, format="fixed"
                                )
                                print(
                                    f"Saved syst_hists split {k_idx} "
                                    f"({len(written)} keys, "
                                    f"{concat_df.memory_usage(deep=True).sum() / (1024**3):.4f} GB)"
                                )
                            else:
                                hdf_pd.put(key=this_key, value=concat_df, format="fixed")
                                print(f"Saved {this_key}: {concat_df.memory_usage(deep=True).sum() / (1024**3):.4f} GB")
                        except Exception as e:
                            print(f"Table {this_key} failed to save, skipping. Exception: {str(e)}")
                        del concat_df
                # Reset counters and buffers
                k_idx += 1
                size_counters = {k: 0 for k in NAMES}
                df_buffers = {k: [] for k in NAMES}

        for k, buffer in df_buffers.items():
            if buffer:
                concat_df = pd.concat(buffer, ignore_index=False)
                this_key = k + "_" + str(k_idx)
                try:
                    if k == "syst_hists":
                        from analysis_village.numucc_1p0pi.syst_histcounts import (
                            put_syst_hists_by_var,
                        )
                        written = put_syst_hists_by_var(
                            hdf_pd, concat_df, split_idx=k_idx, format="fixed"
                        )
                        print(
                            f"Saved syst_hists split {k_idx} "
                            f"({len(written)} keys, "
                            f"{concat_df.memory_usage(deep=True).sum() / (1024**3):.4f} GB)"
                        )
                    else:
                        hdf_pd.put(key=this_key, value=concat_df, format="fixed")
                        print(f"Saved {this_key}: {concat_df.memory_usage(deep=True).sum() / (1024**3):.4f} GB")
                except Exception as e:
                    print(f"Table {this_key} failed to save, skipping. Exception: {str(e)}")
                del concat_df

        # Save the split count metadata
        split_df = pd.DataFrame({"n_split": [k_idx + 1]})  # +1 because k_idx is 0-based
        hdf_pd.put(key="split", value=split_df, format="fixed")
        print(f"Saved split info: {split_df.iloc[0]['n_split']} total splits")

def run_grid(inputfiles):
    # 1) dir/file name style
    JobStartTime = datetime.datetime.now()
    timestamp =  JobStartTime.strftime('%Y_%m_%d_%H%M%S')

    # 2) Define MasterJobDir -- produce grid job submission scripts in $CAFPYANA_GRID_OUT_DIR
    # Optional cut/campaign nesting: dfs/<subdir>/<stamp>__<tag>/  (and same under logs/)
    CAFPYANA_GRID_OUT_DIR = os.environ['CAFPYANA_GRID_OUT_DIR']
    grid_subdir = _resolve_grid_subdir()
    dfs_root = os.path.join(CAFPYANA_GRID_OUT_DIR, "dfs")
    logs_root = os.path.join(CAFPYANA_GRID_OUT_DIR, "logs")
    if grid_subdir:
        dfs_root = os.path.join(dfs_root, grid_subdir)
        logs_root = os.path.join(logs_root, grid_subdir)
        print("[run_df_maker] grid subdir: %s" % grid_subdir)
    MasterJobDir = os.path.join(logs_root, timestamp + '__' + args.output + "_log")
    OutputDir = os.path.join(dfs_root, timestamp + '__' + args.output)
    os.system('mkdir -p ' + MasterJobDir)
    os.system('mkdir -p ' + OutputDir)
    # Freeze VariableConfig at submit time into the campaign output + log dirs.
    _maybe_write_syst_hist_var_config_snapshot(OutputDir)
    _maybe_write_syst_hist_var_config_snapshot(MasterJobDir)

    # 3) grid job is based on number of files
    ngrid = args.NGridJobs
    if(len(inputfiles) <= ngrid):
        ngrid = len(inputfiles)

    NInputfiles = len(inputfiles)
    print("Number of Grid Jobs: %d, number of input caf files: %d" % (ngrid, NInputfiles))

    # 4) prepare bash scripts for each job and make tarball
    flistForEachJob = []
    for i in range(0,ngrid):
        flistForEachJob.append( [] )

    for i_line in range(0,len(inputfiles)):
        flistForEachJob[i_line%ngrid].append(inputfiles[i_line])

    for i_flist in range(0,len(flistForEachJob)):
        flist = flistForEachJob[i_flist]
        out = open(MasterJobDir + '/run_%s.sh'%(i_flist),'w')
        out.write('#!/bin/bash\n')
        # Worker jobs do not inherit the submit-shell environment; configs that branch on e.g.
        # GENIE_KNOB_GROUP must see the same values as the submit host (FLUX_GROUP unused for flux df).
        for _env in ("GENIE_KNOB_GROUP", "FLUX_GROUP", "SYST_HIST_MODE", "SYST_HIST_SAMPLE", "SYST_HIST_EXCLUDE_AR23P", "SYST_HIST_INCLUDE_AR23P", "SYST_HIST_EXCLUDE_SLIM", "SYST_HIST_GENIE_NUNIV", "SYST_HIST_FLUX_NUNIV", "SYST_HIST_G4_NUNIV"):
            _v = os.environ.get(_env, "").strip()
            if _v:
                out.write("export %s=%s\n" % (_env, shlex.quote(_v)))
        # Inject submit-host config (grid clones GitHub; local configs may be uncommitted).
        _cfg_rel = args.config
        _cfg_base = os.path.basename(_cfg_rel)
        _cfg_dir = os.path.dirname(_cfg_rel) or "."
        out.write(
            'if [ -f "${CONDOR_DIR_INPUT}/bin_dir/%s" ]; then\n'
            '  mkdir -p %s\n'
            '  cp -f "${CONDOR_DIR_INPUT}/bin_dir/%s" %s/\n'
            '  echo "[run_%s.sh] injected config %s from submit host"\n'
            'fi\n'
            % (_cfg_base, shlex.quote(_cfg_dir), _cfg_base, shlex.quote(_cfg_dir), i_flist, _cfg_base)
        )
        # Inject pickleable updatecalo makers (and any other local makedf fixes).
        out.write(
            'if [ -f "${CONDOR_DIR_INPUT}/bin_dir/numucc_makedf.py" ]; then\n'
            '  mkdir -p analysis_village/numucc_1p0pi/makedf\n'
            '  cp -f "${CONDOR_DIR_INPUT}/bin_dir/numucc_makedf.py" '
            'analysis_village/numucc_1p0pi/makedf/makedf.py\n'
            '  echo "[run_%s.sh] injected analysis_village/.../makedf.py from submit host"\n'
            'fi\n'
            % i_flist
        )
        # Inject cut thresholds / VERTEX_Z_EXCLUDE (cut campaigns; must match pushed branch).
        out.write(
            'if [ -f "${CONDOR_DIR_INPUT}/bin_dir/numucc_selections.py" ]; then\n'
            '  mkdir -p analysis_village/numucc_1p0pi/makedf\n'
            '  cp -f "${CONDOR_DIR_INPUT}/bin_dir/numucc_selections.py" '
            'analysis_village/numucc_1p0pi/makedf/selections.py\n'
            '  echo "[run_%s.sh] injected analysis_village/.../selections.py from submit host"\n'
            'fi\n'
            % i_flist
        )
        # FSI_compare / slim-throw packs: ship local geniesyst + syst_histcounts (not yet on GitHub).
        out.write(
            'if [ -f "${CONDOR_DIR_INPUT}/bin_dir/numucc_geniesyst.py" ]; then\n'
            '  mkdir -p makedf\n'
            '  cp -f "${CONDOR_DIR_INPUT}/bin_dir/numucc_geniesyst.py" makedf/geniesyst.py\n'
            '  echo "[run_%s.sh] injected makedf/geniesyst.py from submit host"\n'
            'fi\n'
            'if [ -f "${CONDOR_DIR_INPUT}/bin_dir/numucc_syst_histcounts.py" ]; then\n'
            '  mkdir -p analysis_village/numucc_1p0pi\n'
            '  cp -f "${CONDOR_DIR_INPUT}/bin_dir/numucc_syst_histcounts.py" '
            'analysis_village/numucc_1p0pi/syst_histcounts.py\n'
            '  echo "[run_%s.sh] injected analysis_village/.../syst_histcounts.py from submit host"\n'
            'fi\n'
            % (i_flist, i_flist)
        )
        cmd = 'python run_df_maker.py -c ' + args.config + ' -o ' + args.output + '_%d'%i_flist + '.df -ncpu 7 -i'
        for i_f in range(0,len(flist)):
            out.write('echo "[run_%s.sh] input %d : %s"\n'%(i_flist, i_f, flist[i_f]))
            if i_f == 0:
                cmd += ' ' + flist[i_f].split('/')[-1]
            else: 
                cmd += ',' + flist[i_f].split('/')[-1]
            out.write('xrdcp ' + flist[i_f] + ' .\n') ## -- for checking auth
        out.write('ls -alh\n')
        out.write(cmd)
        out.close()

    os.system('cp ./bin/grid_executable.sh %s' %MasterJobDir)
    # Ship the exact config + local makedf (see inject block in run_*.sh).
    try:
        import shutil
        shutil.copy2(os.path.abspath(args.config), os.path.join(MasterJobDir, os.path.basename(args.config)))
        print("[run_df_maker] bundled config into job tarball:", os.path.basename(args.config))
        _wd = os.environ.get("CAFPYANA_WD", os.getcwd())
        _makedf_local = os.path.join(_wd, "analysis_village/numucc_1p0pi/makedf/makedf.py")
        if os.path.isfile(_makedf_local):
            shutil.copy2(_makedf_local, os.path.join(MasterJobDir, "numucc_makedf.py"))
            print("[run_df_maker] bundled numucc makedf.py into job tarball")
        _sel_local = os.path.join(
            _wd, "analysis_village/numucc_1p0pi/makedf/selections.py"
        )
        if os.path.isfile(_sel_local):
            shutil.copy2(_sel_local, os.path.join(MasterJobDir, "numucc_selections.py"))
            print("[run_df_maker] bundled numucc selections.py into job tarball")
        _geniesyst_local = os.path.join(_wd, "makedf/geniesyst.py")
        if os.path.isfile(_geniesyst_local):
            shutil.copy2(_geniesyst_local, os.path.join(MasterJobDir, "numucc_geniesyst.py"))
            print("[run_df_maker] bundled makedf/geniesyst.py into job tarball")
        _syst_hc_local = os.path.join(
            _wd, "analysis_village/numucc_1p0pi/syst_histcounts.py"
        )
        if os.path.isfile(_syst_hc_local):
            shutil.copy2(
                _syst_hc_local, os.path.join(MasterJobDir, "numucc_syst_histcounts.py")
            )
            print("[run_df_maker] bundled syst_histcounts.py into job tarball")
    except Exception as ex:
        print("[run_df_maker] WARNING: could not bundle config/makedf into tarball:", ex)

    # 5) prepare a package for xrootd (prefer installed 5.6.9 build; fall back to legacy 5.6.1 path)
    CAFPYANA_WD = os.environ['CAFPYANA_WD']
    xrd_candidates = [
        CAFPYANA_WD + "/envs/xrootd-5.6.9/build/lib.linux-x86_64-cpython-310",
        CAFPYANA_WD + "/envs/xrootd-5.6.9/build/lib.linux-x86_64-3.9",
        CAFPYANA_WD + "/envs/xrootd-5.6.1/build/lib.linux-x86_64-3.9",
    ]
    xrd_lib = None
    for cand in xrd_candidates:
        if os.path.isdir(os.path.join(cand, "XRootD")) and os.path.isdir(os.path.join(cand, "pyxrootd")):
            xrd_lib = cand
            break
    if xrd_lib is None:
        raise RuntimeError(
            "Could not find XRootD/pyxrootd under envs/xrootd-*; tried:\n  "
            + "\n  ".join(xrd_candidates)
        )
    print("[run_df_maker] bundling XRootD from", xrd_lib)
    os.system("cp -r " + xrd_lib + "/XRootD " + MasterJobDir)
    os.system("cp -r " + xrd_lib + "/pyxrootd " + MasterJobDir)

    os.chdir(MasterJobDir)
    tar_cmd = 'tar cf bin_dir.tar ./'
    os.system(tar_cmd)

    # Resource requests: override via env for heavy histcount / weight jobs.
    job_disk = os.environ.get("JOBSUB_DISK", "10GB").strip() or "10GB"
    job_mem = os.environ.get("JOBSUB_MEMORY", "5GB").strip() or "10GB"
    job_life = os.environ.get("JOBSUB_LIFETIME", "3h").strip() or "3h"
    job_cpu = os.environ.get("JOBSUB_CPU", "7").strip() or "7"
    print(
        "[run_df_maker] jobsub resources: disk=%s memory=%s lifetime=%s cpu=%s"
        % (job_disk, job_mem, job_life, job_cpu)
    )

    submitCMD = '''jobsub_submit \\
-G sbnd \\
--auth-methods="token" \\
-e LC_ALL=C \\
--role=Analysis \\
--resource-provides="usage_model=DEDICATED,OPPORTUNISTIC" \\
--lines '+FERMIHTC_AutoRelease=True' --lines '+FERMIHTC_GraceMemory=5000' --lines '+FERMIHTC_GraceLifetime=3600' \\
--append_condor_requirements='(TARGET.HAS_SINGULARITY=?=true)' \\
--tar_file_name "dropbox://$(pwd)/bin_dir.tar" \\
-N %d \\
--disk %s \\
--cpu %s \\
--memory %s \\
--expected-lifetime %s \\
"file://$(pwd)/grid_executable.sh" \\
"%s" \\
"%s"''' % (ngrid, job_disk, job_cpu, job_mem, job_life, OutputDir, args.output)

    print(submitCMD)
    os.system(submitCMD)
    
    # go back to working dir
    os.chdir(CAFPYANA_WD)
    
if __name__ == "__main__":
    printhelp = ((args.inputfiles == "" and args.inputfilelist == "") or args.config == "" or args.output == "")
    if printhelp:
        parser.print_help()
        print(parser.epilog)
        sys.exit(1)
        
    else:
        ### Organize input list 
        InputSamples = []
        StringForHash = ""
        if args.inputfilelist != "":
            lines = open(args.inputfilelist)
            for line in lines:
                if "#" in line:
                    continue
                line = line.strip('\n')
                InputSamples.append(line)
                StringForHash += line
        else:
            split_inputfiles = args.inputfiles.split(",")
            for split_inputfile in split_inputfiles:
                InputSamples.append(split_inputfile)
                StringForHash += args.inputfiles

        if(args.NFiles > 0 and len(InputSamples) > args.NFiles):
            InputSamples = InputSamples[:args.NFiles]
                
        ### check if it is grid mode for pool mode
        if args.NGridJobs == 0:
            print("Running Pool mode");
            exec(open(args.config).read())
            run_pool(args.output, InputSamples, "auto" if args.NCPU < 0 else args.NCPU)

        elif args.NGridJobs > 0:
            print("Running Grid mode");
            run_grid(InputSamples)
            
        else:
            print("-ngrid must be greater than 0.");
