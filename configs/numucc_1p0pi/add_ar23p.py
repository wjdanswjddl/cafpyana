from preprocess.preprocess import Script

from analysis_village.numucc_1p0pi.makedf.makedf import *

PREPROCESS = [Script("/exp/sbnd/app/users/munjung/xsec/cafpyana_2026Jan17/cafpyana/preprocess/update_reweight_anywhere.sh")]

DFS, ARGS, NAMES = build_genie_knobgroup_config("Ar23p")
