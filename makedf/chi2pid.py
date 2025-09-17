from . import calo
import pandas as pd
import numpy as np
import sqlite3
import uproot

def chi2(hitdf, exprr, expdedx, experr, dedxname="dedx"):
    dedx_exp = pd.cut(hitdf.rr, exprr, labels=expdedx).astype(float)
    dedx_err = pd.cut(hitdf.rr, exprr, labels=experr).astype(float)

    dedx_res = (0.04231 + 0.0001783*hitdf[dedxname]**2)*hitdf[dedxname]

    v_chi2 = (hitdf[dedxname] - dedx_exp)**2 / (dedx_err**2 + dedx_res**2)

    when_chi2 = (hitdf.rr < np.max(exprr)) & ~hitdf.firsthit & ~hitdf.lasthit & (hitdf[dedxname] < 1000.)

    chi2_group = v_chi2[when_chi2].groupby(level=list(range(hitdf.index.nlevels-1)))

    return chi2_group.sum() / chi2_group.size()

def chi2u(hitdf, dedxname="dedx"):
    return chi2(hitdf, muon_rr, muon_dedx, muon_yerr, dedxname)

def chi2p(hitdf, dedxname="dedx"):
    return chi2(hitdf, proton_rr, proton_dedx, proton_yerr, dedxname)

def chi2_ndof(hitdf):
    when_chi2 = (hitdf.rr < np.max(exprr)) & ~hitdf.firsthit & ~hitdf.lasthit & (hitdf.dedx < 1000.)
    chi2_group = v_chi2[when_chi2].groupby(level=list(range(hitdf.index.nlevels-1)))

    return chi2_group.size()

def dqdx(dqdxdf, gain=None, calibrate=None):
    if calibrate == "ICARUS": 
        # get raw dqdx
        dqdx = dqdxdf.integral / dqdxdf.pitch

        # compute y-scale
        ybin = _yz_ybin(dqdxdf.y)
        zbin = _yz_zbin(dqdxdf.z)
        iov = _yz_iov(dqdxdf.run)
        itpc = dqdxdf.tpc // 2 + dqdxdf.cryo*2
        plane = dqdxdf.plane

        yzdf = pd.DataFrame({"ybin": ybin, "zbin": zbin, "itpc": itpc, "plane": plane, "iov": iov})
        yz_scale = yzdf.merge(IC_yz_cal_df, on=["iov", "itpc", "plane", "ybin", "zbin"], how="left", validate="many_to_one").scale
        yz_scale[yz_scale == -999.000000] = 1.
        yz_scale = np.clip(yz_scale, 0.7, 1.3).fillna(1)
        yz_scale.index = dqdxdf.index

        # compute lifetime correction
        iov = _etau_iov(dqdxdf.run)
        etaudf = pd.DataFrame({"itpc": itpc, "iov": iov})
        etau = etaudf.merge(IC_etau_cal_df, on=["iov", "itpc"], how="left", validate="many_to_one").etau
        etau = etau.fillna(np.inf)
        etau.index = dqdxdf.index

        # compute TPC scale
        iov = _tpc_iov(dqdxdf.run)
        tpcdf = pd.DataFrame({"itpc": itpc, "plane": plane, "iov": iov})
        tpc_scale = tpcdf.merge(IC_tpc_cal_df, on=["iov", "itpc", "plane"], how="left", validate="many_to_one").scale
        tpc_scale.index = dqdxdf.index

        # apply the corrections
        t0 = 0 # assume in time
        tick_period = 0.4 # us
        tanode = 850 # ticks
        tdrift = dqdxdf.t*tick_period - t0 - tanode*tick_period

        dqdx = dqdx * tpc_scale * np.exp(tdrift / etau) / yz_scale
    elif calibrate == "SBND": # TODO: add calibrations?
        dqdx  = dqdxdf.integral / dqdxdf.pitch
    else: # if not specified, rely on input calibration
        dqdx = dqdxdf.dqdx

    # apply gain
    if gain == "ICARUS":
        gains = [0.016751, 0.012755, 0.012513]
        gain_perhit = pd.Series(1.0, dqdxdf.index)
        for iplane, g in enumerate(gains):
            gain_perhit[dqdxdf.plane == iplane] = 1.0/g
    elif gain == "SBND": # TODO
        gain_perhit = 1
    else:
        gain_perhit = 1

    return dqdx*gain_perhit

def dedx(dqdxdf, gain=None, calibrate=None, scalegain=1):
    dqdx_v = dqdx(dqdxdf, gain=gain, calibrate=calibrate)

    return calo.recombination_cor(dqdx_v*scalegain, dqdxdf.phi, dqdxdf.efield, dqdxdf.rho)

def _yz_ybin(y):
    return np.searchsorted(yz_ybin, y) - 1

def _yz_zbin(z):
    return np.searchsorted(yz_zbin, z) - 1

def _yz_iov(run): 
    return __iov(run, IC_yz_cal_iovdf)

def _etau_iov(run):
    return __iov(run, IC_etau_cal_iovdf)

def _tpc_iov(run):
    return __iov(run, IC_tpc_cal_iovdf)

def __iov(run, df):
    return pd.cut(run, list(df.run) + [np.inf], labels=df.iov)

# EXPECTED dE/dx FILES
datadir = "/exp/icarus/data/users/gputnam/"
fhist = datadir + "dEdxrestemplates.root"

profp = uproot.open(fhist)["dedx_range_pro"]
profmu = uproot.open(fhist)["dedx_range_mu"]

proton_dedx = profp.values()
proton_rr = profp.axis().edges()
proton_yerr = profp.errors(error_mode="s")
for i in range(len(proton_yerr)):
    if proton_yerr[i] < 1e-6:
        proton_yerr[i] = (proton_yerr[i-1] + proton_yerr[i+1]) / 2
    if proton_dedx[i] < 1e-6:
        proton_dedx[i] = (proton_dedx[i-1] + proton_dedx[i+1]) / 2

muon_dedx = profmu.values()
muon_rr = profmu.axis().edges()
muon_rr_center = profmu.axis().centers()
muon_yerr = profmu.errors(error_mode="s")

# ICARUS CALIBRATION DATABASES
IC_yz_cal_f = "/cvmfs/icarus.opensciencegrid.org/products/icarus/icarus_data/v10_06_00/icarus_data/database/tpc_yz_correction_allplanes_data.db"
IC_yz_cal_db = "tpc_yz_correction_allplanes_data_data"
IC_yz_cal_iov = "tpc_yz_correction_allplanes_data_iovs"

IC_etau_cal_f = "/cvmfs/icarus.opensciencegrid.org/products/icarus/icarus_data/v10_06_00/icarus_data/database/tpc_elifetime_data.db"
IC_etau_cal_db = "tpc_elifetime_data_data"
IC_etau_cal_iov = "tpc_elifetime_data_iovs"

IC_tpc_cal_f = "/cvmfs/icarus.opensciencegrid.org/products/icarus/icarus_data/v10_06_00/icarus_data/database/tpc_dqdxcalibration_allplanes_data.db"
IC_tpc_cal_db = "tpc_dqdxcalibration_allplanes_data_data"
IC_tpc_cal_iov = "tpc_dqdxcalibration_allplanes_data_iovs"

# LOAD THE YZ CALIBRATION
yz_ybin = np.linspace(-180, 130, 32)
yz_ylos = yz_ybin[:-1]
yz_yhis = yz_ybin[1:]
yz_ys = (yz_ylos + yz_yhis) / 2.

yz_zbin = np.linspace(-900, 900, 181)
yz_zlos = yz_zbin[:-1]
yz_zhis = yz_zbin[1:]
yz_zs = (yz_zlos + yz_zhis) / 2.

conn = sqlite3.connect(IC_yz_cal_f)
cursor = conn.cursor()
cursor.execute("SELECT * FROM %s" % IC_yz_cal_db)
rows = cursor.fetchall()
data = list(zip(*rows))
IC_yz_cal_df = pd.DataFrame({
  "iov": data[0],
  "plane": data[2],
  "tpc": data[3],
  "ybin": data[4],
  "zbin": data[5],
  "scale": data[6]
})
IC_yz_cal_df["itpc"] = 0 
IC_yz_cal_df.loc[IC_yz_cal_df.tpc == "EE", "itpc"] = 0
IC_yz_cal_df.loc[IC_yz_cal_df.tpc == "EW", "itpc"] = 1
IC_yz_cal_df.loc[IC_yz_cal_df.tpc == "WE", "itpc"] = 2
IC_yz_cal_df.loc[IC_yz_cal_df.tpc == "WW", "itpc"] = 3
del IC_yz_cal_df["tpc"]

cursor.execute("SELECT * FROM %s WHERE ACTIVE=1" % IC_yz_cal_iov)
rows = cursor.fetchall()
data = list(zip(*rows))
IC_yz_cal_iovdf = pd.DataFrame({
  "iov": data[0],
  "begin_time": data[1],
})
IC_yz_cal_iovdf["run"] = IC_yz_cal_iovdf.begin_time % 1000000000
IC_yz_cal_iovdf.sort_values(by="run", inplace=True)
conn.close()

# LOAD THE LIFETIME CALIBRATION
conn = sqlite3.connect(IC_etau_cal_f)
cursor = conn.cursor()
cursor.execute("SELECT * FROM %s" % IC_etau_cal_db)
rows = cursor.fetchall()
data = list(zip(*rows))
IC_etau_cal_df = pd.DataFrame({
  "iov": data[0],
  "itpc": data[1],
  "etau": data[2]
})

cursor.execute("SELECT * FROM %s WHERE ACTIVE=1" % IC_etau_cal_iov)
rows = cursor.fetchall()
data = list(zip(*rows))
IC_etau_cal_iovdf = pd.DataFrame({
  "iov": data[0],
  "begin_time": data[1],
})
IC_etau_cal_iovdf["run"] = IC_etau_cal_iovdf.begin_time % 1000000000
IC_etau_cal_iovdf.sort_values(by="run", inplace=True)
conn.close()

# LOAD THE TPC SCALE
conn = sqlite3.connect(IC_tpc_cal_f)
cursor = conn.cursor()
cursor.execute("SELECT * FROM %s" % IC_tpc_cal_db)
rows = cursor.fetchall()
data = list(zip(*rows))
IC_tpc_cal_df = pd.DataFrame({
  "iov": data[0],
  "plane": data[2],
  "tpc": data[3],
  "scale": data[4]
})
IC_tpc_cal_df["itpc"] = 0 
IC_tpc_cal_df.loc[IC_tpc_cal_df.tpc == "EE", "itpc"] = 0
IC_tpc_cal_df.loc[IC_tpc_cal_df.tpc == "EW", "itpc"] = 1
IC_tpc_cal_df.loc[IC_tpc_cal_df.tpc == "WE", "itpc"] = 2
IC_tpc_cal_df.loc[IC_tpc_cal_df.tpc == "WW", "itpc"] = 3
del IC_tpc_cal_df["tpc"]

cursor.execute("SELECT * FROM %s WHERE ACTIVE=1" % IC_tpc_cal_iov)
rows = cursor.fetchall()
data = list(zip(*rows))
IC_tpc_cal_iovdf = pd.DataFrame({
  "iov": data[0],
  "begin_time": data[1],
})
IC_tpc_cal_iovdf["run"] = IC_tpc_cal_iovdf.begin_time % 1000000000
IC_tpc_cal_iovdf.sort_values(by="run", inplace=True)
conn.close()

