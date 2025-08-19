# Standard library imports
import sys
import os

# Third-party imports
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Add the head direcoty to sys.path
workspace_root = os.getcwd()
sys.path.insert(0, workspace_root + "/../../")
# Local imports
import kinematics
from makedf.util import *
from sel_tools import *

class PlotObj:
  def __init__(self, title, label):
    self.title = title
    self.label = label

# Colors for plots
HAWKS_COLORS = ["#315031", "#d54c28", "#1e3f54", "#c89648", "#43140b", "#95af8b"]

top_labels = ["Signal",
              "Other numu CC",
              "NC",
              "Out of FV",
              "Cosmic",
              "Other"]

FONTSIZE = 14

def plot_int(df_nd, df_fd, var, title, outfile, label, eff_bool=False):
    var_nd = df_nd[var]
    var_fd = df_fd[var]

    og_nd_sig_ct = df_nd['og_sig_ct'].to_numpy()[0]
    og_fd_sig_ct = df_fd['og_sig_ct'].to_numpy()[0]
    glob_nd_scale = df_nd['glob_scale'].to_numpy()[0]
    glob_fd_scale = df_fd['glob_scale'].to_numpy()[0]

    plt.tight_layout()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))
    pvar_nd = breakdown_mode(var_nd, df_nd)
    pvar_fd = breakdown_mode(var_fd, df_fd)
    ax1.set_title(title+", SBND Only", fontsize=FONTSIZE)
    ax1.set_title("SBND", fontsize=FONTSIZE)
    n, bins, _ = ax1.hist(pvar_nd, bins=np.linspace(0,2,21), stacked=True, label=mode_labels, 
                          color=HAWKS_COLORS, weights=[glob_nd_scale*np.ones_like(p) for p in pvar_nd])

    tot = sum(n[-1])
    nd_ret = [sum(n[0])/tot]
    for i in range(1, len(n)):
        nd_ret.append((sum(n[i])-sum(n[i-1]))/tot)

    nd_sig_ct = len(var_nd[is_signal(df_nd, "SBND")])
    if eff_bool:
        ax1.text(0.5, 0.33, "ND Purity {:.2f}%".format(100*nd_sig_ct/len(df_nd)), 
                 transform=ax1.transAxes, fontsize=FONTSIZE)
        ax1.text(0.5, 0.4, "ND Efficiency {:.2f}%".format(100*nd_sig_ct/og_nd_sig_ct), 
                 transform=ax1.transAxes, fontsize=FONTSIZE)
    ax1.xaxis.set_tick_params(labelsize=FONTSIZE)
    ax1.yaxis.set_tick_params(labelsize=FONTSIZE)
    ax1.set_xlabel(label, fontsize=FONTSIZE)

    ax2.set_title("ICARUS", fontsize=FONTSIZE)
    n, bins, _ = ax2.hist(pvar_fd, bins=np.linspace(0,2,21), stacked=True, label=mode_labels, 
                          color=HAWKS_COLORS, weights=[glob_fd_scale*np.ones_like(p) for p in pvar_fd])
    fd_sig_ct = len(var_fd[is_signal(df_fd, "ICARUS")])
    if eff_bool:
        ax2.text(0.5, 0.33, "FD Purity {:.2f}%".format(100*fd_sig_ct/len(df_fd)), 
                 transform=ax2.transAxes, fontsize=FONTSIZE)
        ax2.text(0.5, 0.4, "FD Efficiency {:.2f}%".format(100*fd_sig_ct/og_fd_sig_ct), 
                 transform=ax2.transAxes, fontsize=FONTSIZE)
    ax2.xaxis.set_tick_params(labelsize=FONTSIZE)
    ax2.yaxis.set_tick_params(labelsize=FONTSIZE)
    ax2.set_xlabel(label, fontsize=FONTSIZE)

    fig.suptitle(title, fontsize=FONTSIZE)
    plt.legend(fontsize=FONTSIZE)
    plt.savefig("EventSelectionPlots/"+outfile, bbox_inches='tight')
    plt.clf()
    plt.close()

    tot = sum(n[-1])
    fd_ret = [sum(n[0])/tot]
    for i in range(1, len(n)):
        fd_ret.append((sum(n[i])-sum(n[i-1]))/tot)

    return [nd_ret, fd_ret]

def plot_fs(df_nd, df_fd, var, title, outfile):
    var_nd = df_nd[var]
    var_fd = df_fd[var]

    og_nd_sig_ct = df_nd['og_sig_ct'].to_numpy()[0]
    glob_nd_scale = df_nd['glob_scale'].to_numpy()[0]
    og_fd_sig_ct = df_fd['og_sig_ct'].to_numpy()[0]
    glob_fd_scale = df_fd['glob_scale'].to_numpy()[0]

    plt.tight_layout()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))
    prim_v = ['np', 'nn', 'nmu', 'npi', 'npi0', 'ng']
    prim_v_labels = ['Protons', 'Neutrons', 'Muons', 'Charged Pions', 'Neutral Pions', 'Gammas']
    b = np.arange(1, 11, 1)
    bin_centers = 0.5 * (b[:-1] + b[1:])
    bin_numbers = range(1, len(bin_centers)+1)

    df_fs = [df_nd[p] for p in prim_v]
    ax1.set_xticks(bin_centers, [str(i) for i in bin_numbers])
    ax1.hist(df_fs, bins=b, stacked=True, label=prim_v_labels, color=HAWKS_COLORS, 
             weights=[glob_nd_scale*np.ones_like(p) for p in df_fs])
    ax1.xaxis.set_tick_params(labelsize=FONTSIZE)
    ax1.yaxis.set_tick_params(labelsize=FONTSIZE)

    df_fs = [df_fd[p] for p in prim_v]
    ax2.set_xticks(bin_centers, [str(i) for i in bin_numbers])
    ax2.hist(df_fs, bins=b, stacked=True, label=prim_v_labels, color=HAWKS_COLORS, 
             weights=[glob_fd_scale*np.ones_like(p) for p in df_fs])
    ax2.xaxis.set_tick_params(labelsize=FONTSIZE)
    ax2.yaxis.set_tick_params(labelsize=FONTSIZE)

    ax1.set_xlabel("Final State Particle Count")
    ax1.set_title("SBND")
    ax2.set_xlabel("Final State Particle Count")
    ax2.set_title("ICARUS", fontsize=FONTSIZE)
    fig.suptitle(title, fontsize=FONTSIZE)
    plt.legend(fontsize=FONTSIZE)
    plt.savefig("EventSelectionPlots/"+outfile, bbox_inches='tight')
    fig.clf()
    plt.close()

def plot_PID_cut(var, cut_vals, title, outfile, xlims=[0, 100], 
                 labels=['SBND', 'ICARUS'], arrow_dir='None', arrow_txt='None'):
    
    var_nd, var_fd = var
    cut_val_nd, cut_val_fd = cut_vals
    b = np.linspace(xlims[0], xlims[1], 40)
    plt.tight_layout()

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 4.5), sharex=True)
    plt.subplots_adjust(hspace=0)

    n1, bins, _ = ax1.hist(var_nd[0], bins=b, histtype='step', color=HAWKS_COLORS[0], 
                           label=labels[0], stacked=False, weights=[1/len(var_nd[0])]*len(var_nd[0]), linestyle='-')
    n2, bins, _ = ax1.hist(var_nd[1], bins=b, histtype='step', color=HAWKS_COLORS[0], 
                           label=labels[1], stacked=False, weights=[1/len(var_nd[1])]*len(var_nd[1]), linestyle='--')
    n3, bins, _ = ax2.hist(var_fd[0], bins=b, histtype='step', color=HAWKS_COLORS[1], 
                           label=labels[2], stacked=False, weights=[1/len(var_fd[0])]*len(var_fd[0]), linestyle='-')
    n4, bins, _ = ax2.hist(var_fd[1], bins=b, histtype='step', color=HAWKS_COLORS[1], 
                           label=labels[3], stacked=False, weights=[1/len(var_fd[1])]*len(var_fd[1]), linestyle='--')

    ax1.legend(fontsize=FONTSIZE)
    ax2.legend(fontsize=FONTSIZE)
    ax2.xaxis.set_tick_params(labelsize=FONTSIZE)
    ax2.yaxis.set_tick_params(labelsize=FONTSIZE)
    max_1 = np.max([n1, n2])*1.2
    max_2 = np.max([n3, n4])*1.2
    ax1.vlines(cut_val_nd, 0, max_1, colors='black', linestyle='--')
    ax2.vlines(cut_val_fd, 0, max_2, colors='black', linestyle='--')
    ax1.set_ylim(0, max_1)
    ax2.set_ylim(0, max_2)

    if arrow_dir == 'right':
       ax1.arrow(cut_val_nd, max_1/2, (xlims[1]-xlims[0])/8, 0, head_width=0.1*max_1, 
                 head_length=(xlims[1]-xlims[0])/40, fc='red', ec='red')
       ax2.arrow(cut_val_fd, max_2/2, (xlims[1]-xlims[0])/8, 0, head_width=0.1*max_2, 
                 head_length=(xlims[1]-xlims[0])/40, fc='red', ec='red')
       ax1.text(cut_val_nd+10, 0.35*max_1, arrow_txt, color='red')
       ax2.text(cut_val_fd+10, 0.35*max_2, arrow_txt, color='red')
    elif arrow_dir == 'left':           
       ax1.arrow(cut_val_nd, max_1/2, -(xlims[1]-xlims[0])/8, 0, head_width=0.1*max_1, 
                 head_length=(xlims[1]-xlims[0])/40, fc='red', ec='red')
       ax2.arrow(cut_val_fd, max_2/2, -(xlims[1]-xlims[0])/8, 0, head_width=0.1*max_2, 
                 head_length=(xlims[1]-xlims[0])/40, fc='red', ec='red')
       ax1.text(cut_val_nd-20, 0.8*max_1, arrow_txt, color='red')
       ax2.text(cut_val_fd-20, 0.8*max_2, arrow_txt, color='red')

    ax1.xaxis.set_tick_params(labelsize=FONTSIZE)
    ax1.yaxis.set_tick_params(labelsize=FONTSIZE)
    fig.suptitle(title, fontsize=FONTSIZE)

    plt.savefig("EventSelectionPlots/"+outfile, bbox_inches='tight')
    plt.clf()
    plt.close()

def plot_nuscore_cut(var_nd, cut_vals_nd, var_fd, cut_vals_fd, title, outfile, xlims=[0, 100], 
                     labels=['SBND', 'ICARUS'], arrow_dir='None', arrow_txt='None'):
    
    plt.tight_layout()
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 4.5), sharex=True)
    plt.subplots_adjust(hspace=0)
    
    b = np.linspace(xlims[0], xlims[1], 40)
    n1, bins, _ = ax1.hist(var_nd[0], bins=b, histtype='step', color=HAWKS_COLORS[0], 
                           label=labels[0], stacked=False, weights=[1/len(var_nd[0])]*len(var_nd[0]), linestyle='-')
    n2, bins, _ = ax1.hist(var_nd[1], bins=b, histtype='step', color=HAWKS_COLORS[0], 
                           label=labels[1], stacked=False, weights=[1/len(var_nd[1])]*len(var_nd[1]), linestyle='--')
    # n3, bins, _ = ax2.hist(var_fd[0], bins=b, histtype='step', color=HAWKS_COLORS[1], 
    #                        label=labels[2], stacked=False, weights=[1/len(var_fd[0])]*len(var_fd[0]), linestyle='-')
    n4, bins, _ = ax2.hist(var_fd[1], bins=b, histtype='step', color=HAWKS_COLORS[1], 
                           label=labels[3], stacked=False, weights=[1/len(var_fd[1])]*len(var_fd[1]), linestyle='--')

    max_1 = np.max([n1, n2])*3.0
    max_2 = np.max(n4)*3.0
    # max_2 = np.max([n3, n4])*1.2
    ax1.vlines(cut_vals_nd, 0, max_1, colors='black', linestyle='--')
    ax2.vlines(cut_vals_fd, 0, max_2, colors='black', linestyle='--')
    ax1.set_ylim(0, max_1)
    ax2.set_ylim(0, max_2)

    ax1.arrow(cut_vals_nd[0], max_1/2, (xlims[1]-xlims[0])/8, 0, head_width=0.1*max_1, 
              head_length=(xlims[1]-xlims[0])/40, fc='red', ec='red')
    ax2.arrow(cut_vals_fd[0], max_2/2, (xlims[1]-xlims[0])/8, 0, head_width=0.1*max_2, 
              head_length=(xlims[1]-xlims[0])/40, fc='red', ec='red')
    ax1.text(cut_vals_nd[0]+(xlims[1]-xlims[0])/20, 0.35*max_1, arrow_txt, color='red')
    ax2.text(cut_vals_fd[0]+(xlims[1]-xlims[0])/20, 0.35*max_2, arrow_txt, color='red')

    fig.suptitle(title, fontsize=FONTSIZE)

    for a in [ax1, ax2]:
        a.legend(fontsize=FONTSIZE)
        a.xaxis.set_tick_params(labelsize=FONTSIZE)
        a.yaxis.set_tick_params(labelsize=FONTSIZE)

    plt.savefig("EventSelectionPlots/"+outfile, bbox_inches='tight')
    plt.clf()
    plt.close()

def plot_cut(var_nd, cut_vals_nd, var_fd, cut_vals_fd, title, outfile, xlims=[0, 100], 
             labels=['SBND', 'ICARUS'], arrow_dir='None', arrow_txt='None'):
    b = np.linspace(xlims[0], xlims[1], 40)
    plt.tight_layout()
    fig, ax1 = plt.subplots(1, 1, figsize=(6, 4.5))
    n, bins, _ = ax1.hist(var_nd, bins=b, histtype='step', color=HAWKS_COLORS[0], label='SBND')
    n, bins, _ = ax1.hist(var_fd, bins=b, histtype='step', color=HAWKS_COLORS[1], label='ICARUS')
    ax1.vlines(cut_vals_fd, 0, np.max(n), colors='black', linestyle='--')
    ax1.xaxis.set_tick_params(labelsize=FONTSIZE)
    ax1.yaxis.set_tick_params(labelsize=FONTSIZE)
    fig.suptitle(title, fontsize=FONTSIZE)

    plt.savefig("EventSelectionPlots/"+outfile, bbox_inches='tight')
    plt.clf()
    plt.close()

def plot_stub_2d(df, cuts, outfile, title='test'):
    plt.figure(figsize=(8, 6), dpi=80)
    dQdx = []
    x = []
    for l in ["0_5", "1", "2", "3"]:
        x.extend(df.stub["l"+l+"cm"].length)
        dQdx.extend(df.stub["l"+l+"cm"].Q/df.stub["l"+l+"cm"].length)

    plt.hlines(cuts, [0, 0.5, 1.0, 2.0], [0.5, 1.0, 2.0, 3.0], color='red', lw=1.0, linestyle='--')
    plt.vlines([0.5, 1.0, 2.0], cuts[1:], cuts[:-1], color='red', lw=1.0, linestyle='--')

    plt.text(0.5, 2e5, 'Select', color='red', fontsize=FONTSIZE*2)
    plt.text(2.0, 6e5, 'Remove', color='red', fontsize=FONTSIZE*2)

    plt.title(title, fontsize=FONTSIZE)
    plt.tick_params(labelsize=FONTSIZE)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.xlabel('Stub Length [cm]', fontsize=FONTSIZE)
    plt.ylabel('dQ/dx [electrons/cm]', fontsize=FONTSIZE)
    plt.hist2d(x, dQdx, bins=[6, 16], range=[[0, 3],[0, 8e5]])
    plt.colorbar()
    plt.savefig("EventSelectionPlots/"+outfile)
    plt.clf()
    plt.close()

def plot_composition(percentages, time_labels=None, components=None, title='', outfile='comp.png'):
    percentages = np.flip(np.array(percentages), axis=0)
    num_time_points, num_components = percentages.shape

    if components is None:
        components = [f"Component {i+1}" for i in range(num_components)]
    if time_labels is None:
        time_labels = [f"Time {i+1}" for i in range(num_time_points)]

    fig, ax = plt.subplots(figsize=(7, 6))
    y = np.arange(num_time_points)

    left = np.zeros(num_time_points)
    for i in range(num_components):
        ax.barh(y, percentages[:, i], left=left, label=components[i])
        left += percentages[:, i]

    ax.set_yticks(y)
    ax.set_yticklabels(reversed(time_labels), fontsize=FONTSIZE)
    plt.tick_params(labelsize=FONTSIZE)
    ax.set_xlabel('Percentage', fontsize=FONTSIZE)
    ax.set_xlim(0, 1)
    ax.legend(fontsize=FONTSIZE)
    ax.set_title('Component Composition Over Time, '+title, fontsize=FONTSIZE)
    plt.tight_layout()
    plt.savefig("EventSelectionPlots/"+outfile)
