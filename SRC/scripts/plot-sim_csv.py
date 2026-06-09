# plot of simulation data only

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import rc
import matplotlib.lines as mlines
import pandas as pd
import csv
import math
import os

numFiles = 250 # number of Files to be loaded, for statistics
case = 'L500' # case name, e.g., L500, L500_OHx25%, L500_OHx50%, L500_OHx2, L500_OHx4

def set_lower_for_zero_offset_mm(ax, offset_mm=1.0):
    """
    Adjust only the lower y-limit so that y=0 is offset_mm above
    the bottom of the axes.

    The current upper y-limit is preserved.
    """

    # Draw first so Matplotlib knows the physical axes size
    ax.figure.canvas.draw()

    # Preserve current upper limit
    ymin_old, ymax_keep = ax.get_ylim()

    # Get axes height in mm
    bbox = ax.get_window_extent().transformed(
        ax.figure.dpi_scale_trans.inverted()
    )
    axes_height_mm = bbox.height * 25.4

    # Fractional position of y=0 above bottom of axes
    f = offset_mm / axes_height_mm

    if not (0 < f < 1):
        raise ValueError("offset_mm must be smaller than the axes height.")

    if ymax_keep <= 0:
        raise ValueError("Current upper y-limit must be positive.")

    # Solve:
    # (0 - ymin_new) / (ymax_keep - ymin_new) = f
    ymin_new = -f * ymax_keep / (1.0 - f)

    # Change lower limit only; keep upper limit fixed
    ax.set_ylim(ymin_new, ymax_keep)

    return ymin_new, ymax_keep


params = {'legend.fontsize': 15,
          'axes.labelsize': 17,
          'axes.titlesize': 17,
          'xtick.labelsize' :16,
          'axes.linewidth' : 2.0,
          'xtick.major.size': 5.5,
          'xtick.major.width': 2.0,
          'xtick.minor.width':3.5,
          'ytick.major.size': 5.5,
          'ytick.major.width': 2.0,
          'xtick.minor.width':3.5,
          'ytick.labelsize': 16,
          'grid.color': 'k',
          'grid.linestyle': ':',
          'grid.linewidth': 1.5,
          'mathtext.fontset' : 'stix',
          'mathtext.rm'      : 'serif',
          'font.family'      : 'serif',
          #'font.serif'       : 'Times New Roman', # or "Times"
          #'font.family'      : 'Times New Roman',
          #'font.weight' : "bold",
          #'axes.labelweight' :"bold",
         }
sizeFigH = 7.5
sizeFigV = 7
fig, ax1 = plt.subplots(figsize=(sizeFigH, sizeFigV))
matplotlib.rcParams.update(params)
plt.subplots_adjust(left=0.15, right=0.85, top=0.76, bottom=0.11)


OHlevel = ['OH x 4']
buildPath = f'../../RESULTS/{case}_results/'

norm_factor = []
outfig = f'{case}_PowdTech.png'
   

#=============================================================================================
# FIRST AXIS - particle volume fraction
#=============================================================================================

####### EXP vol. frac ##########
ax1.set_ylabel('Volume fraction / ($10^{-11}$)',labelpad=8,  fontsize=14, fontname='DejaVu Sans')
#########################################################################
################################

####### SIM vol. frac. #########
# NOW, Cantera density correction (optional at the end)
file_cantera = f"{buildPath}{case}-cantera.csv"
rho_cant = np.genfromtxt(file_cantera, delimiter=',', comments='#', usecols=5).T[1:]
grid_cant = np.genfromtxt(file_cantera, delimiter=',', comments='#', usecols=0).T[1:]
rho_norm = rho_cant / rho_cant[0]

# Seeds number
seeds_file = f'{buildPath}{case}-nuclest.csv'
seeds_nuclest = np.genfromtxt(seeds_file, delimiter=',', comments='#', usecols=1).T[:]

all_x=[]
all_curves=[]
all_T_curves=[]

stats_file = f"../../RESULTS/{case}_results_sens/{case}-compiled_curves.csv"
    
data = np.genfromtxt(stats_file, delimiter=",", names=True)
x_common = data["x_m"]

V_n      = data["V_n"]
V_mean   = data["V_mean"]
V_std    = data["V_std"]
V_median = data["V_median"]
V_q05    = data["V_q05"]
V_q95    = data["V_q95"]

T_n      = data["T_n"]
T_mean   = data["T_mean"]
T_std    = data["T_std"]
T_median = data["T_median"]
T_q05    = data["T_q05"]
T_q95    = data["T_q95"]

x_end_mean  = data["x_end_mean_m"][0]
x_end_std   = data["x_end_std_m"][0]

print(f"Loaded statistics from {stats_file}")

################################################################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111!

# Plot V SIM ============
ax1.plot(x_common*1e3, V_mean, '-', linewidth=2.0, color='brown', label='Vol-frac_MCsim (mean)')
ax1.fill_between(x_common*1e3, V_mean-V_std, V_mean+V_std, color='brown', alpha=0.2, label='MCsim - st. dev.')

# Mark mean final/end position
if np.isfinite(x_end_mean):


    x_mark_mm = x_end_mean * 1e3

    # x in data coordinates, y in axes-fraction coordinates
    trans = ax1.get_xaxis_transform()

    ax1.annotate(
        '',
        xy=(x_mark_mm, 0.01),        # arrow tip near x-axis
        xytext=(x_mark_mm, 0.18),    # arrow tail higher above x-axis
        xycoords=trans,
        textcoords=trans,
        arrowprops=dict(
            arrowstyle='-|>',
            color='black',
            lw=2.4,
            mutation_scale=20
        ),
        annotation_clip=False,
        zorder=30
    )

#=============================================================================================
# SECOND AXIS - temperature
#=============================================================================================
ax2 = ax1.twinx()

####### Gas temp SIM #########
# 2D sim data - extracted from 2D OF (T-of-x_L500.csv)
file_name = f"./2DOF_centerline/T_of_x-{case}.csv"
header = np.genfromtxt(file_name, delimiter=',',  dtype = None, encoding =None , comments = '#', max_rows = 1)
DFB_OF = np.genfromtxt(file_name, delimiter=',', comments='#', usecols=1).T
T_gas  = np.genfromtxt(file_name, delimiter=',', comments='#', usecols=0).T
ax2.set_ylabel('Temperature / ($10^{3}$ K)', labelpad = 6,  fontsize=14, fontname='DejaVu Sans')
ax2.plot(DFB_OF*1e3, T_gas, '-', linewidth = 1.5, color='y', label='Gas phase temp. - SIM')
##############################
ax1.get_yaxis().get_offset_text().set_visible(False)  # hide ×10^n label
ax2.get_yaxis().get_offset_text().set_visible(False)  # hide ×10^n label

# smooth part. temp line
ax2.plot(x_common*1e3, T_mean, linewidth=1.0, zorder = 3, color='slategrey', label='T_PART-MCsim')
#    ax2.fill_between(x_common*1e3, T_mean-T_std, T_mean+T_std, color='darkorchid', alpha=0.2, label='MCsim T - st. dev.')
##############################

####### EXP gas temp #########
########################################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
##############################

####### EXP Particle temp ####
##############################################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
##############################

#    plt.ticklabel_format(style= 'sci', axis= 'y', scilimits= (0,0))
ax1.ticklabel_format(style= 'sci', axis= 'y', scilimits= (0,0))
ax2.ticklabel_format(style= 'sci', axis= 'y', scilimits= (0,0))
ax2.set_ylim([-46,2800])
ax1.set_xlim([0,20])
ax1.tick_params(axis='y', labelsize=12)
ax2.tick_params(axis='y', labelsize=12)
# Then force y=0 to be exactly 1 mm above bottom of axes
set_lower_for_zero_offset_mm(ax1, offset_mm=1.0)
ax1.set_zorder(ax2.get_zorder()+1)
ax1.set_frame_on(False)

ax1.set_xlabel('DFB / mm', labelpad = 6, fontsize=14, fontname='DejaVu Sans')
ax1.tick_params(axis='x', labelsize=12)


c1 = mlines.Line2D([], [], linestyle = 'solid', color='y', linewidth=1.5,label='$T_{gas}$ - 2D sim.')
c2 = mlines.Line2D([], [], linestyle = '-', color='brown',linewidth=1.5, label='Part. vol. frac. - MC sim.')
c4 = mlines.Line2D([], [], linestyle = '-', color='lightslategray',linewidth=1.5, label='$T_{part}$ - MC sim.')

# add arrow in legend
arrow_legend = mlines.Line2D(
    [], [], 
    color='black', 
    marker=r'$\rightarrow$',   # Matplotlib mathtext arrow
    linestyle='None',
    markersize=12,
    label='Mean extinction position'
)

lns = [c1, arrow_legend, c4, c2]
labs = [l.get_label() for l in lns]
legend = plt.legend(lns, labs, loc="upper center", bbox_to_anchor=(-0.05,0.4,1.2,0.8),
                mode="expand", borderaxespad=0, ncol=2,columnspacing=0.12,frameon=False, labelspacing=0.2, handletextpad=0.5, prop={'family': 'DejaVu Sans', 'size': 14})
print(norm_factor)
plt.show()
fig.savefig(outfig,dpi=900)
fig.clf()
plt.close()