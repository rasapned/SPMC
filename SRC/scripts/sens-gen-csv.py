#processing file
import numpy as np
from scipy.signal import savgol_filter
from random import random
from matplotlib import rc
import os

count = 0 # to count the number of curves loaded, for statistics

x_end_all = [] #to store the vanishing place of each curve, for later statistics

#cases = ['L500_OHx25%','L500_OHx50%','L500','L500_OHx2','L500_OHx4']
case = 'L500'
buildPath = f'../../RESULTS/{case}_results/'
numFiles = 2 #250 # number of Files to be loaded, for statistics

# Seeds number
seeds_file = f'{buildPath}{case}-nuclest.csv' 
seeds_nuclest = np.genfromtxt(seeds_file, delimiter=',', comments='#', usecols=1).T[:]

all_x=[] #storing x grids of all curves
all_curves=[] #storing y values of all curves
all_T_curves=[] #storing temperature values of all curves

stats_file = f"../../RESULTS/{case}_results_sens/{case}-compiled_curves.csv"
#stats_file = f"{case}-compiled_curves.csv"

####### SIM vol. frac. #########
# NOW, Cantera density correction (optional at the end)
file_cantera = f'{buildPath}{case}-cantera.csv' 
rho_cant = np.genfromtxt(file_cantera, delimiter=',', comments='#', usecols=5).T[1:]
grid_cant = np.genfromtxt(file_cantera, delimiter=',', comments='#', usecols=0).T[1:]
rho_norm = rho_cant / rho_cant[0]

for i in range(0,numFiles):
    
    #=============================================================================================
    # FIRST AXIS - particle volume fraction
    #=============================================================================================
    
    ####### SIM vol. frac. #########
    base_num = i + 1  # e.g., 1, 15, 250
    possible_dirs = [
        f'../../RESULTS/{case}_results_sens/mech_{base_num:03d}',  # e.g., 001, 015, 250
        f'../../RESULTS/{case}_results_sens/mech_{base_num:02d}',  # e.g., 01, 15, 250
        f'../../RESULTS/{case}_results_sens/mech_{base_num:01d}',  # e.g., 1, 15, 250
    ]

    mc_file = None
    found = False
    for dir_path in possible_dirs:
        candidate_file = f'{dir_path}/{case}-average.csv'
        if os.path.isfile(candidate_file): # Check if the file exists
            mc_file = candidate_file
            found = True
            break

    if not found:
        print(f"Could not find {case}-average.csv for i={i+1} in any format!")
        continue

    # Proceed with mc_file
    print(f"Loaded: {mc_file}")            
    
    DFB_MC, Tp_MC, Dp_MC =  np.genfromtxt(mc_file, delimiter=',', comments='#', usecols=(1,2,3)).T
    count+=1
    # Normalisation of EXP V w.r.t. simulated V
    V_new = Dp_MC**3 * 3.1416 / 6.0 # * 1e-9 # convert from m^3 to mm^3, for better plotting
    for s,j in enumerate(DFB_MC):
        rho_int = np.interp(j,grid_cant,rho_norm)
        seeds_int = np.interp(j,grid_cant,seeds_nuclest)
        V_new[s] *=  seeds_int # * rho_int/rho_init 
    
    all_x.append(DFB_MC)
    all_curves.append(V_new)
    all_T_curves.append(Tp_MC)

    x_end_all.append(DFB_MC[-1])
    
#    ax1.plot(DFB_MC*1e3, V_new, '-', linewidth=0.1, alpha=0.2, color='cornflowerblue', label='Vol-frac_MCsim(diam-norm)')
#    ax1.plot(DFB_MC[sigmaTp_MC>0]*1e3, V_new[sigmaTp_MC>0], '-', linewidth=0.8, color=caseColor[n], label=OHlevel[n])
    ################################
    print(f"Loading file {i} successful!\n")

print(f"Current number of curves: {count}")
# Lengths of grids
lengths = [len(x) for x in all_x]
imax = np.argmax(lengths)

x_common = all_x[imax]   # preserve the longest grid
Ncases = len(all_curves)    # no.  of cases
Nx = len(x_common)          # grid length

# Create one big common matrix, fill with nan
curves = np.full((Ncases, Nx), np.nan)
curves_T = np.full((Ncases, Nx), np.nan)

# Start filling with the actual profiles
for j, (x, y, z) in enumerate(zip(all_x, all_curves, all_T_curves)):
    npts = len(x)

#            # optional safety check: x must match the beginning of x_common
#            if not np.allclose(x, x_common[:npts]):
#                raise ValueError(f"Grid mismatch in case {j}: grid is not a prefix of reference grid")

    curves[j, :npts] = y
    curves_T[j, :npts] = z


def calc_nan_stats(arr):
    n_valid = np.sum(~np.isnan(arr), axis=0)

    mean = np.nanmean(arr, axis=0)
    median = np.nanmedian(arr, axis=0)
    q05 = np.nanpercentile(arr, 5, axis=0)
    q95 = np.nanpercentile(arr, 95, axis=0)

    std = np.full(arr.shape[1], np.nan)
    valid_std = n_valid > 1
    std[valid_std] = np.nanstd(arr[:, valid_std], axis=0, ddof=1)

    return n_valid, mean, std, median, q05, q95

V_n, V_mean, V_std, V_median, V_q05, V_q95 = calc_nan_stats(curves)
T_n, T_mean, T_std, T_median, T_q05, T_q95 = calc_nan_stats(curves_T)

# Mean vanishing place
x_end_all = np.asarray(x_end_all)

x_end_mean = np.mean(x_end_all)
x_end_std = np.std(x_end_all, ddof=1)
Nx = len(x_common)

# =============================================================================================
# Save only statistics
# =============================================================================================

out_stats = np.column_stack(
    (
        x_common,
        V_n,
        V_mean,
        V_std,
        V_median,
        V_q05,
        V_q95,
        T_n,
        T_mean,
        T_std,
        T_median,
        T_q05,
        T_q95,
        np.full(Nx, x_end_mean),
        np.full(Nx, x_end_std),
    )
)

header = ",".join(
    [
        "x_m",
        "V_n",
        "V_mean",
        "V_std",
        "V_median",
        "V_q05",
        "V_q95",
        "T_n",
        "T_mean",
        "T_std",
        "T_median",
        "T_q05",
        "T_q95",
        "x_end_mean_m",
        "x_end_std_m",
    ]
)

np.savetxt(
    stats_file,
    out_stats,
    delimiter=",",
    header=header,
    comments="",
    fmt="%.10e"
)

print(f"Saved statistics to {stats_file}")