import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.gridspec as gridspec
from matplotlib import rc
import matplotlib.lines as mlines
import random
import csv
from scipy.optimize import minimize
from scipy.optimize import basinhopping
from scipy.optimize import differential_evolution
from matplotlib.ticker import MaxNLocator
import os

#### read the data from the Monte Carlo simulations and calculate the reaction constants for each species ####

params = {'legend.fontsize': 16,
          'axes.labelsize': 18,
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
          'grid.linewidth': 0.7,
          'mathtext.fontset' : 'stix',
          'mathtext.rm'      : 'serif',
          'font.family'      : 'serif',
         }

fig2 = plt.figure(figsize=[4.2,6.0])
gs2 = gridspec.GridSpec(3, 2)
gs2.update(wspace=0.08, hspace=0.3) # set the spacing between axes. 
plt.subplots_adjust(left=0.15, right=0.98, wspace=0.1, hspace=0.25, top=0.9, bottom=0.1)

# Define some general constants
N_A = 6.02214076 * 1e23 # 1/mol
PI = 3.1415926
R_g = 8.31446           # J/mol/K
k_B = 1.380649 * 1e-23  # J/K
T_ref = 298          # K
#T_ref = 3000          # K

# define which cases should be read and used for optimisation
cases = ['L500','S500','R500']#,'XL500']
n_cases = len(cases)
phicol = ['m','c','y']#,'sandybrown']

# Define some flame constants
pressure = 3000         # Pa

# Molar masses of species
M_FE = 0.055845
M_H2 = 0.002016
M_O2 = 0.032000
M_O  = 0.016000
M_H2O = 0.01801528
M_H = 0.001008
M_OH = 0.017008    

# Mean free paths * press
# diameter values for the mean free path calculation
d_FE   = 2.52e-10
d_O2   = 2.92e-10
d_H2O  = 2.75e-10
d_H2   = 2.89e-10
d_O   = d_O2/2
d_H   = d_H2/2
d_OH  = 0.97e-10

#### User defined parameters ####
# Include/Exclude temperature steric factor b
b_switch = False
mk_switch = False   
alpha_switch = False
Ek_switch = False

# Release the activation energy fitting
blockE_a = False

# temperature sampling for the final reaction constant curve
T_bf = np.arange(850,1900,5)

############################################

# initialise arrays that store MC read data for all cases (different lengths)
Tp_glob, cov_OFe, k_Fe,k_O2,k_O,k_H2O,k_H2,k_H,k_OH = np.zeros((9,0),dtype = float)

# initialise bestfit array for the coverage (O/Fe ratio)
rat_bf = np.zeros((n_cases,len(T_bf)),dtype=float)

# output figure name
outfig = 'Arrhenius.png'

# assign graphs to axes 
aax2 = plt.subplot(gs2[0,0]) 
aax3 = plt.subplot(gs2[0,1]) 
aax4 = plt.subplot(gs2[1,0]) 
aax5 = plt.subplot(gs2[1,1]) 
aax6 = plt.subplot(gs2[2,0]) 
aax7 = plt.subplot(gs2[2,1]) 

# read all the necessary data for each flame case and get the reaction constants
for n, name in enumerate(cases):
    
    #### READ PART ####

    mc_file = f"../../RESULTS/L500_results/{name}-average.csv" # contains the average values of the Monte Carlo simulations (DFB, particle diameter, particle temperature, temperature before reaction)
    mc_count_file = f"../../RESULTS/L500_results/{name}-count.csv" # contains the collision counts for each species and the total number of collisions (for reaction constant calculation)

    # Cantera 1D sim: read the grid, temperature, density and mole fractions of the gaseous iron containing species
    #sim1D_XFEC5O5, sim1D_XFE2O3, sim1D_XFEO2, sim1D_XFEO, sim1D_XFEOH, sim1D_XFEO2H2, sim1D_XFE2OOOH = np.genfromtxt(name+'-cantera.csv', delimiter=',', comments='#', skip_header=1, usecols=(17,36,10,7,8,9,34)).T
    sim1D_XFEC5O5, sim1D_XFE2O3, sim1D_XFEO2, sim1D_XFEO, sim1D_XFEOH, sim1D_XFEO2H2, sim1D_XFE2OOOH = np.genfromtxt(f"../../RESULTS/L500_results/{name}-cantera.csv", delimiter=',', comments='#', skip_header=1, usecols=(17,36,10,7,8,9,34)).T
    #sim1D_z, sim1D_Tg, sim1D_rho = np.genfromtxt(f"../../RESULTS/{name}_results/{name}-cantera.csv", delimiter=',', comments='#',skip_header=1, usecols=(0,4,5)).T
    sim1D_z, sim1D_Tg, sim1D_rho = np.genfromtxt(f"../../RESULTS/L500_results/{name}-cantera.csv", delimiter=',', comments='#',skip_header=1, usecols=(0,4,5)).T
    # Cantera concentrations of species of interest
    c_O2, c_O, c_H2O, c_H2, c_H, c_OH, c_FEC5O5, c_FE2O3, c_FEO2, c_FEO, c_FEO2H2,c_FE2OOOH = np.genfromtxt(f"../../RESULTS/L500_results/{name}-ConcCant.csv", delimiter=',', comments='#', skip_header=1).T

    # MC-Results
    DFB, Tp, Dp, Tbf = np.genfromtxt(mc_file, delimiter=',', comments='#', usecols=(1,2,3,6)).T # DFB - distance from burner, Tp - particle temperature, Dp - particle diameter, Tbf - temperature before reaction (after thermalisation)
    
    # Reactive collision rates from Monte Carlo
    C_Z, C_Fe, C_O2, C_O, C_H2O, C_H2, C_H, C_OH, C_EV, C_tot, rat = np.genfromtxt(mc_count_file, delimiter=',', comments='#', usecols=(0,1,2,3,4,5,6,7,8,9,10)).T
    
    #### REACTION CONSTANTS CALCULATION PART ####

    # get particle surface area, coverage (O/Fe), particle and gas temperature at all C_Z (collision data grid)
    p_diam = np.interp(C_Z, DFB, Dp) #Unit of Dp is m, so p_diam is in m
    p_Temp = np.interp(C_Z, DFB, Tbf)           # Tbf - temperature before the reaction, but after thermalisation
    g_Temp = np.interp(C_Z, sim1D_z, sim1D_Tg)
    
    # get particle surface 
    p_surf = PI * p_diam**2 #Unit of p_surf is m^2

    # save current case temperature to the global array
    Tp_glob = np.concatenate((Tp_glob,p_Temp))
    # save current O/Fe ratio to the global array (convert from O/Fe ratio to 1 normalised)
    ratnew = rat / 3 * 2
    cov_OFe = np.concatenate((cov_OFe, rat / 3 * 2))
    
    # Handle the iron species concentration in the same way as in the Monte-Carlo simulations
    N_permass = sim1D_XFEC5O5 * pressure / k_B / sim1D_Tg / sim1D_rho       # get IPC concentration in conserved quantities (1/kg)
    Delta_N = N_permass[0] - N_permass                                      # IPC loss w.r.t. initial value
    X_FE = Delta_N * sim1D_rho * sim1D_Tg * k_B / pressure                  # Come back to mole fractions
    X_FE += -sim1D_XFE2O3 * 2 - sim1D_XFEO2 - sim1D_XFEO - sim1D_XFEOH - sim1D_XFEO2H2 - sim1D_XFE2OOOH * 2     # Deduct iron oxidic species 
    X_FE = np.clip(X_FE,0.0,None)                                           # clip negative values
    c_FE = X_FE * 3000 / sim1D_Tg / R_g                                     # convert again to concentration units (kg/m^3)

    print("sim1D_z:", sim1D_z.shape)
    print("sim1D_Tg:", sim1D_Tg.shape)

    print("c_O2:", c_O2.shape)
    print("c_O :", c_O.shape)

    #interpolate concentrations to C_Z grid and convert to mol
    c_FE_int    = np.interp(C_Z, sim1D_z, c_FE) * 1e3 #Unit of c_FE is kmol/m^3, so multiply by 1e3 to get mol/m^3               
    c_O2_int    = np.interp(C_Z, sim1D_z, c_O2) * 1e3 
    c_O_int     = np.interp(C_Z, sim1D_z, c_O)  * 1e3
    c_H2O_int   = np.interp(C_Z, sim1D_z, c_H2O)  * 1e3
    c_H2_int    = np.interp(C_Z, sim1D_z, c_H2)  * 1e3
    c_H_int     = np.interp(C_Z, sim1D_z, c_H)  * 1e3
    c_OH_int    = np.interp(C_Z, sim1D_z, c_OH)  * 1e3
    
    # Mean free paths, Knudsen numbers, transition regime correctors
    gas_Temp    = np.interp(C_Z, sim1D_z, sim1D_Tg)
    constLam = k_B * gas_Temp / 2**0.5 / PI / pressure / p_diam * 2
    kn_FE  = constLam / d_FE**2
    kn_O2  = constLam / d_O2**2
    kn_O   = constLam / d_O**2
    kn_H2O = constLam / d_H2O**2
    kn_H2  = constLam / d_H2**2
    kn_H   = constLam / d_H**2
    kn_OH  = constLam / d_OH**2
    fKn_FE  = (1.333*kn_FE + 1.333*kn_FE**2) / (1 + 1.71*kn_FE + 1.333*kn_FE**2)
    fKn_O2  = (1.333*kn_O2 + 1.333*kn_O2**2) / (1 + 1.71*kn_O2 + 1.333*kn_O2**2)
    fKn_O  = (1.333*kn_O + 1.333*kn_O**2) / (1 + 1.71*kn_O + 1.333*kn_O**2)
    fKn_H2O  = (1.333*kn_H2O + 1.333*kn_H2O**2) / (1 + 1.71*kn_H2O + 1.333*kn_H2O**2)
    fKn_H2  = (1.333*kn_H2 + 1.333*kn_H2**2) / (1 + 1.71*kn_H2 + 1.333*kn_H2**2)
    fKn_H  = (1.333*kn_H + 1.333*kn_H**2) / (1 + 1.71*kn_H + 1.333*kn_H**2)
    fKn_OH  = (1.333*kn_OH + 1.333*kn_OH**2) / (1 + 1.71*kn_OH + 1.333*kn_OH**2)

    # Calculate the reaction constants for each species

    const = p_surf * N_A #* 1e3      # helper (1e3 because Cantera concentrations are in kmols/m^3)
    
    # special case: Fe condensation
    k_Fe_case = C_Fe / c_FE_int / N_A / (PI * k_B * p_Temp / 2 / M_FE * N_A) ** 0.5 / (p_diam)**2 / fKn_FE
    k_Fe = np.concatenate((k_Fe,k_Fe_case))     # add to the global array (storing all cases)

    k_O2_case = C_O2 / c_O2_int / const / (1-ratnew) / fKn_O2
    k_O2 = np.concatenate((k_O2,k_O2_case))
    k_O2_log = np.log(k_O2_case * N_A)

    k_O_case = C_O / c_O_int / const / (1-ratnew) / fKn_O
    k_O = np.concatenate((k_O,k_O_case))
    k_O_log = np.log(k_O_case * N_A)
    
    k_H2O_case = C_H2O / c_H2O_int / const / (1-ratnew) / fKn_H2O
    k_H2O = np.concatenate((k_H2O,k_H2O_case))
    k_H2O_log = np.log(k_H2O_case * N_A)

    k_H2_case = C_H2 / c_H2_int / const / ratnew / fKn_H2
    k_H2 = np.concatenate((k_H2,k_H2_case))
    k_H2_log = np.log(k_H2_case * N_A)
    
    k_H_case = C_H / c_H_int / const / ratnew / fKn_H
    k_H = np.concatenate((k_H,k_H_case))
    k_H_log = np.log(k_H_case * N_A)
    
    k_OH_case = C_OH / c_OH_int / const/ (1-ratnew) / fKn_OH
    k_OH = np.concatenate((k_OH,k_OH_case))
    k_OH_log = np.log(k_OH_case * N_A)
    
        
    #### ARRHENIUS PLOTTING ####
 
    aax2.plot(1000/p_Temp, k_O2_log , '.', markersize = 0.4, color = phicol[n])#,alpha=0.5)
    aax3.plot(1000/p_Temp, k_O_log , '.', markersize = 0.4, color = phicol[n])#,alpha=0.5)
    aax4.plot(1000/p_Temp, k_H2O_log , '.', markersize = 0.4, color = phicol[n])#,alpha=0.5)
    aax5.plot(1000/p_Temp, k_H2_log , '.', markersize = 0.4, color = phicol[n])#,alpha=0.5)
    aax6.plot(1000/p_Temp, k_H_log , '.', markersize = 0.4, color = phicol[n])#,alpha=0.5)
    aax7.plot(1000/p_Temp, k_OH_log , '.', markersize = 0.4, color = phicol[n])#,alpha=0.5)
    
    aax2.xaxis.set_tick_params(labelbottom=False)
    aax3.xaxis.set_tick_params(labelbottom=False)
    aax4.xaxis.set_tick_params(labelbottom=False)
    aax5.xaxis.set_tick_params(labelbottom=False)
  
    #aax2.yaxis.set_tick_params(labelleft=False)
    aax3.yaxis.set_tick_params(labelleft=False)
    aax5.yaxis.set_tick_params(labelleft=False)
    aax7.yaxis.set_tick_params(labelleft=False)

    #### CURVE FITTING ####
    # polynomial fitting (rather USELESS, the first idea was to create a best fit line for species and then find arrh. parameters) 
    order = 3   #polynomial order

    # poly fitting of O/Fe w.r.t. particle temperature (controversial, later used for arrh. curve check)
    param_rat = np.polyfit(p_Temp, rat / 3 * 2, order)

    # create array of points for the best fit line (polynomial)
    for o in range(order+1):
        rat_bf[n][:]   += param_rat[o] *T_bf**(order-o)

# Store the global arrays as numpy arrays
Tp_glob = np.array(Tp_glob)
cov_OFe = np.array(cov_OFe)

min_len = min(len(Tp_glob), len(cov_OFe))

Tp_glob = Tp_glob[:min_len]
cov_OFe = cov_OFe[:min_len]

#### GLOBAL PLOTTING SETTINGS ####

aax2.set_ylim(57.8,62.8)
aax3.set_ylim(57.8,62.8)
aax4.set_ylim(57.8,62.8)
aax5.set_ylim(57.8,62.8)
aax6.set_ylim(57.8,62.8)
aax7.set_ylim(57.8,62.8)

aax6.set_xticks([0.5,1.0])
aax7.set_xticks([0.5,1.0])

c1 = mlines.Line2D([], [], linestyle = '-', color='m', linewidth=4,label='$\phi = 0.5$')# L500')
c2 = mlines.Line2D([], [], linestyle = '-', color='c',linewidth=4, label='$\phi = 1.0$')# S500')
c3 = mlines.Line2D([], [], linestyle = '-', color='y',linewidth=4, label='$\phi = 1.5$')# R500')
lns = [c1, c2, c3]
labs = [l.get_label() for l in lns]

legend = plt.legend(lns, labs, loc="upper center", bbox_to_anchor=(-1.1,3.25,2,0.8),
                mode="expand", borderaxespad=0, ncol=4,columnspacing=0.12,frameon=False, labelspacing=0.1, handletextpad=0.5)

fig2.supylabel(r'ln($k$, ms$^{-1}$mol$^{-1})$', fontsize = 13)
fig2.supxlabel(r'1000/$T$, K$^{-1}$', fontsize = 13)

aax2.set_title(r'R1: O$_{2}$')
aax3.set_title('R2: O')
aax4.set_title(r'R3: H$_{2}$O')
aax5.set_title(r'R4: H$_{2}$')
aax6.set_title('R5: H')
aax7.set_title('R6: OH')

#### ARRHENIUS PARAMETER FITTING ####

# Arrhenius equation with surface coverage parameters
def arrhSurf(Tp_glob, cov_OFe, A, b, E_a, alpha_k, m_k, E_k):
    rateConst = A * (Tp_glob/T_ref)**b * np.exp((-E_a - E_k * cov_OFe) / R_g / Tp_glob) * 10**(cov_OFe*alpha_k) * cov_OFe ** m_k 
    return rateConst

def arrhPlot(Tp_inv, A, b, E_a):
    arrhPlotPoint = np.log(A) - b * (np.log(Tp_inv) - np.log(1/T_ref)) - E_a  / R_g * Tp_inv 
    return arrhPlotPoint

# loss function to be minimised (square difference between arrh. curve and the actual datapoint)
# def loss(para,k_spec, method):
#     if method == 'lin':
#         diff = k_spec - arrhSurf(Tp_glob, cov_OFe, para[0], para[1], para[2], para[3], para[4], para[5])
#     elif method == 'log': 
#         k_fit = arrhSurf(Tp_glob, cov_OFe, para[0], para[1], para[2], para[3], para[4], para[5])
#         eps = 1e-300
#         k_spec = np.maximum(k_spec, eps) ###### von mir eingefügt, damit log nicht negativ wird
#         diff = np.log(np.maximum(k_spec, eps)) - np.log(np.maximum(k_fit, eps))
#     return np.sum(diff**2)

def loss(para, k_spec, method):
    k_fit = arrhSurf(Tp_glob, cov_OFe, para[0], para[1], para[2], para[3], para[4], para[5])
    eps = 1e-300
    # Schutz gegen 0 / negative Werte
    k_spec_safe = np.clip(k_spec, eps, None)
    k_fit_safe  = np.clip(k_fit,  eps, None)
    if method == 'lin':
        diff = k_spec_safe - k_fit_safe
    elif method == 'log':
        diff = np.log(k_spec_safe) - np.log(k_fit_safe)
    else:
        raise ValueError("method must be 'lin' or 'log'")
    return np.sum(diff**2)

# log-RMSE: 
def log_rmse(k_true, k_fit):
    eps = 1e-300
    res = np.log(np.maximum(k_true, eps)) - np.log(np.maximum(k_fit, eps))
    return np.sqrt(np.mean(res**2))

# Arrange the species data in iterable arrays
specArr  = ['O2', 'O', 'H2O', 'H2', 'H', 'OH']
constArr = [k_O2, k_O, k_H2O, k_H2, k_H, k_OH]
axArr    = [aax2,  aax3,  aax4, aax5, aax6, aax7]
method = ['lin', 'lin', 'log', 'log', 'log', 'lin']
method_log = ['log', 'log', 'log', 'log', 'log', 'log']
method_lin = ['lin', 'lin', 'lin', 'lin', 'lin', 'lin']

# Activation temperature for each species as used in SPMC
actTemp =  np.array([1000.0,0.0,1300.0,2000.0,500.0,2200.0])
actEner = actTemp * R_g

# Set minimisation solver method
solver = 'Mixed'
#solvArr = ['Nelder-Mead',  'Nelder-Mead', 'Nelder-Mead', 'Nelder-Mead',  'Nelder-Mead', 'Nelder-Mead']
solvArr = ['L-BFGS-B', 'L-BFGS-B',  'L-BFGS-B', 'L-BFGS-B', 'L-BFGS-B',  'L-BFGS-B', 'L-BFGS-B']
#solvArr = ['SLSQP', 'SLSQP', 'SLSQP', 'SLSQP', 'SLSQP',  'SLSQP', 'SLSQP']

# initialise the parameters output array
opt_Spec = np.zeros((len(specArr),6))
# store minimised loss functions (for comparison between the methods)
storeFun = np.zeros(len(specArr))

# Perform optimisation with scipy.minimize for each species
init_guess =  np.ones((6,))
#init_guess =  np.array([1e10,1e10,1e10,1e10,1e10,1e10])

if not b_switch:
    b_bound = (0.0,0.0)
else:
    b_bound = (None,None)
if not mk_switch:
    mk_bound = (0.0,0.0)
else:
    mk_bound = (None,None)
if not alpha_switch:
    alpha_bound = (0.0,0.0)
else:
    alpha_bound = (None,None)
if not Ek_switch:
    Ek_bound = (0.0,0.0)
else:
    Ek_bound = (None,None)

os.makedirs("./ArrheniusRates", exist_ok=True)
csv_file = open(f"./ArrheniusRates/rates_calc_test15.csv", 'w', newline='')
csv_writer = csv.writer(csv_file)

# original code
for s, spec in enumerate(constArr):
    if blockE_a:
        Ea_bound = (actEner[s], actEner[s])
    else:
        Ea_bound = (None,None)
    result = minimize(loss, init_guess, bounds = ((None,None),b_bound,Ea_bound,alpha_bound,mk_bound,Ek_bound), method=solvArr[s], args = (spec,method_log[s]))
    opt_Spec[s,:] = result.x
    opt_Spec[s,0] *= N_A
    storeFun[s] = result.fun
    for i in range(n_cases):
        checkFun = arrhPlot(1/T_bf, opt_Spec[s,0], opt_Spec[s,1], opt_Spec[s,2])
        axArr[s].plot(np.flip(1000/T_bf), np.flip(checkFun ), '--', linewidth = 0.7, color = 'k')
    print(f'Species {specArr[s]} : {opt_Spec[s]}')

    # Fit evaluation
    k_fit = arrhSurf(Tp_glob, cov_OFe, *result.x)
    q_logrmse = log_rmse(spec, k_fit)
    print(f'{specArr[s]} log-RMSE = {q_logrmse:.4f}')
    csv_writer.writerow(np.append(opt_Spec[s], q_logrmse))

os.makedirs("./ArrheniusPlots", exist_ok=True)
fig2.savefig(f"./ArrheniusPlots/ArrhPlots_test15.png", dpi=600)