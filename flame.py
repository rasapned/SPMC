#!/usr/bin/python3
import sys
from pathlib import Path
import os
import math as mt
import numpy as np
from scipy.special import erf
from scipy.interpolate import interp1d as interp
from scipy.ndimage.interpolation import shift
from random import random
import multiprocessing
import time
import csv

def read_case_file(case_file_name):
    config = {
        "VSdot_H2": 0.0,
        "VSdot_O2": 0.0,
        "VSdot_AR": 0.0,
        "X_IPC": 0.0,
        "yaml_file": None,
        "scale_species": [],
        "scaling_factor": [],
    }

    with open(case_file_name, "r") as f:
        for line in f:
            line = line.split("#", 1)[0].strip()
            if not line:
                continue

            parts = line.split()
            key = parts[0]

            # Flag
            if len(parts) == 1:
                config[key] = True
                continue

            values = parts[1:]

            # Handle lists
            if key == "scale_species":
                config[key] = values

            elif key == "scaling_factor":
                config[key] = list(map(float, values))

            # Single value (your original logic)
            else:
                value = values[0]
                try:
                    if "." in value or "e" in value.lower():
                        value = float(value)
                    else:
                        value = int(value)
                except ValueError:
                    pass
                config[key] = value

    # Optional safety check
    if config["scale_species"] and config["scaling_factor"]:
        if len(config["scale_species"]) != len(config["scaling_factor"]):
            raise ValueError("scale_species and scaling_factor must have same length")

    return config

case_file_name = sys.argv[1]
case_base_name = case_file_name.split('.')[0]
print(f'case name: {case_base_name}')

SPMC_save_file = case_base_name+'-SPMC'
cantera_save_file = "RESULTS/" + case_base_name + "_results/" +  case_base_name+'-cantera.csv'
CFD_T_of_x_file = 'T_of_x-'+case_base_name+'.csv'
cantera_conc_file = "RESULTS/" + case_base_name + "_results/" + case_base_name+'-ConcCant.csv'


# set some constants
PI  = 3.14159265359             # 
N_A = 6.02214129e+23            # Avogadro number           in 1/mol CAUTION: cantera uses kg and kmol!
R_m = 8.3145                    # Univ. gas const.          in J/mol/K
k_B = R_m / N_A                 # Boltzmann const. 1.38e-23 in m^2*kg/s^2/K or J/K
h_P = 6.6261e-34                # Plancks const.            in m^2*kg/s
c_l = 299792458                 # Speed of light            in m/s
SBc = 5.670374419e-8            # Stefan-Boltzmann const.   in J/s/m^2/K^4
SBp = k_B/h_P/c_l * SBc * 3.84  # Modified SBc              in J/s/m^2/K^5 
epsilon_p = 4.e-8               # Particle emissivity       in m
# Some small and big value 
small = 1e-31                   # a sufficiently small number
large = 1e+31                   # a sufficiently large number
maxDt = 1.0e-6                  # maximum time between collisions

###################### CANTERA PART BEGIN ######################################
import cantera as ct
ct.suppress_thermo_warnings()

# Read the input case file
cfg = read_case_file(case_file_name)


VSdot_H2 = cfg["VSdot_H2"]
VSdot_O2 = cfg["VSdot_O2"]
VSdot_AR = cfg["VSdot_AR"]
X_IPC    = cfg["X_IPC"]
yaml_file = cfg["yaml_file"]
# Sensibility study addition
scale_species  = cfg["scale_species"]
scaling_factor = cfg["scaling_factor"]

# chose a reaction mechanism
rxmech = yaml_file
print(rxmech)
# read the temperature profile, preferably from a simulation result
print(CFD_T_of_x_file)
data_file = Path(__file__).parent.joinpath(CFD_T_of_x_file)
T_values, z_loc = np.genfromtxt(str(data_file), delimiter=',', comments='#',usecols=(0,1),skip_header=1).T # np.ndarray.T stands for 'transpose'
z_loc_norm = z_loc / max(z_loc)

# define the case specific constants
D_burner = 36.0e-3       # m
pressure = 3000.0        # Pa
T_burner = T_values[0]   # K
T_stagna = T_values[-1]  # K
width_flame = z_loc[-1]  # m
print(f"T_burner = {T_burner:4.1f} K, T_stagna = {T_stagna:4.1f} K, domain width = {width_flame:1.2e} m")

VSdot_total = VSdot_H2 + VSdot_O2 + VSdot_AR

X_H2 = VSdot_H2 / VSdot_total
X_O2 = VSdot_O2 / VSdot_total
X_AR = 1.0 - (X_H2+X_O2+X_IPC)
# premixed, unburned gas composition
# rewritten for 'new' python string formating methods
comp = f"H2:{X_H2:.5e}, O2:{X_O2:.5e}, FEC5O5:{X_IPC:5e}, AR:{X_AR:.5e}"

print('composition = '+comp)

# create the cantera gas object
gas = ct.Solution(rxmech)

print([s for s in gas.species_names if "Fe" in s])

# compute the mean molecular weight of the unburned mixture
Wm = X_H2*gas.molecular_weights[gas.species_index('H2')] +\
X_O2*gas.molecular_weights[gas.species_index('O2')] +\
X_AR*gas.molecular_weights[gas.species_index('AR')] +\
X_IPC*gas.molecular_weights[gas.species_index('FEC5O5')]

# compute the mass flow rate density 
# Standard state as 1bar and 0°C,  ACHTUNG: Wm of Cantera is [kg/kmol]
rho_s = 1e5/273.15/R_m * (Wm*1e-3) 
print(f"standard density = {rho_s:1.3e} kg/m^3, Wm = {Wm:3.1f} kmol/kg")
mdot = rho_s * VSdot_total/60e6 / (0.25*PI*D_burner*D_burner)
print(f"mdot = {mdot:.5e} kg/s/m^2")

Phi = (X_H2/X_O2) / 2.0
print(f"Phi = {Phi:3.2f}")

S_H2 = VSdot_total * X_H2
S_O2 = VSdot_total * X_O2
S_AR = VSdot_total * X_AR

print('mass flow rates in sccm: S_H2={:3.1f} S_O2={:3.1f} S_AR={:3.1f}'.format(S_H2, S_O2, S_AR))

# set state of gas
gas.TPX = T_burner, pressure, comp

# set up flame object
f = ct.ImpingingJet(gas=gas, width=width_flame)
f.inlet.mdot = mdot
f.surface.T = T_stagna
#f.set_grid_min(5e-5)

#f.set_initial_guess(products='equil')
f.flame.set_fixed_temp_profile(z_loc_norm, T_values)
f.energy_enabled = False
loglevel = 0

f.set_refine_criteria(ratio=3, slope=0.1, curve=0.2, prune=0.06)
f.solve(loglevel, refine_grid=True)

# chemestry rate
print("reaction rates:", f.net_rates_of_progress)
print(f.net_rates_of_progress.shape)

X = f.X.copy()
T = f.T
P = f.P
grid = f.grid
velocity = f.velocity
spread_rate = f.spread_rate
lam = f.L  # lambda (strain rate)
density = f.density

# POWDER TECHNOLOGY (2026) ADDITION: OH SENSIBILITY STUDY
# START ARTIFICIALLY MODIFYING THE MOLE FRAC. (here, of OH)

# Copy the species mole fractions to an array
X_mod = f.X.copy()

# Iterate over species defined in the input file
for i,spec in enumerate(scale_species):
    # Get species index and scale the species
    indx = gas.species_index(spec)
    X_mod[indx, :] *= scaling_factor[i]

# renormalize each grid point
X_mod /= X_mod.sum(axis=0)

for j in range(f.flame.n_points):
    gas.TPX = T[j], P, X_mod[:, j]
    f.set_gas_state(j)

with open(cantera_save_file, 'w') as f_out:

    # header
    header = 'grid,velocity,spread rate,lambda,pressure,T,density,X_' + ',X_'.join(gas.species_names) + '\n'
    f_out.write(header)

    for j in range(len(grid)):
        row = [grid[j], velocity[j], spread_rate[j], lam[j], P,  T[j], density[j], *X_mod[:, j]]
        f_out.write(','.join(f"{val:.9e}" for val in row) + '\n')

# generating a file with the concentrations of the species of interest (for Arrhenius plot)
with open(cantera_conc_file, 'w') as f_out:

    # header
    header = 'c_O2,c_O,c_H2O,c_H2,c_H,c_OH,c_FEC5O5,c_FE2O3(s),c_FEO2,c_FEO,c_FEO2H2,c_FE2OOOH\n#Concentration given in kmol/m^3\n'
    f_out.write(header)
    
    W_mix = 0
    for spec in gas.species_names:
        #print(f"Species: {spec}, Index: {gas.species_index(spec)}")
        W_mix += X_mod[gas.species_index(spec), :] * gas.molecular_weights[gas.species_index(spec)]
    
    #for j in range(len(grid)):
    c_O2 = X_mod[gas.species_index('O2'), :] * density / W_mix #Concentration given in kmol/m^3
    c_O = X_mod[gas.species_index('O'), :] * density / W_mix
    c_H2O = X_mod[gas.species_index('H2O'), :] * density / W_mix
    c_H2 = X_mod[gas.species_index('H2'), :] * density / W_mix
    c_H = X_mod[gas.species_index('H'), :] * density / W_mix
    c_OH = X_mod[gas.species_index('OH'), :] * density / W_mix
    c_FEC5O5 = X_mod[gas.species_index('FEC5O5'), :] * density / W_mix
    c_FE2O3= X_mod[gas.species_index('FE2O3(s)'), :] * density / W_mix
    c_FEO2 = X_mod[gas.species_index('FEO2'), :] * density / W_mix
    c_FEO = X_mod[gas.species_index('FEO'), :] * density / W_mix
    c_FEO2H2 = X_mod[gas.species_index('FEO2H2'), :] * density / W_mix
    c_FE2OOOH = X_mod[gas.species_index('FE2OOOH'), :] * density / W_mix
    print("H2 conc =", c_OH)
    #print("All concentrations =", f.concentrations)
    print("f.H2 = ", f.concentrations[gas.species_index('OH'), :])
    #print("H2 conc formula =", c_H2[gas.species_index('H2')])
    #print("H2 conc cantera =", f.concentrations[gas.species_index('H2')])
    #print(X_mod.shape)
    #print(f.concentrations.shape)
    #print(gas.species_index('H2'))
    #print(f.X[gas.species_index('H2'),200])
    #print(X_mod[gas.species_index('H2'),200])
    #row = [c_O2, c_O, c_H2O, c_H2, c_H, c_OH, c_FEC5O5, c_FE2O3, c_FEO2, c_FEO, c_FEO2H2, c_FE2OOOH]
    #f_out.write(','.join(f"{val:.9e}" for val in row) + '\n') 
    data = np.column_stack((c_O2,c_O,c_H2O,c_H2,c_H,c_OH,c_FEC5O5,c_FE2O3,c_FEO2,c_FEO,c_FEO2H2,c_FE2OOOH))

    np.savetxt(f_out, data, delimiter=',',fmt='%.9e')

    # calculate chemistry rates


    # for j in range(len(grid)):
    #     c_O2 = X_mod[gas.species_index('O2'), j] * density[j] / gas.molecular_weights[gas.species_index('O2')] #Concentration given in kmol/m^3
    #     c_O = X_mod[gas.species_index('O'), j] * density[j] / gas.molecular_weights[gas.species_index('O')]
    #     c_H2O = X_mod[gas.species_index('H2O'), j] * density[j] / gas.molecular_weights[gas.species_index('H2O')]
    #     c_H2 = X_mod[gas.species_index('H2'), j] * density[j] / gas.molecular_weights[gas.species_index('H2')]
    #     c_H = X_mod[gas.species_index('H'), j] * density[j] / gas.molecular_weights[gas.species_index('H')]
    #     c_OH = X_mod[gas.species_index('OH'), j] * density[j] / gas.molecular_weights[gas.species_index('OH')]
    #     c_FEC5O5 = X_mod[gas.species_index('FEC5O5'), j] * density[j] / gas.molecular_weights[gas.species_index('FEC5O5')]
    #     c_FE2O3= X_mod[gas.species_index('FE2O3(s)'), j] * density[j] / gas.molecular_weights[gas.species_index('FE2O3(s)')]
    #     c_FEO2 = X_mod[gas.species_index('FEO2'), j] * density[j] / gas.molecular_weights[gas.species_index('FEO2')]
    #     c_FEO = X_mod[gas.species_index('FEO'), j] * density[j] / gas.molecular_weights[gas.species_index('FEO')]
    #     c_FEO2H2 = X_mod[gas.species_index('FEO2H2'), j] * density[j] / gas.molecular_weights[gas.species_index('FEO2H2')]
    #     c_FE2OOOH = X_mod[gas.species_index('FE2OOOH'), j] * density[j] / gas.molecular_weights[gas.species_index('FE2OOOH')]
    #     print(f.X[gas.species_index('H2'),10])
    #     print(X_mod[gas.species_index('H2'),10])
    #     print(f.X[gas.species_index('H'),10])
    #     print(X_mod[gas.species_index('H'),10])
    #     print(f.X[gas.species_index('OH'),10])
    #     print(X_mod[gas.species_index('OH'),10])
    #     print(f.X[gas.species_index('H2O'),10])
    #     print(X_mod[gas.species_index('H2O'),10])
    #     print("rho =", density[j])
    #     print("MW =", gas.mean_molecular_weight)
    #     print("H2 X =", X_mod[gas.species_index('H2'),j])
    #     print("H2 conc formula =", c_H2)
    #     print("H2 conc cantera =", gas.concentrations[gas.species_index('H2')])
    #     row = [c_O2, c_O, c_H2O, c_H2, c_H, c_OH, c_FEC5O5, c_FE2O3, c_FEO2, c_FEO, c_FEO2H2, c_FE2OOOH]
    #     f_out.write(','.join(f"{val:.9e}" for val in row) + '\n') 

###################### CANTERA PART END ########################################
