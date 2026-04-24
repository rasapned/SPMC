#!/usr/bin/python3
import sys
import math as mt
import numpy as np
from scipy.special import erf
from scipy.interpolate import interp1d as interp
from random import random
import time
import csv
import yaml
import cantera as ct

# Define file names
case_file_name = sys.argv[1]
output_file = sys.argv[2]
case_base_name = case_file_name.split('.')[0]
SPMC_save_file = case_base_name + "_results/" + case_base_name+'-SPMC-'+ output_file+ '.csv'
CFD_T_of_x_file = 'T_of_x-'+case_base_name+'.csv'
cantera_save_file = case_base_name + "_results/" + case_base_name+'-cantera.csv'
nuclest_file = case_base_name + "_results/" + case_base_name+'-nuclest.csv'
with open(case_file_name) as f:
    for line in f:
        line = line.split("#", 1)[0].strip()
        if line.startswith("yaml_file"):
            rxmech = line.split()[1]
            break

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

################### READ FLAME DATA FROM THE CSV FILE (GENERATED THROUGH CANTERA)
data = np.genfromtxt(cantera_save_file, delimiter=",", names=True)

# Read basic flame data
f_z   = data["grid"]
f_T   = data["T"]
T_burner = f_T[0]
f_rho = data["density"]
f_vel = data["velocity"]
f_P   = data["pressure"]

# Flame mesh no. of cells
f_points = len(f_z)
# Initialise cantera time array (derive from grid and velocity)
f_time = np.zeros(f_points)
for i in range(1, f_points):
    f_time[i] = f_time[i-1] + (f_z[i] - f_z[i-1]) / (0.5 * (f_vel[i] + f_vel[i-1]))

#Some thermo species data have to be taken from Cantera
gas = ct.Solution(rxmech)
    
################### READ FLAME DATA END #######################################

################### READ PARTICLE PROPERTIES ##################################

# Open input file
with open('part_params.yaml', 'r') as f:
    conf = yaml.safe_load(f)

# Particle parameters
ni_FE_pMC    = int(conf['particle']['initial']['n_Fe'])         # initial number of metal atoms
ni_O_pMC     = int(conf['particle']['initial']['n_O'])          # initial number of oxygen
if conf['particle']['initial']['T_p'] == 'auto':               # Determine initial particle temperature:
    T_pMC    = f_T[0]                                                  # Same at gas (automatic)
else:
    T_pMC   = float(conf['particle']['initial']['T_p'])               # Read the value
metal_spec  = str(conf['particle']['initial']['metal_spec'])   # name of the particle metallic species
rho_p_fe    = float(conf['particle']['material']['rho_fe'])    # Metal density
rho_p_fe2o3 = float(conf['particle']['material']['rho_fe2o3']) # Maximum metal oxidation species density
pST_max     = float(conf['particle']['material']['pST_max'])   # maximum oxygen/metal ratio
C_pMC       = float(conf['particle']['material']['C_p'])       # particle specific heat capacity 
epsilon_p   = float(conf['particle']['material']['epsilon_p']) # particle emissivity 

# Gas metal reservoir properties (nucleation properties)
nucl_mode   = str(conf['nucleation']['mode'])
# Remaining variables read further in the code

# Evaporation parameters
H_evap     = float(conf['evaporation']['H_evap']) / N_A        # thermal dissociation of Fe and O expressed as an Arrhenius fit to literature values
A_evap     = float(conf['evaporation']['A_evap'])              # Arrhenius fit preexponential factor
T_Aevap    = float(conf['evaporation']['T_Aevap'])             # Arrhenius fit activation temperature
Evap_relax = float(conf['evaporation']['relax'])               # evaporation relaxation factor

# Collision/Material parameters
elast_accomodation = float(conf['collision']['elastic_accommodation'])

# Collision parameters and definitions
class proto_SOI:        # Define the collision/reaction class
    def __init__(self, name, shape, T_a, Hr, type, star, addO, addFE, sc1, sc2):
        self.name = name
        self.shape = shape
        self.T_a = T_a
        self.Hr = Hr / N_A
        self.type = type
        self.star = star
        self.addO = addO
        self.addFE = addFE
        self.sc1 = sc1
        self.sc2 = sc2
        self.coll_count = 0
        self.coll_comp = 0.0
SOI = []

for c in conf['collisions']:        # Read collision/reaction types
    # Logic for stick/oxi/red/etch coefficients
    # interaction if (random() < (sc1[s] + sc2[s]*pST))
    #                 sc1=+1, sc2=-1 : interaction with surface-Fe --> oxydation or etching
    #                 sc1= 0, sc2=+1 : interaction with surface-O  --> reduction
    #                 sc1=+1, sc2=+1 : stick always
    #                 sc1= 0, sc2= 0 : stick never
   if c["type"] == "stick":
        sc1, sc2 = 1, 1
    elif c["type"] == "oxi":
        sc1, sc2 = 1, -1
    elif c["type"] == "red":
        sc1, sc2 = 0, 1
    elif c["type"] == "etch":
        sc1, sc2 = 1, -1
    elif c["type"] == "elastic":
        sc1, sc2 = 0, 0
    else:
        raise ValueError(f"Unknown collision type: {c['type']}")
    
    SOI.append(proto_SOI(
        c['name'], c['shape'], c['T_a'], c['Hr'], c['type'], c['star']
        c['addO'], c['addFE'], sc1, sc2
    ))

# Read MC parameters
if conf['MC_params']['z_max'] == 'auto':                    # MC sim. until z_max
    z_max    = f_z[-1]                                          # Same as last Cantera flame point
else:
    T_pMC   = float(conf['particle']['initial']['T_p'])         # or read the value from input file
MBD_points  = int(conf['MC_params']['MBD_points'])      # Maxwell-Boltzmann func. distribution resolution
MBD_width   = float(conf['MC_params']['MBD_width'])     # Distribution width (multiple of mean molecule speed)
MBD_deltaT  = float(conf['MC_params']['MBD_deltaT'])    # max. temp. deviation from current before recalculating Maxwell-Botzmann

################### READ PARTICLE PROPERTIES END ##############################

################### READ NUCLEATION FILE ######################################

with open(nuclest_file) as f:
    first_line = f.readline().strip()

# Start simulation after a certain amount of precursor decomposes
# (prevents particle extinction in unfavourable growth conditions)
start_idx = int(first_line.split(":")[1].strip())

# Read estimated number of monomers (nuclei)
f_z_nucl, ni_mono = np.genfromtxt(str(nuclest_file), delimiter=',', comments='#').T # np.ndarray.T stands for 'transpose'

################### READ NUCLEATION FILE END ##################################

################### NOW DO THE MONTE CARLO MAGIC ##############################

#---------------------------------
# Initialise all necessary vectors
#---------------------------------

# Open the results file
MC_file = open(SPMC_save_file,'w')

# For each collision species, read the gas mole fractions from Cantera file
n_species = len(SOI)
f_X = np.zeros((n_species, f_points))
for s in range(n_species):
    col_name = "X_" + SOI[s].name.upper()
    f_X[s] = data[col_name]
    # Get the particle metal species index
    if col_name == ('X_' + metal_spec):
        metal_idx = s
    if col_name == 'X_O':
        oxi_idx = s

# Define the gas metal reservoir (nucleation source)
# Two modes available:
#   - add: add all defined metallic species
#   - decomp: precursor decomposed minus selected removed species 
if nucl_mode == 'decomp':
    # Read precursor from the input file
    prec_name = conf['nucleation']['precursor']
    X_prec = data[prec_name]
    # Read removed species from input file
    removed_species = conf['nucleation'].get('removed_species',[])
    factor = conf['nucleation'].get('factor',[])
    
    # Calculate number of molecules
    N_permass = X_prec * f_P / k_B / f_T / f_rho
    Delta_N = N_permass[0] - N_permass

    # Recalculate the particle metal gas pool
    f_X[metal_idx] = Delta_N * f_rho * f_T * k_B / f_P
    for spec, fac in zip(removed_species, factor):
        print(spec+' removed from iron concentration')
        col_name = "X_" + spec
        if col_name not in data.dtype.names:         # Check if the species exists in the CSV cantera file
            raise ValueError(f"{col_name} not found in input file")
        remove = data[col_name] * fac
        f_X[metal_idx] -= remove
        f_X[metal_idx] = np.clip(f_X[metal_idx],0.0,None)

elif nucl_mode == 'add':
    # Read added species from input file
    added_species = conf['nucleation'].get('added_species',[])
    factor = conf['nucleation'].get('factor',[])
    
    # Create temp metal array
    f_X_metal = np.zeros_like(f_z)

    for spec, fac in zip(added_species, factor):
        col_name = "X_" + spec  # column name in CSV

        if col_name not in data.dtype.names:         # Check if the species exists in the CSV cantera file
            raise ValueError(f"{col_name} not found in input file")

        print(f"{spec} added to the iron concentration")

        f_X_metal += data[col_name] * fac

    f_X[metal_idx] = f_XFE

# Initialise vectors for thermo species data (taken usually from Cantera)
W_g     = np.zeros(n_species)       # molar weight
M_g     = np.zeros(n_species)       # molecule weight    
s_g     = np.zeros(n_species)       # sigma (transport diameter)
geo_g   = [None]*n_species          # molecule geometry
rel_g   = np.zeros(n_species)       # rotational relaxation

# fetch SOI data from Cantera mechanism
for s in range(len(SOI)):
    s_SOI      = gas.species_index(SOI[s].name)
    W_g[s]     = gas.molecular_weights[s_SOI]
    M_g[s]     = W_g[s] / N_A * 1e-3
    s_g[s]     = gas.species(s_SOI).transport.diameter * 0.5        #WHY *0.5???? TO BE CHECKED
    geo_g[s]   = gas.species(s_SOI).transport.geometry
    rel_g[s]   = gas.species(s_SOI).transport.rotational_relaxation

# Initialise some incident gas molecules instantaneous arrays
X_g   = np.zeros(n_species)       # incident gas spec. mole fraction
n_g   = np.zeros(n_species)       # incident gas spec. number of molecules
cp_g  = np.zeros(n_species)       # incident gas spec. heat capacity

# Initialise particle-gas (pg) collision arrays
f_pg    = np.zeros(n_species)          # collision frequency
M_red   = np.zeros(n_species)          # particle-gas molecule reduced mass 

# Initialise MonteCarlo particle properties variables
M_pMC   = ni_FE_pMC * M_g[metal_idx] + ni_O_pMC * M_g[oxi_idx]  # total mass of the particle based on number of atoms
E_pMC   = C_pMC * T_pMC * M_pMC                                 # particle internal energy
pST_pMC = ni_O_pMC/ni_FE_pMC                                    # particle stoichiometry O/FE
rho_pMC = rho_p_fe - (rho_p_fe - rho_p_fe2o3)*pST_pMC/pST_max   # stoichiometry weighted particle density
R_pMC   = (3.0*(M_pMC/rho_pMC)/4.0/PI)**(1.0/3.0)               # particle radius

# Collision counters
n_ell = 0               # total number of elastic collisions
n_ine = 0               # total number of inelastic collisions
count_ell = 0           # currently not in use, but later for reduced output

# Record computation real time arrays (after progressing by 1 mm along the centreline)
t_hist = []
z_ref = f_z[start_idx]
t_ref = time.time()
CompTime_file = open(case_base_name + '_results/CompTime.dat','w')

#---------------------------------
# Maxwell Boltzmann cum. function
#---------------------------------

def CMBD(v, m, T):
    ''' Cumulative Distribution function of the Maxwell-Boltzmann speed distribution '''
    k_B = 1.38065e-23
    a = np.sqrt(k_B*T/m)
    # expression could be simplified but is used for conformity with text books
    return erf(v/(np.sqrt(2)*a)) - np.sqrt(2/np.pi)* v* np.exp(-v**2/(2*a**2))/a
    
MBD_v   = np.zeros((n_species,MBD_points))      # CMBD x-axis (molecule velocity)
MBD_cdf = np.zeros((n_species,MBD_points))      # CMBD y-axis (cum. probability)
MBD_T   = 100                                   # CMBD temperature (if it deviates too much from current, recalc. CMBD)

#---------------------------------
# Start MC loop
#---------------------------------

z = f_z[start_idx]              # Start point
t = f_time[start_idx]           # Start time
icount = 1                      # Iteration counter
# Results file 
MC_head_str = '#time, z, Tg, Tp, Dp, Mp,S p, n_ine, n_ell, n_coll, n_FE, n_O, SOI, coll_type'
MC_file.write(MC_head_str+'\n')

# Loop start
while (z<z_max):
   
    #''''''''''''''''''''''''''''''''''''
    # Interpolate current position & flame conditions (Cantera resolution much larger than MC)
    #''''''''''''''''''''''''''''''''''''
    T_g = np.interp(t, f_time, f_T)             # Gas temp
    P_g = np.interp(t, f_time, f_P)             # Pressure
    z = np.interp(t, f_time, f_z)               # Grid position
    n_FE_seeds = np.interp(t, f_time, ni_mono)  # Seeds
    
    # Set current gas state (to get correct thermo data)
    gas.TP = T_g, P_g
    
    # Update gas metal concentration (Fe depletion/coupling)
    Fe_depl = n_FE_seeds * MC_ni_FE
    
    for s in range(len(SOI)):
        s_SOI    = gas.species_index(SOI[s].name)               # species name
        X_g[s]   = np.interp(t, f_time, f_X[s])                 # species mol. frac.
        n_g[s]   = P_g / T_g / k_B * X_g[s] - (s==0)*Fe_depl    # convert to no. of molecules
        cp_g[s]  = gas.species(s_SOI).thermo.cp(T_g) / W_g[s]   # specific heat 

    #''''''''''''''''''''''''''''''''''''
    # (re-) Compute the Maxwell-Boltznamm distribution if T_g deviates too much from MBD_T
    #''''''''''''''''''''''''''''''''''''
    if abs(T_g-MBD_T) > MBD_deltaT:
        for s in range(len(SOI)):                       # Build MB distribution for each incident species

            # Build the velocity axis (x-axis)
            v_mean = np.sqrt(3.0 * k_B * T_g / M_g[s])  # kinetic th. mean mol. speed
            v_max = v_mean * MBD_width                  # max. distribution width (speed)     
            v_inc = v_max / (len(MBD_v[s]) + 1)         # distribution (speed) increment  
            for j in range(1, len(MBD_v[s])):           # Build the distribution x axis for each MBD point
                MBD_v[s][j] = MBD_v[s][j-1] + v_inc
            # Build the probability axis (y-axis)
            MBD_cdf[s] = CMBD(vs_g[s], M_g[s], T_g) 
        
        # Update the temperature for which the CMBD was calculated
        MBD_T = T_g                                     
        
    #''''''''''''''''''''''''''''''''''''
    # kMC algorithm: collision frequencies and timestepping
    #''''''''''''''''''''''''''''''''''''
    for s in range(len(SOI)):
        M_red[s] = (MC_M * M_g[s]) / (MC_M + M_g[s])    # Particle-gas molec. reduced mass
        f_pg[s]   = n_g[s] * PI * (MC_R + s_g[s])**2 * (8.0*k_B*T_g/PI/M_red[s])**0.5   # species collision frequency
    f_cum = sum(f_pg)                                   # Cumulative collision frequency (of any species molecule)        
    dt_av = 1.0/f_cum                                   # Collision characteristic time (average)
    dt = -np.log(max(random(), 1.0e-7)) * dt_av         # MC time step (derived from Poisson distribution)

    # computed averaged collision event: how many theoretical collision events of each species would happen within current timestep?
    for s in range(len(SOI)):
        SOI[s].coll_comp = SOI[s].coll_comp + dt*f_pg[s]

    #''''''''''''''''''''''''''''''''''''
    # kMC algorithm: determine the collision partner (random process)
    #''''''''''''''''''''''''''''''''''''
    rnd = random()
    prob = 0.
    for s in range(len(SOI)):
        prob += f_pg[s]*dt_av
        if rnd < prob:          # if rnd smaller than cumulative prob => select species
            break
    
    # actual collision event with species s: increase the count
    SOI[s].coll_count = SOI[s].coll_count + 1
    
    #''''''''''''''''''''''''''''''''''''
    # Compute the kinetic energy of incident gas molecule sampled from Maxwell-Boltzmann distribution (random process)
    #''''''''''''''''''''''''''''''''''''
    v_gm = interp(cdf[s], vs_g[s])(random())            # gas molecule velocity sampled from MBD
    E_gm = 0.5 * M_g[s] * v_gm**2                       # gas molecule kinetic energy
    T_gm = M_g[s] * v_gm**2 / 3.0 / k_B                 # convert to gas molecule temperature
    
    #''''''''''''''''''''''''''''''''''''
    # Single particle just before collision
    #''''''''''''''''''''''''''''''''''''
    # Particle stoichiometry ratio: 0: metal only, 1: fully oxidised
    pST_pMC = ni_O_pMC/ni_FE_pMC
    pST_rat = pST_pMC / pST_max
    # Save particle temperature before the collision and potential reaction (this temp. will be considered for Arrh. constants)
    Tbefore_pMC = T_pMC
    
    #''''''''''''''''''''''''''''''''''''
    # Single particle - gas molecule collision
    #''''''''''''''''''''''''''''''''''''
    E_pMC += -(4+SOI[s].shape)*0.5*k_B*(T_pMC - T_gm) * elast_accomodation  # Particle energy change (equilibration with gas molecule)
    T_pMC = E_pMC / M_pMC / C_pMC                                           # Updated particle temperature
 
    #''''''''''''''''''''''''''''''''''''
    # Determine collision partner at the particle surface (random process)
    #''''''''''''''''''''''''''''''''''''
    if random() < (SOI[s].sc1 + SOI[s].sc2*pST_rat):        # True: matching reaction partner; else, elastic collision
        
        #''''''''''''''''''''''''''''''''''''
        # Again MBD for the gas molecule (only applicable if stoichiometry condition for reactive collision is met) (random process) 
        #''''''''''''''''''''''''''''''''''''
        v_mean = np.sqrt(3.0 * k_B * T_pMC / M_g[s])
        v_max = v_mean * MBD_width
        v_arr = np.linspace(0.0,v_max,num=MBD_points)
        
        cdf_arr = CMBD(v_arr, M_g[s], T_pMC) 

        v_cur = interp(cdf_arr, v_arr)(random())
        T_cur = M_g[s] * v_cur**2 / 3.0 / k_B

        #''''''''''''''''''''''''''''''''''''
        # Potential inelastic collision (if molecule temperature > activation temperature)
        #''''''''''''''''''''''''''''''''''''
        if (T_cur > SOI[s].T_a):                    # else, elastic col.
            print(f'Part. temp:  {T_pMC:6.1f}')                   
            print(f'Molec. temp: {T_cur:6.1f}')
            ni_FE_pMC += SOI[s].addFE               # add metal atoms
            ni_O_pMC  += SOI[s].addO                # add oxygen atoms
            star = SOI[s].star
            n_ine += 1                              # inelastic coll. increase
            count_ell = 0
            n_coll = n_ine+n_ell
            M_pMC  = ni_FE_pMC * M_g[metal_idx] + ni_O_pMC * M_g[oxi_idx]
            pST_pMC = ni_O_pMC/ni_FE_pMC
            rho_pMC = rho_p_fe - (rho_p_fe - rho_p_fe2o3)*pST_pMC/pST_max
            R_pMC = (3.0*(M_pMC/rho_pMC)/4.0/PI)**(1.0/3.0)
            if s==0 and random() < pST_rat:
                E_pMC += 289.0 * 1e3 / N_A
            if (s==1 or s==2):
                E_pMC += SOI[s].Hr - (pST_rat * pST_max / (1 + pST_rat*pST_max) * 127.42 * 1e3 / N_A)
            else:
                E_pMC += SOI[s].Hr
            T_pMC = E_pMC / M_pMC / C_pMC
            #
            outstr = f"{t:1.8e}, {z:1.8e}, {T_g:6.1f}, {T_pMC:6.1f}, {R_pMC*2:1.5e}, {M_pMC:1.5e}, {pST_pMC:4.3f}, {n_ine:6d}, {n_ell:6d}, {n_coll:6d}, {ni_FE_pMC:5d}, {ni_O_pMC:5d}, {SOI[s].name:3s}, {star:3s}, {Tbefore_pMC:6.1f}"
            print(outstr)
            MC_file.write(outstr+'\n')
        
        # elastic
        else:
            n_ell += 1
            n_coll = n_ine+n_ell
            count_ell += 1
            star = 'el'
            # this part is actually not needed anymore, as after the second MDB roll we have plenty of reactive reactions at all stages
            if count_ell == 200:
                outstr = f"{t:1.8e}, {z:1.8e}, {T_g:6.1f}, {T_pMC:6.1f}, {R_pMC*2:1.5e}, {M_pMC:1.5e}, {pST_pMC:4.3f}, {n_ine:6d}, {n_ell:6d}, {n_coll:6d}, {ni_FE_pMC:5d}, {ni_O_pMC:5d}, {SOI[s].name:3s}, {star:3s}, {Tbefore_pMC:6.1f}"
                count_ell = 0
                print(outstr)
                MC_file.write(outstr+'\n')

    # same, elastic
    else:
        n_ell += 1
        n_coll = n_ine+n_ell
        count_ell += 1
        star = 'el'
        # this part is actually not needed anymore, as after the second MDB roll we have plenty of reactive reactions at all stages
        if count_ell == 200:
            outstr = f"{t:1.8e}, {z:1.8e}, {T_g:6.1f}, {MC_Tp:6.1f}, {MC_R*2:1.5e}, {MC_M:1.5e}, {MC_ST:4.3f}, {n_ine:6d}, {n_ell:6d}, {n_coll:6d}, {MC_ni_FE:5d}, {MC_ni_O:5d}, {SOI[s].name:3s}, {star:3s}, {MC_Tp_before:6.1f}"
            count_ell = 0
            print(outstr)
            MC_file.write(outstr+'\n')

    # Evaporation losses
    if (ni_O_pMC > 5 and ni_FE_pMC > 5):
        T_pMC = E_pMC / M_pMC / C_pMC
        ndot_FEO = A_evap * mt.exp(-T_Aevap/MC_Tp) * N_A * 4.0*PI*R_pMC**2.0 * Evap_relax
        n_FEO = mt.trunc(ndot_FEO * dt)
        # loosing Fe and O one by one until particle is too cold
        while (n_FEO>0 and ni_O_pMC > 8 and ni_FE_pMC > 8):
            print(E_pMC/M_pMC/C_pMC, ndot_FEO * dt, n_FEO)
            ni_O_pMC  -= 1   # loosing FE+O atoms
            ni_FE_pMC -= 1
            star = 'ev'
            M_pMC  = ni_FE_pMC * M_g[metal_idx] + ni_O_pMC * M_g[oxi_idx] # compute new mass of particle
            pST_pMC = ni_O_pMC / ni_FE_pMC
            rho_pMC = rho_p_fe - (rho_p_fe - rho_p_fe2o3)*pST_pMC/pST_max
            R_pMC = (3.0*(M_pMC/rho_pMC)/4.0/PI)**(1.0/3.0)
            E_pMC += H_evap        # loosing a lot of heat
            T_pMC = max(T_g, E_pMC/M_pMC/C_pMC)
            outstr = f"{t:1.8e}, {z:1.8e}, {T_g:6.1f}, {MC_Tp:6.1f}, {MC_R*2:1.5e}, {MC_M:1.5e}, {MC_ST:4.3f}, {n_ine:6d}, {n_ell:6d}, {n_coll:6d}, {MC_ni_FE:5d}, {MC_ni_O:5d}, {SOI[s].name:3s}, {star:3s}, {MC_Tp_before:6.1f}"
            print(outstr)
            MC_file.write(outstr+'\n')
            ndot_FEO = A_evap * mt.exp(-T_Aevap/T_pMC) * N_A * 4.0*PI* R_pMC**2.0 * Evap_relax
            n_FEO = mt.trunc(ndot_FEO * dt)
    
    # Radiation losses/gain according to Landstroem et al. (2005) https://doi.org/10.1007/s00339-005-3284-3
    # But at least, they are integrated at the end of the large dt
    dErdt_R =  -4.0*PI*(R_pMC)**2 * epsilon_p * SBp * (T_pMC-T_g)**5
    E_pMC +=  dErdt_R * dt

    icount +=1

    # record computation time after progressing by 1 mm forward
    if z >= x_ref+0.001:
        x_ref = z
        time_elapsed = time.time() - t_ref
        CompTime_file.write(f"{time_elapsed:1.8e}"+'\n')
        t_hist.append(time_elapsed)
        t_ref = time.time()

    t += dt
    if (ni_FE_pMC < 6) or (t > t_max):
        print('# particle burned or evaporated or time up - finalized')
        MC_file.close()
        CompTime_file.close()
        print('computed vs counted collisions')
        for s in range(len(SOI)):
            print(int(SOI[s].coll_comp), SOI[s].coll_count)
        exit()

print('computed vs counted collisions')
for s in range(len(SOI)):
    print(int(SOI[s].coll_comp), SOI[s].coll_count)
print('# maximum number of iterations reached')
MC_file.close()
CompTime_file.close()

