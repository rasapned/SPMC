import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from random import random
import matplotlib
import sys
import math as mt

######################################################################################################
# THIS CODE ESTIMATES THE NUMBER OF PARTICLE PRECURSOR SEEDS
# Necessary for the one way coupling. Rationale: IPC is added in small quantities, creation of particles
# depletes the Fe reservoir in the gas such that the particles can't grow infinitely. To account for this
# effect in a single particle simulation, we need to estimate the total partciles population.
######################################################################################################

# Set basic constants
PI  = 3.14159265359             # Pi
N_A = 6.02214129e+23            # Avogadro number           in 1/mol CAUTION: cantera uses kg and kmol!
R_m = 8.3145                    # Univ. gas const.          in J/mol/K
k_B = R_m / N_A                 # Boltzmann const. 1.38e-23 in m^2*kg/s^2/K or J/K

# Material constants - Fe
m_a_Fe = 55.845e-3              # molar mass
rho_Fe = 7874.0                 # density
r1 = 1.6e-10                    # Molecule radius

# Determine molecule properties like diameter, area, volume...
# Option 1: using predefined molecule radius r1
m1 = m_a_Fe / N_A
v1=4/3*PI*r1**3
s1=4*PI*r1**2
d1=2*r1
# Option 2: using bulk properties
#v1 = m1 / rho_Fe 
#r1 = (v1 * 3.0/4.0 / PI)**(1.0/3.0)
#d1 = r1 * 2
#s1 = 4*PI*r1**2

# Read the current flame case
case_base_name = sys.argv[1]
case = case_base_name.split('.')[0]
colour = 'royalblue'

# Read Cantera file: basic flame quantities
head_x = ['grid','T','density','velocity']
# Read Cantera file: Set which species should be considered as Fe particle monomers?
head_y = ['X_FE7','X_FE6','X_FE5','X_FE4','X_FE3','X_FE2','X_FE']#,'X_FEC4O4','X_FEC3O3','X_FEC2O2','X_FECO','X_FE2CO','X_FE3CO','X_FEO','X_FEOH','X_FEO2H2','X_FEO2','X_FEH','X_FEOOH','X_FE4CO','X_FE5CO','X_FE6CO','X_FE7CO','X_FE2C2O2','X_FE3C2O2','X_FE2C3O3','X_FE2OOOH','X_FE2OO2H2']
# How much Fe in a species?
factor = [7, 6, 5, 4, 3, 2, 1]#,1,1,1,1,2,3,1,1,1,1,1,1,4,5,6,7,2,3,2,2,2] 

# Create visualisation figures
# Create 2x2 grid
fig, axes = plt.subplots(2, 2, figsize=(16, 12), dpi=100)
# Assign each subplot
ax1, bx1, cx1, dx1 = axes.flatten()
# Tick params
for ax in [ax1, bx1, cx1, dx1]:
    ax.tick_params(axis='both', which='major', labelsize=24)
    ax.set_xlabel('DFB / mm', color='black', fontsize=18)
    ax.set_xlim([0, 8])
# Y labels
ax1.set_ylabel('Critical cluster size / -', fontsize=18)
ax1.set_ylim([0, 20])
bx1.set_ylabel(r'Number density of Fe monomers / #/m$^3$', fontsize=18)
cx1.set_ylabel(r'Instantaneous nascent nuclei / 1/m$^3$', fontsize=18)
dx1.set_ylabel(r'Cumulative nuclei / 1/m$^3$', fontsize=18)
# Logarithmic scales
dx1.set_yscale('log')

# Read the grid, velocity, species concentration and temperature data from cantera file
infile = case + "_results/" + case + '-cantera.csv'
# 1D sim data
header = np.genfromtxt(infile, delimiter=',',  dtype = None, encoding =None , comments = '#', max_rows = 1)
x_indx = np.in1d(header,head_x).nonzero()[0]
x,v,T,rho = np.genfromtxt(infile, delimiter=',', comments='#', skip_header = 1, usecols=x_indx).T

# Set the nucleation species array 
tot_y = np.zeros((len(T)))
for i, val in enumerate(head_y):
    y_indx = np.in1d(header,val).nonzero()[0]
    y = np.genfromtxt(infile, delimiter=',', comments='#',skip_header = 1, usecols=y_indx).T
    if val == 'X_FEC5O5':
        tot_y = y
        mode = 0
    else:
        tot_y += y*factor[i]
        mode = 1

# Calculate saturation equilibrium concentration of monomers
p_sat = 10**(11.353 - 19574.0 / T)      # vapour pressure
n_sat = p_sat / k_B / T                 # monomers concentration from ideal gas 

# Actual monomers concentration
if mode == 0:
    N_permass = tot_y * 3000 / k_B / T / rho
    Delta_N = N_permass[0] - N_permass
    n_mon = Delta_N * rho

elif mode == 1:
    n_mon = tot_y * 3000 / k_B / T

# Saturation ratio
Sat = n_mon/n_sat
# Surface energy (tension) dependence on the temperature 
sigma = 2.40328 - 2.85 * 1e-4 * T
# Dimensionless surface tension
theta = sigma * s1 / k_B / T
# Number of molecules in a critical cluster
g_star = (2.0*theta/3.0/np.log(Sat))**3
# Critical diameter
d_star = (g_star * v1 * 3.0/4.0 / PI) ** (1.0/3.0) * 2

# Saturation ratio - alternative: correction of the saturation equilibrium concentration of monomers:
#   - both surface tension (sigma) and saturation concentration (n_sat) can depend on the size of the critical cluster
#   - since they depend on each other, they have to be calculated iteratively
# Surface tension correction for small sizes of droplets (Tolman correction)
#    sigma_tolman = sigma / (1 + 4 * 0.1 * 1e-9 / d_star)
#    theta_tolman = sigma_tolman * s1 / k_B / T
#    g_star = (2.0*theta_tolman/3.0/np.log(Sat))**3
#    d_star = (g_star * v1 * 3.0/4.0 / PI) ** (1.0/3.0) * 2
#    r_star = d_star / 2.
#    Sat = np.exp(2 * sigma_tolman * m1 / (k_B * T * rho_Fe * r_star))
#    g_star = (2.0*theta_tolman/3.0/np.log(Sat))**3
#    n_sat = n_mon / Sat

# Number of molecules in a critical cluster - ceiling
g_star_ceil = np.int_(np.ceil(g_star))

# Initialise nucleation arrays
Nk_star = np.zeros((50))        # Number of seeds per cluster size (capped at 50 atoms per nucleus)
Jk_star = [0.,0.]               # Instantaneous nucleation rate (two zeros needed for the velocity integration)
Nk_star_corr = [0.]             # Number of seeds along the axis (normalised to clusters of size 8)
Nk_star_cum = [0.]              # Cumulative number of seeds 

# Save output
lines = []
lines.append(f"{x[0]},0.0\n")

#########################
# Start actual nucleation estimation
# Open the nucleation output file
with open(case + "-nuclest.csv", "w") as file_c:
    
    # Write the x coordinate
    file_c.write(f"{x[0]},0.0\n")

    i=0
    j=0
    # Iterate over x axis:
    for x_i in x[1:]:
        
        i+=1
        # Stop the simulation if the array is decreasing (last three terms decreasing) and term is smaller than 1 promile of the maximum value
        if Jk_star[-2] > Jk_star[-1]: 
            if Jk_star[-3] > Jk_star[-2]:
                if Jk_star[-1] < (0.001* max(Jk_star)):
                   lines.append(f"{x_i},{Nk_star_cum[-1]}\n")
                   continue
        
        j+=1
        # Instantaneous nucleation rate (classical nucleation)
        Jk_star_cur = n_sat[i] * n_mon[i] * v1 * (2.0 * sigma[i] / PI / m1)**0.5 * np.exp(theta[i] - (4.0*theta[i]**3.0/27.0/(np.log(Sat[i]))**2))
#            Jk_star_cur = n_sat[i] * n_mon[i] * v1 * (2.0 * sigma_tolman[i] / PI / m1)**0.5 * np.exp(theta_tolman[i] - (4.0*theta_tolman[i]**3.0/27.0/(np.log(Sat[i]))**2)) # alteranative using tolman correction
        # Add to the array
        Jk_star.append(Jk_star_cur)

        # Derive the time change from delta x and velocity
        time = (x_i - x[i-1])/(0.5*(v[i]+v[i-1]))
        
        # Calculate nuclei if size between 1 and 50 atoms
        if g_star[i] <= 50 and g_star[i] > 1.0:

            # Number of nuclei (and append to the array)
            Nk_star[g_star_ceil[i]-1] += Jk_star_cur * g_star[i]/g_star_ceil[i] * time
            Nk_star_corr.append(Jk_star_cur * g_star[i]/g_star_ceil[i]/8. * time)
        
            # Correct the number of monomers and nucleation size in the next timestep (depletion of gas Fe)
            n_mon[i+1] -= Jk_star_cur * g_star[i] * time 
            Sat[i+1] = n_mon[i+1]/n_sat[i+1]
            g_star[i+1] = (2.0*theta[i+1]/3.0/np.log(Sat[i+1]))**3
            g_star_ceil[i+1] = np.int_(np.ceil(g_star[i+1]))
        
        else:
            # Append no nucleation
            Nk_star_corr.append(0)
        
        # Calculate cumulative nucleation seeds
        Nk_star_cum.append(Nk_star_cum[-1] + Nk_star_corr[-1])
        
        # Write results
        lines.append(f"{x[i]},{Nk_star_cum[-1]}\n")

# Remove the first artificial 0
Jk_star.pop(0)

# Determine the starting position (at % of maximum nucleation rate)
with open(case_base_name) as f:
    for line in f:
        line = line.split("#", 1)[0].strip()
        if line.startswith("nucl_start_pos"):
            start_pos_scale = float(line.split()[1])
            break
threshold = start_pos_scale * max(Jk_star)
for idx,val in enumerate(Jk_star):
    if val > threshold:
        print(f"\nSimulation start position at {start_pos_scale*100}% of maximum nucleation rate:\n\t\t\t{x[idx]} mm (Index: {idx})\n")
        break

# Open the nucleation output file and write results
with open(case + "_results/" + case + "-nuclest.csv", "w") as file_c:
    file_c.write(f"# Start position index: {idx} \n")
    file_c.writelines(lines)
 
# Plot nucleation estimate
ax1.plot(x*1e3, g_star_ceil, '-', color=colour,linewidth=6)
bx1.plot(x[:j+1]*1e3, n_mon[:j+1], '-', color=colour, linewidth=6)
cx1.plot(x[:j+1]*1e3, Nk_star_corr, '-', color=colour, linewidth=6)
dx1.plot(x[:j+1]*1e3, Nk_star_cum, '-', color=colour, linewidth=6)

fig.tight_layout()
plt.show()
plt.close()
fig.savefig(f'{case}_results/nuclest_{case}.png')
