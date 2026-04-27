# SPMC main

SPMC legacy code + 2026 Powder Technology paper extension: allows to artificially scale species concentration for sensitivity studies

## Description

This repository contains a variant of the SPMC (Single Particle Monte Carlo) code, specifically modified to allow artificial scaling of species concentrations. This is useful for sensitivity studies and analysis of combustion or chemical reaction models, particularly in the context of nanoparticle formation and growth in flames.

The code simulates the growth of single iron oxide particles through Monte Carlo methods, coupled with Cantera for gas-phase chemistry and flame structure.

## Features

- Monte Carlo simulation of particle growth
- Coupling with Cantera for detailed gas-phase kinetics
- Support for different nucleation modes (add/decomp)
- Evaporation and radiation models
- Output in CSV format for analysis

## Files

- `SRC/`: Source code directory
  - `nuclest.py`: Nucleation estimation script
  - `SPMC-0D-main-200326.py`: Main SPMC 0D code
  - `runningAver.py`: Running average calculations
- `prerun.sh`: Runs the flame and nucleation estimation scripts
- `exe.sh`: SPMC execution (parallel)
- `L500_OHx25%.case`: Case file for main flame parameters (here: lean H2/O2/Ar, 500 ppm IPC impinging flame)
- `part_params.yaml`: Particle and Monte-Carlo parameters
- `flame.py`: Cantera flame simulation file 
- `Iron_elte_Syngas-newTransp.yaml`: Chemistry mechanism
- `T_of_x-L500_OHx25%.csv`: Temperature profile data from OF 2D sim., fed to cantera (optional, depending of the cantera flame file)

## Installation

1. Ensure Python 3.x is installed.
2. Install required packages:
   ```bash
   pip install cantera numpy scipy pyyaml
   ```
3. Clone or download this repository.

## Usage

1. Adjust/replace your flame file acc. to your desired case. This file serves as a template.
2. Define the `<name>.case` file with basic flame parameters. The <name> of this file defines the whole simulation name. Set if any species should be artficially scaled. If flame requires some external data (T_of_x-)
3. Run `./prerun.sh` to execute your `flame.py` and `nuclest.py` to estimate the number of nuclei. This also creates <name>_results folder, where the results files (Cantera and nucleation estimation) are written. All case files are copied there. 
4. Define particle and Monte-Carlo parameters in `part_params.yaml`.
5. Specify number of cores in `exe.sh` for parallel runs.
6. Run `./exe.sh` to start (parallel) SPMC simulation(s). `Runningaverage.py` is run automatically afterwards. Results are written to the <name>_results folder.

## Requirements

- Python 3.x
- Cantera (for gas-phase chemistry)
- NumPy, SciPy, PyYAML
- Bash shell for scripts

## Contributing

Please create issues for bugs or feature requests. For contributions, create a pull request with a clear description of changes.

## License

[Add license information here]