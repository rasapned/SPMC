# SPMC-varGasSpec

SPMC branch: allows to artificially scale species concentration for sensitivity studies

## Description

This repository contains a variant of the SPMC (likely Single Particle Monte Carlo) code, specifically modified to allow artificial scaling of species concentrations. This is useful for sensitivity studies and analysis of combustion or chemical reaction models.

## Files

- `exe.sh`: Execution script
- `flame.py`: Main Python script for flame simulation
- `Iron_elte_Syngas-newTransp.yaml`: Configuration file for iron element with syngas and new transport properties
- `L500_OHx25%.case`: Case file for L500 with 25% OH
- `part_params.yaml`: Particle parameters configuration
- `prerun.sh`: Pre-run setup script
- `T_of_x-L500_OHx25%.csv`: Temperature profile data
- `SRC/`: Source code directory
  - `nuclest.py`: Nucleation script
  - `runningAver.py`: Running average calculations
  - `SPMC-0D-main-200326.py`: Main SPMC 0D code

## Usage

1. Run `prerun.sh` to set up the environment
2. Execute `exe.sh` to run the simulation
3. Analyze results from output files

## Requirements

- Python 3.x
- Required Python packages (install via pip or conda)
- Bash shell for scripts

## Contributing

Please create issues for bugs or feature requests. For contributions, create a pull request with a clear description of changes.

## License

[Add license information here]