#!/usr/bin/env bash

set -euo pipefail

# Step 1: Find *.case files
shopt -s nullglob
case_files=( *.case )
shopt -u nullglob

# Step 2: Check how many were found
count=${#case_files[@]}

if [[ "$count" -eq 0 ]]; then
    echo "Error: No .case files found in the directory." >&2
    exit 1
elif [[ "$count" -gt 1 ]]; then
    echo "Error: Multiple .case files found. Ensure only one file exists." >&2
    printf '%s\n' "${case_files[@]}"
    exit 1
fi

# Step 3: Exactly one file
case_file="${case_files[0]}"
echo "Found case file: $case_file"

# Step 4: create results folder based on filename (without extension)
case_name="${case_file%.case}"
results_dir="${case_name}_results"

# Check if directory exists
if [[ -d "RESULTS/$results_dir" ]]; then
    read -rp "Directory 'RESULTS/$results_dir' already exists. Continue and possibly overwrite contents? [y/N]: " answer
    case "$answer" in
        [yY]|[yY][eE][sS])
            echo "Continuing..."
            ;;
        *)
            echo "Aborting."
            exit 1
            ;;
    esac
else
    mkdir -p "RESULTS/$results_dir"
    echo "Created results directory: RESULTS/$results_dir"
fi

# Step 5: Run Python scripts: solve Cantera flame and estimate number of seeds
python flame.py $case_file
echo "Flame results file saved!"
python SRC/nuclest.py "$case_file"
echo "Number of seeds estimated and saved to file!"

# Step 6: Copy the case file, T-of-x and so on
cp $case_file "RESULTS/$results_dir"
cp T_of_x* "RESULTS/$results_dir"
