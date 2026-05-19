#!/usr/bin/env bash

set -euo pipefail

# Run the main code on X number of cores
ITERATIONS=30

# Number of different part_params files
FILES_NO=250

# Path to the main file
PYTHON_SCRIPT="SRC/SPMC-0D-main-200326.py"

########################################
# Step 1: Find *.case files
########################################
shopt -s nullglob
case_files=( *.case )
shopt -u nullglob

########################################
# Step 2: Validate number of case files
########################################
count=${#case_files[@]}

if [[ "$count" -eq 0 ]]; then
    echo "Error: No .case files found in the directory." >&2
    exit 1
elif [[ "$count" -gt 1 ]]; then
    echo "Error: Multiple .case files found. Ensure only one file exists." >&2
    printf '%s\n' "${case_files[@]}"
    exit 1
fi

########################################
# Step 3: Use the single case file
########################################
case_file="${case_files[0]}"
echo "Found case file: $case_file"

########################################
# Step 4: Determine results directory
########################################
case_name="${case_file%.case}"
results_dir="${case_name}_results"

# Loop over each part_params file
for ((j=246; j<=FILES_NO; j++))
do
    current_file=$(printf "sample-T_a-gaussian/sampled_yaml/collisions_sampled_%03d.yaml" "${j}")
    ln -sf "$current_file" part_params.yaml

    echo "Using part_params file no. ${j}"

    for ((i=31; i<=ITERATIONS+30; i++))
    do
	    taskset -c "$i" python "$PYTHON_SCRIPT" "$case_file" "$((i-30))" &
    done

    wait

    echo "All runs completed successfully."

    python SRC/runningAver.py $case_file $ITERATIONS

    cp -r "$results_dir" "res_part-params_${j}"

    shopt -s nullglob 
    rm "$results_dir"/*SPMC-[0-9].csv "$results_dir"/*SPMC-[1-3][0-9].csv "$results_dir"/*average.csv "$results_dir"/*count.csv 
    shopt -u nullglob 
done
