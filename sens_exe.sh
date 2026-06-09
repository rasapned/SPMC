#!/usr/bin/env bash

set -euo pipefail

# Run the main code on X number of cores
ITERATIONS=1

# Number of different part_params files
FILES_NO=2

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

mkdir -p "RESULTS/${case_file%.case}_results_sens"

# Loop over each part_params file
for ((j=1; j<=FILES_NO; j++))
do
    current_file=$(printf "SRC/scripts/sample-T_a-gaussian/sampled_yaml/collisions_sampled_%03d.yaml" "${j}")
    ln -sf "$current_file" part_params.yaml

    echo "Using part_params file no. ${j}"

    for ((i=0; i<=ITERATIONS; i++))
    do
	    taskset -c "$i" python "$PYTHON_SCRIPT" "$case_file" "$((i))" &
    done

    wait

    echo "All runs completed successfully."

    python SRC/runningAver.py $case_file $ITERATIONS

    cp -r "RESULTS/$results_dir" "RESULTS/${case_file%.case}_results_sens/mech_${j}"
    #mv "RESULTS/${case_file%.case}_res_part-params/$results_dir" "RESULTS/${case_file%.case}_res_part-params/res_part-params_${j}"

    shopt -s nullglob
    # Remove RESULTS > case-name_results
    rm "RESULTS/$results_dir"/*SPMC-[0-9].csv "RESULTS/$results_dir"/*SPMC-[1-3][0-9].csv "RESULTS/$results_dir"/*average.csv "RESULTS/$results_dir"/*count.csv
  
    # Remove RESULTS > case-name_res_part-params > res_part-params_*
    rm "RESULTS/${case_name}_results_sens/mech_${j}"/*SPMC-[0-9].csv "RESULTS/${case_name}_results_sens/mech_${j}"/*SPMC-[1-3][0-9].csv 
    
    shopt -u nullglob
done
