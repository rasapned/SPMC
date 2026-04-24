#!/usr/bin/env bash

set -euo pipefail

# Run the main code on X number of cores
ITERATIONS=1

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

########################################
# Step 5: Ensure results directory exists
########################################
if [[ ! -d "$results_dir" ]]; then
    echo "Error: Results directory '$results_dir' does not exist. Aborting." >&2
    exit 1
fi

echo "Using existing results directory: $results_dir"

########################################
# Step 6: Run iterations in parallel
########################################
pids=()
for ((i=1; i<=ITERATIONS; i++))
do
    taskset -c "$i" python "$PYTHON_SCRIPT" "$case_file" "$i" &
    pids+=($!)
done

# Wait for all processes and check exit codes
all_success=true
for pid in "${pids[@]}"; do
    if ! wait "$pid"; then
        echo "Error: Process $pid failed" >&2
        all_success=false
    fi
done

if $all_success; then
    echo "All runs completed successfully."
else
    echo "Some runs failed. Check output above." >&2
    exit 1
fi

########################################
# Step 7: Averaging the results and statistics
########################################
#python SRC/runningAverage.py $case_file $ITERATIONS
