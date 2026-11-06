#!/usr/bin/env bash

set -u

first_dir="${1:-.}"
second_dir="${2:-.}"

find "$first_dir" -mindepth 1 -maxdepth 1 -type d -print0 |
while IFS= read -r -d '' subdir; do
    # Get the basename of the subdirectory
    subdir_name=$(basename "$subdir")
    
    input_file="$subdir/true_net.nwk"
    output_dir="$second_dir/$subdir_name"
    output_file="$output_dir/true_net-rooted.nwk"

    if [[ -f "$input_file" ]]; then
        # Create output directory if it doesn't exist
        mkdir -p "$output_dir"
        echo "Processing: $input_file"
        julia /u/quadros3/research/phylogeny/project_01/experiment_01/1a/scripts/julia/reroot_outgroup.jl "$input_file" OUT > "$output_file"
    else
        echo "Skipping $subdir: $input_file not found" >&2
    fi
done

