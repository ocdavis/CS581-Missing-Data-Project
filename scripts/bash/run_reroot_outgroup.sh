#!/usr/bin/env bash

set -u

parent_dir="${1:-.}"

find "$parent_dir" -mindepth 1 -maxdepth 1 -type d -print0 |
while IFS= read -r -d '' subdir; do
    input_file="$subdir/snaq-astral-iqtree-1.nwk"
    output_file="$subdir/reroot-snaq-astral-iqtree-1.nwk"

    if [[ -f "$input_file" ]]; then
        echo "Processing: $input_file"
        julia /u/quadros3/research/phylogeny/project_01/experiment_01/1a/scripts/julia/reroot_outgroup.jl "$input_file" OUT > "$output_file"
    else
        echo "Skipping $subdir: $input_file not found" >&2
    fi
done
