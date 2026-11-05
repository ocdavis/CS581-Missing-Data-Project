#!/bin/bash

if [ $# -ne 7 ]; then
    echo "usage: run_snaq_wrapper.sh [gene_trees_dir] [network_dir] [output_dir] [num_ret] [threads] [gene_tree_method] [network_method]"
    exit 1
fi

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

gene_tree_path="$1"
network_path="$2"
output_path="$3"
n_ret="$4"
threads="$5"
gene_tree_method="$6"
network_method="$7"

stats_path="${output_path}/snaq_stats.log"
mkdir -p "$output_path"

find "$network_path" -type f -name "$network_method" | while read network_file; do
    dir_path=$(dirname "$network_file")
    sub_path=$(realpath --relative-to="$network_path" "$dir_path")
    
    replicate_gene_tree_path="${gene_tree_path}/${sub_path}/${gene_tree_method}"
    
    if [ ! -f "$replicate_gene_tree_path" ]; then
        echo "$replicate_gene_tree_path does not exist, continuing."
        continue
    fi
    
    replicate_output_path="${output_path}/${sub_path}"
    mkdir -p "$replicate_output_path"
    
    file_name=$(basename "$network_file")
    f_name="${file_name%.*}"
    
    snaq_output="${replicate_output_path}/snaq_${f_name}.nwk"
    
    /usr/bin/time -a -f "${sub_path},${f_name},%e,%M" -o "$stats_path" \
        julia "${script_dir}/../julia/run_snaq.jl" "$replicate_gene_tree_path" "$network_file" "$snaq_output" "$n_ret" "$threads"
done

echo "Stats saved to $stats_path"
