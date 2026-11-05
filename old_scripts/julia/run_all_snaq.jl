#!/usr/bin/env julia
using Distributed
using PhyloNetworks, SNaQ, DataFrames

if length(ARGS) != 7
    print("usage: run_snaq.jl [gene_trees_dir] [network_dir] [output_dir] [num_ret] [threads] [gene_tree_method] [network_method]\n")
    exit(1)
end

gene_tree_path = ARGS[1]
network_path = ARGS[2]
output_path = ARGS[3]
n_ret = parse(Int, ARGS[4])
threads = parse(Int, ARGS[5])
gene_tree_method = ARGS[6]
network_method = ARGS[7]

addprocs(threads)
@everywhere using PhyloNetworks, SNaQ, DataFrames

function infer_all_snaq(network_path, gene_tree_path, output_path, n_ret, threads, gene_tree_method, network_method)
    for (root, dirs, files) in walkdir(network_path)
        sub_path = relpath(root, network_path)
        replicate_gene_tree_path = joinpath(gene_tree_path, sub_path, gene_tree_method)
        
        if !isfile(replicate_gene_tree_path)
            println("$replicate_gene_tree_path does not exist, continuing.")
            continue
        end
        
        replicate_output_path = joinpath(output_path, sub_path)
        mkpath(replicate_output_path)
        
        for file in files
            if file == network_method
                network_file = joinpath(root, file)
                
                # Read gene trees and network
                gtrees = readmultinewick(replicate_gene_tree_path)
                start_tree = readnewick(network_file)
                
                # Calculate quartet counts
                q, t = countquartetsintrees(gtrees)
                cf = readtableCF(DataFrame(tablequartetCF(q, t), copycols=false))
                
                f_name = splitext(file)[1]
                output_file = joinpath(replicate_output_path, "snaq_$(f_name).nwk")
                
                # Run SNaQ
                result_net = snaq!(start_tree, cf, hmax=n_ret, filename=output_file, seed=0, runs=threads)
                writenewick(result_net, output_file)
                
                break
            end
        end
    end
end

infer_all_snaq(network_path, gene_tree_path, output_path, n_ret, threads, gene_tree_method, network_method)
