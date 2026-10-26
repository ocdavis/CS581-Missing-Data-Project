#!/usr/bin/env julia
# Reroot an extended Newick network/tree at a given outgroup taxon and
# print the rerooted topology to stdout. The input file is never modified.
#
# Usage: julia reroot_outgroup.jl <newick_file> <outgroup_label>

using PhyloNetworks

function main(args)
    if length(args) != 2
        println(stderr, "Usage: julia reroot_outgroup.jl <newick_file> <outgroup_label>")
        exit(1)
    end
    newick_path, outgroup = args

    if !isfile(newick_path)
        println(stderr, "Error: file not found: $newick_path")
        exit(1)
    end

    net = readnewick(newick_path)

    try
        rootatnode!(net, outgroup)
    catch e
        println(stderr, "Error: could not root network at outgroup \"$outgroup\": $e")
        exit(1)
    end

    println(writenewick(net))
end

main(ARGS)
