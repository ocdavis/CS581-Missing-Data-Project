#!/usr/bin/env Rscript
# Infer a tree of blobs with ECToBlob (MSCquartets R package).
#
# ECToBlob needs two inputs: a collection of gene trees (one Newick tree per
# line) and a tree assumed to be a resolution of the true tree of blobs
# (e.g. a TREE-QMC species tree or an ASTRAL tree). It tests every quartet
# for fit to a star tree and to the resolution displayed on the starting
# tree, contracts unsupported edges, and returns a sequence of trees from
# which one is selected via its p-value (indexEarly or indexLate).
#
# Usage:
#   Rscript ectoblob.R --genes <gene trees> --tree <start tree> --out <file> \
#       [--alpha 0.05] [--beta 0.8] [--qtype bi|quad|mul] \
#       [--correction Bon|Cauchy|BBC|CBC] [--index early|late] \
#       [--outgroup taxon] [--save-table]

suppressPackageStartupMessages(library(MSCquartets))

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
    i <- match(flag, args)
    if (!is.na(i) && i < length(args)) return(args[i + 1])
    default
}

has_flag <- function(flag) !is.na(match(flag, args))

genes_file <- get_arg("--genes")
tree_file <- get_arg("--tree")
out_file <- get_arg("--out")
if (is.null(genes_file) || is.null(tree_file) || is.null(out_file)) {
    stop("usage: ectoblob.R --genes <gene trees> --tree <start tree> ",
         "--out <output tree> [--alpha a] [--beta b] [--qtype t] ",
         "[--correction c] [--index early|late] [--outgroup taxon] ",
         "[--save-table]")
}

alpha <- as.numeric(get_arg("--alpha", "0.05"))
beta <- as.numeric(get_arg("--beta", "0.8"))
qtype <- get_arg("--qtype", "bi")
correction <- get_arg("--correction", "Bon")
index <- get_arg("--index", "early")
outgroup <- get_arg("--outgroup")
save_table <- has_flag("--save-table")

gts <- read.tree(genes_file)
start_tree <- read.tree(tree_file)

out <- ECToBlob(gts, start_tree, alpha = alpha, beta = beta,
                qType = qtype, testCorrection = correction, plot = 0)

idx_name <- if (index == "early") "indexEarly" else "indexLate"
newick <- out$treeList[[out[[idx_name]]]]$tree

tree <- read.tree(text = newick)
if (!is.null(outgroup)) {
    if (!is.na(match(outgroup, tree$tip.label))) {
        tree <- root(tree, outgroup, resolve.root = FALSE)
    } else {
        warning("outgroup '", outgroup, "' not in tree; leaving unrooted")
    }
}
write.tree(tree, file = out_file)

if (save_table) {
    saveRDS(out$pTable, sub("\\.tre$", ".pTable.rds", out_file))
}
