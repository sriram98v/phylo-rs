#!/usr/bin/env Rscript

# Load the necessary library
library(ape)
options(warn=-1)


# Get input tree file
# Check if the correct number of arguments is provided
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: ape-rfs.R <tree1_file> <tree2_file>")
}

# Start the clock
start <- Sys.time()

# Read the input files from command-line arguments
tree1_file <- args[1]
tree2_file <- args[2]

# Read the trees
tree1 <- read.tree(tree1_file)
tree2 <- read.tree(tree2_file)

# Compute the unnormalized Robinson-Foulds distance
rf_distance <- dist.topo(tree1, tree2)

end <- Sys.time()

duration <- difftime(end, start, units="secs")
print(paste("", rf_distance))
print(paste("Internal time:", duration))

# Output the result
# cat("Unnormalized Robinson-Foulds distance between the trees:", rf_distance, "\n")
