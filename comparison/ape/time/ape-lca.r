#!/usr/bin/env Rscript

# Load the necessary library
library(ape)
options(warn=-1)

# Get input tree file
# Check if the correct number of arguments is provided
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: ape-mrca.R <tree1_file>")
}

treefile <- args[1]

# Read the trees
tree1 <- read.tree(treefile)
# tree2 <- read.tree(tree2_file)

# Start the clock
start <- Sys.time()

tips <- c("Tip60", "Tip126", "Tip118", "Tip161", "Tip25", "Tip127", "Tip38", "Tip186", "Tip21", "Tip56", "Tip78", "Tip70", "Tip115")

mrca(tree1, full = FALSE)
out <- getMRCA(tree1, tips)

# Stop the clock
end <- Sys.time()
duration <- difftime(end, start, units="secs")
print(paste("Internal time:", duration))