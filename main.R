#########################################################################################

# Geospatial Constrained Clustering as a Boolean Satisfiability (SAT) Problem

# Nahian Ahmed (ahmedna@oregonstate.edu) and Kazi Ahmed Asif Fuad (fuadk@oregonstate.edu)
# School of EECS, Oregon State University

# May 22, 2024

# This file runs experiments on geospatial clustering using various constraints and parameters.
# It sets up the working directory, seeds for reproducibility, and sources the necessary 
# R scripts for running experiments and plotting results. The experiments include generating 
# and verifying clustering results based on different types of constraints (must-link, 
# cannot-link, and mixed) and evaluating these results using various metrics.

#########################################################################################

# Set working directory to current working directory
setwd("~/Documents/Research/code/GeoConClustSAT/")

# Set the seed for reproducibility
r_seed <- 1
set.seed(r_seed)

# Run all experiments
source("experiments.R")

# Plot results
source("plot.R")
