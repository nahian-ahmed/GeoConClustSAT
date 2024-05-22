#########################################################################################

# Geospatial Constrained Clustering as a Boolean Satisfiability (SAT) Problem
#
# May 22, 2024
#
# This file contains functions for performing geospatial constrained clustering using 
# SAT solvers. It includes functions to generate constraints in DIMACS CNF format, 
# call the Z3 solver, and parse the solver output. 
# The main function (geo_con_clust_sat()) integrates these components
# to perform the clustering.

#########################################################################################

# Load necessary libraries
library (dplyr, quietly = T)

# Source necessary R scripts
source("helpers.R")


# Function to generate partitioning constraints
#
# Description:
#   This function generates SAT constraints for partitioning a data frame into k partitions.
#   It ensures that each row is assigned to exactly one partition (coverage constraint),
#   that no row is assigned to more than one partition (no-overlap constraint),
#   and that each cluster has at least one instance assigned to it (non-empty cluster constraint).
#
# Inputs:
#   df : A data frame containing the data to be partitioned.
#   k  : An integer representing the number of partitions.
#
# Outputs:
#   A list containing:
#     - clauses : A list of clauses representing the partitioning constraints.
#     - num_vars: An integer representing the total number of variables (nrow(df) * k).
add_partitioning_constraints <- function(df, k) {
    num_vars <- nrow(df) * k  # Total number of variables (X_{i,j})
    total_vars <- num_vars  # Total number of variables (no auxiliary variables needed)

    clauses <- list()
    clause_idx <- 1

    # Coverage constraint: Ensure each row is assigned to exactly one partition.
    for (i in seq_len(nrow(df))) {
        clauses[[clause_idx]] <- (1:k) + (i - 1) * k
        clause_idx <- clause_idx + 1
    }

    # No-overlap constraint: Ensure no row is assigned to more than one partition.
    for (i in seq_len(nrow(df))) {
        for (j1 in 1:(k-1)) {
            for (j2 in (j1+1):k) {
                clauses[[clause_idx]] <- c(-((i - 1) * k + j1), -((i - 1) * k + j2))
                clause_idx <- clause_idx + 1
            }
        }
    }

    # Non-empty cluster constraint: Ensure each cluster contains at least one instance.
    for (j in 1:k) {
        clause <- integer(nrow(df))
        for (i in seq_len(nrow(df))) {
            clause[i] <- (i - 1) * k + j
        }
        clauses[[clause_idx]] <- clause
        clause_idx <- clause_idx + 1
    }

    return(list(clauses = clauses, num_vars = total_vars))
}




# Function to add instance-level constraints
#
# Description:
#   This function generates a list of clauses representing instance-level 
#   must-link and cannot-link constraints for clustering. It ensures that 
#   specified pairs of instances are either in the same cluster (must-link) 
#   or in different clusters (cannot-link) according to the constraints provided.
#
# Inputs:
#   df            : A data frame containing the data.
#   k             : An integer representing the number of partitions.
#   constraints_df: A data frame containing must-link and cannot-link constraints.
#   id_col        : A string representing the name of the column in df that contains unique identifiers.
#
# Outputs:
#   A list of clauses representing the instance-level constraints.
add_instance_level_constraints <- function(df, k, constraints_df, id_col) {
    clauses <- list()
    if (nrow(constraints_df) == 0) {
        return(clauses)
    }
    
    for (i in seq_len(nrow(constraints_df))) {
        instance1 <- constraints_df[i, "id_col_1"]
        instance2 <- constraints_df[i, "id_col_2"]
        idx1 <- which(df[[id_col]] == instance1)
        idx2 <- which(df[[id_col]] == instance2)
        
        if (length(idx1) == 0 || length(idx2) == 0) {
            stop(paste("Instance IDs", instance1, "or", instance2, "not found in df"))
        }
        
        if (constraints_df[i, "constraint"] == "must-link") {
            for (j in 1:k) {
                # Must-Link constraints
                # If two instances \(i_r\) and \(i_s\) are must-linked, they must be in the same cluster.
                # Therefore, for each cluster \(j\):
                # Clause 1: If instance1 is not in cluster \(j\), then instance2 must be in cluster \(j\).
                # \((\neg X_{i_r, j} \lor X_{i_s, j})\)
                clauses <- append(clauses, list(c(-((idx1 - 1) * k + j), (idx2 - 1) * k + j)))
                # Clause 2: If instance2 is not in cluster \(j\), then instance1 must be in cluster \(j\).
                # \((\neg X_{i_s, j} \lor X_{i_r, j})\)
                clauses <- append(clauses, list(c(-((idx2 - 1) * k + j), (idx1 - 1) * k + j)))
            }
        } else if (constraints_df[i, "constraint"] == "cannot-link") {
            for (j in 1:k) {
                # Cannot-Link constraints
                # If two instances \(i_r\) and \(i_s\) are cannot-linked, they must not be in the same cluster.
                # Therefore, for each cluster \(j\):
                # Clause 1: If instance1 is in cluster \(j\), then instance2 must not be in cluster \(j\).
                # \((X_{i_r, j} \rightarrow \neg X_{i_s, j})\)
                # \((\neg X_{i_r, j} \lor \neg X_{i_s, j})\)
                clauses <- append(clauses, list(c(-((idx1 - 1) * k + j), -((idx2 - 1) * k + j))))
            }
        }
    }
    
    return(clauses)
}


# Function to add cluster-level constraints
#
# Description:
#   This function generates a list of clauses representing cluster-level 
#   constraints for clustering based on pairwise distances. It ensures that 
#   pairs of instances either are not in the same cluster (max constraint) 
#   or must be in the same cluster (min constraint) according to the given 
#   threshold and constraint type.
#
# Inputs:
#   df            : A data frame containing the data.
#   k             : An integer representing the number of partitions.
#   pairwise_dist : A numeric matrix representing the pairwise distances between instances.
#   threshold     : A numeric value representing the distance threshold for constraints.
#   max_constraint: A boolean indicating whether to use max constraint (TRUE) or min constraint (FALSE).
#
# Outputs:
#   A list of clauses representing the cluster-level constraints.
add_cluster_level_constraints <- function(df, k, pairwise_dist, threshold, max_constraint = TRUE) {
    clauses <- list()
    
    # Determine indices based on the threshold and max_constraint
    if (max_constraint) {
        indices <- which(pairwise_dist > threshold, arr.ind = TRUE)
    } else {
        indices <- which(pairwise_dist < threshold, arr.ind = TRUE)
    }
    
    # Filter out pairs where i >= j to avoid redundancy, selecting the lower triangle of the matrix
    indices <- indices[indices[, 1] < indices[, 2], , drop = FALSE]
    
    # Preallocate list for clauses
    num_pairs <- nrow(indices)
    
    # If there are no pairs, return an empty list
    if (num_pairs == 0) {
        return(clauses)
    }
    
    num_clauses <- num_pairs * 2 * k  # Two clauses per pair per cluster
    clauses <- vector("list", num_clauses)
    
    clause_idx <- 1
    for (idx in 1:num_pairs) {
        i <- indices[idx, 1]
        j <- indices[idx, 2]
        for (cluster in 1:k) {
            if (max_constraint) {
                # Maximum Diameter Constraints (Cannot-Link)
                # If two instances \(i_r\) and \(i_s\) are cannot-linked, they must not be in the same cluster.
                # Therefore, for each cluster \(j\):
                # Clause: If instance1 is in cluster \(j\), then instance2 must not be in cluster \(j\).
                # \((\neg X_{i_r, j} \lor \neg X_{i_s, j})\)
                clauses[[clause_idx]] <- c(-((i - 1) * k + cluster), -((j - 1) * k + cluster))
                clause_idx <- clause_idx + 1
            } else {
                # Minimum Separation Constraints (Must-Link)
                # If two instances \(i_r\) and \(i_s\) are must-linked, they must be in the same cluster.
                # Therefore, for each cluster \(j\):
                # Clause 1: If instance1 is not in cluster \(j\), then instance2 must be in cluster \(j\).
                # \((\neg X_{i_r, j} \lor X_{i_s, j})\)
                clauses[[clause_idx]] <- c(-((i - 1) * k + cluster), (j - 1) * k + cluster)
                clause_idx <- clause_idx + 1
                # Clause 2: If instance1 is in cluster \(j\), then instance2 must be in cluster \(j\).
                # \((X_{i_r, j} \lor \neg X_{i_s, j})\)
                clauses[[clause_idx]] <- c(((i - 1) * k + cluster), -((j - 1) * k + cluster))
                clause_idx <- clause_idx + 1
            }
        }
    }
    
    return(clauses)
}


# Function to generate DIMACS CNF format
#
# Description:
#   This function generates a DIMACS CNF formatted string from a list of clauses 
#   representing the partitioning constraints.
#
# Inputs:
#   clauses  : A list of clauses representing the partitioning constraints.
#   num_vars : An integer representing the total number of variables.
#
# Outputs:
#   A string in DIMACS CNF format representing the clauses.
generate_dimacs <- function(clauses, num_vars) {
    
    # Preallocate space for all the clauses and the header
    cnf <- vector("character", length(clauses) + 1)
    cnf[1] <- paste("p cnf", num_vars, length(clauses))

    
    
    # Efficiently concatenate each clause
    cnf[2:(length(clauses) + 1)] <- sapply(clauses, function(clause) {
        paste(paste(clause, collapse = " "), "0")
    })
    
    # Combine all parts into a single string
    return(paste(cnf, collapse = "\n"))
}

# Function to write CNF file
#
# Description:
#   This function writes a given DIMACS CNF formatted string to a specified file path.
#
# Inputs:
#   cnf       : A string in DIMACS CNF format representing the clauses.
#   file_path : A string representing the path where the CNF file will be written.
#
# Outputs:
#   None
write_cnf_file <- function(cnf, file_path) {
    writeLines(cnf, file_path)
}

# Function to call Z3 solver with parallel processing
#
# Description:
#   This function calls the Z3 solver on a specified CNF file and captures the output.
#
# Inputs:
#   cnf_file : A string representing the path to the CNF file to be solved.
#
# Outputs:
#   A string representing the output from the Z3 solver.
call_z3_solver <- function(cnf_file) {
     
    # Construct the Z3 command with parallel processing
    z3_command <- paste("z3 -dimacs", cnf_file)
    
    # Execute the Z3 command
    result <- system(z3_command, intern = TRUE)
    
   
    return(result)
}

# Function to parse Z3 output and generate clustering
#
# Description:
#   This function parses the output from the Z3 solver to generate a clustering solution.
#   It checks if the problem is satisfiable, extracts the solution, and assigns each instance
#   in the data frame to a cluster based on the solution.
#
# Inputs:
#   output: A character vector containing the output from the Z3 solver.
#   df     : A data frame containing the instances to be clustered.
#   k      : An integer representing the number of clusters.
#
# Outputs:
#   A numeric vector of length nrow(df) where each element represents the cluster assigned to the corresponding instance.
parse_z3_output <- function(output, df, k) {
    # Check if the first line indicates a satisfiable result
    if (output[1] != "sat") {
        stop("The problem is unsatisfiable")
    }
    
    # Extract the solution line and convert it to numeric
    solution <- unlist(strsplit(output[2], "\\s+"))
    solution <- as.numeric(solution[solution != "" & solution != "0"])
    
    # Initialize clustering with NA values
    clustering <- rep(NA, nrow(df))
    
    # Process each solution value
    for (sol in solution) {
        if (sol > 0) {
            # Adjust for 1-based indexing in R
            var <- sol - 1
            instance <- (var %/% k) + 1
            cluster <- (var %% k) + 1
            
            # Assign cluster to instance
            if (instance <= nrow(df)) {
                clustering[instance] <- cluster
            } else {
                warning(paste("Instance index out of bounds:", instance))
            }
        }
    }
    
    # Check for unassigned instances
    if (any(is.na(clustering))) {
        stop("Some instances are not assigned to any cluster")
    }
    
    return(clustering)
}


# Main function for geospatial constrained clustering as a SAT problem
#
# Description:
#   This function performs geospatial constrained clustering using a SAT solver. 
#   It generates constraints based on the given data, objective function, and parameters, 
#   and then uses a SAT solver to find the optimal clustering solution.
#
# Inputs:
#   df            : DataFrame, dataset of instances.
#   k             : Integer, number of clusters.
#   sp_features   : Vector of strings, names of columns used for calculating geospatial distances.
#   env_features  : Vector of strings, names of columns used for calculating attribute distances.
#   rho           : Numeric, weighting parameter.
#   constraints_df: DataFrame, contains the must-link and cannot-link constraints.
#   f             : String, objective function to be optimized ("min_diameter" or "max_separation").
#   id_col        : String, name of the column containing instance IDs.
#
# Outputs:
#   List containing:
#     - result       : DataFrame with two columns (id and cluster) indicating the cluster assignments.
#     - optimal_value: Numeric, the optimal value of the objective function.
#     - call_times   : DataFrame with information about each SAT solver call, including iteration number,
#                      number of clauses, SAT/UNSAT output, and time taken.
geo_con_clust_sat <- function(df, k, sp_features, env_features, rho, constraints_df, f, id_col) {
    # Create pairwise distance matrix
    pairwise_dist <- create_pairwise_distance_matrix(df, sp_features, env_features, rho)
    
    # Create a sorted array of distinct pairwise distances
    distinct_distances <- sort(unique(as.vector(pairwise_dist)))
    
    # Initialize binary search bounds
    low <- 1
    high <- length(distinct_distances)
    
    optimal_value <- NULL
    optimal_clustering <- NULL
    
    # Generate partitioning constraints
    partitioning_constraints <- add_partitioning_constraints(df, k)
    clauses <- partitioning_constraints$clauses
    num_vars <- partitioning_constraints$num_vars

    # Generate instance-level constraints
    instance_constraints <- add_instance_level_constraints(df, k, constraints_df, id_col)
    clauses <- append(clauses, instance_constraints)
    
    iter <- 1
    call_times <- list() # Initialize the list to store call times
    
    while (low <= high) {
        mid <- floor((low + high) / 2)
        threshold <- distinct_distances[mid]
        
        cat("\tIteration ", iter,". high: ", high, ", low: ", low, ", L[mid]: ", threshold, "\n")
        
        clauses_iter <- clauses

        # Generate cluster-level constraints
        if (f == "min_diameter") {
            cluster_clauses <- add_cluster_level_constraints(df, k, pairwise_dist, threshold, max_constraint = TRUE)
        } else if (f == "max_separation") {
            cluster_clauses <- add_cluster_level_constraints(df, k, pairwise_dist, threshold, max_constraint = FALSE)
        }
        
        clauses_iter <- append(clauses_iter, cluster_clauses)
        
        cat("\t\tall constraints generated\n")
        cat("\t\tno. of clauses: ", length(clauses_iter),"\n")
        
        # Generate DIMACS CNF
        cnf <- generate_dimacs(clauses_iter, num_vars)
        cat("\t\tDIMACS CNF generated\n")

        # Write CNF to file
        cnf_file <- tempfile(fileext = ".cnf")
        write_cnf_file(cnf, cnf_file)
        cat("\t\tCNF file written\n")

        cat("\t\tcalling solver now..\n")
        
        # Call Z3 solver and measure time taken
        start_time <- Sys.time()
        output <- call_z3_solver(cnf_file)
        end_time <- Sys.time()
        time_taken <- as.numeric(difftime(end_time, start_time, units = "secs"))
        cat("\t\tsolver returned output\n")
        
        # Check satisfiability
        if (grepl("^sat", output[1])) {
            cat("\t\tSATISFIABLE\n\n")
            sat_output <- "SAT"
            optimal_value <- threshold
            optimal_clustering <- parse_z3_output(output, df, k)
            
            
            if (f == "min_diameter") {
                high <- mid - 1
            } else if (f == "max_separation") {
                low <- mid + 1
            }
        } else {
            cat("\t\tUNSATISFIABLE\n\n")
            sat_output <- "UNSAT"
            if (f == "min_diameter") {
                low <- mid + 1
            } else if (f == "max_separation") {
                high <- mid - 1
            }
        }
        
        # Record the iteration number, SAT output, and time taken
        call_times[[iter]] <- list(iteration = iter, num_clauses = length(clauses_iter), sat_output = sat_output, time_taken = time_taken)
        
        iter <- iter + 1
    }
    
    result <- data.frame(id = df[, id_col], cluster = optimal_clustering)
    
    # Convert call_times list to dataframe
    call_times_df <- do.call(rbind, lapply(call_times, as.data.frame))
    
    return(list(result = result, optimal_value = optimal_value, call_times = call_times_df))
}


