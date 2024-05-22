#########################################################################################

# Geospatial Constrained Clustering as a Boolean Satisfiability (SAT) Problem

# May 22, 2024

# This file contains functions for:
# - Generating must-link and cannot-link constraints for geospatial constrained clustering
# - Sampling data while maintaining cluster integrity
# - Creating pairwise distance matrices combining geospatial and attribute distances
# - Verifying the validity of clustering results against various constraints
# - Evaluating clustering results using metrics such as Adjusted Rand Index (ARI),
#   Adjusted Mutual Information (AMI), and Normalized Information Distance (NID)

#########################################################################################

library (dplyr, quietly = T)
library (mclust, quietly = T)
library (aricode, quietly = T)


# Function to generate must-link constraints
#
# Description:
#   This function generates a specified number of must-link constraints for a dataset, 
#   ensuring that pairs of instances from the same cluster are linked together.
#
# Inputs:
#   df     : DataFrame, the dataset to analyze.
#   id_col : String, the name of the ID column.
#   gt_col : String, the name of the ground truth cluster column.
#   n_con  : Integer, the number of must-link constraints to generate.
#
# Outputs:
#   DataFrame, a dataframe with the must-link constraints.
generate_must_link <- function(df, id_col, gt_col, n_con) {
    # Initialize an empty data frame to store the must-link constraints
    must_link <- data.frame(id_col_1 = character(),
                            id_col_2 = character(),
                            constraint = character(),
                            stringsAsFactors = FALSE)
    
    # Split the dataset into clusters based on the ground truth column
    clusters <- split(df, df[[gt_col]])
    
    # Calculate the maximum possible number of must-link constraints
    max_possible_con <- sum(sapply(clusters, function(cluster) choose(nrow(cluster), 2)))
    
    # Check if the requested number of must-link constraints is feasible
    if (n_con > max_possible_con) {
        stop(paste("Error: Cannot generate", n_con, "must-link constraints. Maximum possible is", max_possible_con))
    }
    
    # Generate must-link constraints until the desired number is reached
    while (nrow(must_link) < n_con) {
        for (cluster in clusters) {
            if (nrow(cluster) > 1) {
                # Generate all possible pairs of instances within the cluster
                pairs <- t(combn(cluster[[id_col]], 2))
                
                # Create a data frame with the pairs and the "must-link" label
                pairs <- data.frame(pairs, constraint = "must-link")
                colnames(pairs) <- c("id_col_1", "id_col_2", "constraint")
                
                # Add the new pairs to the must-link constraints data frame
                must_link <- bind_rows(must_link, pairs)
                
                # Remove duplicate constraints
                must_link <- must_link[!duplicated(must_link), ]
                
                # Stop if the desired number of constraints is reached
                if (nrow(must_link) >= n_con) break
            }
        }
    }
    
    # Return the specified number of must-link constraints
    return(must_link[1:n_con, ])
}


# Function to generate cannot-link constraints
#
# Description:
#   This function generates a specified number of cannot-link constraints for a dataset, 
#   ensuring that pairs of instances from different clusters are not linked together.
#
# Inputs:
#   df     : DataFrame, the dataset to analyze.
#   id_col : String, the name of the ID column.
#   gt_col : String, the name of the ground truth cluster column.
#   n_con  : Integer, the number of cannot-link constraints to generate.
#
# Outputs:
#   DataFrame, a dataframe with the cannot-link constraints.
generate_cannot_link <- function(df, id_col, gt_col, n_con) {
    # Initialize an empty data frame to store the cannot-link constraints
    cannot_link <- data.frame(id_col_1 = character(),
                              id_col_2 = character(),
                              constraint = character(),
                              stringsAsFactors = FALSE)
    
    # Get the unique clusters in the dataset
    unique_clusters <- unique(df[[gt_col]])
    
    # Calculate the maximum possible number of cannot-link constraints
    max_possible_con <- choose(length(unique_clusters), 2) * choose(nrow(df) / length(unique_clusters), 2)
    
    # Check if the requested number of cannot-link constraints is feasible
    if (n_con > max_possible_con) {
        stop(paste("Error: Cannot generate", n_con, "cannot-link constraints. Maximum possible is", max_possible_con))
    }
    
    # Generate cannot-link constraints until the desired number is reached
    while (nrow(cannot_link) < n_con) {
        # Generate random pairs of instances
        pairs <- data.frame(id_col_1 = sample(df[[id_col]], n_con, replace = TRUE),
                            id_col_2 = sample(df[[id_col]], n_con, replace = TRUE),
                            stringsAsFactors = FALSE)
        
        # Remove pairs where the instances are the same
        pairs <- pairs[pairs$id_col_1 != pairs$id_col_2, ]
        
        # Remove duplicate pairs
        pairs <- pairs[!duplicated(pairs), ]
        
        # Label the pairs as "cannot-link"
        pairs$constraint <- "cannot-link"
        
        # Check if the instances in each pair belong to different clusters
        for (i in seq_len(nrow(pairs))) { 
            id1 <- pairs$id_col_1[i]
            id2 <- pairs$id_col_2[i]
            if (df[df[[id_col]] == id1, gt_col] != df[df[[id_col]] == id2, gt_col]) {
                # Add the valid cannot-link pair to the constraints data frame
                cannot_link <- bind_rows(cannot_link, pairs[i, ])
                
                # Stop if the desired number of constraints is reached
                if (nrow(cannot_link) >= n_con) break
            }
        }
    }
    
    # Return the specified number of cannot-link constraints
    return(cannot_link[1:n_con, ])
}


# Function to generate a constraints data frame for constrained clustering
#
# Description:
#   This function generates a specified number of must-link or cannot-link constraints, 
#   or a mixture of both, based on the given constraint type.
#
# Inputs:
#   df      : Data frame of instances.
#   id_col  : Column name for the unique identifier of each instance.
#   gt_col  : Column name for the ground truth cluster of each instance.
#   con_type: Type of constraint ("none", "must-link", "cannot-link", "mixed").
#   n_con   : Total number of constraints to generate.
#
# Outputs:
#   DataFrame, a dataframe with n_con rows and 3 columns:
#     id_col_1: ID of instance 1.
#     id_col_2: ID of instance 2.
#     constraint: Either "must-link" or "cannot-link".
get_constraints_dataframe <- function(df, id_col, gt_col, con_type, n_con) {
    # Generate constraints based on con_type
    if (con_type == "none") {
        constraints_df <- data.frame(id_col_1 = character(),
                                     id_col_2 = character(),
                                     constraint = character(),
                                     stringsAsFactors = FALSE)
    } else if (con_type == "must-link") {
        constraints_df <- generate_must_link(df, id_col, gt_col, n_con)
    } else if (con_type == "cannot-link") {
        constraints_df <- generate_cannot_link(df, id_col, gt_col, n_con)
    } else if (con_type == "mixed") {
        half_con <- n_con %/% 2
        must_link <- generate_must_link(df, id_col, gt_col, half_con)
        cannot_link <- generate_cannot_link(df, id_col, gt_col, n_con - half_con)
        constraints_df <- bind_rows(must_link, cannot_link)
    } else {
        stop("Invalid constraint type. Choose from 'none', 'must-link', 'cannot-link', 'mixed'.")
    }
    
    # Remove any rows with all NA values
    constraints_df <- constraints_df[!apply(is.na(constraints_df), 1, all), ]
    
    return(constraints_df)
}


# Function to sample a dataset while maintaining cluster integrity
#
# Description:
#   This function samples a specified percentage of instances from a dataset, ensuring 
#   that the integrity of the clusters is maintained.
#
# Inputs:
#   df         : DataFrame, the dataset to be sampled.
#   percent    : Numeric, the approximate percentage of the dataset to sample (between 0 and 100).
#   id_col     : String, the name of the ID column.
#   cluster_col: String, the name of the cluster column.
#
# Outputs:
#   DataFrame, the sampled dataset.
sample_df <- function(df, percent, id_col, cluster_col) {
    # Calculate the target sample size in terms of number of instances
    target_sample_size <- ceiling((percent / 100) * nrow(df))
    
    # Get unique clusters
    unique_clusters <- unique(df[[cluster_col]])
    
    # Initialize an empty dataframe to store the sampled data
    sampled_df <- data.frame()
    
    # Randomly shuffle the clusters
    shuffled_clusters <- sample(unique_clusters)
    
    # Initialize a counter for the number of sampled instances
    sampled_instance_count <- 0
    
    # Sample clusters until the target sample size is reached or exceeded
    for (cluster in shuffled_clusters) {
        cluster_data <- subset(df, df[[cluster_col]] == cluster)
        sampled_df <- rbind(sampled_df, cluster_data)
        sampled_instance_count <- sampled_instance_count + nrow(cluster_data)
        
        if (sampled_instance_count >= target_sample_size) {
            break
        }
    }
    
    # Return the sampled dataframe
    return(sampled_df)
}

# Function to get the number of unique clusters in a dataframe
#
# Description:
#   This function calculates and returns the number of unique clusters in a given dataset.
#
# Inputs:
#   df         : DataFrame, the dataset to analyze.
#   cluster_col: String, the name of the cluster column.
#
# Outputs:
#   Integer, the number of unique clusters.
get_true_k <- function(df, cluster_col) {
    # Calculate the number of unique clusters
    num_clusters <- length(unique(df[[cluster_col]]))
    
    return(num_clusters)
}

# Function to normalize a matrix
#
# Description:
#   This function normalizes the values in a given numeric matrix to the range [0, 1].
#
# Inputs:
#   mat : A numeric matrix to be normalized.
#
# Outputs:
#   A numeric matrix with values normalized to the range [0, 1].
normalize_matrix <- function(mat) {
    return((mat - min(mat)) / (max(mat) - min(mat)))
}

# Function to calculate geospatial distances
#
# Description:
#   This function calculates the pairwise Euclidean distances between instances in a data 
#   frame based on specified geospatial features.
#
# Inputs:
#   df          : A data frame containing the data.
#   sp_features : A vector of strings representing the names of columns used for geospatial distance calculation.
#
# Outputs:
#   A numeric matrix of pairwise Euclidean distances based on geospatial features.
calculate_geospatial_distance <- function(df, sp_features) {
    dist_matrix <- as.matrix(dist(df[, sp_features], method = "euclidean"))
    return(dist_matrix)
}

# Function to calculate attribute distances
#
# Description:
#   This function calculates the pairwise Euclidean distances between instances in a data 
#   frame based on specified attribute features.
#
# Inputs:
#   df           : A data frame containing the data.
#   env_features : A vector of strings representing the names of columns used for attribute distance calculation.
#
# Outputs:
#   A numeric matrix of pairwise Euclidean distances based on attribute features.
calculate_attribute_distance <- function(df, env_features) {
    dist_matrix <- as.matrix(dist(df[, env_features], method = "euclidean"))
    return(dist_matrix)
}

# Function to create pairwise distance matrix
#
# Description:
#   This function creates a weighted pairwise distance matrix by combining geospatial and 
#   attribute distances, with an optional normalization step.
#
# Inputs:
#   df          : A data frame containing the data.
#   sp_features : A vector of strings representing the names of columns used for geospatial distance calculation.
#   env_features: A vector of strings representing the names of columns used for attribute distance calculation.
#   rho         : A numeric value representing the weighting parameter for the geospatial distance.
#   normalize   : A boolean value indicating whether to normalize the distance matrices (default is TRUE).
#
# Outputs:
#   A numeric matrix representing the weighted pairwise distances.
create_pairwise_distance_matrix <- function(df, sp_features, env_features, rho, normalize = TRUE) {
    sp_dist <- calculate_geospatial_distance(df, sp_features)
    env_dist <- calculate_attribute_distance(df, env_features)
    
    if (normalize) {
        sp_dist <- normalize_matrix(sp_dist)
        env_dist <- normalize_matrix(env_dist)
    }
    
    pairwise_dist <- rho * sp_dist + (1 - rho) * env_dist
    return(pairwise_dist)
}

# Function to verify clustering
#
# Description:
#   This function verifies the validity of a clustering result against various constraints and criteria.
#
# Inputs:
#   clustering_df: A data frame containing the clustering results with columns for instance IDs and cluster assignments.
#   df           : A data frame containing the original data to be clustered.
#   id_col       : A string specifying the name of the column in 'df' that contains instance IDs.
#   sp_features  : A vector of strings specifying the names of spatial features in 'df'.
#   env_features : A vector of strings specifying the names of environmental features in 'df'.
#   constraints_df: A data frame containing instance-level constraints with columns for instance IDs and constraints.
#   rho          : A numeric value specifying the weight of spatial features relative to environmental features in distance calculation.
#   threshold    : A numeric value specifying the threshold distance for cluster-level constraints.
#   f            : A string specifying the type of cluster-level constraint, either "min_diameter" or "max_separation".
#
# Outputs:
#   A boolean value indicating whether the clustering result is valid (TRUE) or not (FALSE).
verify_clustering <- function(clustering_df, df, id_col, sp_features, env_features, constraints_df, rho, threshold, f) {
    
    id <- "id"
    
    # Ensure that every instance is assigned to a cluster and no instance is assigned to more than one cluster
    if (anyDuplicated(clustering_df[[id]]) > 0) {
        return(FALSE)
    }
    
    if (any(is.na(clustering_df$cluster))) {
        return(FALSE)
    }
    
    # Ensure that every cluster has at least one instance
    if (any(table(clustering_df$cluster) == 0)) {
        return(FALSE)
    }

    # Check instance-level constraints
    if (nrow(constraints_df) > 0) {
        for (i in seq_len(nrow(constraints_df))) {
            instance1 <- constraints_df[i, "id_col_1"]
            instance2 <- constraints_df[i, "id_col_2"]
            idx1 <- clustering_df$cluster[clustering_df[[id]] == instance1]
            idx2 <- clustering_df$cluster[clustering_df[[id]] == instance2]
            
            if (length(idx1) == 0 || length(idx2) == 0) {
                return(FALSE)
            }
            
            if (constraints_df[i, "constraint"] == "must-link" && idx1 != idx2) {
                return(FALSE)
            } else if (constraints_df[i, "constraint"] == "cannot-link" && idx1 == idx2) {
                return(FALSE)
            }
        }
    }

    # Calculate the pairwise distance matrix
    pairwise_dist <- create_pairwise_distance_matrix(df, sp_features, env_features, rho)

    # Check cluster-level constraints
    if (f == "min_diameter") {
        indices <- which(pairwise_dist > threshold, arr.ind = TRUE)
        for (idx in seq_len(nrow(indices))) {
            i <- indices[idx, 1]
            j <- indices[idx, 2]
            instance1 <- df[[id_col]][i]
            instance2 <- df[[id_col]][j]
            cluster1 <- clustering_df$cluster[clustering_df[[id]] == instance1]
            cluster2 <- clustering_df$cluster[clustering_df[[id]] == instance2]
            if (cluster1 == cluster2) {
                return(FALSE)
            }
        }
    } else if (f == "max_separation") {
        indices <- which(pairwise_dist < threshold, arr.ind = TRUE)
        for (idx in seq_len(nrow(indices))) {
            i <- indices[idx, 1]
            j <- indices[idx, 2]
            instance1 <- df[[id_col]][i]
            instance2 <- df[[id_col]][j]
            cluster1 <- clustering_df$cluster[clustering_df[[id]] == instance1]
            cluster2 <- clustering_df$cluster[clustering_df[[id]] == instance2]
            if (cluster1 != cluster2) {
                return(FALSE)
            }
        }
    }

    return(TRUE)
}

# Function to evaluate clustering
#
# Description:
#   This function evaluates the clustering result by comparing predicted clusters against ground truth clusters using 
#   various metrics including Adjusted Rand Index (ARI), Adjusted Mutual Information (AMI), and Normalized Information Distance (NID).
#
# Inputs:
#   pred_df : A data frame containing the predicted clustering results with an 'id' column and a 'cluster' column.
#   og_df   : A data frame containing the original data with ground truth cluster assignments.
#   id_col  : A string specifying the name of the column in 'og_df' that contains instance IDs.
#   gt_col  : A string specifying the name of the column in 'og_df' that contains ground truth cluster assignments.
#
# Outputs:
#   A list containing three elements:
#     - ari: Adjusted Rand Index
#     - ami: Adjusted Mutual Information
#     - nid: Normalized Information Distance
evaluate_clustering <- function(pred_df, og_df, id_col, gt_col) {
    # Ensure 'id' exists in both data frames
    if (!"id" %in% colnames(pred_df)) {
        stop("Column 'id_col' not found in pred_df")
    }
    if (!id_col %in% colnames(og_df)) {
        stop(paste("Column", id_col, "not found in og_df"))
    }
    
    # Ensure 'gt_col' exists in og_df
    if (!gt_col %in% colnames(og_df)) {
        stop(paste("Column", gt_col, "not found in og_df"))
    }
    
    # Merge pred_df with og_df using id_col
    merged_df <- merge(pred_df, og_df, by.x = "id", by.y = id_col)
    
    # Extract the true and predicted cluster assignments
    pred_sites <- as.factor(merged_df$cluster)
    og_sites <- as.factor(merged_df[[gt_col]])
    
    # Calculate ARI, AMI, and NID
    ari <- mclust::adjustedRandIndex(og_sites, pred_sites)
    ami <- aricode::AMI(og_sites, pred_sites)
    nid <- aricode::NID(og_sites, pred_sites)
    
    return(list(ari = ari, ami = ami, nid = nid))
}


