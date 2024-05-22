#########################################################################################

# Geospatial Constrained Clustering as a Boolean Satisfiability (SAT) Problem

# May 22, 2024

# This file runs experiments on geospatial clustering using various constraints and parameters

#########################################################################################

# Source necessary R scripts
source("clust.R")
source("helpers.R")


# Function to run experiments over different combinations of constraints, number of constraints,
# selection methods, lambda values, and rho values
#
# Description:
#   This function runs experiments on geospatial clustering using various combinations of 
#   constraint types, number of constraints, lambda values, and rho values. It generates 
#   constraints, runs the geospatial constrained clustering algorithm, and records the results.
#   The function saves the constraints and clustering results to CSV files, verifies the clustering,
#   and evaluates the clustering using metrics such as ARI, AMI, and NID.
#
# Inputs:
#   df             : Data frame of instances.
#   id_col         : String, the name of the ID column.
#   gt_col         : String, the name of the ground truth cluster column.
#   sp_features    : Vector of strings, names of the spatial features in the data frame.
#   env_features   : Vector of strings, names of the environmental features in the data frame.
#   constraint_type: Vector of strings, types of constraints (e.g., "none", "must-link", "cannot-link", "mixed").
#   num_constraint : Vector of integers, number of constraints to be tested.
#   lambda_values  : Vector of integers, lambda values (percentage of instances to be used for clustering).
#   rho_values     : Vector of numeric, rho values for weighted distance calculation.
#   objectives     : Vector of strings, names of the objective functions to be tested (e.g., "min_diameter", "max_separation").
#
# Outputs:
#   None
#   The function saves multiple CSV files containing the generated constraints, clustering results,
#   call times, and a summary of the experiment results, including verification and evaluation metrics.
run_experiments <- function(df, id_col, gt_col, sp_features, env_features, constraint_type, num_constraint, lambda_values, rho_values, objectives) {
    
    # Initialize a dataframe to store the output values
    output_df <- data.frame(con_type = character(),
                            n_con = integer(),
                            k = integer(),
                            rho = numeric(),
                            f = character(),
                            gamma = numeric(),
                            solver_calls = numeric(),
                            sat_calls = numeric(),
                            unsat_calls = numeric(),
                            sat_avg_nclauses = numeric(),
                            unsat_avg_nclauses = numeric(),
                            sat_avg_time = numeric(),
                            unsat_avg_time = numeric(),
                            verification = character(),
                            ari = character(),
                            ami = character(),
                            nid = character(),
                            stringsAsFactors = FALSE)

    
    # Iterate over each constraint type
    for (con_type in constraint_type) {
        
        num_constraint_t <- num_constraint
        if (con_type == "none"){
            num_constraint_t <- c(0)
        }
        # Iterate over each number of constraints
        for (n_con in num_constraint_t) {
            
            # Generate constraints dataframe
            constraints_df <- get_constraints_dataframe(df, id_col, gt_col, con_type, n_con)
            
            # Save the generated constraints to a CSV file
            write.csv(constraints_df, paste0("results/constraints/constraints_type=", con_type, "_num=", n_con, ".csv"), row.names = FALSE)

            # Iterate over each lambda value
            for (lambda in lambda_values) {

                # Calculate the number of clusters based on lambda
                k <- ceiling((lambda / 100) * nrow(df))
                
                # Iterate over each rho value
                for (rho in rho_values) {
                    
                    # Iterate over each objective function
                    for (f in objectives) {
                        
                        # Print the current combination of indexing variables
                        info_str <- paste0("Processing: Constraint type = ", con_type, ", No. of constraints = ", n_con, ", k = ", k, ", rho = ", rho, ", f = ", f)
                        cat(strrep("-",nchar(info_str)),"\n")
                        cat(info_str, "\n")
                        cat(strrep("-",nchar(info_str)),"\n")

                        # Measure the time taken to run geo_con_clust_sat
                        result <- geo_con_clust_sat(df, k, sp_features, env_features, rho, constraints_df, f, id_col)
                  

                        # Print the first few rows of the result dataframe and the optimal value
                        # print(head(result$result, 5))
                        cat("clusters created: ", length(unique(result$result$cluster)), "\n")
                        cat("optimal value: ", result$optimal_value, "\n")
                        call_times_df <- result$call_times

                        # Save the result dataframe to a CSV file
                        result_path <- paste0("results/clusterings/clustering_type=", con_type, "_num=", n_con, "_k=", k, "_rho=", rho, "_f=", f, ".csv")
                        write.csv(result$result, result_path, row.names = FALSE)

                        # Save the result dataframe to a CSV file
                        result_path <- paste0("results/call_times/call_times_type=", con_type, "_num=", n_con, "_k=", k, "_rho=", rho, "_f=", f, ".csv")
                        write.csv(call_times_df, result_path, row.names = FALSE)

                        # Verify clustering
                        verification <- verify_clustering(result$result, df, id_col, sp_features, env_features, constraints_df, rho, result$optimal_value, f)
                        # cat("verification: ", verification, "\n")

                        # Evaluate clustering
                        evaluation <- evaluate_clustering(result$result, df, id_col, gt_col)
                        # cat("evaluation: ARI = ", evaluation$ari,", AMI = ", evaluation$ari,", NID = ", evaluation$nid, "\n")
                                                
                        
                        # print(call_times_df)

                        # Add the optimal value and elapsed time to the output_df
                        output_df <- rbind(output_df, data.frame(con_type = con_type, 
                                                                n_con = n_con, 
                                                                k = k, 
                                                                rho = rho, 
                                                                f = f, 
                                                                gamma = round(result$optimal_value, 6),
                                                                solver_calls = nrow(call_times_df),
                                                                sat_calls = nrow(call_times_df[call_times_df$sat_output=="SAT",]),
                                                                unsat_calls = nrow(call_times_df[call_times_df$sat_output=="UNSAT",]),
                                                                sat_avg_nclauses = round(mean(as.numeric(call_times_df[call_times_df$sat_output=="SAT",]$num_clauses)), 4),
                                                                unsat_avg_nclauses = round(mean(as.numeric(call_times_df[call_times_df$sat_output=="UNSAT",]$num_clauses)), 4),
                                                                sat_avg_time = round(mean(call_times_df[call_times_df$sat_output=="SAT",]$time_taken), 4),
                                                                unsat_avg_time = round(mean(call_times_df[call_times_df$sat_output=="UNSAT",]$time_taken), 4),
                                                                verification = verification,
                                                                ari = round(evaluation$ari, 6),
                                                                ami = round(evaluation$ami, 6),
                                                                nid = round(evaluation$nid, 6)))
                        cat(strrep("-",nchar(info_str)),"\n\n")
                    }
                }
            }
        }
    }

    # Save the output values dataframe to a CSV file
    write.csv(output_df, "results/output.csv", row.names = FALSE)
}

title_str <- "GeoConClustSAT: Geospatial Constrained Clustering as a Boolean Satisfiability (SAT) Problem"
title_border <- strrep("#", nchar(title_str))

cat("\n\n\n\n")
cat(title_border, "\n")
cat(title_str, "\n")
cat(title_border, "\n\n\n")
                        
# Define spatial and environmental features
sp_features <- c("longitude", "latitude")
env_features <- c("summer_nbr_TCB_mean_2400", "spring_nbr_TCG_mean_75", "fall_nbr_B2_stdDev_300", "slope_mean_300", "summer_b5_B4_stdDev_1200")

# Define columns for unique identifiers and ground truth clusters
id_col <- "checklist_id"
gt_col <- "site"

# Load the dataset
df <- read.csv("data.csv")

df_perc <- 15

# Sample the dataset
df <- sample_df(df, df_perc, id_col, gt_col)

write.csv(df, paste0("data_",df_perc,".csv"), row.names = FALSE)

# Convert spatial and environmental features columns to numeric if they are not already
df[sp_features] <- lapply(df[sp_features], as.numeric)
df[env_features] <- lapply(df[env_features], as.numeric)

# Calculate the number of unique clusters and lambda value
true_k <- get_true_k(df, gt_col)
true_lambda <- (true_k/nrow(df)) * 100

cat("Clusters (k): ", true_k, ", Instances (n): ", nrow(df), ", Lambda (k/n): ", true_lambda,"\n\n")

# Define the parameters for the experiments
constraint_type <- c("none", "must-link", "cannot-link", "mixed")
num_constraint <- c(1,2,3,4,5)
lambda_values <- c(true_lambda) # lambda = ((No. of clusters)/(No. of instances)) * 100 
rho_values <- c(0.25, 0.5, 0.75)
objectives <- c("min_diameter")



# Temporary override for testing
# constraint_type <- c("none")
# num_constraint <- c(3)
# lambda_values <- c(true_lambda)
# rho_values <- c(0.5)
# objectives <- c("max_separation")




# Run the experiments
run_experiments(
    df = df, 
    id_col = id_col, 
    gt_col = gt_col, 
    sp_features = sp_features, 
    env_features = env_features, 
    constraint_type = constraint_type, 
    num_constraint = num_constraint, 
    lambda_values = lambda_values, 
    rho_values = rho_values,
    objectives = objectives
)