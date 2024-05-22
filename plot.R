#########################################################################################

# Geospatial Constrained Clustering as a Boolean Satisfiability (SAT) Problem

# May 22, 2024

# This file contains functions for:
# - Generating and plotting the results of geospatial constrained clustering
# - Creating individual and combined plots to visualize the number of clauses, 
#   time taken by the SAT solver, and various evaluation metrics
# - Saving the generated plots to specified file paths

#########################################################################################


library(ggplot2, quietly = T)
library(gridExtra, quietly = T)
library(grid, quietly = T)

# Function to plot results from a CSV file
#
# Description:
#   This function reads a CSV file containing results of geospatial constrained clustering,
#   filters the data, creates individual and combined plots for various metrics, and saves
#   the plots to specified file paths.
#
# Inputs:
#   df_path : String, the file path to the CSV file containing the results.
#
# Outputs:
#   Saves the generated plots as PNG files in the specified directory.
plot_results <- function(df_path) {
  
  # Read the CSV file
  df <- read.csv(df_path)
  
  # Filter the dataframe to keep rows where rho == 0.5
  df <- df[df$rho == 0.5, ]
  
  # Helper function to create individual plots
  create_plot <- function(data, y, y_label, y_transformation = NULL) {
    if (!is.null(y_transformation)) {
      data[[y]] <- y_transformation(data[[y]])
    }
    
    # Create the plot with custom aesthetic mappings
    p <- ggplot(data, aes(x = n_con, y = get(y), color = con_type, shape = con_type)) +
      geom_line(aes(linetype = con_type, group = con_type), linewidth = 1) +
      geom_point(data = data[data$con_type != "none", ], size = 3) +
      geom_hline(aes(yintercept = mean(get(y)[data$con_type == "none"])), color = "black", linetype = "dashed") +
      scale_color_manual(values = c("none" = "black", "must-link" = "blue", "cannot-link" = "red", "mixed" = "purple"),
                         breaks = c("none", "must-link", "cannot-link", "mixed")) +
      scale_shape_manual(values = c("none" = NA, "must-link" = 19, "cannot-link" = 4, "mixed" = 17),
                         breaks = c("none", "must-link", "cannot-link", "mixed")) +
      scale_linetype_manual(values = c("none" = "dotdash", "must-link" = "solid", "cannot-link" = "solid", "mixed" = "solid"),
                            breaks = c("none", "must-link", "cannot-link", "mixed")) +
      scale_x_continuous(breaks = 1:5, limits = c(1, 5)) +
      labs(title = y_label, x = NULL, y = NULL) +
      theme_linedraw() +
      theme(legend.title = element_text(size = 14, hjust = 0.5, margin = margin(b = 10)), 
            legend.text = element_text(size = 12),
            legend.key.width = unit(1.5, "lines"),
            plot.title = element_text(size = 12, hjust = 0.5), # Increased title size
            axis.text = element_text(size = 10),
            legend.position = "bottom",
            plot.margin = unit(c(1, 1, 1, 1), "lines"),
            plot.tag.position = "bottomright") +
      guides(color = guide_legend(title = "Constraint Type", nrow = 1), shape = guide_legend(title = "Constraint Type"), linetype = guide_legend(title = "Constraint Type"))
    
    return(p)
  }
  
  # Suppress plotting to avoid windows popping open
  pdf(NULL)
  
  # Subplot creation for "No. of Clauses and Time Taken by SAT Solver"
  p1 <- create_plot(df, "unsat_avg_nclauses", "Avg. No. of Clauses when Unsatisfiable")
  p2 <- create_plot(df, "sat_avg_nclauses", "Avg. No. of Clauses when Satisfiable")
  p3 <- create_plot(df, "unsat_avg_time", "Avg. Time Taken when Unsatisfiable\n(seconds)")
  p4 <- create_plot(df, "sat_avg_time", "Avg. Time Taken when Satisfiable\n(minutes)", function(x) x / 60)
  
  # Combine plots for the first set with common x-axis label and single legend
  legend1 <- get_legend(p1)
  p_combined1 <- arrangeGrob(p1 + theme(legend.position = "none"), 
                             p2 + theme(legend.position = "none"),
                             p3 + theme(legend.position = "none"), 
                             p4 + theme(legend.position = "none"), 
                             nrow = 2, ncol = 2,
                             top = textGrob("No. of Clauses and Time Taken by SAT Solver", gp = gpar(fontsize = 14)), # Increased top title size
                             bottom = textGrob("Number of Constraints", gp = gpar(fontsize = 12)))
  p_combined1 <- grid.arrange(p_combined1, legend1, ncol = 1, heights = c(10, 1))
  
  # Save the first plot without displaying it
  ggsave("results/plots/solver_times.png", plot = p_combined1, width = 10, height = 8, dpi = 300, units = "in", device = "png")
  
  # Subplot creation for "Evaluation Metrics"
  p5 <- create_plot(df, "gamma", "Optimal Value (gamma)")
  p6 <- create_plot(df, "ari", "Adjusted Rand Index")
  p7 <- create_plot(df, "ami", "Adjusted Mutual Information")
  p8 <- create_plot(df, "nid", "Normalized Information Distance")
  
  # Combine plots for the second set with common x-axis label and single legend
  legend2 <- get_legend(p5)
  p_combined2 <- arrangeGrob(p5 + theme(legend.position = "none"), 
                             p6 + theme(legend.position = "none"),
                             p7 + theme(legend.position = "none"), 
                             p8 + theme(legend.position = "none"), 
                             nrow = 2, ncol = 2,
                             top = textGrob("Evaluation Metrics", gp = gpar(fontsize = 14)), # Increased top title size
                             bottom = textGrob("Number of Constraints", gp = gpar(fontsize = 12)))
  p_combined2 <- grid.arrange(p_combined2, legend2, ncol = 1, heights = c(10, 1))
  
  # Save the second plot without displaying it
  ggsave("results/plots/evaluation_metrics.png", plot = p_combined2, width = 10, height = 8, dpi = 300, units = "in", device = "png")
  
  # Close the suppressed plotting device
  dev.off()
}

# Helper function to extract the legend
#
# Description:
#   This function extracts the legend from a ggplot object.
#
# Inputs:
#   myggplot : A ggplot object from which to extract the legend.
#
# Outputs:
#   A grob object representing the legend.
get_legend <- function(myggplot) {
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# Example usage of plot_results function
plot_results("results/output.csv")
