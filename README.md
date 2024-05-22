
# GeoConClustSAT

## Geospatial Constrained Clustering as a Boolean Satisfiability (SAT) Problem


### Overview

GeoConClustSAT leverages the capabilities of SAT solvers to handle various clustering constraints, such as must-link and cannot-link constraints, ensuring that the clustering results adhere to the specified constraints.

Clone the project repository from GitHub:
```bash
git clone https://github.com/nahian-ahmed/GeoConClustSAT.git
```

### Project Structure

The project consists of the following main components:

1. **main.R**: The main script that sets up the environment, sources necessary scripts, and runs the experiments.
2. **experiments.R**: Manages and executes experiments on geospatial clustering using various constraints and parameters.
3. **helpers.R**: Contains helper functions for generating constraints, sampling data, creating distance matrices, and verifying clustering results.
4. **clust.R**: Contains functions for performing geospatial constrained clustering using SAT solvers, including generating constraints in DIMACS CNF format, calling the Z3 solver, and parsing the solver output.
5. **plot.R**: Provides functions for generating and plotting the results of the clustering experiments.

### Getting Started

#### Prerequisites

- R (version 3.6 or higher)
- Required R packages: `ggplot2`, `grid`, `gridExtra`, `dplyr`, `mclust`, `aricode`

#### Installation

Navigate to the project directory:

```bash
cd GeoConClustSAT
```

Install the required R packages:

```r
install.packages(c("ggplot2", "grid", "gridExtra", "dplyr", "mclust", "aricode"))
```

### Usage

#### Running Experiments

To run the experiments, execute the `main.R` script:

```r
source("main.R")
```

This will set up the working directory, source the necessary scripts, and run the experiments with the specified parameters.

#### Generating Plots

To generate and visualize the results, execute the `plot.R` script:

```r
source("plot.R")
```

This will create individual and combined plots to visualize the number of clauses, time taken by the SAT solver, and various evaluation metrics.

### Using Z3 Solver

This project uses command-line calls to the Z3 solver. Ensure that the Z3 solver is installed and accessible from the command line. You can download the Z3 solver from [https://github.com/Z3Prover/z3](https://github.com/Z3Prover/z3).

### Functions

#### experiments.R

- **run_experiments**: Runs experiments over different combinations of constraints and parameters.

#### helpers.R

- **generate_constraints**: Generates must-link and cannot-link constraints for geospatial constrained clustering.
- **sample_data**: Samples data while maintaining cluster integrity.
- **create_distance_matrix**: Creates pairwise distance matrices combining geospatial and attribute distances.
- **verify_clustering**: Verifies the validity of clustering results.
- **evaluate_clustering**: Evaluates the clustering results using metrics such as Adjusted Rand Index (ARI), Adjusted Mutual Information (AMI), and Normalized Information Distance (NID).

#### clust.R

- **geo_con_clust_sat**: Performs geospatial constrained clustering using SAT solvers by generating constraints, calling the Z3 solver, and parsing the solver output.

#### plot.R

- **plot_results**: Generates and plots the results of geospatial constrained clustering.

### Contact

For any questions or issues, please contact:

- Nahian Ahmed: ahmedna@oregonstate.edu
