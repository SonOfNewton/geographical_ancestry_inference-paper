# geographical_ancestry_inference-paper project master script
# Copyright (c) 2026 Zhejiang University 
# This master file sources all other R scripts shell and script in the project.
# It is intended to be run from the root directory of the project and allow for replication of all analyses and figures.

# Before running this script, please set the working directory (setwd()) to the root directory of the project 
# Ensure all dependencies are installed
pkgs <- c("gaia","igraph","ggplot2","tidyverse","ggpubr","parallel","sf","terra","igraph","viridis","here") 

# install if not installed (nothing to do with updates)
new.packages <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages,repos = "http://cran.us.r-project.org")
lapply(pkgs, library, character.only = TRUE)

# SLiM application needs to be installed
# required input files: 
# in data/geo: 2020_walking_only_friction_surface.geotiff, landgrid_adjmat_afro-eurasia.csv landgrid_afro-eurasia.gpkg
# in data/pop: popc_5000BC.asc (and other files including 10000BC, 0AD, 1500AD, 2000AD)
# in data/genetics: hgdp_tgp_sgdp_high_cov_ancients_chr18_p.dated.trees

message("\n Master code started: all librairies loaded\n")

#define number of cores to use
maxcores <- parallel::detectCores()    
usecores <- maxcores - 1
if (usecores < 1) usecores <- 1

#generate folders if not existent
folders <- c(here('output'),here('output','figures'),here('output','tables'),
             here('data'),here('data','trees'),here('data','math'),here('data','geo'),here('data','pop'))
lapply(folders, dir.create, recursive = TRUE, showWarnings = FALSE)

# load all functions required
source(here("code","functions.R"),verbose=FALSE)
message("\n all functions loaded\n")


# ##################### for figure 1 ##############################################
# carry out the plot on China coverage and overlap ratio/intensity over time
source(here("code", "visualization", "tree_demo.R"),verbose=FALSE)
message("\n successfully generated abstract illustration plot\n\n")

system2(command = "bash",args = c(here("code", "generation", "slim_math.sh"), "line", "naive", usecores))
message("\n successfully generated tree data for figure 1\n\n")

source(here("code", "visualization", "visualize_line.R"),verbose=FALSE)
message("\n successfully generated actual and estimated tree sequence plot\n\n")

# tip 1: math.slim can also run on SLiM gui, simply by adding absolute path in function initialize():
# setwd("/Users/jianjunlian/Desktop/work/研二/HbS_v2/real_github/geographical_ancestry_inference-paper");
# and parameter will not be passed by bash but use default instead:
# if (!exists("TOPOLOGY"))
#   defineConstant("TOPOLOGY", "line");
# if (!exists("MODEL"))
#   defineConstant("MODEL", "naive");

# tip 2: .sh files are sensitive to space, so edit with care.
# (a redundant space at the end of the line can fail the whole pipeline)

# ##################### for figure 2 ##############################################
system2(command = "bash",args = c(here("code", "generation", "slim_math.sh"), "line", "friction", usecores))
system2(command = "bash",args = c(here("code", "generation", "slim_math.sh"), "square", "friction", usecores))
system2(command = "bash",args = c(here("code", "generation", "slim_math.sh"), "cube", "friction", usecores))
message("\n successfully generated tree sequence data for figure 2\n\n")

source(here("code", "visualization", "compare_friction_with_naive.R"),verbose=FALSE)
message("\n successfully generated actual and estimated tree sequence plot\n\n")

# ##################### for figure 3 ##############################################
my_world = "afro-eurasia"
source(here("code", "generation", "data_preparation.R"),verbose=FALSE)
message("\n data preparation complete\n\n")

system2(command = "bash",args = c(here("code", "generation", "slim_empirical.sh"), usecores))
message("\n successfully generated tree sequence data for figure 3\n\n")

source(here("code", "simulation", "select_worlds.R"),verbose=FALSE)
system2(command = "bash",args = c(here("code", "simulation", "run_gaia.sh"), my_world, "friction", usecores))
system2(command = "bash",args = c(here("code", "simulation", "run_gaia.sh"), my_world, "naive", usecores))
message("\n simulation complete\n\n")

source(here("code", "visualization", "compare_flux.R"),verbose=FALSE)
source(here("code", "visualization", "compare_ancestor_estimates.R"),verbose=FALSE)
message("\n visualization complete\n\n")

# ##################### for figure 4 ##############################################
my_world == "asia-americas"
source(here("code", "generation", "data_preparation.R"),verbose=FALSE)
message("\n data preparation complete\n\n")



# End
message("\n End of master.R. Thanks and have a good day! \n\n")


