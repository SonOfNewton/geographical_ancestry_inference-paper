# geographical_ancestry_inference-paper project master script
# Copyright (c) 2026 Zhejiang University 
# This master file sources all other R scripts shell and script in the project.
# It is intended to be run from the root directory of the project and allow for replication of all analyses and figures.

setwd("/Users/jianjunlian/Desktop/work/研二/HbS_v2/real_github/geographical_ancestry_inference-paper")
# Before running this script, please set the working directory (setwd()) to the root directory of the project 
# Ensure all dependencies are installed
pkgs <- c("gaia","igraph","ggplot2","tidyverse","ggpubr") 

# remotes::install_github("hrbrmstr/waffle") #waffle not on CRAN anymore 
# install if not installed (nothing to do with updates)
new.packages <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages,repos = "http://cran.us.r-project.org")
lapply(pkgs, library, character.only = TRUE)

message("\n Master code started: all librairies loaded\n")

#necessary for polygon issues in topology
sf_use_s2(FALSE)

# start year of the spatial analysis
start_year <- 2000
# end year of the spatial analysis
end_year <- 2024

#define max cores that R can run
maxRcores <- parallel::detectCores()    
freecores <- 1

#generate folders if not existent
folders <- c(here('output'),here('output','raster'),here('output','tables'),
here('output','figsi'),here('output','rdata'),here('output','figsmanuscript'))
lapply(folders, dir.create, recursive = TRUE, showWarnings = FALSE)

# load all functions required
source(here("code","functions.R"),verbose=FALSE)
message("\n all functions loaded\n")

# do a minimum data cleaning used for all processes
source(here("code","dataprep.R"),verbose=FALSE)
message("\n data preparatation completed\n\n")

# carry out the spatial analysis 
# configuration for traveltime.R (computationally heavy)
#quick and dirty
quick <- TRUE
# check results during the pipeline with saving maps etc
check <- FALSE 
# write raster maps for China (SSD/2013) for illustration
wraster <- TRUE
# put number of hours for the max deployments as threshold
samplesize <- 10 #number of samples to compute spatial coverage and overlap metrics. Default 500
#minimum of deployments to compute overlap (default should be 2)
mindeploy <- 2
#generate all spatial metrics
source(here("code","traveltime.R"),verbose=FALSE)
message("\n successfully generated spatial coverage and travel time metrics\n\n")

# #####################for figure 2 ##############################################
# carry out the plot on China coverage and overlap ratio/intensity over time
source(here("code","coverageplot.R"),verbose=FALSE)
message("\n successfully generated coverage and overlap plots\n\n")

# generate a map of pko in countries with Chinese deployment (fig.2 panel)
source(here("code","mappko.R"),verbose=FALSE)
message("\n successfully generated deployment and spatial coverage map by country\n\n")

# # generate case study in SSD with Chinese deployments, access to strategic locations (fig.2 panel)
#select sample used for case study (it should be a sample that has about 6h restraint)
casestudysample <- 53 # sample index (choose map with about 6h travel time restraint89)
source(here("code","casestudyplot.R"),verbose=FALSE)
message("\n generated a map for figure 2 showing deployments by country \n\n")

# #####################for figure 3 ##############################################

# generate plots related to time to nearest road and nearest city (sup.fig)
source(here("code","time2objectsplot.R"),verbose=FALSE)
message("\n successfully generated time to road comparison plots\n\n")

# #####################for figure 1 ##############################################
# waffle plots troops and fatalities
source(here("code","waffletroop_plot.R"),verbose=FALSE)
message("\n successfully waffle plot on troop size and fatalities by TCC\n\n")

# waffle plots civilians and military
source(here("code","waffleciv_nonciv_plot.R"),verbose=FALSE)
message("\n successfully waffle plots of nb. of Chinese military and civilian by mission \n\n")

# End
message("\n End of master.R. Thanks and have a good day! \n\n")


