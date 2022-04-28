################
## YeastMain.R
##
## @authors Martyna Lukaszewicz, Erkan Buzbas, Paul Hohenlohe
## @contact martyna@uidaho.edu
##
## @license This is free software distributed without any warranty, without
##          even the implied warranty of merchantability or fitness for a particular
##          purpose. You are welcome to use, modify, or distribute it 
##          under the terms of the GNU General Public License,
##          version 3, as published by the Free Software Foundation.
##          See http://www.r-project.org/Licenses/GPL-3
##
################
## @description 
# Type: Main for analyzing SNP data 
# Model: a population genetics model of divergent selection 
#        with migration and recombination (described as in Author et al. 20xx)
# Analysis: using Approximate Bayesian Computation
##
## @param None
##
## @return None
##
## @lastChange 2021-04-12
##
## @changes
##  added scaling.R and LocLinRegCorABC.R paths
################
# begin all
#-------------------------------------------------------------
rm(list=ls())

# library(data.table)
# library(permute)
# library(matrixStats)
# library(MCMCpack)
#############


## PATHS
#############
# Edit full path to the base directory of your computer 
baseDir <- "~/" ## specify working directory
scriptDir <- paste0(baseDir, "/Dev/Ver0") # code directory
outputDir <- paste0(baseDir, "/Data") # data directory
#############
## FUNCTIONS
#############
source(paste0(scriptDir, "/Simulator.R"))
source(paste0(scriptDir, "/DoABCrej.R"))
source(paste0(scriptDir, "/LocLinRegCorABC.R"))
source(paste0(scriptDir, "/scaling.R"))
###################
## INPUT PARAMETERS
###################
## All parameters that are necessary to simulate data
## and run the analysis
nsim <- 1e5

###################
## OUTPUT 
###################
## Write to files in appropriate directories