#-------------------------------------------------------------------------------
#          MODELLING GENERALISATION GRADIENTS AS AUGMENTED GAUSSIANS
#-------------------------------------------------------------------------------

seed_num <- 9999
set.seed(seed_num)
file_name_root <- "output/"
dir.create(file_name_root)

# load packages 
library(tidyverse)
library(patchwork)
library(rstan)
library(bayestestR)
library(loo)

# user-defined parameters ------------------------------------------------------

# experiment-specific parameters
dim_vals <- seq(-.5, +.5, .05) # dimension values, CS+ should be at 0

# mcmc parameters
n_chains <- 4
n_iter <- 20000
n_burnin <- 1000
n_thin <- 1
n_samp <- 50
hdi_limit <- .95
rope_low <- c(-.05, .1, .1, 70) # ROPE limits for raw parameters (M, W-, W+, H)
rope_high <- c(+.05, .2, .2, 80)
# rope_low_diffs <- c(-.05, -.05, -.05, -2.5) # ROPE limits for group diffs (M, W-, W+, H)
# rope_high_diffs <- abs(rope_low_diffs)
augG_params <- c("M", "SDPlus", "SDMinus", "height", "noise", "predR", "log_lik",
                 "M_group", "SDPlus_group", "SDMinus_group", "height_group")
augG_group_params <- c("M_group", "SDPlus_group", "SDMinus_group", "height_group")

# figure parameters
graph_file_type <- ".jpeg"
fig_cols <- c("black", "red")
scat_shape <- 16
scat_size <- 1
scat_col <- alpha(fig_cols[2], .1)
gg_height <- 10
gg_width <- 10
dpi <- 600
n_facet_rows <- 5
x_breaks <- 1:11 # number of breaks on stimulus dimension
x_labs <- c(min(dim_vals), rep("", 4), "CS+", rep("", 4), max(dim_vals)) # labels on stimulus dimension (length should match x_breaks)
group_names <- c("easy", "medium", "hard")
fig_shapes <- c(21, 24, 22) # length should match group_names
density_cols <- c("#301A4B", "#6DB1BF", "#D4775E") # length should match group_names

# ------------------------------------------------------------------------------

# build augmented Gaussian model string
source("R/build_model.R")

# source functions
source("R/functions.R")

# ------------------------------------------------------------------------------
# example
file_name_root <- paste0("output/", "NSW57", "-")

nsw57out <- Fit_Aug_Gaussian(fileName = "data/NSW57-Data-Sim.csv",
                             modelFile = "models/AugGaus.stan",
                             modelName = "similarity",
                             groupNames = group_names,
                             dimVals = dim_vals,
                             params = augG_params, groupParams = augG_group_params,
                             ropeLow = rope_low, ropeHigh = rope_high, nSamp = n_samp,
                             ropeLowDiffs = rope_low_diffs, ropeHighDiffs = rope_high_diffs,
                             hdiLim = hdi_limit, nRow = n_facet_rows, figMult = 2)
