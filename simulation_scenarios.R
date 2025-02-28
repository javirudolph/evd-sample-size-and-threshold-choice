# Get results for this on the allinone script

# Libraries -----------------------------------------------
set.seed(271367)

library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyr)
library(purrr)
library(fitdistrplus)
library(extRemes)


source("functions.R")

points_for_boxplots <- 30
B <- 100

# ORIGINAL SCENARIO --------------------------------------
dir_scenario <- "allinone_original"

if(fs::dir_exists(paste0("data/", dir_scenario)) == FALSE){
  fs::dir_create(paste0("data/", dir_scenario))}

desired_means <- c(28, 32, 40, 50)
desired_sds <- c(49.7, 39.9, 33.3, 31.1)
pars <- desired_mean_sd(mu_x = desired_means, sd_x = desired_sds)
source("simulation_process.R")

######################
# Supplementary scenarios

# SAME MEAN INCREASING SD ------------------------------
dir_scenario <- "ain1_samemean_upsd"

if(fs::dir_exists(paste0("data/", dir_scenario)) == FALSE){
   fs::dir_create(paste0("data/", dir_scenario))}

desired_means <- c(30, 30, 30, 30)
desired_sds <- c(20, 30, 40, 50)
pars <- desired_mean_sd(mu_x = desired_means, sd_x = desired_sds)
source("simulation_process.R")


#INCREASING MEANS DECRESING SD --------------------------------------
dir_scenario <- "ain1_upmean_downsd"

if(fs::dir_exists(paste0("data/", dir_scenario)) == FALSE){
   fs::dir_create(paste0("data/", dir_scenario))}

desired_means <- c(30, 35, 40, 50)
desired_sds <- c(50, 40, 35, 30)
pars <- desired_mean_sd(mu_x = desired_means, sd_x = desired_sds)
source("simulation_process.R")


# SAME SD INCREASING MEAN ------------------------------
dir_scenario <- "ain1_samesd_upmean"

if(fs::dir_exists(paste0("data/", dir_scenario)) == FALSE){
   fs::dir_create(paste0("data/", dir_scenario))}

desired_means <- c(30, 40, 50, 60)
desired_sds <- c(25, 25, 25, 25)
pars <- desired_mean_sd(mu_x = desired_means, sd_x = desired_sds)
source("simulation_process.R")

# INCREASING SD INCREASING MEAN ------------------------------
dir_scenario <- "ain1_upsd_upmean"

if(fs::dir_exists(paste0("data/", dir_scenario)) == FALSE){
   fs::dir_create(paste0("data/", dir_scenario))}

desired_means <- c(30, 40, 50, 60)
desired_sds <- c(20, 30, 40, 50)
pars <- desired_mean_sd(mu_x = desired_means, sd_x = desired_sds)
source("simulation_process.R")


# DECREASING BOTH ------------------------------
dir_scenario <- "ain1_downsd_downmean"

if(fs::dir_exists(paste0("data/", dir_scenario)) == FALSE){
   fs::dir_create(paste0("data/", dir_scenario))}

desired_means <- c(60, 50, 40, 30)
desired_sds <- c(50, 40, 30, 20)
pars <- desired_mean_sd(mu_x = desired_means, sd_x = desired_sds)
source("simulation_process.R")



