# collect best lhoods from different fastsimcoal runs

library(readr)
library(purrr)
library(stringr)
library(parallel)
options(scipen=999)
setwd(paste0("/home/martin/nes/nes_rad"))
# how many bootstraps
# nboot <- 100
# num_sim <- 50


# check folders in fsc_run/
sim_name <- "lgp_bot"
path_to_runs <- paste0("fastsimcoal_analyses/fsc_run_", sim_name, "/")
folders <- list.files(path_to_runs , pattern = "fscrun[0-9]")

# collect likelihoods 
collect_vars <- function(folder){
  fsc_pars <- read_delim(paste0(path_to_runs, folder, "/nes/nes.bestlhoods"), delim = "\t")
}

all_vars <- folders %>% 
  map(safely(collect_vars)) %>% 
  purrr::transpose()

all_vars_working <- do.call(rbind, all_vars[[1]]) %>% 
  dplyr::mutate(lh_diff = MaxEstLhood - MaxObsLhood)

all_vars_working[which.max(all_vars_working$MaxEstLhood), ]

hist(all_vars_working$NBOT)
