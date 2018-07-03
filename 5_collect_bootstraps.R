# collect vars from bootstrap replicates

library(readr)
library(purrr)
library(dplyr)
library(stringr)
library(parallel)
library(tibble)
options(scipen=999)
# how many bootstraps
nboot <- 100
num_sims <- 10
setwd("/home/martin/nes/nes_rad")

sim_name <- "lgp_bot"
# folders to collect from
# bootstrap folders
folders_boot <- list.files(paste0("fastsimcoal_analyses/nonpar_bootstrap_", sim_name, "/"), pattern = "boot")
# simulations per bootstrap replicate folders
# folders_sim <- list.files("fastsimcoal_analyses/analyse_bootstrap_sfs/boot_1/", pattern = "fscrun[0-9]")
folders_sim <- paste0("fscrun", 1:num_sims)
# combine
folders_df <- data.frame("boot" = rep(folders_boot, each = length(folders_sim)), "sims" = rep(folders_sim, times = length(folders_boot)))

# collect likelihoods 
collect_vars <- function(folder_boot, folder_sim){
  fsc_pars <- read_delim(paste0("fastsimcoal_analyses/nonpar_bootstrap_", sim_name, "/", folder_boot, "/", folder_sim, "/nes/nes.bestlhoods"), delim = "\t")
}
# collect everything (usually just take maximum likelihood per bootstrap)
all_vars <- map2(folders_df$boot, folders_df$sims, possibly(collect_vars, NA_real_)) 
# put all in df
all_vars <- all_vars[!is.na(all_vars)]
all_vars_working <- bind_rows(all_vars)


library(readr)
write_delim(all_vars_working, path = "data/fsc_out/nonpar_boot_lpg_talk.txt", delim = " ")
hist(all_vars_working$NBOT, breaks = 100)

par_ests <- all_vars_working %>% 
  mutate(boots = as.factor(MaxObsLhood)) %>% 
  group_by(boots) %>% 
  arrange(desc(MaxEstLhood), .by_group = TRUE) %>% 
  top_n(n = 2, MaxEstLhood)


hist(par_ests$NANC, breaks = 10)

skimr::skim(all_vars_working)
quantile(all_vars_working$NCUR)




df <- data.frame(x = c(10, 4, 1, 6, 3, 1, 1))
df %>% top_n(2)