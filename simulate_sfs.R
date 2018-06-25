# test SFS under different scenarios
library(readr)
library(stringr)

sim_name <- "nes_sim_test"
sample_size <- 94

# dirs
current_dir <- getwd()
simulate_folder <- "/home/martin/nes/nes_RAD/fastsimcoal_analyses/simulate_sfs/"
if(!dir.exists(simulate_folder)) dir.create(simulate_folder)

Ncur <- 10000
Nbot <- 10
Nhist <- 20000
Nlg <- 2000
T_bot_start <- 28
T_bot_end <- 20
T_lg <- 1000

Nbotprop <- Nbot/Ncur
Nhistprop <- Nhist/Nbot
Nlgprop <- Nlg/Nhist
growth_rate <- log(Nlg/Nhist) / (T_lg - T_bot_start)
  
  
# create par file
sink(paste0(simulate_folder, sim_name, ".par"))
cat("//Number of population samples (demes)");  cat("\n")
cat("1"); cat("\n")
cat("//Population effective sizes (number of genes)"); cat("\n")
cat(Ncur); cat("\n")
cat("//Sample sizes"); cat("\n")
cat("188"); cat("\n")
cat("//Growth rates : negative growth implies population expansion"); cat("\n")
cat("0"); cat("\n")
cat("//Number of migration matrices : 0 implies no migration between demes"); cat("\n")
cat("0"); cat("\n")
cat("//historical event: time, source, sink, migrants, new size, new growth rate, migr. matrix"); cat("\n")
cat("3 historical event"); cat("\n")
cat(paste0(T_bot_end, " 0 0 0 ", Nbotprop, " 0 0")); cat("\n")
cat(paste0(T_bot_start, " 0 0 0 ", Nhistprop, " ", growth_rate, " 0")); cat("\n")
cat(paste0(T_lg, " 0 0 0 ", Nlgprop));cat("\n")
cat("//Number of independent loci [chromosome]"); cat("\n")
cat("200000 0"); cat("\n")
cat("//Per chromosome: Number of linkage blocks"); cat("\n")
cat("1"); cat("\n")
cat("//per Block: data type, num loci, rec. rate and mut rate + optional parameters"); cat("\n")
cat("DNA 100 0  0.000000025 OUTEXP")
sink()

# number of simulations
num_sim <- 1

setwd(simulate_folder)
system(paste("/home/martin/bin/fsc26 -i ", sim_name, ".par -n ", num_sim, " -j -m -s0 -x -q", sep = ""))
setwd(current_dir)

# readLines(paste0(simulate_folder, "nes_sim_iceage/nes_sim_iceage_1/nes_sim_iceage_MAFpop0.obs"))[3] %>% 
nes_sfs <- readLines(paste0(simulate_folder, sim_name, "/",sim_name, "_1/", sim_name, "_MAFpop0.obs"))[3] %>% 
            str_split(pattern = "\t") %>% 
            unlist() %>% 
            as.numeric() %>% 
            .[1:c(sample_size + 1)]

sum(nes_sfs[-1]) / sum(nes_sfs)          
sfs_prop_sim <- nes_sfs/sum(nes_sfs)
barplot(sfs_prop_sim[-c(1)])
barplot(sfs_prop_sim[-c(1,2)])
barplot(nes_sfs[-c(1)])

# proportions
# polymorphic
sum(nes_sfs[-1]) / sum(nes_sfs) # 0.0022
# singletons
sum(nes_sfs[2]) / sum(nes_sfs) # 0.0016
# doubletons
sum(nes_sfs[3]) / sum(nes_sfs) #0.00005
# singletons among polymorphic
sum(nes_sfs[2]) / sum(nes_sfs[-1]) # 0.72
# doubletons among polymorphic
sum(nes_sfs[3]) / sum(nes_sfs[-1]) # 0.026


