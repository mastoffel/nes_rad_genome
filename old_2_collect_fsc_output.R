# collect best lhoods from different fastsimcoal runs

library(readr)
library(purrr)
library(stringr)
library(parallel)
options(scipen=999)
# how many bootstraps
nboot <- 100

# check folders in fsc_run/
folders <- list.files("fastsimcoal_analyses/fsc_run/", pattern = "fscrun[0-9]")

# collect likelihoods 
collect_vars <- function(folder){
  fsc_pars <- read_delim(paste0("fastsimcoal_analyses/fsc_run/", folder, "/nes/nes.bestlhoods"), delim = "\t")
}

all_vars <- folders %>% 
  map(safely(collect_vars)) %>% 
  purrr::transpose()
all_vars_working <- do.call(rbind, all_vars[[1]])

# determine highest likelihood
min_index <- which.min(all_vars_working$MaxEstLhood)
# folder with best lh 
folder_name <- folders[min_index]

### Parametric bootstrapping
bootstrap_folder <- "/home/martin/nes/nes_RAD/fastsimcoal_analyses/bootstrap_test/"

# (1) manipulate par_file
par_file <- readLines(paste0("fastsimcoal_analyses/fsc_run/", folder_name, "/nes/nes_maxL.par"))
# Number of sites in observed SFS: Number of lines in fasta file / 2 (roughly)
num_sites <- 2000000
# new par_file
loci_row <- which(str_detect(par_file, "Number of independent loci")) + 1
datatype_row <- which(str_detect(par_file, "FREQ"))
# replace loci
par_file[loci_row] <- str_c(num_sites, " 0")
# replace datatype
par_file[datatype_row] <- "DNA 1000 0  2.5e-8 OUTEXP"
# write to new folder
if (!dir.exists(bootstrap_folder)) {
  system(paste0("mkdir ", bootstrap_folder)) 
}
writeLines(par_file, paste0(bootstrap_folder, "nes.par"))

# (2) create 100 SFS in the bootstrap folder
setwd(bootstrap_folder)
system(paste("/home/martin/bin/fsc26 -i nes.par -n ", nboot, " -j -m -s0 -x -q", sep = ""))
setwd(paste0("/home/martin/nes/nes_RAD"))


#### Second part

setwd(paste0("/home/martin/nes/nes_RAD"))
# (3) Run Fsc 50 times in each folder
# create new folder for likelihood maximization of all 100 bootstrap SFS samples
if (!dir.exists("/home/martin/nes/nes_RAD/fastsimcoal_analyses/analyse_bootstrap_sfs")) {
  system("mkdir /home/martin/nes/nes_RAD/fastsimcoal_analyses/analyse_bootstrap_sfs")
}
setwd("/home/martin/nes/nes_RAD/fastsimcoal_analyses/analyse_bootstrap_sfs")

run_fsc <- function(run_num, num_boot){
  # create directory for fsc run
  system(paste0("mkdir fscrun", run_num))
  # paste all relevant files into directory
  system(paste0("cp ../../fastsimcoal_files/nes* fscrun", run_num))
  system(paste0("cp ../../bootstrap/nes/nes_", num_boot, "/nes_MAFpop0.obs fscrun", run_num))
  # change to directory
  setwd(paste0("/home/martin/nes/nes_RAD/fastsimcoal_analyses/analyse_bootstrap_sfs/boot_", num_boot, "/fscrun", run_num))
  # run fsc
  system("~/bin/fsc26 -t nes.tpl -n 100000 -m -e nes.est -M -L 40 -q -w 0.01 -x --foldedSFS -C 5 --nosingleton")
  # change back
  setwd(paste0("/home/martin/nes/nes_RAD/fastsimcoal_analyses/analyse_bootstrap_sfs/boot_", num_boot))
}

# run all
for (num_boot in 1:nboot) {
  
  if (!dir.exists(paste0("/home/martin/nes/nes_RAD/fastsimcoal_analyses/analyse_bootstrap_sfs/boot_", num_boot))) {
    system(paste0("mkdir /home/martin/nes/nes_RAD/fastsimcoal_analyses/analyse_bootstrap_sfs/boot_", num_boot))
  }
  setwd(paste0("/home/martin/nes/nes_RAD/fastsimcoal_analyses/analyse_bootstrap_sfs/boot_", num_boot))
  
  cl <- makeCluster(getOption("cl.cores", 10))
  parLapply(cl, 1:10, run_fsc, num_boot)
  stopCluster(cl)
  
  setwd("/home/martin/nes/nes_RAD/fastsimcoal_analyses/analyse_bootstrap_sfs")
}

setwd("/home/martin/nes/nes_RAD/")

