# use simulated parametric bootstraps SFS to generate CI
library(parallel)
nboot <- 100

setwd(paste0("/home/martin/nes/nes_rad"))
# (3) Run Fsc 50 times in each folder
# create new folder for likelihood maximization of all 100 bootstrap SFS samples
if (!dir.exists("/home/martin/nes/nes_rad/fastsimcoal_analyses/analyse_bootstrap_sfs")) {
  system("mkdir /home/martin/nes/nes_rad/fastsimcoal_analyses/analyse_bootstrap_sfs")
}
setwd("/home/martin/nes/nes_rad/fastsimcoal_analyses/analyse_bootstrap_sfs")

run_fsc <- function(run_num, num_boot){
  # create directory for fsc run
  system(paste0("mkdir fscrun", run_num))
  # paste all relevant files into directory
  system(paste0("cp ../../fastsimcoal_files/nes* fscrun", run_num))
  system(paste0("cp ../../bootstrap/nes/nes_", num_boot, "/nes_MAFpop0.obs fscrun", run_num))
  # change to directory
  setwd(paste0("/home/martin/nes/nes_rad/fastsimcoal_analyses/analyse_bootstrap_sfs/boot_", num_boot, "/fscrun", run_num))
  # run fsc
  system("~/bin/fsc26 -t nes.tpl -n 100000 -m -e nes.est -M -L 40 -q -w 0.01 -x --foldedSFS -C 5 --nosingleton")
  # change back
  setwd(paste0("/home/martin/nes/nes_rad/fastsimcoal_analyses/analyse_bootstrap_sfs/boot_", num_boot))
}

# run all
for (num_boot in 1:nboot) {
  
  if (!dir.exists(paste0("/home/martin/nes/nes_rad/fastsimcoal_analyses/analyse_bootstrap_sfs/boot_", num_boot))) {
    system(paste0("mkdir /home/martin/nes/nes_rad/fastsimcoal_analyses/analyse_bootstrap_sfs/boot_", num_boot))
  }
  setwd(paste0("/home/martin/nes/nes_rad/fastsimcoal_analyses/analyse_bootstrap_sfs/boot_", num_boot))
  
  cl <- makeCluster(getOption("cl.cores", 12))
  parLapply(cl, 1:10, run_fsc, num_boot)
  stopCluster(cl)
  
  setwd("/home/martin/nes/nes_rad/fastsimcoal_analyses/analyse_bootstrap_sfs")
}

setwd("/home/martin/nes/nes_rad/")