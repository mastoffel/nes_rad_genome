# run several iterations of fastsimcoal
library(parallel)

# run name
fsc_name <- "lgp_bot" #"bot"

folder_name_fsc <- paste0("/home/martin/nes/nes_rad/fastsimcoal_analyses/fsc_run_", fsc_name)
if (!dir.exists(folder_name_fsc)) {
  dir.create(folder_name_fsc)
}
setwd(folder_name_fsc)

# if (length(list.files("./")) != 0) { 
#     system("rm -r *")
#   }
# system("~/bin/fsc26")
#folders <- list.files("fsc_run/", pattern = "fsc_run[0-9]")

run_fsc <- function(run_num, fsc_name, folder_name_fsc){
  # create directory for fsc run
  if (!dir.exists(paste0("fscrun", run_num))) system(paste0("mkdir fscrun", run_num))
  # paste all relevant files into directory
  # obs file is same
  system(paste0("cp ../fastsimcoal_files/", fsc_name, "/nes_MAFpop0.obs fscrun", run_num))
  system(paste0("cp ../fastsimcoal_files/", fsc_name, "/nes.est fscrun", run_num))
  system(paste0("cp ../fastsimcoal_files/", fsc_name, "/nes.tpl fscrun", run_num))
  
  # change to directory
  setwd(paste0(folder_name_fsc,"/fscrun", run_num))
  # run fsc
  system("~/bin/fsc26 -t nes.tpl -n 100000 -m -e nes.est -M -L 40 -q -w 0.01 --nosingleton -0 --nosingleton  --foldedSFS -x -C 10") # --nosingleton-0
   
  # change back
  setwd(folder_name_fsc)
}

# run all
cl <- makeCluster(getOption("cl.cores", 10))
parLapply(cl, 1:10, run_fsc, fsc_name, folder_name_fsc)
stopCluster(cl)

setwd("/home/martin/nes/nes_rad/")
