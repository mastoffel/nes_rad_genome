# imports an sfs from angsd to the fastsimcoal_files folder
sfs <- scan("../angsd_analysis/nes37.sfs")
# add zeros for compatability with fsc
sfs <- c(sfs, rep(0, length(sfs)-1))
length(sfs)
# create names
sfs_names <- sapply(0:(length(sfs)-1), function(x) paste0("d0_", x))

sink("fastsimcoal_analyses/fastsimcoal_files/nes_MAFpop0.obs")
cat("1 observations")
cat("\n")
cat(sfs_names)
cat("\n")
cat(sfs)
sink()