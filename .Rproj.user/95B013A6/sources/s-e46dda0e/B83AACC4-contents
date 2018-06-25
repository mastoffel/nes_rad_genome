# Create all files to run fastsimcoal
library(stringr)

# parameters

Ncur <- "unif 2000 10000"
Nhist <- "unif 5000 50000"
Nbot <- "unif 1 50"
Tbotend <- "20"
Tbotstart <- "30"

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


##### BOTTLENECK MODEL ##########

# define folder for fastsimcoal files
name_run <- "bot"
folder_run <- paste0("/home/martin/nes/nes_rad/fastsimcoal_analyses/fastsimcoal_files/", name_run)
if (!dir.exists(folder_run)) dir.create(folder_run)

# create obs file from ANGSD SFS

# imports an sfs from angsd to the fastsimcoal_files folder
sfs <- scan("../angsd_analysis/SFS12/nes.sfs")
# add zeros for compatability with fsc
sfs <- c(sfs, rep(0, length(sfs) - 1))
length(sfs)
# create names
sfs_names <- sapply(0:(length(sfs) - 1), function(x) paste0("d0_", x))

sink(paste0(folder_run, "/nes_MAFpop0.obs"))
cat("1 observations")
cat("\n")
cat(sfs_names)
cat("\n")
cat(sfs)
sink()

# create tpl file
sink(file = paste0(folder_run, "/nes.tpl"))
cat(c("//Number of population samples (demes)",
      "1",
      "//Population effective sizes (number of genes)", 
      "NCUR", 
      "//Sample sizes", 
      "188",
      "//Growth rates : negative growth implies population expansion", 
      "0", 
      "//Number of migration matrices : 0 implies no migration between demes",  
      "0", 
      "//historical event: time, source, sink, migrants, new size, new growth rate, migr. matrix", 
      "2 historical event",  
      "20 0 0 0 RESBOT 0 0",    
      "30 0 0 0 RESENDBOT 0 0",
      "//Number of independent loci [chromosome]",
      "1 0",
      "//Per chromosome: Number of linkage blocks", 
      "1", 
      "//per Block: data type, num loci, rec. rate and mut rate + optional parameters",
      "FREQ 1 0 2.5e-8"),
    sep = "\n"
)
sink()


# est file

sink(file = paste0(folder_run, "/nes.est"))
cat(c("// Priors and rules file",
      "// *********************",
      "[PARAMETERS]", 
      "1 NCUR unif 5000 20000 output", 
      "1 NANC unif 5000 40000 output", 
      "1 NBOT unif 1 50 output",
      "[RULES]", 
      "NCUR < NANC",
      "[COMPLEX PARAMETERS]",  
      "0 RESBOT = NBOT/NCUR hide", 
      "0 RESENDBOT = NANC/NBOT hide"),
    sep = "\n"
)
sink()





##### LGP + BOTTLENECK MODEL ##########

name_run <- "lgp_bot"
folder_run <- paste0("/home/martin/nes/nes_rad/fastsimcoal_analyses/fastsimcoal_files/", name_run)
if (!dir.exists(folder_run)) dir.create(folder_run)

# create obs file from ANGSD SFS

# imports an sfs from angsd to the fastsimcoal_files folder
sfs <- scan("../angsd_analysis/SFS12/nes.sfs")
# add zeros for compatability with fsc
sfs <- c(sfs, rep(0, length(sfs) - 1))
length(sfs)
# create names
sfs_names <- sapply(0:(length(sfs) - 1), function(x) paste0("d0_", x))

sink(paste0(folder_run, "/nes_MAFpop0.obs"))
cat("1 observations")
cat("\n")
cat(sfs_names)
cat("\n")
cat(sfs)
sink()


# create tpl file
sink(file = paste0(folder_run, "/nes.tpl"))
cat(c("//Number of population samples (demes)",
      "1",
      "//Population effective sizes (number of genes)", 
      "NCUR", 
      "//Sample sizes", 
      "188",
      "//Growth rates : negative growth implies population expansion", 
      "0", 
      "//Number of migration matrices : 0 implies no migration between demes",  
      "0", 
      "//historical event: time, source, sink, migrants, new size, new growth rate, migr. matrix", 
      "3 historical event",  
      "20 0 0 0 RESBOT 0 0",    
      "30 0 0 0 RESENDBOT 0 0",
      "1500 0 0 0 RESENDLGP GROWTH 0",
      "//Number of independent loci [chromosome]",
      "1 0",
      "//Per chromosome: Number of linkage blocks", 
      "1", 
      "//per Block: data type, num loci, rec. rate and mut rate + optional parameters",
      "FREQ 1 0 2.5e-8"),
    sep = "\n"
)
sink()


# est file

sink(file = paste0(folder_run, "/nes.est"))
cat(c("// Priors and rules file",
      "// *********************",
      "[PARAMETERS]", 
      "1 NCUR unif 5000 20000 output", 
      "1 NANC unif 5000 40000 output", 
      "1 NBOT unif 1 50 output",
      "1 NLGP unif 500 4000 output",
      "[RULES]", 
      "[COMPLEX PARAMETERS]",  
      "0 RESBOT = NBOT/NCUR hide", 
      "0 RESENDBOT = NANC/NBOT hide",
      "0 RESENDLGP = NLPG/NBOT hide",
      "0 SIZERATIO  = NLPG/NANC hide",
      "0 TEMP = log(SIZERATIO) hide",
      "0 GROWTH = TEMP/1470"),
    sep = "\n"
)
sink()

