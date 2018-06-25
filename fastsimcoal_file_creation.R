# test run with fastsimcoal

name_run <- "bot"
folder_run <- paste0("/home/martin/nes/nes_RAD/fastsimcoal_analyses/fastsimcoal_files/", name_run)
if (!dir.exists(folder_run)) dir.create(folder_run)

library(stringr)

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
      "1 NCUR logunif 1000 50000 output", 
      "1 NANC logunif 1000 50000 output", 
      "1 NBOT unif 1 30 output",
      "[RULES]", 
      "[COMPLEX PARAMETERS]",  
      "0 RESBOT = NBOT/NCUR hide", 
      "0 RESENDBOT = NANC/NBOT hide"),
    sep = "\n"
)
sink()

