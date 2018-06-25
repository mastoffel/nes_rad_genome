# Quantifying inbreeding
# Filter in a variety of ways and visualise inbreeding ~ missingness
# sMLH, IBCS, relatedness (IBD)

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
options(scipen=999)
library(reshape2)
library(inbreedR)
#source("scripts/mkTrend.R")

system("mkdir data/inbreeding")

#~~ Filter raw vcf file

#~~~~~~~~~~~~~~~~~~~~~~#
#  maf 0.05, 60% geno  #
#~~~~~~~~~~~~~~~~~~~~~~# 

system("vcftools --vcf data/ArcGaz_biallelic_sub.recode.vcf --out data/inbreeding/maf05_g60 --maf 0.05 --max-missing 0.60 --recode --recode-INFO-all")
system("vcftools --vcf data/inbreeding/maf05_g60.recode.vcf --plink --out data/inbreeding/maf05_g60")
system("plink --file data/inbreeding/maf05_g60 --make-bed --recodeAD --out data/inbreeding/maf05_g60")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  maf 0.05, 60% geno, GQ 5, DP 8  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

system("vcftools --vcf data/ArcGaz_biallelic_sub.recode.vcf --geno-depth --out data/inbreeding/raw")

depth <- fread("data/inbreeding/raw.gdepth", colClasses = "numeric") %>%
  unite(SNP, CHROM, POS)
depth <- select(depth, -SNP)

depth <- depth %>%
  mutate(mean = rowMeans(.), 
         sum = rowSums(.), 
         logmeandepth = log10(mean))

CI <- 0.90
CI_meanDP <- stats::quantile(depth$mean, c((1-CI)/2,1-(1-CI)/2), na.rm=TRUE)
CI_sumDP <- stats::quantile(depth$sum, c((1-CI)/2,1-(1-CI)/2), na.rm=TRUE)
CI_meanDP
CI_sumDP

hist(depth$logmeandepth)

# ** See filters below: incorrect 

system("vcftools --vcf data/ArcGaz_biallelic_sub.recode.vcf --out data/inbreeding/GQ_DP --minGQ 5 --minDP 8 --maxDP 30 --recode --recode-INFO-all")

system("vcftools --vcf data/inbreeding/GQ_DP.recode.vcf --out data/inbreeding/GQ_DP_maf05_g30 --maf 0.05 --max-missing 0.30 --recode --recode-INFO-all")

system("vcftools --vcf data/inbreeding/GQ_DP.recode.vcf --plink --out data/inbreeding/GQ_DP_maf05_g30")

#

# recode and reorder bim
system("scripts/recode_full_map.sh data/inbreeding/GQ_DP_maf05_g30.map data/inbreeding/GQ_DP_maf05_g30.map")
system("plink --file data/inbreeding/GQ_DP_maf05_g30 --make-bed --out data/inbreeding/GQ_DP_maf05_g30 --allow-extra-chr")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#            Remove individuals with high amounts of missing data          #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Identify individuals with lots of missing data

system("vcftools --vcf data/inbreeding/maf05_g60.recode.vcf --missing-indv --out data/inbreeding/maf05_g60")

system("vcftools --vcf data/inbreeding/GQ_DP_maf05_g30.recode.vcf --missing-indv --out data/inbreeding/GQ_DP_maf05_g30")
#system("vcftools --vcf data/inbreeding/GQ_DP.recode.vcf --missing-indv --out data/inbreeding/GQ_DP_maf05_g30")


# make histogram

import_files <- function(file){
  fread(file)
}

IM_files <- paste("data/inbreeding/", list.files(path = "data/inbreeding", pattern="*.imiss"), sep = "")
IM <- lapply(IM_files, import_files)
names(IM) <- IM_files

hist_fun <- function(x){
  hist(x$F_MISS)
}

par(mfrow = c(1, 2))
histograms <- lapply(IM, hist_fun)
par(mfrow = c(1, 1))

# 40% maf05_g60
# 70% GQ_DP_maf05_g30

system("mawk '$5 > 0.4' data/inbreeding/maf05_g60.imiss | cut -f1 > data/inbreeding/maf05_g60_lowDP.indv")
system("mawk '$5 > 0.7' data/inbreeding/GQ_DP_maf05_g30.imiss | cut -f1 > data/inbreeding/GQ_DP_maf05_g30_lowDP.indv")

# feed to vcftools
system("vcftools --vcf data/inbreeding/maf05_g60.recode.vcf --remove data/inbreeding/maf05_g60_lowDP.indv --recode --recode-INFO-all --out data/inbreeding/maf05_g60_ldi")
system("vcftools --vcf data/inbreeding/GQ_DP_maf05_g30.recode.vcf --remove data/inbreeding/GQ_DP_maf05_g30_lowDP.indv --recode --recode-INFO-all --out data/inbreeding/GQ_DP_maf05_g30_ldi")

#system("vcftools --vcf data/inbreeding/GQ_DP.recode.vcf --remove data/inbreeding/GQ_DP_maf05_g30_lowDP.indv --recode --recode-INFO-all --out data/inbreeding/GQ_DP_maf05_g30_ldi")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#       Restrict Dataset        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Now that poor coverage individuals have been removed:
# Remove mendel errors
# Keep variants called in a high percentage of individuals


#~~ Mendel errors

system("vcftools --vcf data/inbreeding/maf05_g60_ldi.recode.vcf --exclude-positions data/mendel/MendelErrorStrict_SNPs_tab.txt --recode --recode-INFO-all --out data/inbreeding/maf05_g60_ldi_ME")
system("vcftools --vcf data/inbreeding/GQ_DP_maf05_g30_ldi.recode.vcf --exclude-positions data/mendel/MendelErrorStrict_SNPs_tab.txt --recode --recode-INFO-all --out data/inbreeding/GQ_DP_maf05_g30_ldi_ME")


#~~ Filter for SG individuals here

system("vcftools --vcf data/inbreeding/GQ_DP_maf05_g30_ldi_ME.recode.vcf --keep data/raw/SSB_IDs.txt --recode --recode-INFO-all --out data/inbreeding/GQ_DP_maf05_g30_ldi_SSB")
system("bcftools query -l data/inbreeding/GQ_DP_maf05_g30_ldi_SSB.recode.vcf | wc -l ")

#~~ Genotyping rate

system("vcftools --vcf data/inbreeding/maf05_g60_ldi_ME.recode.vcf --max-missing 0.9 --maf 0.05 --recode --recode-INFO-all --out data/inbreeding/maf05_g90_ldi_ME")
system("vcftools --vcf data/inbreeding/maf05_g90_ldi_ME.recode.vcf --plink --out data/inbreeding/maf05_g90_ldi_ME")
system("plink --file data/inbreeding/maf05_g90_ldi_ME --make-bed --recodeAD --out data/inbreeding/maf05_g90_ldi_ME")

# for this dataset, 30 % is necessary to keep a high number of SNPs
system("vcftools --vcf data/inbreeding/GQ_DP_maf05_g30_ldi_SSB.recode.vcf --max-missing 0.75 --maf 0.02 --recode --recode-INFO-all --out data/inbreeding/GQ_DP_maf05_G_ldi_ME")
system("bcftools query -l data/inbreeding/GQ_DP_maf05_G_ldi_ME.recode.vcf | wc -l ")

system("vcftools --vcf data/inbreeding/GQ_DP_maf05_G_ldi_ME.recode.vcf --plink --out data/inbreeding/GQ_DP_maf05_G_ldi_ME")
system("wc -l data/inbreeding/GQ_DP_maf05_G_ldi_ME.ped")
system("plink --file data/inbreeding/GQ_DP_maf05_G_ldi_ME --recodeAD --make-bed --out data/inbreeding/GQ_DP_maf05_G_ldi_ME")
system("wc -l data/inbreeding/GQ_DP_maf05_G_ldi_ME.fam")


#~~ Prepare file for measuring relatedness
# High overlap, and high maf -- informative

system("vcftools --vcf data/inbreeding/GQ_DP_maf05_g30_ldi_ME.recode.vcf --max-missing 0.8 --maf 0.2 --recode --recode-INFO-all --out data/inbreeding/relatedness")
system("vcftools --vcf data/inbreeding/relatedness.recode.vcf --plink --out data/inbreeding/relatedness")
system("plink --file data/inbreeding/relatedness --recodeAD --make-bed --out data/inbreeding/relatedness")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#         Linkage pruning          #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# lapply these functions

# recode map file for LD pruning
system("scripts/recode_full_map.sh data/inbreeding/maf05_g90_ldi_ME.map data/inbreeding/maf05_g90_ldi_ME.map")
system("scripts/recode_full_map.sh data/inbreeding/GQ_DP_maf05_G_ldi_ME.map data/inbreeding/GQ_DP_maf05_G_ldi_ME.map")
system("scripts/recode_full_map.sh data/inbreeding/relatedness.map data/inbreeding/relatedness.map")


# Make bed
system("plink --file data/inbreeding/maf05_g90_ldi_ME --make-bed --out data/inbreeding/maf05_g90_ldi_ME --allow-extra-chr --debug")
system("plink --file data/inbreeding/GQ_DP_maf05_G_ldi_ME --make-bed --out data/inbreeding/GQ_DP_maf05_G_ldi_ME --allow-extra-chr --debug")
system("plink --file data/inbreeding/relatedness --make-bed --out data/inbreeding/relatedness --allow-extra-chr --debug")

# Run PLINK --indep
system("plink --bfile data/inbreeding/maf05_g90_ldi_ME --indep 50 5 2 --nonfounders --out data/inbreeding/maf05_g90_ldi_ME --allow-extra-chr --debug")
system("wc -l data/inbreeding/maf05_g90_ldi_ME.prune.in")

# Run PLINK --indep
system("plink --bfile data/inbreeding/GQ_DP_maf05_G_ldi_ME --indep 50 5 2 --nonfounders --out data/inbreeding/GQ_DP_maf05_G_ldi_ME --allow-extra-chr --debug")
system("wc -l data/inbreeding/GQ_DP_maf05_G_ldi_ME.prune.in")

# Run PLINK --indep
system("plink --bfile data/inbreeding/relatedness --indep 50 5 2 --nonfounders --out data/inbreeding/relatedness --allow-extra-chr --debug")
system("wc -l data/inbreeding/relatedness.prune.in")

# LD filter
# Make plink raw file
system("plink --bfile data/inbreeding/maf05_g90_ldi_ME --extract data/inbreeding/maf05_g90_ldi_ME.prune.in --out data/inbreeding/maf05_g90_ldi_ME_LD --recodeAD --allow-extra-chr --debug")
# Make standard plink files
system("plink --bfile data/inbreeding/maf05_g90_ldi_ME --extract data/inbreeding/maf05_g90_ldi_ME.prune.in --out data/inbreeding/maf05_g90_ldi_ME_LD --make-bed --recode --allow-extra-chr --debug")


# Make plink raw file
system("plink --bfile data/inbreeding/GQ_DP_maf05_G_ldi_ME --extract data/inbreeding/GQ_DP_maf05_G_ldi_ME.prune.in --out data/inbreeding/GQ_DP_maf05_G_ldi_ME_LD --recodeAD --allow-extra-chr --debug")
# Make standard plink files
system("plink --bfile data/inbreeding/GQ_DP_maf05_G_ldi_ME --extract data/inbreeding/GQ_DP_maf05_G_ldi_ME.prune.in --out data/inbreeding/GQ_DP_maf05_G_ldi_ME_LD --make-bed --recode --allow-extra-chr --debug")



# Make plink raw file
system("plink --bfile data/inbreeding/relatedness --extract data/inbreeding/relatedness.prune.in --out data/inbreeding/relatedness_lD --recodeAD --allow-extra-chr --debug")
# Make standard plink files
system("plink --bfile data/inbreeding/relatedness --extract data/inbreeding/relatedness.prune.in --out data/inbreeding/relatedness_lD --make-bed --recode --allow-extra-chr --debug")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#      Inbreeding coefficients     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Get sMLH values

get_sMLH_from_plinkraw <- function(file) {
  
  library(data.table)
  library(dplyr)
  x <- fread(file, colClasses = "character")
  ids <- x$IID
  x <- dplyr::select(x, contains("HET"))
  row.names(x) <- ids
  NAs <- apply(x, 1, function(x) sum(is.na(x)))
  
  sMLH <- as.data.frame(MLH(x))
  sMLH$ANIMAL <- ids
  sMLH$NAS <- NAs
  sMLH <- dplyr::filter(sMLH, grepl("^AG", ANIMAL))
  colnames(sMLH) <- c("sMLH", "ANIMAL", "NAs")
  sMLH <- sMLH
  
}


raw_files <- paste("data/inbreeding/", list.files(path = "data/inbreeding", pattern="*.raw"), sep = "")
sMLH <- lapply(raw_files[c(1,4)], get_sMLH_from_plinkraw)
names(sMLH) <- lapply(raw_files[c(1,4)], function(x) gsub("data/inbreeding/|\\.raw", "", x))

sMLH <- rbindlist(sMLH, idcol = "Run") %>%
  dplyr::rename(IID = ANIMAL)

#~~ Load microsatellite data

ms <- fread("data/raw/Rack61_microsatellites_Jul2017.csv", header = T) %>%
  right_join(filter(sMLH, Run == "maf05_g90_ldi_ME_LD"), by = c("ID" = "IID"))
IDs <- ms$ID
rownames(ms) <- ms$ID

ms <- select(ms, -Run, -ID, -sMLH, -NAs)
raw_ms <- convert_raw(ms)
ms_sMLH <- data.frame(MLH(raw_ms))
ms_sMLH$IID <- IDs
colnames(ms_sMLH) <- c("ms_sMLH", "IID")
rm(ms)
rm(raw_ms)

#~~ g2 (boot over loci)

get_g2_from_plinkraw <- function(file) {
  
  x <- fread("data/inbreeding/GQ_DP_maf05_G_ldi_ME_LD.raw", colClasses = "numeric") %>%
    filter(grepl("^AG", IID))
  
  ids <- x$IID
  x <- dplyr::select(x, contains("HET"))
  row.names(x) <- ids
  
  library(inbreedR)
  g2 <- g2_snps(x, nperm = 1000, nboot = 1000, CI = 0.95)
  g2
  g2 <- g2
  
}

g2 <- lapply(raw_files[1], get_g2_from_plinkraw)
g2 <- g2[[1]]
plot(g2, col = "grey")
g2
var(filter(sMLH, Run == "GQ_DP_maf05_G_ldi_ME_LD")$sMLH)



#~~ Get Fhats

# recode bim files for GCTA

recode_bim_1chr <- function(file){
  file <- fread(file) %>%
    mutate(V1 = 1) %>%
    fwrite(file, quote = F, row.names = F,
           col.names = F, sep = " ")
}

bim_files <- paste("data/inbreeding/", list.files(path = "data/inbreeding", pattern="*.bim"), sep = "")
lapply(bim_files, recode_bim_1chr)

# get fhats using gcta

plink_files <- paste("data/inbreeding/", list.files(path = "data/inbreeding", pattern="*.ped"), sep = "")
plink_files <- lapply(plink_files, function(x) gsub(".ped", "", x))

for (i in 1:length(plink_files)){
  
  #system(paste0("~/programs/gcta64 --bfile ", plink_files[i]," --autosome --ibc --out ", plink_files[i]," --thread-num 10"))
  
  system(paste0("gcta64 --bfile ", plink_files[i]," --autosome --ibc --out ", plink_files[i]," --thread-num 10"))
  
}

# load fhats

load_fhats <- function(file) {
  
  fhats <- fread(file, header = T) %>%
    filter(grepl("^AG", IID)) %>%
    mutate(Animal = ifelse(grepl("^AGP", IID), "Pup", "Adult"))
  # filter(grepl("^AG", IID)) %>%
  # mutate(Animal = ifelse(grepl("^AGP", IID), "Pup", "Adult"))
  
  #%>% mutate(Animal = ifelse(grepl("^AGP", IID), "Pup",
  #                       ifelse(grepl("^AGM", IID), "Father", "Mother")))
  
}


ibc_files <- paste("data/inbreeding/", list.files(path = "data/inbreeding", pattern="*.ibc"), sep = "")
fhats <- lapply(ibc_files[c(1,5)], load_fhats)
names(fhats) <- lapply(ibc_files[c(1,5)], function(x) gsub("data/inbreeding/|\\.ibc", "", x))

ibcs <- rbindlist(fhats, idcol = "Run") %>%
  left_join(sMLH, by = c("IID", "Run")) %>%
  left_join(ms_sMLH, by  = "IID")

# ~~ Plot g2

library(ggthemr)

ggthemr(palette = "pale", layout = "clean", 
        line_weight = 0.7, text_size = 20, type = "outer")
swatch()
to_swap <- swatch()[3:4]

g2$g2 / var(filter(ibcs, Run == "GQ_DP_maf05_G_ldi_ME_LD")$sMLH)
g2$g2 / var(filter(ibcs, Run == "GQ_DP_maf05_G_ldi_ME_LD")$Fhat1)
g2$g2 / var(filter(ibcs, Run == "GQ_DP_maf05_G_ldi_ME_LD")$Fhat2)
g2$g2 / var(filter(ibcs, Run == "GQ_DP_maf05_G_ldi_ME_LD")$Fhat3)

View(filter(ibcs, Run == "GQ_DP_maf05_G_ldi_ME_LD"))


g2_plot <- data.frame(g2$g2_boot)
lcl <- g2$CI_boot[1]
ucl <- g2$CI_boot[2]
g2_boot_summary <- data.frame(lcl, ucl)

ibcs_vars_summary <- data.frame(c(g2$g2,
                                  var(filter(ibcs, Run == "GQ_DP_maf05_G_ldi_ME_LD")$sMLH),
                                  var(filter(ibcs, Run == "GQ_DP_maf05_G_ldi_ME_LD")$Fhat1),
                                  var(filter(ibcs, Run == "GQ_DP_maf05_G_ldi_ME_LD")$Fhat2),
                                  var(filter(ibcs, Run == "GQ_DP_maf05_G_ldi_ME_LD")$Fhat3)),
                                c("g2", "sMLH", "Fhat1", "Fhat2", "Fhat3"))

colnames(ibcs_vars_summary) <- c("val", "var")

plot(g2)


# g2 bootstrapping distribution showing empirical g2 with CIs

require(gridExtra)
library(sitools)

png("figs/g2_boot.png", units = "in", res = 300, width = 8, height = 7)

cbPalette <- c( "#1B9E77", "#66A61E", "#E6AB02", "black", "#7570B3", "#D95F02", "#E7298A")

g2_CI_plot <- ggplot(g2_plot, aes(g2$g2_boot)) + 
  geom_histogram(colour = "grey45", fill = "grey45") +
  geom_errorbarh(aes(xmin = g2_boot_summary$lcl , xmax = g2_boot_summary$ucl , y = 90),
                 size = 0.8, color = "black", linetype = "solid", height = 0) +
  geom_vline(data = ibcs_vars_summary, aes(xintercept = val, colour = var), size = 0.8,
             linetype = c("dashed", "solid", "solid", "solid", "solid")) +
  scale_colour_manual(values = cbPalette, name = "",
                      breaks = c("g2", "sMLH", "Fhat1", "Fhat2", "Fhat3"),
                      labels = c(expression(italic(g[2])), 
                                 expression("var"(italic(hat(F)["I"]))), 
                                 expression("var"(italic(hat(F)["II"]))), 
                                 expression("var"(italic(hat(F)["III"]))), 
                                 expression("var"("sMLH")))) +
  labs(y = "Counts", x = expression(italic(g[2]))) +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain")) +
  theme(legend.position = "right")

dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#          Adult Pup Distribution          #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Adult pup distribution

#library(RColorBrewer)
#cbPalette <- c(brewer.pal(6, "Dark2"))

png("figs/ibc_pup_adult.png", units = "in", res = 300, width = 11, height = 8)

adult_pup_fhat3 <- ggplot(filter(ibcs, Run == "GQ_DP_maf05_G_ldi_ME_LD"), aes(x=Fhat3, fill=Animal, color = Animal)) +
  #geom_histogram(aes(y=..density..),  position="identity", alpha=0.9, col = "#3262AB", fill = "#3262AB") +
  geom_density(alpha=0.8) +
  labs(y = "Density", x = expression(italic(hat(F)["III"]))) +
  xlim(c(-0.2,0.2)) +
  #xlim(c(0.75, 1.3)) +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain")) +
  scale_color_manual(values=c("black", "black")) +
  # scale_color_manual(values=c("grey33", "#EB9050")) +
  scale_fill_manual(values=c("grey33", "#EB9050"))

dev.off()



png("figs/ibc_and_g2.png", units = "in", res = 300, width = 13, height = 6)
grid.arrange(adult_pup_fhat3, g2_CI_plot,
             ncol=2, nrow=1)
dev.off()


#~~ Permutation test

pups <- filter(filter(ibcs, Run == "GQ_DP_maf05_G_ldi_ME_LD"), Animal == "Pup")$Fhat3
adults <- filter(filter(ibcs, Run == "GQ_DP_maf05_G_ldi_ME_LD"), Animal == "Adult")$Fhat3
obsdiff <- mean(pups)-mean(adults)

N <- 22

avgdiff <- replicate(10000, {
  all <- sample(c(pups,adults)) # scrambles data
  newpups <- all[1:N] # first half of data
  newadults <- all[(N+1):(2*N)] # second half of data
  return(mean(newpups) - mean(newadults)) # new mean
})

hist(avgdiff)
abline(v=obsdiff, col="red", lwd=2)
(sum(abs(avgdiff) > abs(obsdiff)) + 1) / (length(avgdiff) + 1)

median(pups)
median(adults)

#~~~~~~~~~~~~~~~~~~~~~~#
#     Joy plots        #
#~~~~~~~~~~~~~~~~~~~~~~#

melt_ibcs <- melt(filter(ibcs, Run == "GQ_DP_maf05_G_ldi_ME_LD"), id = c("Run", "FID", "IID", "NOMISS", "Animal", "NAs"))
library(ggjoy)

png("figs/joy_plot_ibcs.png", units = "in", res = 300, width = 11, height = 8)
ggplot(melt_ibcs, aes(x = value, y = variable, fill = variable)) + geom_joy(scale = 5, alpha = 0.5) +
  scale_fill_manual(values=c('#EB9050', 'grey33', '#EB9050', 'grey33')) +
  labs(y = "Inbreeding level", x = "Inbreeding coefficient") +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain")) 
dev.off()


#~~~~~~~~~~~~~~~~~~~~~~#
#    NA bias plots     #
#~~~~~~~~~~~~~~~~~~~~~~#


# Plots for presentations and Supplementary
# Presentation colour = #3262AB

plot1 <- ggplot(filter(ibcs, Run == "maf05_g90_ldi_ME_LD"), aes(x=Fhat3, color=.id, fill=.id), col = "grey33") +
  geom_histogram(aes(y=..density..),  position="identity", alpha=0.9, col = "grey33", fill = "grey33") +
  geom_density(alpha=0.6, col = "grey33", fill = "grey33") +
  labs(y = "Density", x = expression(italic(hat(F)["III"]))) +
  xlim(c(-0.5,0.5)) +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain"))
# axis.title = element_text(color = "black"))

plot2 <- ggplot(filter(ibcs, Run == "maf05_g90_ldi_ME_LD"), aes(factor(Animal), Fhat1)) + 
  geom_violin(aes(fill = factor(Animal))) +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain"))

plot3 <- ggplot(filter(ibcs, Run == "maf05_g90_ldi_ME_LD"), aes(x=NOMISS, y = Fhat3)) +
  geom_point(size = 2, alpha = 1/1.5, col = "grey33") +
  labs(x = "SNPs Genotyped", y = expression(italic(hat(F)["III"]))) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  scale_x_continuous(labels = f2si) +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain"))

# Depth filtered 
library(scales)

plot4 <- ggplot(filter(ibcs, Run == "GQ_DP_maf05_G_ldi_ME_LD"), aes(x=Fhat3, color=.id, fill=.id), col = "grey33") +
  geom_histogram(aes(y=..density..),  position="identity", alpha=0.9, col = "grey33", fill = "grey33") +
  geom_density(alpha=0.6, col = "grey33", fill = "grey33") + #3262AB
  labs(y = "Density", x = expression(italic(hat(F)["III"]))) +
  #xlim(c(-0.5,0.5)) +
  xlim(c(-0.3,0.3)) +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain")) +
  ggtitle('A') + theme(plot.title=element_text(hjust=0, size = 12))

plot5 <- ggplot(filter(ibcs, Run == "GQ_DP_maf05_G_ldi_ME_LD"), aes(factor(Animal), Fhat3)) + 
  geom_violin(aes(fill = factor(Animal))) +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain"))

plot6 <- ggplot(filter(ibcs, Run == "GQ_DP_maf05_G_ldi_ME_LD"), aes(x=NOMISS, y = Fhat3)) +
  geom_point(size = 2, alpha = 1/1.5, col = "grey33") +
  labs(x = "SNPs Genotyped", y = expression(italic(hat(F)["III"]))) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  scale_x_continuous(labels = f2si) +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain"))

#~~~~~~~~~~~~~~~~~~~~~~~~~#
#   Correlation plots     #
#~~~~~~~~~~~~~~~~~~~~~~~~~#

# For Supplementary ?

#lm_eqn = function(m) {
#  eq <- substitute(italic(r)^2~"="~r2, 
#                   list(r2 = format(summary(m)$r.squared, digits = 3)))
#  as.character(as.expression(eq));   
#}

corr_eqn <- function(x,y, digits = 2) {
  corr_coef <- round(cor(x, y), digits = digits)
  paste("italic(r) == ", corr_coef)
}


#m1 <- lm(Fhat1 ~ Fhat3, filter(ibcs, Run == "GQ_DP_maf05_G_ldi_ME_LD"))
labels1 = data.frame(x = -0.04, y = 0.15, label = corr_eqn(filter(ibcs, Run == "GQ_DP_maf05_G_ldi_ME_LD")$Fhat1, 
                                                           filter(ibcs, Run == "GQ_DP_maf05_G_ldi_ME_LD")$Fhat3))

plot7 <- ggplot(filter(ibcs, Run == "GQ_DP_maf05_G_ldi_ME_LD"), aes(x=Fhat3, y = Fhat1)) +
  geom_point(size = 2, alpha = 1/1.5, col = "grey33") + #3262AB
  labs(x = expression(italic(hat(F)["III"])), y = expression(italic(hat(F)["I"]))) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain")) +
  # geom_text(aes(x = -0.03, y = 0.15, label = lm_eqn(m1)), parse = TRUE,
  #           size = 5, fontface = "plain") +
  geom_text(data = labels1, aes(x = x, y = y, label = label), parse = TRUE,
            size = 5, fontface = "plain") +
  ggtitle('B') + theme(plot.title=element_text(hjust=0, size = 12))


labels2 = data.frame(x = -0.04, y = 0.15, label = corr_eqn(filter(ibcs, Run == "GQ_DP_maf05_G_ldi_ME_LD")$Fhat2, 
                                                           filter(ibcs, Run == "GQ_DP_maf05_G_ldi_ME_LD")$Fhat3))

plot8 <- ggplot(filter(ibcs, Run == "GQ_DP_maf05_G_ldi_ME_LD"), aes(x=Fhat3, y = Fhat2)) +
  geom_point(size = 2, alpha = 1/1.5, col = "grey33") +
  labs(x = expression(italic(hat(F)["III"])), y = expression(italic(hat(F)["II"]))) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain")) +
  geom_text(data = labels2, aes(x = x, y = y, label = label), parse = TRUE,
            size = 5, fontface = "plain") +
  ggtitle('C') + theme(plot.title=element_text(hjust=0, size = 12))


labels3 = data.frame(x = 0.1, y = 0.29, label = corr_eqn(filter(ibcs, Run == "GQ_DP_maf05_G_ldi_ME_LD")$sMLH, 
                                                         filter(ibcs, Run == "GQ_DP_maf05_G_ldi_ME_LD")$Fhat3))
plot9 <- ggplot(filter(ibcs, Run == "GQ_DP_maf05_G_ldi_ME_LD"), aes(x=Fhat3, y = sMLH)) +
  geom_point(size = 2, alpha = 1/1.5, col = "grey33") +
  labs(x = expression(italic(hat(F)["III"])), y = "sMLH") +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain")) +
  geom_text(data = labels3, aes(x = x, y = y, label = label), parse = TRUE,
            size = 5, fontface = "plain") +
  ggtitle('D') + theme(plot.title=element_text(hjust=0, size = 12))


labels4 = data.frame(x = 0.1, y = 0.29, label = corr_eqn(filter(ibcs, Run == "GQ_DP_maf05_G_ldi_ME_LD")$sMLH, 
                                                         filter(ibcs, Run == "GQ_DP_maf05_G_ldi_ME_LD")$Fhat2))
plot10 <- ggplot(filter(ibcs, Run == "GQ_DP_maf05_G_ldi_ME_LD"), aes(x=Fhat2, y = sMLH)) +
  geom_point(size = 2, alpha = 1/1.5, col = "grey33") +
  labs(x = expression(italic(hat(F)["II"])), y = "sMLH") +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain")) +
  geom_text(data = labels4, aes(x = x, y = y, label = label), parse = TRUE,
            size = 5, fontface = "plain")

labels5 = data.frame(x = 0.285, y = 0.9, label = corr_eqn(filter(ibcs, Run == "GQ_DP_maf05_G_ldi_ME_LD")$sMLH, 
                                                          filter(ibcs, Run == "GQ_DP_maf05_G_ldi_ME_LD")$ms_sMLH))
plot11 <- ggplot(filter(ibcs, Run == "GQ_DP_maf05_G_ldi_ME_LD"), aes(x=sMLH, y = ms_sMLH)) +
  geom_point(size = 2, alpha = 1/1.5, col = "grey33") +
  labs(x = "RAD sMLH", y = "Microsatellite sMLH") +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain")) +
  geom_text(data = labels5, aes(x = x, y = y, label = label), parse = TRUE,
            size = 5, fontface = "plain")


labels6 = data.frame(x = 0.3, y = 0.15, label = corr_eqn(filter(ibcs, Run == "GQ_DP_maf05_G_ldi_ME_LD")$Fhat3, 
                                                         filter(ibcs, Run == "GQ_DP_maf05_G_ldi_ME_LD")$ms_sMLH))
plot12 <- ggplot(filter(ibcs, Run == "GQ_DP_maf05_G_ldi_ME_LD"), aes(x=ms_sMLH, y = Fhat3)) +
  geom_point(size = 2, alpha = 1/1.5, col = "grey33") +
  labs(y = expression(italic(hat(F)["III"])), x = "Microsatellite sMLH") +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain")) +
  geom_text(data = labels6, aes(x = x, y = y, label = label), parse = TRUE,
            size = 5, fontface = "plain")


png("figs/ibc_miss.png", units = "in", res = 300, width = 12, height = 8)
grid.arrange(plot1, plot3,
             plot4, plot6, ncol=2, nrow=2)
dev.off()


png("figs/ibc_correlations.png", units = "in", res = 300, width = 10, height = 9) # 22, 6
grid.arrange(plot4, plot7,
             plot8, plot9, ncol=2, nrow=2)
dev.off()








#~~ Change in inbreeding through time


cbPalette <- c("grey33", "#EB9050")

test <- filter(ibcs, Run == "GQ_DP_maf05_G_ldi_ME_LD") %>%
  mutate(Year = substr(IID, 4, 5)) %>%
  mutate(Year_pre = case_when(grepl("^0", .$Year) ~ "20",
                              grepl("^9", .$Year) ~ "19")) %>%
  unite("YEAR", c("Year_pre", "Year"), sep = "") %>%
  mutate(YEAR = as.numeric(YEAR))



png("figs/ibc_time.png", units = "in", res = 300, width = 11, height = 8)

ggplot(test, aes(YEAR, Fhat3, col = Animal)) +
  geom_point(aes(colour = factor(Animal)), size = 6, alpha = 0.6) +
  #geom_smooth(method = "lm", se = FALSE, col = "grey33") +
  #  geom_smooth(method = "lm", se = FALSE, col = "grey33", size = 2) +
  geom_smooth(method = "lm", se = FALSE, size = 2) +
  # geom_smooth(span = 1, colour = "grey30", fill = "grey") +
  # geom_smooth(span = 1) +
  scale_colour_manual(values = cbPalette,
                      name = "Age") +
  labs(y = expression(italic(hat(F)["III"])), x = "Birth Year") +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain")) 

dev.off()


# Mean of year
mean_year <- test %>%
  dplyr::group_by(YEAR) %>%
  dplyr::summarise(mean = mean(Fhat3),
                   CI_lower = stats::quantile(Fhat3, c((1-CI)/2,1-(1-CI)/2), na.rm=TRUE)[[1]],
                   CI_upper = stats::quantile(Fhat3, c((1-CI)/2,1-(1-CI)/2), na.rm=TRUE)[[2]])


ggplot(mean_year, aes(YEAR, mean)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_manual(values = cbPalette,
                      name = "Age") +
  geom_ribbon(data=mean_year,aes(ymin=CI_lower,ymax=CI_upper),alpha=0.3) +
  labs(y = expression(italic(hat(F)["III"])), x = "Year") +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain")) +
  theme(axis.text.x = element_text(angle = 70, hjust = 1))




# Analysis

gm1 <- glm(mean ~ YEAR, data = mean_year)
gm0 <- glm(mean ~ 1, data = mean_year)
anova(gm0, gm1)
AIC(gm0, gm1)
coef(gm1)

# investigate residual correlations
acf(resid(gm1), main = "Residuals of linear model")

layout(matrix(1:4, ncol = 2))
plot(gm1)
layout(1)

library(nlme)

gg0 <- gls(mean ~ 1, data = mean_year, method = "ML")
gg1 <- gls(mean ~ YEAR, data = mean_year, method = "ML")
anova(gg0, gg1)


gg1 <- update(gg1, method = "REML")
gg2 <- gls(mean ~ YEAR, data = mean_year,
           correlation = corARMA(form = ~ YEAR, p = 1), method = "REML")
gg3 <- gls(mean ~ YEAR, data = mean_year,
           correlation = corARMA(form = ~ YEAR, p = 2), method = "REML")

anova(gg1, gg2, gg3)
confint(gg1)

# Visualise vitted trend
pred <- data.frame(YEAR = 1994:2003)
pred <- transform(pred, yhat = predict(gg1, newdata = pred))
with(pred, yhat)

# Plot

ggplot(mean_year, aes(YEAR, mean)) +
  geom_point(size = 3) +
  # geom_smooth(method = "lm", se = FALSE) +
  geom_line(data = pred, 
            aes(x = YEAR, y = yhat), size = 1.5) +
  geom_line(data = pred, 
            aes(x = YEAR, y = yhat), size = 1.5, col  = "red") +
  scale_colour_manual(values = cbPalette,
                      name = "Age") +
  labs(y = expression(italic(hat(F)["III"])), x = "Year") +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain")) +
  theme(axis.text.x = element_text(angle = 70, hjust = 1))

# Heritability of heterozygosity

head(ibcs)

pedigree <- read.table("data/raw/new_pedigree.txt", header = T) %>%
  filter(grepl("^AGP", ANIMAL))

mother_ibcs <- filter(ibcs, Run == "GQ_DP_maf05_G_ldi_ME_LD") %>%
  filter(grepl("^AGF", IID))

father_ibcs <- filter(ibcs, Run == "GQ_DP_maf05_G_ldi_ME_LD") %>%
  filter(grepl("^AGM", IID))

her_het <- filter(ibcs, Run == "GQ_DP_maf05_G_ldi_ME_LD") %>%
  left_join(pedigree, by = c("IID" = "ANIMAL")) %>%
  filter(grepl("^AGP", IID)) %>%
  left_join(mother_ibcs,by = c("MOTHER" = "IID")) %>%
  left_join(father_ibcs, by = c("FATHER" = "IID")) %>%
  select(IID, Fhat3.x, sMLH.x, MOTHER, Fhat3.y, sMLH.y, FATHER, Fhat3, sMLH)

colnames(her_het) <- c("PUP", "PUP_Fhat3", "PUP_sMLH", 
                       "MOTHER", "MOTHER_Fhat3", "MOTHER_sMLH",
                       "FATHER", "FATHER_Fhat3", "FATHER_sMLH")

her_het$MID_Fhat3 <- rowMeans(her_het[,c("MOTHER_Fhat3", "FATHER_Fhat3")], na.rm=TRUE)
her_het$MID_sMLH <- rowMeans(her_het[,c("MOTHER_sMLH", "FATHER_sMLH")], na.rm=TRUE)

library(stringr)

df <- her_het %>%
  mutate(Year = substr(PUP, 4, 5)) %>%
  mutate(Year_pre = case_when(grepl("^0", .$Year) ~ "20",
                              grepl("^9", .$Year) ~ "19")) %>%
  unite("YEAR", c("Year_pre", "Year"), sep = "") %>%
  mutate(YEAR = as.numeric(YEAR))


plot(df$MOTHER_Fhat3 ~ df$YEAR) +
  abline(lm(df$MOTHER_Fhat3 ~ df$YEAR))
plot(df$FATHER_Fhat3 ~ df$YEAR) +
  abline(lm(df$FATHER_Fhat3 ~ df$YEAR))
plot(df$MOTHER_sMLH ~ df$YEAR) +
  abline(lm(df$MOTHER_sMLH ~ df$YEAR))
plot(df$FATHER_sMLH ~ df$YEAR) +
  abline(lm(df$FATHER_sMLH ~ df$YEAR))

plot(her_het$PUP_Fhat3 ~ her_het$MOTHER_Fhat3)
plot(her_het$PUP_Fhat3 ~ her_het$FATHER_Fhat3)
plot(her_het$PUP_sMLH ~ her_het$MOTHER_sMLH)
plot(her_het$PUP_sMLH ~ her_het$FATHER_sMLH)

plot(her_het$PUP_Fhat3 ~ her_het$MID_Fhat3) +
  abline(lm(her_het$PUP_Fhat3 ~ her_het$MID_Fhat3))

plot(her_het$PUP_sMLH ~ her_het$MID_sMLH) +
  abline(lm(her_het$PUP_sMLH ~ her_het$MID_sMLH))







# Method 2






mean_year_all <- test %>%
  group_by(YEAR) %>%
  summarise(mean = mean(Fhat3))

ggplot(mean_year_all, aes(YEAR, mean)) +
  geom_point(size = 3) +
  #geom_line() +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_manual(values = cbPalette,
                      name = "Age") +
  labs(y = expression(italic(hat(F)["III"])), x = "Year") +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain")) +
  theme(axis.text.x = element_text(angle = 70, hjust = 1))

#~~ Time series analysis

summary(lm(test$Fhat3 ~ test$YEAR))

install.packages("Kendall")
library(Kendall)
MannKendall(test$Fhat3)
MannKendall(filter(test, Animal == "Pup")$Fhat3)
MannKendall(filter(test, Animal == "Adult")$Fhat3)
MannKendall(filter(mean_year, Animal == "Adult")$mean)

library(pastecs)
trend.test(mean_year_all$mean, R=1)
MannKendall(mean_year_all$mean)
mkTrend(mean_year_all$mean)

mean_year_pup <- test %>%
  group_by(YEAR) %>%
  filter(Animal == "Pup") %>%
  summarise(mean = mean(Fhat3))

acf(mean_year_pup$mean)

MannKendall(mean_year_pup$mean)
MannKendall(mean_year_pup$mean)
mkTrend(mean_year_pup$mean)

library(boot)
MKtau<-function(z) MannKendall(z)$tau
tsboot(mean_year_pup$mean, MKtau, R=500, l=5, sim="fixed")


mean_year_adult <- test %>%
  group_by(YEAR) %>%
  filter(Animal == "Adult") %>%
  summarise(mean = mean(Fhat3))

MannKendall(mean_year_adult$mean)
mkTrend(mean_year_adult$mean)















#~~ Other plots

# Melt ibcs

ibcs_melt <- melt(ibcs, id = c("Run", "FID", "IID",
                               "Animal"))

unique(ibcs_melt$Run)

# Subset runs

df <- ibcs %>%
  filter(Run == "GQ_DP_maf05_G_ldi_ME_LD" | Run == "GQ_DP_maf05_G_ldi_ME" | 
           Run == "maf05_g90_ldi_ME_LD" | Run == "maf05_g60")


# violin
ggplot(fhats, aes(factor(Animal), Fhat1)) + 
  geom_violin(aes(fill = factor(Animal))) +
  facet_grid(. ~ .id)

#dens
ggplot(df, aes(x=Fhat3, color=Animal, fill=Animal)) +
  geom_histogram(aes(y=..density..),  position="identity", alpha=0.5) +
  geom_density(alpha=0.6) +
  facet_grid(.~ Run)

ggplot(filter(fhats, .id == "GQ_DP_maf05_g90_ldi_ld.ibc"), aes(x=Fhat3, color=Animal, fill=Animal)) +
  geom_histogram(aes(y=..density..),  position="identity", alpha=0.5) +
  geom_density(alpha=0.6)

pairs(filter(fhats, .id == "GQ_DP_maf05_g90_ldi.ibc")[c(4:7)])

# nomiss
ggplot(fhats, aes(x=NOMISS, y = Fhat1)) +
  geom_point() +
  facet_grid(. ~ .id, scales = "free")


# Pairs plot

df <- ibcs %>%
  filter(Run == "GQ_DP_maf05_G_ldi_ME_LD")
pairs(df[c(5,6,7,9)])

library(GGally)
ggpairs(df[c(5,6,7,9)])






#~~ Relatedness


system("plink --bfile data/inbreeding/relatedness_ld --genome --out data/inbreeding/relatedness_ld --allow-extra-chr")

ibds <- fread("data/inbreeding/relatedness_ld.genome", header = T)


ped <- fread("data/raw/full_pedigree.txt") %>%
  dplyr::filter(grepl("^AG", V1))
ped[ped == 0] <- NA
colnames(ped) <- c("id", "dam", "sire")

#install.packages("pedantics")
library(pedantics)
library(reshape2)

stats.g<-pedigreeStats(ped, graphicalReport='n', includeA = T)

mat <- stats.g$Amatrix


relate <- setNames(melt(mat), c('rows', 'vars', 'values')) %>%
  dplyr::filter(values != 1)
colnames(relate) <- c("IID1", "IID2", "R")

# Add pedigree relatedness to SNP relatedness 

test <- dplyr::left_join(ibds, relate, by = c("IID1", "IID2")) %>%
  filter(grepl("^AG", FID1))


library(RColorBrewer)
cbPalette <- c(brewer.pal(8, "Dark2")[c(1)], brewer.pal(8, "Dark2")[c(2)])


ggthemr(palette = "pale", layout = "clean", 
        line_weight = 0.7, text_size = 20, type = "outer")


png("figs/snp_relatedness.png", units = "in", res = 300, width = 10, height = 8)

ggplot(test, aes(x=R, y = PI_HAT)) +
  geom_point(size = 4, alpha = 0.5, col = "#3262AB") +
  ylab("SNP based relatedness") +
  xlab("Pedigree based relatedness")

dev.off()

# FOR EXTRACTING PARENTAL RELATEDNESS TO CORRELATE WITH INBREEDING

# probably very bad code to get mum dad combos

arranged <- test %>%
  mutate(MOTHER1 = case_when(grepl("AGF", .$IID1) & grepl("AGM", .$IID2) ~ .$IID1),
         FATHER1 = case_when(grepl("AGF", .$IID1) & grepl("AGM", .$IID2) ~ .$IID2),
         MOTHER2 = case_when(grepl("AGF", .$IID2) & grepl("AGM", .$IID1) ~ .$IID2),
         FATHER2 = case_when(grepl("AGF", .$IID2) & grepl("AGM", .$IID1) ~ .$IID1),
         MOTHER = case_when(is.na(MOTHER1) ~ MOTHER2, is.na(MOTHER2) ~ MOTHER1),
         FATHER = case_when(is.na(FATHER1) ~ FATHER2, is.na(FATHER2) ~ FATHER1)) %>%
  dplyr::select(-MOTHER1, -MOTHER2, -FATHER1, -FATHER2) %>%
  filter(!is.na(MOTHER))

# select triads

triads <- read.table("data/raw/new_pedigree.txt", header = T) %>%
  filter(grepl("^AG", ANIMAL))

x <- triads %>%
  left_join(arranged, by = c("MOTHER", "FATHER"))


