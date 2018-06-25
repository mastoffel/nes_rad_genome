# inbreeding analysis

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
library(vcfR)
library(purrr)
source("martin.R")
# system("mkdir data/inbreeding")


# Count number of SNPs in vcf file
system("grep -v '#' data/nes_rad.vcf | wc -l")

#~~ Filter raw vcf file
#~~~~~~~~~~~~~~~~~~~~~~#
# 5% geno  # maf 0.01 # minDP 10
#~~~~~~~~~~~~~~~~~~~~~~#  
#  #--min-meanDP 10  --maf 0.01 --max-missing-count --max-missing 0.05
system("/home/martin/bin/vcftools --vcf data/nes_rad.vcf --out data/inbreeding/nes_filtered --remove-indels --min-alleles 2 --max-alleles 2 --minDP 7 --maxDP 20 --min-meanDP 10 --max-meanDP 20 --max-missing 0.50  --maf 0.01  --recode --recode-INFO-all")
system("grep -v '#' data/inbreeding/nes_filtered.recode.vcf | wc -l")

# # filter 1
# system("/home/martin/bin/vcftools --vcf data/nes_rad.vcf --out data/inbreeding/nes_filtered1 --remove-indels --minGQ 20 --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all")
# system("grep -v '#' data/inbreeding/nes_filtered1.recode.vcf | wc -l")
# 
# # filter 2
# system("/home/martin/bin/vcftools --vcf data/inbreeding/nes_filtered1.recode.vcf --out data/inbreeding/nes_filtered2 --minDP 10 --maxDP 20 --recode --recode-INFO-all")
# 
# # filter 3
# system("/home/martin/bin/vcftools --vcf data/inbreeding/nes_filtered2.recode.vcf --out data/inbreeding/nes_filtered3  --max-missing-count 93 --recode --recode-INFO-all")
# 
# # filter 4
# system("/home/martin/bin/vcftools --vcf data/inbreeding/nes_filtered3.recode.vcf --out data/inbreeding/nes_filtered4 --min-meanDP 10 --max-meanDP 20 --maf 0.01 --recode --recode-INFO-all")

#--min-meanDP 10 --max-meanDP 20 --minDP 10 --maxDP 20 
#       --max-missing-count 90 --remove-indels --maf 0.01 --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all")


# 
# load vcf_file
#vcf_file <- "data/inbreeding/nes_filtered.recode.vcf"
# read vcf
#nes_vcf <- read.vcfR(vcf_file, verbose = FALSE )

# create plink files from vcf
system("/home/martin/bin/vcftools --vcf data/inbreeding/nes_filtered.recode.vcf --plink --out data/inbreeding/nes_filtered")

# add names to chromosomes in map file for plink to work
nes_map <- read_delim("data/inbreeding/nes_filtered.map", delim = "\t", col_names = FALSE) %>% 
  mutate(X1 = paste0("chr", X1)) %>% 
  mutate(X2 = paste0("chr", X2))
write_delim(nes_map, path = "data/inbreeding/nes_filtered.map", delim = "\t", col_names = FALSE)

# create plink raw file
system("/home/martin/bin/plink --file data/inbreeding/nes_filtered --make-bed --allow-extra-chr --recodeAD --out data/inbreeding/nes_filtered")


get_sMLH_from_plinkraw <- function(file) {
  
  x <- fread(file, colClasses = "character")
  ids <- x$IID
  x <- dplyr::select(x, contains("HET"))
  row.names(x) <- ids
  NAs <- apply(x, 1, function(x) sum(is.na(x)))
  
  sMLH <- as.data.frame(sMLH(x))
  sMLH$ANIMAL <- ids
  sMLH$NAS <- NAs
  colnames(sMLH) <- c("sMLH", "ANIMAL", "NAs")
  sMLH 
  
}

raw_files <- paste("data/inbreeding/", list.files(path = "data/inbreeding", pattern="*.raw"), sep = "")
# sMLH dataframe
nes_sMLH <- get_sMLH_from_plinkraw(raw_files)

plot(nes_sMLH$NAs,  nes_sMLH$sMLH)
summary(lm(nes_sMLH$NAs ~  nes_sMLH$sMLH))
#~~ g2 (boot over loci)

get_g2_from_plinkraw <- function(file) {
  
  x <- fread(file, colClasses = "numeric")
  
  ids <- x$IID
  x <- dplyr::select(x, contains("HET"))
  row.names(x) <- ids
  
  g2 <- g2_snps(x, nboot = 100, nperm = 100, CI = 0.95)
  g2
  
}
# load(file = "data/inbreeding/nes_g2.rda")
set.seed(102)
nes_g2 <- get_g2_from_plinkraw(raw_files)
# save(nes_g2, file = "data/inbreeding/nes_g2.rda") # 7407 snps
plot(nes_g2, col = "grey")
nes_g2
# var(filter(sMLH, Run == "GQ_DP_maf05_G_ldi_ME_SSB_hwe_LD")$sMLH)

get_r2hf_from_plinkraw <- function(file) {
  
  x <- fread(file, colClasses = "numeric")
  
  ids <- x$IID
  x <- dplyr::select(x, contains("HET"))
  row.names(x) <- ids
  
  r2 <- r2_hf(x, nboot = 100, CI = 0.95, type = "snps")
  r2
  
}
nes_r2hf <- get_r2hf_from_plinkraw(raw_files)
plot(nes_r2hf)
#~~ Get Fhats

# recode bim files for GCTA

recode_bim_1chr <- function(file){
  file <- fread(file) %>%
    mutate(V1 = 1) %>%
    fwrite(file, quote = F, row.names = F,
           col.names = F, sep = " ")
}

bim_files <- paste("data/inbreeding/", list.files(path = "data/inbreeding", pattern="*.bim"), sep = "")
recode_bim_1chr(bim_files)


# get fhats using gcta
plink_files <- paste("data/inbreeding/", list.files(path = "data/inbreeding", pattern="*.ped"), sep = "")
plink_files <- lapply(plink_files, function(x) gsub(".ped", "", x))

for (i in 1:length(plink_files)){
  system(paste0("/home/martin/bin/gcta64 --bfile ", plink_files[i]," --autosome --ibc --out ", plink_files[i]," --thread-num 10"))
}

# load fhats

load_fhats <- function(file) {
  fhats <- fread(file, header = T)
}

ibc_file <- paste("data/inbreeding/", list.files(path = "data/inbreeding", pattern="*.ibc"), sep = "")
fhats <- load_fhats(ibc_file)
names(fhats)[1] <- "ANIMAL"

ibcs <- fhats %>% left_join(nes_sMLH, by = c("ANIMAL")) %>% select(-IID)

plot(ibcs$NAs, ibcs$Fhat3)
summary(lm(NAs~Fhat3, data = ibcs))
# ~~ Plot g2

library(ggplot2)

g2_plot <- data.frame(nes_g2$g2_boot)
lcl <- nes_g2$CI_boot[1]
ucl <- nes_g2$CI_boot[2]
g2_boot_summary <- data.frame(lcl, ucl)

ibcs_vars_summary <- data.frame(c(nes_g2$g2,
                                  var(ibcs$sMLH),
                                  var(ibcs$Fhat1),
                                  var(ibcs$Fhat2),
                                  var(ibcs$Fhat3)),
                                c("g2", "sMLH", "Fhat1", "Fhat2", "Fhat3"))

colnames(ibcs_vars_summary) <- c("val", "var")

plot(nes_g2)
# save(nes_g2, file = "data/inbreeding/nes_g2.rda")

hist(ibcs$Fhat3)
hist(ibcs$sMLH)
# g2 bootstrapping distribution showing empirical g2 with CIs

require(gridExtra)
library(sitools)
cbPalette <- c( "#1B9E77", "#66A61E", "#E6AB02", "black", "#7570B3", "#D95F02", "#E7298A")
library(forcats)
png("figs/g2_boot_comp.png", units = "in", res = 300, width = 6, height = 3)
# remove F2 and F1
ibcs_vars_summary <- ibcs_vars_summary %>% 
  filter(!(var %in% c("Fhat1", "Fhat2", "Fhat3")))

g2s_df <- data.frame("val" = c( 0.005737434, 0.0014, 0.0011, 0.004, 0.0002), 
                     "var" = fct_inorder(c("Northern elephant seal", "Soay sheep", "Red deer", "Antarctic fur seal", "Blue tit")))


cbPalette <- c("black", "#7570B3", "#E6AB02", "#66A61E", "cornflowerblue")
g2_CI_plot <-
  ggplot(g2_plot, aes(nes_g2$g2_boot)) +
  geom_histogram(colour = "grey45", fill = "grey45") +
  geom_errorbarh(aes(xmin = g2_boot_summary$lcl , xmax = g2_boot_summary$ucl , y = 100),
                 size = 0.5, color = "black", linetype = "solid", height = 0) +
  geom_vline(data = g2s_df, aes(xintercept = val, colour = var), size = 0.8,
             linetype = c("dashed", rep("solid", 4)), show.legend = T) +
  # without lines
  #geom_vline(data = g2s_df[1, ], aes(xintercept = val, colour = var), size = 0.8,
  #           linetype = c("dashed"), show.legend = T) +
  scale_colour_manual(values = cbPalette, name = " ") +
  #                     breaks = c("g2", "sMLH"),
  #                     labels = c(expression(italic(g[2])),
  #                                #expression("var"(italic(hat(F)["III"]))),
  #                                expression("var"("sMLH")))) +
  
  labs(y = "Counts", x = expression(italic(g[2]))) +
  theme_martin() +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain"),
        axis.title.x = element_text(face = "plain"),
        axis.title.y = element_text(face = "plain"),
        plot.title=element_text(hjust=0, size = 18, face = "plain"))  +
  geom_vline(xintercept = 0, color = "darkgrey", size = 0.5)

g2_CI_plot
dev.off()


# get individual data
library(readxl)
library(stringr)
library(skimr)
# for translating sequencing barcodes to real ids
sample_ids <- read_excel("data/sample_ids.xlsx") %>% 
  mutate(id_seq = paste(library, barcode, sep = "_"))

# read in data from marine mammal center
nes_mmc <- read_excel("data/nes_rad_data.xlsx")
names(nes_mmc)[1] <- "sample_id"

### todo
# prepare data
nes_ibcs <- ibcs %>% select(ANIMAL, NOMISS, Fhat3, sMLH, NAs) %>% 
  rename(id_seq = ANIMAL) %>% 
  # mutate(id_seq = stringr::str_sub(ANIMAL, str_length("merged_sample_") + 1, end = -1)) %>% 
  dplyr::select(id_seq, NOMISS, Fhat3, sMLH, NAs)
# join and sort
nes <- left_join(nes_ibcs, sample_ids, by = "id_seq") %>% 
  select(sample_id, sMLH, Fhat3, NOMISS, NAs, id_seq)

# join with necropsy data
nes_full <- left_join(nes, nes_mmc, "sample_id") %>% 
  mutate(admit_weight = as.numeric(`admit weight`)) %>% 
  mutate(congenital_defect = as.factor(`Congenital defect`)) %>% 
  mutate(malnutrition = as.factor(Malnutrition))
#skim(nes_full)

table(nes_full$group)
table(nes_full$group_broad)

# some check
str(nes_full)

point_size <- 4
ggplot(data = nes_full, aes(x = group_broad, y = sMLH)) +
  geom_boxplot(alpha = 1, col = "darkgrey",  size = 0.5, width = 0.7,  outlier.shape = NA) + #aes(fill = BreedingType),
  geom_jitter(size = point_size, alpha = 0.6, shape = 21, col = "black", width = 0.2, fill = "lightgrey") +
  theme_martin(base_family = "Lato", highlight_family = "Lato")

write_delim(nes_full, "data/inbreeding/all_inbr_coef.txt", delim = " ")

names(nes_full)

ggplot(data = nes_full, aes(x = group_broad, y = sMLH)) +
  geom_point() +
  geom_boxplot()

ggplot(data = nes_full, aes(x = Fhat3, y = sMLH)) +
  geom_point() 


ggplot(data = nes_full, aes(x = Fhat3)) +
  geom_histogram() 

# seals <- seals %>% dplyr::mutate(groups = seal_data$group) %>% 
#   dplyr::mutate(groups_sum = if_else(groups %in% c("control1", "control2", "control3", "control4"), "control",
#                                      if_else(groups == "bacteria", "bacteria", if_else(groups == "worms", "worms", "")))) %>% 
#   dplyr::mutate(blubber = seal_data$`Relative Blubber depth`) %>% 
#   dplyr::mutate(weight = as.numeric(seal_data$`death weight`)) %>% 
#   dplyr::mutate(msat_het = as.numeric(seal_data$HL)) 








