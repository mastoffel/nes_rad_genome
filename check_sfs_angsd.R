# min ind 77
sfs <- scan("../angsd_analysis/SFS14/nes.sfs")

# proportions
# polymorphic
sum(sfs[-1]) / sum(sfs) # 0.0022
# singletons
sum(sfs[2]) / sum(sfs) # 0.0016
# doubletons
sum(sfs[3]) / sum(sfs) #0.00005
# singletons among polymorphic
sum(sfs[2]) / sum(sfs[-1]) # 0.72
# doubletons among polymorphic
sum(sfs[3]) / sum(sfs[-1]) # 0.026



sum(sfs[-c(1,2)]) / sum(sfs)
sum(sfs[-c(1,2)]) /sum(sfs[-1])

sfs_prop <- sfs/sum(sfs)
barplot(sfs_prop[-c(1)])
barplot(sfs_prop[-c(1,2)])

sum(sfs[-1])


sfs_sim <- readLines("fastsimcoal_analyses/bootstrap_test/nes/nes_1/nes_MAFpop0.obs")
sfs_sim <- as.numeric(unlist(str_split(sfs_sim[3], pattern = "\\\t")))[1:97]
sfs2 <- sfs_sim
plot(sfs[-c(1,2,3,4)], sfs_sim[-c(1,2,3,4)])
cor(sfs[-c(1,2,3)], sfs_sim[-c(1,2,3)])

barplot(sqrt(sfs[-1]))
barplot(sqrt(sfs2[-1]))

barplot(sfs[-c(1,2)])
barplot(sfs2[-c(1,2)])
barplot(sfs3[-c(1,2)])

plot(sfs[-c(1,2,3)], sfs2[-c(1,2,3)])

test <- unlist(sfs_sim[2, ])
sfs2 <- scan("../angsd_analysis/nes39.sfs")
plot(sfs[-1], sfs2[-1])
plot(sfs[-c(1,2)], sfs2[-c(1,2)])
plot(sfs[-c(1,2,3)], sfs2[-c(1,2,3)])
cor(sfs[-1], sfs2[-1])
cor(sfs[-c(1,2)], sfs2[-c(1,2)])
cor(sfs[-c(1,2,3)], sfs2[-c(1,2,3)])

barplot(sfs[-c(1,2)])
barplot(sfs2[-c(1,2)])
# barplot(norm(sfs2[-c(1,2,3)]))
barplot(sfs[-c(1,2,3)])
barplot(sfs2[-c(1,2,3)])

norm <- function(x) x/sum(x)
sfs <- norm(sfs[-1])
barplot(sfs)

sfs <- norm(sfs[-c(1,length(sfs))])
sfs <- norm(sfs)
barplot(sfs, xlab="Number Chromosomes with Derived Allele", names=1:length(sfs),
    ylab="Proportion of Variable Sites", main=" Pop2 SFS",col='blue')
barplot(sfs[-1])
barplot(sfs[-c(1,2,3)])
hist(sfs[-c(1,2,3)], 10)

sfs
options(scipen=999)
cor(sfs2, sfs)
cor(sfs2[-c(1,2)], sfs[-c(1,2)])

sum(sfs3)
plot(sfs, sfs2, xlim = c(0,170), ylim = c(0,400))
cor(sfs3, sfsold)

plot(sfs3, sfsold)

barplot(sfs2)

barplot(sfs2)
barplot(sfs[-c(1)])
sum(sfs)


sfs3 <- scan("../../demography/angsd/nes25.sfs")
# create names
sfs_names <- sapply(1:length(sfs7), function(x) paste0("d0_", x))

sink("../../demography/fastsimcoal/nes_MAFpop0.obs")
cat("1 observations")
cat("\n")
cat(sfs_names)
cat("\n")
cat(sfs7)
sink()

# fastsimcoal parameter estimates 
# min ind 48
sfs8 <- scan("../../demography/angsd/nesSFS8.sfs")
barplot(sfs[-1])
sum(sfs8)

cor(sfs7, sfs8)

# min ind 86
sfs9 <- scan("../../demography/angsd/nesSFS9.sfs")
barplot(sfs9[-1])
sum(sfs9)

cor(sfs8, sfs9)
# min ind 77 geno_min 5
sfs5 <- scan("../../demography/angsd/nesSFS5.sfs")
barplot(sfs5[-1])
sum(sfs5)

plot(sfs5, sfs7)
