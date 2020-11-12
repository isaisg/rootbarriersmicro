library(ohchibi)
library(palettesPM)
library(Rmisc)
library(ggpubr)
library(DESeq2)
library(dendextend)

set.seed(130816)
setwd('/home/isai/Documents/results/rootbarriersmicro/scripts/')
dir.create("../cleandata")
dir.create("../figures")


Dat <- readRDS(file = "../cleandata/dat_censusmutants_16s.RDS")


## Determine usable genotypes per fraction and rep
Dat_rar <- Dat$Rarefied

#Define Bray Curtis as the dissimilarity
distfun <- function(x,method) vegan::vegdist(x = x, method = "bray")

Dat_sub <- Dat_rar %>% subset.Dataset(Fraction == "Root",drop = T,clean = T)

Tab_bray <- distfun(t(Dat_sub$Tab))
mypermanova <- adonis(Tab_bray ~   Genotype,
                      data = Dat_sub$Map,
                      strata = Dat_sub$Map$Rep,
                      permutations = 9999)


Dat_sub <- Dat_rar %>% subset.Dataset(Fraction == "Shoot",drop = T,clean = T)

Tab_bray <- distfun(t(Dat_sub$Tab))
mypermanova <- adonis(Tab_bray ~   Genotype,
                      data = Dat_sub$Map,
                      strata = Dat_sub$Map$Rep,
                      permutations = 9999)
mypermanova


Dat_sub <- Dat_rar %>% subset.Dataset(Fraction == "Soil",drop = T,clean = T)

Tab_bray <- distfun(t(Dat_sub$Tab))
mypermanova <- adonis(Tab_bray ~   Genotype,
                      data = Dat_sub$Map,
                      strata = Dat_sub$Map$Rep,
                      permutations = 9999)
mypermanova
