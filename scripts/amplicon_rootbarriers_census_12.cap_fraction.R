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

##PCoa ###
mcap <- oh.cap(Tab = Dat_rar$Tab ,Map = Dat_rar$Map,formula = "Fraction + Genotype + Condition(Rep)",
               perms = 9999,distfun = distfun)
p <- chibi.cap(list_ohpco = mcap,col_val = "Fraction",comp_a = "CAP1",comp_b = "CAP2",size = 8,size_panel_border = 2) +
  scale_fill_fraction() + theme(legend.position = "none")
oh.save.pdf(p = p,outname = "cap_fraction_census_microbiome.pdf",outdir = "../figures/",width = 8,height = 8)
Tab_bray <- distfun(t(Dat_rar$Tab))
mypermanova <- adonis(Tab_bray ~ Fraction + Genotype,
                      data = Dat_rar$Map,
                      strata = Dat_rar$Map$Rep,
                      permutations = 9999)
mypermanova
