library(ohchibi)
library(phytools)
library(RColorBrewer)

set.seed(130816)
setwd('/Users/isaisalasgonzalez/Documents/rootbarriersmicro/scripts/')


Dat <- readRDS(file = "../cleandata/dat_screeningbarriers.RDS")
Map <- Dat$Map
tree <- Dat$tree

significant_strains <- read.table(file = "../rawdata/41_isolates.tsv") %$% V3

Map$Chosen <- "No"
Map$Chosen[which(Map$taxon_oid %in% significant_strains)] <- "Yes"

toremove <- Map %>% subset(Chosen == "No") %$% taxon_oid %>% as.character
tree <- ape::drop.tip(phy = tree,tip = toremove)
write.tree(phy = tree,file = "../cleandata/tree_41syncom.newick")
