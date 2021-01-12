library(ohchibi)
library(DESeq2)

set.seed(130816)
setwd('/home/isai/Documents/results/rootbarriersmicro/scripts/')


Dat <- readRDS(file = "../cleandata/dat_rootbarriers_amplicon_bacteria_syncom_mutants.RDS")
Dat_raw <- Dat$RawCounts


Dat_sub <- Dat_raw %>% 
  subset.Dataset(subset = Genotype != "blank",drop = T,clean = T) %>%
  subset.Dataset(subset = Genotype != "noPlant",drop = T,clean = T) %>%
  subset.Dataset(subset = Treatment != "NB",drop = T,clean = T) %>%
  subset.Dataset(subset = Fraction != "Inoculum",drop = T,clean = T) 


### Fraction model ###
dds <- DESeqDataSetFromMatrix(countData = Dat_sub$Tab,
                              colData = Dat_sub$Map,
                              design = ~ Replica + Fraction + Genotype )

dds <- DESeq(dds)
dds_main <- dds


#Look for the interaction effect
#Remove NaCl given it is not present in Agar
Dat_sub <- Dat_raw %>% 
  subset.Dataset(subset = Genotype != "blank",drop = T,clean = T) %>%
  subset.Dataset(subset = Genotype != "noPlant",drop = T,clean = T) %>%
  subset.Dataset(subset = Treatment != "NB",drop = T,clean = T) %>%
  subset.Dataset(subset = Fraction != "Inoculum",drop = T,clean = T) 


dds <- DESeqDataSetFromMatrix(countData = Dat_sub$Tab,
                              colData = Dat_sub$Map,
                              design = ~ Replica + Fraction + Genotype + Fraction:Genotype  )

dds <- DESeq(dds,test = "LRT",reduced = ~Replica + Fraction + Genotype )
plotDispEsts(dds)
dds_int <- dds

##### Test stress effect inside each  fraction
### Root ###
Dat_sub <- Dat_raw %>% 
  subset.Dataset(subset = Genotype != "blank",drop = T,clean = T) %>%
  subset.Dataset(subset = Genotype != "noPlant",drop = T,clean = T) %>%
  subset.Dataset(subset = Treatment != "NB",drop = T,clean = T) %>%
  subset.Dataset(subset = Fraction != "Inoculum",drop = T,clean = T) %>%
  subset.Dataset(subset = Fraction == "Root",drop = T,clean = T)

dds <- DESeqDataSetFromMatrix(countData = Dat_sub$Tab,
                              colData = Dat_sub$Map,
                              design = ~Replica +  Genotype)

dds <- DESeq(dds)
plotDispEsts(dds)
dds_root <- dds



### Agar ###
Dat_sub <- Dat_raw %>% 
  subset.Dataset(subset = Genotype != "blank",drop = T,clean = T) %>%
  subset.Dataset(subset = Genotype != "noPlant",drop = T,clean = T) %>%
  subset.Dataset(subset = Treatment != "NB",drop = T,clean = T) %>%
  subset.Dataset(subset = Fraction != "Inoculum",drop = T,clean = T) %>% 
  subset.Dataset(subset = Fraction == "Agar",drop = T,clean = T) 

dds <- DESeqDataSetFromMatrix(countData = Dat_sub$Tab,
                              colData = Dat_sub$Map,
                              design = ~ Replica + Genotype )

dds <- DESeq(dds)
plotDispEsts(dds)
dds_agar <- dds


### Shoot ###
Dat_sub <- Dat_raw %>% 
  subset.Dataset(subset = Genotype != "blank",drop = T,clean = T) %>%
  subset.Dataset(subset = Genotype != "noPlant",drop = T,clean = T) %>%
  subset.Dataset(subset = Treatment != "NB",drop = T,clean = T) %>%
  subset.Dataset(subset = Fraction != "Inoculum",drop = T,clean = T) %>% 
  subset.Dataset(subset = Fraction == "Shoot",drop = T,clean = T) 

dds <- DESeqDataSetFromMatrix(countData = Dat_sub$Tab,
                              colData = Dat_sub$Map,
                              design = ~ Replica + Genotype )

dds <- DESeq(dds)
plotDispEsts(dds)
dds_shoot <- dds
rm(dds)

ResList <- list(MainEffects = dds_main,
                Interaction = dds_int,
                AgarMutants = dds_agar,     
                RootMutants = dds_root,
                ShootMutants = dds_shoot
)

saveRDS(object = ResList,
        file = "../cleandata/deseq2_models_rootbarriers_amplicon_bacteria_syncom_mutants.RDS")

rm(list=ls())
dev.off()
