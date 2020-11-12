library(ohchibi)
library(DESeq2)

set.seed(130816)
setwd('/home/isai/Documents/results/rootbarriersmicro/scripts/')


Dat <- readRDS(file = "../cleandata/dat_rootbarriers_amplicon_bacteria_syncom_stresses.RDS")
Dat_raw <- Dat$RawCounts


Dat_sub <- Dat_raw %>% 
  subset.Dataset(subset = Fraction != "Blank",drop = T,clean = T) %>%
  subset.Dataset(subset = Bacteria != "NB",drop = T,clean = T)  %>%
  subset.Dataset(subset = Fraction != "Inoculum",drop = T,clean = T) %>%
  subset.Dataset(subset = Plant != "NoPlant",drop = T,clean = T)

Dat_sub$Map$Fraction <- Dat_sub$Map$Fraction %>%
  gsub(pattern = "Agar",replacement = "AgarPlant") %>%
  factor


### Fraction model ###
dds <- DESeqDataSetFromMatrix(countData = Dat_sub$Tab,
                              colData = Dat_sub$Map,
                              design = ~ Experiment + Fraction + Stress )

dds <- DESeq(dds)
dds_main <- dds


#Look for the interaction effect
#Remove NaCl given it is not present in Agar
Dat_sub <- Dat_raw %>% 
  subset.Dataset(subset = Fraction != "Blank",drop = T,clean = T) %>%
  subset.Dataset(subset = Bacteria != "NB",drop = T,clean = T)  %>%
  subset.Dataset(subset = Fraction != "Inoculum",drop = T,clean = T) %>%
  subset.Dataset(subset = Plant != "NoPlant",drop = T,clean = T) 
Dat_sub$Map$Fraction <- Dat_sub$Map$Fraction %>%
  gsub(pattern = "Agar",replacement = "AgarPlant") %>%
  factor

dds <- DESeqDataSetFromMatrix(countData = Dat_sub$Tab,
                              colData = Dat_sub$Map,
                              design = ~ Experiment + Fraction + Stress + Fraction:Stress  )

dds <- DESeq(dds,test = "LRT",reduced = ~Experiment + Fraction + Stress )
plotDispEsts(dds)
dds_int <- dds

##### Test stress effect inside each  fraction
### Root ###
Dat_sub <- Dat_raw %>% 
  subset.Dataset(subset = Fraction == "Root",drop = T,clean = T) %>%
  subset.Dataset(subset = Bacteria != "NB",drop = T,clean = T)  %>%
  subset.Dataset(subset = Fraction != "Inoculum",drop = T,clean = T) %>%
  subset.Dataset(subset = Plant != "NoPlant",drop = T,clean = T) 

dds <- DESeqDataSetFromMatrix(countData = Dat_sub$Tab,
                              colData = Dat_sub$Map,
                              design = ~Experiment +  Stress)

dds <- DESeq(dds)
plotDispEsts(dds)
dds_root <- dds



### Agar ###
Dat_sub <- Dat_raw %>% 
  subset.Dataset(subset = Fraction == "Agar",drop = T,clean = T) %>%
  subset.Dataset(subset = Bacteria != "NB",drop = T,clean = T)  %>%
  subset.Dataset(subset = Fraction != "Inoculum",drop = T,clean = T) %>%
  subset.Dataset(subset = Plant != "NoPlant",drop = T,clean = T) 

dds <- DESeqDataSetFromMatrix(countData = Dat_sub$Tab,
                              colData = Dat_sub$Map,
                              design = ~ Experiment + Stress )

dds <- DESeq(dds)
plotDispEsts(dds)
dds_agar <- dds


### Shoot ###
Dat_sub <- Dat_raw %>% 
  subset.Dataset(subset = Fraction == "Shoot",drop = T,clean = T) %>%
  subset.Dataset(subset = Bacteria != "NB",drop = T,clean = T)  %>%
  subset.Dataset(subset = Fraction != "Inoculum",drop = T,clean = T) %>%
  subset.Dataset(subset = Plant != "NoPlant",drop = T,clean = T) 

dds <- DESeqDataSetFromMatrix(countData = Dat_sub$Tab,
                              colData = Dat_sub$Map,
                              design = ~ Experiment + Stress )

dds <- DESeq(dds)
plotDispEsts(dds)
dds_shoot <- dds
rm(dds)

ResList <- list(MainEffects = dds_main,
                Interaction = dds_int,
                AgarStress = dds_agar,     
                RootStress = dds_root,
                ShootStress = dds_shoot
)

saveRDS(object = ResList,
        file = "../cleandata/deseq2_models_rootbarriers_amplicon_bacteria_syncom_stresses.RDS")
