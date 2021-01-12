library(ohchibi)
library(egg)
library(emmeans)

set.seed(130816)
setwd('/home/isai/Documents/results/rootbarriersmicro/scripts/')


###### Suberin ##########
############# Suberin ####################
df_suberin <- read.table(file = "../rawdata/dat_syncom_mutants_suberin.csv",header = T,sep = "\t")
colnames(df_suberin)[7] <- "Batch"


toremove <- c("pCASP1::abi1","pCASP1::abi2","pCASP1::MYB41","pELTP::abi1")

df_suberin <- which(!(df_suberin$Genotype %in% toremove)) %>%
  df_suberin[.,] %>% droplevels
df_suberin$Genotype <- df_suberin$Genotype %>% gsub(pattern = " ",replacement = "")

df <- df_suberin %>% subset(Batch == "Batch5" | Batch == "Batch6")

m1 <- aov(formula = Suberin~Genotype+Bacteria+Genotype:Bacteria+Batch,data = df) 
res <- emmeans(m1,pairwise~Genotype:Bacteria,adjust = "none")
dfc <- res$contrasts %>% as.data.frame  
dfc <- dfc$contrast %>% grep(pattern = "NB.*SC") %>%
  dfc[.,]
dfc[c(1,7,13,19,25),]
