library(ohchibi)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(harrietr)
library(egg)
library(paletteer)
library(ggtree)
library(scales)
library(car)
library(Rmisc)
library(ggpubr)

set.seed(130816)
setwd('/home/isai/Documents/results/rootbarriersmicro/scripts/')


Dat <- readRDS("../cleandata/dat_41isolatesmonoassociation_cfu_suberin_primroot_dryweight_ionome_zscoreRaw.distance.RDS")
df <- Dat$Merged$Raw

#select 41 tested
st <- read.table(file = "../rawdata/41_isolates.tsv") %$% V2 %>% as.character

df <- which(df$Strains %in% c(st,"NB")) %>% df[.,] %>% droplevels
colnames(df)

df <- df[,c(1,12,13,24:57)]


#add the group information
df_groups <- read.table(file = "../cleandata/df_41_suberin_results.distance.tsv",header = T,sep = "\t")
merged <- merge(df,df_groups, by = "Strains",all.x = T)
merged$Group <- merged$Group %>% as.character
merged$Group[which(is.na(merged$Group))] <- "HK"
merged$Group[merged$Strains %>% grep(pattern = "NB")] <- "NB"

merged$Group <- merged$Group %>% 
  factor(levels = c("NB","HK","Group1","Group2","Group3","Group4","Group5"))



paleta_alive <- c("#414141","#D9D9D9",paletteer_d(package = "LaCroixColoR",palette = "PassionFruit",n = 5))
names(paleta_alive) <-  c("NB","HK","Group1","Group2","Group3","Group4","Group5")

df <- merged

mions <- colnames(df) %>% grep(pattern = "Strains|Dist|Group",invert = T,value = T) %>%
  grep(pattern = "Mean",value = T) %>% gsub(pattern = "Mean",replacement = "")

mmodels <- list()
Res <- NULL
for(ion in mions){
  mm <- paste0("Mean",ion)
  mc <- paste0("CI",ion)
  targets <- c("Strains","MeanDistSub","CIDistSub",mm,mc,"Group")
  df_temp <- which(colnames(df) %in% targets) %>%
    df[,.]
  df_temp$Down <- df_temp[,4]-df_temp[,5]
  df_temp$Up<- df_temp[,4]+df_temp[,5]
  m1 <- lm(df_temp[,4]~ MeanDistSub,df_temp)
  mmodels[[ion]] <- m1
}

saveRDS(object = mmodels,file = "../cleandata/models_lm_41.distance.RDS")
