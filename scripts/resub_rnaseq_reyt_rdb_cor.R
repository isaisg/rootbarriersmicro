library(ohchibi)
library(DESeq2)
library(clusterProfiler)
library(org.At.tair.db)
library(scales)
library(egg)
library(ggpubr)

setwd('/Users/isaisalasgonzalez/Documents/rootbarriersmicro/scripts/')

set.seed(130816)

Dat <- readRDS("../cleandata/design1_dat_deseq2.RDS")
Dat_z <- Dat$Dat_z %>%
  subset.Dataset(Genotype == "WT" | Genotype == "myb36",drop = T,clean = T)

melted_av <- Dat_z$Tab %>% melt
colnames(melted_av) <- c("Gene","Name","value")
melted_av <- merge(melted_av,Dat_z$Map, by = "Name")

Tab_sub <- dcast(data = melted_av,formula = Gene~Genotype,value.var = "value",
                 fun.aggregate = mean)
melted_av <- Tab_sub %>% melt
colnames(melted_av) <- c("Gene","Genotype","value")
melted_av_reyt <- melted_av


#Load the expression for rootbariers data
melted_av_rb <-readRDS(file = "../cleandata/dat_rnaseq_col_myb36.RDS") %$% melted_av


melted_av_reyt$Dataset <- "Reyt"
melted_av_rb$Dataset <- "RDB"
melted_av_rb <- melted_av_rb %>% subset(Bacteria == "NB") %>% droplevels
melted_av_rb <- melted_av_rb[,c("Gene","Genotype","value","Dataset")]


melted_av_reyt$Genotype <- melted_av_reyt$Genotype %>% gsub(pattern = "WT",replacement = "Col_0")


melted_av <- rbind(melted_av_reyt,melted_av_rb)

df <- dcast(data = melted_av,formula = Gene~Genotype + Dataset,value.var = "value")

toremove <- c(which(is.na(df$Col_0_RDB)),which(is.na(df$Col_0_Reyt))) %>% unique

toremove <- c(which(is.na(df$myb36_RDB)),which(is.na(df$myb36_Reyt))) %>% unique

df[-toremove,] %>%
  ggplot(data = .,mapping = aes(myb36_Reyt,myb36_RDB)) +
  stat_cor() +
  geom_point()



#What about fold changes
res_reyt <- Dat$Res_DESeq2 %>% subset(Genotype == "myb36")
dds <-readRDS(file = "../cleandata/dat_rnaseq_col_myb36.RDS") %$% dds
res_rdb <- results(dds,contrast = c("group","myb36_NB","Col_0_NB")) %>% as.data.frame
res_rdb$Gene <- rownames(res_rdb)

merged <- merge(res_reyt[,c("Gene","log2FoldChange")],res_rdb[c("Gene","log2FoldChange")], by = "Gene") 
ggplot(data = merged,aes(log2FoldChange.y,log2FoldChange.x)) + 
  stat_cor() +
  geom_point()


a <- res_reyt %>% subset(padj < 0.1) %$% Gene
b <- res_rdb %>% subset(padj < 0.1) %$% Gene

c <- c(a,b) %>% unique

which(merged$Gene %in% c) %>% merged[.,] %>%
  droplevels %>%
  ggplot(data = .,aes(log2FoldChange.y,log2FoldChange.x)) + 
  stat_cor() +
  geom_point()
  
