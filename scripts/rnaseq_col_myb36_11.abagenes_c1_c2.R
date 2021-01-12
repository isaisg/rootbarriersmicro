library(ohchibi)
library(DESeq2)
library(clusterProfiler)
library(org.At.tair.db)
library(scales)
library(egg)
library(ggpubr)
library(ggvenn)

setwd('/home/isai/Documents/results/rootbarriersmicro/scripts/')

set.seed(130816)
Dat <-readRDS(file = "../cleandata/dat_rnaseq_col_myb36.RDS")
Res <- readRDS(file = "../cleandata/res_rnaseq_col_myb36.RDS")

melted_av <- Dat$melted_av
clean_list <- Res$SetsList
df_class <- Res$df_classification

mgenes <- df_class %>% subset(Classification == "bacteria_up_unique" | Classification == "bacteria_down_unique") %$% Gene %>%
  as.character %>% unique

## Read aba  genes ###
Res_aba <- readRDS(file = "../cleandata/aba_robust_genes.RDS")
chosen_up <- Res_aba$aba_genes_up
chosen_down <- Res_aba$aba_genes_down


maba <- union(chosen_up,chosen_down)

chosengenes <- intersect(mgenes,maba)

## Load the description of this genes and see which ons is key component of the aba
mdf <- read.table(file = "../cleandata/mgene_description_20131231.txt",header = F,sep = "\t",fill = NA,quote = "",comment.char = "")
mdf$V1 <- mdf$V1 %>% gsub(pattern = "\\..*",replacement = "")
df <- data.frame(Gene = chosengenes,Description = match(chosengenes,mdf$V1) %>%mdf$V2[.])

df_aba <- read.table(file = "../cleandata/aag1550_Table_S1.txt",skip = 3,header = F,sep = "\t")
colnames(df_aba) <- c("Gene","h1","h4","h8","h12","h24","h36","h60")


rownames(df_aba) <- df_aba$Gene
mat_aba <- df_aba[,-1] %>% as.matrix

df_aba <- rowMeans(mat_aba) %>% as.data.frame

df$Av <- match(df$Gene,rownames(df_aba)) %>%
  df_aba[.,]


df <- with(df,order(-Av)) %>%
  df[.,]

df$Classification <- match(df$Gene,df_class$Gene) %>%
  df_class$Classification[.] %>% droplevels


which(df$Gene %in% c("AT1G16540","AT1G30100","AT1G52340",
                     "AT1G67080","AT1G78390","AT2G27150","AT2G47130","AT3G14440",
                     "AT3G24220","AT4G18350","AT5G67030"))
      
colnames(df)[3] <- "AverageVFoldChangeEligeAbaPaper"

df$Classification <- df$Classification %>% gsub(pattern = "bacteria_down_unique",replacement = "C2") %>%
  gsub(pattern = "bacteria_up_unique",replacement = "C1")

write.table(x = df,file = "../cleandata/df_aba_c1c2.csv",append = F,quote = F,sep = ",",row.names = F,col.names = T)
