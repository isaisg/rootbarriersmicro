library(ohchibi)
library(DESeq2)
library(clusterProfiler)
library(org.At.tair.db)
library(scales)
library(egg)
library(ggpubr)

setwd('/Users/isaisalasgonzalez/Documents/rootbarriersmicro/scripts/')

## Load the number scheme
df_num <- read.table(file = "../rawdata/df_genotypes_to_classification.csv",header = T,sep = ",")
df_num$Number <- 1:20


set.seed(130816)

Dat <- readRDS("../cleandata/design1_dat_deseq2.RDS")
Dat_z <- Dat$Dat_z %>%
  subset.Dataset(Genotype != "ucc1") %>%
  subset.Dataset(Genotype != "Wtcif") %>%
  subset.Dataset(Genotype != "sgn3cif") %>%
  subset.Dataset(Genotype != "ucc1.38ucc2GK",drop = T,clean = T)

melted_av <- Dat_z$Tab %>% melt
colnames(melted_av) <- c("Gene","Name","value")
melted_av <- merge(melted_av,Dat_z$Map, by = "Name")

Tab_sub <- dcast(data = melted_av,formula = Gene~Genotype,value.var = "value",
      fun.aggregate = mean)
melted_av <- Tab_sub %>% melt
colnames(melted_av) <- c("Gene","Genotype","value")


res <- Dat$Res_DESeq2
melted_av <- merge(melted_av,
                   res[,c("log2FoldChange","Gene","Genotype","padjG")], 
                   by = c("Gene","Genotype"),all.x = TRUE)
melted_av$Significance <- NA
melted_av$Significance[which(melted_av$padjG < 0.1)] <- "q < 0.1"

## Read the ABA genes ###
Res_aba <- readRDS(file = "../cleandata/aba_robust_genes.RDS")
chosen_up <- Res_aba$aba_genes_up
chosen_down <- Res_aba$aba_genes_down

#Adjust the names so they match the numbers
melted_av$Genotype <- melted_av$Genotype %>% 
  gsub(pattern = "^WT$",replacement = "Col_0") %>%
  gsub(pattern = "^esb1$",replacement = "esb1_1") %>%
  gsub(pattern = "^esb1sgn3$",replacement = "esb1_1sng3_3") %>%
  gsub(pattern = "^myb36$",replacement = "myb36_2") %>%
  gsub(pattern = "^sgn3$",replacement = "sgn3_3") %>%
  gsub(pattern = "^sm$",replacement = "myb36_2_sgn3_3")

melted_av$Number <- match(melted_av$Genotype,df_num$Genotype) %>%
  df_num$Number[.] %>% factor


####### ABA ##########
#Upregulated
melted_sub <- which(melted_av$Gene %in% chosen_up ) %>%
  melted_av[.,] %>% droplevels

lims_col <- melted_sub %>%
  subset(Genotype == "WT") %$% value %>%
  quantile



p1 <- ggplot(data = melted_sub,aes(Number,value)) +
  geom_rect(mapping = aes(ymin = lims_col[2],ymax = lims_col[4],xmin = -Inf,xmax =Inf),
            color = NA,fill = "#A6CEE3",alpha = 0.1) +   geom_jitter(aes(color = Significance),alpha = 0.6,size = 2,stroke = 0)+
  scale_color_manual(values = c("red"),na.value = "grey",name = "Significance vs Col-0")+
  geom_boxplot(outlier.size = 0,outlier.stroke = 0,fill = NA)  +
  theme_ohchibi(font_family = "Helvetica") +
  theme(
    
  ) +
  ylab(label = "Standardised gene expression") +
  scale_y_continuous(limits = c(-1,2))


reyt_up_genes <- melted_sub %>%
  subset((Genotype == "esb1" | Genotype == "myb36") & Significance == "q < 0.1") %$%
  Gene %>% as.character %>% unique
saveRDS(object = reyt_up_genes,file = "../cleandata/reyt_up_genes.RDS")

#Downregulated
  melted_sub <- which(melted_av$Gene %in% chosen_down ) %>%
    melted_av[.,] %>% droplevels
  
  lims_col <- melted_sub %>%
    subset(Genotype == "WT") %$% value %>%
    quantile
  
  
 p2 <-  ggplot(data = melted_sub,aes(Number,value)) +
    geom_rect(mapping = aes(ymin = lims_col[2],ymax = lims_col[4],xmin = -Inf,xmax =Inf),
              color = NA,fill = "#A6CEE3",alpha = 0.1) +   geom_jitter(aes(color = Significance),alpha = 0.6,size = 2,stroke = 0)+
    scale_color_manual(values = c("red"),na.value = "grey",name = "Significance vs Col-0")+
    geom_boxplot(outlier.size = 0,outlier.stroke = 0,fill = NA)  +
    theme_ohchibi(font_family = "Helvetica") +
    theme(
      
    ) +
    ylab(label = "Standardised gene expression") +
   scale_y_continuous(limits = c(-1,2))
 

#### Ethylene ###
mregs <- readRDS(file = "../cleandata/regulons_rnaseq.RDS")
mgenes <- mregs$ethylene
 
melted_sub <- which(melted_av$Gene %in% mgenes ) %>%
  melted_av[.,] %>% droplevels

lims_col <- melted_sub %>%
  subset(Genotype == "WT") %$% value %>%
  quantile


p3 <- ggplot(data = melted_sub,aes(Number,value)) +
  geom_rect(mapping = aes(ymin = lims_col[2],ymax = lims_col[4],xmin = -Inf,xmax =Inf),
            color = NA,fill = "#A6CEE3",alpha = 0.1) +   geom_jitter(aes(color = Significance),alpha = 0.6,size = 2,stroke = 0)+
  scale_color_manual(values = c("red"),na.value = "grey",name = "Significance vs Col-0")+
  geom_boxplot(outlier.size = 0,outlier.stroke = 0,fill = NA)  +
  theme_ohchibi(font_family = "Helvetica") +
  theme(
    
  ) +
  ylab(label = "Standardised gene expression") +
  scale_y_continuous(limits = c(-1,2))

#### Check the patterns
Res <- readRDS(file = "../cleandata/res_rnaseq_col_myb36.RDS")
df_classification <- Res$df_classification

match(reyt_up_genes,df_classification$Gene) %>%
  df_classification[.,] %$% Classification %>% table

composition <- egg::ggarrange(p1 + theme(legend.position = "none"),p2 + theme(legend.position = "none"),p3 + theme(legend.position = "none"),nrow = 1)

oh.save.pdf(p = composition,outname = "resub_composition_reyt_hormones.pdf",outdir = "../figures/",width = 18,height = 4)
rm(list=ls())
dev.off()
gc()
