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

######## Heatmap ##########
mgenes <- clean_list %>% unlist %>% as.character %>% unique
melted_sub <- which(melted_av$Gene %in% mgenes) %>% melted_av[.,] %>% droplevels
Tab <- acast(melted_sub,Gene ~group,value.var  ="value")
mdist_genes <- as.dist(1-cor(Tab %>% t))
mclust_genes <- hclust(d = mdist_genes,method ="ward.D2")
mclust_genes %>% plot


order_genes <- mclust_genes$order %>% mclust_genes$labels[.]


mdist_samples <- dist(Tab %>% t)
mclust_samples <- hclust(d = mdist_samples,method ="ward.D2")
mclust_samples %>% plot

order_samples <- mclust_samples$order %>% mclust_samples$labels[.]


melted_sub <- merge(melted_sub,df_class, by = "Gene") 


melted_sub$Gene <- melted_sub$Gene %>% factor(levels = order_genes)
melted_sub$group <- melted_sub$group %>% factor(levels = order_samples)

melted_sub$value %>% sort %>% plot
melted_sub$group <- melted_sub$group %>% factor(levels = c("Col_0_NB","Col_0_SC","myb36_NB","myb36_SC"))

p <- ggplot(data = melted_sub,mapping = aes(group,Gene)) + 
  geom_tile(aes(fill = value),color = "#00000000") + 
  facet_grid(Classification~.,scales = "free",space = "free") +
  theme_ohchibi() + 
  scale_x_discrete(expand = c(0,0)) +
  scale_fill_paletteer_c(package = "pals"
                         ,palette = "kovesi.diverging_bwr_55_98_c37",
                         limits = c(-0.5,0.5),oob = squish) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.background.y = element_blank(),
    strip.text.y = element_text(size = 15,family = "Arial",face = "bold",angle = 90)
  ) 

#add a bar with gene placement 

Res_aba <- readRDS(file = "../cleandata/aba_robust_genes.RDS")
chosen_up <- Res_aba$aba_genes_up
chosen_down <- Res_aba$aba_genes_down


melted_sub$ABAUP <- NA
melted_sub$ABAUP[which(melted_sub$Gene %in% chosen_up)] <- "Selected"
melted_sub$ABADOWN <- NA
melted_sub$ABADOWN[which(melted_sub$Gene %in% chosen_down)] <- "Selected"



melted_sub$Bar <- "Bar"

p2 <- ggplot(data = melted_sub,aes(Bar,Gene))+
  geom_tile(aes(fill = ABAUP,color = ABAUP)) +
  facet_grid(Classification~.,scales = "free",space = "free") +
  theme_ohchibi() + 
  scale_x_discrete(expand = c(0,0)) +
  scale_fill_manual(values = c("black"),na.value = "#00000000")+
  scale_color_manual(values = c("black"),na.value = "#00000000")+
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.background.y = element_blank(),
    strip.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  )  +
  xlab(label = "ABA upregulated")

p3 <- ggplot(data = melted_sub,aes(Bar,Gene))+
  geom_tile(aes(fill = ABADOWN, color =ABADOWN )) +
  facet_grid(Classification~.,scales = "free",space = "free") +
  theme_ohchibi() + 
  scale_x_discrete(expand = c(0,0)) +
  scale_fill_manual(values = c("black"),na.value = "#00000000")+
  scale_color_manual(values = c("black"),na.value = "#00000000")+
  
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.background.y = element_blank(),
    strip.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  )  +
  xlab(label = "ABA downregulated")

## add bar with ethylene genes
mregs <- readRDS(file = "../../variovorax/rnaseq/cleandata/regulons_rnaseq.RDS")
mgenes <- mregs$ethylene
melted_sub$ETH <- NA
melted_sub$ETH[which(melted_sub$Gene %in% mgenes)] <- "Selected"


p4 <- ggplot(data = melted_sub,aes(Bar,Gene))+
  geom_tile(aes(fill = ETH, color =ETH )) +
  facet_grid(Classification~.,scales = "free",space = "free") +
  theme_ohchibi() + 
  scale_x_discrete(expand = c(0,0)) +
  scale_fill_manual(values = c("black"),na.value = "#00000000")+
  scale_color_manual(values = c("black"),na.value = "#00000000")+
  
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.background.y = element_blank(),
    strip.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  )  +
  xlab(label = "Ethylene")
#### add a bar with embo genes ######

melted_sub$EMBO <- NA
df_embo <- read.table(file = "../rawdata/embo_930genes.csv",header = T,sep = "\t")
#mgenes <- df_embo %>% subset(cluster2 == 1) %$% rowname %>% as.character
#melted_sub$EMBO[which(melted_sub$Gene %in% mgenes)] <- "EMBO1"
mgenes <- df_embo %>% subset(cluster2 == 2) %$% rowname %>% as.character
melted_sub$EMBO[which(melted_sub$Gene %in% mgenes)] <- "EMBO2"
#mgenes <- df_embo %>% subset(cluster2 == 3) %$% rowname %>% as.character
#melted_sub$EMBO[which(melted_sub$Gene %in% mgenes)] <- "EMBO3"
#mgenes <- df_embo %>% subset(cluster2 == 4) %$% rowname %>% as.character
#melted_sub$EMBO[which(melted_sub$Gene %in% mgenes)] <- "EMBO4"


p5 <- ggplot(data = melted_sub,aes(Bar,Gene))+
  geom_tile(aes(fill = EMBO, color =EMBO )) +
  facet_grid(Classification~.,scales = "free",space = "free") +
  theme_ohchibi() + 
  scale_x_discrete(expand = c(0,0)) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.background.y = element_blank(),
    strip.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  )  +
  xlab(label = "EMBO clusters") +
  scale_fill_manual(values = c("black"),na.value = "#00000000")+
  scale_color_manual(values = c("black"),na.value = "#00000000")


mtemp <- melted_sub %>% subset(EMBO == "EMBO2") %>% droplevels
mtemp[,c(1,6)] %>% unique %$% Classification %>% table

melted_sub[,c(1,6)] %>% unique %$% Classification %>% table

composition <- egg::ggarrange(p,p2,p3,p4,p5,nrow =1 ,widths = c(1,0.03,0.03,0.03,0.03))
oh.save.pdf(p = composition,outname = "rnaseq_heatmap_sets.other.pdf",outdir = "../figures/",width = 16,height = 24)


rm(list=ls())
dev.off()
