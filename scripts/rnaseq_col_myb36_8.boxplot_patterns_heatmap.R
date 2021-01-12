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
melted_sub <- merge(melted_sub,df_class, by = "Gene") 

ggplot(data =melted_sub,mapping = aes(group,value)) +
  geom_boxplot() +
  facet_grid(.~Classification)

melted_sub$Bacteria <- melted_sub$Bacteria %>% factor(levels = c("NB","SC"))


p <- chibi.boxplot(Map = melted_sub,x_val = "Genotype",
              y_val = "value",col_val = "Bacteria",
              facet_formula = "Classification",strip_text_size = 10,
              median_colored_as_points = T,size_median = 2,alpha_point = 1,
              size_boxplot = 1,mpalette = c("#525252","#b48c36"))



oh.save.pdf(p = p + theme(legend.position = "none"),outname = "boxplot_rnaseq_patterns.pdf",outdir = "../figures/",width = 24,height = 8)
