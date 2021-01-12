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

dds <- Dat$dds
melted_av <- Dat$melted_av
### Compute specific contrasts my pooling
res_bacteria_main <- results(dds, 
                             contrast = list(c("groupmyb36_SC", "groupCol_0_SC"), 
                                             c("groupmyb36_NB", "groupCol_0_NB")), 
                             listValues = c(1/2, -1/2)) %>% as.data.frame

res_interaction <- results(dds, 
                           contrast = list(c("groupmyb36_SC", "groupCol_0_NB"), 
                                           c("groupCol_0_SC", "groupmyb36_NB")), 
                           listValues = c(1/2, -1/2)) %>% as.data.frame
res_genotype_main <- results(dds, 
                             contrast = list(c("groupmyb36_SC", "groupmyb36_NB"), 
                                             c("groupCol_0_SC", "groupCol_0_NB")), 
                             listValues = c(1/2, -1/2)) %>% as.data.frame

### Genotype effect 
down_genotype <- res_genotype_main %>% subset(padj < 0.05 & log2FoldChange < 0) %>% rownames
up_genotype <- res_genotype_main %>% subset(padj < 0.05 & log2FoldChange > 0) %>% rownames


which(melted_av$Gene %in% down_genotype) %>%
  melted_av[.,] %>%
  chibi.boxplot(x_val = "Genotype",y_val = "value",style = "open"
                ,facet_formula = "Bacteria")

which(melted_av$Gene %in% up_genotype) %>%
  melted_av[.,] %>%
  chibi.boxplot(x_val = "Genotype",y_val = "value",style = "open"
                ,facet_formula = "Bacteria")

### Bacteria effect
down_bacteria <- res_bacteria_main %>% subset(padj < 0.05 & log2FoldChange < 0) %>% rownames
up_bacteria <- res_bacteria_main %>% subset(padj < 0.05 & log2FoldChange > 0) %>% rownames


which(melted_av$Gene %in% down_bacteria) %>%
  melted_av[.,] %>%
  chibi.boxplot(x_val = "Genotype",y_val = "value",style = "open"
                ,facet_formula = "Bacteria")

which(melted_av$Gene %in% up_bacteria) %>%
  melted_av[.,] %>%
  chibi.boxplot(x_val = "Genotype",y_val = "value",style = "open"
                ,facet_formula = "Bacteria")

### Interaction ###

#Almost no interaction!
interaction <- res_interaction %>% subset(padj < 0.1) %>% rownames


which(melted_av$Gene %in% interaction) %>%
  melted_av[.,] %>%
  chibi.boxplot(x_val = "Genotype",y_val = "value",style = "open"
                ,facet_formula = "Bacteria")



### Create list ###
mlista <- list(
  up_bacteria = up_bacteria,
  down_bacteria = down_bacteria,
  up_genotype = up_genotype,
  down_genotype = down_genotype
)



########### Create sets ########
clean_genotype <- union(mlista$up_genotype,mlista$down_genotype)
bacteria_up_unique <- which(!(mlista$up_bacteria %in% clean_genotype)) %>% mlista$up_bacteria[.]
bacteria_up_up <- intersect(mlista$up_bacteria,mlista$up_genotype)
bacteria_up_down <- intersect(mlista$up_bacteria,mlista$down_genotype)

bacteria_down_unique <- which(!(mlista$down_bacteria %in% clean_genotype)) %>% mlista$down_bacteria[.]
bacteria_down_up <- intersect(mlista$down_bacteria,mlista$up_genotype)
bacteria_down_down <- intersect(mlista$down_bacteria,mlista$down_genotype)


clean_bacteria <- union(mlista$up_bacteria,mlista$down_bacteria)
genotype_up_unique <- which(!(mlista$up_genotype %in% clean_bacteria)) %>% mlista$up_genotype[.]
genotype_up_up <- intersect(mlista$up_genotype,mlista$up_bacteria)
genotype_up_down <- intersect(mlista$up_genotype,mlista$down_bacteria)

genotype_down_unique <- which(!(mlista$down_genotype %in% clean_bacteria)) %>% mlista$down_genotype[.]
genotype_down_up <- intersect(mlista$down_genotype,mlista$up_bacteria)
genotype_down_down <- intersect(mlista$down_genotype,mlista$down_bacteria)


clean_list <- list(
  bacteria_up_unique = bacteria_up_unique,
  bacteria_down_unique = bacteria_down_unique,
  genotype_up_unique = genotype_up_unique,
  genotype_down_unique = genotype_down_unique,
  bacteria_up_genotype_up = bacteria_up_up,
  bacteria_up_genotype_down = bacteria_up_down,
  bacteria_down_genotype_up = bacteria_down_up,
  bacteria_down_genotype_down = bacteria_down_down
)



## Build the heatmap
df_class <- NULL
for(nom in names(clean_list)){
  df_class <- data.frame(Classification = nom,Gene =
                           clean_list[[nom]],row.names = NULL) %>%
    rbind(df_class,.)
}


### Save structures
res <- list(
  RawContrasts = mlista,
  SetsList = clean_list,
  df_classification = df_class
)
saveRDS(object = res,file = "../cleandata/res_rnaseq_col_myb36.RDS")