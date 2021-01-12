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


######### Gene ontology enrichment ############
####### Enrichment pattern

mlist <- list(
  C1= clean_list$bacteria_up_unique,
  C2 = clean_list$bacteria_down_unique,
  C3 = clean_list$genotype_up_unique,
  C4 = clean_list$genotype_down_unique,
  C5 = clean_list$bacteria_up_genotype_up,
  C6 = clean_list$bacteria_up_genotype_down,
  C7 = clean_list$bacteria_down_genotype_up,
  C8 = clean_list$bacteria_down_genotype_down
)

cg <- compareCluster(geneCluster=mlist,
                     fun="enrichGO",
                     keyType       = "TAIR",
                     OrgDb         = org.At.tair.db,
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.1)

p <- dotplot(cg, showCategory=10, includeAll=TRUE)+
  scale_color_paletteer_c(package = "viridis",palette = "plasma",na.value = "#BFBFBF") 

oh.save.pdf(p = p,outname = "rnaseq_go_enrichments_top10.pdf",outdir = "../figures/",width = 10,height = 12)

bp2 <- simplify(cg, cutoff=0.7, by="p.adjust", select_fun=min)


p <- dotplot(bp2, showCategory=15, includeAll=TRUE)+
  scale_color_paletteer_c(package = "viridis",palette = "plasma",na.value = "#BFBFBF") 

oh.save.pdf(p = p,outname = "rnaseq_go_enrichments_top15.simplified.pdf",outdir = "../figures/",width = 10,height = 14)

rm(list=ls())
dev.off()
