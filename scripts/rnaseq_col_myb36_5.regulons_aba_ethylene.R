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

## read aba genes ##
Res_aba <- readRDS(file = "../cleandata/aba_robust_genes.RDS")
chosen_up <- Res_aba$aba_genes_up
chosen_down <- Res_aba$aba_genes_down

melted_av$Type <- "ABA downregulated"

a <- which(melted_av$Gene %in% chosen_down) %>%
  melted_av[.,] %>%
  subset(value < 1  & value > -1) %>%
chibi.boxplot(Map = .,x_val = "Genotype",
              y_val = "value",col_val = "Bacteria",facet_formula = "Type",
              strip_text_size = 20,mpalette = c("#525252","#b48c36"),
              median_colored_as_points = T,size_median = 2,alpha_point = 1,size_boxplot = 1) 






b <- which(melted_av$Gene %in% chosen_down) %>%
  melted_av[.,] %>% 
  subset(value < 1  & value > -1) %>%
  chibi.boxplot(x_val = "Bacteria",y_val = "value",style = "open")

melted_av$Type <- "ABA upregulated"


c <- which(melted_av$Gene %in% chosen_up) %>%
  melted_av[.,] %>%
  subset(value < 1  & value > -1) %>%
  chibi.boxplot(Map = .,x_val = "Genotype",
                y_val = "value",col_val = "Bacteria",facet_formula = "Type",
                strip_text_size = 20,mpalette = c("#525252","#b48c36"),
                median_colored_as_points = T,size_median = 2,alpha_point = 1,size_boxplot = 1) 





d <- which(melted_av$Gene %in% chosen_up) %>%
  melted_av[.,] %>% 
  subset(value < 1  & value > -1) %>%
  chibi.boxplot(x_val = "Bacteria",y_val = "value",style = "open")


oh.save.pdf(p = a,outname = "rnaseq_boxplot_aba_down.pdf",outdir = "../figures/",width = 10,height = 16)
oh.save.pdf(p = c,outname = "rnaseq_boxplot_aba_up.pdf",outdir = "../figures/",width = 10,height = 16)



###### Ethylene #####
mregs <- readRDS(file = "../../variovorax/rnaseq/cleandata/regulons_rnaseq.RDS")
mgenes <- mregs$ethylene

melted_av$Type <- "Ethylene upregulated"

e <- which(melted_av$Gene %in% mgenes) %>%
  melted_av[.,] %>%
  subset(value < 1  & value > -1) %>%
  chibi.boxplot(Map = .,x_val = "Genotype",
                y_val = "value",col_val = "Bacteria",facet_formula = "Type",
                strip_text_size = 20,mpalette = c("#525252","#b48c36"),
                median_colored_as_points = T,size_median = 2,alpha_point = 1,size_boxplot = 1) 



oh.save.pdf(p = e,outname = "rnaseq_boxplot_ethylene.pdf",outdir = "../figures/",width = 10,height = 16)


f <- which(melted_av$Gene %in% mgenes) %>%
  melted_av[.,] %>% 
  subset(value < 1  & value > -1) %>%
  chibi.boxplot(x_val = "Bacteria",y_val = "value",style = "open")


which(melted_av$Gene %in% mgenes) %>%
  melted_av[.,] %>%
  ggplot(data = .,aes(x=Bacteria,y = value))+
  geom_line(aes(group =Gene )) +
  geom_point() +
  facet_grid(.~Genotype) +
  theme_ohchibi()


Dat <-readRDS(file = "../cleandata/dat_rnaseq_col_myb36.RDS")
Res <- readRDS(file = "../cleandata/res_rnaseq_col_myb36.RDS")
df_embo <- read.table(file = "../rawdata/embo_930genes.csv",header = T,sep = "\t")

melted_av <- Dat$melted_av 
melted_av$EMBO <- NA

for(i in (df_embo$cluster2 %>% unique)){
  mgenes <- df_embo %>% subset(cluster2 == i) %$% rowname %>% as.character
  melted_av$EMBO[which(melted_av$Gene %in% mgenes)] <- paste0("EMBO",i)
}


melted_avembo <- melted_av[which(!is.na(melted_av$EMBO)),]
melted_avembo$Type <- "Cluster2 EMBO"

p <- melted_avembo %>%
  subset(EMBO == "EMBO2") %>% 
  #subset(value < 1  & value > -1) %>%
  chibi.boxplot(Map = .,x_val = "Genotype",
                y_val = "value",col_val = "Bacteria",facet_formula = "Type",
                strip_text_size = 20,mpalette = c("#525252","#b48c36"),
                median_colored_as_points = T,size_median = 2,alpha_point = 1,size_boxplot = 1) 


oh.save.pdf(p = p,outname = "embo2.pdf",outdir = "../figures/",width = 10,height = 16)

### create composition ##
composition <- egg::ggarrange(a  + theme(legend.position = "none"),
               c + theme(legend.position = "none"),
               e + theme(legend.position = "none"),
               p + theme(legend.position = "none"),
               nrow = 1)
oh.save.pdf(p = composition,outname = "boxplots_rnaseq_all.pdf",outdir = "../figures/",width = 24,height = 8)
