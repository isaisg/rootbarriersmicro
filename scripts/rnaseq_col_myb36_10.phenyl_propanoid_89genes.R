library(ohchibi)
library(DESeq2)
library(clusterProfiler)
library(org.At.tair.db)
library(scales)
library(egg)
library(ggpubr)

setwd('/Users/isaisalasgonzalez/Documents/rootbarriersmicro/scripts/')

set.seed(130816)
Dat <-readRDS(file = "../cleandata/dat_rnaseq_col_myb36.RDS")
Res <- readRDS(file = "../cleandata/res_rnaseq_col_myb36.RDS")

melted_av <- Dat$melted_av

df_mgenes <- readRDS(file = "../rawdata/df_88genes_aracyc_regulators_pp_pathway.RDS")
#Append AT1G07890

df_mgenes <- df_mgenes %>% 
  rbind(data.frame(Gene = "AT1G07890",Symbol = "AT1G07890",SymbolFirst = "AT1G07890"))

genes <- df_mgenes$Gene %>% as.character

which(melted_av$Gene %in% genes) %>%
  melted_av[.,] %>% droplevels %>%
  chibi.boxplot(Map = .,x_val = "Genotype",y_val = "value",col_val = "Bacteria")

melted_av <- which(melted_av$Gene %in% genes) %>%
  melted_av[.,] %>% droplevels

melted_av$Name <- match(melted_av$Gene,df_mgenes$Gene) %>%
  df_mgenes$SymbolFirst[.] %>% as.character



##### Cluster
Tab <- acast(data = melted_av,formula =Name~group,
             fun.aggregate = mean,value.var = "value" )

mclust_genes <- hclust(d = as.dist(1-cor(Tab %>% t)),method = "ward.D2") 
mclust_genes %>% plot
df_clust <- cutree(tree = mclust_genes,k =11) %>%
  data.frame(Name = names(.),Cluster = paste0("C",.),row.names = NULL)
df_clust <- df_clust[,-1]
order_genes <- mclust_genes$order %>% mclust_genes$labels[.]


melted_av <- merge(melted_av,df_clust, by = "Name")

melted_av$Name <- melted_av$Name %>% factor(levels = order_genes)

order_clusters <- with(melted_av,order(Name)) %>%
  melted_av[.,] %$% Cluster %>% as.character %>% unique

melted_av$Cluster <- melted_av$Cluster %>% factor(levels = order_clusters %>% rev)
melted_av$value %>% sort %>% plot

melted_av$Label <- round(melted_av$value,3)

order_genes <- with(melted_av,order(Name)) %>%
  melted_av$Gene[.] %>% as.character %>% unique
melted_av$Gene <- melted_av$Gene %>% factor(levels = order_genes)


p <- ggplot(data = melted_av,mapping = aes(Bacteria,Gene)) + 
  geom_tile(aes(fill = value),color = "#00000000") + 
  facet_grid(Cluster~Genotype,scales = "free",space = "free") +
  theme_ohchibi() + 
  geom_text(aes(label = Label)) +
  scale_x_discrete(expand = c(0,0)) +
  scale_fill_paletteer_c("pals::kovesi.diverging_bwr_55_98_c37",
                         limits = c(-0.5,0.5),oob = squish) +
  theme(
    axis.text.y = element_text(size = 13),
    #axis.ticks.y = element_blank(),
    strip.background.x = element_blank(),
    strip.text.x = element_text(size = 15,family = "Arial",face = "bold",angle = 0),
    strip.background.y = element_blank(),
    strip.text.y =element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) 

oh.save.pdf(p = p,outname = "phenylpropanoid_heatmap.89genes_locus.pdf",outdir = "../figures/",width = 8,height = 16)

######## Classify each gene to cluster #######
df_class <- Res$df_classification
df_class$ClassCluster <- df_class$Classification %>% 
  gsub(pattern = "bacteria_up_unique",replacement = "C1") %>%
  gsub(pattern = "bacteria_down_unique",replacement = "C2") %>%
  gsub(pattern = "genotype_up_unique",replacement = "C3") %>%
  gsub(pattern = "genotype_down_unique",replacement = "C4") %>%
  gsub(pattern = "bacteria_up_genotype_up",replacement = "C5") %>%
  gsub(pattern = "bacteria_up_genotype_down",replacement = "C6") %>%
  gsub(pattern = "bacteria_down_genotype_up",replacement = "C7") %>%
  gsub(pattern = "bacteria_down_genotype_down",replacement = "C8") 

df_noms <- melted_av[,c(2,1)] %>% unique
df_noms$ClassCluster <- match(df_noms$Gene,df_class$Gene) %>%
  df_class$ClassCluster[.]

df_noms <- with(df_noms,order(Name) %>% rev) %>%
  df_noms[.,]

write.table(x = df_noms,file = "../cleandata/df_names_phenylpropanoid2cluster_89genes.csv",
            append = F,quote = F,sep = "\t",row.names = F,col.names = T)

rm(list=ls())
dev.off()
