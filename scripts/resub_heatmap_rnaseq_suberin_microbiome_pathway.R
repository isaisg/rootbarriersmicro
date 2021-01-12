library(ohchibi)
library(DESeq2)
library(clusterProfiler)
library(org.At.tair.db)
library(scales)
library(egg)
library(ggpubr)
library(ggvenn)

setwd('/Users/isaisalasgonzalez/Documents/rootbarriersmicro/scripts/')

set.seed(130816)
Dat <-readRDS(file = "../cleandata/dat_rnaseq_col_myb36.RDS")
Res <- readRDS(file = "../cleandata/res_rnaseq_col_myb36.RDS")

melted_av <- Dat$melted_av


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

dim(df_class)
Res$SetsList %>% unlist %>% length

genes <- which(df_class$ClassCluster  %in% c("C5","C6","C7","C8")) %>%
  df_class[.,] %>% droplevels %$% Gene %>% as.character


melted_av <- which(melted_av$Gene %in% genes) %>%
  melted_av[.,] %>% droplevels

melted_av$ClassCluster <- match(melted_av$Gene,df_class$Gene) %>%
  df_class$ClassCluster[.] %>% as.character



##### Cluster
Tab <- acast(data = melted_av,formula =Gene~group,
             fun.aggregate = mean,value.var = "value" )

## Approach where we identify contrast ###
mc7 <- df_class %>% subset(ClassCluster == "C7") %$% Gene %>% as.character
Tab <- which(rownames(Tab) %in% mc7) %>%
  Tab[.,]
#chosen_genes <- which((Tab[,1] > Tab[,2]) & (Tab[,3] > Tab[,1]) & (Tab[,2] < Tab[,4]) & (Tab[,4]>Tab[,1]) & (Tab[,3] > Tab[,2])) %>%  names


mclust_genes <- hclust(d = as.dist(1-cor(Tab %>% t)),method = "ward.D2")
mclust_genes %>% plot
df_clust <- cutree(tree = mclust_genes,k =6) %>%
  data.frame(Gene = names(.),Cluster = paste0("C",.),row.names = NULL)
df_clust <- df_clust[,-1]
order_genes <- mclust_genes$order %>% mclust_genes$labels[.]


melted_av <- merge(melted_av,df_clust, by = "Gene")

melted_av$Gene <- melted_av$Gene %>% factor(levels = order_genes)

order_clusters <- with(melted_av,order(Gene)) %>%
  melted_av[.,] %$% Cluster %>% as.character %>% unique

melted_av$Cluster <- melted_av$Cluster %>% factor(levels = order_clusters %>% rev)
melted_av$value %>% sort %>% plot

melted_av$Label <- round(melted_av$value,3)



p <- ggplot(data = melted_av,mapping = aes(Bacteria,Gene)) +
  geom_tile(aes(fill = value),color = "00000000") +
  facet_grid(Cluster~Genotype,scales = "free",space = "free") +
  theme_ohchibi() +
  scale_x_discrete(expand = c(0,0)) +
  scale_fill_paletteer_c("pals::kovesi.diverging_bwr_55_98_c37",
                         limits = c(-0.5,0.5),oob = squish) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.background.x = element_blank(),
    strip.text.x = element_text(size = 15,family = "Arial",face = "bold",angle = 0),
    strip.background.y = element_blank(),
    strip.text.y  =  element_text(size = 15,family = "Arial",face = "bold",angle = 270),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )


mgenes <- df_clust %>% subset(Cluster == "C2") %$% Gene %>% as.character
#write(x = mgenes,file = "../cleandata/48_genes_overwrite.txt")
melted_sub <- which(melted_av$Gene %in% mgenes) %>%
  melted_av[.,] %>% droplevels

p_boxplot <- chibi.boxplot(Map = melted_sub,x_val = "Genotype",
                           y_val = "value",col_val = "Bacteria",
                           strip_text_size = 10,
                           median_colored_as_points = T,size_median = 2,alpha_point = 1,
                           size_boxplot = 1,mpalette = c("#525252","#b48c36")) +
  theme(axis.title.x = element_blank(),
        legend.position = "none") 

p <- ggplot(data = melted_sub,mapping = aes(Bacteria,Gene)) +
  geom_tile(aes(fill = value),color = "00000000") +
  facet_grid(Cluster~Genotype,scales = "free",space = "free") +
  geom_text(aes(label = round(value,2))) +
  theme_ohchibi() +
  scale_x_discrete(expand = c(0,0)) +
  scale_fill_paletteer_c("pals::kovesi.diverging_bwr_55_98_c37",
                         limits = c(-0.5,0.5),oob = squish) +
  theme(
    #axis.text.y = element_blank(),
    #axis.ticks.y = element_blank(),
    strip.background.x = element_blank(),
    strip.text.x = element_text(family = "Helvetica",face = "plain",size = 15),
    strip.background.y = element_blank(),
    strip.text.y  =  element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(family = "Helvetica",face = "plain",size = 10)
  )
oh.save.pdf(p = p,outname = "resub_heatmap_suberin_microbiome_pathway.pdf",outdir = "../figures/",width = 8,height = 10)
