library(ohchibi)
library(scales)
library(ggtree)

##This is the script of the enrichment in the soil

set.seed(130816)
setwd('/Users/isaisalasgonzalez/Documents/rootbarriersmicro/scripts/')

### Load the projection in capspace
mconst <- readRDS(file = "../cleandata/plots_constellations_correlations_amplicon_wild.RDS")


#Read the deseq2 results
df_deseq <- read.table(file = "../cleandata/revision_deseq2_results_censusmutants_genotypes_inside_fractions_all.tsv",header = T,sep = "\t")

#Define threshold for enrichment 
threshold <- 0.1
hclust_method <- "ward.D2"


#### Root ##########
#Constellation of the cap
data_const <- mconst$constellation_root_amplicon$data %>% 
  subset(Genotype != "Col_0") %>% droplevels
dist_const <- data_const[,c("CAP1","CAP2")]
rownames(dist_const) <- data_const$Genotype
dist_const <- dist(dist_const)

df_deseq_sub <- df_deseq %>% 
  subset(Fraction == "Root") %>% droplevels

#Deseq results
level <- "ASV"
df_temp <- df_deseq_sub %>% subset(Level == level) %>%
  droplevels

#Ad just pvalue
df_temp$padjSub <- df_temp$pvalue %>% p.adjust(method = "fdr")


df_sig <- df_temp %>% subset(padjSub < threshold) %>% droplevels
df_sig$Direction <- -1
df_sig$Direction[which(df_sig$log2FoldChange > 0)] <- 1
mids <- df_sig$Id %>% unique

#Subset the original matrix to obtain fold changes
df_lfc <- which(df_temp$Id %in% mids) %>%
  df_temp[.,] %>% droplevels

#Create matrices of enrichment and lfc
Tab_man <- acast(data = df_sig,formula = Id~Genotype,value.var = "Direction",fill = 0)
Tab_lfc <- acast(data = df_lfc,formula = Id~Genotype,value.var = "log2FoldChange")


mclust_genotype <- hclust(d = dist_const,method = "ward.D2")
df_mclust_genotype <- cutree(tree = mclust_genotype,k = 6) %>% 
  data.frame(Genotype = names(.),CG = paste0("CGen",.),row.names = NULL)
df_mclust_genotype <- df_mclust_genotype[,-1]


#Cluster the asv patterns
mclust_asv <- hclust(d = (1-cor(Tab_lfc  %>%t)) %>% as.dist,method = "ward.D2") 
mclust_asv %>%plot

#Merge the df_clust structures with the structure used to plot
melted <- Tab_lfc %>% melt
colnames(melted) <- c("Id","Genotype","value")

melted <- merge(melted,df_mclust_genotype, by = "Genotype")



#Append the taxonomy information
Dat <- readRDS(file = "../cleandata/dat_censusmutants_16s.RDS")
df_tax <- Dat$df_tax
df_tax <- match(melted$Id,df_tax$Id) %>%
  df_tax[.,]
melted <- cbind(melted,df_tax[,-c(1,2)])

#Order the genotype and asvs
order_genotype <- mclust_genotype$order %>% mclust_genotype$labels[.]
order_asv <- mclust_asv$order %>% mclust_asv$labels[.]

melted$Genotype <- melted$Genotype %>% factor(levels = order_genotype)
melted$Id <- melted$Id %>% factor(levels = order_asv)

#Determine the order of the clusters
order_cg <- with(melted,order(Genotype)) %>%
  melted[.,] %$% CG %>% unique

melted$CG <- melted$CG %>% factor(levels = order_cg %>% rev)


## Add the number identifier
df_num <- read.table(file = "../rawdata/df_genotypes_to_classification.csv",header = T,sep = ",")
df_num$Number <- 1:20
melted$NumberGenotype <- match(melted$Genotype,df_num$Genotype) %>%
  df_num$Number[.]


order_numgeno <- with(melted,order(Genotype)) %>%
  melted$NumberGenotype[.] %>% unique
melted$NumberGenotype <- melted$NumberGenotype %>% factor(levels = order_numgeno)

#Append significance 
df_sig$Significance <- "q < 0.1"
melted <- merge(melted,df_sig[,c("Genotype","Id","Significance")],by = c("Genotype","Id"), all.x = TRUE)


### Plot 
p_heatmap <- ggplot(data = melted,aes(Id,NumberGenotype)) +
  geom_raster(aes(fill =value ), stat = "identity") +
  geom_tile(aes(color = Significance),fill = '#00000000', size = 1,width = 0.85,height = 0.85) + 
  facet_grid(CG~.,space = "free",scales = "free")  +
  scale_color_manual(values = c("black"))+
  scale_fill_paletteer_c("pals::kovesi.diverging_bwr_55_98_c37",
                         limits = c(-5,5),oob = squish,na.value = "white") +
  theme_ohchibi(size_panel_border = 2) +
  theme(
    strip.background.x = element_blank(),
    strip.background.y = element_blank(),
    strip.text.x = element_text(family = "Helvetica",face = "plain",size = 15),
    #strip.text.y = element_text(family = "Helvetica",face = "plain",size = 20),
    strip.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(family = "Helvetica",face = "plain",
                               size = 11,angle = 90,vjust = 0.5,hjust = 1),
    #axis.text.x = element_blank(),
    #axis.ticks.x = element_blank(),
    axis.text.y = element_text(family = "Helvetica",face = "plain",size = 18),
    #legend.position = "none"
  ) +
  scale_x_discrete(expand = c(0,0))  +
  scale_y_discrete(expand = c(0,0))


#Prepare the dendrogram to put in the same graph
tree <- mclust_genotype %>% as.phylo
p_tree_genotype <- ggtree(tree,ladderize = F,size = 1.2) + 
  #The expand parth is fundamental to make the tree fit the composition with the other plots
  scale_y_continuous(expand = c(0.026,0.026))

melted_asv <- melted

#Create composition
composition_asv  <- egg::ggarrange(p_tree_genotype,p_heatmap,nrow = 1,widths = c(0.1,1))
oh.save.pdf(p = composition_asv,outname = "revision_amplicon_heatmap_deseq2_asv.pdf",
            outdir = "../figures/",width = 16,height = 12)
p_heatmap_asv <- p_heatmap

### Append the total RA
Dat_sub <- Dat$RelativeAbundance %>% subset.Dataset(Fraction == "Root",drop = T,clean = T)
df_ra <- Dat_sub$Tab %>% melt
colnames(df_ra) <- c("Id","DADA2_Header","RA")
df_ra <- merge(df_ra,Dat_sub$Map, by = "DADA2_Header")
df_sum <- aggregate(RA~Id,df_ra,mean) 
df_sum <- which(df_sum$Id %in% mids) %>%
  df_sum[.,] %>% droplevels

melted_sub <- melted %>% subset(Significance == "q < 0.1") %>% droplevels
mgenos <- melted_sub$Genotype %>% levels
df_cumra <- NULL
for(geno in mgenos){
  chosen_ids <- melted_sub %>% subset(Genotype == geno) %>% droplevels %$% Id
  mval <- which(df_sum$Id %in% chosen_ids) %>%
    df_sum$RA[.] %>% sum
  df_cumra <- data.frame(Genotype = geno,CumRA = mval) %>%
    rbind(df_cumra,.)
}


#### Genus level #######
level <- "Genus"
df_temp <- df_deseq_sub %>% subset(Level == level) %>%
  droplevels

#Ad just pvalue
df_temp$padjSub <- df_temp$pvalue %>% p.adjust(method = "fdr")


df_sig <- df_temp %>% subset(padjSub < threshold) %>% droplevels
df_sig$Direction <- -1
df_sig$Direction[which(df_sig$log2FoldChange > 0)] <- 1
mids <- df_sig$Id %>% unique

#Subset the original matrix to obtain fold changes
df_lfc <- which(df_temp$Id %in% mids) %>%
  df_temp[.,] %>% droplevels

#Create matrices of enrichment and lfc
Tab_man <- acast(data = df_sig,formula = Id~Genotype,value.var = "Direction",fill = 0)
Tab_lfc <- acast(data = df_lfc,formula = Id~Genotype,value.var = "log2FoldChange")


#Cluster the asv patterns
mclust_asv <- hclust(d = (1-cor(Tab_lfc  %>%t)) %>% as.dist,method = "ward.D2") 
mclust_asv %>%plot

#Merge the df_clust structures with the structure used to plot
melted <- Tab_lfc %>% melt
colnames(melted) <- c("Id","Genotype","value")


#Order the genotype and asvs
order_asv <- mclust_asv$order %>% mclust_asv$labels[.]

melted$Genotype <- melted$Genotype %>% factor(levels = melted_asv$Genotype %>% levels)
melted$Id <- melted$Id %>% factor(levels = order_asv)


melted$CG <- match(melted$Genotype,melted_asv$Genotype) %>%
  melted_asv$CG[.]

melted$NumberGenotype <- match(melted$Genotype,melted_asv$Genotype) %>%
  melted_asv$NumberGenotype[.]


#Append significance 
df_sig$Significance <- "q < 0.1"
melted <- merge(melted,df_sig[,c("Genotype","Id","Significance")],by = c("Genotype","Id"), all.x = TRUE)

mat <- melted$Id %>% as.character %>% strsplit(split = "\\;") %>% unlist %>%
  matrix(data = .,byrow = T,ncol = 7)
nid <- paste0(mat[,3]," ",mat[,7])  %>% gsub(pattern = " p__",replacement = "p__") %>%
  gsub(pattern = "__",replacement = "_")
melted$NId <- nid

order_nids <- with(melted,order(Id)) %>%
  melted$NId[.] %>% unique
melted$NId <- melted$NId %>% factor(levels = order_nids)

### Plot 
p_heatmap <- ggplot(data = melted,aes(NId,NumberGenotype)) +
  geom_raster(aes(fill =value ), stat = "identity") +
  geom_tile(aes(color = Significance),fill = '#00000000', size = 1,width = 0.85,height = 0.85) + 
  facet_grid(CG~.,space = "free",scales = "free")  +
  scale_color_manual(values = c("black"))+
  scale_fill_paletteer_c("pals::kovesi.diverging_bwr_55_98_c37",
                         limits = c(-5,5),oob = squish,na.value = "white") +
  theme_ohchibi(size_panel_border = 2) +
  theme(
    strip.background.x = element_blank(),
    strip.background.y = element_blank(),
    strip.text.x = element_text(family = "Helvetica",face = "plain",size = 15),
    #strip.text.y = element_text(family = "Helvetica",face = "plain",size = 20),
    strip.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(family = "Helvetica",face = "plain",
                               size = 11,angle = 90,vjust = 0.5,hjust = 1),
    #axis.text.x = element_blank(),
    #axis.ticks.x = element_blank(),
    axis.text.y = element_text(family = "Helvetica",face = "plain",size = 18),
    legend.position = "none"
  ) +
  scale_x_discrete(expand = c(0,0))  +
  scale_y_discrete(expand = c(0,0))


#Prepare the dendrogram to put in the same graph
tree <- mclust_genotype %>% as.phylo
p_tree_genotype <- ggtree(tree,ladderize = F,size = 1.2) + 
  #The expand parth is fundamental to make the tree fit the composition with the other plots
  scale_y_continuous(expand = c(0.026,0.026))

melted_asv <- melted

#Create composition
composition_genus <- egg::ggarrange(p_tree_genotype,p_heatmap,nrow = 1,widths = c(0.1,1))
oh.save.pdf(p = composition_genus,outname = "revision_amplicon_heatmap_deseq2_genus.pdf",
            outdir = "../figures/",width = 16,height = 12)

p_heatmap_genus <- p_heatmap

##### Family level ########
level <- "Family"
df_temp <- df_deseq_sub %>% subset(Level == level) %>%
  droplevels

#Ad just pvalue
df_temp$padjSub <- df_temp$pvalue %>% p.adjust(method = "fdr")


df_sig <- df_temp %>% subset(padjSub < threshold) %>% droplevels
df_sig$Direction <- -1
df_sig$Direction[which(df_sig$log2FoldChange > 0)] <- 1
mids <- df_sig$Id %>% unique

#Subset the original matrix to obtain fold changes
df_lfc <- which(df_temp$Id %in% mids) %>%
  df_temp[.,] %>% droplevels

#Create matrices of enrichment and lfc
Tab_man <- acast(data = df_sig,formula = Id~Genotype,value.var = "Direction",fill = 0)
Tab_lfc <- acast(data = df_lfc,formula = Id~Genotype,value.var = "log2FoldChange")


#Cluster the asv patterns
mclust_asv <- hclust(d = (1-cor(Tab_lfc  %>%t)) %>% as.dist,method = "ward.D2") 
mclust_asv %>%plot

#Merge the df_clust structures with the structure used to plot
melted <- Tab_lfc %>% melt
colnames(melted) <- c("Id","Genotype","value")


#Order the genotype and asvs
order_asv <- mclust_asv$order %>% mclust_asv$labels[.]

melted$Genotype <- melted$Genotype %>% factor(levels = melted_asv$Genotype %>% levels)
melted$Id <- melted$Id %>% factor(levels = order_asv)


melted$CG <- match(melted$Genotype,melted_asv$Genotype) %>%
  melted_asv$CG[.]

melted$NumberGenotype <- match(melted$Genotype,melted_asv$Genotype) %>%
  melted_asv$NumberGenotype[.]


#Append significance 
df_sig$Significance <- "q < 0.1"
melted <- merge(melted,df_sig[,c("Genotype","Id","Significance")],by = c("Genotype","Id"), all.x = TRUE)

mat <- melted$Id %>% as.character %>% strsplit(split = "\\;") %>% unlist %>%
  matrix(data = .,byrow = T,ncol = 6)
nid <- paste0(mat[,3]," ",mat[,6])  %>% gsub(pattern = " p__",replacement = "p__") %>%
  gsub(pattern = "__",replacement = "_")
melted$NId <- nid

order_nids <- with(melted,order(Id)) %>%
  melted$NId[.] %>% unique
melted$NId <- melted$NId %>% factor(levels = order_nids)


### Plot 
p_heatmap <- ggplot(data = melted,aes(NId,NumberGenotype)) +
  geom_raster(aes(fill =value ), stat = "identity") +
  geom_tile(aes(color = Significance),fill = '#00000000', size = 1,width = 0.85,height = 0.85) + 
  facet_grid(CG~.,space = "free",scales = "free")  +
  scale_color_manual(values = c("black"))+
  scale_fill_paletteer_c("pals::kovesi.diverging_bwr_55_98_c37",
                         limits = c(-5,5),oob = squish,na.value = "white") +
  theme_ohchibi(size_panel_border = 2) +
  theme(
    strip.background.x = element_blank(),
    strip.background.y = element_blank(),
    strip.text.x = element_text(family = "Helvetica",face = "plain",size = 15),
    #strip.text.y = element_text(family = "Helvetica",face = "plain",size = 20),
    strip.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(family = "Helvetica",face = "plain",
                               size = 11,angle = 90,vjust = 0.5,hjust = 1),
    #axis.text.x = element_blank(),
    #axis.ticks.x = element_blank(),
    axis.text.y = element_text(family = "Helvetica",face = "plain",size = 18),
    legend.position = "none"
  ) +
  scale_x_discrete(expand = c(0,0))  +
  scale_y_discrete(expand = c(0,0))


#Prepare the dendrogram to put in the same graph
tree <- mclust_genotype %>% as.phylo
p_tree_genotype <- ggtree(tree,ladderize = F,size = 1.2) + 
  #The expand parth is fundamental to make the tree fit the composition with the other plots
  scale_y_continuous(expand = c(0.026,0.026))

melted_asv <- melted

#Create composition
composition_family <- egg::ggarrange(p_tree_genotype,p_heatmap,nrow = 1,widths = c(0.1,1))
oh.save.pdf(p = composition_family,outname = "revision_amplicon_heatmap_deseq2_family.pdf",
            outdir = "../figures/",width = 16,height = 12)
p_heatmap_family <- p_heatmap

##### Order level ########
level <- "Order"
df_temp <- df_deseq_sub %>% subset(Level == level) %>%
  droplevels

#Ad just pvalue
df_temp$padjSub <- df_temp$pvalue %>% p.adjust(method = "fdr")


df_sig <- df_temp %>% subset(padjSub < threshold) %>% droplevels
df_sig$Direction <- -1
df_sig$Direction[which(df_sig$log2FoldChange > 0)] <- 1
mids <- df_sig$Id %>% unique

#Subset the original matrix to obtain fold changes
df_lfc <- which(df_temp$Id %in% mids) %>%
  df_temp[.,] %>% droplevels

#Create matrices of enrichment and lfc
Tab_man <- acast(data = df_sig,formula = Id~Genotype,value.var = "Direction",fill = 0)
Tab_lfc <- acast(data = df_lfc,formula = Id~Genotype,value.var = "log2FoldChange")


#Cluster the asv patterns
mclust_asv <- hclust(d = (1-cor(Tab_lfc  %>%t)) %>% as.dist,method = "ward.D2") 
mclust_asv %>%plot

#Merge the df_clust structures with the structure used to plot
melted <- Tab_lfc %>% melt
colnames(melted) <- c("Id","Genotype","value")


#Order the genotype and asvs
order_asv <- mclust_asv$order %>% mclust_asv$labels[.]

melted$Genotype <- melted$Genotype %>% factor(levels = melted_asv$Genotype %>% levels)
melted$Id <- melted$Id %>% factor(levels = order_asv)


melted$CG <- match(melted$Genotype,melted_asv$Genotype) %>%
  melted_asv$CG[.]

melted$NumberGenotype <- match(melted$Genotype,melted_asv$Genotype) %>%
  melted_asv$NumberGenotype[.]


#Append significance 
df_sig$Significance <- "q < 0.1"
melted <- merge(melted,df_sig[,c("Genotype","Id","Significance")],by = c("Genotype","Id"), all.x = TRUE)


mat <- melted$Id %>% as.character %>% strsplit(split = "\\;") %>% unlist %>%
  matrix(data = .,byrow = T,ncol = 5)
nid <- paste0(mat[,3]," ",mat[,5])  %>% gsub(pattern = " p__",replacement = "p__") %>%
  gsub(pattern = "__",replacement = "_")
melted$NId <- nid

order_nids <- with(melted,order(Id)) %>%
  melted$NId[.] %>% unique
melted$NId <- melted$NId %>% factor(levels = order_nids)


### Plot 
p_heatmap <- ggplot(data = melted,aes(NId,NumberGenotype)) +
  geom_raster(aes(fill =value ), stat = "identity") +
  geom_tile(aes(color = Significance),fill = '#00000000', size = 1,width = 0.85,height = 0.85) + 
  facet_grid(CG~.,space = "free",scales = "free")  +
  scale_color_manual(values = c("black"))+
  scale_fill_paletteer_c("pals::kovesi.diverging_bwr_55_98_c37",
                         limits = c(-5,5),oob = squish,na.value = "white") +
  theme_ohchibi(size_panel_border = 2) +
  theme(
    strip.background.x = element_blank(),
    strip.background.y = element_blank(),
    strip.text.x = element_text(family = "Helvetica",face = "plain",size = 15),
    #strip.text.y = element_text(family = "Helvetica",face = "plain",size = 20),
    strip.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(family = "Helvetica",face = "plain",
                               size = 11,angle = 90,vjust = 0.5,hjust = 1),
    #axis.text.x = element_blank(),
    #axis.ticks.x = element_blank(),
    axis.text.y = element_text(family = "Helvetica",face = "plain",size = 18),
    legend.position = "none"
  ) +
  scale_x_discrete(expand = c(0,0))  +
  scale_y_discrete(expand = c(0,0))


#Prepare the dendrogram to put in the same graph
tree <- mclust_genotype %>% as.phylo
p_tree_genotype <- ggtree(tree,ladderize = F,size = 1.2) + 
  #The expand parth is fundamental to make the tree fit the composition with the other plots
  scale_y_continuous(expand = c(0.026,0.026))

melted_asv <- melted

#Create composition
composition_order <- egg::ggarrange(p_tree_genotype,p_heatmap,nrow = 1,widths = c(0.1,1))
oh.save.pdf(p = composition_order,outname = "revision_amplicon_heatmap_deseq2_order.pdf",
            outdir = "../figures/",width = 16,height = 12)
p_heatmap_order <- p_heatmap

##### Class level ########
level <- "Class"
df_temp <- df_deseq_sub %>% subset(Level == level) %>%
  droplevels

#Ad just pvalue
df_temp$padjSub <- df_temp$pvalue %>% p.adjust(method = "fdr")


df_sig <- df_temp %>% subset(padjSub < threshold) %>% droplevels
df_sig$Direction <- -1
df_sig$Direction[which(df_sig$log2FoldChange > 0)] <- 1
mids <- df_sig$Id %>% unique

#Subset the original matrix to obtain fold changes
df_lfc <- which(df_temp$Id %in% mids) %>%
  df_temp[.,] %>% droplevels

#Create matrices of enrichment and lfc
Tab_man <- acast(data = df_sig,formula = Id~Genotype,value.var = "Direction",fill = 0)
Tab_lfc <- acast(data = df_lfc,formula = Id~Genotype,value.var = "log2FoldChange")


#Cluster the asv patterns
mclust_asv <- hclust(d = (1-cor(Tab_lfc  %>%t)) %>% as.dist,method = "ward.D2") 
mclust_asv %>%plot

#Merge the df_clust structures with the structure used to plot
melted <- Tab_lfc %>% melt
colnames(melted) <- c("Id","Genotype","value")


#Order the genotype and asvs
order_asv <- mclust_asv$order %>% mclust_asv$labels[.]

melted$Genotype <- melted$Genotype %>% factor(levels = melted_asv$Genotype %>% levels)
melted$Id <- melted$Id %>% factor(levels = order_asv)


melted$CG <- match(melted$Genotype,melted_asv$Genotype) %>%
  melted_asv$CG[.]

melted$NumberGenotype <- match(melted$Genotype,melted_asv$Genotype) %>%
  melted_asv$NumberGenotype[.]


#Append significance 
df_sig$Significance <- "q < 0.1"
melted <- merge(melted,df_sig[,c("Genotype","Id","Significance")],by = c("Genotype","Id"), all.x = TRUE)

mat <- melted$Id %>% as.character %>% strsplit(split = "\\;") %>% unlist %>%
  matrix(data = .,byrow = T,ncol = 4)
nid <- paste0(mat[,3]," ",mat[,4])  %>% gsub(pattern = " p__",replacement = "p__") %>%
  gsub(pattern = "__",replacement = "_")
melted$NId <- nid

order_nids <- with(melted,order(Id)) %>%
  melted$NId[.] %>% unique
melted$NId <- melted$NId %>% factor(levels = order_nids)


### Plot 
p_heatmap <- ggplot(data = melted,aes(NId,NumberGenotype)) +
  geom_raster(aes(fill =value ), stat = "identity") +
  geom_tile(aes(color = Significance),fill = '#00000000', size = 1,width = 0.85,height = 0.85) + 
  facet_grid(CG~.,space = "free",scales = "free")  +
  scale_color_manual(values = c("black"))+
  scale_fill_paletteer_c("pals::kovesi.diverging_bwr_55_98_c37",
                         limits = c(-5,5),oob = squish,na.value = "white") +
  theme_ohchibi(size_panel_border = 2) +
  theme(
    strip.background.x = element_blank(),
    strip.background.y = element_blank(),
    strip.text.x = element_text(family = "Helvetica",face = "plain",size = 15),
    #strip.text.y = element_text(family = "Helvetica",face = "plain",size = 20),
    strip.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(family = "Helvetica",face = "plain",
                               size = 11,angle = 90,vjust = 0.5,hjust = 1),
    #axis.text.x = element_blank(),
    #axis.ticks.x = element_blank(),
    axis.text.y = element_text(family = "Helvetica",face = "plain",size = 18),
    legend.position = "none"
  ) +
  scale_x_discrete(expand = c(0,0))  +
  scale_y_discrete(expand = c(0,0))


#Prepare the dendrogram to put in the same graph
tree <- mclust_genotype %>% as.phylo
p_tree_genotype <- ggtree(tree,ladderize = F,size = 1.2) + 
  #The expand parth is fundamental to make the tree fit the composition with the other plots
  scale_y_continuous(expand = c(0.026,0.026))

melted_asv <- melted

#Create composition
composition_class <- egg::ggarrange(p_tree_genotype,p_heatmap,nrow = 1,widths = c(0.1,1))
oh.save.pdf(p = composition_class,outname = "revision_amplicon_heatmap_deseq2_class.pdf",
            outdir = "../figures/",width = 16,height = 12)
p_heatmap_class <- p_heatmap

##### Phylum level ########
level <- "Phylum"
df_temp <- df_deseq_sub %>% subset(Level == level) %>%
  droplevels

#Ad just pvalue
df_temp$padjSub <- df_temp$pvalue %>% p.adjust(method = "fdr")


df_sig <- df_temp %>% subset(padjSub < threshold) %>% droplevels
df_sig$Direction <- -1
df_sig$Direction[which(df_sig$log2FoldChange > 0)] <- 1
mids <- df_sig$Id %>% unique

#Subset the original matrix to obtain fold changes
df_lfc <- which(df_temp$Id %in% mids) %>%
  df_temp[.,] %>% droplevels

#Create matrices of enrichment and lfc
Tab_man <- acast(data = df_sig,formula = Id~Genotype,value.var = "Direction",fill = 0)
Tab_lfc <- acast(data = df_lfc,formula = Id~Genotype,value.var = "log2FoldChange")


#Cluster the asv patterns
mclust_asv <- hclust(d = (1-cor(Tab_lfc  %>%t)) %>% as.dist,method = "ward.D2") 
mclust_asv %>%plot

#Merge the df_clust structures with the structure used to plot
melted <- Tab_lfc %>% melt
colnames(melted) <- c("Id","Genotype","value")


#Order the genotype and asvs
order_asv <- mclust_asv$order %>% mclust_asv$labels[.]

melted$Genotype <- melted$Genotype %>% factor(levels = melted_asv$Genotype %>% levels)
melted$Id <- melted$Id %>% factor(levels = order_asv)


melted$CG <- match(melted$Genotype,melted_asv$Genotype) %>%
  melted_asv$CG[.]

melted$NumberGenotype <- match(melted$Genotype,melted_asv$Genotype) %>%
  melted_asv$NumberGenotype[.]


#Append significance 
df_sig$Significance <- "q < 0.1"
melted <- merge(melted,df_sig[,c("Genotype","Id","Significance")],by = c("Genotype","Id"), all.x = TRUE)


nid <- melted$Id %>% gsub(pattern = ".* p__",replacement = "p_") 
melted$NId <- nid

order_nids <- with(melted,order(Id)) %>%
  melted$NId[.] %>% unique
melted$NId <- melted$NId %>% factor(levels = order_nids)


### Plot 
p_heatmap <- ggplot(data = melted,aes(NId,NumberGenotype)) +
  geom_raster(aes(fill =value ), stat = "identity") +
  geom_tile(aes(color = Significance),fill = '#00000000', size = 1,width = 0.85,height = 0.85) + 
  facet_grid(CG~.,space = "free",scales = "free")  +
  scale_color_manual(values = c("black"))+
  scale_fill_paletteer_c("pals::kovesi.diverging_bwr_55_98_c37",
                         limits = c(-5,5),oob = squish,na.value = "white") +
  theme_ohchibi(size_panel_border = 2) +
  theme(
    strip.background.x = element_blank(),
    strip.background.y = element_blank(),
    strip.text.x = element_text(family = "Helvetica",face = "plain",size = 15),
    #strip.text.y = element_text(family = "Helvetica",face = "plain",size = 20),
    strip.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(family = "Helvetica",face = "plain",
                               size = 11,angle = 90,vjust = 0.5,hjust = 1),
    #axis.text.x = element_blank(),
    #axis.ticks.x = element_blank(),
    axis.text.y = element_text(family = "Helvetica",face = "plain",size = 18),
    legend.position = "none"
  ) +
  scale_x_discrete(expand = c(0,0))  +
  scale_y_discrete(expand = c(0,0))


#Prepare the dendrogram to put in the same graph
tree <- mclust_genotype %>% as.phylo
p_tree_genotype <- ggtree(tree,ladderize = F,size = 1.2) + 
  #The expand parth is fundamental to make the tree fit the composition with the other plots
  scale_y_continuous(expand = c(0.026,0.026))

melted_asv <- melted




### Create composite figure
all_composition <- egg::ggarrange(p_tree_genotype,p_heatmap_asv,p_heatmap_genus,p_heatmap_family,p_heatmap_order,p_heatmap_class,
                                  nrow = 1,widths = c(0.1,1,1,1,1,1))
oh.save.pdf(p = all_composition,outname = "revision_amplicon_heatmap_deseq2_root_all.pdf",
            outdir = "../figures/",width = 48,height = 12)

### Create composite figure
all_composition <- egg::ggarrange(p_tree_genotype,p_heatmap_asv,p_heatmap_family,p_heatmap_class,
                                  nrow = 1,widths = c(0.1,0.8,1,0.5))
oh.save.pdf(p = all_composition,outname = "revision_amplicon_heatmap_deseq2_root_chosen.pdf",
            outdir = "../figures/",width = 29,height = 12)

#rm(list=ls())
#dev.off()
#gc()
