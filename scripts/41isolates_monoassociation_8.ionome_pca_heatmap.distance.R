library(ohchibi)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(harrietr)
library(egg)
library(paletteer)
library(ggtree)
library(scales)
library(car)
library(Rmisc)

set.seed(130816)
setwd('/home/isai/Documents/results/rootbarriersmicro/scripts/')


Dat <- readRDS("../cleandata/dat_41isolatesmonoassociation_cfu_suberin_primroot_dryweight_ionome.RDS")

##### Ionome #######
df_ionome <- Dat$Ionome %>% subset(Strains != "CDEF") %>% droplevels

#Aggregate the dataframe to display as heatmap
melted <- df_ionome  %>%
  acast(formula = Index~Ion,fill = 0,
        value.var = "NormValue") %>%
  scale 

#Remove outliers using quantile
maxvalue <- quantile(melted,0.99)
minvalue <- quantile(melted,0.01)
melted[melted > maxvalue] <- maxvalue
melted[melted < minvalue] <- minvalue

#Perform pca
mpca <- prcomp(x = melted,center = F,scale. = F)
summary(mpca)
scores <- mpca$x %>% as.data.frame
scores$Index <- rownames(scores)

Map <- df_ionome[,c(1,4,5,6,7,8,9,11)] %>% unique
scores <- merge(scores,Map, by = "Index")


dfpc1 <- summarySE(data = scores,measurevar = "PC1",groupvars = "Strains")
dfpc2 <- summarySE(data = scores,measurevar = "PC2",groupvars = "Strains")


merged <- merge(dfpc1,dfpc2,by = "Strains")
merged <- merged[,c(1,3,6,8,11)]
colnames(merged) <- c("Strains","PC1","ci1","PC2","ci2")


Map <- df_ionome[,c(1,5,11)] %>% unique %>% droplevels

merged$Type <- match(merged$Strains,Map$Strains) %>% Map$Type[.]
summary(mpca)


paleta_alive <- c("#FFB919","#00CC7A","#C200FF")

#load the groups structure
df_groups <- read.table(file = "../cleandata/df_41_suberin_results.distance.tsv",header = T,sep = "\t")
merged <- merge(merged,df_groups, by = "Strains",all.x = T)
merged$Group <- merged$Group %>% as.character
merged$Group[which(is.na(merged$Group))] <- "HK"
merged$Group[merged$Strains %>% grep(pattern = "NB")] <- "NB"

merged$Group <- merged$Group %>% 
  factor(levels = c("NB","HK","Group1","Group2","Group3","Group4","Group5"))



paleta_alive <- c("#414141","#D9D9D9",paletteer_d(package = "LaCroixColoR",palette = "PassionFruit",n = 5))

p <- ggplot(data = merged,aes(PC1,PC2)) +
  geom_vline(xintercept = 0,size = 2,color = "#D9D9D9",linetype = "longdash")+
  geom_hline(yintercept = 0,size = 2,color = "#D9D9D9",linetype = "longdash")+
  #stat_ellipse(aes(group = Type,color = Type),size = 1)+
  geom_errorbarh(mapping = aes(xmin =PC1 - ci1,xmax = PC1 + ci1, ),size = 0.01,alpha = 0.1) +
  geom_errorbar(mapping = aes(ymin =PC2 - ci2,ymax = PC2 + ci2,),size = 0.01,alpha = 0.1)  + 
  geom_point(shape = 21,size = 10,aes(fill = Group),stroke = 1) + 
  theme_ohchibi(size_panel_border = 2) + 
  xlab(label = "PC1 28.22%") + ylab(label = "PC2 19.20%") +
  scale_fill_manual(values = paleta_alive) +
  scale_color_manual(values = paleta_alive)


oh.save.pdf(p = p,outname = "pca_ionome_41_hk_alive.distance.pdf",outdir = "../figures/",width = 10,height = 8)

##### Test differences inside each ion against NB
Res <- NULL
for(ion in df_ionome$Ion %>% unique){
  temp <- df_ionome %>% subset(Ion == ion) %>% droplevels
  temp$Strains <- temp$Strains %>% factor %>% relevel(ref = "NB")
  temp$UId <- paste(temp$Batch,temp$Replicates,sep = "_")
  temp_nb <- temp %>% subset(Strains == "NB") %>% droplevels
  temp_other <- temp %>% subset(Strains != "NB") %>% droplevels
  for(strain in temp_other$Strains %>% unique){
    df_subject <- temp_other %>% subset(Strains == strain) %>% droplevels
    exps <- df_subject$Batch %>% unique
    temp_nb_inner <- temp_nb[which(temp_nb$Batch %in% exps),]
    #Levene test
    m1 <- rbind(df_subject,temp_nb_inner) %>% 
      leveneTest(NormValue ~Strains,data = .) 
    pval <- m1$`Pr(>F)`[1]
    if(pval < 0.05){
      mtest <- t.test(x = temp_nb_inner$NormValue,df_subject$NormValue,
                      paired = F,exact = F,var.equal = F)  
    }else{
      mtest <- t.test(x = temp_nb_inner$NormValue,df_subject$NormValue,
                      paired = F,exact = F,var.equal = T)  
    }
    y <- df_subject$NormValue %>% mean
    x <- temp_nb_inner$NormValue %>% mean
    mfold <- y/x
    mestimate <- (df_subject$NormValue %>% mean)-(temp_nb_inner$NormValue %>%mean) 
    #This estimate is equal to an estimate calculated by anova given there are no other predictos
    mtest <- data.frame(Ion = ion, Strains = strain,
                        p.value = mtest$p.value, estimate = mestimate,fold = mfold)
    Res <- rbind(Res,mtest)
    
  }
}

#Adjust for multiple testing
Res$p.adj <- Res$p.value %>% p.adjust(method = "fdr")
Res$p.adj %>% hist


#Create a matrix to cluster based on the z score 
Tab <-  df_ionome %>% acast(formula = Index~Ion,fill = 0,
                            value.var = "NormValue") %>%
  scale 
melted <- Tab
#Remove outliers using quantile
maxvalue <- quantile(melted,0.99)
minvalue <- quantile(melted,0.01)
melted[melted > maxvalue] <- maxvalue
melted[melted < minvalue] <- minvalue

melted <- melted %>% melt

Map <- df_ionome[,c(1,5,11)] %>% unique %>% droplevels
colnames(melted) <- c("Index","Ion","zscore")
melted <- merge(melted,Map, by = "Index")                    

Tab <- acast(data = melted,formula = Strains~Ion,fun.aggregate = mean,
             fill = 0,value.var = "zscore")


#Cluster the ions and strains 
mclust_strains <- dist(x = Tab) %>% hclust(method = "ward.D")
#mclust_strains <- hclust(d = 1-cor(Tab %>% t) %>% as.dist,method = "complete")
mclust_ions <- hclust(d = 1-cor(Tab) %>% as.dist,method = "ward.D")
#mclust_ions <- dist(x = Tab %>% t) %>% hclust(method ="complete")
order_strains <- mclust_strains$order %>% mclust_strains$labels[.]
order_ions <- mclust_ions$order %>% mclust_ions$labels[.]

#Divide in clusters the ions and strains
df_clust_ions <- mclust_ions %>% cutree(k = 4) %>% 
  data.frame(Ion = names(.),row.names = NULL, ClusterIon = paste0("CI",.))
df_clust_ions <- df_clust_ions[,-1]

df_clust_strains <- mclust_strains %>% cutree(k = 3)%>% 
  data.frame(Strains = names(.),row.names = NULL, ClusterStrains = paste0("CS",.))
df_clust_strains <- df_clust_strains[,-1]


melted <- Tab %>% melt
colnames(melted) <- c("Strains","Ion","zscore")

#Merge with the clusters of ion
melted <- merge(melted,df_clust_ions, by = "Ion")
melted <- merge(melted,df_clust_strains, by = "Strains")


#Load screening barriers dat
Dat <- readRDS(file = "../cleandata/dat_screeningbarriers.RDS")
meta <- Dat$Map[,c(16,7)]
colnames(meta)[1] <- "UStrains"
melted$UStrains <- melted$Strains %>% gsub(pattern = "HK",replacement = "")
melted <- merge(melted,meta, by = "UStrains", all.x = TRUE)
melted$Nom <- paste0(melted$Genus," ",melted$Strains)
melted$Nom <- melted$Nom %>% gsub(pattern = "NA ",replacement = "")



melted$UId <- paste(melted$Strains,melted$Ion,sep = "-")

Res$UId <- paste(Res$Strains,Res$Ion,sep = "-")
melted <- merge(melted,Res[,-c(1,2)],by = "UId",all.x = TRUE)
melted$Significance <- rep("NoSignificant")
melted$Significance[which(melted$p.adj < 0.05)] <- "Significant"

melted$Strains <- melted$Strains %>% factor(levels = order_strains)
melted$Ion <- melted$Ion %>% factor(levels = order_ions)

#Determine the order of the clusters of ions
melted <- with(melted,order(Ion)) %>% melted[.,]
order_cluster_ion <- melted$ClusterIon %>% unique %>% as.character
melted$ClusterIon <- melted$ClusterIon %>% factor(levels = order_cluster_ion %>% rev)

#Determine the order of the strain clsuters
melted <- with(melted,order(Strains)) %>% melted[.,]
order_cluster_strains <- melted$ClusterStrains %>% unique %>% as.character
melted$ClusterStrains <- melted$ClusterStrains %>% factor(levels = order_cluster_strains )

order_nom <- melted$Nom %>% unique
melted$Nom <- melted$Nom %>% factor(levels = order_nom)


p_heatmap <- ggplot(data = melted,mapping = aes(x = Ion, y = Nom)) + 
  facet_grid(ClusterIon~ClusterStrains,scales = "free",space = "free") +
  geom_raster(aes(fill = zscore)) +  
  geom_tile(aes(color = Significance),fill = '#00000000', size = 0.5,width = 0.85,height = 0.85) + 
  scale_fill_paletteer_c(package = "pals",palette = "kovesi.diverging_bwr_55_98_c37",limits = c(-1.5,1.5),oob = squish) +
  scale_color_manual(values = c("#00000000","black"))+
  theme_ohchibi(size_axis_text.x = 10,
                angle_text.x = 90,
                size_axis_text.y = 10,
                size_axis_title.x = 0,
                size_axis_title.y = 0,
                legend_proportion_size = 0,
                size_title_text = 0,size_legend_text = 0,
                size_panel_border = 1,
                size_lines_panel = 0) +
  theme(aspect.ratio = 0,axis.ticks = element_blank(),
        panel.spacing = unit(0.075, "lines"),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text.y = element_text(family = "Arial",face = "bold",hjust = 1,vjust = 0.5),
        axis.text.x = element_text(hjust = 1,vjust = 0.5,family = "Arial",face = "bold")
  ) +
  coord_flip(expand = F)




#Create a bar for the type of each strain treatment
df_un <- melted[,c(3,7)] %>% unique
df_un$Bar <- rep("Bar")
df_un <- merge(df_un,df_ionome[,c(5,11)], by = "Strains")

df_un <- merge(df_un,df_groups, by = "Strains",all.x = T) 
df_un$Group <- df_un$Group %>% as.character
df_un$Group[which(is.na(df_un$Group))] <- "HK"
df_un$Group[df_un$Strains %>% grep(pattern = "NB")] <- "NB"

df_un$Group <- df_un$Group %>% 
  factor(levels = c("NB","HK","Group1","Group2","Group3","Group4","Group5"))




p_strip <- ggplot(data = df_un,aes(Strains,Bar)) +
  facet_grid(.~ClusterStrains,scales = "free",space = "free") +
  geom_tile(aes(fill = Group)) + 
  theme_ohchibi(size_panel_border = 1) +
  scale_y_discrete(expand = c(0,0)) +
  theme(
    aspect.ratio = 0,axis.ticks = element_blank(),
    panel.spacing = unit(0.075, "lines"),
    legend.position = "none",
    strip.background = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    strip.text = element_blank()
  )+
  scale_fill_manual(values = paleta_alive)



tree <- mclust_strains %>% as.phylo
p_tree_strains  <- ggtree(tree,ladderize = F,size = 0.7) + 
  #The expand parth is fundamental to make the tree fit the composition with the other plots
  scale_y_continuous(expand = c(0.009,0.009)) +
  coord_flip()  + scale_x_reverse()

tree <- mclust_ions %>% as.phylo
p_tree_ions <- ggtree(tree,ladderize = F,size = 0.7) + 
  #The expand parth is fundamental to make the tree fit the composition with the other plots
  scale_y_continuous(expand = c(0.026,0.026))

#Create blank plot
p_blank <- ggplot() +theme_void()



#Create composition
composition <- egg::ggarrange(p_blank,p_tree_strains,
                              p_blank,p_strip,
                              p_tree_ions,p_heatmap,
                              ncol = 2,nrow = 3,byrow = T,
                              heights = c(0.25,0.05,1),
                              widths = c(0.05,1),padding = unit(20,"line"))


oh.save.pdf(p = composition,outname = "heatmap_composition_ionome_nb__hk_alive.distance.pdf",
            outdir = "../figures/",height = 5.5,width = 12)


####### Suberin heatmap
melted_sub <- melted$Strains %>% grep(pattern = "HK",invert = T) %>%
  melted[.,] %>% droplevels


df_suberin <- read.table(file = "../cleandata/df_41_suberin_results.distance.tsv",header = T)
order_strains <- df_suberin$Strains %>% as.character
melted_sub$Strains <- melted_sub$Strains %>% factor(levels = order_strains)

order_nom <- with(melted_sub,order(Strains)) %>%
  melted_sub$Nom[.] %>% as.character %>% unique

melted_sub$Nom <- melted_sub$Nom %>% factor(levels = order_nom)

p_heatmap <- ggplot(data = melted_sub,mapping = aes(x = Ion, y = Nom)) + 
  facet_grid(ClusterIon~.,scales = "free",space = "free") +
  geom_raster(aes(fill = zscore)) +  
  geom_tile(aes(color = Significance),fill = '#00000000', size = 0.5,width = 0.85,height = 0.85) + 
  scale_fill_paletteer_c(package = "pals",palette = "kovesi.diverging_bwr_55_98_c37",limits = c(-1.5,1.5),oob = squish) +
  scale_color_manual(values = c("#00000000","black"))+
  theme_ohchibi(size_axis_text.x = 10,
                angle_text.x = 90,
                size_axis_text.y = 10,
                size_axis_title.x = 0,
                size_axis_title.y = 0,
                #legend_proportion_size = 0,
                size_title_text = 0,size_legend_text = 0,
                size_panel_border = 1,
                size_lines_panel = 0) +
  theme(aspect.ratio = 0,axis.ticks = element_blank(),
        panel.spacing = unit(0.075, "lines"),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text.y = element_text(family = "Arial",face = "bold",hjust = 1,vjust = 0.5),
        axis.text.x = element_text(hjust = 1,vjust = 0.5,family = "Arial",face = "bold")
  ) +
  coord_flip(expand = F)

oh.save.pdf(p = p_heatmap,outname = "heatmap_suberin_ordered_ionome.distance.pdf",
            outdir = "../figures/",height = 5.5,width = 12)


### heatmap ###
p_heatmap <- ggplot(data = melted_sub,mapping = aes(x = Ion, y = Nom)) + 
  facet_grid(ClusterIon~.,scales = "free",space = "free") +
  geom_raster(aes(fill = zscore)) +  
  geom_tile(aes(color = Significance),fill = '#00000000', size = 0.5,width = 0.85,height = 0.85) + 
  scale_fill_paletteer_c(package = "pals",palette = "kovesi.diverging_bwr_55_98_c37",limits = c(-1.5,1.5),oob = squish) +
  scale_color_manual(values = c("#00000000","black"))+
  theme_ohchibi(size_axis_text.x = 10,
                angle_text.x = 90,
                size_axis_text.y = 10,
                size_axis_title.x = 0,
                size_axis_title.y = 0,
                #legend_proportion_size = 0,
                size_title_text = 0,
                size_panel_border = 1,
                size_lines_panel = 0) +
  theme(aspect.ratio = 0,axis.ticks = element_blank(),
        panel.spacing = unit(0.075, "lines"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text.y = element_text(family = "Arial",face = "bold",hjust = 1,vjust = 0.5),
        axis.text.x = element_text(hjust = 1,vjust = 0.5,family = "Arial",face = "bold")
  ) +
  coord_flip(expand = F)

oh.save.pdf(p = p_heatmap,outname = "heatmap_suberin_ordered_ionome.distance.legend.pdf",
            outdir = "../figures/",height = 5.5,width = 12)




rm(list=ls())
dev.off()