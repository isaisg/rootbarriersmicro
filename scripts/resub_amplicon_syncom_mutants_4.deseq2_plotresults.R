library(ohchibi)
library(DESeq2)
library(scales)
library(paletteer)

set.seed(130816)
setwd('/Users/isaisalasgonzalez/Documents/rootbarriersmicro/scripts/')

Dat <- readRDS(file = "../cleandata/deseq2_models_rootbarriers_amplicon_bacteria_syncom_mutants.RDS")

###### Root #########
dds <- Dat$RootMutants

#Everything against Full
genos <- dds$Genotype %>% as.character %>% unique %>% 
  grep(pattern = "Col-0",invert = T,value = T)

Res <- NULL
for(geno in genos){
  cat("Working on ",geno,"\n")
  cont <- c("Genotype",geno,"Col-0")
  res <- dds %>% results(contrast = cont) %>%
    as.data.frame
  #res <- match(usable_asvs,rownames(res)) %>%
  #  res[.,]
  res$Contrast <- cont[2:3]%>%paste0(collapse = "_")
  res$Id <- rownames(res)
  rownames(res) <- NULL
  Res <- rbind(Res,res)
}


Res$Genotype <- Res$Contrast %>% gsub(pattern = "_.*",replacement = "")



Res$Significance <- "NoSignificant"
Res$Significance[which(Res$padj < 0.05)] <- "Significant"

#Adjust the Stresses
Res$Genotype <- Res$Genotype %>% 
  factor()

## Split the ones that have more than one id for visualization purpose
Res_un <- Res$Id %>% grep(pattern = "_",invert = T) %>% 
  Res[.,] %>% droplevels

Res_do <- Res$Id %>% grep(pattern = "_",invert = F) %>% 
  Res[.,] %>% droplevels


#Loop over creating duplicate rows
Res_dup <- NULL
for(i in 1:nrow(Res_do)){
  temp <- Res_do[i,]
  mst <- temp$Id %>% strsplit(split = "_") %>% unlist
  for(st in mst){
    mend <- temp
    mend$Id <- st
    Res_dup <- rbind(Res_dup,mend)
  }
}


Res <- rbind(Res_un,Res_dup)


df_sub <- read.table(file = "../cleandata/df_41_suberin_results.distance.tsv",header = T,sep = "\t")
Res$Group <- match(Res$Id,df_sub$Strains) %>%
  df_sub$Group[.] %>% as.character

#Order the strains based on the suberin pattern
df_sub <- with(df_sub,order(Norm_Distance_tip_cz)) %>%
  df_sub[.,]
order_id <- df_sub$Strains %>% as.character %>% unique


Res$Id <- Res$Id %>% factor(levels = order_id)
Res$Fraction <- "Root"

##Create a matrix based on the enrichment patterns
Res$Direction <- 0
Res$Direction[which((Res$padj < 0.1) & (Res$log2FoldChange > 0))] <- 1
Res$Direction[which((Res$padj < 0.1) & (Res$log2FoldChange < 0))] <- -1

#Determine clustering and direction
Tab_man <- acast(data = Res,formula = Id~Genotype,value.var = "Direction")
dist_man <- dist(x = Tab_man %>% t,method = "manhattan")
mclust_man <- hclust(d = dist_man,method = "complete")
order_genotype <- mclust_man $order %>% mclust_man$labels[.]

Res$Genotype <- Res$Genotype %>% factor(levels = order_genotype)

dist_man <- dist(x = Tab_man,method = "manhattan")
mclust_man <- hclust(d = dist_man,method = "ward.D2")
order_ids <- mclust_man $order %>% mclust_man$labels[.]

Res$Id <- Res$Id %>% factor(levels = c(order_ids))

#Inverted version
p_root <- ggplot(data = Res,aes(Id,Genotype)) +
  geom_raster(aes(fill =log2FoldChange ), stat = "identity") +
  geom_tile(aes(color = Significance),fill = '#00000000', size = 1,width = 0.85,height = 0.85) + 
  scale_color_manual(values = c("#00000000","black"))+
  facet_grid(Fraction~.,space = "free",scales = "free")  +
  scale_fill_paletteer_c("pals::kovesi.diverging_bwr_55_98_c37",
                         limits = c(-5,5),oob = squish,na.value = "white") +
  theme_ohchibi(size_panel_border = 2) +
  theme(
    strip.background.x = element_blank(),
    strip.background.y = element_blank(),
    strip.text.x = element_text(family = "Helvetica",face = "plain",size = 15),
    strip.text.y = element_text(family = "Helvetica",face = "plain",size = 20),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(family = "Helvetica",face = "plain",
                               size = 22,angle = 90,vjust = 0.5,hjust = 1),
    axis.text.y = element_text(family = "Helvetica",face = "plain",size = 18),
    legend.position = "none"
    
  ) +
  scale_x_discrete(expand = c(0,0))  +
  scale_y_discrete(expand = c(0,0))
Res_root <- Res

###### Shoot #####
dds <- Dat$ShootMutants

#Everything against Full
genos <- dds$Genotype %>% as.character %>% unique %>% 
  grep(pattern = "Col-0",invert = T,value = T)

Res <- NULL
for(geno in genos){
  cat("Working on ",geno,"\n")
  cont <- c("Genotype",geno,"Col-0")
  res <- dds %>% results(contrast = cont) %>%
    as.data.frame
  #res <- match(usable_asvs,rownames(res)) %>%
  #  res[.,]
  res$Contrast <- cont[2:3]%>%paste0(collapse = "_")
  res$Id <- rownames(res)
  rownames(res) <- NULL
  Res <- rbind(Res,res)
}


Res$Genotype <- Res$Contrast %>% gsub(pattern = "_.*",replacement = "")



Res$Significance <- "NoSignificant"
Res$Significance[which(Res$padj < 0.05)] <- "Significant"

#Adjust the Stresses
Res$Genotype <- Res$Genotype %>% 
  factor()

## Split the ones that have more than one id for visualization purpose
Res_un <- Res$Id %>% grep(pattern = "_",invert = T) %>% 
  Res[.,] %>% droplevels

Res_do <- Res$Id %>% grep(pattern = "_",invert = F) %>% 
  Res[.,] %>% droplevels


#Loop over creating duplicate rows
Res_dup <- NULL
for(i in 1:nrow(Res_do)){
  temp <- Res_do[i,]
  mst <- temp$Id %>% strsplit(split = "_") %>% unlist
  for(st in mst){
    mend <- temp
    mend$Id <- st
    Res_dup <- rbind(Res_dup,mend)
  }
}


Res <- rbind(Res_un,Res_dup)


df_sub <- read.table(file = "../cleandata/df_41_suberin_results.distance.tsv",header = T,sep = "\t")
Res$Group <- match(Res$Id,df_sub$Strains) %>%
  df_sub$Group[.] %>% as.character

#Order the strains based on the suberin pattern
df_sub <- with(df_sub,order(Norm_Distance_tip_cz)) %>%
  df_sub[.,]
order_id <- df_sub$Strains %>% as.character %>% unique


Res$Id <- Res$Id %>% factor(levels = order_id)
Res$Fraction <- "Shoot"

##Create a matrix based on the enrichment patterns
Res$Direction <- 0
Res$Direction[which((Res$padj < 0.1) & (Res$log2FoldChange > 0))] <- 1
Res$Direction[which((Res$padj < 0.1) & (Res$log2FoldChange < 0))] <- -1
Res$Genotype <- Res$Genotype %>% factor(levels = order_genotype)

#Inverted version
p_shoot <- ggplot(data = Res,aes(Id,Genotype)) +
  geom_raster(aes(fill =log2FoldChange ), stat = "identity") +
  geom_tile(aes(color = Significance),fill = '#00000000', size = 1,width = 0.85,height = 0.85) + 
  scale_color_manual(values = c("#00000000","black"))+
  facet_grid(Fraction~Group,space = "free",scales = "free")  +
  scale_fill_paletteer_c("pals::kovesi.diverging_bwr_55_98_c37",
                         limits = c(-5,5),oob = squish,na.value = "white") +
  theme_ohchibi(size_panel_border = 2) +
  theme(
    strip.background.x = element_blank(),
    strip.background.y = element_blank(),
    strip.text.x = element_text(family = "Helvetica",face = "plain",size = 15),
    strip.text.y = element_text(family = "Helvetica",face = "plain",size = 20),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(family = "Helvetica",face = "plain",
                               size = 22,angle = 90,vjust = 0.5,hjust = 1),
    axis.text.y = element_text(family = "Helvetica",face = "plain",size = 18),
    legend.position = "none"
    
  ) +
  scale_x_discrete(expand = c(0,0))  +
  scale_y_discrete(expand = c(0,0))

Res_shoot <- Res

###### Agar #####
dds <- Dat$AgarMutants

#Everything against Full
genos <- dds$Genotype %>% as.character %>% unique %>% 
  grep(pattern = "Col-0",invert = T,value = T)

Res <- NULL
for(geno in genos){
  cat("Working on ",geno,"\n")
  cont <- c("Genotype",geno,"Col-0")
  res <- dds %>% results(contrast = cont) %>%
    as.data.frame
  #res <- match(usable_asvs,rownames(res)) %>%
  #  res[.,]
  res$Contrast <- cont[2:3]%>%paste0(collapse = "_")
  res$Id <- rownames(res)
  rownames(res) <- NULL
  Res <- rbind(Res,res)
}


Res$Genotype <- Res$Contrast %>% gsub(pattern = "_.*",replacement = "")



Res$Significance <- "NoSignificant"
Res$Significance[which(Res$padj < 0.05)] <- "Significant"

#Adjust the Stresses
Res$Genotype <- Res$Genotype %>% 
  factor()

## Split the ones that have more than one id for visualization purpose
Res_un <- Res$Id %>% grep(pattern = "_",invert = T) %>% 
  Res[.,] %>% droplevels

Res_do <- Res$Id %>% grep(pattern = "_",invert = F) %>% 
  Res[.,] %>% droplevels


#Loop over creating duplicate rows
Res_dup <- NULL
for(i in 1:nrow(Res_do)){
  temp <- Res_do[i,]
  mst <- temp$Id %>% strsplit(split = "_") %>% unlist
  for(st in mst){
    mend <- temp
    mend$Id <- st
    Res_dup <- rbind(Res_dup,mend)
  }
}


Res <- rbind(Res_un,Res_dup)


df_sub <- read.table(file = "../cleandata/df_41_suberin_results.distance.tsv",header = T,sep = "\t")
Res$Group <- match(Res$Id,df_sub$Strains) %>%
  df_sub$Group[.] %>% as.character

#Order the strains based on the suberin pattern
df_sub <- with(df_sub,order(Norm_Distance_tip_cz)) %>%
  df_sub[.,]
order_id <- df_sub$Strains %>% as.character %>% unique


Res$Id <- Res$Id %>% factor(levels = order_id)
Res$Fraction <- "Agar"

##Create a matrix based on the enrichment patterns
Res$Direction <- 0
Res$Direction[which((Res$padj < 0.1) & (Res$log2FoldChange > 0))] <- 1
Res$Direction[which((Res$padj < 0.1) & (Res$log2FoldChange < 0))] <- -1
Res$Genotype <- Res$Genotype %>% factor(levels = order_genotype)

#Inverted version
p_agar <- ggplot(data = Res,aes(Id,Genotype)) +
  geom_raster(aes(fill =log2FoldChange ), stat = "identity") +
  geom_tile(aes(color = Significance),fill = '#00000000', size = 1,width = 0.85,height = 0.85) + 
  scale_color_manual(values = c("#00000000","black"))+
  facet_grid(Fraction~Group,space = "free",scales = "free")  +
  scale_fill_paletteer_c("pals::kovesi.diverging_bwr_55_98_c37",
                         limits = c(-5,5),oob = squish,na.value = "white") +
  theme_ohchibi(size_panel_border = 2) +
  theme(
    strip.background.x = element_blank(),
    strip.background.y = element_blank(),
    strip.text.x = element_text(family = "Helvetica",face = "plain",size = 15),
    strip.text.y = element_text(family = "Helvetica",face = "plain",size = 20),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(family = "Helvetica",face = "plain",
                               size = 22,angle = 90,vjust = 0.5,hjust = 1),
    axis.text.y = element_text(family = "Helvetica",face = "plain",size = 18),
    legend.position = "none"
    
  ) +
  scale_x_discrete(expand = c(0,0))  +
  scale_y_discrete(expand = c(0,0))

Res_agar <- Res

Res <- rbind(Res_root,Res_shoot,Res_agar)

Res$Fraction <- Res$Fraction %>%
  factor(levels = c("Root","Shoot","Agar"))

Res$Id <- Res$Id %>% factor(levels = c(order_ids,"L288"))
p <- ggplot(data = Res,aes(Id,Genotype)) +
  geom_raster(aes(fill =log2FoldChange ), stat = "identity") +
  geom_tile(aes(color = Significance),fill = '#00000000', size = 1,width = 0.85,height = 0.85) + 
  scale_color_manual(values = c("#00000000","black"))+
  facet_grid(Fraction~.,space = "free",scales = "free")  +
  scale_fill_paletteer_c("pals::kovesi.diverging_bwr_55_98_c37",
                         limits = c(-5,5),oob = squish,na.value = "white") +
  theme_ohchibi(size_panel_border = 2) +
  theme(
    strip.background.x = element_blank(),
    strip.background.y = element_blank(),
    strip.text.x = element_text(family = "Arial",face = "bold",size = 15),
    strip.text.y = element_text(family = "Arial",face = "bold",size = 20),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(family = "Arial",face = "bold",
                               size = 22,angle = 90,vjust = 0.5,hjust = 1),
    axis.text.y = element_text(family = "Arial",face = "bold",size = 18)
    
  ) +
  scale_x_discrete(expand = c(0,0)) 

oh.save.pdf(p = p,outname = "resub_clustered_heatmap_genotype_syncom_by_fraction.pdf",
            outdir = "../figures/",width = 18,height =8)


### Clustered version ####

p <- ggplot(data = Res,aes(Id,Genotype)) +
  geom_raster(aes(fill =log2FoldChange ), stat = "identity") +
  geom_tile(aes(color = Significance),fill = '#00000000', size = 1,width = 0.85,height = 0.85) + 
  scale_color_manual(values = c("#00000000","black"))+
  facet_grid(Fraction~.,space = "free",scales = "free")  +
  scale_fill_paletteer_c("pals::kovesi.diverging_bwr_55_98_c37",
                         limits = c(-5,5),oob = squish,na.value = "white") +
  theme_ohchibi(size_panel_border = 2) +
  theme(
    strip.background.x = element_blank(),
    strip.background.y = element_blank(),
    strip.text.x = element_text(family = "Arial",face = "bold",size = 15),
    strip.text.y = element_text(family = "Arial",face = "bold",size = 20),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(family = "Arial",face = "bold",
                               size = 22,angle = 90,vjust = 0.5,hjust = 1),
    axis.text.y = element_text(family = "Arial",face = "bold",size = 18)
    
  ) +
  scale_x_discrete(expand = c(0,0)) 

oh.save.pdf(p = p,outname = "resub_clustered_heatmap_genotype_syncom_by_fraction.pdf",
            outdir = "../figures/",width = 18,height =8)


# Group version
Res$Id <- Res$Id %>% factor(levels = order_id)

p <- ggplot(data = Res,aes(Id,Genotype)) +
  geom_raster(aes(fill =log2FoldChange ), stat = "identity") +
  geom_tile(aes(color = Significance),fill = '#00000000', size = 1,width = 0.85,height = 0.85) + 
  scale_color_manual(values = c("#00000000","black"))+
  facet_grid(Fraction~Group,space = "free",scales = "free")  +
  scale_fill_paletteer_c("pals::kovesi.diverging_bwr_55_98_c37",
                         limits = c(-5,5),oob = squish,na.value = "white") +
  theme_ohchibi(size_panel_border = 2) +
  theme(
    strip.background.x = element_blank(),
    strip.background.y = element_blank(),
    strip.text.x = element_text(family = "Arial",face = "bold",size = 15),
    strip.text.y = element_text(family = "Arial",face = "bold",size = 20),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(family = "Arial",face = "bold",
                               size = 22,angle = 90,vjust = 0.5,hjust = 1),
    axis.text.y = element_text(family = "Arial",face = "bold",size = 18)
    
  ) +
  scale_x_discrete(expand = c(0,0)) 

oh.save.pdf(p = p,outname = "resub_group_heatmap_genotype_syncom_by_fraction.pdf",
            outdir = "../figures/",width = 18,height =8)

rm(list=ls())
dev.off()
gc()