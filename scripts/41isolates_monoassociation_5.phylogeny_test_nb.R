library(ohchibi)
library(ggtree)
library(palettesPM)
library(paletteer)
library(car)



#Define palettes 
set.seed(130816)
paleta_enrichment <- c("#B3B3B3","#2480FF","#E65037")
names(paleta_enrichment) <- c("NoSignificant","Down","Up")
shade_nb <- "#408078"

size_point <- 10
size_axis_title_y <- 20
method_adjust <- "bonferroni"
#Horizontal version of the distribution of the traits
#Merge datasets into a single structure
set.seed(130816)
setwd('/home/isai/Documents/results/rootbarriersmicro/scripts/')

Dat_mono <- readRDS(file = "../cleandata/dat_41isolatesmonoassociation_cfu_suberin_primroot_dryweight_ionome.RDS")
Dat_screen <- readRDS(file = "../cleandata/dat_screeningbarriers.RDS")

#Load the phylogenetic tree
tree <- Dat_screen$tree
meta <- Dat_mono$Map
chosen_strains <- meta$Strains_original %>% unique %>% 
  grep(pattern = "CDEF|NB",invert = T,value = T)
chosen_taxonoid <- meta$taxon_oid %>% unique %>% na.omit %>% as.character
todrop <- which(!(tree$tip.label %in% chosen_taxonoid)) %>% tree$tip.label[.]
tree <- ape::drop.tip(phy = tree,tip = todrop)

df_un <- meta[,-2] %>% unique
df_un <- match(chosen_strains,df_un$Strains_original) %>%
  df_un[.,] %>% droplevels
rownames(df_un) <- df_un$taxon_oid %>% as.character
df_un <- df_un[,c(2,1,3:ncol(df_un))]

#Create palette
df_colors <- pm.colors.phyla() %>% 
  data.frame(Phyla = names(.), Color = .,row.names = NULL)
df_deino <- data.frame(Phyla = "Deinococcus-Thermus", Color = "#7FE85D")
df_colors <- rbind(df_colors,df_deino)

paleta_phyla <- df_colors$Color %>% as.character
names(paleta_phyla) <- df_colors$Phyla

##### Phylogenetic tree ####
p <- ggtree(tree,ladderize = F,size = 2) + 
  #The expand parth is fundamental to make the tree fit the composition with the other plots
  scale_y_reverse(expand = c(0.02,0.02))
p + geom_text(aes(label =node))
dev.off()
p  <- ggtree::flip(p,46,57) 

p_tree <- p  %<+% df_un + aes(color=(Phylum)) +
  geom_tiplab(align = T,size = 0,linetype = "solid") +
  scale_color_manual(values = paleta_phyla,na.value = "black") +
  theme(legend.position = "none") 
p_tree <- p_tree + coord_flip()  + scale_x_reverse()

#Determine the order of the tips
d <- p$data %>% subset(isTip)
mtips <- with(d, label[order(y, decreasing=T)])
mstrains <- match(mtips,meta$taxon_oid) %>% meta$Strains_original[.] %>% as.character

##### Suberin ####
df_sub <- Dat_mono$Suberin
df_nb <- df_sub %>% subset(Strains == "NB") %>% droplevels
df_sub <- df_sub %>% subset(Strains != "NB") %>% droplevels
mean_nb <- df_nb$Norm_Distance_tip_cz %>% mean

#Test differences 
Res_sub <- NULL
for(st in df_sub$Strains %>% unique){
  df_temp <- df_sub %>% subset(Strains == st) %>% droplevels
  df_temp_nb <- which(df_nb$Batch %in% (df_temp$Batch %>% unique)) %>%
    df_nb[.,]
  
  x <- df_temp_nb$Norm_Distance_tip_cz
  y <- df_temp$Norm_Distance_tip_cz 
  df_temp <- rbind(df_temp_nb,df_temp) 
  df_temp$Strains <- df_temp$Strains   %>% factor
  m <- car::leveneTest(Norm_Distance_tip_cz ~ Strains,data = df_temp) 
  pval <- m$`Pr(>F)`[1]
  if(pval < 0.05){
    m1 <- t.test(x = x,y = y,exact = F,var.equal = F)
    
  }else{
    m1 <- t.test(x = x,y = y,exact = F,var.equal = T)
    
  }
  Res_sub <- data.frame(Strains = st,Mean = mean(y),p.value = m1$p.value,
                        est = mean(y)-mean_nb) %>%
    rbind(Res_sub,.)
}
Res_sub$p.adj <- Res_sub$p.value %>% p.adjust(method = method_adjust)
Res_sub$Significance <- rep("NoSignificant")
Res_sub$Significance[which(Res_sub$Mean > mean_nb & Res_sub$p.adj < 0.1)] <- "Up"
Res_sub$Significance[which(Res_sub$Mean < mean_nb & Res_sub$p.adj < 0.1)] <- "Down"

Res_suberin_raw <- Res_sub
#Order according to the phylogeny
df_sub <- which(df_sub$Strains %in% mstrains) %>% df_sub[.,] %>%
  droplevels
df_sub$Strains <- factor(x = df_sub$Strains,levels = mstrains)
Res_sub <- match(mstrains,Res_sub$Strains) %>%
  Res_sub[.,]
Res_sub$Strains <- Res_sub$Strains %>% factor(levels = mstrains)

Res_suberin <- Res_sub


df_sub$Significance <- rep("NoSignificant")
ch <- Res_sub %>% subset(Significance == "Up") %$% Strains %>%
  as.character
df_sub$Significance[which(df_sub$Strains %in% ch)] <- "Up"

ch <- Res_sub %>% subset(Significance == "Down") %$% Strains %>%
  as.character
df_sub$Significance[which(df_sub$Strains %in% ch)] <- "Down"

quant_nb <- df_nb$Norm_Distance_tip_cz %>% quantile 
df_seg <- data.frame(one = quant_nb[1],two = quant_nb[2],three = quant_nb[3],
                     four = quant_nb[4],five = quant_nb[5])



p_sub <- ggplot(data = df_sub,aes(Strains,Norm_Distance_tip_cz)) +
  geom_point(size = NA,stroke = NA) +
  geom_rect(data = df_seg,
            mapping = aes(ymin = two,ymax = four,xmin = 0,xmax = Inf),
            inherit.aes = F,color = NA,fill = shade_nb,alpha = 0.5) +
  geom_boxplot(aes (fill = Significance),size = 1, color = "black",outlier.colour = NA,outlier.fill = NA,outlier.size = NA) + 
  geom_point(shape = 20,color = "black",size = 3,alpha = 0.5) + 
  theme_ohchibi(size_panel_border = 2) +
  #scale_x_discrete(expand = c(0,0)) +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = size_axis_title_y)
    
  ) + 
  ylab(label = "Suberine (Distance to tip)")+ 
  scale_fill_manual(values = paleta_enrichment)



### Primary root elongation ###
df_sub <- Dat_mono$PrimaryRoot
df_nb <- df_sub %>% subset(Strains == "NB") %>% droplevels
df_sub <- df_sub %>% subset(Strains != "NB") %>% droplevels
mean_nb <- df_nb$Norm_Main_root_elongation %>% mean

#Test differences 
Res_sub <- NULL
for(st in df_sub$Strains %>% unique){
  df_temp <- df_sub %>% subset(Strains == st) %>% droplevels
  df_temp_nb <- which(df_nb$Batch %in% (df_temp$Batch %>% unique)) %>%
    df_nb[.,]
  x <- df_temp_nb$Norm_Main_root_elongation
  y <- df_temp$Norm_Main_root_elongation 
  df_temp <- rbind(df_temp_nb,df_temp) 
  df_temp$Strains <- df_temp$Strains   %>% factor
  m <- car::leveneTest(Norm_Main_root_elongation ~ Strains,data = df_temp) 
  pval <- m$`Pr(>F)`[1]
  if(pval < 0.05){
    m1 <- t.test(x = x,y = y,exact = F,var.equal = F)
    
  }else{
    m1 <- t.test(x = x,y = y,exact = F,var.equal = T)
    
  }
  Res_sub <- data.frame(Strains = st,Mean = mean(y),p.value = m1$p.value,
                        est = mean(y) -mean_nb) %>%
    rbind(Res_sub,.)
}
Res_sub$p.adj <- Res_sub$p.value %>% p.adjust(method = method_adjust)
Res_sub$Significance <- rep("NoSignificant")
Res_sub$Significance[which(Res_sub$Mean > mean_nb & Res_sub$p.adj < 0.1)] <- "Up"
Res_sub$Significance[which(Res_sub$Mean < mean_nb & Res_sub$p.adj < 0.1)] <- "Down"

Res_root_raw <- Res_sub

#Order according to the phylogeny
df_sub <- which(df_sub$Strains %in% mstrains) %>% df_sub[.,] %>%
  droplevels
df_sub$Strains <- factor(x = df_sub$Strains,levels = mstrains)
Res_sub <- match(mstrains,Res_sub$Strains) %>%
  Res_sub[.,]
Res_sub$Strains <- Res_sub$Strains %>% factor(levels = mstrains)

Res_primroot <- Res_sub 

df_sub$Significance <- rep("NoSignificant")
ch <- Res_sub %>% subset(Significance == "Up") %$% Strains %>%
  as.character
df_sub$Significance[which(df_sub$Strains %in% ch)] <- "Up"

ch <- Res_sub %>% subset(Significance == "Down") %$% Strains %>%
  as.character
df_sub$Significance[which(df_sub$Strains %in% ch)] <- "Down"

quant_nb <- df_nb$Norm_Main_root_elongation %>% quantile 
df_seg <- data.frame(one = quant_nb[1],two = quant_nb[2],three = quant_nb[3],
                     four = quant_nb[4],five = quant_nb[5])


p_root <- ggplot(data = df_sub,aes(Strains,Norm_Main_root_elongation)) +
  geom_point(size = NA,stroke = NA) +
  geom_rect(data = df_seg,
            mapping = aes(ymin = two,ymax = four,xmin = 0,xmax = Inf),
            inherit.aes = F,color = NA,fill = shade_nb,alpha = 0.5) +
  geom_boxplot(aes (fill = Significance),size = 1, color = "black",outlier.colour = NA,outlier.fill = NA,outlier.size = NA) + 
  #geom_point(shape = 20,color = "black",size = 3,alpha = 0.5) + 
  geom_sina(size = 2, alpha = 0.1,color = "black") +
  theme_ohchibi(size_panel_border = 2) +
  #scale_x_discrete(expand = c(0,0)) +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = size_axis_title_y)
    
  ) + 
  ylab(label = "Primary root elongation (cm)")+ 
  scale_fill_manual(values = paleta_enrichment)


#### Dry weight #########
df_sub <- Dat_mono$DryWeight
df_nb <- df_sub %>% subset(Strains == "NB") %>% droplevels
df_sub <- df_sub %>% subset(Strains != "NB") %>% droplevels
mean_nb <- df_nb$Norm_WeightbyPlant %>% mean
#Remove clear outlier
df_sub <- df_sub[-which(df_sub$Norm_WeightbyPlant > 3),] %>% droplevels

#Test differences 
Res_sub <- NULL
for(st in df_sub$Strains %>% unique){
  df_temp <- df_sub %>% subset(Strains == st) %>% droplevels
  df_temp_nb <- which(df_nb$Batch %in% (df_temp$Batch %>% unique)) %>%
    df_nb[.,]
  x <- df_temp_nb$Norm_WeightbyPlant
  y <- df_temp$Norm_WeightbyPlant 
  df_temp <- rbind(df_temp_nb,df_temp) 
  df_temp$Strains <- df_temp$Strains   %>% factor
  m <- car::leveneTest(Norm_WeightbyPlant ~ Strains,data = df_temp) 
  pval <- m$`Pr(>F)`[1]
  if(pval < 0.05){
    m1 <- t.test(x = x,y = y,exact = F,var.equal = F)
    
  }else{
    m1 <- t.test(x = x,y = y,exact = F,var.equal = T)
    
  }
  Res_sub <- data.frame(Strains = st,Mean = mean(y),p.value = m1$p.value,
                        est = mean(y)-mean_nb) %>%
    rbind(Res_sub,.)
}
Res_sub$p.adj <- Res_sub$p.value %>% p.adjust(method = method_adjust)
Res_sub$Significance <- rep("NoSignificant")
Res_sub$Significance[which(Res_sub$Mean > mean_nb & Res_sub$p.adj < 0.1)] <- "Up"
Res_sub$Significance[which(Res_sub$Mean < mean_nb & Res_sub$p.adj < 0.1)] <- "Down"

Res_weight_raw <- Res_sub

#Order according to the phylogeny
df_sub <- which(df_sub$Strains %in% mstrains) %>% df_sub[.,] %>%
  droplevels
df_sub$Strains <- factor(x = df_sub$Strains,levels = mstrains)
Res_sub <- match(mstrains,Res_sub$Strains) %>%
  Res_sub[.,]
Res_sub$Strains <- Res_sub$Strains %>% factor(levels = mstrains)

Res_weight <- Res_sub

df_sub$Significance <- rep("NoSignificant")
ch <- Res_sub %>% subset(Significance == "Up") %$% Strains %>%
  as.character
df_sub$Significance[which(df_sub$Strains %in% ch)] <- "Up"

ch <- Res_sub %>% subset(Significance == "Down") %$% Strains %>%
  as.character
df_sub$Significance[which(df_sub$Strains %in% ch)] <- "Down"

quant_nb <- df_nb$Norm_WeightbyPlant %>% quantile 
df_seg <- data.frame(one = quant_nb[1],two = quant_nb[2],three = quant_nb[3],
                     four = quant_nb[4],five = quant_nb[5])

p_dry <- ggplot(data = df_sub,aes(Strains,Norm_WeightbyPlant)) +
  geom_point(size = NA,stroke = NA) +
  geom_rect(data = df_seg,
            mapping = aes(ymin = two,ymax = four,xmin = 0,xmax = Inf),
            inherit.aes = F,color = NA,fill = shade_nb,alpha = 0.5) +
  geom_boxplot(aes (fill = Significance),size = 1, color = "black",outlier.colour = NA,outlier.fill = NA,outlier.size = NA) + 
  geom_point(shape = 20,color = "black",size = 3,alpha = 0.5) + 
  theme_ohchibi(size_panel_border = 2) +
  #scale_x_discrete(expand = c(0,0)) +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = size_axis_title_y)
    
  ) + 
  ylab(label = "Dry Weight")+ 
  scale_fill_manual(values = paleta_enrichment)


#### CFU data ####
df_sub <- Dat_mono$CFU
df_nb <- df_sub %>% subset(Strains == "NB") %>% droplevels
df_sub <- df_sub %>% subset(Strains != "NB") %>% droplevels

#Load screening barriers dat
Dat <- readRDS(file = "../cleandata/dat_screeningbarriers.RDS")
meta <- Dat$Map[,c(16,7)]
colnames(meta)[1] <- "Strains"
df_sub <- merge(df_sub,meta, by = "Strains", all.x = TRUE)
df_sub$Nom <- paste0(df_sub$Genus," ",df_sub$Strains)

df_sub$Strains <- df_sub$Strains %>% factor(levels = mstrains)

order_nom <- with(df_sub,order(Strains)) %>%
  df_sub[.,] %$% Nom %>% as.character %>%
  unique
df_sub$Nom <- df_sub$Nom %>% factor(levels = order_nom)

paleta_fractions <- c(pm.colors.fractions(),"#FFC6CD")
names(paleta_fractions)[11] <- "Agar"

p_cfu_agar <- df_sub %>% subset(Fraction == "Agar") %>% droplevels %>%
  ggplot(data = .,aes(Strains,Log_Norm_cfu_weight)) +
  geom_point(size = NA,stroke = NA) +
  geom_rect(data = df_seg,
            mapping = aes(ymin = 0,ymax = 0.05,xmin = 0,xmax = Inf),
            inherit.aes = F,color = NA,fill = shade_nb,alpha = 0.5) +
  geom_boxplot(aes (fill = Fraction),size = 1, color = "black",outlier.colour = NA,outlier.fill = NA,outlier.size = NA) + 
  geom_point(shape = 20,color = "black",size = 3,alpha = 0.5) + 
  theme_ohchibi(size_panel_border = 2) +
  #scale_x_discrete(expand = c(0,0)) +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = size_axis_title_y)
    
  ) + 
  ylab(label = "CFU Agar")+ 
  scale_fill_manual(values = paleta_fractions)


p_cfu_root <- df_sub %>% subset(Fraction == "Root") %>% droplevels %>%
  ggplot(data = .,aes(Strains,Log_Norm_cfu_weight)) +
  geom_point(size = NA,stroke = NA) +
  geom_rect(data = df_seg,
            mapping = aes(ymin = 0,ymax = 0.05,xmin = 0,xmax = Inf),
            inherit.aes = F,color = NA,fill = shade_nb,alpha = 0.5) +
  geom_boxplot(aes (fill = Fraction),size = 1, color = "black",outlier.colour = NA,outlier.fill = NA,outlier.size = NA) + 
  geom_point(shape = 20,color = "black",size = 3,alpha = 0.5) + 
  theme_ohchibi(size_panel_border = 2) +
  #scale_x_discrete(expand = c(0,0)) +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = size_axis_title_y)
    
  ) + 
  ylab(label = "CFU Root")+ 
  scale_fill_manual(values = paleta_fractions)




p_cfu_shoot <- df_sub %>% subset(Fraction == "Shoot") %>% droplevels %>%
  ggplot(data = .,aes(Nom,Log_Norm_cfu_weight)) +
  geom_point(size = NA,stroke = NA) +
  geom_rect(data = df_seg,
            mapping = aes(ymin = 0,ymax = 0.05,xmin = 0,xmax = Inf),
            inherit.aes = F,color = NA,fill = shade_nb,alpha = 0.5) +
  geom_boxplot(aes (fill = Fraction),size = 1, color = "black",outlier.colour = NA,outlier.fill = NA,outlier.size = NA) + 
  geom_point(shape = 20,color = "black",size = 3,alpha = 0.5) + 
  theme_ohchibi(size_panel_border = 2) +
  #scale_x_discrete(expand = c(0,0)) +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_text(angle = 90,family = "Arial",
                               face = "bold",size = 30,hjust = 1,vjust = 0.5),
    #axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = size_axis_title_y)
    
  ) + 
  ylab(label = "CFU Shoot")+ 
  scale_fill_manual(values = paleta_fractions)


composition <- egg::ggarrange(p_tree,p_sub ,p_root, p_dry,
                              p_cfu_agar,p_cfu_root,p_cfu_shoot,ncol = 1,heights = c(0.4,1,1,1,1,1,1))


oh.save.pdf(p = composition,outname = "41isolates_monoassociation_phenotypes_phylogeny.pdf",outdir = "../figures",width = 25,height = 25)



### Test the phylogenetic signal of the predictors
df <- merge(Res_suberin,Res_primroot, by = "Strains") %>%
  merge(Res_weight, by = "Strains")
colnames(df) <- c("Strains",
                  "Suberin_Mean","Suberin_p.value","Suberin_p.adj","Suberin_est","Suberin_Significance",
                  "PrimRoot_Mean","PrimRoot_p.value","PrimRoot_p.adj","PrimmtRoot_est","PrimRoot_Signifiance",
                  "DryWeight_Mean","DryWeight_p.value","DryWeight_p.adj","DryWeight_est","DryWeight_Significance")

#Merge here with the taxon_oids
Dat <- readRDS(file = "../cleandata/dat_screeningbarriers.RDS")
meta <- Dat$Map[,c(16,7,1)]
colnames(meta)[1] <- "Strains"

df <- merge(df,meta, by = "Strains", all.x = TRUE)

Tab <- data.frame(Suberin = df$Suberin_Mean,PrimRoot = df$PrimRoot_Mean,
                  DryWeight = df$DryWeight_Mean) %>% as.matrix
rownames(Tab) <- df$taxon_oid

Tab <- tree$tip.label %>% match(rownames(Tab)) %>%
  Tab[.,]

#Test phylogenetic signal
df_ps <- apply(X = Tab,MARGIN = 2,
               FUN = function(x)phylosig(tree = tree,x,method = "lambda",test = TRUE) %>% unlist) %>%
  melt


## Combine all raw structures and put them in the same data.frame ###
Res_root_raw$Variable <- "Root"
Res_suberin_raw$Variable <- "Suberin"
Res_weight_raw$Variable <- "Weight"

allRes <- rbind(Res_root_raw,Res_suberin_raw,Res_weight_raw)
write.table(x = allRes,file = "../cleandata/res_41_suberinraw_primroot_weight_vsnb.tsv",
            append = F,quote = F,sep = "\t",row.names = F,col.names = T)
rm(list=ls())
dev.off()

