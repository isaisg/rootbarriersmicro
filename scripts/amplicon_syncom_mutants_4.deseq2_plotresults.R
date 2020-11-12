library(ohchibi)
library(DESeq2)
library(scales)
library(paletteer)

set.seed(130816)
setwd('/home/isai/Documents/results/rootbarriersmicro/scripts/')

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
Res_root <- Res

###### Agar #########
dds <- Dat$AgarMutants

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
Res_agar <- Res

###### Shoot #########
dds <- Dat$ShootMutants

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
Res_shoot <- Res

#Create structure to plot
Res_root$Fraction <- "Root"
Res_agar$Fraction <- "Agar"
Res_shoot$Fraction <- "Shoot"

Res <- rbind(Res_agar,Res_root,Res_shoot)

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

#Now do the plotting
p <- ggplot(data = Res,aes(Genotype,Id)) +
  geom_raster(aes(fill =log2FoldChange ), stat = "identity") +
  geom_tile(aes(color = Significance),fill = '#00000000', size = 1,width = 0.85,height = 0.85) + 
  scale_color_manual(values = c("#00000000","black"))+
  facet_grid(Group~Fraction,space = "free",scales = "free")  +
  scale_fill_paletteer_c(package = "pals",palette = "kovesi.diverging_bwr_55_98_c37",
                         limits = c(-5,5),oob = squish,na.value = "white") +
  theme_ohchibi(size_panel_border = 2) +
  theme(
    strip.background.x = element_blank(),
    strip.background.y = element_blank(),
    strip.text.x = element_text(family = "Arial",face = "bold",size = 25),
    strip.text.y = element_text(family = "Arial",face = "bold",size = 20),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(family = "Arial",face = "bold",
                               size = 22,angle = 90,vjust = 0.5,hjust = 1),
    axis.text.y = element_text(family = "Arial",face = "bold",size = 18)
    
  ) +
  scale_x_discrete(expand = c(0,0)) 


#Inverted version
p <- ggplot(data = Res,aes(Id,Genotype)) +
  geom_raster(aes(fill =log2FoldChange ), stat = "identity") +
  geom_tile(aes(color = Significance),fill = '#00000000', size = 1,width = 0.85,height = 0.85) + 
  scale_color_manual(values = c("#00000000","black"))+
  facet_grid(Fraction~Group,space = "free",scales = "free")  +
  scale_fill_paletteer_c(package = "pals",palette = "kovesi.diverging_bwr_55_98_c37",
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

oh.save.pdf(p = p,outname = "heatmap_genotype_syncom_by_fraction.pdf",
            outdir = "../figures/",width = 18,height =8)




### Check patterns across groups
p <- ggplot(data = Res,aes(Group,log2FoldChange)) +
  geom_hline(yintercept = 0,size = 1) +
  geom_point(aes(color = Group),alpha = 0.4,size = 1) +
  stat_summary(fun.y = mean, geom = "point",size =4,aes(color = Group)) +
  facet_grid(.~Fraction + Genotype) +
  theme_ohchibi(size_panel_border = 2) +
  scale_y_continuous(limits = c(-10,10),oob = rescale_none) +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background.x = element_blank(),
    strip.text.x = element_text(size = 10,family = "Arial",face = "bold") ,
    axis.title.x = element_blank()
  ) +
  scale_color_paletteer_d(package ="LaCroixColoR",palette = "PassionFruit" ) 

oh.save.pdf(p = p,outname = "bubuplot_log2foldchange_genotype_syncom_by_fraction.pdf",
            outdir = "../figures/",width = 21,height =7)



#### quantify number of enrichments
#and use an statistic to say we are seeing something random
Res$CountUp <- 0
Res$CountUp[which(Res$Significance == "Significant" & Res$log2FoldChange > 0)] <- 1

Res$CountDown <- 0
Res$CountDown[which(Res$Significance == "Significant" & Res$log2FoldChange <0)] <- 1



acast(data = Res,formula =Group~Fraction + Genotype,fun.aggregate = sum,value.var = "CountUp" )  %>%
  rowSums()
acast(data = Res,formula =Group~Fraction + Genotype,fun.aggregate = sum,value.var = "CountDown" )  %>%
  rowSums()

melt_up <- dcast(data = Res,formula =Group~Fraction + Genotype,fun.aggregate = sum,value.var = "CountUp" )  %>%
  melt

melt_down <- dcast(data = Res,formula =Group~Fraction + Genotype,fun.aggregate = sum,value.var = "CountDown" )  %>%
  melt


melt_up$Direction <- "Enriched"
melt_down$Direction <- "Depleted"

melted <- rbind(melt_up,melt_down)

melted$Fraction <- melted$variable %>% gsub(pattern = "_.*",replacement = "")
melted$Genotype <- melted$variable %>% gsub(pattern = ".*_",replacement = "")

melted$Direction <- melted$Direction %>% factor(levels = c("Enriched","Depleted"))

p <- ggplot(data = melted,aes(Group,value)) +
  geom_bar(stat = "identity",color = "#414141",fill = "#414141") + 
  facet_grid(Fraction~Direction) + 
  theme_ohchibi(size_panel_border = 2) +
  ylab(label = "Count") + 
  theme(
    strip.background.x = element_blank(),
    strip.background.y = element_blank(),
    strip.text.x = element_text(size = 20,family = "Arial",face = "bold"),
    strip.text.y = element_text(size = 20,family = "Arial",face = "bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1)
  )

oh.save.pdf(p = p,outname = "barplot_amplicon_mutants_groups.pdf",
            outdir = "../figures/",width = 10,height =7)

rm(list=ls())