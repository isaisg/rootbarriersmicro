library(ohchibi)
library(palettesPM)
library(Rmisc)
library(ggpubr)
library(DESeq2)
library(dendextend)
library(multcomp)
library(limma)
library(scales)

set.seed(130816)
setwd('/home/isai/Documents/results/rootbarriersmicro/scripts/')
dir.create("../cleandata")
dir.create("../figures")

### Do a first round to determine the most enriched orders ###
df_deseq_temp <- read.table(file = "../cleandata/deseq2_results_censusmutants_genotypes_inside_fractions_all.tsv",header = T,sep = "\t")
df_deseq_temp <- df_deseq_temp %>% subset(Level == "Order") %>% droplevels
df_deseq_temp$Id <- df_deseq_temp$Id %>% gsub(pattern = "Root; k__Bacteria; ",replacement = "")
df_deseq_temp$padj <- df_deseq_temp$pvalue %>% p.adjust(method = "fdr")
df_sup_sig <- df_deseq_temp %>% subset(padj < 0.05) %$% Id %>% table %>% data.frame()
colnames(df_sup_sig)[1] <- "Id"
with(df_sup_sig,order(-Freq)) %>% df_sup_sig[.,]
by_num <- df_sup_sig %>% subset(Freq >=3) %$% Id %>% as.character

Dat <- readRDS(file = "../cleandata/dat_censusmutants_16s.RDS")

Dat_rab <- Dat$RelativeAbundance

Tab <- Dat_rab$Tab
Tab_batch <- removeBatchEffect(x = Tab,Dat_rab$Map$Rep)

Dat_batch <- create_dataset(Tab = Tab_batch,Map = Dat_rab$Map,Tax = Dat_rab$Tax)



Dat_class <- Dat_batch %>% collapse_by_taxonomy.Dataset(Dat = .,level = 5)
rownames(Dat_class$Tab) <- rownames(Dat_class$Tab) %>% gsub(pattern = "Root; k__Bacteria; ",replacement = "") 

df_freq <- Dat_class$Tab %>% rowMeans %>% data.frame(Id = names(.),RA = .,row.names = NULL)
df_freq <- with(df_freq,order(-RA)) %>% df_freq[.,]

df_freq$RA[1:40] %>% sum
df_freq$Id[1:25]


mchosen <- df_freq$Id[1:40] %>% as.character
mchosen <- c(mchosen ,by_num) %>% as.character %>% unique
Tab_chosen <- which(rownames(Dat_class$Tab) %in% mchosen) %>% Dat_class$Tab[.,] %>% as.matrix

Tab_other <- which(!(rownames(Dat_class$Tab) %in% mchosen)) %>% Dat_class$Tab[.,] %>% as.matrix
Other <- colSums(Tab_other)
Tab_chosen <- rbind(Tab_chosen,Other) 

melted <- Tab_chosen %>% melt
melted$Var1 <- melted$Var1 %>% gsub(pattern = "Root; k__Bacteria; ",replacement = "") 
colnames(melted) <- c("Id","DADA2_Header","RA")

melted <- merge(melted,Dat_class$Map, by = "DADA2_Header")
melted <- which(melted$Id %in% mchosen) %>% melted[.,] %>% droplevels
melted$Id <- melted$Id %>% factor(levels = mchosen)



## Read the  deseq 2 results ############
df_deseq <- read.table(file = "../cleandata/deseq2_results_censusmutants_genotypes_inside_fractions_all.tsv",header = T,sep = "\t")
df_deseq <- df_deseq %>% subset(Level == "Order") %>% droplevels
df_deseq$Id <- df_deseq$Id %>% gsub(pattern = "Root; k__Bacteria; ",replacement = "")
df_deseq <- df_deseq[which(df_deseq$Id %in% mchosen),] %>% droplevels
df_deseq$padjClass <- df_deseq$pvalue %>% p.adjust(method = "fdr")
df_deseq$SignificanceClass <- NA
df_deseq$SignificanceClass[which(df_deseq$padjClass < 0.1)] <- "S"



df_sum <- aggregate(RA ~Id+Genotype+Fraction,melted,mean)

df_sum <- merge(df_sum,df_deseq[,c(5,7,8,10,11,12,13,14)], by = c("Id","Fraction","Genotype"),all.x = TRUE) 

df_sum$RAAbs <- abs(df_sum$RA)


df_sum %>% subset(SignificanceClass == "S")
df_sum$RAAbs %>% sort %>% plot



df_sum$Id <- df_sum$Id %>% factor(levels = mchosen %>% rev)


#Here we should use RA instead
#melted_plot <- acast(data =df_sum,formula = Id~Genotype+Fraction,value.var = "RAAbs") %>% t %>%
#  scale %>% t  %>% melt

melted_plot <- acast(data =df_sum,formula = Id~Genotype+Fraction,value.var = "RA") %>% t %>%
  scale %>% t  %>% melt


melted_plot$Genotype <- melted_plot$Var2 %>% gsub(pattern = "_Root|_Shoot|_Soil",replacement = "")
melted_plot$Fraction <- melted_plot$Var2 %>% gsub(pattern = ".*_",replacement = "")
colnames(melted_plot)[1] <- "Id"

melted_plot <- merge(melted_plot,df_deseq[,c(5,7,8,10,11,12,13,14)], by = c("Id","Fraction","Genotype"),all.x = TRUE) 
melted_plot$Id <- melted_plot$Id %>% factor(levels = mchosen %>% rev)

### Adjust the panels
df_groups <- read.table(file = "../rawdata/df_genotypes_to_classification.csv",header = T,sep = ",")
df_groups$Num <- 1:20

melted_plot <- merge(melted_plot,df_groups, by = "Genotype")
melted_plot$Num <- melted_plot$Num %>% factor
melted_plot$Phylum <- melted_plot$Id %>% gsub(pattern = "\\;.*",replacement = "") %>% 
  gsub(pattern = "p__",replacement = "") %>% factor

melted_plot$Order <- melted_plot$Id %>% gsub(pattern = ".*o__",replacement = "") %>% factor
order_orders <- rev(melted_plot$Order %>% levels)
melted_plot$Order <- melted_plot$Order %>% factor(levels = order_orders)
p <- ggplot(data = melted_plot,aes(Num,Order )) +
  geom_raster(aes(fill = value)) +  
  geom_tile(aes(color = SignificanceClass),fill = '#00000000', size = 1.5,width = 0.85,height = 0.85) + 
  scale_fill_paletteer_c(package = "pals",palette = "kovesi.diverging_bwr_55_98_c37",
                         limits = c(-2,2),oob = squish) + 
  facet_grid(Phylum~Fraction,space = "free",scales = "free") +
  scale_color_manual(values = c("black"),na.value = "#00000000")+
  theme_ohchibi(size_axis_text.x = 20,
                angle_text.x = 90,
                size_axis_text.y = 20,
                size_axis_title.x = 22,
                size_axis_title.y = 0,
                legend_proportion_size = 1,
                size_title_text = 12,
                size_legend_text = 12,
                size_panel_border = 1.5,
                size_lines_panel = 0) +
  theme(axis.ticks = element_blank(),
        panel.spacing = unit(0.075, "lines"),
        strip.background = element_blank(),
        strip.text.x = element_text(family = "Arial",face = "bold",size = 20),
        strip.text.y = element_text(family = "Arial",face = "bold",size = 10),
        axis.text.y = element_text(family = "Arial",face = "bold",hjust = 1,vjust = 0.5,size = 10),
        axis.text.x = element_text(hjust = 1,vjust = 0.5,family = "Arial",face = "bold",size = 15),
        axis.title.x = element_blank()
  )     +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) 
melted_plot$Num <- melted_plot$Num %>% factor(levels = melted_plot$Num %>% levels %>% rev)
### Inverted plot ####
p <- ggplot(data = melted_plot,aes(Order,Num )) +
  geom_raster(aes(fill = value)) +  
  geom_tile(aes(color = SignificanceClass),fill = '#00000000', size = 1.5,width = 0.85,height = 0.85) + 
  scale_fill_paletteer_c(package = "pals",palette = "kovesi.diverging_bwr_55_98_c37",
                         limits = c(-2,2),oob = squish) + 
  facet_grid(Fraction~Phylum,space = "free",scales = "free") +
  scale_color_manual(values = c("black"),na.value = "#00000000")+
  theme_ohchibi(size_axis_text.x = 20,
                angle_text.x = 90,
                size_axis_text.y = 20,
                size_axis_title.x = 22,
                size_axis_title.y = 0,
                legend_proportion_size = 1,
                size_title_text = 12,
                size_legend_text = 12,
                size_panel_border = 1.5,
                size_lines_panel = 0) +
  theme(axis.ticks = element_blank(),
        panel.spacing = unit(0.075, "lines"),
        strip.background = element_blank(),
        strip.text.x = element_text(family = "Arial",face = "bold",size = 10),
        strip.text.y = element_text(family = "Arial",face = "bold",size = 20),
        axis.text.y = element_text(family = "Arial",face = "bold",hjust = 1,vjust = 0.5,size = 10),
        axis.text.x = element_text(hjust = 1,vjust = 0.5,family = "Arial",face = "bold",size = 15),
        axis.title.x = element_blank()
  )     +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) 
oh.save.pdf(p = p,outname = "heatmap_census_genotype_effect_per_fraction_vertical.pdf",outdir = "../figures/",width = 18,height = 14)

