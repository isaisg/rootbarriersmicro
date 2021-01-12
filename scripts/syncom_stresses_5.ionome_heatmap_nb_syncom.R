library(ohchibi)
library(emmeans)
library(paletteer)
library(scales)
library(ggtree)

set.seed(130816)
setwd('/home/isai/Documents/results/rootbarriersmicro/scripts/')

Dat <- readRDS(file = "../cleandata/dat_syncom_stresses.RDS") %$%
  Ionome

paleta_enrichment <- readRDS(file = "../cleandata/dat_syncom_stresses.RDS") %$%
  paleta_enrichment

#Perform testing inside each ion
#We wanna do pairwise comparison between  StressNB vs Full NB
#StressSC vs Full SC and #Intra Stress SC vs Stress NB
Tab_z <- Dat$Tab %>% t %>% scale %>% t 
melted <- Tab_z %>% melt
colnames(melted) <- c("Ion","Index","zscore")

melted <- merge(melted,Dat$Map, by = "Index") 
melted$Stress <- melted$Stress %>% factor %>% relevel(ref = "Full")
melted$Type <- melted$Type %>% factor(levels = c("NB","HK","SynCom"))

#Do a model per ion
Res_em <- NULL
Res_p <- NULL
for(ion in melted$Ion %>% unique){
  melted_sub <- melted %>% subset(Ion  == ion) %>% droplevels
  m1 <- lm(data = melted_sub,formula = zscore ~ Type + Stress + Type:Stress)
  m1_res <- emmeans(m1,pairwise~Type:Stress,adjust = "none")
  m1_em <- m1_res$emmeans %>% as.data.frame
  m1_p <- m1_res$contrasts %>% as.data.frame
  m1_em$Ion <- ion
  m1_p$Ion <- ion
  Res_em <- rbind(Res_em,m1_em)
  Res_p <- rbind(Res_p,m1_p)
}

#Arrange the pvalues dataframe
Res_p$contrast <- Res_p$contrast %>% gsub(pattern = ",",replacement = "_") %>%
  gsub(pattern = " ",replacement = "")

#Determine the pairwise contrasts
c1 <- Res_p$contrast %>% grep(pattern = "Full.*Full",value = T)  %>% unique
c2 <- Res_p$contrast %>% grep(pattern = "Mn.*Mn",value = T)  %>% unique
c3 <- Res_p$contrast %>% grep(pattern = "Zn.*Zn",value = T)  %>% unique
c4 <- Res_p$contrast %>% grep(pattern = "Fe.*Fe",value = T)  %>% unique
c5 <- Res_p$contrast %>% grep(pattern = "Mo.*Mo",value = T)  %>% unique
c6 <- Res_p$contrast %>% grep(pattern = "_-K.*_-K",value = T)  %>% unique
c7 <- Res_p$contrast %>% grep(pattern = "P.*P",value = T)  %>% unique
c8 <- Res_p$contrast %>% grep(pattern = "NaCl.*NaCl",value = T)  %>% unique

#Now the contrast against Full in each type of syncom
c9 <- Res_p$contrast %>% grep(pattern = "NB.*NB",value = T) %>%
  grep(pattern = "Full",value = T) %>% unique

c10 <- Res_p$contrast %>% grep(pattern = "HK.*HK",value = T) %>%
  grep(pattern = "Full",value = T) %>% unique


c11 <- Res_p$contrast %>% grep(pattern = "SynCom.*SynCom",value = T) %>%
  grep(pattern = "Full",value = T) %>% unique

allc <- c(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11)
allc <- allc %>% unique

Res_p_sub <- which(Res_p$contrast %in% allc) %>%
  Res_p[.,] %>% droplevels

Res_p_sub$p.adj <- Res_p_sub$p.value %>% p.adjust(method = "fdr")
Res_p_sub$Significance <- "NoSignificant"
pval_thres <- 0.05
Res_p_sub$Significance[which(Res_p_sub$p.adj < pval_thres)] <- "Significant"


#Against the full 
agfull <- c(c9,c10,c11)
Res_agfull <- which(Res_p_sub$contrast %in% agfull) %>%
  Res_p_sub[.,] %>% droplevels

df_temp <- Res_agfull$contrast %>% gsub(pattern = "NB_Full|HK_Full|SynCom_Full",replacement = "") %>%
  gsub(pattern = "^-",replacement = "") %>% as.data.frame
colnames(df_temp) <- c("UId")
Res_agfull <- cbind(Res_agfull,df_temp)
Res_agfull$UId <- paste(Res_agfull$UId,Res_agfull$Ion,sep = "_")
df_temp <- data.frame(UId = Res_agfull$UId,
                      p.adj_agfull =Res_agfull$p.adj, estimate_agfull = Res_agfull$estimate)
df_temp$Significance_agfull <- "NoSignificant"
df_temp$Significance_agfull[which(df_temp$p.adj_agfull < pval_thres & df_temp$estimate_agfull > 0)] <- "Down"
df_temp$Significance_agfull[which(df_temp$p.adj_agfull < pval_thres & df_temp$estimate_agfull < 0)] <- "Up"


#Structure to plot
Res_em$UId <- paste(Res_em$Type,Res_em$Stress,Res_em$Ion, sep = "_")
Res_em <- merge(Res_em,df_temp, by = "UId",all.x = TRUE)


#Analysis for the Syncoms vs NB comparison
Res_em_sub <- Res_em %>% subset(Type != "HK") %>% droplevels

#Order the stresses based on level of extremeness
Res_em_sub$Stress <- Res_em_sub$Stress %>% factor(levels = c("Full","-P","-Mn","-Mo","-Zn","-K","+NaCl","-Fe"))


#Append the intra stress comparison
intrast <- c(c1,c2,c3,c4,c5,c6,c7,c8) %>% unique
Res_intrast <- which(Res_p_sub$contrast %in% intrast) %>%
  Res_p_sub[.,] %>% droplevels
df_temp <- Res_intrast$contrast %>% grep(pattern = "NB.*SynCom",value = F) %>%
  Res_intrast[.,] %>% droplevels
df_temp$Stress <- df_temp$contrast %>% gsub(pattern = ".*_",replacement = "")
df_temp$SI <- paste(df_temp$Stress,df_temp$Ion,sep = "_")



df_temp <- data.frame(SI = df_temp$SI,p.value_intra = df_temp$p.value,
                      p.adj_intra =df_temp$p.adj, estimate_intra = df_temp$estimate)
df_temp$Significance_intra <- "NoSignificant"
df_temp$Significance_intra[which(df_temp$p.adj_intra < pval_thres)] <- "Significant"


Res_em_sub$SI <- paste(Res_em_sub$Stress,Res_em_sub$Ion,sep = "_")

Res_em_sub <- merge(Res_em_sub,df_temp, by = "SI", all.x = TRUE)


Res_em_sub$Significance_agfull <- Res_em_sub$Significance_agfull %>%
  factor(levels = c("NoSignificant","Down","Up"))

Res_em_sub$Significance_intra_s <- Res_em_sub$Significance_intra 
Res_em_sub$Significance_intra_s <- Res_em_sub$Significance_intra_s %>%
  gsub(pattern = "NoSignificant",replacement = NA) %>%
  gsub(pattern = "Significant",replacement = "*")

#Determine the order of the ions clustering them based on their z-score pattern
Tab <- acast(data = Res_em_sub,formula =Ion~Type+Stress,
             value.var = "emmean")
#mclust_ions <- dist(Tab) %>% hclust(method = "ward.D") 
mclust_ions <- hclust(d = 1-cor(Tab %>% t) %>% as.dist,method = "complete")

order_ions <- mclust_ions$order %>% mclust_ions$labels[.]

Res_em_sub$Ion <- Res_em_sub$Ion %>% factor(levels = order_ions)

p_heatmap <- ggplot(data = Res_em_sub,aes(Type,Ion)) +
  geom_raster(aes(fill = emmean)) +  
  geom_tile(aes(color = Significance_agfull),fill = '#00000000', size = 1.5,width = 0.85,height = 0.85) + 
  scale_fill_paletteer_c(package = "pals",palette = "kovesi.diverging_bwr_55_98_c37",
                         limits = c(-2.5,2.5),oob = squish) +
  geom_text(aes(label = Significance_intra_s),size = 8) + 
  geom_text(aes(label = Significance_intra_s),size = 8.2) + 
  geom_text(aes(label = Significance_intra_s),size = 8.4) + 
  
  facet_grid(.~Stress) +
  scale_color_manual(values = c("#00000000","black","black"),na.value = "#00000000")+
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
        #legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(family = "Arial",face = "bold",size = 20),
        axis.text.y = element_text(family = "Arial",face = "bold",hjust = 1,vjust = 0.5),
        axis.text.x = element_text(hjust = 1,vjust = 0.5,family = "Arial",face = "bold"),
        axis.title.x = element_blank()
  )  

#Prepare the dendrogram to put in the same graph
tree <- mclust_ions %>% as.phylo
p_tree_ions <- ggtree(tree,ladderize = F,size = 1.2) + 
  #The expand parth is fundamental to make the tree fit the composition with the other plots
  scale_y_continuous(expand = c(0.026,0.026))


#Create arranged composition
composition <- egg::ggarrange(
  p_tree_ions,p_heatmap,
  ncol = 2,nrow = 1,byrow = T,
  widths = c(0.1,1))


oh.save.pdf(p = composition,outname = "heatmap_ionome_syncom_stresses.pdf",
            outdir = "../figures/",height = 12,width = 14)

write.table(x = Res_em_sub,file = "../cleandata/res_syncomstresses_ionome.tsv",
            append = F,quote = F,sep = "\t",row.names = F,col.names = T)
