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
library(ggpubr)

set.seed(130816)
setwd('/home/isai/Documents/results/rootbarriersmicro/scripts/')


Dat <- readRDS("../cleandata/dat_41isolatesmonoassociation_cfu_suberin_primroot_dryweight_ionome_zscoreRaw.distance.RDS")
df <- Dat$Merged$zscore

#select 41 tested
st <- read.table(file = "../rawdata/41_isolates.tsv") %$% V2 %>% as.character

df <- which(df$Strains %in% c(st,"NB")) %>% df[.,] %>% droplevels
colnames(df)

df <- df[,c(1,12,13,24:57)]


#add the group information
df_groups <- read.table(file = "../cleandata/df_41_suberin_results.distance.tsv",header = T,sep = "\t")
merged <- merge(df,df_groups, by = "Strains",all.x = T)
merged$Group <- merged$Group %>% as.character
merged$Group[which(is.na(merged$Group))] <- "HK"
merged$Group[merged$Strains %>% grep(pattern = "NB")] <- "NB"

merged$Group <- merged$Group %>% 
  factor(levels = c("NB","HK","Group1","Group2","Group3","Group4","Group5"))



paleta_alive <- c("#414141","#D9D9D9",paletteer_d(package = "LaCroixColoR",palette = "PassionFruit",n = 5))
names(paleta_alive) <-  c("NB","HK","Group1","Group2","Group3","Group4","Group5")

df <- merged

mions <- colnames(df) %>% grep(pattern = "Strains|Dist|Group",invert = T,value = T) %>%
  grep(pattern = "Mean",value = T) %>% gsub(pattern = "Mean",replacement = "")

mplots <- list()
Res <- NULL
for(ion in mions){
  mm <- paste0("Mean",ion)
  mc <- paste0("CI",ion)
  targets <- c("Strains","MeanDistSub","CIDistSub",mm,mc,"Group")
  df_temp <- which(colnames(df) %in% targets) %>%
    df[,.]
  df_temp$Down <- df_temp[,4]-df_temp[,5]
  df_temp$Up<- df_temp[,4]+df_temp[,5]
  
  p <- ggplot(data =df_temp,aes_string(x = "MeanDistSub",mm) ) +
    geom_smooth(method = "lm",se = T,size = 2,color = "black") +
    geom_errorbarh(mapping = aes(xmin =MeanDistSub - CIDistSub,xmax = MeanDistSub + CIDistSub, ),size = 0.5,alpha = 0.1) +
    geom_errorbar(mapping = aes(ymin =Down,ymax = Up),size = 0.5,alpha = 0.1)  + 
    geom_point(shape = 21,size = 7,aes(fill = Group),stroke = 1) + 
    theme_ohchibi(size_panel_border = 2) + 
    scale_fill_manual(values = paleta_alive) +
    scale_color_manual(values = paleta_alive) +
    xlab(label = "Suberin z-score") +
    ylab(label = paste0(ion," z-score")) +
    theme(legend.position = "none") +
    stat_cor() 
  mplots[[ion]] <- p
  m1 <- cor.test(df_temp$MeanDistSub,df_temp[,4])
  Res <- data.frame(Ion = ion,m1$estimate ,m1$p.value,row.names = NULL) %>%
    rbind(Res,.)
}

Res <- with(Res,order(m1.estimate)) %>%
  Res[.,]
Res$p.adj <- Res$m1.p.value %>% p.adjust(method = "bonferroni")
Res
write.table(x = Res,file = "../cleandata/res_41_correlations_ion_suberin.distance.tsv",
            append = F,quote = F,sep = "\t",row.names = F,col.names = T)

composition <- egg::ggarrange(mplots$B11,mplots$Na23,mplots$Co59,mplots$Sr88,
                              mplots$Mg24,mplots$Zn66,mplots$Mo98,mplots$P31,
                              mplots$As75,mplots$Fe56,mplots$Cd111,mplots$Rb85,
                              mplots$Cu63,mplots$Ca43,mplots$S34,mplots$Mn55,
                              mplots$K39,
                              nrow = 5,ncol = 4)
oh.save.pdf(p = composition,outname = "corr_suberin_ions.distance.pdf",outdir = "../figures/",width = 32,height = 32)

############ Summarized version of the figure #################
Res$Ion <- Res$Ion %>% gsub(pattern = "[0-9]+",replacement = "")
morder <- Res$Ion %>% as.character
Res$Significance <- "No"
Res$Significance[which(Res$p.adj < 0.05)] <- "Yes"
Res$Ion <- Res$Ion %>% factor(levels = morder)


p <- ggplot(data = Res,aes(Ion,m1.estimate)) +
  geom_bar(stat = "identity",aes(fill = Significance),color = "black", size = 1) +
  theme_ohchibi(size_panel_border = 2) +
  ylab(label = "Pearson r") +
  scale_y_continuous(breaks = seq(-1,1,0.2)) +
  scale_fill_manual(values = c("#D9D9D9","red")) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    panel.grid.major.y = element_line(size = 0.5,color = "#D9D9D9"),
    #panel.grid.major. = element_line(size = 0.5,color = "#D9D9D9")
    
  ) 
oh.save.pdf(p = p,outname = "corr_suberin_ions.distance.minimal.pdf",outdir = "../figures/",width = 8,height = 4)
