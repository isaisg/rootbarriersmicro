library(ohchibi)
library(palettesPM)
library(Rmisc)
library(ggpubr)
library(DESeq2)
library(dendextend)
library(scales)
library(ggrepel)
set.seed(130816)
setwd('/Users/isaisalasgonzalez/Documents/rootbarriersmicro/scripts/')
dir.create("../cleandata")
dir.create("../figures")

#### Read LM results ####
Res <- readRDS(file = "../cleandata/res_lm_cap_amplicon_census.RDS")

### Correlations ######
merged_cap <- Res$merged_cap
paleta_geno <- paletteer_d(package = "ggthemes",palette = "Tableau 20")
names(paleta_geno) <- merged_cap$Genotype %>% levels

merged_cap <- merged_cap %>% subset(Genotype != "dir9_dir18_esb1_1") %>% droplevels




df <- merged_cap %>% subset(Fraction == "Root")
Tab_root <- df[,c(2,4)] %>% as.matrix
rownames(Tab_root) <- df$Genotype


df <- merged_cap %>% subset(Fraction == "Shoot")
Tab_shoot <- df[,c(2,4)] %>% as.matrix
rownames(Tab_shoot) <- df$Genotype


df <- merged_cap %>% subset(Fraction == "Soil")
Tab_soil <- df[,c(2,4)] %>% as.matrix
rownames(Tab_soil) <- df$Genotype


mdist_root <- Tab_root %>% dist
mdist_shoot <- Tab_shoot %>% dist
mdist_soil <- Tab_soil %>% dist


mclust_root <-mdist_root %>% hclust
mclust_shoot <- mdist_shoot %>% hclust
mclust_soil <- mdist_soil %>% hclust


mmantel_root_shoot <- mantel(mdist_root,mdist_shoot,permutations = 9999)
mmantel_root_soil <- mantel(mdist_root,mdist_soil,permutations = 9999)
mmantel_shoot_soil <- mantel(mdist_shoot,mdist_soil,permutations = 9999)

### Plot correlations ###
merged_mantel <- mdist_root %>% as.matrix %>% melt_dist() %>%
  merge(
    mdist_shoot %>% as.matrix %>% melt_dist(),
    by = c("iso1","iso2")
      
  )%>%
  merge(
    mdist_soil %>% as.matrix %>% melt_dist(),
    by = c("iso1","iso2")
    
  )
colnames(merged_mantel)[3:5] <- c("Root","Shoot","Soil")

mtext <- paste0("r = ",round(mmantel_root_soil$statistic,3),"\np = ",format.pval(mmantel_root_soil$signif))
p1 <- ggplot(merged_mantel,aes(Soil,Root)) +
  geom_smooth(method = "lm",color = "red",size = 2,fill = "#D9D9D9") +
  geom_point(size = 3) +
  theme_ohchibi(size_panel_border = 2) + 
  annotate(geom = "text",x = 3,y = 0.6,label = mtext)+
  scale_y_continuous(limits = c(0,3.5)) +
  scale_x_continuous(limits = c(0,3.5)) +
  coord_equal() 


mtext <- paste0("r = ",round(mmantel_shoot_soil$statistic,3),"\np = ",format.pval(mmantel_shoot_soil$signif))
p2 <- ggplot(merged_mantel,aes(Soil,Shoot)) +
  geom_smooth(method = "lm",color = "red",size = 2,fill = "#D9D9D9") +
  geom_point(size = 3) +
  theme_ohchibi(size_panel_border = 2) + 
  coord_equal() +
  annotate(geom = "text",x = 3,y = 0.6,label = mtext)+
  scale_y_continuous(limits = c(0,3.5)) +
  scale_x_continuous(limits = c(0,3.5)) +
  coord_equal() 


mtext <- paste0("r = ",round(mmantel_root_shoot$statistic,3),"\np = ",format.pval(mmantel_root_shoot$signif))
p3 <- ggplot(merged_mantel,aes(Root,Shoot)) +
  geom_smooth(method = "lm",color = "red",size = 2,fill = "#D9D9D9") +
  geom_point(size =3) +
  theme_ohchibi(size_panel_border = 2) + 
  coord_equal() +
  annotate(geom = "text",x = 3,y = 0.6,label = mtext) +
  scale_y_continuous(limits = c(0,3.5)) +
  scale_x_continuous(limits = c(0,3.5)) +
  coord_equal() 



#composition <- egg::ggarrange(p1,p2,p3,nrow = 1)


################ Do cap constellations ################
merged_cap <- Res$merged_cap
Res_lm <- Res$Res_lm
df_sum <- dcast(data = Res_lm,formula = Fraction+Genotype~Axis,value.var = "Significance") 
df_sum$Pattern <- "NS"
df_sum$Pattern[which((df_sum$CAP1 == "Significant") & (df_sum$CAP2 == "Significant"))] <- "Both"
df_sum$Pattern[which((df_sum$CAP1 == "Significant") & (df_sum$CAP2 != "Significant"))] <- "First"
df_sum$Pattern[which((df_sum$CAP1 != "Significant") & (df_sum$CAP2 == "Significant"))] <- "Second"

df_sum$Pattern %>% table
head(df_sum)


merged_cap <- merge(merged_cap,df_sum[,c(1,2,5)], by = c("Genotype","Fraction"),all.x = TRUE) 

df_num <- read.table(file = "../rawdata/df_genotypes_to_classification.csv",header = T,sep = ",")
df_num$Number <- 1:20
merged_cap <- merge(merged_cap,df_num, by = "Genotype")


merged_cap$PatternEasy <- merged_cap$Pattern %>% 
  gsub(pattern = "Second|First|Both",replacement = "Significant") 
merged_cap$PatternEasy[which(is.na(merged_cap$PatternEasy))] <- "Col_0"
merged_cap$PatternEasy <- merged_cap$PatternEasy %>% factor(levels = c("Col_0","NS","Significant"))


p_all <- ggplot(data = merged_cap,aes(CAP1,CAP2)) + 
geom_vline(xintercept = 0,linetype = "longdash",size =2,color = "#D9D9D9") +
  geom_hline(yintercept = 0,linetype = "longdash",size =2,color = "#D9D9D9") +
  geom_linerange(aes(ymin = CAP2 - ci.y,ymax = CAP2 + ci.y),size = 0.1)+
  geom_linerange(aes(xmin = CAP1 - ci.x,xmax = CAP1 + ci.x),size = 0.1) +
  geom_text_repel(aes(color = PatternEasy,label = Number),size = 8) +
  theme_ohchibi() + 
  facet_grid(.~Fraction) +
  scale_x_continuous(limits = c(-2.75,2.75),oob = rescale_none)+
  scale_y_continuous(limits = c(-3,3),oob = rescale_none) +
  scale_shape_manual(values = 21:24)+
  theme(
    legend.position = "none",
    strip.background.x = element_blank(),
    strip.text.x = element_text(family = "Arial",face = "bold",size = 25) 
  )  +
  scale_color_manual(values = c("darkblue","black","red"))


p_root <- ggplot(data = merged_cap %>% subset(Fraction == "Root"),aes(CAP1,CAP2)) + 
  geom_vline(xintercept = 0,linetype = "longdash",size =2,color = "#D9D9D9") +
  geom_hline(yintercept = 0,linetype = "longdash",size =2,color = "#D9D9D9") +
  geom_linerange(aes(ymin = CAP2 - ci.y,ymax = CAP2 + ci.y),size = 0.1,alpha = 0.2)+
  geom_linerange(aes(xmin = CAP1 - ci.x,xmax = CAP1 + ci.x),size = 0.1,alpha = 0.2) +
  geom_text_repel(aes(color = PatternEasy,label = Number),size = 8) +
  theme_ohchibi(size_panel_border = 2) + 
  facet_grid(.~Fraction) +
  scale_x_continuous(limits = c(-2.75,2.75),oob = rescale_none)+
  scale_y_continuous(limits = c(-3,3),oob = rescale_none) +
  scale_shape_manual(values = 21:24)+
  theme(
    legend.position = "none",
    strip.background.x = element_blank(),
    strip.text.x = element_text(family = "Arial",face = "bold",size = 25) 
  )  +
  scale_color_manual(values = c("darkblue","black","red"))
#coord_equal()

p_shoot <- ggplot(data = merged_cap %>% subset(Fraction == "Shoot"),aes(CAP1,CAP2)) + 
  geom_vline(xintercept = 0,linetype = "longdash",size =2,color = "#D9D9D9") +
  geom_hline(yintercept = 0,linetype = "longdash",size =2,color = "#D9D9D9") +
  geom_linerange(aes(ymin = CAP2 - ci.y,ymax = CAP2 + ci.y),size = 0.1,alpha = 0.2)+
  geom_linerange(aes(xmin = CAP1 - ci.x,xmax = CAP1 + ci.x),size = 0.1,alpha = 0.2) +
  geom_text_repel(aes(color = PatternEasy,label = Number),size = 8) +
  theme_ohchibi(size_panel_border = 2) + 
  facet_grid(.~Fraction) +
  scale_x_continuous(limits = c(-2.75,2.75),oob = rescale_none)+
  scale_y_continuous(limits = c(-3,3),oob = rescale_none) +
  scale_shape_manual(values = 21:24)+
  theme(
    legend.position = "none",
    strip.background.x = element_blank(),
    strip.text.x = element_text(family = "Arial",face = "bold",size = 25) 
  )  +
  scale_color_manual(values = c("darkblue","black","red"))
  #coord_equal()


p_soil <- ggplot(data = merged_cap %>% subset(Fraction == "Soil"),aes(CAP1,CAP2)) + 
  geom_vline(xintercept = 0,linetype = "longdash",size =2,color = "#D9D9D9") +
  geom_hline(yintercept = 0,linetype = "longdash",size =2,color = "#D9D9D9") +
  geom_linerange(aes(ymin = CAP2 - ci.y,ymax = CAP2 + ci.y),size = 0.1,alpha = 0.2)+
  geom_linerange(aes(xmin = CAP1 - ci.x,xmax = CAP1 + ci.x),size = 0.1,alpha = 0.2) +
  geom_text_repel(aes(color = PatternEasy,label = Number),size = 8) +
  theme_ohchibi(size_panel_border = 2) + 
  facet_grid(.~Fraction) +
  scale_x_continuous(limits = c(-2.75,2.75),oob = rescale_none)+
  scale_y_continuous(limits = c(-3,3),oob = rescale_none) +
  scale_shape_manual(values = 21:24)+
  theme(
    legend.position = "none",
    strip.background.x = element_blank(),
    strip.text.x = element_text(family = "Arial",face = "bold",size = 25) 
  )  +
  scale_color_manual(values = c("darkblue","black","red"))
  #coord_equal()

composition <- egg::ggarrange(p_root,p_shoot,p_soil,
  p1,p2,p3,nrow = 2,ncol = 3)
oh.save.pdf(p = composition,outname = "composition_correlation_amplicon_census_between_fractions_caps.pdf",
            outdir = "../figures/",width = 20,height = 14)

### save plots
mplots <- list(
  constellation_root_amplicon = p_root,
  constellation_shoot_amplicon = p_shoot,
  constellation_soil_amplicon = p_soil,
  correlation_amplicon_soil_root = p1,
  correlation_amplicon_soil_shoot = p2,
  correlation_amplicon_root_shoot = p3
  
  
)
saveRDS(object = mplots,
        file = "../cleandata/plots_constellations_correlations_amplicon_wild.RDS")

df_num

rm(list=ls())
dev.off()
gc()
