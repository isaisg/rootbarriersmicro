library(ohchibi)
library(emmeans)
library(egg)
library(car)
library(ggrepel)
library(paletteer)
library(ggpubr)
#Merge datasets into a single structure
set.seed(130816)
setwd('/home/isai/Documents/results/rootbarriersmicro/scripts/')

Dat <- readRDS("../cleandata/dat_syncom_stresses.RDS")




########## Suberin ################
df_suberin_norm <- Dat$Suberin

df_suberin_norm$Rep <- df_suberin_norm$Rep %>% factor
##Test with anova
m1 <- aov(formula = NormSuberin ~ Stress + Bacteria + Stress:Bacteria,data = df_suberin_norm)

m1_em <- emmeans(m1,pairwise ~Stress:Bacteria,adjust = "none")
df_em <- m1_em$emmeans %>% as.data.frame
df_pval <- m1_em$contrasts %>% as.data.frame
df_pval$contrast <- df_pval$contrast %>%
  gsub(pattern = " ",replacement = "") %>%
  gsub(pattern = ",",replacement = "_")
chosen <- c("-Fe_NB--Fe_SynCom","-Mn_NB--Mn_SynCom","-K_NB--K_SynCom",
            "-Mo_NB--Mo_SynCom","-Zn_NB--Zn_SynCom","-P_NB--P_SynCom",
            "Full_NB-Full_SynCom","+NaCl_NB-+NaCl_SynCom")
df_pval_sub <- which(df_pval$contrast  %in% chosen) %>%
  df_pval[.,] %>% droplevels


df_em$Stress <- df_em$Stress %>% factor %>% relevel(ref = "Full")
df_em$Significance <- "NoSignificant"

df_pval_sub$p.adj <- df_pval_sub$p.value %>% p.adjust(method = "fdr")
st_diff <- df_pval_sub %>% subset(p.adj < 0.05) %$% contrast %>% 
  gsub(pattern = "_.*",replacement = "")
df_em$Significance[which(df_em$Stress %in% st_diff)] <- "Significant"


x <- df_em %>% subset(Bacteria == "NB") %$% emmean
y <- df_em %>% subset(Bacteria == "SynCom") %$% emmean
t.test(x,y)
p1 <- ggplot(df_em,aes(Bacteria,emmean)) +
  #stat_summary(fun.data = "mean_cl_normal",geom = "rect", xmin = 1,xmax = 2,
  #             ,alpha = 0.2,fill = Dat$paleta_type[c(1,3)]) +
  #stat_summary(fun.y = "mean",geom = "rect", xmin = 1,xmax = 2,
  #             aes(ymax = ..y.., ymin = ..y..+0.005),
  #             ,alpha = 0.8,fill = Dat$paleta_type[c(1,3)]) +
  stat_summary(fun.y = "mean",geom = "errorbar",aes(ymax = ..y.., ymin = ..y..),
               width = 1, linetype = "longdash",color = "red")+
  geom_line(aes(group = Stress,color = Significance),size = 2) +
  geom_point(shape = 21,aes(fill = Stress),size = 12) +
  theme_ohchibi(size_panel_border = 2) +
  scale_fill_manual(values = Dat$paleta_stress %>% as.character) +
  scale_color_manual(values = c("#D9D9D9","black")) +
  xlab(label = NULL) + ylab(label = "Suberin") +
  theme(
   # legend.position = "none"
  ) +

  theme(
   # panel.grid.major.y = element_line(size = 0.5,color = "#D9D9D9"),
  #  panel.grid.minor.y = element_line(size = 0.5,color = "#D9D9D9"),
    panel.background = element_rect(fill = "white"),
    strip.text.x = element_text(family = "Arial",face = "bold",size = 20),
    axis.title.x = element_blank()
  )   



df_em %>% subset(Bacteria == "NB") %$% emmean %>% var
df_em %>% subset(Bacteria == "SynCom") %$% emmean %>% var

## sup figure ##
p <- chibi.boxplot(Map = df_suberin_norm,x_val = "Bacteria",
              y_val = "NormSuberin",style = "open",facet_formula = "Stress")

oh.save.pdf(p = p,outname = "suberin_syncom_stress_nb_syncom.sup.pdf",outdir = "../figures/",width = 18,height = 9)

####### Try another implementation ######3
p <- chibi.boxplot(Map = df_suberin_norm,x_val = "Bacteria",
                   y_val = "NormSuberin",style = "open",facet_formula = "Stress")

paleta_type <- Dat$paleta_type[c(1,3)]

p_col <- chibi.boxplot(Map = df_suberin_norm,x_val = "Bacteria",col_val = "Bacteria",
              y_val = "NormSuberin",facet_formula = "Stress",
              mpalette = paleta_type,median_colored_as_points = T,
              size_boxplot = 1,size_median = 4,size_axis_text.x = 25,size_panel_border = 2) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1),
    panel.grid.minor.y = element_blank()
  ) 

oh.save.pdf(p = p_col,outname = "suberin_syncom_stress_nb_syncom.supcolor.pdf",outdir = "../figures/",width = 18,height = 9)





########## Dry weight ############
df_suberin_norm <- Dat$DryWeight
df_suberin_norm$Replicates <- df_suberin_norm$Replicates %>% factor
m1 <- aov(formula = Normweight_mg_plant ~ Stress + Bacteria + Stress:Bacteria,
          data = df_suberin_norm)
m1_em <- emmeans(m1,pairwise ~Stress:Bacteria,adjust = "none")
df_em <- m1_em$emmeans %>% as.data.frame
df_pval <- m1_em$contrasts %>% as.data.frame
df_pval$contrast <- df_pval$contrast %>%
  gsub(pattern = " ",replacement = "") %>%
  gsub(pattern = ",",replacement = "_")
chosen <- c("-Fe_NB--Fe_SynCom","-Mn_NB--Mn_SynCom","-K_NB--K_SynCom",
            "-Mo_NB--Mo_SynCom","-Zn_NB--Zn_SynCom","-Pi_NB--Pi_SynCom",
            "Full_NB-Full_SynCom","+NaCl_NB-+NaCl_SynCom")
df_pval_sub <- which(df_pval$contrast  %in% chosen) %>%
  df_pval[.,] %>% droplevels

df_em$Stress <- df_em$Stress %>% factor %>% relevel(ref = "Full")
df_em$Significance <- "NoSignificant"

df_pval_sub$p.adj <- df_pval_sub$p.value %>% p.adjust(method = "fdr")
st_diff <- df_pval_sub %>% subset(p.adj < 0.05) %$% contrast %>% 
  gsub(pattern = "_.*",replacement = "")
df_em$Significance[which(df_em$Stress %in% st_diff)] <- "Significant"

x <- df_em %>% subset(Bacteria == "NB") %$% emmean
y <- df_em %>% subset(Bacteria == "SynCom") %$% emmean
t.test(x,y)


p2 <- ggplot(df_em,aes(Bacteria,emmean)) +
  stat_summary(fun.y = "mean",geom = "errorbar",aes(ymax = ..y.., ymin = ..y..),
               width = 1, linetype = "longdash",color = "red")+
  geom_line(aes(group = Stress,color = Significance),size = 2) +
  geom_point(shape = 21,aes(fill = Stress),size = 12) +
  theme_ohchibi(size_panel_border = 2) +
  scale_fill_manual(values = Dat$paleta_stress %>% as.character) +
  scale_color_manual(values = c("#D9D9D9","black")) +
  xlab(label = NULL) + ylab(label = "Dry weight")+
  theme(
   # panel.grid.major.y = element_line(size = 0.5,color = "#D9D9D9"),
    #panel.grid.minor.y = element_line(size = 0.5,color = "#D9D9D9"),
    panel.background = element_rect(fill = "white"),
    strip.text.x = element_text(family = "Arial",face = "bold",size = 20),
    axis.title.x = element_blank(),
    #legend.position = "none"
  )   + scale_y_continuous(breaks = c(0.2,0.4,0.6,0.8,1,1.2))

df_em %>% subset(Bacteria == "NB") %$% emmean %>% var
df_em %>% subset(Bacteria == "SynCom") %$% emmean %>% var


##### sup figure ###
p <- chibi.boxplot(Map = df_suberin_norm,x_val = "Bacteria",
                   y_val = "Normweight_mg_plant",style = "open",facet_formula = "Stress")

oh.save.pdf(p = p,outname = "dryweight_syncom_stress_nb_syncom.sup.pdf",outdir = "../figures/",width = 18,height = 9)



p_col <- chibi.boxplot(Map = df_suberin_norm,x_val = "Bacteria",col_val = "Bacteria",
                       y_val = "Normweight_mg_plant",facet_formula = "Stress",
                       mpalette = paleta_type,median_colored_as_points = T,
                       size_boxplot = 1,size_median = 4,size_axis_text.x = 25,size_panel_border = 2) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1),
    panel.grid.minor.y = element_blank()
  ) 
oh.save.pdf(p = p_col,outname = "dryweight_syncom_stress_nb_syncom.supcol.pdf",outdir = "../figures/",width = 18,height = 9)



#Correlatoin
df_sub <- Dat$Suberin

df_dw <- Dat$DryWeight %>% subset(Plates != "5") %>%
  droplevels


df_sub$Plate <- rep(c(1:4),40)
df_sub <- df_sub[,c(1,2,6,8,11)]

df_sub$Stress <- df_sub$Stress %>% gsub(pattern = " ",replacement = "")
df_sub$UId <- paste0(df_sub$Bacteria,"|",df_sub$Stress,"|",
                     df_sub$Rep,"|",df_sub$Plate)

df_dw$Stress <- df_dw$Stress %>% gsub(pattern = "Pi",replacement = "P")
df_dw$UId <- paste0(df_dw$Bacteria,"|",df_dw$Stress,"|",df_dw$Replicates,"|",df_dw$Plates)


mint <- intersect(df_dw$UId %>% unique,df_sub$UId %>% unique)

df_dw_merged <- match(df_sub$UId,df_dw$UId) %>%
  df_dw[.,]

df_sub$Normweight_mg_plant <- df_dw_merged$Normweight_mg_plant
df_sub <- df_sub %>% na.omit 

#df_sub <- df_sub %>% subset(Stress != "-P") %>% droplevels

df_sub$Stress <- df_sub$Stress %>% factor %>% relevel(ref = "Full")

#calculate means to add lines
df_mean_sub <- aggregate(NormSuberin ~Bacteria,df_sub,mean) 
df_mean_dw <-  aggregate(Normweight_mg_plant ~Bacteria,df_sub,mean) 

p3 <- ggplot(data =df_sub,aes(NormSuberin,Normweight_mg_plant)) +
  geom_hline(data = df_mean_dw,aes(yintercept = Normweight_mg_plant),
             color = "red",size = 1,linetype = "longdash",alpha = 0.5) +
  geom_vline(data = df_mean_sub,aes(xintercept = NormSuberin),
             color = "red",size = 1,linetype = "longdash",alpha = 0.5) +
  geom_smooth(method = "lm",se = FALSE,color = "black",size = 4) +
  geom_point(shape = 21,size = 10,aes(fill = Stress)) +
  scale_fill_manual(values = Dat$paleta_stress) +
  facet_grid(.~Bacteria) +
  stat_cor(size = 6) + 
  theme_ohchibi(size_panel_border = 2) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_text(family = "Arial",face = "bold",size = 20)
  ) +
  xlab(label = "Suberin") + ylab(label = "Dry Weight") +
  scale_y_continuous(breaks = c(0.2,0.4,0.6,0.8,1,1.2,1.4,1.6)) +
  theme(
    panel.grid.major.y = element_line(size = 0.5,color = "#D9D9D9"),
  #  panel.grid.minor.y = element_line(size = 0.5,color = "#D9D9D9"),
    panel.grid.major.x = element_line(size = 0.5,color = "#D9D9D9"),

    panel.background = element_rect(fill = "white"),
  )
#composition <- egg::ggarrange(p1,p2,nrow = 1)



oh.save.pdf(p = p1,outname = "suberin_syncom_stress_nb_syncom.pdf",outdir = "../figures/",width = 11,height = 12)

oh.save.pdf(p = p2,outname = "dryweight_syncom_stress_nb_syncom.pdf",outdir = "../figures/",width = 11,height = 12)

oh.save.pdf(p = p3,outname = "correlation_suberin_dryweight_syncom_stresses.pdf",
            outdir = "../figures/",width = 26,height = 10)
