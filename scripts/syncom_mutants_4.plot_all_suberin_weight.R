library(ohchibi)
library(egg)
library(emmeans)

set.seed(130816)
setwd('/Users/isaisalasgonzalez/Documents/rootbarriersmicro/scripts/')

mDat <- readRDS(file = "../cleandata/dat_syncom_genotypes.RDS") 
paleta_type <- mDat$paleta_type[c(1,3)]
res <- readRDS(file = "../cleandata/res_syncom_mutants_suberin_root_weight.RDS")


######### Suberin ###########
df_sub <- mDat$Suberin
Res <- res$Res_suberin
morder_genos <- c("Col-0","esb1","myb36","sgn3","sgn3myb36",
                  "pCASP1::CDEF","esb1.1pCASP1::CDEF","sgn3.3myb36pELTP::CEDF",
                  "etr1.1","ein3.1","aba2.1","abi4.1"
                  )
df_sub$Genotype <- df_sub$Genotype %>%
  factor(levels =morder_genos)
Res$Genotype <- Res$Genotype %>%
  factor(levels =morder_genos)

##
p1 <- df_sub %>%
  chibi.boxplot(Map = .,x_val = "Genotype",col_val = "Type",
                y_val = "NormSuberin",
                mpalette = paleta_type,median_colored_as_points = T,
                size_boxplot = 1,size_median = 4,size_axis_text.x = 25,size_panel_border = 2) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1),
    panel.grid.minor.y = element_blank(),
  )  +
  geom_text(data = Res,,mapping = aes(x = Genotype,y =  4,label = Letter),size =3)  +
  ylab(label = "Suberin") +
  scale_y_continuous(breaks = seq(0,6,by = 0.5))+
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5),color = "black",size = 0.5, linetype = "longdash")



######### Dey weight ###########
df_sub <- mDat$DryWeight
Res <- res$Res_weight
df_sub$Genotype <- df_sub$Genotype %>%
  factor(levels = morder_genos)
Res$Genotype <- Res$Genotype %>%
  factor(levels = morder_genos)

##
p2 <- df_sub %>%
  chibi.boxplot(Map = .,x_val = "Genotype",col_val = "Type",
                y_val = "Normweight_mg_plant",
                mpalette = paleta_type,median_colored_as_points = T,
                size_boxplot = 1,size_median = 4,size_axis_text.x = 25,size_panel_border = 2) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1),
    panel.grid.minor.y = element_blank(),
  )  +
  geom_text(data = Res,,mapping = aes(x = Genotype,y =  1.1,label = Letter),size =3)  +
  ylab(label = "Dry weight") +
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2))+
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5),color = "black",size = 0.5, linetype = "longdash")



######### Primary root ###########
df_sub <- mDat$Root
Res <- res$Res_root
df_sub$Genotype <- df_sub$Genotype %>%
  factor(levels = morder_genos)
Res$Genotype <- Res$Genotype %>%
  factor(levels = morder_genos)

##
p3 <- df_sub %>%
  chibi.boxplot(Map = .,x_val = "Genotype",col_val = "Type",
                y_val = "NormRoot",
                mpalette = paleta_type,median_colored_as_points = T,
                size_boxplot = 1,size_median = 4,size_axis_text.x = 25,size_panel_border = 2) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1),
    panel.grid.minor.y = element_blank(),
  )  +
  geom_text(data = Res,,mapping = aes(x = Genotype,y =  1.1,label = Letter),size =3)  +
  ylab(label = "Primary root elongation") +
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5),color = "black",size = 0.5, linetype = "longdash")



oh.save.pdf(p = p1,outname = "boxplot_suberin_all.pdf",outdir = "../figures/",
            width = 24,height = 12)
oh.save.pdf(p = p2,outname = "boxplot_weight_all.pdf",outdir = "../figures/",
            width = 24,height = 12)
oh.save.pdf(p = p3,outname = "boxplot_root_al.pdf",outdir = "../figures/",
            width = 24,height = 12)
#rm(list=ls())
#dev.off()

composition <- egg::ggarrange(p1 + theme(axis.text.x = element_blank(),legend.position = "none"),
                              p2 + theme(axis.text.x = element_blank(),legend.position = "none"),
                              p3 + theme(axis.text.x = element_text(size = 10),legend.position = "none"),ncol = 1)
oh.save.pdf(p = composition,outname = "boxplot_all.pdf",outdir = "../figures/",width = 24,height = 12)
