library(ohchibi)
library(egg)
library(emmeans)

set.seed(130816)
setwd('/home/isai/Documents/results/rootbarriersmicro/scripts/')

mDat <- readRDS(file = "../cleandata/dat_syncom_genotypes.RDS") 
paleta_type <- mDat$paleta_type[c(1,3)]
res <- readRDS(file = "../cleandata/res_syncom_mutants_suberin_root_weight.RDS")


######### Suberin ######
#####
df_sub <- mDat$Suberin
Res <- res$Res_suberin
df_sub$Genotype <- df_sub$Genotype %>%
  factor(levels = c("Col-0","pCASP1::CDEF","esb1.1pCASP1::CDEF","sgn3.3myb36pELTP::CEDF",
                    "esb1","sgn3","myb36","sgn3myb36",
                    "aba2.1","abi4.1","etr1.1","ein3.1"))
Res$Genotype <- Res$Genotype %>%
  factor(levels = c("Col-0","pCASP1::CDEF","esb1.1pCASP1::CDEF","sgn3.3myb36pELTP::CEDF",
                    "esb1","sgn3","myb36","sgn3myb36",
                    "aba2.1","abi4.1","etr1.1","ein3.1"))

#
##
p1 <- df_sub$Genotype %>% grep(pattern = "CDEF|Col-0|CEDF") %>%
  df_sub[.,] %>% droplevels %>%
  subset(NormSuberin < 4.5) %>%
  chibi.boxplot(Map = .,x_val = "Genotype",col_val = "Type",
                y_val = "NormSuberin",
                mpalette = paleta_type,median_colored_as_points = T,
                size_boxplot = 1,size_median = 4,size_axis_text.x = 25,size_panel_border = 2) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1),
    panel.grid.minor.y = element_blank()
  )  +
  geom_text(data = Res$Genotype %>%  grep(pattern = "CDEF|Col-0|CEDF") %>%
              Res[.,] %>% droplevels,,mapping = aes(x = Genotype,y =  4,label = Letter),size =10)  +
  ylab(label = "Suberin") 


##Weight

df_sub <- mDat$DryWeight
Res <- res$Res_weight
df_sub$Genotype <- df_sub$Genotype %>%
  factor(levels = c("Col-0","pCASP1::CDEF","esb1.1pCASP1::CDEF","sgn3.3myb36pELTP::CEDF",
                    "esb1","sgn3","myb36","sgn3myb36",
                    "aba2.1","abi4.1","etr1.1","ein3.1"))
Res$Genotype <- Res$Genotype %>%
  factor(levels = c("Col-0","pCASP1::CDEF","esb1.1pCASP1::CDEF","sgn3.3myb36pELTP::CEDF",
                    "esb1","sgn3","myb36","sgn3myb36",
                    "aba2.1","abi4.1","etr1.1","ein3.1"))


p2<- df_sub$Genotype %>% grep(pattern = "CDEF|Col-0|CEDF") %>%
  df_sub[.,] %>% droplevels %>%
  chibi.boxplot(Map = .,x_val = "Genotype",col_val = "Type",
                y_val = "Normweight_mg_plant",
                mpalette = paleta_type,median_colored_as_points = T,
                size_boxplot = 1,size_median = 4,size_axis_text.x = 25,size_panel_border = 2) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1),
    panel.grid.minor.y = element_blank()
  )  +
  geom_text(data = Res$Genotype %>%  grep(pattern = "CDEF|Col-0|CEDF") %>%
              Res[.,] %>% droplevels,,mapping = aes(x = Genotype,y =  1.1,label = Letter),size =10)  +
  ylab(label = "Dry weight")  +
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1,1.25))

###### Primary root######
df_sub <- mDat$Root
Res <- res$Res_root

df_sub$Genotype <- df_sub$Genotype %>%
  factor(levels = c("Col-0","pCASP1::CDEF","esb1.1pCASP1::CDEF","sgn3.3myb36pELTP::CEDF",
                    "esb1","sgn3","myb36","sgn3myb36",
                    "aba2.1","abi4.1","etr1.1","ein3.1"))
Res$Genotype <- Res$Genotype %>%
  factor(levels = c("Col-0","pCASP1::CDEF","esb1.1pCASP1::CDEF","sgn3.3myb36pELTP::CEDF",
                    "esb1","sgn3","myb36","sgn3myb36",
                    "aba2.1","abi4.1","etr1.1","ein3.1"))

p3<- df_sub$Genotype %>% grep(pattern = "CDEF|Col-0|CEDF") %>%
  df_sub[.,] %>% droplevels %>%
  chibi.boxplot(Map = .,x_val = "Genotype",col_val = "Type",
                y_val = "NormRoot",
                mpalette = paleta_type,median_colored_as_points = T,
                size_boxplot = 1,size_median = 4,size_axis_text.x = 25,size_panel_border = 2) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1),
    panel.grid.minor.y = element_blank()
  )  +
  geom_text(data = Res$Genotype %>%  grep(pattern = "CDEF|Col-0|CEDF") %>%
              Res[.,] %>% droplevels,,mapping = aes(x = Genotype,y =  1.1,label = Letter),size =10)  +
  ylab(label = "Primary root elongation")  
  



oh.save.pdf(p = p1,outname = "boxplot_suberin_cdef.pdf",outdir = "../figures/",
            width = 18,height = 12)
oh.save.pdf(p = p2,outname = "boxplot_weight_cdef.pdf",outdir = "../figures/",
            width = 18,height = 12)
oh.save.pdf(p = p3,outname = "boxplot_root_cdef.pdf",outdir = "../figures/",
            width = 18,height = 12)
rm(list=ls())
dev.off()
