library(ohchibi)
library(plyr)

#Merge datasets into a single structure
set.seed(130816)
setwd('/home/isai/Documents/results/rootbarriersmicro/scripts/')


#Check the distributions of the raw data
Dat <- readRDS("../cleandata/dat_41isolatesmonoassociation_cfu_suberin_primroot_dryweight_ionome_zscoreRaw.RDS")


### Now lets do a version with  only 41
df_a <- Dat$Raw$Raw$PrimaryRoot %>%
  subset(Strains!= "NB") %>% subset(Strains != "CDEF") %>% 
  aggregate(zscore ~Strains,.,FUN = mean) 

df_b <- Dat$Raw$Raw$Suberin %>%
  subset(Strains!= "NB") %>% subset(Strains != "CDEF")  %>%
  aggregate(zscore ~Strains,.,FUN = mean) 

df_c <- Dat$Raw$zscore$DryWeight %>%
  subset(Strains!= "NB") %>% subset(Strains != "CDEF") %>%
  aggregate(zscore ~Strains,.,FUN = mean) 


df_d <- Dat$Raw$zscore$CFU  %>%
  subset(Fraction == "Root") %>%
  subset(Strains!= "NB") %>% subset(Strains != "CDEF")  %>%
  aggregate(zscore ~Strains,.,FUN = mean) 


#Remove heat killed from the correlations
df_a <- df_a$Strains %>% grep(pattern = "HK",invert = T) %>% 
  df_a[.,] %>% droplevels
df_b <- df_b$Strains %>% grep(pattern = "HK",invert = T) %>% 
  df_b[.,] %>% droplevels
df_c <- df_c$Strains %>% grep(pattern = "HK",invert = T) %>% 
  df_c[.,] %>% droplevels
df_d <- df_d$Strains %>% grep(pattern = "HK",invert = T) %>% 
  df_d[.,] %>% droplevels



df_a$Variable <- "Root"
df_b$Variable <- "Suberin"
df_c$Variable <- "Dryweight"
df_d$Variable <- "CFU"
colnames(df_d)[2] <- "V1"

merged <- rbind(df_a,df_b) %>% rbind(.,df_c) %>% rbind(.,df_d)

colnames(merged)[2] <- "zscore"
merged$Variable <- merged$Variable %>% factor(levels = c("Root","Dryweight","CFU","Suberin"))

#Compute the variance and create a data frame to annotate the final figure

root_var <- var(df_a$V1)
suberin_var <- var(df_b$V1)
weight_var <- var(df_c$V1)
cfu_var <- var(df_d$V1)

df_labels <- data.frame(Variable = c("Root","Dryweight","CFU","Suberin"),
           label = c(root_var,weight_var,cfu_var,suberin_var))

df_labels$label <- round(df_labels$label,digits = 2)
p <- ggplot(data = merged,aes(  zscore)) +
  geom_density(aes(fill = Variable,color = Variable),alpha = 0.1,size = 1) +
  geom_text(data = df_labels,aes(x = -2,y = 1,label = label),size = 10) +
  theme_ohchibi(size_panel_border = 2) +
  facet_grid(.~Variable) + 
  theme(
    strip.background.x = element_blank(),
    strip.text.x = element_text(family = "Arial",face = "bold",size = 20)
  ) +
  scale_y_continuous(expand = c(0,0)) + 
  ylab(label = "Density") + xlab(label = "z-score") +
  scale_fill_paletteer_d(package = "jcolors",palette = "pal3") +
  scale_color_paletteer_d(package = "jcolors",palette = "pal3") 
  
oh.save.pdf(p = p,outname = "41_variance_variables.pdf",outdir = "../figures/",width = 12,height = 8)

#correlation with cfu
mdf <- Dat$Raw$Raw

ggplot(data = )