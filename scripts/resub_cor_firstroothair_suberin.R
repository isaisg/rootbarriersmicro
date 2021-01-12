library(ohchibi)
library(ggpubr)
library(Rmisc)
library(ggrepel)
library(scales)
library(multcomp)

#
set.seed(130816)
setwd('/Users/isaisalasgonzalez/Documents/rootbarriersmicro/scripts/')

Dat_a <- readRDS(file = "../cleandata/dat_root_length_vs_number_of_cells_41.RDS")
Dat_b <- readRDS(file = "../cleandata/dat_lignin_xylem_firstroot_diameter_meristem_41.RDS")

df <- Dat_b$Dat_z$Tab %>% t
Map <- match(rownames(df),Dat_b$Dat_z$Map$Id) %>%
  Dat_b$Dat_z$Map[.,]
df_p <- cbind(df,Map[,-4])

#Read the screnning barrier dataset
Dat <- readRDS(file = "../cleandata/dat_screeningbarriers.RDS")
df <- Dat$Suberin
df <- dcast(data = df,formula = Strains + Plant + Batch ~Localization,value.var = "Normvalue_ra")
df$ZoneSum <- df$No_expression + df$Discrete_expression

df$ZoneSum <- df$ZoneSum %>% scale
df_sum_sub <- aggregate(ZoneSum~Strains,df,mean)
colnames(df_sum_sub) <- c("Strain","ZoneSum")
df_sum_frh <- aggregate(DistanceFirstRootHair~Strain,df_p,mean)

mdf <- merge(df_sum_frh,df_sum_sub, by = "Strain")

df <- Dat$Propidium
df$Normvalue <- df$Normvalue %>% scale
df_sum_sub <- aggregate(Normvalue~Strains,df,mean)
mdf$PIScr <- match(mdf$Strain,df_sum_sub$Strains) %>%
  df_sum_sub$V1[.]



Dat <- readRDS(file = "../cleandata/dat_41isolatesmonoassociation_cfu_suberin_primroot_dryweight_ionome.RDS")
df <- Dat$Suberin
df$Norm_Distance_tip_cz <- df$Norm_Distance_tip_cz %>% scale
mtemp <- aggregate(Norm_Distance_tip_cz~Strains,df,mean)
mdf$Dist41 <- match(mdf$Strain,mtemp$Strains) %>%
  mtemp$V1[.]


#Append new PI
Dat_a <- readRDS(file = "../cleandata/dat_root_length_vs_number_of_cells_41.RDS")
df_pi <- Dat_a$df_pi
df_pi$Norm_PI_permeability <- df_pi$Norm_PI_permeability %>% scale
mtemp <- aggregate(Norm_PI_permeability~Strain,df_pi,mean)
mnames <- mtemp$Strain
#Check consistency of the names
#Check consitency of names
df_names <- read.table(file = "../rawdata/41_strains_selected.txt") %$% V1
df_names <- data.frame(Name = df_names,Id = df_names %>%
                         gsub(pattern = "^L|^R",replacement = ""))
df_names <- df_names %>%
  rbind(data.frame(Name = "NB", Id = "NB"))

# 
mtemp$Strain <- match(mtemp$Strain,df_names$Id) %>%
   df_names$Name[.]

mdf$PINew <- match(mdf$Strain,mtemp$Strain) %>%
  mtemp$V1[.]

plot(mdf$PIScr,mdf$PINew)
cor.test(mdf$PIScr,mdf$PINew)

cor.test(mdf$DistanceFirstRootHair,mdf$PIScr)
cor.test(mdf$DistanceFirstRootHair,mdf$PINew)
plot(mdf$DistanceFirstRootHair,mdf$PINew)

plot(mdf$DistanceFirstRootHair,mdf$PIScr)
cor.test(mdf$DistanceFirstRootHair,mdf$PIScr)

#Suberin
cor.test(mdf$ZoneSum,mdf$Dist41)
cor.test(mdf$DistanceFirstRootHair,mdf$ZoneSum)
cor.test(mdf$DistanceFirstRootHair,mdf$Dist41)


#Prepare graph showing correlation
p <- ggplot(data = mdf %>% subset(Strain != "NB"),aes(DistanceFirstRootHair,ZoneSum)) +
  geom_smooth(method = "lm",size = 1,color = "red") +
  geom_point(size = 3) +
  theme_ohchibi(font_family = "Helvetica") +
  scale_fill_manual(values = c("#E68A65","#E6C54E","#59B0E6") %>% rev) +
  scale_color_manual(values = c("#E68A65","#E6C54E","#59B0E6") %>% rev) +
  ylab(label = "Suberin")  + xlab(label = "First root hair") +
  stat_cor(label.x  = 0.5,label.y = -2)
oh.save.pdf(p = p,outname = "res_cor_firstroothair_suberin.pdf",outdir = "../figures/",width = 6,height = 6)
rm(list=ls())
dev.off()
gc()