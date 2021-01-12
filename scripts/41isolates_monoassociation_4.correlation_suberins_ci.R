library(ohchibi)
library(plyr)
library(ggpubr)

#Merge datasets into a single structure
set.seed(130816)
setwd('/home/isai/Documents/results/rootbarriersmicro/scripts/')

Dat <- readRDS("../cleandata/dat_41isolatesmonoassociation_cfu_suberin_primroot_dryweight_ionome.RDS")

df_suberin <- Dat$Suberin
df_suberin$zscore <- df_suberin$NormProp

maxvalue <- df_suberin$zscore %>% quantile(0.99)
minvalue <- df_suberin$zscore %>% quantile(0.01)
df_suberin$zscore[which(df_suberin$zscore > maxvalue)] <- maxvalue
df_suberin$zscore[which(df_suberin$zscore < minvalue)] <- minvalue
sum_suberin <- ddply(df_suberin, c("Strains"), summarise,
                     N    = length(zscore),
                     mean = mean(zscore),
                     sd   = sd(zscore),
                     se   = sd / sqrt(N),
                     ci = (qt(0.95/2 + 0.5, N- 1))*se
                     
)

#Original Suberin
Dat <- readRDS(file = "../cleandata/dat_screeningbarriers.RDS")

df_sub <- Dat$Suberin
melted <- acast(data = df_sub,formula = UId ~Localization,value.var = "Normvalue_ra") 


suma <- melted[,1] + melted[,2]

melted <- cbind(melted,suma)
colnames(melted)[4] <- "SumNoDiscrete"

melted <- melted %>% melt
colnames(melted) <- c("UId","Zone","zscore")
melted$Strains <- match(melted$UId,df_sub$UId) %>% df_sub$Strains[.] %>% as.character
sum_sub_original <- ddply(melted, c("Strains","Zone"), summarise,
                          N    = length(zscore),
                          mean = mean(zscore),
                          sd   = sd(zscore),
                          se   = sd / sqrt(N),
                          ci = (qt(0.95/2 + 0.5, N- 1))*se
                          
)
sum_sub_original_mean <- dcast(data = sum_sub_original,Strains~Zone,value.var = "mean")
sum_sub_original_se <- dcast(data = sum_sub_original,Strains~Zone,value.var = "ci")
sum_sub_original <- merge(sum_sub_original_mean,sum_sub_original_se,by ="Strains")


merged <- merge(sum_sub_original,sum_suberin,by = "Strains")
colnames(merged)
merged <- merged[,c(1,5,9,11,14)]

colnames(merged) <- c("Strains","MeanOriginal","CIOriginal","Mean41","CI41")


cor.test(merged$MeanOriginal,merged$Mean41)
cor.test(merged$MeanOriginal,merged$Mean41,method = "spearman",exact = F)

#Here we need consisten aesthetic as before

#read information about the groups
df_groups <- read.table(file = "../cleandata/df_41_suberin_results.distance.tsv",header = T,sep = "\t")

merged <- merge(merged,df_groups, by = "Strains") 


p1 <- ggplot(data =merged,aes(MeanOriginal,Mean41)) +
  geom_smooth(method = "lm",se = T,color = "black",size = 2) +
  geom_point(size = 10,aes(fill = Group),shape = 21) + theme_ohchibi(size_panel_border = 2) +
  stat_cor() +
  scale_fill_paletteer_d(package ="LaCroixColoR",palette = "PassionFruit" )  +
  scale_x_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1))



##Repeat using the 
Dat <- readRDS("../cleandata/dat_41isolatesmonoassociation_cfu_suberin_primroot_dryweight_ionome.RDS")

df_suberin <- Dat$Suberin
df_suberin$zscore <- df_suberin$Norm_Distance_tip_cz

maxvalue <- df_suberin$zscore %>% quantile(0.99)
minvalue <- df_suberin$zscore %>% quantile(0.01)
df_suberin$zscore[which(df_suberin$zscore > maxvalue)] <- maxvalue
df_suberin$zscore[which(df_suberin$zscore < minvalue)] <- minvalue
sum_suberin <- ddply(df_suberin, c("Strains"), summarise,
                     N    = length(zscore),
                     mean = mean(zscore),
                     sd   = sd(zscore),
                     se   = sd / sqrt(N),
                     ci = (qt(0.95/2 + 0.5, N- 1))*se
                     
)

merged2 <- merge(sum_sub_original,sum_suberin,by = "Strains")
colnames(merged2)
merged2 <- merged2[,c(1,5,9,11,14)]

colnames(merged2) <- c("Strains","MeanOriginal","CIOriginal","Mean41","CI41")

cor.test(merged2$MeanOriginal,merged2$Mean41)
cor.test(merged2$MeanOriginal,merged2$Mean41,method = "spearman",exact = F)
merged2 <- merge(merged2,df_groups, by = "Strains") 


p2 <- ggplot(data =merged2,aes(MeanOriginal,Mean41)) +
  geom_errorbarh(mapping = aes(xmin = MeanOriginal - CIOriginal,
                               xmax = MeanOriginal + CIOriginal),size = 0.1,alpha = 0.2)+
  geom_errorbar(mapping = aes(ymin = Mean41 - CI41,ymax = Mean41 + CI41),size = 0.1,alpha = 0.2)+
  geom_smooth(method = "lm",se = T,color = "black",size = 2) +
  geom_point(size = 10,aes(fill = Group),shape = 21) + theme_ohchibi(size_panel_border = 2) +
  stat_cor() +
  scale_fill_paletteer_d(package ="LaCroixColoR",palette = "PassionFruit" )  +
  scale_x_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1)) 



composition <- egg::ggarrange(p1 + theme(legend.position = "none"),p2 + theme(legend.position = "none"),nrow = 1)

oh.save.pdf(p = p2,outname = "correlation_suberin_old_new.ci.pdf",
            outdir = "../figures/",width = 10,height = 8)
rm(list=ls())
dev.off()
