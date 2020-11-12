library(ohchibi)
library(palettesPM)
library(Rmisc)
library(ggpubr)
library(DESeq2)
library(dendextend)
library(multcomp)

set.seed(130816)
setwd('/home/isai/Documents/results/rootbarriersmicro/scripts/')
dir.create("../cleandata")
dir.create("../figures")


Dat <- readRDS(file = "../cleandata/dat_censusmutants_16s.RDS")
#Dat <- readRDS(file = "../cleandata/dat_censusmutants_16s.clean.RDS")


## Determine usable genotypes per fraction and rep
Dat_rar <- Dat$Rarefied

Map_count <- Dat_rar$Map
Map_count$Count <- 1
df_count <- aggregate(Count ~ Genotype+Fraction,Map_count,sum)

#Define Bray Curtis as the dissimilarity
distfun <- function(x,method) vegan::vegdist(x = x, method = "bray")
paleta_geno <- paletteer_d(package = "ggthemes",palette = "Tableau 20")
names(paleta_geno) <- Dat_rar$Map$Genotype %>% levels

mfractions <- Dat_rar$Map$Fraction %>% levels
Res_lm <- NULL
Res_merged <- NULL
Res_p <- list()
mgenos <- Dat_rar$Map$Genotype %>% levels %>% grep(pattern = "Col_0",invert = T,value = T)
for(fraction in mfractions){

  ###### CAP approach ####
  Dat_sub <- Dat_rar %>% subset.Dataset(Fraction == fraction,drop = T,clean = T)
  
  mcap <- oh.cap(Tab = Dat_sub$Tab ,Map = Dat_sub$Map ,
                 formula = "Genotype + Condition(Rep)",
                 distfun = distfun,perms = 9999)
  Mapa_cap <- mcap$Map_cap
  
  df_cap1 <- summarySE(data = Mapa_cap,measurevar = c("CAP1"),groupvars = "Genotype")
  df_cap2 <- summarySE(data = Mapa_cap,measurevar = c("CAP2"),groupvars = "Genotype")
  
  merged_cap <- merge(df_cap1,df_cap2, by = "Genotype")
  merged_cap <- merged_cap[,c(1,3,6,8,11)]
  
  pcap <- ggplot(data = merged_cap,aes(CAP1,CAP2)) + 
    geom_errorbar(aes(ymin = CAP2 - ci.y,ymax = CAP2 + ci.y))+
    geom_errorbarh(aes(xmin = CAP1 - ci.x,xmax = CAP1 + ci.x)) +
    geom_point(shape = 21,aes(fill = Genotype),size = 8) +
    theme_ohchibi() + scale_fill_manual(values = paleta_geno) +
    theme(legend.position = "none")
  Res_p[[fraction]] <- pcap
  maxis <- c("CAP1","CAP2")
  df <- NULL
  for(ax in maxis){
    for(geno in mgenos){
      temp <- Mapa_cap %>% subset(Genotype == geno )
      y <- temp[,colnames(temp) %in% ax]
      if(length(y) ==0){
        next
      }
      temp <- Mapa_cap %>% subset(Genotype == "Col_0" )
      x <- temp[,colnames(temp) %in% ax]
      estimate <- mean(y)-mean(x)
      mlevene <- rbind(data.frame(value = y,geno = "other"),data.frame(value =x,geno = "col")) %>%
        car::leveneTest(value ~geno, data = .) 
      if(mlevene$`Pr(>F)`[1] < 0.05){
        pval <- tryCatch(expr =t.test(x,y,var.equal = F) %$% p.value,error = function(e)p.value = NA )
        df <- data.frame(Fraction = fraction,Genotype = geno,Axis = ax,Estimate = estimate,p.value = pval) %>%
          rbind(df,.)
      }else{
        pval <- tryCatch(expr =t.test(x,y,var.equal = T) %$% p.value,error = function(e)p.value = NA )
        df <- data.frame(Fraction = fraction,Genotype = geno,Axis = ax,Estimate = estimate,p.value = pval) %>%
          rbind(df,.)
      }

    }
  }
  Res_lm <- rbind(Res_lm,df)
  merged_cap$Fraction <- fraction
  Res_merged <- rbind(Res_merged,merged_cap)
}
Res_lm <- merge(Res_lm,df_count,by = c("Genotype","Fraction")) 

res_lm <- Res_lm
res_lm$Direction  <- "Down"
res_lm$Direction[which(res_lm$Estimate > 0)] <- "Up"
res_lm$EstimateaAbs <- abs(res_lm$Estimate)


res_lm$p.adj <- res_lm$p.value %>% p.adjust(method = "bonferroni")
res_lm$Significance <- "NS"
res_lm$Significance[which(res_lm$p.adj < 0.1)] <- "Significant"
res_lm$Significance <- res_lm$Significance %>% factor



res_lm$Direction <- res_lm$Direction %>% factor(levels = c("Down","Up"))
Res_lm <- res_lm


aggregate(Significance ~Fraction,Res_lm,table)

mres <- list(merged_cap = Res_merged,
             Res_lm = Res_lm,
             Res_p = Res_p)

saveRDS(object = mres,file = "../cleandata/res_lm_cap_amplicon_census.RDS")
rm(list=ls())

