library(ohchibi)
library(palettesPM)
library(paletteer)
library(emmeans)
library(extrafont)
loadfonts(device = "pdf")
library(egg)
library(multcomp)

#Set random seed
set.seed(130816)


#Read bacterial  dat
setwd('/home/isai/Documents/results/rootbarriersmicro/scripts/')

Dat <- readRDS(file = "../cleandata/dat_censusmutants_16s.RDS")
Dat_rar <- Dat$Rarefied


#Calculate Shannon and Richness
Dat_rar$Map$Shannon <- vegan::diversity(x = Dat_rar$Tab, index = "shannon", MARGIN = 2 )
Dat_rar$Map$Richness <- colSums(Dat_rar$Tab > 0)
Dat_rar$Map$LDepth <- log(Dat_rar$Map$Depth)

########## Root ###############

#Remove bulksoil samples
Dat_sub <- subset.Dataset(x = Dat_rar,subset = Fraction =="Root",drop = T,clean = T) 

###Shannon### anova framework 
m1 <- aov(formula = Shannon ~ Genotype + Rep ,
          data = Dat_sub$Map)

msum <- glht(m1, linfct=mcp(Genotype="Dunnett")) %>% summary 
df_shannon <- emmeans(m1,specs = "Genotype")  %>% as.data.frame
df_shannon$p.value <-  c(1,msum$test$pvalues)



###Richness### anova framework 
m1 <- aov(formula = Richness ~ Genotype + Rep ,
          data = Dat_sub$Map)

msum <- glht(m1, linfct=mcp(Genotype="Dunnett")) %>% summary 
df_richness <- emmeans(m1,specs = "Genotype")  %>% as.data.frame
df_richness$p.value <-  c(1,msum$test$pvalues)


df_root <- data.frame(Genotype = df_shannon$Genotype,
           Shannon = df_shannon$emmean,
           p.value_Shannon = df_shannon$p.value,
           Richness = df_richness$emmean,
           p.value_Richness = df_richness$p.value,
           Fraction = "Root"
           )

#### Shoot ###########

#Remove bulksoil samples
Dat_sub <- subset.Dataset(x = Dat_rar,subset = Fraction =="Shoot",drop = T,clean = T) 

###Shannon### anova framework 
m1 <- aov(formula = Shannon ~ Genotype + Rep ,
          data = Dat_sub$Map)

msum <- glht(m1, linfct=mcp(Genotype="Dunnett")) %>% summary 
df_shannon <- emmeans(m1,specs = "Genotype")  %>% as.data.frame
df_shannon$p.value <-  c(1,msum$test$pvalues)



###Shannon### anova framework 
m1 <- aov(formula = Richness ~ Genotype + Rep ,
          data = Dat_sub$Map)

msum <- glht(m1, linfct=mcp(Genotype="Dunnett")) %>% summary 
df_richness <- emmeans(m1,specs = "Genotype")  %>% as.data.frame
df_richness$p.value <-  c(1,msum$test$pvalues)


df_shoot <- data.frame(Genotype = df_shannon$Genotype,
                      Shannon = df_shannon$emmean,
                      p.value_Shannon = df_shannon$p.value,
                      Richness = df_richness$emmean,
                      p.value_Richness = df_richness$p.value,
                      Fraction = "Shoot"
)

#### Soil ###########

#Remove bulksoil samples
Dat_sub <- subset.Dataset(x = Dat_rar,subset = Fraction =="Soil",drop = T,clean = T) 

###Shannon### anova framework 
m1 <- aov(formula = Shannon ~ Genotype + Rep ,
          data = Dat_sub$Map)

msum <- glht(m1, linfct=mcp(Genotype="Dunnett")) %>% summary 
df_shannon <- emmeans(m1,specs = "Genotype")  %>% as.data.frame
df_shannon$p.value <-  c(1,msum$test$pvalues)



###Shannon### anova framework 
m1 <- aov(formula = Richness ~ Genotype + Rep ,
          data = Dat_sub$Map)

msum <- glht(m1, linfct=mcp(Genotype="Dunnett")) %>% summary 
df_richness <- emmeans(m1,specs = "Genotype")  %>% as.data.frame
df_richness$p.value <-  c(1,msum$test$pvalues)


df_soil <- data.frame(Genotype = df_shannon$Genotype,
                       Shannon = df_shannon$emmean,
                       p.value_Shannon = df_shannon$p.value,
                       Richness = df_richness$emmean,
                       p.value_Richness = df_richness$p.value,
                       Fraction = "Shoot"
)


df <- rbind(df_root,df_shoot,df_soil)

### Put the anova framework ###
Map <- Dat_rar$Map
Map$Fraction <- Map$Fraction %>% factor(levels = c("Soil","Root","Shoot"))
df_num <- read.table(file = "../rawdata/df_genotypes_to_classification.csv",header = T,sep = ",")
df_num$Number <- 1:20
Map <- merge(Map,df_num, by = "Genotype")
Map$Number <- Map$Number %>% factor

p <- chibi.boxplot(Map,x_val = "Number",y_val = "Shannon",facet_formula = "Fraction",style = "open",size_axis_text.x = 15) +
  theme(axis.title.x = element_blank()) 
oh.save.pdf(p = p,outname = "alpha_diversity_amplicon.pdf",outdir = "../figures/",width = 16,height = 10)
