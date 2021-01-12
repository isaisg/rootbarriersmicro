library(ohchibi)
library(egg)
library(emmeans)
library(multcomp)
library(car)

set.seed(130816)
setwd('/Users/isaisalasgonzalez/Documents/rootbarriersmicro/scripts/')

## Load the number scheme
df_num <- read.table(file = "../rawdata/df_genotypes_to_classification.csv",header = T,sep = ",")
df_num$Number <- 1:20
df_num$Genotype

#### Natural soil data 11 weeks ####
df_wild11 <- read.table(file = "../rawdata/shootarea_wildsoilds_week11.csv",
                        header = T,sep = ",",quote = "",comment.char = "")
df_wild11$genotype <- df_wild11$genotype %>%
  gsub(pattern = "^C4H::F5H$",replacement = "C4H_FSH") %>%
  gsub(pattern = "^casp1.1 casp3.1$",replacement = "casp1_1_casp3_1") %>%
  gsub(pattern = "Col-0",replacement = "Col_0") %>%
  gsub(pattern = "^dir9-1 dir18-1 esb1-1$",replacement = "dir9_dir18_esb1_1") %>%
  gsub(pattern = "^erk1-3$",replacement = "Kinase_salk_060966") %>%
  gsub(pattern = "erk1-3 rbk1-1",replacement = "Kinase_salk_060966_Kinase_salk_043441") %>%
  gsub(pattern = "^esb1-1$",replacement = "esb1_1") %>%
  gsub(pattern = "^esb1-1 sgn3-3$",replacement = "esb1_1sng3_3") %>%
  gsub(pattern = "^horst-1$",replacement = "horst_1") %>%
  gsub(pattern = "^myb36-2$",replacement = "myb36_2") %>%
  gsub(pattern = "^myb36-2 sgn3-3$",replacement = "myb36_2_sgn3_3") %>%
  gsub(pattern = "pCASP1::CDEF1\\(esb1-1\\)",replacement = "pCASP1_CDEFesb1_1") %>%
  gsub(pattern = "pCASP1::CDEF1\\(wild-type\\)",replacement = "pCASP1_CDEFwt") %>%
  gsub(pattern = "pELTP::CDEF1\\(myb36-2 sgn3-3\\)",replacement = "pELTP_CDEFsgn3_3_myb36_1") %>%
  gsub(pattern = "pELTP::CDEF1\\(sgn3-3\\)",replacement = "pELTP_CDEFsgn3_3") %>%
  gsub(pattern = "^ralph-1$",replacement = "ralph_1") %>%
  gsub(pattern = "ralph-1 horst-1",replacement = "ralph_1_horst_1") %>%
  gsub(pattern = "^rbk1-1$",replacement = "Kinase_salk_043441") %>%
  gsub(pattern = "^sgn3-3$",replacement = "sgn3_3") %>%
  gsub(pattern = "tic-2",replacement = "tic")


# Append the information of number and grouping
df_wild11$NumGenotype <- match(df_wild11$genotype,df_num$Genotype) %>%
  df_num$Number[.]
df_wild11$GroupGenotype <- match(df_wild11$genotype,df_num$Genotype) %>%
  df_num$Group[.]

df_wild11$NumGenotype <- df_wild11$NumGenotype %>% factor
df_wild11$GroupGenotype <- df_wild11$GroupGenotype %>% factor


#Scale the area so the values can be compared to agar
mgenos <- df_wild11$genotype %>% unique %>%
  grep(pattern = "Col_0",invert = T,value = T)
df_col0 <- df_wild11 %>% subset(genotype == "Col_0")

mcol <- df_col0$Area_plant %>% mean
Res_11 <- NULL
for(mg in mgenos){
  mtemp_11 <- df_wild11 %>% subset(genotype == mg) %>% droplevels
  #Levene test
  df_temp <- rbind(df_col0,mtemp_11)
  df_temp$genotype <- df_temp$genotype %>% factor %>% relevel(ref = "Col_0")
  m <- car::leveneTest(as.numeric(Area_plant) ~ genotype,data = df_temp)
  pval <- m$`Pr(>F)`[1]
  x <- df_col0$Area_plant
  y <- mtemp_11$Area_plant 
  mmean <- y %>% mean
  mest <- mmean - mcol
  if(pval < 0.05){
    m0 <- t.test(x = x,y = y,exact = F,var.equal = F)
    
  }else{
    m0 <- t.test(x = x,y = y,exact = F,var.equal = T)
    
  }
  m0 <- t.test(x = df_col0$Area_plant,y = mtemp_11$Area_plant)
  mpv <- m0$p.value
  Res_11 <- data.frame(genotype = mg,coefficients = mest,pvalue = mpv) %>%
    rbind(Res_11,.)
}
Res_11$padj <- Res_11$pvalue %>% p.adjust(method = "fdr")

Res_11$NumGenotype <- match(Res_11$genotype,df_num$Genotype) %>%
  df_num$Number[.]
Res_11$GroupGenotype <- match(Res_11$genotype,df_num$Genotype) %>%
  df_num$Group[.]
Res_11$Significance <- NA
Res_11$Significance[which(Res_11$padj < 0.1)] <- "q < 0.1"

df_wild11$Significance <- match(df_wild11$genotype,Res_11$genotype) %>%
  Res_11$Significance[.]
df_wild11$Significance <- df_wild11$Significance %>% factor

### Plot for the soil data
lims_col <- df_wild11 %>% subset(genotype == "Col_0") %$% Area_plant %>%
  quantile()

p2 <- ggplot(data = df_wild11,aes(NumGenotype,Area_plant)) + 
  geom_sina(alpha = 0.3,shape = 21,size = 4) +
  stat_summary(mapping = aes(color  = Significance),
               fun.ymin = function(z) { quantile(z,0.25) },
               fun.ymax = function(z) { quantile(z,0.75) },
               fun.y = median,size = 0.75) +
  #facet_grid(.~GroupGenotype,space = "free",scales = "free") +
  theme_ohchibi() + 
  theme(
    strip.background.x = element_blank(),
    strip.text.x = element_text(family = "Helvetica",size = 15)
  ) +
  scale_color_manual(values = c("red"),na.value = "black") +
  theme(
    axis.title.x = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(family = "Helvetica",face = "plain"),
    axis.text.y = element_text(family = "Helvetica",face = "plain"),
    axis.title.y = element_text(family = "Helvetica",face = "plain",size = 20)
  ) +
  ylab(label = "Shoot circumference (cm)") +
  geom_rect(mapping = aes(ymin = lims_col[2],ymax = lims_col[4],xmin = -Inf,xmax =Inf),
            color = NA,fill = "#A6CEE3",alpha = 0.1) + 
  scale_y_continuous(breaks = seq(0,30,5))


# oh.save.pdf(p = p2,outname = "resub_soil_shootcircum_week11.pdf",outdir = "../figures/",
#             width = 10,height = 10)


p2 <- ggplot(data = df_wild11,aes(NumGenotype,Area_plant)) + 
  geom_sina(alpha = 0.3,shape = 21,size = 4) +
  stat_summary(fun.data = mean_cl_normal,aes(color = Significance)) +
  #facet_grid(.~GroupGenotype,space = "free",scales = "free") +
  theme_ohchibi() + 
  theme(
    strip.background.x = element_blank(),
    strip.text.x = element_text(family = "Helvetica",size = 15)
  ) +
  scale_color_manual(values = c("red"),na.value = "black") +
  theme(
    axis.title.x = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(family = "Helvetica",face = "plain"),
    axis.text.y = element_text(family = "Helvetica",face = "plain"),
    axis.title.y = element_text(family = "Helvetica",face = "plain",size = 20)
  ) +
  ylab(label = "Shoot circumference (cm)") +
  scale_y_continuous(breaks = seq(0,30,5))
oh.save.pdf(p = p2,outname = "resub_soil_shootcircum_week11.pdf",outdir = "../figures/",
            width = 10,height = 10)


rm(list=ls())
dev.off()
gc()
