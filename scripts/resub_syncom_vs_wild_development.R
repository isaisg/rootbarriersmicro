library(ohchibi)
library(egg)
library(emmeans)
library(multcomp)
library(ggrepel)

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

df_wild11$Area_plantZ <- df_wild11$Area_plant %>% scale


#### Shoot area in the syncom ###
mDat <- readRDS(file = "../cleandata/dat_syncom_genotypes.RDS") 

morder_genos <- c("Col-0","esb1","myb36","sgn3","sgn3myb36","pCASP1::CDEF","esb1.1pCASP1::CDEF","sgn3.3myb36pELTP::CEDF")

######### Dey weight ###########
df_sub <- mDat$DryWeight
df_sub <- which(df_sub$Genotype %in% morder_genos) %>%
  df_sub[.,] %>% droplevels
df_sub$Genotype <- df_sub$Genotype %>%
  factor(levels = morder_genos)

df_a <- df_sub %>% subset(Type == "NB") %>% droplevels
df_b <- df_sub %>% subset(Type != "NB") %>% droplevels

df_a$Normweight_mg_plantZ <- df_a$Normweight_mg_plant %>% scale
df_b$Normweight_mg_plantZ <- df_b$Normweight_mg_plant %>% scale

df_sub <- rbind(df_a,df_b)

merged <- df_sub[,c("Genotype","Type","Normweight_mg_plantZ")]

merged$Genotype <- merged$Genotype %>% 
  gsub(pattern = "^Col-0$",replacement = "Col_0") %>%
  gsub(pattern = "^myb36$",replacement = "myb36_2") %>%
  gsub(pattern = "^esb1$",replacement = "esb1_1") %>%
  gsub(pattern = "^sgn3$",replacement = "sgn3_3") %>%
  gsub(pattern = "^sgn3myb36$",replacement = "myb36_2_sgn3_3") %>%
  gsub(pattern = "^pCASP1::CDEF$",replacement = "pCASP1_CDEFwt") %>%
  gsub(pattern = "^esb1.1pCASP1::CDEF$",replacement = "pCASP1_CDEFesb1_1") %>%
  gsub(pattern = "^sgn3.3myb36pELTP::CEDF$",replacement = "pELTP_CDEFsgn3_3_myb36_1") %>%
  factor() %>% relevel(ref = "Col_0")

merged$Number <- match(merged$Genotype,df_num$Genotype) %>%
  df_num$Number[.]
merged$Group <- match(merged$Genotype,df_num$Genotype) %>%
  df_num$Group[.]
merged$Number <- merged$Number %>% factor
merged$Group <- merged$Group %>% factor


df_a <- df_wild11[,c("NumGenotype","Area_plantZ")]
df_a$Type <- "WildSoil"

df_b <- merged[,c("Number","Normweight_mg_plantZ","Type")]

colnames(df_a) <- c("Number","value","Type")
colnames(df_b) <- c("Number","value","Type")

merged <- rbind(df_a,df_b)
merged <- merged %>% droplevels

## Create graph
merged$Significance <- NA
merged$Type <- merged$Type %>% gsub(pattern = "NB",replacement = "NB (agar system)") %>%
  gsub(pattern = "SynCom",replacement = "SynCom (agar system)") %>%
  factor

chosen_nums <- merged %>% subset(Type == "SynCom (agar system)") %$% Number %>% as.character %>%  unique

merged <- which(merged$Number %in% chosen_nums) %>%
  merged[.,] %>% droplevels

mtypes <- merged$Type %>% as.character %>% unique
Res <- NULL
for (ty in mtypes){
  a <- merged %>% subset(Number == 1 & Type == ty) %$% value %>% mean
  mnums <- chosen_nums %>% grep(pattern = "^1$",value = T,invert = T)
  for (num in mnums){
      b <- merged %>% subset(Number == num & Type == ty) %$% value %>% mean
      diff <- b-a
      Res <- data.frame(Number = num, Estimate = diff,Type = ty) %>%
        rbind(Res,.)
  }
}


mpaleta <- mDat$paleta_type[c(1,3)]
mpaleta <- c(mpaleta,"#950095")
names(mpaleta) <- c("NB (agar system)","SynCom (agar system)","WildSoil")

Res$Number <- Res$Number %>% factor(levels = Res$Number %>% as.numeric %>% unique %>% sort)

df_sig_a <- data.frame(Type = "WildSoil",Number = c(2,11,16,19) , Significance = "q < 0.1")

df_sig_b <- data.frame(Type = c("NB (agar system)","NB (agar system)","SynCom (agar system)"),
                      Number  = c(16,19,19),Significance = "q < 0.1")

df_sig <- rbind(df_sig_a,df_sig_b)
Res <- merge(Res,df_sig, by = c("Number","Type"),all.x = TRUE)

Res$Type <- Res$Type %>% factor(levels = c("WildSoil","NB (agar system)","SynCom (agar system)"))
p1 <- ggplot(data = Res,aes(Type,Estimate)) + 
  geom_linerange(aes(x=Type, ymax=Estimate, ymin=0),size = 0.6) +
  geom_hline(yintercept = 0,color = "#D9D9D9",size = 1) +
  geom_point(size = 5,aes(fill = Type),shape = 21) +
  facet_grid(.~Number) +
  theme_ohchibi(legend_proportion_size = 0.5) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background.x = element_blank(),
    strip.text.x = element_text(family = "Helvetica",face = "plain",size = 20),
    axis.text.y = element_text(family = "Helvetica",size = 15,face = "plain"),
    axis.title.y = element_text(family = "Helvetica",size = 20,face = "plain"),
    legend.text = element_text(family = "Helvetica",size = 15,face = "plain"),
    legend.title = element_text(family = "Helvetica",size = 20,face = "plain"),
    panel.spacing = unit(0.1, "lines")
  ) +
  ylab(label = "Standarised estimate in relation to Col-0") +
  scale_fill_manual(values = mpaleta) 


p2 <- ggplot(data = Res,aes(Type,Estimate)) + 
  geom_hline(yintercept = 0,size = 0.7,color = "#D9D9D9") +
  stat_summary(fun.data  = mean_cl_normal,geom = "pointrange",color = "red",size = 1) +
  geom_text_repel(aes(label = Number,color = Significance),size = 5, position = position_jitter(width = 0.2)) +
  theme_ohchibi() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    #axis.ticks.x = element_blank(),
    strip.background.x = element_blank(),
    strip.text.x = element_text(family = "Helvetica",face = "plain",size = 20),
    axis.text.y = element_text(family = "Helvetica",size = 15,face = "plain"),
    axis.title.y = element_text(family = "Helvetica",size = 20,face = "plain"),
    legend.text = element_text(family = "Helvetica",size = 15,face = "plain"),
    legend.title = element_text(family = "Helvetica",size = 20,face = "plain"),
    panel.spacing = unit(0.1, "lines")
  ) +
  ylab(label = "Standarised estimate in relation to Col-0")

composition <- egg::ggarrange(p = p1 + theme(legend.position = "none"),
               p2 + theme(legend.position = "noen"),nrow = 1,widths = c(1,0.4))
oh.save.pdf(p = composition,outname = "resub_syncom_vs_wild_development.pdf",outdir = "../figures/",
            4,width = 10)


### Absolute version of the estimate 
Res$AbsEstimate <- abs(Res$Estimate)

ggplot(data = Res,aes(Type,AbsEstimate)) + 
  geom_hline(yintercept = 0,size = 0.7,color = "#D9D9D9") +
  stat_summary(fun.data  = mean_cl_normal,geom = "pointrange",color = "red",size = 1) +
  geom_text_repel(aes(label = Number,color = Significance),size = 5, position = position_jitter(width = 0.2)) +
  theme_ohchibi() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    #axis.ticks.x = element_blank(),
    strip.background.x = element_blank(),
    strip.text.x = element_text(family = "Helvetica",face = "plain",size = 20),
    axis.text.y = element_text(family = "Helvetica",size = 15,face = "plain"),
    axis.title.y = element_text(family = "Helvetica",size = 20,face = "plain"),
    legend.text = element_text(family = "Helvetica",size = 15,face = "plain"),
    legend.title = element_text(family = "Helvetica",size = 20,face = "plain"),
    panel.spacing = unit(0.1, "lines")
  ) +
  ylab(label = "Standarised estimate in relation to Col-0")

aov(formula = AbsEstimate ~ Type,data = Res) %>%
  emmeans(object = .,specs = "Type") %>% CLD

rm(list=ls())
dev.off()
gc()
