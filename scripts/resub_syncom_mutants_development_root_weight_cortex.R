library(ohchibi)
library(egg)
library(emmeans)
library(multcomp)

set.seed(130816)
setwd('/Users/isaisalasgonzalez/Documents/rootbarriersmicro/scripts/')


### Syncom data
mDat <- readRDS(file = "../cleandata/dat_syncom_genotypes.RDS") 
paleta_type <- mDat$paleta_type[c(1,3)]

morder_genos <- c("Col-0","esb1","myb36","sgn3","sgn3myb36","pCASP1::CDEF","esb1.1pCASP1::CDEF","sgn3.3myb36pELTP::CEDF")


######### Dey weight ###########
df_sub <- mDat$DryWeight
df_sub <- which(df_sub$Genotype %in% morder_genos) %>%
  df_sub[.,] %>% droplevels
df_sub$Genotype <- df_sub$Genotype %>%
  factor(levels = morder_genos)

#Test
m1_sum <- df_sub %>% subset(Type == "NB") %>%
  aov(formula = Normweight_mg_plant~Genotype) %>%
  glht(., linfct = mcp( Genotype = "Dunnett")) %>%
  summary
res_dw_nb <- data.frame(genotype = names(m1_sum$test$coefficients) %>% gsub(pattern = " - .*",replacement = ""),
                     coefficients = m1_sum$test$coefficients %>% as.numeric,
                     pvalue = m1_sum$test$pvalues %>% as.numeric)

m1_sum <- df_sub %>% subset(Type == "SynCom") %>%
  aov(formula = Normweight_mg_plant~Genotype) %>%
  glht(., linfct = mcp( Genotype = "Dunnett")) %>%
  summary
res_dw_sc <- data.frame(genotype = names(m1_sum$test$coefficients) %>% gsub(pattern = " - .*",replacement = ""),
                        coefficients = m1_sum$test$coefficients %>% as.numeric,
                        pvalue = m1_sum$test$pvalues %>% as.numeric)


#Create structue for plotting
df_min <- df_sub[,c("Genotype","Type","Normweight_mg_plant")]
res_dw_nb$Type <- "NB"
res_dw_sc$Type <- "SynCom"

res <- rbind(res_dw_nb,res_dw_sc)
colnames(res)[1] <- c("Genotype")
res$padj <- p.adjust(p = res$pvalue,method = "fdr")

df_min <- merge(df_min,res, by = c("Genotype","Type"),all.x = TRUE)
df_dw <- df_min





######### Primary root ###########
df_sub <- mDat$Root
df_sub <- which(df_sub$Genotype %in% morder_genos) %>%
  df_sub[.,] %>% droplevels
df_sub$Genotype <- df_sub$Genotype %>%
  factor(levels = morder_genos)

m1_sum <- df_sub %>% subset(Type == "NB") %>%
  aov(formula = NormRoot~Genotype) %>%
  glht(., linfct = mcp( Genotype = "Dunnett")) %>%
  summary
res_root_nb <- data.frame(genotype = names(m1_sum$test$coefficients) %>% gsub(pattern = " - .*",replacement = ""),
                       coefficients = m1_sum$test$coefficients %>% as.numeric,
                       pvalue = m1_sum$test$pvalues %>% as.numeric)



m1_sum <- df_sub %>% subset(Type == "SynCom") %>%
  aov(formula = NormRoot~Genotype) %>%
  glht(., linfct = mcp( Genotype = "Dunnett")) %>%
  summary
res_root_sc<- data.frame(genotype = names(m1_sum$test$coefficients) %>% gsub(pattern = " - .*",replacement = ""),
                       coefficients = m1_sum$test$coefficients %>% as.numeric,
                       pvalue = m1_sum$test$pvalues %>% as.numeric)


#Create structue for plotting
df_min <- df_sub[,c("Genotype","Type","NormRoot")]
res_root_nb$Type <- "NB"
res_root_sc$Type <- "SynCom"

res <- rbind(res_root_nb,res_root_sc)
colnames(res)[1] <- c("Genotype")
res$padj <- p.adjust(p = res$pvalue,method = "fdr")


df_min <- merge(df_min,res, by = c("Genotype","Type"),all.x = TRUE)
df_root <- df_min


###### Cortex #########
df <- read.table(file = "../rawdata/data_volume_cortex.csv",sep = ",",header = T)
df$Genotype <- df$Genotype %>% gsub(pattern = "WT",replacement = "Col-0") %>%
  factor
df$Type <- df$Type %>% gsub(pattern = "SC",replacement = "SynCom") %>%
  factor

m1_sum <- df %>% subset(Type == "NB") %>%
  aov(formula = Volume~Genotype) %>%
  glht(., linfct = mcp( Genotype = "Dunnett")) %>%
  summary
res_cortex_nb <- data.frame(genotype = names(m1_sum$test$coefficients) %>% gsub(pattern = " - .*",replacement = ""),
                          coefficients = m1_sum$test$coefficients %>% as.numeric,
                          pvalue = m1_sum$test$pvalues %>% as.numeric)



m1_sum <- df %>% subset(Type == "SynCom") %>%
  aov(formula = Volume~Genotype) %>%
  glht(., linfct = mcp( Genotype = "Dunnett")) %>%
  summary
res_cortex_sc <- data.frame(genotype = names(m1_sum$test$coefficients) %>% gsub(pattern = " - .*",replacement = ""),
                            coefficients = m1_sum$test$coefficients %>% as.numeric,
                            pvalue = m1_sum$test$pvalues %>% as.numeric)

#Create structue for plotting
df_min <- df[,c("Genotype","Type","Volume")]
res_cortex_nb$Type <- "NB"
res_cortex_sc$Type <- "SynCom"

res <- rbind(res_cortex_nb,res_cortex_sc)
colnames(res)[1] <- c("Genotype")
res$padj <- p.adjust(p = res$pvalue,method = "fdr")


df_min <- merge(df_min,res, by = c("Genotype","Type"),all.x = TRUE)
df_volume <- df_min


### Perepare merged structure for plotting ###
colnames(df_root)[3] <- "Variable"
colnames(df_dw)[3] <- "Variable"
colnames(df_volume)[3] <- "Variable"

df_root$Var <- "Root"
df_dw$Var <- "DryWeight"
df_volume$Var <- "Volume"


df_volume$Genotype <- df_volume$Genotype %>%
  gsub(pattern = "^Col-0 CDEF$",replacement = "pCASP1::CDEF") %>%
  gsub(pattern = "^esb1CDEF$",replacement = "esb1.1pCASP1::CDEF") %>%
  gsub(pattern = "^sgn3myb36 CDEF$",replacement = "sgn3.3myb36pELTP::CEDF")

merged <- rbind(df_root,df_dw,df_volume)

merged$Significance <- NA
merged$Significance[which(merged$padj < 0.1)] <- "q < 0.1"

## Load the number scheme
df_num <- read.table(file = "../rawdata/df_genotypes_to_classification.csv",header = T,sep = ",")
df_num$Number <- 1:20
df_num$Genotype


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

merged$Significance <- merged$Significance %>% factor


p1 <- merged %>% subset(Var == "Root") %>% droplevels %>%
  ggplot(data = .,aes(Number,Variable)) + 
  geom_sina(alpha = 0.1,shape = 21) +
  stat_summary(fun.data = mean_cl_normal,aes(color = Significance)) +
  facet_grid(.~Type,space = "free",scales = "free") +
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
  ylab(label = "Primary root elongation (cm)")  +
  scale_y_continuous(breaks = seq(0,7,1),limits = c(0,7))


p2 <- merged %>% subset(Var == "DryWeight") %>% droplevels %>%
  ggplot(data = .,aes(Number,Variable)) + 
  geom_sina(alpha = 0.1,shape = 21) +
  stat_summary(fun.data = mean_cl_normal,aes(color = Significance)) +
  facet_grid(.~Type,space = "free",scales = "free") +
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
  ylab(label = "Dry weight")  +
  scale_y_continuous(breaks = seq(0,1,0.1),limits = c(0,1))



p3 <- merged %>% subset(Var == "Volume") %>% droplevels %>%
  ggplot(data = .,aes(Number,Variable)) + 
  geom_sina(alpha = 0.1,shape = 21) +
  stat_summary(fun.data = mean_cl_normal,aes(color = Significance)) +
  facet_grid(.~Type,space = "free",scales = "free") +
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
  ylab(label = "Cell volumne")  +
  scale_y_continuous(breaks = seq(24000,58000,5000),limits = c(24000,58000))


### Save the graphs as a composition 
composition <- egg::ggarrange(p1,p2,p3,nrow = 1)
oh.save.pdf(p = composition,outname = "resub_syncom_mutants_development_root_weight_cortex.pdf",
            outdir = "../figures/",width = 20,height = 5)
rm(list=ls())
dev.off()
gc()
