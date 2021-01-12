library(ohchibi)
library(emmeans)
library(egg)
library(car)
library(ggrepel)
library(paletteer)
library(ggpubr)
#Merge datasets into a single structure
set.seed(130816)
setwd('/home/isai/Documents/results/rootbarriersmicro/scripts/')

Dat <- readRDS("../cleandata/dat_syncom_stresses.RDS")
paleta <- Dat$paleta_type[c(1,3)]

df <- read.table(file = "../rawdata/data_syncomstress_monomers.csv",header = T,sep = "\t")
df$Treatment <- df$Treatment %>% gsub(pattern = "SC",replacement = "SynCom") %>% factor
df <- df %>% subset(Genotype == "WT") %>% droplevels
colnames(df)

### Test each individual compounts 
Res <- NULL
for(class in (df$Monomer.class %>% as.character %>% unique)){
  df_temp <- subset(df,Monomer.class == class) %>% droplevels
  for(mono in df_temp$Monomer  %>% as.character %>% unique){
    cat("Working on ",class," ",mono,"\n")
    df_sub <- df_temp %>% subset(Monomer == mono)
    mres <- leveneTest(Amount~Treatment,df_sub)
    x <- df_sub %>% subset(Treatment == "NB") %$% Amount
    y <- df_sub %>% subset(Treatment == "SynCom") %$% Amount
    
    if(mres$`Pr(>F)`[1]< 0.05){
      mres_t <- t.test(x,y,var.equal = F)
    }else{
      mres_t <- t.test(x,y,var.equal = T)
    }
    Res <- data.frame(Monomer.class = class,Monomer = mono,p.value.levene =mres$`Pr(>F)`[1],
             p.value.ttest = mres_t$p.value) %>%
      rbind(Res,.)
  }
}

Res$p.adj.ttest <- Res$p.value.ttest  %>%p.adjust(method = "fdr")
Res$Significance <- NA
Res$Significance[which(Res$p.adj.ttest < 0.05)] <- "*"



p <- df %>% subset(Monomer != "Total") %>% droplevels %>%
  chibi.boxplot(Map = .,x_val = "Monomer",
                   y_val = "Amount",col_val = "Treatment",
                facet_formula = "Monomer.class",mpalette = paleta,
                alpha_point = 1,median_colored_as_points = T) +
  scale_y_continuous(breaks = seq(0,6,by = 1))

p1 <- df %>% subset(Monomer == "Total") %>% droplevels %>%
  chibi.boxplot(Map = .,x_val = "Monomer",
                y_val = "Amount",col_val = "Treatment",
                facet_formula = "Monomer.class",mpalette = paleta,
                alpha_point = 1,median_colored_as_points = T) +
  scale_y_continuous(breaks = c(15,16,17,18,19,20,21,22,23,24,25,26,27,28))
composition <- egg::ggarrange(p + theme(legend.position = "none"),p1 + theme(legend.position =  "none"),
               nrow = 1,widths = c(1,0.075))
oh.save.pdf(p = composition,outname = "monomer.quantification.sup.pdf",
            outdir = "../figures/",width = 18,height = 9)

