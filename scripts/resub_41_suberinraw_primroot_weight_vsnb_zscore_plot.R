library(ohchibi)
library(emmeans)

set.seed(130816)
setwd('/Users/isaisalasgonzalez/Documents/rootbarriersmicro/scripts/')

Res <- read.table(file = "../cleandata/res_41_suberinraw_primroot_weight_vsnb.zscore.updated.tsv",header = T)

Res <- Res$Strains %>% grep(pattern = "HK",invert = T) %>%
  Res[.,] %>% droplevels
Res$absest <- abs(Res$est)
Res$Variable <- Res$Variable %>% factor(levels = c("Weight","Root","Suberin"))

Res %>% subset(Significance != "NoSignificant") %$% Variable %>% table


p1 <- chibi.boxplot(Map = Res,x_val = "Variable",y_val = "absest",style = "open")  +
  theme(axis.title.x = element_blank()) + ylab(label = "Effect size vs NB")

Res  %>%
  aov(formula =absest~Variable+Strains ) %>%
  emmeans(specs = "Variable") %>% CLD


torem <- Res %>% subset(absest >3) %$% Strains %>% as.character

p2 <- which(!(Res$Strains %in% torem)) %>% Res[.,] %>%
  chibi.boxplot(Map = .,x_val = "Variable",y_val = "absest",style = "open")  +
  theme(axis.title.x = element_blank()) + ylab(label = "Effect size vs NB")

which(!(Res$Strains %in% torem)) %>% Res[.,] %>%
  aov(formula =absest~Variable+Strains ) %>%
  emmeans(specs = "Variable") %>% CLD

composition <- egg::ggarrange(p1,p2,nrow = 1)
oh.save.pdf(p = composition,outname = "resub_41_signal_suberin_root_weight.zscore.pdf",outdir = "../figures/",width = 16,height = 12)
