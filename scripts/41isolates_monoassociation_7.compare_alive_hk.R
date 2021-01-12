library(ohchibi)

set.seed(130816)
setwd('/home/isai/Documents/results/rootbarriersmicro/scripts/')

Res <- read.table(file = "../cleandata/res_41_suberinraw_primroot_weight_vsnb.tsv",header = T)
Res <- Res %>% subset(Variable != "Suberin") %>% droplevels

Res$Type <- "Alive"
Res$Type[Res$Strains %>% grep(pattern = "HK")] <- "HK"
Res$Type <- Res$Type %>% factor

## Root 
df <- Res %>% subset(Variable == "Root")
m <- car::leveneTest(est ~ Type,data = df) 
m$`Pr(>F)`[1]
x <- df %>% subset(Type == "Alive" & Strains != "CDEF") %$% est
y <- df %>% subset(Type == "HK") %$% est
m1 <- t.test(x = x,y = y,exact = F,var.equal = F)
pval_root <- m1$p.value

## Dry weight
df <- Res %>% subset(Variable == "Weight")
m <- car::leveneTest(est ~ Type,data = df) 
m$`Pr(>F)`[1]
x <- df %>% subset(Type == "Alive" & Strains != "CDEF") %$% est
y <- df %>% subset(Type == "HK") %$% est
m1 <- t.test(x = x,y = y,exact = F,var.equal = F)
pval_weight <- m1$p.value


Res$StrainsUnique <- Res$Strains %>% gsub(pattern = "HK",replacement = "")

#create pval data.frame
df_pval <- data.frame(Root = pval_root,Weight = pval_weight) %>% melt
colnames(df_pval) <- c("Variable","Label")
df_pval$Label <- format.pval(df_pval$Label)
## Do a graph 
p <- Res %>% 
  #subset(est > - 1) %>% 
  ggplot(data = .,aes(Type,est)) +
  geom_hline(yintercept = 0,linetype = "longdash",size = 2,color = "#D9D9D9")+
  geom_point() +
  geom_line(aes(group = StrainsUnique)) +
  stat_summary(fun.y = mean,geom = "point",size =8,color = "red") +
  facet_grid(.~Variable) +
  theme_ohchibi(size_panel_border = 2) +
  theme(strip.background.x = element_blank(),
        strip.text.x = element_text(family = "Arial",face = "bold",size = 20),
        axis.title.x = element_blank()
        )  +
  ylab(label = "Estimate") +
  geom_text(data = df_pval,aes(1.5,0.55,label =Label )) 

oh.save.pdf(p = p,outname = "41_root_weight_alive_hk.pdf",outdir = "../figures/",width = 14,height = 10)

### absolute estimate version 
rm(list=ls())


Res <- read.table(file = "../cleandata/res_41_suberinraw_primroot_weight_vsnb.tsv",header = T)
Res <- Res %>% subset(Variable != "Suberin") %>% droplevels

Res$Type <- "Alive"
Res$Type[Res$Strains %>% grep(pattern = "HK")] <- "HK"
Res$Type <- Res$Type %>% factor
Res$est <- abs(Res$est)

## Root 
df <- Res %>% subset(Variable == "Root")
m <- car::leveneTest(est ~ Type,data = df) 
m$`Pr(>F)`[1]
x <- df %>% subset(Type == "Alive" & Strains != "CDEF") %$% est
y <- df %>% subset(Type == "HK") %$% est
m1 <- t.test(x = x,y = y,exact = F,var.equal = F)
pval_root <- m1$p.value

## Dry weight
df <- Res %>% subset(Variable == "Weight")
m <- car::leveneTest(est ~ Type,data = df) 
m$`Pr(>F)`[1]
x <- df %>% subset(Type == "Alive" & Strains != "CDEF") %$% est
y <- df %>% subset(Type == "HK") %$% est
m1 <- t.test(x = x,y = y,exact = F,var.equal = F)
pval_weight <- m1$p.value


Res$StrainsUnique <- Res$Strains %>% gsub(pattern = "HK",replacement = "")

#create pval data.frame
df_pval <- data.frame(Root = pval_root,Weight = pval_weight) %>% melt
colnames(df_pval) <- c("Variable","Label")
#df_pval$Label <- format.pval(df_pval$Label)
## Do a graph 
p <- Res %>% 
  subset(Variable == "Root") %>% droplevels %>%
  ggplot(data = .,aes(Type,est)) +
  geom_hline(yintercept = 0,linetype = "longdash",size = 2,color = "#D9D9D9")+
  geom_point() +
  geom_line(aes(group = StrainsUnique)) +
  stat_summary(fun.y = mean,geom = "point",size =8,color = "red") +
  facet_grid(.~Variable) +
  theme_ohchibi(size_panel_border = 2) +
  theme(strip.background.x = element_blank(),
        strip.text.x = element_text(family = "Arial",face = "bold",size = 20),
        axis.title.x = element_blank()
  )  +
  ylab(label = "Estimate") +
  geom_text(data = df_pval %>% subset(Variable == "Root"),aes(2,3,label =Label ))

## Do a graph 
p2 <- Res %>% 
  subset(Variable == "Weight") %>% droplevels %>%
  ggplot(data = .,aes(Type,est)) +
  geom_hline(yintercept = 0,linetype = "longdash",size = 2,color = "#D9D9D9")+
  geom_point() +
  geom_line(aes(group = StrainsUnique)) +
  stat_summary(fun.y = mean,geom = "point",size =8,color = "red") +
  facet_grid(.~Variable) +
  theme_ohchibi(size_panel_border = 2) +
  theme(strip.background.x = element_blank(),
        strip.text.x = element_text(family = "Arial",face = "bold",size = 20),
        axis.title.x = element_blank()
  )  +
  ylab(label = "Estimate") +
  geom_text(data = df_pval %>% subset(Variable == "Weight"),aes(1.5,0.55,label =Label )) +
  scale_y_continuous(breaks = c(0,0.1,0.2,0.3,0.4,0.5))

p <- egg::ggarrange(p,p2,nrow = 1)

oh.save.pdf(p = p,outname = "41_root_weight_alive_hk_abs.pdf",outdir = "../figures/",width = 15,height = 10)
