library(ohchibi)
library(emmeans)
library(egg)
library(car)
library(ggpubr)
library(paletteer)

#Merge datasets into a single structure
set.seed(130816)
setwd('/home/isai/Documents/results/rootbarriersmicro/scripts/')


Dat <- readRDS(file = "../cleandata/dat_syncom_stresses.RDS")
df_sub <- Dat$Suberin

df_sub$Plate <- rep(c(1:4),40)
#df_sub <- df_sub[,c(3,4,5,8,9)]
df_sub <- df_sub[,c(1,2,6,8,11)]

df_sub$Stress <- df_sub$Stress %>% gsub(pattern = " ",replacement = "")
df_sub$UId <- paste0(df_sub$Bacteria,"|",df_sub$Stress,"|",
                     df_sub$Rep,"|",df_sub$Plate)


df_sub <- df_sub %>% na.omit 

df_sub$Stress <- df_sub$Stress %>% factor %>% relevel(ref = "Full")


#IOonome
Dat_ionome <- Dat$Ionome
#Dat_ionome <- Dat_ionome %>% subset.Dataset(Type != "HK",drop = T,clean = T)
Tab_z <- Dat_ionome$Tab %>% t %>% scale %>% t 
melted <- Tab_z %>% melt
colnames(melted) <- c("Ion","Index","zscore")

melted <- merge(melted,Dat_ionome$Map, by = "Index") 
melted$Stress <- melted$Stress %>% factor %>% relevel(ref = "Full")
melted$Type <- melted$Type %>% factor(levels = c("NB","HK","SynCom"))

melted$Rep <- melted$Replicates %>% gsub(pattern = "Rep",replacement = "")
melted$Plate <- melted$Plates %>% gsub(pattern = "Plate",replacement = "")
melted$Bacteria <- melted$Type

melted$UId <- paste0(melted$Bacteria,"|",melted$Stress,"|",
                     melted$Rep,"|",melted$Plate)

dim(melted)
dim(df_sub)

merged <- merge(melted[,c(2,3,15)],df_sub, by = "UId")


a <-melted$UId %>% unique
b <- df_sub$UId %>% unique

length(a)
length(b)
length(intersect(a,b))

paleta_stress <- Dat$paleta_stress

p <- ggplot(data = merged,aes(NormSuberin,zscore)) +
  geom_smooth(method = "lm",se = FALSE,color = "black",size = 2) +
  geom_point(shape = 16,size = 4,aes(color = Stress),alpha = 1) +
  scale_color_manual(values = paleta_stress) +
  facet_grid(Ion~Bacteria,scales = "free_y") +
  stat_cor(size = 7) + 
  theme_ohchibi(size_panel_border = 2,size_axis_text.y = 10) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_text(family = "Arial",face = "bold",size = 20),
    strip.text.y = element_text(family = "Arial",face = "bold",size = 20)
  ) +
  xlab(label = "Suberin") + ylab(label = "z-score")
oh.save.pdf(p = p,outname = "correlation_ion_nb_syncom_stresses.pdf",outdir = "../figures/",width =20,height = 30 )


#Calculate the p.values
Res_ion <- NULL
for(ion in merged$Ion %>% unique){
  merged_sub <- merged %>% subset(Ion == ion) %>% droplevels
  for(type in merged_sub$Bacteria %>% unique){
    temp <- merged_sub %>% subset(Bacteria == type) %>% droplevels
    mres <- cor.test(temp$NormSuberin,temp$zscore)
    Res_ion <- data.frame(r = mres$estimate,p.value = mres$p.value,Bacteria = type,Ion = ion) %>%
      rbind(Res_ion,.)
  }
}

Res_ion$p.adj <- Res_ion$p.value %>% p.adjust(method = "fdr")
Res_ion$Significance <- "NotSignificant"
Res_ion$Significance[which(Res_ion$p.adj < 0.05)] <- "Significant"


a <- Res_ion %>% subset(Bacteria ==  "NB") %$% r
mean(a)
mean(a %>% abs)

b <- Res_ion %>% subset(Bacteria ==  "SynCom") %$% r
mean(b)
mean(b %>% abs)

m1 <- wilcox.test((a),(b))
mpval <- m1$p.value %>% format.pval
t.test((a),(b))

p <- ggplot(Res_ion,aes(Bacteria,r)) +
  geom_hline(yintercept = 0,size = 2,color = "#D9D9D9",linetype = "longdash") +
  stat_summary(fun.y = "mean",geom = "errorbar",aes(ymax = ..y.., ymin = ..y..),
               width = 1, linetype = "longdash",color = "red")+
  geom_line(aes(group = Ion),size = 2) +
  geom_point(shape = 21,aes(fill = Ion),size = 12) +
  theme_ohchibi(size_panel_border = 2) +
  xlab(label = NULL) + ylab(label = "Pearson Correlation Coefficient") +
  scale_fill_paletteer_d(package = "ggthemes",palette = "Tableau 20")+
  scale_color_paletteer_d(package = "ggthemes",palette = "Tableau 20") +
  scale_y_continuous(breaks = seq(-1,1,0.2),limits = c(-1,1))

oh.save.pdf(p = p,outname = "cumulative_correlation_raw.pdf",outdir = "../figures/",width = 10,height = 12)

m1 <- wilcox.test((a %>% abs),(b %>% abs))
mpval <- m1$p.value %>% format.pval
t.test((a %>% abs),(b %>% abs))


p <- ggplot(Res_ion,aes(Bacteria,abs(r))) +
  stat_summary(fun.y = "mean",geom = "errorbar",aes(ymax = ..y.., ymin = ..y..),
               width = 1, linetype = "longdash",color = "red")+
  geom_line(aes(group = Ion),size = 2) +
  geom_point(shape = 21,aes(fill = Ion),size = 12) +
  geom_text(aes(1.5,0.8,label = paste0("pval = ",mpval)),size = 10)+
  theme_ohchibi(size_panel_border = 2) +
  xlab(label = NULL) + ylab(label = "Absolute Pearson Correlation Coefficient") +
  scale_fill_paletteer_d(package = "ggthemes",palette = "Tableau 20")+
  scale_color_paletteer_d(package = "ggthemes",palette = "Tableau 20") +
  scale_y_continuous(breaks = seq(0,1,0.2),limits = c(0,1))
oh.save.pdf(p = p,outname = "cumulative_correlation_abs.pdf",outdir = "../figures/",width = 10,height = 12)


Res_ion$absr <- abs(Res_ion$r)
leveneTest(absr ~Bacteria,data = Res_ion)
bartlett.test(absr ~Bacteria,data = Res_ion)
