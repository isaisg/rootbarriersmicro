library(ohchibi)
library(emmeans)
library(egg)

set.seed(130816)
setwd('/home/isai/Documents/results/rootbarriersmicro/scripts/')

Dat <- readRDS(file = "../cleandata/dat_syncom_stresses.RDS") %$%
  Ionome
paleta_stress <- readRDS(file = "../cleandata/dat_syncom_stresses.RDS") %$%
  paleta_stress
paleta_type <- readRDS(file = "../cleandata/dat_syncom_stresses.RDS") %$%
  paleta_type

#Perform testing inside each ion
#We wanna do pairwise comparison between  StressNB vs Full NB
#StressSC vs Full SC and #Intra Stress SC vs Stress NB
Tab_z <- Dat$Tab %>% t %>% scale %>% t 

mpca <- prcomp(x = Tab_z %>% t,center = F,scale. = F)
scores <- mpca$x %>% 
  as.data.frame
scores$Index <- rownames(scores)
scores <- merge(scores,Dat$Map, by = "Index")



 p1 <- ggplot(scores,aes(PC1,PC2)) +
   geom_vline(xintercept = 0,linetype = "dashed",size = 4 , color = "#D9D9D9") +
   geom_hline(yintercept = 0,linetype = "dashed",size = 4 , color = "#D9D9D9") +
   geom_point(shape = 21,aes(fill = Stress),size = 12,stroke = 1) +
   theme_ohchibi(size_ticks = NA,size_panel_border = 2)  +
   theme(
     axis.text.x = element_text(hjust = 0.5),
     legend.position = "none"
   ) +
   xlab(label = "PC1 (30.01% Variance Explained)") +
   ylab(label = "PC2 (13.73% Variance Explained)")  +
   scale_fill_manual(values = paleta_stress) + coord_fixed()

 p2 <- ggplot(scores,aes(PC1,PC2)) +
   geom_vline(xintercept = 0,linetype = "dashed",size = 4 , color = "#D9D9D9") +
   geom_hline(yintercept = 0,linetype = "dashed",size = 4 , color = "#D9D9D9") +
   geom_point(shape = 21,aes(fill = Type),size = 12,stroke = 1) +
   theme_ohchibi(size_ticks = NA,size_panel_border = 2)  +
   theme(
     axis.text.x = element_text(hjust = 0.5),
     legend.position = "none"
   ) +
   xlab(label = "PC1 (30.01% Variance Explained)") +
   ylab(label = "PC2 (13.73% Variance Explained)")  +
   scale_fill_manual(values = paleta_type) + coord_fixed()

 composition <- ggarrange(p1 +theme(legend.position = "none"),p2 + theme(legend.position = "none"),nrow = 1)
 oh.save.pdf(p = composition,outname = "composition_ionome_syncom_stresses_all.pc.pdf",
             outdir = "../figures/",width = 25,height = 15)

 
 
 p1 <- ggplot(scores,aes(PC1,PC2)) +
   geom_vline(xintercept = 0,linetype = "dashed",size = 4 , color = "#D9D9D9") +
   geom_hline(yintercept = 0,linetype = "dashed",size = 4 , color = "#D9D9D9") +
   geom_point(aes(fill = Stress,shape = Type),size = 12,stroke = 1) +
   theme_ohchibi(size_ticks = NA,size_panel_border = 2)  +
   theme(
     axis.text.x = element_text(hjust = 0.5),
     legend.position = "none"
   ) +
   xlab(label = "PC1 (30.01% Variance Explained)") +
   ylab(label = "PC2 (13.73% Variance Explained)")  +
   scale_fill_manual(values = paleta_stress) + coord_fixed() +
   scale_shape_manual(values = c(21,22,24) %>% rev)
 
 oh.save.pdf(p = p1,outname = "composition_ionome_syncom_stresses_all.pc.shape.pdf",
             outdir = "../figures/",width = 25,height = 15)
 
 

### Perform cap constraining for the rep effect
size_legend_text <- 45
size_axis_title <- 20
size_axis_text <- 15
size_legend_title <- 55
legend_proportion_size <- 4
size_points <- 10
strip_text_size <- 25
size_median <- 4

distfun <- function(x,method) vegan::vegdist(x = x, method = "euclidean")
Dat_z <- create_dataset(Tab = Tab_z,Map = Dat$Map)
mcap <- oh.cap(Tab = Dat_z$Tab,Map = Dat_z$Map,
               formula = "Stress  + Type ",
               distfun = distfun,perms = 9999,sqrt = T)
mcap$perm_anova_terms



p1 <- chibi.cap(list_ohpco = mcap,col_val = "Stress",
          comp_a = "CAP1",comp_b = "CAP2" ,size =size_points ,alpha=1,
          size_legend_text = size_legend_text,size_axis_title = size_axis_title,size_axis_line = 2,
          size_title_text = size_legend_title,size_axis_text =  size_axis_text,
          font_family = "Arial",legend_proportion_size = 4,lines_zero = T,size_panel_border = 2) +
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())  +
  scale_fill_manual(values = paleta_stress) + coord_fixed()
p2 <- chibi.cap(list_ohpco = mcap,col_val = "Type",
                comp_a = "CAP1",comp_b = "CAP2" ,size =size_points ,alpha=1,
                size_legend_text = size_legend_text,size_axis_title = size_axis_title,size_axis_line = 2,
                size_title_text = size_legend_title,size_axis_text =  size_axis_text,
                font_family = "Arial",legend_proportion_size = 4,lines_zero = T,size_panel_border = 2) +
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())  +
  scale_fill_manual(values = paleta_type) + coord_fixed()


###Perform mean projection ####
Mapa <- mcap$Map_cap
df_sum <- aggregate(CAP1~Stress+Type,Mapa,mean) %>%
  merge(.,aggregate(CAP2~Stress+Type,Mapa,mean), by = c("Stress","Type"))

df_temp <- df_sum %>% subset(Type == "NB") 
tab <- df_temp[,3:4] %>% as.matrix
rownames(tab) <- df_temp$Stress
Tab_nb <- tab

df_temp <- df_sum %>% subset(Type == "HK") 
tab <- df_temp[,3:4] %>% as.matrix
rownames(tab) <- df_temp$Stress
Tab_hk <- tab

df_temp <- df_sum %>% subset(Type == "SynCom") 
tab <- df_temp[,3:4] %>% as.matrix
rownames(tab) <- df_temp$Stress
Tab_syncom <- tab

mantel(Tab_nb %>% dist,Tab_syncom %>% dist,permutations = 9999)
mantel(Tab_nb %>% dist,Tab_hk %>% dist,permutations = 9999)

mantel(Tab_hk %>% dist,Tab_syncom %>% dist,permutations = 9999)

melt_nb <- Tab_nb %>% dist %>% as.matrix %>% melt_dist
melt_hk <- Tab_hk %>% dist %>% as.matrix %>% melt_dist
melt_sc <- Tab_syncom %>% dist %>% as.matrix %>% melt_dist


melt_nb$Type <- "NB"
melt_hk$Type <- "HK"
melt_sc$Type <- "SynCom"

melted <- rbind(melt_nb,melt_hk,melt_sc)

melted$iso1 %>% grep(pattern = "Full")
melted$iso2 %>% grep(pattern = "Full")

melted$iso1 %>% grep(pattern = "Full") %>% melted[.,] %>%
  chibi.boxplot(Map = .,x_val = "Type",y_val = "dist",style = "open")

melted$Type <- melted$Type %>% factor(levels = c("NB","HK","SynCom"))


melted %>% subset(iso1 == "Full") %>% droplevels %>%
  subset(Type != "HK") %>%
  ggplot(data = .,aes(Type,dist)) +
  stat_summary(fun.y = "mean",geom = "errorbar",aes(ymax = ..y.., ymin = ..y..),
               width = 1, linetype = "longdash",color = "red")+
  geom_point(aes(fill = iso2),shape = 21,size = 12) +
  geom_line(aes(group = iso2)) + 
  theme_ohchibi() +
  scale_fill_manual(values = paleta_stress) +
  scale_y_continuous(breaks = seq(0,1.1,by = 0.1)) +
  theme(
    # panel.grid.major.y = element_line(size = 0.5,color = "#D9D9D9"),
    #  panel.grid.minor.y = element_line(size = 0.5,color = "#D9D9D9"),
    panel.background = element_rect(fill = "white"),
    strip.text.x = element_text(family = "Arial",face = "bold",size = 20),
    axis.title.x = element_blank()
  )    +
  ylab(label = "Distance to Full")


melted %>% subset(iso1 == "Full") %>% droplevels %>%
  subset(Type != "SynCom") %>%
  ggplot(data = .,aes(Type,dist)) +
  stat_summary(fun.y = "mean",geom = "errorbar",aes(ymax = ..y.., ymin = ..y..),
               width = 1, linetype = "longdash",color = "red")+
  geom_point(aes(fill = iso2),shape = 21,size = 12) +
  geom_line(aes(group = iso2)) + 
  theme_ohchibi() +
  scale_fill_manual(values = paleta_stress) +
  scale_y_continuous(breaks = seq(0,1.1,by = 0.1)) +
  theme(
    # panel.grid.major.y = element_line(size = 0.5,color = "#D9D9D9"),
    #  panel.grid.minor.y = element_line(size = 0.5,color = "#D9D9D9"),
    panel.background = element_rect(fill = "white"),
    strip.text.x = element_text(family = "Arial",face = "bold",size = 20),
    axis.title.x = element_blank()
  )    +
  ylab(label = "Distance to Full")

melted %>% subset(iso1 == "Full") %>% droplevels %>%
  subset(Type != "NB") %>%
  ggplot(data = .,aes(Type,dist)) +
  stat_summary(fun.y = "mean",geom = "errorbar",aes(ymax = ..y.., ymin = ..y..),
               width = 1, linetype = "longdash",color = "red")+
  geom_point(aes(fill = iso2),shape = 21,size = 12) +
  geom_line(aes(group = iso2)) + 
  theme_ohchibi() +
  scale_fill_manual(values = paleta_stress) +
  scale_y_continuous(breaks = seq(0,1.1,by = 0.1)) +
  theme(
    # panel.grid.major.y = element_line(size = 0.5,color = "#D9D9D9"),
    #  panel.grid.minor.y = element_line(size = 0.5,color = "#D9D9D9"),
    panel.background = element_rect(fill = "white"),
    strip.text.x = element_text(family = "Arial",face = "bold",size = 20),
    axis.title.x = element_blank()
  )    +
  ylab(label = "Distance to Full")


melted %>% subset(iso1 == "Full") %>% droplevels %>%
  ggplot(data = .,aes(Type,dist)) +
  stat_summary(fun.y = "mean",geom = "errorbar",aes(ymax = ..y.., ymin = ..y..),
               width = 1, linetype = "longdash",color = "red")+
  geom_point(aes(fill = iso2),shape = 21,size = 12) +
  geom_line(aes(group = iso2)) + 
  theme_ohchibi() +
  scale_fill_manual(values = paleta_stress) +
  scale_y_continuous(breaks = seq(0,1.1,by = 0.1)) +
  theme(
    # panel.grid.major.y = element_line(size = 0.5,color = "#D9D9D9"),
    #  panel.grid.minor.y = element_line(size = 0.5,color = "#D9D9D9"),
    panel.background = element_rect(fill = "white"),
    strip.text.x = element_text(family = "Arial",face = "bold",size = 20),
    axis.title.x = element_blank()
  )    +
  ylab(label = "Distance to Full")

### Do the formal test ### 
df_res <- lm(formula = CAP1 ~ Stress+Type + Stress:Type,data = Mapa) %>% 
  emmeans(pairwise~Stress:Type,adjust = "none") 

res <- df_res$contrasts  %>% as.data.frame
res <- res$contrast %>% grep(pattern = "Full,NB.*NB|Full,SynCom.*SynCom|Full,HK.*HK") %>% res[.,]
res$p.adj <- res$p.value %>% p.adjust(method = "bonferroni") 
res$Significance <- "NS"
res$Significance[which(res$p.adj < 0.05)] <- "Significant"

res$Type <- "NB"
res$Type[res$contrast %>% grep(pattern = "SynCom")] <- "SynCom"
res$Type[res$contrast %>% grep(pattern = "HK")] <- "HK"
res$Type <- res$Type %>% factor(levels = c("NB","HK","SynCom"))

res$Stress <- res$contrast %>% gsub(pattern = ".*- ",replacement = "") %>%
  gsub(pattern = ",.*",replacement = "")

res$Axis <- "CAP1"
res_cap1 <- res


df_res <- lm(formula = CAP2 ~ Stress+Type + Stress:Type,data = Mapa) %>% 
  emmeans(pairwise~Stress:Type,adjust = "none") 

res <- df_res$contrasts  %>% as.data.frame
res <- res$contrast %>% grep(pattern = "Full,NB.*NB|Full,SynCom.*SynCom|Full,HK.*HK") %>% res[.,]
res$p.adj <- res$p.value %>% p.adjust(method = "bonferroni") 
res$Significance <- "NS"
res$Significance[which(res$p.adj < 0.05)] <- "Significant"

res$Type <- "NB"
res$Type[res$contrast %>% grep(pattern = "SynCom")] <- "SynCom"
res$Type[res$contrast %>% grep(pattern = "HK")] <- "HK"
res$Type <- res$Type %>% factor(levels = c("NB","HK","SynCom"))

res$Stress <- res$contrast %>% gsub(pattern = ".*- ",replacement = "") %>%
  gsub(pattern = ",.*",replacement = "")

res$Axis <- "CAP2"
res_cap2 <- res

res <- rbind(res_cap1,res_cap2)

res %>% subset(Type != "HK") %>% droplevels %>%
  ggplot(data = .,aes(Type,abs(estimate))) +
  geom_point(aes(fill = Stress,color = Significance),shape = 21, size = 12) +
  facet_grid(.~Axis) +
  stat_summary(fun.y = "mean",geom = "errorbar",aes(ymax = ..y.., ymin = ..y..),
               width = 1, linetype = "longdash",color = "red")+
  theme_ohchibi()  +
  scale_fill_manual(values = paleta_stress) +
  scale_y_continuous(breaks = seq(0,1.1,by = 0.1)) +
  theme(
    # panel.grid.major.y = element_line(size = 0.5,color = "#D9D9D9"),
    #  panel.grid.minor.y = element_line(size = 0.5,color = "#D9D9D9"),
    panel.background = element_rect(fill = "white"),
    strip.text.x = element_text(family = "Arial",face = "bold",size = 20),
    axis.title.x = element_blank()
  )    +
  ylab(label = "Distance to Full")



### Formula constraining
mcap <- oh.cap(Tab = Dat_z$Tab,Map = Dat_z$Map,
               formula = "Type + Condition(Stress)",
               distfun = distfun,perms = 9999,sqrt = T)
mcap$perm_anova_terms

p3 <- chibi.cap(list_ohpco = mcap,col_val = "Type",
                comp_a = "CAP1",comp_b = "CAP2" ,size =size_points ,alpha=1,
                size_legend_text = size_legend_text,size_axis_title = size_axis_title,size_axis_line = 2,
                size_title_text = size_legend_title,size_axis_text =  size_axis_text,
                font_family = "Arial",legend_proportion_size = 4,lines_zero = T,size_panel_border = 2) +
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())  +
  scale_fill_manual(values = paleta_type) 


composition <- egg::ggarrange(p1,p2,nrow = 1)
oh.save.pdf(p = composition,outname = "composition_ionome_cap_syncom_stresses_all.pdf",
            outdir = "../figures/",width = 18,height = 15)


oh.save.pdf(p = p3,outname = "ionome_cap_constrained_syncom_stresses_all.pdf",
            outdir = "../figures/",width = 9,height = 7)


