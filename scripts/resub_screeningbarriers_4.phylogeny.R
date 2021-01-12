library(ohchibi)
library(egg)
library(paletteer)
library(ggtree)
library(car)
library(palettesPM)
library(scales)
library(Rmisc)

#Graphiac parameters
color_line <- "red"
#color_shadow <- "#097a60"
color_shadow <- "#eeeeee"
color_border_circle <- "#222831"
#color_fill_circle <- "#07614C"
color_fill_circle <- "#222831"
color_error <- "#00adb5"
size_point <- 3
size_text <- 11
stroke_point <- 0
width_line <- 2.5


setwd('/Users/isaisalasgonzalez/Documents/rootbarriersmicro/scripts/')



merged <- read.table(file = "../cleandata/res_merged_propidium_suberin.tsv",header = T,sep = "\t")
Dat <- readRDS(file = "../cleandata/dat_screeningbarriers.RDS")


merged_cor <- merged %>% subset(Propidium_Mean < 0.2) %>% droplevels
mcor <- cor.test(merged_cor$Propidium_Mean,merged_cor$Suberin_Mean)


dfpval <- data.frame(label = paste0("r = ",round(mcor$estimate,3),
                                    "\n p-value = ",base::format.pval(pv = mcor$p.value,digits = 2))
)


p <- ggplot(data = merged_cor,aes(Propidium_Mean,Suberin_Mean)) +
  geom_smooth(method = "lm",se = F,color = color_line,fill = color_shadow,size = width_line) +
  geom_text(dfpval,mapping = aes(0.1,0.15,label =label),size = size_text, color = "black",fontface = "bold") +
  #geom_errorbar(mapping = aes(ymin = Suberin_Mean - Suberin_ci,
  #                            ymax = Suberin_Mean +Suberin_ci ),apha = 0.1)+
  #geom_errorbarh(mapping = aes(xmin = Propidium_Mean - Propidium_ci,
  #                          xmax = Propidium_Mean +Propidium_ci ,alpha = 0.1))
  geom_point(size = size_point,fill = color_fill_circle,shape = 21, color = color_border_circle,stroke = stroke_point) +
  theme_ohchibi(size_legend_text = 30,size_panel_border = 2)  +
  xlab(label = "Propidium") +
  ylab(label = "No Expression + Discrete Expression Zones")  



#Signal plot
merged_sub <-   merged %>% subset(Strains != "NB") %>% droplevels
mdf_dist <-melt(data = merged,id.vars = "Strains",
                measure.vars = c("Propidium_est","Suberin_est"))
mdf_dist$variable <- mdf_dist$variable %>% gsub(pattern = "_fc",replacement = "")
colnames(mdf_dist)[2] <- "Variable"

ks.test(merged_sub$Propidium_est,merged_sub$Suberin_est)



##### Phylogenetic figure ########
tree <- Dat$tree
Map <- Dat$Map

tree$tip.label <- match(tree$tip.label,Map$taxon_oid) %>%
  Map$Strains[.] 


rownames(Map) <- Map$Strains
Map <- Map[,c(16,1:15)]
#Define palette
df_colors <- pm.colors.phyla() %>% 
  data.frame(Phyla = names(.), Color = .,row.names = NULL)
df_deino <- data.frame(Phyla = "Deinococcus-Thermus", Color = "#7FE85D")
df_colors <- rbind(df_colors,df_deino)

paleta_phyla <- df_colors$Color %>% as.character
names(paleta_phyla) <- df_colors$Phyla


p <- ggtree(tree,ladderize = F,size = 0.65) + 
  #The expand parth is fundamental to make the tree fit the composition with the other plots
  scale_y_reverse(expand = c(0,0))
p  <- ggtree::flip(p,719,596) 
p_tree <- p  %<+% Map  +
  geom_tiplab(align = T,size = 0,linetype = "solid",aes(color=(Phylum))) +
  scale_color_manual(values = paleta_phyla,na.value = "black") +
  theme(legend.position = "none")
p_tree <- p_tree + coord_flip()  + scale_x_reverse()

#Determine the other of the tips
d <- p$data %>% subset(isTip)
mtips <- with(d, label[order(y, decreasing=T)])

merged_tree <- match(mtips,merged$Strains) %>%
  merged[.,] %>% droplevels


merged_tree$Strains <- merged_tree$Strains %>% factor(levels = mtips)
merged_tree$SignificanceProp <- "NotSignificant"
merged_tree$SignificanceProp[which((merged_tree$Propidium_p.adj < 0.1) & (merged_tree$Propidium_est > 0 ))] <- "Up"
merged_tree$SignificanceProp[which((merged_tree$Propidium_p.adj < 0.1) & (merged_tree$Propidium_est < 0 ))] <- "Down"


paleta_prop <- c("#CCCCCC","#E65037","#2480FF")
names(paleta_prop) <- c("NotSignificant","Up","Down")

mean_nb_prop <- merged %>% subset(Strains == "NB") %$% Propidium_Mean
p_prop <- ggplot(data = merged_tree,aes(Strains,Propidium_Mean)) +
  geom_bar(stat = "identity",aes(fill = SignificanceProp),width = 1,color = "#00000000") +
  geom_hline(yintercept = mean_nb_prop,size =1, color = "black") +
  scale_y_continuous(expand = c(0,0),limits = c(0,0.3),oob = rescale_none,breaks = c(0,0.1,0.2,0.3)) +
  theme_ohchibi(y_vjust = 0) +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(family = "Arial",
                                face="bold",size=12,colour="black",
                                angle = 90,vjust = 1),
    axis.text.y =element_text(family = "Arial",
                              face="bold",size=10,colour="black",
                              angle = 0,vjust = 0.5),
    strip.text.x = element_blank(),
    panel.spacing.x = unit(0.2, "lines")
  )   +
  scale_fill_manual(values =paleta_prop) +
  ylab(label = "Propidium")

#Now do the suberin figure
melted_suberin <- Dat$Suberin %>% dcast(data = .,formula =Strains~Localization,
                                        fun.aggregate = mean,value.var = "Normvalue_ra" )  %>%
  melt
nb_sub <- melted_suberin %>% subset(Strains == "NB") %>% droplevels
melted_sub <- which(melted_suberin$Strains %in% mtips) %>%
  melted_suberin[.,]
melted_sub$Strains <- melted_sub$Strains %>% factor(levels = mtips)
melted_sub$variable <- melted_sub$variable %>% 
  factor(levels = c("Continous_expression","Discrete_expression","No_expression"))

#Suberin plot 
paleta_suberin <- c("#C4AA88","#697A55","#8399B3")


names(paleta_suberin) <- c("No_expression","Discrete_expression","Continous_expression")

p_sub <- ggplot(data = melted_sub,aes(Strains,value, fill = variable)) +
  geom_bar(stat = "identity",width = 2) + 
  geom_hline(yintercept = nb_sub %>% subset(variable == "No_expression") %$% value,size =1, 
             color = "black",linetype = "solid") +
  geom_hline(yintercept = 1-(nb_sub %>% subset(variable == "Continous_expression") %$% value),size =1, 
             color = "black",linetype = "solid")+
  scale_y_continuous(expand = c(0,0)) +
  theme_ohchibi(y_vjust = 0)  +
  theme(
    panel.grid.major.x   = element_blank(),
    panel.grid.minor.x  = element_blank(),
    panel.grid.major.y   = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.spacing.x = unit(0.2, "lines"),
    axis.ticks.x =element_blank(),
    axis.text.y =element_text(family = "Arial",
                              face="bold",size=10,colour="black",
                              angle = 0,vjust = 0.5),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(family = "Arial",face="bold",size=12,
                                colour="black",
                                angle = 90,vjust = 1),
    legend.background = element_blank(),
    legend.position = "none",
    strip.text.x = element_blank()
  ) +
  scale_fill_manual(values = paleta_suberin) +
  ylab(label = "Suberin")



#Create a colorbar for order
df_un <- Map
df_un$Sig <- "Bar"
df_un$Strains <- df_un$Strains %>% factor(levels = mtips)
df_un <- with(df_un,order(Strains)) %>%
  df_un[.,]
df_order <- data.frame(Order = df_un$Order %>% as.character %>% unique ,Color = rep(c("black","white"),9))
paleta_order <- df_order$Color %>% as.character
names(paleta_order) <- df_order$Order %>% as.character


order_order <- df_un$Order %>% as.character %>% unique
df_un$Order <- df_un$Order %>% factor(levels = order_order)


p_colorbar <- ggplot(df_un,aes(Sig,Strains, fill = Order,color = Order)) +
  geom_tile(size = 0)+
  theme_ohchibi(size_axis_text.y = 0,size_axis_title.x = 0,
                size_axis_title.y = 0,size_axis_text.x = 0,
  ) +
  theme(axis.ticks = element_blank(), 
        axis.text.x = element_blank(),
        panel.spacing = unit(0.25, "lines"),
        axis.text.y = element_blank(),
        legend.position = "none",
        panel.grid.major.x = element_line(size = 0.4, color = "#D9D9D9"),
        panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank())+
  coord_flip(expand = FALSE) +
  scale_fill_manual(values = paleta_order) +
  scale_color_manual(values = paleta_order) 

### Create bar for significance in suberin pattern
df_sub_sig <- data.frame(Strains = merged$Strains,Estimate = merged$Suberin_est,padj  = merged$Suberin_p.adj)
df_sub_sig$Significance <- "NotSignificant"
df_sub_sig$Significance[which((df_sub_sig$padj < 0.1) & (df_sub_sig$Estimate < 0))] <- "Down"
df_sub_sig$Significance[which((df_sub_sig$padj < 0.1) & (df_sub_sig$Estimate >0))] <- "Up"
df_sub_sig <- match(mtips,df_sub_sig$Strains) %>% df_sub_sig[.,] %>% droplevels
df_sub_sig$Strains <- df_sub_sig$Strains %>% factor(levels = mtips) 
df_sub_sig$Sig <- "Bar"

p_sub_sig <- ggplot(df_sub_sig,aes(Sig,Strains, fill = Significance,color = Significance)) +
  geom_tile(size = 0)+
  theme_ohchibi(size_axis_text.y = 0,size_axis_title.x = 0,
                size_axis_title.y = 0,size_axis_text.x = 0,
  ) +
  theme(axis.ticks = element_blank(), 
        axis.text.x = element_blank(),
        panel.spacing = unit(0.25, "lines"),
        axis.text.y = element_blank(),
        legend.position = "none",
        panel.grid.major.x = element_line(size = 0.4, color = "#D9D9D9"),
        panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank())+
  coord_flip(expand = FALSE) +
  scale_fill_manual(values = paleta_prop) +
  scale_color_manual(values = paleta_prop) 

##Append the total cell information
Mapa <- Dat$Suberin
df <- Mapa[,c("Strains","Plant","Batch","Total_cells_Normvalue")] %>% unique

#Test each strain versus NB in the corresponding batch
mstrains <- df$Strains %>% levels %>% grep(pattern = "NB",invert = T,value = T)
Res <- NULL
for(st in mstrains){
  df_sub <- df %>% subset(Strains == st) %>% droplevels
  mbatches <- df_sub$Batch %>% as.character %>% unique
  df_nb <- which(df$Batch %in% mbatches) %>%
    df[.,] %>% subset(Strains == "NB") %>% droplevels
  df_both <- rbind(df_sub,df_nb)  
  m1 <- car::leveneTest(Total_cells_Normvalue ~ Strains,df_both)
  pvalue <- m1$`Pr(>F)`[1]
  x <- df_sub$Total_cells_Normvalue
  y <- df_nb$Total_cells_Normvalue
  if(pvalue < 0.05){
    m1 <- t.test(x,y,exact = F,var.equal = F)
    m_w <- wilcox.test(x,y,exact = F,var.equal = F)
  }else{
    m1 <- t.test(x,y,exact = F,var.equal = T)
    m_w <- wilcox.test(x,y,exact = F,var.equal = T)
    
  }
  pval <- m1$p.value
  temp <- data.frame(Strains = st,lavene = pvalue,
                     p.value = pval,
                     p.valueW = m_w$p.value)
  Res <- rbind(Res,temp)
}
plot(Res$p.value,Res$p.valueW)
cor.test(Res$p.value,Res$p.valueW)
cor.test(Res$p.value,Res$p.valueW,method = "spearman",exact = F)

#Stay with t.test p.value
Res$p.adj <- Res$p.value %>% p.adjust(method = "fdr")

Res$Significance <- NA
Res$Significance[which(Res$p.adj < 0.1)] <- "q < 0.1"
#Determine quantiles
df_quant <- df %>% subset(Strains == "NB") %$% Total_cells_Normvalue %>%
  quantile() 

mean_nb <- df_quant[3] %>% as.numeric

#We want to show effect of number of ceclls
#Calculate summary for the number of cells
df_plot <- df %>%
  summarySE(data = .,measurevar = "Total_cells_Normvalue",groupvars = "Strains")

df_plot <- which(df_plot$Strains %in% mtips) %>%
  df_plot[.,] %>% droplevels

df_plot$Strains <- df_plot$Strains %>% factor(levels = mtips)

#Append the test information
df_plot$Significance <- match(df_plot$Strains,Res$Strains) %>%
  Res$Significance[.]

#Plot

p_tocells <- ggplot(data = df_plot %>% subset(Strains != "NB"),aes(Strains,Total_cells_Normvalue)) +
  geom_point(size =0,stroke = NA,color = NA) +
  geom_hline(yintercept = mean_nb) +
  geom_bar(stat = "identity",aes(fill = Significance),width = 1,color = "#00000000") +
  #geom_pointrange(mapping = aes(
  #  ymin = Total_cells_Normvalue -ci,
  #  ymax = Total_cells_Normvalue + ci,color = Significance))  +
  ylab(label = "Total number of cells") +
  theme_ohchibi(font_family = "Helvetica") +
  theme(
    axis.text.y =element_text(family = "Helvetica",
                              face="bold",size=10,colour="black",
                              angle = 0,vjust = 0.5),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(family = "Helvetica",face="bold",size=12,
                                colour="black",
                                angle = 90,vjust = 1),
    axis.ticks.x = element_blank()
    
  ) +
  ylab(label = "Normalized total number of cells") +
  scale_y_continuous(limits = c(0,120),oob = rescale_none,expand = c(0,0)) +
  scale_fill_manual(values = "red",na.value = "#D9D9D9")



#Create composition 
composition <- egg::ggarrange(p_tree,p_colorbar,p_prop,p_sub,p_sub_sig,p_tocells,ncol = 1,heights = c(0.15,0.04,1,1,0.075,1))


oh.save.pdf(p = composition,outname = "screening_barriers_phylogeny.pdf",outdir = "../figures/",
            height = 10,width = 15)

#Do the phylogenetic singal
merged_tree <- match(tree$tip.label,df_plot$Strains) %>%
  df_plot[.,]
sig_pi <- phylosig(tree,merged_tree$Total_cells_Normvalue,method="lambda",test=TRUE)

sig_pi$lambda
sig_pi$P

rm(list=ls())
dev.off()
gc()
