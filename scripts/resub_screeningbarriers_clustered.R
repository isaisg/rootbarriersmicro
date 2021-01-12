library(ohchibi)
library(Rmisc)
library(scales)

setwd('/Users/isaisalasgonzalez/Documents/rootbarriersmicro/scripts/')

set.seed(130816)

Dat <- readRDS(file = "../cleandata/dat_screeningbarriers.RDS")

Mapa <- Dat$Suberin
df <- Mapa[,c("Strains","Plant","Batch","Total_cells_Normvalue")] %>% unique

df_quant <- df %>% subset(Strains == "NB") %$% Total_cells_Normvalue %>%
  quantile() 

#We want to show effect of number of ceclls
#Calculate summary for the number of cells
df_plot <- df %>%
  summarySE(data = .,measurevar = "Total_cells_Normvalue",groupvars = "Strains")


#Order the cells according to both PI and suberin
df_cont <- aggregate(Normvalue_ra ~ Strains,
                     Dat$Suberin %>% subset(Localization  =="Continous_expression" ),mean)

mtips <- with(df_cont,order(-Normvalue_ra)) %>% df_cont[.,] %$% Strains %>% as.character 
df_plot$Strains <- df_plot$Strains %>% factor(levels = mtips)

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

#Read the result from 
merged <- read.table(file = "../../rootbarriersmicro/cleandata/res_merged_propidium_suberin.tsv",header = T,sep = "\t")

#Palettes for paleta de sisifo

paleta_prop <- c("#CCCCCC","red","red")
names(paleta_prop) <- c("NotSignificant","Up","Down")

mean_nb_prop <- merged %>% subset(Strains == "NB") %$% Propidium_Mean
merged_tree <- merged %>% subset(Strains != "NB")


merged_tree$SignificanceProp <- "NotSignificant"
merged_tree$SignificanceProp[which((merged_tree$Propidium_p.adj < 0.1) & (merged_tree$Propidium_est > 0 ))] <- "Up"
merged_tree$SignificanceProp[which((merged_tree$Propidium_p.adj < 0.1) & (merged_tree$Propidium_est < 0 ))] <- "Down"

merged_tree$Strains <- merged_tree$Strains %>% factor(levels = mtips)

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
    axis.title.y = element_text(family = "Helvetica",
                                face="bold",size=12,colour="black",
                                angle = 90,vjust = 1),
    axis.text.y =element_text(family = "Helvetica",
                              face="bold",size=10,colour="black",
                              angle = 0,vjust = 0.5),
    strip.text.x = element_blank(),
    panel.spacing.x = unit(0.2, "lines")
  )   +
  scale_fill_manual(values =paleta_prop) +
  ylab(label = "Propidium")


#Suberin plot 
### suberin

melted_suberin <- Dat$Suberin %>% dcast(data = .,formula =Strains~Localization,
                                        fun.aggregate = mean,value.var = "Normvalue_ra" )  %>%
  melt
nb_sub <- melted_suberin %>% subset(Strains == "NB") %>% droplevels
melted_sub <- melted_suberin %>% subset(Strains != "NB") %>% droplevels
melted_sub$Strains <- melted_sub$Strains %>% factor(levels = mtips)
melted_sub$variable <- melted_sub$variable %>% 
  factor(levels = c("Continous_expression","Discrete_expression","No_expression"))


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
    axis.text.y =element_text(family = "Helvetica",
                              face="bold",size=10,colour="black",
                              angle = 0,vjust = 0.5),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(family = "Helvetica",face="bold",size=12,
                                colour="black",
                                angle = 90,vjust = 1),
    legend.background = element_blank(),
    legend.position = "none",
    strip.text.x = element_blank()
  ) +
  scale_fill_manual(values = paleta_suberin) +
  ylab(label = "Suberin")


### Create bar for significance in suberin pattern
df_sub_sig <- data.frame(Strains = merged$Strains,Estimate = merged$Suberin_est,padj  = merged$Suberin_p.adj)
df_sub_sig$Significance <- "NotSignificant"
df_sub_sig$Significance[which((df_sub_sig$padj < 0.1) & (df_sub_sig$Estimate < 0))] <- "Down"
df_sub_sig$Significance[which((df_sub_sig$padj < 0.1) & (df_sub_sig$Estimate >0))] <- "Up"
df_sub_sig <- match(mtips,df_sub_sig$Strains) %>% df_sub_sig[.,] %>% droplevels
df_sub_sig <- df_sub_sig %>% subset(Strains != "NB")
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



##add a bar of chosen strains
st <- read.table(file = "../../rootbarriersmicro/rawdata//41_strains_selected.txt",stringsAsFactors = F) %$% V1

merged_tree$Chosen <- "No"
merged_tree$Chosen[which(merged_tree$Strains %in% st)] <- "Yes"
merged_tree$Bar <- "Bar"
p_bar <- ggplot(data = merged_tree,aes(Strains,Bar)) +
  geom_tile(aes(fill = Chosen,color = Chosen)) +
  scale_color_manual(values = c("black","red"))+
  scale_fill_manual(values = c("black","red")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none"
  ) +
  scale_y_discrete(expand = c(0,0))

#Create composition adding the number of cells distribution
composition <- ggarrange(p_bar,p_prop,p_sub,p_sub_sig,p_tocells + theme(legend.position = "none"),ncol = 1,heights = c(0.05,1,1,0.075,1))

oh.save.pdf(p = composition,outname = "revision_screening_barriers_sisifo.pdf",outdir = "../figures/",
            height = 12,width = 15)

rm(list=ls())
dev.off()
