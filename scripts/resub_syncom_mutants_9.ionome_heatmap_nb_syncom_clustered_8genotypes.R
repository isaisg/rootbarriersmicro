library(ohchibi)
library(emmeans)
library(paletteer)
library(scales)
library(ggtree)

set.seed(130816)
setwd('/Users/isaisalasgonzalez/Documents/rootbarriersmicro/scripts/')

Dat <- readRDS(file = "../cleandata/dat_syncom_genotypes.RDS") %$%
  Ionome

paleta_enrichment <- readRDS(file = "../cleandata/dat_syncom_stresses.RDS") %$%
  paleta_enrichment

#Perform testing inside each ion
#We wanna do pairwise comparison between  StressNB vs Full NB
#StressSC vs Full SC and #Intra Stress SC vs Stress NB
Tab_z <- Dat$Tab %>% t %>% scale %>% t 
melted <- Tab_z %>% melt
colnames(melted) <- c("Ion","Index","zscore")

melted <- merge(melted,Dat$Map, by = "Index") 
melted$Genotype <- melted$Genotype %>%
  factor(levels = c("Col-0","pCASP1::CDEF","esb1.1pCASP1::CDEF","sgn3.3myb36pELTP::CEDF",
                    "esb1","sgn3","myb36","sgn3myb36",
                    "aba2.1","abi4.1","etr1.1","ein3.1"))

melted$Type <- melted$Type %>% factor(levels = c("NB","HK","SynCom"))

#Do a model per ion
Res_em <- NULL
Res_p <- NULL
for(ion in melted$Ion %>% unique){
  melted_sub <- melted %>% subset(Ion  == ion) %>% droplevels
  m1 <- lm(data = melted_sub,formula = zscore ~ Type + Genotype + Type:Genotype)
  m1_res <- emmeans(m1,pairwise~Type:Genotype,adjust = "none")
  m1_em <- m1_res$emmeans %>% as.data.frame
  m1_p <- m1_res$contrasts %>% as.data.frame
  m1_em$Ion <- ion
  m1_p$Ion <- ion
  Res_em <- rbind(Res_em,m1_em)
  Res_p <- rbind(Res_p,m1_p)
}

#Arrange the pvalues dataframe
Res_p$contrast <- Res_p$contrast %>% gsub(pattern = ",",replacement = "_") %>%
  gsub(pattern = " ",replacement = "")

Res_p$contrast %>% as.character %>% unique
#Determine the pairwise contrasts
c1 <- c("(NBCol-0)-(SynComCol-0)",
        "NBpCASP1::CDEF-SynCompCASP1::CDEF",
        "NBesb1.1pCASP1::CDEF-SynComesb1.1pCASP1::CDEF",
        "NBsgn3.3myb36pELTP::CEDF-SynComsgn3.3myb36pELTP::CEDF",
        "NBesb1-SynComesb1",
        "NBsgn3-SynComsgn3",
        "NBmyb36-SynCommyb36",
        "NBsgn3myb36-SynComsgn3myb36",
        "NBaba2.1-SynComaba2.1",
        "NBabi4.1-SynComabi4.1",
        "NBetr1.1-SynCometr1.1",
        "NBein3.1-SynComein3.1"
)

#Now the contrast against Full in each type of syncom
c9 <- Res_p$contrast %>% grep(pattern = "NB.*NB",value = T) %>%
  grep(pattern = "Col-0",value = T) %>% unique

c11 <- Res_p$contrast %>% grep(pattern = "SynCom.*SynCom",value = T) %>%
  grep(pattern = "Col-0",value = T) %>% unique

allc <- c(c1,c9,c11)
allc <- allc %>% unique

Res_p_sub <- which(Res_p$contrast %in% allc) %>%
  Res_p[.,] %>% droplevels

Res_p_sub$p.adj <- Res_p_sub$p.value %>% p.adjust(method = "fdr")
Res_p_sub$Significance <- "NoSignificant"
pval_thres <- 0.05
Res_p_sub$Significance[which(Res_p_sub$p.adj < pval_thres)] <- "Significant"


#Against the full 
agfull <- c(c9,c11)
Res_agfull <- which(Res_p_sub$contrast %in% agfull) %>%
  Res_p_sub[.,] %>% droplevels

df_temp <- Res_agfull$contrast %>% gsub(pattern = "NBCol-0|HKCol-0|SynComCol-0",replacement = "") %>%
  gsub(pattern = "\\(\\)-",replacement = "") %>% as.data.frame
colnames(df_temp) <- c("UId")
Res_agfull <- cbind(Res_agfull,df_temp)
Res_agfull$UId <- paste(Res_agfull$UId,Res_agfull$Ion,sep = "_")
df_temp <- data.frame(UId = Res_agfull$UId,
                      p.adj_agfull =Res_agfull$p.adj, estimate_agfull = Res_agfull$estimate)
df_temp$Significance_agfull <- "NoSignificant"
df_temp$Significance_agfull[which(df_temp$p.adj_agfull < pval_thres & df_temp$estimate_agfull > 0)] <- "Down"
df_temp$Significance_agfull[which(df_temp$p.adj_agfull < pval_thres & df_temp$estimate_agfull < 0)] <- "Up"


#Structure to plot
Res_em$UId <- paste0(Res_em$Type,Res_em$Genotype,"_",Res_em$Ion)
Res_em <- merge(Res_em,df_temp, by = "UId",all.x = TRUE)


#Analysis for the Syncoms vs NB comparison
Res_em_sub <- Res_em %>% subset(Type != "HK") %>% droplevels

#Order the stresses based on level of extremeness
Res_em_sub$Genotype <- Res_em_sub$Genotype %>%
  factor(levels = c("Col-0","pCASP1::CDEF","esb1.1pCASP1::CDEF","sgn3.3myb36pELTP::CEDF",
                    "esb1","sgn3","myb36","sgn3myb36",
                    "aba2.1","abi4.1","etr1.1","ein3.1"))

#Append the intra stress comparison
intrast <- c(c1) %>% unique
Res_intrast <- which(Res_p_sub$contrast %in% intrast) %>%
  Res_p_sub[.,] %>% droplevels
df_temp <- Res_intrast$contrast %>% grep(pattern = "NB.*SynCom",value = F) %>%
  Res_intrast[.,] %>% droplevels

df_temp$Genotype <- df_temp$contrast %>% gsub(pattern = ".*SynCom",replacement = "") %>%
  gsub(pattern = "\\)$",replacement = "")
df_temp$SI <- paste(df_temp$Genotype,df_temp$Ion,sep = "_")



df_temp <- data.frame(SI = df_temp$SI,p.value_intra = df_temp$p.value,
                      p.adj_intra =df_temp$p.adj, estimate_intra = df_temp$estimate)
df_temp$Significance_intra <- "NoSignificant"
df_temp$Significance_intra[which(df_temp$p.adj_intra < pval_thres)] <- "Significant"


Res_em_sub$SI <- paste(Res_em_sub$Genotype,Res_em_sub$Ion,sep = "_")

Res_em_sub <- merge(Res_em_sub,df_temp, by = "SI", all.x = TRUE)


Res_em_sub$Significance_agfull <- Res_em_sub$Significance_agfull %>%
  factor(levels = c("NoSignificant","Down","Up"))

Res_em_sub$Significance_intra_s <- Res_em_sub$Significance_intra 
Res_em_sub$Significance_intra_s <- Res_em_sub$Significance_intra_s %>%
  gsub(pattern = "NoSignificant",replacement = NA) %>%
  gsub(pattern = "Significant",replacement = "*")


##### Subset only the genotypes used for the 8 plots
chosen_genos <- c("Col-0","esb1","esb1.1pCASP1::CDEF","myb36",
  "pCASP1::CDEF","sgn3","sgn3.3myb36pELTP::CEDF","sgn3myb36")
Res_em_sub <- which(Res_em_sub$Genotype %in% chosen_genos) %>%
  Res_em_sub[.,] %>% droplevels


#Determine the order of the ions clustering them based on their z-score pattern
Tab <- acast(data = Res_em_sub,formula =Ion~Type+Genotype,
             value.var = "emmean")
#mclust_ions <- dist(Tab) %>% hclust(method = "ward.D") 
mclust_ions <- hclust(d = 1-cor(Tab %>% t) %>% as.dist,method = "complete")

order_ions <- mclust_ions$order %>% mclust_ions$labels[.]

Res_em_sub$Ion <- Res_em_sub$Ion %>% factor(levels = order_ions)


### Order the genotypes ###
mclust_geno <- hclust(d = 1-cor(Tab ) %>% as.dist,method = "complete")
order_geno <- mclust_geno$order %>% mclust_geno$labels[.]

Res_em_sub$Ion <- Res_em_sub$Ion %>% factor(levels = order_ions)


Res_em_sub$UId <- paste0(Res_em_sub$Type,"_",Res_em_sub$Genotype)

Res_em_sub$UId <- Res_em_sub$UId %>%
  factor(levels = order_geno)

p_heatmap <- ggplot(data = Res_em_sub,aes(UId,Ion)) +
  geom_raster(aes(fill = emmean)) +  
  geom_tile(aes(color = Significance_agfull),fill = '#00000000', size = 1.5,width = 0.85,height = 0.85) + 
  scale_fill_paletteer_c("pals::kovesi.diverging_bwr_55_98_c37",
                         limits = c(-2.5,2.5),oob = squish) +
  #geom_text(aes(label = Significance_intra_s),size = 8) + 
  #geom_text(aes(label = Significance_intra_s),size = 8.2) + 
  #geom_text(aes(label = Significance_intra_s),size = 8.4) + 
  scale_color_manual(values = c("#00000000","black","black"),na.value = "#00000000")+
  theme_ohchibi(size_axis_text.x = 20,
                angle_text.x = 90,
                size_axis_text.y = 20,
                size_axis_title.x = 22,
                size_axis_title.y = 0,
                legend_proportion_size = 1,
                size_title_text = 12,
                size_legend_text = 12,
                size_panel_border = 1.5,
                size_lines_panel = 0) +
  theme(axis.ticks = element_blank(),
        panel.spacing = unit(0.075, "lines"),
        #legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(family = "Arial",face = "bold",size = 20),
        axis.text.y = element_text(family = "Arial",face = "bold",hjust = 1,vjust = 0.5),
        axis.text.x = element_text(hjust = 1,vjust = 0.5,family = "Arial",face = "bold"),
        axis.title.x = element_blank()
  )   +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0))

#Prepare the dendrogram to put in the same graph
tree <- mclust_ions %>% as.phylo
p_tree_ions <- ggtree(tree,ladderize = F,size = 1.2) + 
  #The expand parth is fundamental to make the tree fit the composition with the other plots
  scale_y_continuous(expand = c(0.026,0.026))

tree <- mclust_geno %>% as.phylo 

p_tree_geno <- ggtree(tree,ladderize = F,size = 1.2)  +
  coord_flip() + scale_x_reverse() +
  scale_y_continuous(expand = c(0.026,0.026))


p_blank <- ggplot() + theme_void()

#Create arranged composition
composition <- egg::ggarrange(
  p_blank,p_tree_geno,
  p_tree_ions,p_heatmap,
  ncol = 2,nrow = 2,byrow = T,
  widths = c(0.1,1),
  heights = c(0.1,1)
)


oh.save.pdf(p = composition,outname = "heatmap_ionome_syncom_mutants.clustered.8genotypes.pdf",
            outdir = "../figures/",height = 12,width = 14)

### Version ordered by number an dgroups 
## Load the number scheme
df_num <- read.table(file = "../rawdata/df_genotypes_to_classification.csv",header = T,sep = ",")
df_num$Number <- 1:20
df_num$Genotype
Res_em_sub$Genotype <- Res_em_sub$Genotype %>% 
  gsub(pattern = "^Col-0$",replacement = "Col_0") %>%
  gsub(pattern = "^esb1$",replacement = "esb1_1") %>%
  gsub(pattern = "^esb1.1pCASP1::CDEF$",replacement = "pCASP1_CDEFesb1_1") %>%
  gsub(pattern = "^myb36$",replacement = "myb36_2") %>%
  gsub(pattern = "^pCASP1::CDEF$",replacement = "pCASP1_CDEFwt") %>%
  gsub(pattern = "^sgn3$",replacement = "sgn3_3") %>%
  gsub(pattern = "^sgn3.3myb36pELTP::CEDF$",replacement = "pELTP_CDEFsgn3_3_myb36_1") %>%
  gsub(pattern = "^sgn3myb36$",replacement = "myb36_2_sgn3_3")
Res_em_sub$Number <- match(Res_em_sub$Genotype,df_num$Genotype) %>%
  df_num$Number[.]
Res_em_sub$Group <- match(Res_em_sub$Genotype,df_num$Genotype) %>%
  df_num$Group[.]


Res_em_sub$Number <- Res_em_sub$Number %>% factor
Res_em_sub$Group <- Res_em_sub$Group %>% factor

p_heatmap_two <- ggplot(data = Res_em_sub,aes(Number,Ion)) +
  geom_raster(aes(fill = emmean)) +  
  geom_tile(aes(color = Significance_agfull),fill = '#00000000', size = 1.5,width = 0.85,height = 0.85) + 
  scale_fill_paletteer_c("pals::kovesi.diverging_bwr_55_98_c37",
                         limits = c(-2.5,2.5),oob = squish) +
  facet_grid(.~Type,space = "free",scales = "free") +
  scale_color_manual(values = c("#00000000","black","black"),na.value = "#00000000")+
  theme_ohchibi(size_axis_text.x = 20,
                angle_text.x = 90,
                size_axis_text.y = 20,
                size_axis_title.x = 22,
                size_axis_title.y = 0,
                legend_proportion_size = 1,
                size_title_text = 12,
                size_legend_text = 12,
                size_panel_border = 1.5,
                size_lines_panel = 0) +
  theme(axis.ticks = element_blank(),
        panel.spacing = unit(0.075, "lines"),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(family = "Arial",face = "bold",size = 20),
        axis.text.y = element_text(family = "Arial",face = "bold",hjust = 1,vjust = 0.5),
        axis.text.x = element_text(hjust = 1,vjust = 0.5,family = "Arial",face = "bold"),
        axis.title.x = element_blank()
  )   +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0))


composition <- egg::ggarrange(
  p_tree_ions,p_heatmap_two,
  ncol = 2,nrow = 1,byrow = T,
  widths = c(0.1,1)
)


oh.save.pdf(p = composition,outname = "heatmap_ionome_syncom_mutants.parition.8genotypes.pdf",
            outdir = "../figures/",height = 12,width = 14)



rm(list=ls())
dev.off()
gc()
