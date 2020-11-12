library(ohchibi)
library(palettesPM)
library(Rmisc)
library(ggpubr)
library(DESeq2)
library(dendextend)
library(data.tree)
library(limma)

set.seed(130816)
setwd('/home/isai/Documents/results/rootbarriersmicro/scripts/')
dir.create("../cleandata")
dir.create("../figures")

Res <- readRDS(file = "../cleandata/deseq2_results_censusmutants_genotypes_inside_fractions.RDS")

mlevels <- c("Phylum","Class","Order","Family","Genus","ASV")


#### Soil ######
Res_sub <- Res$Soil
df_wide <- NULL
for(i in 1:6){
  df_sub <- Res_sub[i] %>% as.data.frame
  df_sub$Level <- mlevels[i]
  df_wide <- rbind(df_wide,df_sub)
}

df_wide$Fraction <- "Soil"
res_soil <- df_wide
rm(df_wide)


##### Root ############
Res_sub <- Res$Root
df_wide <- NULL
for(i in 1:6){
  df_sub <- Res_sub[i] %>% as.data.frame
  df_sub$Level <- mlevels[i]
  df_wide <- rbind(df_wide,df_sub)
}

df_wide$Fraction <- "Root"
res_root <- df_wide
rm(df_wide)

##### Shoot ############
Res_sub <- Res$Shoot
df_wide <- NULL
for(i in 1:6){
  df_sub <- Res_sub[i] %>% as.data.frame
  df_sub$Level <- mlevels[i]
  df_wide <- rbind(df_wide,df_sub)
}

df_wide$Fraction <- "Shoot"
res_shoot <- df_wide
rm(df_wide)

res_all <- rbind(res_soil,res_root,res_shoot)
res_all$Level <- res_all$Level %>% 
  factor(levels = c("ASV","Genus","Family","Order","Class","Phylum"))
res_all$padjG <- p.adjust(res_all$pvalue,method = "fdr")
res_all$Significance <- 0
res_all$Significance[which(res_all$padjG < 0.05)] <- 1

write.table(x = res_all,file = "../cleandata/deseq2_results_censusmutants_genotypes_inside_fractions_all.tsv",
            append = F,quote = F,sep = "\t",row.names = F,col.names = T)
# 
# ### Check estimate effects ########
# res_lm <- readRDS("../cleandata/res_lm_cap_amplicon_census.RDS") %$% Res_lm
# 
# 
# 
# ## correlate ###
# ## root ###
# df_plot <- res_all %>% subset(Fraction == "Root" & Level == "ASV") %>%
#   aggregate(abs(log2FoldChange )~Genotype,.,mean)
# 
# df_lm <- res_lm %>% subset(Fraction == "Root" & Axis == "CAP1")
# df_plot <- merge(df_plot,df_lm, by = "Genotype")
# 
# df_plot$EstimateAbs <- abs(df_plot$Estimate)
# colnames(df_plot)[2] <- "log2FoldChangeAbs"
# ggplot(data = df_plot,aes(EstimateaAbs,log2FoldChangeAbs)) +
#   geom_smooth(method = "lm") +
#   geom_point() +
#   theme_ohchibi() +
#   stat_cor()
# 
# ## shoot ###
# df_plot <- res_all %>% subset(Fraction == "Shoot" & Level == "ASV") %>%
#   aggregate(abs(log2FoldChange )~Genotype,.,mean)
# 
# df_lm <- res_lm %>% subset(Fraction == "Shoot" & Axis == "CAP2")
# df_plot <- merge(df_plot,df_lm, by = "Genotype")
# 
# df_plot$EstimateAbs <- abs(df_plot$Estimate)
# colnames(df_plot)[2] <- "log2FoldChangeAbs"
# ggplot(data = df_plot,aes(EstimateaAbs,log2FoldChangeAbs)) +
#   geom_smooth(method = "lm") +
#   geom_point() +
#   theme_ohchibi() +
#   stat_cor()
# 
# ## soil ###
# df_plot <- res_all %>% subset(Fraction == "Soil" & Level == "ASV") %>%
#   aggregate(abs(log2FoldChange )~Genotype,.,mean)
# 
# df_lm <- res_lm %>% subset(Fraction == "Soil" & Axis == "CAP1")
# df_plot <- merge(df_plot,df_lm, by = "Genotype")
# 
# df_plot$EstimateAbs <- abs(df_plot$Estimate)
# colnames(df_plot)[2] <- "log2FoldChangeAbs"
# ggplot(data = df_plot,aes(EstimateaAbs,log2FoldChangeAbs)) +
#   geom_smooth(method = "lm") +
#   geom_point() +
#   theme_ohchibi() +
#   stat_cor()
# 
# 
# 
# 
# ###### Create dataframe of relative abundance
# Dat <- readRDS(file = "../cleandata/dat_censusmutants_16s.clean.RDS")
# Dat_rab <- Dat$RelativeAbundance
# 
# Tab <- Dat_rab$Tab
# Tab_batch <- removeBatchEffect(x = Tab,Dat_rab$Map$Rep)
# 
# 
# Dat_batch <- create_dataset(Tab = Tab_batch,Map = Dat_rab$Map,Tax = Dat_rab$Tax)
# #remove batch effect
# 
# 
# mDats <- list()
# counter <- 0
# for(i in 3:7){
#   Dat_temp <- collapse_by_taxonomy.Dataset(Dat = Dat_batch,level = i)
#   counter <- counter+1
#   mDats[[counter]] <- Dat_temp
# }
# mDats[[6]] <- Dat_batch
# 
# ### Loop over the dats and compute average asv abundance by fraction+genotype
# taxonomic_level <- c("Phylum","Class","Order","Family","Genus","ASV")
# df_ra <- NULL
# for(i in 1:6){
#  temp <- mDats[[i]] 
#  melt_temp <- temp$Tab %>% melt
#  colnames(melt_temp) <- c("Id","DADA2_Header","RA")
#  melt_temp <- merge(melt_temp,temp$Map, by = "DADA2_Header")
#  melt_temp <- melt_temp %>% subset(Genotype == "Col_0")  %>% droplevels
#  melt_sum <- dcast(data = melt_temp,formula = Id~Fraction,fun.aggregate = mean,value.var = "RA")
#  melt_sum <- melt_sum %>% melt
#  melt_sum$Level <- taxonomic_level[i]
#  colnames(melt_sum)[2:3] <- c("Fraction","RA")
#  df_ra <- rbind(df_ra,melt_sum)
# }
# 
# 
# df_sum <- aggregate(Significance ~ Level+Genotype+Fraction,res_all,sum)
# 
# ### compute statistic per 
# Res_cumra <- NULL
# for(lev in res_all$Level %>% as.character %>% unique){
#   for(fraction in res_all$Fraction %>% as.character %>% unique){
#     for(geno in res_all$Genotype %>% as.character %>% unique){
#       mchosen <- res_all %>%
#         subset(Level == lev & Fraction == fraction & Genotype == geno &   Significance ==1)  %$%
#         Id
#       mtemp <- df_ra %>% subset(Level == lev & Fraction == fraction  ) %>% droplevels
#       mval <- mtemp$RA[which(mtemp$Id %in% mchosen)] %>% sum
#       Res_cumra <- data.frame(Level = lev,Fraction = fraction,Genotype = geno, CumRA = mval) %>%
#         rbind(Res_cumra,.)
#     }
#   }
# }
# 
# 
# df_sum$Genotype <- df_sum$Genotype %>% factor
# Res_cumra$Genotype <- Res_cumra$Genotype %>% factor(levels = df_sum$Genotype %>% levels)
# 
# 
# ggplot(data = df_sum,aes(Fraction,Significance)) +
#   geom_bar(stat = "identity",aes(fill = Fraction)) + 
#   facet_grid(.~Level) + 
#   theme_ohchibi(size_panel_border = 1,size_axis_text.x = 12) + 
#   ylab(label = "Significant tests versus Col-0") +
#   scale_fill_fraction() +
#   scale_y_continuous(expand = c(0,0)) +
#   theme(
#     strip.background.x = element_blank(),
#     strip.text.x = element_text(family = "Arial",face = "bold",size = 25),
#     axis.title.x = element_blank(),
#     legend.position = "none"
#   )
# 
# 
# 
# 
# ## Count genotypes in x axis
# ggplot(data = df_sum,aes(Genotype,Significance)) +
#   geom_bar(stat = "identity",aes(fill = Fraction)) + 
#   facet_grid(.~Level) + 
#   theme_ohchibi(size_panel_border = 1,size_axis_text.x = 12) + 
#   ylab(label = "Significant tests versus Col-0") +
#   scale_fill_fraction() +
#   scale_y_continuous(expand = c(0,0)) +
#   theme(
#     strip.background.x = element_blank(),
#     strip.text.x = element_text(family = "Arial",face = "bold",size = 25),
#     axis.title.x = element_blank(),
#     legend.position = "none",
#     axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1)
#   )  
# 
# 
# p1 <- df_sum %>% subset(Fraction == "Root" & Level == "ASV") %>%
#   ggplot(data =.,aes(Genotype,Significance)) +
#   geom_bar(stat = "identity",aes(fill = Fraction)) + 
#   facet_grid(.~Level) + 
#   theme_ohchibi(size_panel_border = 1,size_axis_text.x = 12) + 
#   ylab(label = "Significant tests versus Col-0") +
#   scale_fill_fraction() +
#   scale_y_continuous(expand = c(0,0)) +
#   theme(
#     strip.background.x = element_blank(),
#     strip.text.x = element_text(family = "Arial",face = "bold",size = 25),
#     axis.title.x = element_blank(),
#     legend.position = "none",
#     axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1)
#   )  
# 
# 
# 
# p2 <- Res_cumra %>% subset(Fraction == "Root" & Level == "ASV") %>%
#   ggplot(data = .,aes(Genotype,CumRA)) +
#   geom_bar(stat = "identity",aes(fill = Fraction)) +
#   facet_grid(.~Level) + 
#   theme_ohchibi(size_panel_border = 1,size_axis_text.x = 12) + 
#   ylab(label = "Cumulative Relative Abundance") +
#   scale_fill_fraction() +
#   scale_y_continuous(expand = c(0,0)) +
#   theme(
#     strip.background.x = element_blank(),
#     strip.text.x = element_text(family = "Arial",face = "bold",size = 25),
#     axis.title.x = element_blank(),
#     legend.position = "none",
#     axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1)
#   )  
#   
# 
# 
# 
# Res_cumra %>% subset(Fraction == "Shoot" & Level == "ASV") %>%
#   ggplot(data = .,aes(Genotype,CumRA)) +
#   geom_bar(stat = "identity",aes(fill = Fraction)) +
#   facet_grid(.~Level) + 
#   theme_ohchibi(size_panel_border = 1,size_axis_text.x = 12) + 
#   ylab(label = "Cumulative Relative Abundance") +
#   scale_fill_fraction() +
#   scale_y_continuous(expand = c(0,0)) +
#   theme(
#     strip.background.x = element_blank(),
#     strip.text.x = element_text(family = "Arial",face = "bold",size = 25),
#     axis.title.x = element_blank(),
#     legend.position = "none",
#     axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1)
#   )  
# 
# ####### Phylogenetic signal test #########
# 
# ####### Integrate relative abundance information #######
# Dat <- readRDS(file = "../cleandata/dat_censusmutants_16s.clean.RDS")
# Dat_rab <- Dat$RelativeAbundance
# df_tax <- Dat$df_tax
# 
# df_res <- df_tax
# #Clean for weird characters so the tree costruction works
# df_res$Genus <- df_res$Genus %>%
#   gsub(pattern = "\\/",replacement = "") 
# 
# 
# df_res$Phylum <- paste0("p_",df_res$Phylum)
# df_res$Class <- paste0(df_res$Phylum,"|c_",df_res$Class)
# df_res$Order <- paste0(df_res$Class,"|o_",df_res$Order)
# df_res$Family <- paste0(df_res$Order,"|f_",df_res$Family)
# df_res$Genus <- paste0(df_res$Family,"|g_",df_res$Genus)
# df_res$ASV_Id_Long <- paste0(df_res$Genus,"|",df_res$Id)
# df_res$ASV_Id_Long <- df_res$ASV_Id_Long %>% 
#   gsub(pattern = "\\(|\\)",replacement = "")
# df_res$pathString <- paste0("Root/",df_res$Phylum,"/",df_res$Class,"/",
#                             df_res$Order,"/",df_res$Family,"/",df_res$Genus,"/",
#                             df_res$ASV_Id_Long)
# df_res$pathString <- df_res$pathString %>% 
#   gsub(pattern = "\\(|\\)",replacement = "")
# 
# df_res <- which(df_res$Id %in% rownames(Dat_rab$Tab )) %>%
#   df_res[.,] %>% droplevels
# 
# population <- as.Node(df_res)
# 
# tree <- data.tree::as.phylo.Node(population)
# 
# melted_rab <- Dat_rab$Tab %>% melt
# colnames(melted_rab) <- c("Id","DADA2_Header","RA")
# melted_rab <- merge(melted_rab,Dat_rab$Map, by = "DADA2_Header")
# 
# 
# 
# tree$tip.label <- tree$tip.label %>% gsub(pattern = ".*\\|",replacement = "") %>% as.character
# 
# res_sub <- res_all %>% subset(Level == "ASV") %>% droplevels
# Tab_psig <- acast(data = res_sub,formula = Id ~Fraction+Genotype,
#                   value.var = "log2FoldChange")
# 
# Res_psignal <- NULL
# for(i in 1:ncol(Tab_psig)){
#   x <- na.omit(Tab_psig[,i])   
#   todrop <- which(!(tree$tip.label %in% names(x))) %>% tree$tip.label[.]
#   tree_sub <- ape::drop.tip(phy = tree,tip = todrop)
#   rp <- phylosig(tree = tree_sub,x,method = "lambda",test = T)
#   Res_psignal <- data.frame(lambda = rp$lambda,p.value = rp$P,contrast = colnames(Tab_psig)[i]) %>%
#     rbind(Res_psignal,.)
#   
# }
# Res_psignal$Fraction <- Res_psignal$contrast %>% gsub(pattern = "_.*",replacement = "")
# Res_psignal$Genotype <- Res_psignal$contrast %>% gsub(pattern = "Root_|Shoot_|Soil_",replacement = "")
# Res_psignal$p.adj <- p.adjust(Res_psignal$p.value,method = "fdr")
# Res_psignal$Significance <- "NS"
# Res_psignal$Significance[which(Res_psignal$p.adj < 0.1)] <- "Significant"
# 
# ggplot(data = Res_psignal,aes(Genotype,lambda)) +
#   geom_bar(stat = "identity",aes(fill  = Significance)) +
#   facet_grid(.~Fraction) 
# with(Res_psignal,order(-lambda)) %>% Res_psignal[.,]
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
