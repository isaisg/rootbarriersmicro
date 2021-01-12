library(ohchibi)
library(emmeans)
library(egg)
library(car)
library(ggrepel)
library(paletteer)
#Merge datasets into a single structure
set.seed(130816)
setwd('/home/isai/Documents/results/rootbarriersmicro/scripts/')

############# Suberin ####################
df_suberin <- read.table(file = "../rawdata/data_syncomstress_suberin.tsv",header = T,sep = "\t")

df_suberin$Batch <- paste0(df_suberin$Batch,"_",df_suberin$Rep)
#df_suberin <- na.omit(df_suberin)

#Calculane normalization factor  for the weight of the plants
#Probably worth having two normalization factors one for NB and one for SC
#Use mean to normalize inside each batch
df_suberin %>% subset(Bacteria == "NB" & Stress == "Full") %>%
  chibi.boxplot(Map = .,x_val = "Batch", y_val = "Suberin",style = "open")

#df_sub <- df_suberin %>% subset(Bacteria == "NB") %>% droplevels
df_sub <- df_suberin
overall_me <- mean(subset(df_sub,Bacteria == "NB" &  Stress == "Full")$Suberin)
df_sub_norm <- NULL
for(batch in unique(df_sub$Batch)){
  local_me <- subset(df_sub,Bacteria == "NB" & Batch == batch & Stress == "Full") %$% Suberin %>% mean 
  nfactor <- overall_me/local_me
  temp <- subset(df_sub,Batch==batch)
  temp$NormSuberin <- temp$Suberin*nfactor
  df_sub_norm <- rbind(df_sub_norm,temp)
}

df_sub_norm %>% subset(Bacteria == "NB" & Stress == "Full") %>% 
  chibi.boxplot(Map = .,x_val = "Batch",y_val = "NormSuberin",
                style = "open")  + 
  theme(axis.text.x = element_text(angle = 90))

df_suberin_norm <- df_sub_norm

############# Primary root elongation ####################
df_suberin <- read.table(file = "../rawdata/data_syncomstress_primaryroot.csv",header = T,sep = "\t")

df_suberin$Stress <- df_suberin$Media %>%
  gsub(pattern = "MS",replacement =  "")  %>%
  gsub(pattern = "-Full",replacement =  "Full")  %>%
  factor %>% relevel(ref = "Full")

df_suberin$Batch <- paste0(df_suberin$Batch,"_",df_suberin$Rep)

#Calculane normalization factor  for the weight of the plants
#Probably worth having two normalization factors one for NB and one for SC
#Use mean to normalize inside each batch
df_suberin %>% subset(Bacteria == "NB" & Stress == "Full") %>%
  chibi.boxplot(Map = .,x_val = "Batch", y_val = "main.root.elongation",style = "open")

#df_sub <- df_suberin %>% subset(Bacteria == "NB") %>% droplevels
df_sub <- df_suberin
overall_me <- mean(subset(df_sub,Bacteria == "NB" &  Stress == "Full")$main.root.elongation)
df_sub_norm <- NULL
for(batch in unique(df_sub$Batch)){
  local_me <- subset(df_sub,Bacteria == "NB" & Batch == batch & Stress == "Full") %$% main.root.elongation %>% mean 
  nfactor <- overall_me/local_me
  temp <- subset(df_sub,Batch==batch)
  temp$NormRoot <- temp$main.root.elongation*nfactor
  df_sub_norm <- rbind(df_sub_norm,temp)
}

df_sub_norm %>% subset(Bacteria == "NB" & Stress == "Full") %>% 
  chibi.boxplot(Map = .,x_val = "Batch",y_val = "NormRoot",
                style = "open")  + 
  theme(axis.text.x = element_text(angle = 90))

df_root_norm <- df_sub_norm


############### Dry weight ##############
df <- read.table(file = "../rawdata/data_syncomstress_dryweight.csv",header = T,sep = "\t")
df$Type <- "SC"
df$Type[df$Sample_name %>% grep(pattern = "No Bacteria")] <- "NB"
df$Type[df$Sample_name %>% grep(pattern = "HK")] <- "HK"

df$Stress <- "Full"
df$Stress[df$Sample_name %>% grep(pattern = "NaCl")] <- "+NaCl"
df$Stress[df$Sample_name %>% grep(pattern = "-Fe")] <- "-Fe"
df$Stress[df$Sample_name %>% grep(pattern = "-Mn")] <- "-Mn"
df$Stress[df$Sample_name %>% grep(pattern = "-Zn")] <- "-Zn"
df$Stress[df$Sample_name %>% grep(pattern = "-Mo")] <- "-Mo"
df$Stress[df$Sample_name %>% grep(pattern = "-P")] <- "-Pi"
df$Stress[df$Sample_name %>% grep(pattern = "-K")] <- "-K"

df$Type <- df$Type %>% factor
df$Stress <- df$Stress %>% factor %>% relevel(ref = "Full")


#Calculane normalization factor  for the weight of the plants
#Probably worth having two normalization factors one for NB and one for SC
#Use mean to normalize inside each batch
df_suberin <- df
colnames(df_suberin)[10] <- "Bacteria"
colnames(df_suberin)[3] <- "Batch"
df_suberin$Batch <- paste0("Batch",df_suberin$Batch)

df_suberin <- df_suberin %>%subset(Bacteria != "HK") %>% droplevels
df_suberin %>% subset(Bacteria == "NB" & Stress == "Full") %>%
  chibi.boxplot(Map = .,x_val = "Batch", y_val = "weight_mg_plant",style = "open")

#df_sub <- df_suberin %>% subset(Bacteria == "NB") %>% droplevels
df_sub <- df_suberin
overall_me <- mean(subset(df_sub,Bacteria == "NB" &  Stress == "Full")$weight_mg_plant)
df_sub_norm <- NULL
for(batch in unique(df_sub$Batch)){
  local_me <- subset(df_sub,Bacteria == "NB" & Batch == batch & Stress == "Full") %$% weight_mg_plant %>% mean 
  nfactor <- overall_me/local_me
  temp <- subset(df_sub,Batch==batch)
  temp$Normweight_mg_plant <- temp$weight_mg_plant*nfactor
  df_sub_norm <- rbind(df_sub_norm,temp)
}

df_sub_norm %>% subset(Bacteria == "NB" & Stress == "Full") %>%
  chibi.boxplot(Map = .,x_val = "Batch", y_val = "Normweight_mg_plant",style = "open")


## models ###
df_dw_norm <- df_sub_norm

df_dw_norm$Stress <- df_dw_norm$Stress %>% gsub(pattern = " ",replacement = "")
df_dw_norm$Bacteria <- df_dw_norm$Bacteria %>% gsub(pattern = " ",replacement = "")

########### Ionome ################
#Read the ionome matrix
df_ionome <- read.table("../rawdata/data_syncomstress_ionome.csv",header = T,sep = ",",quote = "",comment.char = "")
Map <- df_ionome[,1:6]
Tab <- df_ionome[,7:ncol(df_ionome)]

Map$Sample.Id <- Map$Sample.Id %>% gsub(pattern = " ",replacement = "_")
colnames(Map)[2] <- "Sample_name"
Map$Sample_name <- Map$Sample_name %>% gsub(pattern = " ",replacement = "") %>%
  gsub(pattern = "NoBacteria",replacement = "NB")
colnames(Map)[3] <- "Batch"
Map$Batch <- paste0("Batch",Map$Batch)
Map$Plates <- paste0("Plate",Map$Plates)
Map$Replicates <- paste0("Rep",Map$Replicates )
Map$Run <- paste0("Run",Map$Run)
colnames(Map)[1] <- c("Sample_Id")
Map$Index <- rownames(Tab)


Tab <- match(Map$Index,rownames(Tab)) %>%
  Tab[.,]


rownames(Tab) <- Map$Index
colnames(Tab) <- colnames(Tab) %>% gsub(pattern = "\\.",replacement = "")
Tab <- as.matrix(Tab)


#
melted <- Map  %$% Index %>%
  match(rownames(Tab)) %>% Tab[.,]  %>% melt
colnames(melted) <- c("Index","Ion","value")
melted <- merge(melted,Map, by = "Index")
colnames(melted)[5] <- c("Strains")

melted %>% subset( Ion == "Na23" & Strains == "FullNB") %>% droplevels %>%
  ggplot(.,aes(Batch,value)) + 
  facet_wrap(facets = "Ion") + geom_boxplot()  + geom_point()

#There exists variation associated with biological variation besides the measurement of the ionome
#We should adjust for this variation
melted <- melted %>% subset(Strains != "NB" & Strains != "L22" & Strains != "L22HK") %>%
  droplevels

#Do it across all Ions
melted_normalized <- NULL
for(ion in melted$Ion %>% unique){
  temp <- melted %>% subset(Ion == ion) %>% droplevels
  temp_nb <- temp %>% subset(Strains == "FullNB") %>% droplevels
  overall_me <- temp_nb$value %>% mean
  for(exp in temp$Batch %>% unique){
    local_me <- subset(temp, Batch == exp & Strains == "FullNB") %$% value %>% mean 
    nfactor <- overall_me/local_me
    ftemp <- temp %>% subset(Batch == exp) %>% droplevels
    ftemp$NormValue <- ftemp$value*nfactor
    melted_normalized <- rbind(melted_normalized,ftemp)
  }
}

melted_normalized %>% subset( Ion == "Na23" & Strains == "FullNB") %>% droplevels %>%
  ggplot(.,aes(Batch,NormValue)) + 
  facet_wrap(facets = "Ion") + geom_boxplot()  + geom_point()

Tab <- acast(data = melted_normalized,formula =Index ~ Ion,value.var = "NormValue" ,fill = 0)
Map <- match(rownames(Tab),rownames(Map)) %>%
  Map[.,]
Dat <- create_dataset(Tab = Tab %>% t,Map = Map)

Dat$Map$Stress <- Dat$Map$Sample_name %>% gsub(pattern = "NB",replacement = "") %>%
  gsub(pattern = "SynCom",replacement = "") %>%
  gsub(pattern = "HK",replacement = "") %>%
  gsub(pattern = "-",replacement = "")

Dat$Map$Type <- rep("NB")
Dat$Map$Type[Dat$Map$Sample_name %>% grep(pattern = "HK")] <- "HK"
Dat$Map$Type[Dat$Map$Sample_name %>% grep(pattern = "HK|NB",value = F,invert = T)] <- "SynCom"


Dat$Map$Stress <- Dat$Map$Stress %>% 
  gsub(pattern = "NaCl",replacement = "+NaCl") %>%
  gsub(pattern = "Zn",replacement = "-Zn") %>%
  gsub(pattern = "Mn",replacement = "-Mn") %>%
  gsub(pattern = "Fe",replacement = "-Fe") %>% 
  gsub(pattern = "Mo",replacement = "-Mo") %>%
  gsub(pattern = "P",replacement = "-P") %>%
  gsub(pattern = "K",replacement = "-K")  %>%
  factor %>% relevel(ref = "Full")

#Define palettes
paleta_stress <- paletteer_d(package = "ochRe",palette = "namatjira_qual" )[1:8]
names(paleta_stress) <- levels(Dat$Map$Stress)

Dat$Map$Type <- Dat$Map$Type %>% factor(levels = c("NB","HK","SynCom"))
paleta_type <- c("#525252","#D9D9D9","#b48c36")
names(paleta_type) <- levels(Dat$Map$Type)

paleta_enrichment <- c("#B3B3B3","#2480FF","#E65037")
names(paleta_enrichment) <- c("NoSignificant","Down","Up")

#Remove spaces
df_suberin_norm$Stress <- df_suberin_norm$Stress %>% 
  gsub(pattern = " ",replacement = "") %>%
  factor %>% relevel(ref = "Full")
df_suberin_norm$Bacteria <- df_suberin_norm$Bacteria %>% gsub(pattern = "SC",replacement = "SynCom") %>% factor

df_dw_norm$Stress <- df_dw_norm$Stress %>% factor %>%relevel(ref = "Full")
df_dw_norm$Bacteria <- df_dw_norm$Bacteria %>% 
  gsub(pattern = "SC",replacement = "SynCom") %>% factor

## normalize suberin by root
df_root_norm$Stress <- df_root_norm$Stress %>% 
  gsub(pattern = " ",replacement = "") %>%
  factor %>% relevel(ref = "Full")

df_root_norm$Bacteria <- df_root_norm$Bacteria %>% 
  gsub(pattern = "SC",replacement = "SynCom") %>%
  factor

df_av_root <- aggregate(NormRoot~Bacteria+Stress+Batch,
          df_root_norm,mean)


df_suberin_norm <- merge(df_suberin_norm,df_av_root,
      by = c("Bacteria","Stress","Batch"),all.x= TRUE) 
df_suberin_norm$PropNorm <- df_suberin_norm$NormSuberin/df_suberin_norm$NormRoot



#Save the structure with all the datasets
mDat <- list(Suberin = df_suberin_norm,
             Root = df_root_norm,
             DryWeight = df_dw_norm,Ionome = Dat,
             paleta_stress = paleta_stress,
             paleta_type = paleta_type,
             paleta_enrichment = paleta_enrichment)
saveRDS(object = mDat,file = "../cleandata/dat_syncom_stresses.RDS")
rm(list=ls())
dev.off()
