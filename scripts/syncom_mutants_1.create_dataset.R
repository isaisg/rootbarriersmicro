library(ohchibi)

#Merge datasets into a single structure
set.seed(130816)
setwd('/home/isai/Documents/results/rootbarriersmicro/scripts/')

############# Suberin ####################
df_suberin <- read.table(file = "../rawdata/dat_syncom_mutants_suberin.csv",header = T,sep = "\t")
colnames(df_suberin)[7] <- "Batch"


toremove <- c("pCASP1::abi1","pCASP1::abi2","pCASP1::MYB41","pELTP::abi1")

df_suberin <- which(!(df_suberin$Genotype %in% toremove)) %>%
  df_suberin[.,] %>% droplevels
df_suberin$Genotype <- df_suberin$Genotype %>% gsub(pattern = " ",replacement = "")

#Calculane normalization factor  for the weight of the plants
#Probably worth having two normalization factors one for NB and one for SC
#Use mean to normalize inside each batch
df_suberin %>% subset(Bacteria == "NB" & Genotype == "Col-0") %>%
  chibi.boxplot(Map = .,x_val = "Batch", y_val = "Suberin",style = "open")



#df_sub <- df_suberin %>% subset(Bacteria == "NB") %>% droplevels
df_sub <- df_suberin
overall_me <- mean(subset(df_sub,Bacteria == "NB" &  Genotype == "Col-0")$Suberin)
df_sub_norm <- NULL
for(batch in unique(df_sub$Batch)){
  local_me <- subset(df_sub,Bacteria == "NB" & Batch == batch & Genotype == "Col-0") %$% Suberin %>% mean 
  nfactor <- overall_me/local_me
  temp <- subset(df_sub,Batch==batch)
  temp$NormSuberin <- temp$Suberin*nfactor
  df_sub_norm <- rbind(df_sub_norm,temp)
}


df_sub_norm %>% subset(Bacteria == "NB" & Genotype == "Col-0") %>% 
  chibi.boxplot(Map = .,x_val = "Batch",y_val = "NormSuberin",
                style = "open")  + 
  theme(axis.text.x = element_text(angle = 90))



df_suberin_norm <- df_sub_norm

###### Primary root elongation ##########
df_root <- read.table(file = "../rawdata/dat_syncom_mutants_primaryroot.csv",header = T,sep = "\t")

df_root$Genotype <- df_root$Genotype %>%
  gsub(pattern = "-",replacement = "\\.") %>%
  gsub(pattern = "esb.1",replacement = "esb1") %>%
gsub(pattern = "Col\\.0",replacement = "Col-0") %>%
  gsub(pattern = "abi1\\.1",replacement = "abi1") %>%
gsub(pattern = "\\/",replacement = "") %>%
  gsub(pattern = "sgn3\\.3",replacement = "sgn3") %>%
  gsub(pattern = "sgn3myb36 pELTP",replacement = "sgn3.3myb36 pELTP") %>%
  gsub(pattern = " ",replacement = "")

df_root <- which(!(df_root$Genotype %in% toremove)) %>% df_root[.,] %>% droplevels
df_root$Genotype %>% table
colnames(df_root)[8] <- "Batch"

#Calculane normalization factor  for the weight of the plants
#Probably worth having two normalization factors one for NB and one for SC
#Use mean to normalize inside each batch
df_root %>% subset(Bacteria == "NB" & Genotype == "Col-0") %>%
  chibi.boxplot(Map = .,x_val = "Batch", y_val = "RootElongation",style = "open")



#df_sub <- df_suberin %>% subset(Bacteria == "NB") %>% droplevels
df_sub <- df_root
overall_me <- mean(subset(df_sub,Bacteria == "NB" &  Genotype == "Col-0")$RootElongation)
df_root_norm <- NULL
for(batch in unique(df_sub$Batch)){
  local_me <- subset(df_sub,Bacteria == "NB" & Batch == batch & Genotype == "Col-0") %$% RootElongation %>% mean 
  nfactor <- overall_me/local_me
  temp <- subset(df_sub,Batch==batch)
  temp$NormRoot <- temp$RootElongation*nfactor
  df_root_norm <- rbind(df_root_norm,temp)
}


df_root_norm %>% subset(Bacteria == "NB" & Genotype == "Col-0") %>% 
  chibi.boxplot(Map = .,x_val = "Batch",y_val = "NormRoot",
                style = "open")  + 
  theme(axis.text.x = element_text(angle = 90))



df_root_norm <- df_root_norm


################### Dry weight #################
df <- read.table(file = "../rawdata/dat_syncom_mutants_weight.csv",header = T,sep = "\t")
df$Type <- "SC"
df$Type[df$Sample_name %>% grep(pattern = "No Bacteria|NB")] <- "NB"
df$Type[df$Sample_name %>% grep(pattern = "HK")] <- "HK"

df$Type <- df$Type %>% factor


#Calculane normalization factor  for the weight of the plants
#Probably worth having two normalization factors one for NB and one for SC
#Use mean to normalize inside each batch
df_suberin <- df
colnames(df_suberin)[10] <- "Bacteria"
colnames(df_suberin)[3] <- "Batch"
df_suberin$Batch <- paste0("Batch",df_suberin$Batch)

df_suberin <- df_suberin %>%subset(Bacteria != "HK") %>% droplevels

df_suberin$Genotype <- df_suberin$Sample_name %>% 
  gsub(pattern = " .*",replacement = "")  %>%
gsub(pattern = "abi2-1",replacement = "aba2.1") %>%
  gsub(pattern = "\\-",replacement = "\\.") %>%
  gsub(pattern = "Col\\.0",replacement = "Col-0") %>%
  gsub(pattern = "\\/",replacement = "")  %>%
  gsub(pattern = "pCASP1::CDEFesb1",replacement = "esb1.1pCASP1::CDEF") %>%
  gsub(pattern = "pCASP1::CDEFmyb36sgn3.3",replacement = "sgn3.3myb36pELTP::CEDF") %>%
  gsub(pattern = "^sgn3.3$",replacement = "sgn3") %>%
  gsub(pattern = "sgn3.3myb36",replacement = "sgn3myb36") %>%
  gsub(pattern = "sgn3myb36pELTP::CEDF",replacement = "sgn3.3myb36pELTP::CEDF")

toremove <- c(toremove,"pCASP1::abi1.1")
  

df_suberin <- which(!(df_suberin$Genotype %in% toremove)) %>%
  df_suberin[.,] %>% droplevels
df_suberin$Genotype %>% table


df_suberin %>% subset(Bacteria == "NB" & Genotype == "Col-0") %>%
  chibi.boxplot(Map = .,x_val = "Batch", y_val = "weight_mg_plant",style = "open")

#df_sub <- df_suberin %>% subset(Bacteria == "NB") %>% droplevels
df_sub <- df_suberin
overall_me <- mean(subset(df_sub,Bacteria == "NB" &  Genotype == "Col-0")$weight_mg_plant)
df_sub_norm <- NULL
for(batch in unique(df_sub$Batch)){
  local_me <- subset(df_sub,Bacteria == "NB" & Batch == batch & Genotype == "Col-0") %$% weight_mg_plant %>% mean 
  nfactor <- overall_me/local_me
  temp <- subset(df_sub,Batch==batch)
  temp$Normweight_mg_plant <- temp$weight_mg_plant*nfactor
  df_sub_norm <- rbind(df_sub_norm,temp)
}

df_sub_norm %>% subset(Bacteria == "NB" & Genotype == "Col-0") %>%
  chibi.boxplot(Map = .,x_val = "Batch", y_val = "Normweight_mg_plant",style = "open")


## models ###
df_dw_norm <- df_sub_norm


#Read the ionome matrix
df_ionome <- read.table("../rawdata/dat_syncom_mutants_ionome.csv",header = T,sep = "\t",
                        quote = "",comment.char = "")
Map <- df_ionome[,1:6]
Tab <- df_ionome[,7:ncol(df_ionome)]

Map$Sample.Id <- Map$Sample.Id %>% gsub(pattern = " ",replacement = "_")
colnames(Map)[2] <- "Sample_name"
colnames(Map)[3] <- "Batch"
Map$Batch <- paste0("Batch",Map$Batch)
Map$Plates <- paste0("Plate",Map$Plates)
Map$Replicates <- paste0("Rep",Map$Replicates )
Map$Run <- paste0("Run",Map$Run)
colnames(Map)[1] <- c("Sample_Id")
Map$Bacteria <- Map$Sample_name %>% gsub(pattern = ".* ",replacement = "") %>%
  factor
Map$Genotype <- Map$Sample_name %>% gsub(pattern = " .*",replacement = "") 

Map$Index <- rownames(Tab)


## Match genotypes to all the other datasets
Map$Genotype <- Map$Genotype %>%
  gsub(pattern = "sgn3-3$",replacement = "sgn3",perl = T) %>%
  gsub(pattern = "sgn3-3/myb36",replacement = "sgn3myb36") %>%
  gsub(pattern = "pCASP1::CDEF/esb1",replacement = "esb1.1pCASP1::CDEF") %>%
  gsub(pattern = "pCASP1::CDEF/myb36/sgn3",replacement = "sgn3.3myb36pELTP::CEDF") %>%
  gsub(pattern = "abi2-1",replacement = "aba2.1") %>%
  gsub(pattern = "abi4-1",replacement = "abi4.1") %>%
  gsub(pattern = "etr1-1",replacement = "etr1.1") %>%
  gsub(pattern = "ein3-1",replacement = "ein3.1") 


Map <- which(!(Map$Genotype %in% toremove)) %>%
  Map[.,] %>% droplevels


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

melted %>% subset( Ion == "Na23" & Genotype == "Col-0" & Bacteria == "NB") %>% droplevels %>%
  ggplot(.,aes(Batch,value)) + 
  facet_wrap(facets = "Ion") + geom_boxplot()  + geom_point()

#There exists variation associated with biological variation besides the measurement of the ionome
#We should adjust for this variation


#Do it across all Ions
melted_normalized <- NULL
for(ion in melted$Ion %>% unique){
  temp <- melted %>% subset(Ion == ion) %>% droplevels
  temp_nb <- temp %>% subset(Bacteria == "NB" & Genotype == "Col-0") %>% droplevels
  overall_me <- temp_nb$value %>% mean
  for(exp in temp$Batch %>% unique){
    local_me <- subset(temp, Batch == exp & Bacteria == "NB" & Genotype == "Col-0") %$% value %>% mean 
    nfactor <- overall_me/local_me
    ftemp <- temp %>% subset(Batch == exp) %>% droplevels
    ftemp$NormValue <- ftemp$value*nfactor
    melted_normalized <- rbind(melted_normalized,ftemp)
  }
}

melted_normalized %>% subset( Ion == "Na23" & Bacteria == "NB" & Genotype == "Col-0") %>% droplevels %>%
  ggplot(.,aes(Batch,NormValue)) + 
  facet_wrap(facets = "Ion") + geom_boxplot()  + geom_point()

Tab <- acast(data = melted_normalized,formula =Index ~ Ion,value.var = "NormValue" ,fill = 0)
Map <- match(rownames(Tab),rownames(Map)) %>%
  Map[.,]
Dat <- create_dataset(Tab = Tab %>% t,Map = Map)

paleta_enrichment <- c("#B3B3B3","#2480FF","#E65037")
names(paleta_enrichment) <- c("NoSignificant","Down","Up")

#Save the structure with all the datasets
df_suberin_norm$Genotype <- df_suberin_norm$Genotype %>%
  gsub(pattern = " ",replacement = "") %>%
  factor %>% relevel(ref = "Col-0")
df_root_norm$Genotype <- df_root_norm$Genotype %>%
  gsub(pattern = " ",replacement = "") %>%
  factor %>% relevel(ref = "Col-0")


  
df_dw_norm$Genotype <- df_dw_norm$Genotype %>%
  gsub(pattern = " ",replacement = "") %>%
  factor %>% relevel(ref = "Col-0")

Dat$Map$Genotype <- Dat$Map$Genotype %>%
  gsub(pattern = " ",replacement = "") %>%
  factor %>% relevel(ref = "Col-0")


Dat$Map$Bacteria <- Dat$Map$Bacteria %>% gsub(pattern = "SC",replacement = "SynCom") %>%
  factor

df_dw_norm$Bacteria <- df_dw_norm$Bacteria %>% gsub(pattern = "SC",replacement = "SynCom") %>%
  factor
df_suberin_norm$Bacteria <- df_suberin_norm$Bacteria%>% gsub(pattern = "SC",replacement = "SynCom") %>%
  factor

colnames(Dat$Map)[7] <- "Type"
colnames(df_dw_norm)[10] <- "Type"
colnames(df_suberin_norm)[3] <- "Type"
colnames(df_root_norm)[4] <- "Type"


df_root_norm$Type <- df_root_norm$Type %>%
  gsub(pattern = "SC",replacement = "SynCom") %>%
  factor


paleta_type <- c("#525252","#D9D9D9","#b48c36")
names(paleta_type) <- c("NB","HK","SynCom")



mDat <- list(Suberin = df_suberin_norm,
             Root = df_root_norm,
             DryWeight = df_dw_norm,
             Ionome = Dat,
             paleta_type = paleta_type,
             paleta_enrichment = paleta_enrichment)
saveRDS(object = mDat,file = "../cleandata/dat_syncom_genotypes.RDS")
rm(list=ls())
dev.off()


