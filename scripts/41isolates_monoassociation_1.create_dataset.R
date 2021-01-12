library(ohchibi)

#Merge datasets into a single structure
set.seed(130816)
setwd('/home/isai/Documents/results/rootbarriersmicro/scripts/')

################## DRY weight ##############

### Load the dryweight data as template ###
df_wei <- read.table(file = "../rawdata/data_41isolatesmonoassociation_dryweight.csv",header = T,sep = "\t")
colnames(df_wei)[2] <- "Strains"
colnames(df_wei)[3] <- "Batch"
colnames(df_wei)[9] <- "WeightbyPlant"
df_wei$Batch <- paste0("Batch",df_wei$Batch)
df_wei$Strains <- df_wei$Strains %>% gsub(pattern = "No Bacteria",replacement = "NB") %>%
  gsub(pattern = " ",replacement = "")

df_wei %>% subset(Strains == "NB") %>%
  chibi.boxplot(Map = .,x_val = "Batch",y_val = "WeightbyPlant",style = "open")

strains_wei <- df_wei$Strains%>% unique %>% as.character

#Use mean to normalize inside each batch
overall_me <- mean(subset(df_wei,Strains=="NB")$WeightbyPlant)
df_wei_norm <- NULL
for(batch in unique(df_wei$Batch)){
  local_me <- subset(df_wei, Batch == batch & Strains == "NB") %$% WeightbyPlant %>% mean 
  nfactor <- overall_me/local_me
  temp <- subset(df_wei,Batch==batch)
  temp$Norm_WeightbyPlant <- temp$WeightbyPlant*nfactor
  df_wei_norm <- rbind(df_wei_norm,temp)
}

df_wei_norm %>% subset(Strains == "NB") %>%
  chibi.boxplot(Map = .,x_val = "Batch",y_val = "Norm_WeightbyPlant",style = "open")



############# Suberin ####################
df_suberin <- read.table(file = "../rawdata/data_41isolatesmonoassociation_suberin.csv",header = T,sep = "\t")
colnames(df_suberin)[2] <- "Genotype"
colnames(df_suberin)[3] <- "Strains"
df_suberin$Media <- df_suberin$Media %>% 
  sub(pattern = " ",replacement = "") %>% factor


strains_suberin <- df_suberin$Strains %>% unique %>% as.character

#Split to add the synthetic batches
df_new <- NULL
counter <- 1
df_temp <- NULL
flag = 1
for(i in 1:nrow(df_suberin)){
  temp <- df_suberin[i,]
  if(temp[,1] == 1 & flag != 1){
    df_temp$Batch <- paste0("Batch",counter)
    df_new <- rbind(df_new,df_temp)
    counter <- counter +1
    df_temp <- NULL
  }else{
    flag <- 0
    df_temp <- rbind(df_temp,temp)
  }
}
df_temp$Batch <- paste0("Batch",counter)
df_new <- rbind(df_new,df_temp)

#Loop to add the batches as in the dryweight
df_suberin_withbatches <- NULL
for(batch in df_new$Batch %>% unique){
  mstrains <- df_new %>% subset(Batch == batch) %$% Strains %>% unique %>%
    grep(pattern = "NB|CDEF",value = T,invert = T)
  fbatches <- which(df_wei_norm$Strains %in% mstrains ) %>%
    df_wei_norm[.,] %$% Batch %>% unique %>%
    as.character %>%  sort
  st_dw <- df_wei_norm %>% subset(Batch == fbatches) %$% Strains %>% unique %>%
    grep(pattern = "NB|HK|CDEF",value = T,invert = T)
  int <- intersect(mstrains,st_dw) %>% length
  if((length(mstrains) == int) & (length(st_dw) == int)){
    fbatches <- fbatches %>% gsub(pattern = "Batch",replacement = "") %>%
      as.numeric %>% sort
    cat(fbatches,"\n")
    df_one <- df_new %>% subset(Batch == batch & Rep == 1) %>% droplevels
    df_one$RBatch <- rep(paste0("Batch",fbatches[1]),nrow(df_one))
    df_two <- df_new %>% subset(Batch == batch & Rep == 2) %>% droplevels
    df_two$RBatch <- rep(paste0("Batch",fbatches[2]),nrow(df_two))
    df_both <- rbind(df_one,df_two)
    df_suberin_withbatches <- rbind(df_suberin_withbatches,df_both)
  }else{
    cat("Failed\n")
  }
}

df_suberin <- df_suberin_withbatches[,-8]
colnames(df_suberin)[8] <- "Batch"
df_suberin$Rep <- paste0("Rep",df_suberin$Rep)



#Calculane normalization factor  for the weight of the plants
df_suberin %>% subset(Strains == "NB") %>% chibi.boxplot(Map = .,x_val = "Batch",
                                                         y_val = "Distance_tip_cz",style = "open")

#Use mean to normalize inside each batch
overall_me <- mean(subset(df_suberin,Strains=="NB")$Distance_tip_cz)
df_suberin_norm <- NULL
for(batch in unique(df_suberin$Batch)){
  local_me <- subset(df_suberin, Batch == batch & Strains == "NB") %$% Distance_tip_cz %>% mean 
  nfactor <- overall_me/local_me
  temp <- subset(df_suberin,Batch==batch)
  temp$Norm_Distance_tip_cz <- temp$Distance_tip_cz*nfactor
  df_suberin_norm <- rbind(df_suberin_norm,temp)
}


############# Primary root elongation ####################
df_mroot <- read.table(file = "../rawdata/data_41isolatesmonoassociation_primaryroot.csv",header = T,sep = "\t")
colnames(df_mroot)[3] <- "Strains"
df_mroot$Strains  <- df_mroot$Strains %>% 
  gsub(pattern = "^CL69$",replacement = "RCL69",perl = T) %>%
  gsub(pattern = "-",replacement = "")


strains_root <- df_mroot$Strains %>% as.character %>%
  gsub(pattern = "HK",replacement = "") %>% unique
which(!(strains_root %in% strains_wei)) %>%
  strains_root[.]
which(!(strains_root %in% strains_suberin)) %>%
  strains_root[.]

#RM32 is not measured in CFU and suberin

#Split to add the synthetic batches
df_new <- NULL
counter <- 1
df_temp <- NULL
flag = 1
for(i in 1:nrow(df_mroot)){
  temp <- df_mroot[i,]
  if(temp[,1] == 1 & flag != 1){
    df_temp$Batch <- paste0("Batch",counter)
    df_new <- rbind(df_new,df_temp)
    counter <- counter +1
    df_temp <- NULL
  }else{
    flag <- 0
    df_temp <- rbind(df_temp,temp)
  }
}
df_temp$Batch <- paste0("Batch",counter)
df_new <- rbind(df_new,df_temp)


#Loop to add the batches as in the dryweight
df_mroot_withbatches <- NULL
for(batch in df_new$Batch %>% unique){
  mstrains <- df_new %>% subset(Batch == batch) %$% Strains %>% unique %>%
    grep(pattern = "NB|CDEF",value = T,invert = T)
  fbatches <- which(df_wei_norm$Strains %in% mstrains ) %>%
    df_wei_norm[.,] %$% Batch %>% unique %>%
    as.character %>%  sort
  st_dw <- df_wei_norm %>% subset(Batch == fbatches) %$% Strains %>% unique %>%
    grep(pattern = "NB|CDEF",value = T,invert = T)
  int <- intersect(mstrains,st_dw) %>% length
  if((length(mstrains) == int) & (length(st_dw) == int)){
    fbatches <- fbatches %>% gsub(pattern = "Batch",replacement = "") %>%
      as.numeric %>% sort
    cat(fbatches,"\n")
    df_one <- df_new %>% subset(Batch == batch & Rep == 1) %>% droplevels
    df_one$RBatch <- rep(paste0("Batch",fbatches[1]),nrow(df_one))
    df_two <- df_new %>% subset(Batch == batch & Rep == 2) %>% droplevels
    df_two$RBatch <- rep(paste0("Batch",fbatches[2]),nrow(df_two))
    df_both <- rbind(df_one,df_two)
    df_mroot_withbatches <- rbind(df_mroot_withbatches,df_both)
  }else{
    cat(batch,"\tFailed\n")
    fbatches <- fbatches %>% gsub(pattern = "Batch",replacement = "") %>%
      as.numeric %>% sort
    cat(fbatches,"\n")
    df_one <- df_new %>% subset(Batch == batch & Rep == 1) %>% droplevels
    df_one$RBatch <- rep(paste0("Batch",fbatches[1]),nrow(df_one))
    df_two <- df_new %>% subset(Batch == batch & Rep == 2) %>% droplevels
    df_two$RBatch <- rep(paste0("Batch",fbatches[2]),nrow(df_two))
    df_both <- rbind(df_one,df_two)
    df_mroot_withbatches <- rbind(df_mroot_withbatches,df_both)
  }
}


df_mroot <- df_mroot_withbatches[,-8]
colnames(df_mroot)[8] <- "Batch"
df_mroot$Rep <- paste0("Rep",df_mroot$Rep)



#Use median to normalize inside each batch
overall_me <- mean(subset(df_mroot,Strains=="NB")$Main_root_elongation)
df_mroot_norm <- NULL
for(batch in unique(df_mroot$Batch)){
  local_me <- subset(df_mroot, Batch == batch & Strains == "NB") %$% Main_root_elongation %>% mean 
  nfactor <- overall_me/local_me
  temp <- subset(df_mroot,Batch==batch)
  temp$Norm_Main_root_elongation <- temp$Main_root_elongation*nfactor
  df_mroot_norm <- rbind(df_mroot_norm,temp)
}

############### IONOME ###############

#Read the ionome matrix
df_ionome <- read.table("../rawdata/data_41isolatesmonoassociation_ionome.csv",header = T,sep = "\t",quote = "",comment.char = "")
Map <- df_ionome[,1:6]
Tab <- df_ionome[,7:ncol(df_ionome)]
colnames(Map)[1] <- "Sample.Id"
colnames(Map)[2] <- "Sample_name"

Map$Sample.Id <- Map$Sample.Id %>% gsub(pattern = " ",replacement = "_")
Map$Sample_name <- Map$Sample_name %>% gsub(pattern = " ",replacement = "") %>%
  gsub(pattern = "NoBacteria",replacement = "NB")
colnames(Map)[3] <- "Batch"
Map$Batch <- paste0("Batch",Map$Batch)
Map$Plates <- paste0("Plate",Map$Plates)
Map$Replicates <- paste0("Rep",Map$Replicates )
Map$Run <- paste0("Run",Map$Run)
colnames(Map)[1:2] <- c("Sample_Id","Strains")
Map$Index <- rownames(Tab)


Tab <- match(Map$Index,rownames(Tab)) %>%
  Tab[.,]


rownames(Tab) <- Map$Index
colnames(Tab) <- colnames(Tab) %>% gsub(pattern = "\\.",replacement = "")
Tab <- as.matrix(Tab)



melted <- Map  %$% Index %>%
  match(rownames(Tab)) %>% Tab[.,]  %>% melt
colnames(melted) <- c("Index","Ion","value")
melted <- merge(melted,Map, by = "Index")
melted %>% subset( Ion == "Na23" & Strains == "NB") %>% droplevels %>%
  ggplot(.,aes(Batch,value)) + 
  facet_wrap(facets = "Ion") + geom_boxplot()  + geom_point()

#There exists variation associated with biological variation besides the measurement of the ionome
#We should adjust for this variation
#Do it across all Ions
melted_normalized <- NULL
for(ion in melted$Ion %>% unique){
  temp <- melted %>% subset(Ion == ion) %>% droplevels
  temp_nb <- temp %>% subset(Strains == "NB") %>% droplevels
  overall_me <- temp_nb$value %>% mean
  for(exp in temp$Batch %>% unique){
    local_me <- subset(temp, Batch == exp & Strains == "NB") %$% value %>% mean 
    nfactor <- overall_me/local_me
    ftemp <- temp %>% subset(Batch == exp) %>% droplevels
    ftemp$NormValue <- ftemp$value*nfactor
    melted_normalized <- rbind(melted_normalized,ftemp)
  }
}

melted_normalized %>% subset( Ion == "Na23" & Strains == "NB") %>% droplevels %>%
  ggplot(.,aes(Batch,NormValue)) + 
  facet_wrap(facets = "Ion") + geom_boxplot()  + geom_point()

melted_normalized$Type <- rep("HK",nrow(melted_normalized))
melted_normalized$Type[melted_normalized$Strains %>% grep(pattern = "HK",invert = T)] <- "Alive"
melted_normalized$Type[which(melted_normalized$Strains == "NB")] <- "NB"
melted_normalized$Type[which(melted_normalized$Strains == "CDEF")] <- "CDEF"

melted_normalized$Type <- melted_normalized$Type %>%
  factor(levels = c("NB","CDEF","HK","Alive"))

##### CFU  data ####
############# CFU ####################
df_cfu <- read.table(file = "../rawdata/data_41isolatesmonoassociation_cfu.csv",header = T,sep = "\t")
colnames(df_cfu)[2] <- "Strains"
df_cfu$cfu[which(is.na(df_cfu$cfu))] <- 0
colnames(df_cfu)[3] <- "sample"

#Split to add the synthetic batches
df_new <- NULL
counter <- 1
df_temp <- NULL
flag = 1
for(i in 1:nrow(df_cfu)){
  temp <- df_cfu[i,]
  if(temp[,1] == 1 & flag != 1){
    df_temp$Batch <- paste0("Batch",counter)
    df_new <- rbind(df_new,df_temp)
    counter <- counter +1
    df_temp <- NULL
  }else{
    flag <- 0
    df_temp <- rbind(df_temp,temp)
  }
}
df_temp$Batch <- paste0("Batch",counter)
df_new <- rbind(df_new,df_temp)


#Loop to add the batches as in the dryweight
df_cfu_withbatches <- NULL
for(batch in df_new$Batch %>% unique){
  mstrains <- df_new %>% subset(Batch == batch) %$% Strains %>% unique %>%
    grep(pattern = "NB|CDEF",value = T,invert = T)
  fbatches <- which(df_wei_norm$Strains %in% mstrains ) %>%
    df_wei_norm[.,] %$% Batch %>% unique %>%
    as.character %>%  sort
  st_dw <- df_wei_norm %>% subset(Batch == fbatches) %$% Strains %>% unique %>%
    grep(pattern = "NB|HK|CDEF",value = T,invert = T)
  int <- intersect(mstrains,st_dw) %>% length
  if((length(mstrains) == int) & (length(st_dw) == int)){
    fbatches <- fbatches %>% gsub(pattern = "Batch",replacement = "") %>%
      as.numeric %>% sort
    cat(fbatches,"\n")
    df_one <- df_new %>% subset(Batch == batch & Rep == "Rep1") %>% droplevels
    df_one$RBatch <- rep(paste0("Batch",fbatches[1]),nrow(df_one))
    df_two <- df_new %>% subset(Batch == batch & Rep == "Rep2") %>% droplevels
    df_two$RBatch <- rep(paste0("Batch",fbatches[2]),nrow(df_two))
    df_both <- rbind(df_one,df_two)
    df_cfu_withbatches <- rbind(df_cfu_withbatches,df_both)
  }else{
    cat(batch,"\tFailed\n")
  }
}

df_cfu <- df_cfu_withbatches[,-12]
colnames(df_cfu)[12] <- "Batch"



#Use median to normalize inside fractions
df_cfu_norm <- NULL
for(fraction in unique(df_cfu$sample)){
  temp <- df_cfu %>% subset(sample == fraction) %>% droplevels
  overall_me <- mean(subset(temp,Strains=="NB")$mg.FW)
  #Normalize batch 
  for(batch in unique(temp$Batch)){
    local_me <- subset(temp, Batch == batch & Strains == "NB") %$% mg.FW %>% mean 
    nfactor <- overall_me/local_me
    intemp <- subset(temp,Batch==batch)
    intemp$Norm_mg.FW <- intemp$mg.FW*nfactor
    df_cfu_norm <- rbind(df_cfu_norm,intemp)
  }
  
}


#Change the 
df_cfu_norm$Norm_cfu_weight <- df_cfu_norm$cfu / df_cfu_norm$Norm_mg.FW

#Add pseudocount
df_cfu_norm$Log_Norm_cfu_weight <- log10(df_cfu_norm$Norm_cfu_weight+1)
df_cfu_norm$Log_Norm_cfu_weight[which(is.infinite(x = df_cfu_norm$Log_Norm_cfu_weight))] <- 0
df_cfu_norm$Log_Norm_cfu_weight[which(is.na(df_cfu_norm$Log_Norm_cfu_weight))] <- 0
colnames(df_cfu_norm)[3] <- "Fraction"
df_cfu_norm$Fraction <- df_cfu_norm$Fraction %>% gsub(pattern = "shoot",replacement = "Shoot") %>%
  gsub(pattern = "root",replacement = "Root") %>%
  gsub(pattern = "agar",replacement = "Agar") %>%
  factor(levels = c("Agar","Root","Shoot"))




#Add unique sample ids
df_cfu_norm$Sample_Id <- paste0("Sample",1:nrow(df_cfu_norm))
df_suberin_norm$Sample_Id <- paste0("Sample",1:nrow(df_suberin_norm))
df_mroot_norm$Sample_Id <- paste0("Sample",1:nrow(df_mroot_norm))
df_wei_norm$Sample_Id <- paste0("Sample",1:nrow(df_wei_norm))

#Clean spaces from the strains variable
df_cfu_norm$Strains <- df_cfu_norm$Strains %>% gsub(pattern = " ",replacement = "")
df_suberin_norm$Strains <- df_suberin_norm$Strains %>% gsub(pattern = " ",replacement = "")
df_mroot_norm$Strains <- df_mroot_norm$Strains %>% gsub(pattern = " ",replacement = "")
melted_normalized$Strains <- melted_normalized$Strains %>% gsub(pattern = " ",replacement = "")

Dat_init <- readRDS("../cleandata/dat_screeningbarriers.RDS")
Map_original <- Dat_init$Map

Map$Strains_original <- Map$Strains %>% gsub(pattern = "HK",replacement = "")
Map_original <- Map_original[,c(16,1,3:7)]
colnames(Map_original)[1] <- "Strains_original"

Map <- merge(Map,Map_original,by = "Strains_original",all.x = TRUE) 
Map <- Map[,-c(2,4,5,6,7,8)] %>% unique


#Add the normalize suberin by primary root length
df_cf_root <- aggregate(Norm_Main_root_elongation ~ Strains+Batch,df_mroot_norm,FUN = "mean")
df_cf_root$UId <- paste0(df_cf_root$Strains,"_",df_cf_root$Batch)

#Merge in the context of suberin
df_suberin_norm$UId <- paste0(df_suberin_norm$Strains,"_",df_suberin_norm$Batch)


nroot <- match(df_suberin_norm$UId,df_cf_root$UId) %>%
  df_cf_root$Norm_Main_root_elongation[.]


df_suberin_norm$Norm_Main_root_elongation <- nroot

#For two of them we only measured main root in one rep so we will use that as proxy
#for the second rep
which(is.na(df_suberin_norm$Norm_Main_root_elongation))
df_suberin_norm[which(is.na(df_suberin_norm$Norm_Main_root_elongation)),]

toaddmissing <- df_suberin_norm %>% subset(Strains == "CDEF") %$%
  Norm_Main_root_elongation %>% na.omit %>% mean

indices <- which(is.na(df_suberin_norm$Norm_Main_root_elongation))
df_suberin_norm$Norm_Main_root_elongation[indices[1:8]] <- toaddmissing

df_suberin_norm %>% subset(Strains == "RMF329")

df_suberin_norm$Norm_Main_root_elongation[indices[9:12]] <- 5.230773

#Do the normalization
df_suberin_norm$NormProp <- df_suberin_norm$Norm_Distance_tip_cz / df_suberin_norm$Norm_Main_root_elongation

#### Now add the raw primary root elongation just in case

df_cf_root <- aggregate(Main_root_elongation ~ Strains+Batch,df_mroot_norm,FUN = "mean")
df_cf_root$UId <- paste0(df_cf_root$Strains,"_",df_cf_root$Batch)


nroot <- match(df_suberin_norm$UId,df_cf_root$UId) %>%
  df_cf_root$Main_root_elongation[.]


df_suberin_norm$Main_root_elongation <- nroot

#For two of them we only measured main root in one rep so we will use that as proxy
#for the second rep
which(is.na(df_suberin_norm$Main_root_elongation))
df_suberin_norm[which(is.na(df_suberin_norm$Main_root_elongation)),]

toaddmissing <- df_suberin_norm %>% subset(Strains == "CDEF") %$%
  Main_root_elongation %>% na.omit %>% mean

indices <- which(is.na(df_suberin_norm$Main_root_elongation))
df_suberin_norm$Main_root_elongation[indices[1:8]] <- toaddmissing

df_suberin_norm %>% subset(Strains == "RMF329")

df_suberin_norm$Main_root_elongation[indices[9:12]] <- 5.721091


df_suberin_norm$Prop <- df_suberin_norm$Distance_tip_cz / df_suberin_norm$Main_root_elongation

#Create raw dataset
Raw <- list(CFU = df_cfu_norm, Suberin = df_suberin_norm, 
            PrimaryRoot = df_mroot_norm,DryWeight = df_wei_norm, 
            Ionome = melted_normalized,Map = Map
)


############################
saveRDS(object = Raw,
        file = "../cleandata/dat_41isolatesmonoassociation_cfu_suberin_primroot_dryweight_ionome.RDS")
#Clean workspace
dev.off()
gc()
rm(list=ls())

