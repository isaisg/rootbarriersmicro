library(ohchibi)
library(plyr)

#Merge datasets into a single structure
set.seed(130816)
setwd('/home/isai/Documents/results/rootbarriersmicro/scripts/')


Dat <- readRDS("../cleandata/dat_41isolatesmonoassociation_cfu_suberin_primroot_dryweight_ionome.updated.RDS")


################# Raw values ########################
##### Ionome #######
df_ionome <- Dat$Ionome

#Aggregate the dataframe to display as heatmap
melted <- df_ionome %>% acast(formula = Index~Ion,fill = 0,
                              value.var = "NormValue") 


melted_quant <- matrix(data = NA,nrow = nrow(melted),ncol = ncol(melted))
for(i in 1:ncol(melted)){
  x <- melted[,i]
  maxvalue <- quantile(x,0.99)
  minvalue <- quantile(x,0.01)
  x[which(x > maxvalue)] <- maxvalue
  x[which((x < minvalue))] <- minvalue
  melted_quant[,i] <- x
}
colnames(melted_quant) <- colnames(melted)
melted <- melted_quant 
melted <- melted %>% melt

Map <- df_ionome[,c(1,5,11)] %>% unique %>% droplevels
colnames(melted) <- c("Index","Ion","zscore")
melted <- merge(melted,Map, by = "Index")                    
melted_ionome <- melted
sum_ionome <- ddply(melted, c("Strains","Ion"), summarise,
                    N    = length(zscore),
                    mean = mean(zscore),
                    sd   = sd(zscore),
                    se   = sd / sqrt(N),
                    ci = (qt(0.95/2 + 0.5, N- 1))*se
)

sum_ionome_mean <- dcast(data = sum_ionome,Strains~Ion,value.var = "mean")
sum_ionome_se <- dcast(data = sum_ionome,Strains~Ion,value.var = "ci")
colnames(sum_ionome_mean)[2:ncol(sum_ionome_mean)] <- paste0("Mean",colnames(sum_ionome_mean)[2:ncol(sum_ionome_mean)])
colnames(sum_ionome_se)[2:ncol(sum_ionome_se)] <- paste0("CI",colnames(sum_ionome_se)[2:ncol(sum_ionome_se)])

sum_ionome <- merge(sum_ionome_mean,sum_ionome_se,by ="Strains")



Tab_ionome <- acast(data = melted,formula = Strains~Ion,fun.aggregate = mean,
                    fill = 0,value.var = "zscore")
Tab_ionome <- Tab_ionome %>% as.data.frame
Tab_ionome$Strains <- rownames(Tab_ionome)

rownames(Tab_ionome) <- NULL
Tab_ionome <- Tab_ionome[,c(18,1:17)]

##### CFU #######
df_cfu <- Dat$CFU
df_temp <- df_cfu %>% subset(Fraction == "Agar") %>% droplevels
Tab_cfu_agar <- acast(df_temp,formula = Sample_Id~Fraction,value.var = "Log_Norm_cfu_weight",fill = 0)
df_temp <- df_cfu %>% subset(Fraction == "Root") %>% droplevels
Tab_cfu_root <- acast(df_temp,formula = Sample_Id~Fraction,value.var = "Log_Norm_cfu_weight",fill = 0)
df_temp <- df_cfu %>% subset(Fraction == "Shoot") %>% droplevels
Tab_cfu_shoot <- acast(df_temp,formula = Sample_Id~Fraction,value.var = "Log_Norm_cfu_weight",fill = 0)

#Quantile normalizations
maxvalue <- Tab_cfu_agar %>% melt %$% value  %>% quantile(0.99)
minvalue <- Tab_cfu_agar %>% melt %$% value  %>% quantile(0.01)
Tab_cfu_agar[Tab_cfu_agar > maxvalue]  <- maxvalue
Tab_cfu_agar[Tab_cfu_agar < minvalue]  <- minvalue

maxvalue <- Tab_cfu_root %>% melt %$% value  %>% quantile(0.99)
minvalue <- Tab_cfu_root %>% melt %$% value  %>% quantile(0.01)
Tab_cfu_root[Tab_cfu_root > maxvalue]  <- maxvalue
Tab_cfu_root[Tab_cfu_root < minvalue]  <- minvalue

maxvalue <- Tab_cfu_shoot %>% melt %$% value  %>% quantile(0.99)
minvalue <- Tab_cfu_shoot %>% melt %$% value  %>% quantile(0.01)
Tab_cfu_shoot[Tab_cfu_shoot > maxvalue]  <- maxvalue
Tab_cfu_shoot[Tab_cfu_shoot < minvalue]  <- minvalue



Tab_cfu <- rbind(Tab_cfu_agar %>% melt,Tab_cfu_root %>% melt) %>%
  rbind(Tab_cfu_shoot %>% melt)

colnames(Tab_cfu) <- c("Sample_Id","Fraction","zscore")
Tab_cfu <- merge(Tab_cfu,df_cfu[,-3], by = "Sample_Id")


sum_cfu <- ddply(Tab_cfu, c("Strains","Fraction"), summarise,
                 N    = length(zscore),
                 mean = mean(zscore),
                 sd   = sd(zscore),
                 se   = sd / sqrt(N),
                 ci = (qt(0.95/2 + 0.5, N- 1))*se
                 
)
sum_cfu_mean <- dcast(data = sum_cfu,Strains~Fraction,value.var = "mean")
sum_cfu_se <- dcast(data = sum_cfu,Strains~Fraction,value.var = "ci")
sum_cfu <- merge(sum_cfu_mean,sum_cfu_se,by ="Strains")


##### Suberin #####
df_suberin <- Dat$Suberin
#df_suberin$zscore <- df_suberin$Norm_Distance_tip_cz
df_suberin$zscore <- df_suberin$Norm_Distance_tip_cz



maxvalue <- df_suberin$zscore %>% quantile(0.99)
minvalue <- df_suberin$zscore %>% quantile(0.01)
df_suberin$zscore[which(df_suberin$zscore > maxvalue)] <- maxvalue
df_suberin$zscore[which(df_suberin$zscore < minvalue)] <- minvalue
sum_suberin <- ddply(df_suberin, c("Strains"), summarise,
                     N    = length(zscore),
                     mean = mean(zscore),
                     sd   = sd(zscore),
                     se   = sd / sqrt(N),
                     ci = (qt(0.95/2 + 0.5, N- 1))*se
                     
)


### Primary Root ####
df_mroot <- Dat$PrimaryRoot
df_mroot$zscore <- df_mroot$Norm_Main_root_elongation 

maxvalue <- df_mroot$zscore %>% quantile(0.99)
minvalue <- df_mroot$zscore %>% quantile(0.01)
df_mroot$zscore[which(df_mroot$zscore > maxvalue)] <- maxvalue
df_mroot$zscore[which(df_mroot$zscore < minvalue)] <- minvalue
sum_root <- ddply(df_mroot, c("Strains"), summarise,
                  N    = length(zscore),
                  mean = mean(zscore),
                  sd   = sd(zscore),
                  se   = sd / sqrt(N),
                  ci = (qt(0.95/2 + 0.5, N- 1))*se
                  
)


### Dry weight
df_weight <- Dat$DryWeight
df_weight$zscore <- df_weight$Norm_WeightbyPlant 

maxvalue <- df_weight$zscore %>% quantile(0.99)
minvalue <- df_weight$zscore %>% quantile(0.01)
df_weight$zscore[which(df_weight$zscore > maxvalue)] <- maxvalue
df_weight$zscore[which(df_weight$zscore < minvalue)] <- minvalue
sum_weight <- ddply(df_weight, c("Strains"), summarise,
                    N    = length(zscore),
                    mean = mean(zscore),
                    sd   = sd(zscore),
                    se   = sd / sqrt(N),
                    ci = (qt(0.95/2 + 0.5, N- 1))*se
                    
)



#Load the original Pi data
Dat <- readRDS(file = "../cleandata/dat_screeningbarriers.RDS")
df_prop <- Dat$Propidium
df_prop$zscore <- df_prop$Normvalue_ra 
maxvalue <- df_prop$zscore %>% quantile(0.99)
minvalue <- df_prop$zscore %>% quantile(0.01)
df_prop$zscore[which(df_prop$zscore > maxvalue)] <- maxvalue
df_prop$zscore[which(df_prop$zscore < minvalue)] <- minvalue

sum_prop <- ddply(df_prop, c("Strains"), summarise,
                  N    = length(zscore),
                  mean = mean(zscore),
                  sd   = sd(zscore),
                  se   = sd / sqrt(N),
                  ci = (qt(0.95/2 + 0.5, N- 1))*se
                  
)



#Original Suberin
df_sub <- Dat$Suberin
melted <- acast(data = df_sub,formula = UId ~Localization,value.var = "Normvalue_ra") 


suma <- melted[,1] + melted[,2]

melted <- cbind(melted,suma)
colnames(melted)[4] <- "SumNoDiscrete"

melted <- melted %>% melt
colnames(melted) <- c("UId","Zone","zscore")
melted$Strains <- match(melted$UId,df_sub$UId) %>% df_sub$Strains[.] %>% as.character
sum_sub_original <- ddply(melted, c("Strains","Zone"), summarise,
                          N    = length(zscore),
                          mean = mean(zscore),
                          sd   = sd(zscore),
                          se   = sd / sqrt(N),
                          ci = (qt(0.95/2 + 0.5, N- 1))*se
                          
)
sum_sub_original_mean <- dcast(data = sum_sub_original,Strains~Zone,value.var = "mean")
sum_sub_original_se <- dcast(data = sum_sub_original,Strains~Zone,value.var = "se")
sum_sub_original <- merge(sum_sub_original_mean,sum_sub_original_se,by ="Strains")



#Create the structure with the summarized data
merged <- merge(sum_prop, sum_sub_original, by = "Strains", all = TRUE) %>%
  merge(sum_suberin, by = "Strains", all = TRUE) %>%
  merge(sum_cfu, by = "Strains", all = TRUE) %>%
  merge(sum_root, by = "Strains", all = TRUE) %>%
  merge(sum_weight, by = "Strains", all = TRUE)



merged <- merged[,c(1,3,6:14,16,19,20:25,27,30,32,35)]

colnames(merged) <- c("Strains","MeanPropidium","CIPropidium",
                      "MeanNoExpSub","MeanDiscExpSub","MeanContExpSub","MeanSumNoDiscExpSub",
                      "CINoExpSub","CIDiscExpSub","CIContExpSub","CISumNoDiscExpSub",
                      "MeanDistSub","CIDistSub",
                      "MeanCFUAgar","MeanCFURoot","MeanCFUShoot",
                      "CICFUAgar","CICFURoot","CIFUShoot",
                      "MeanPrimRoot","CIPrimRoot",
                      "MeanDryWeight","CIDryWeight")

merged <- merge(merged,sum_ionome, by = "Strains", all = TRUE)
colnames(merged)

#Reorder the columns so they follow the mean +se structure
merged <- merged[,c(1,2,3,4,8,5,9,6,10,7,11,12,13,14,17,15,18,16,19,20:23,
                    24,41,25,42,26,43,27,44,28,45,29,46,30,47,31,48,32,49,33,50,34,51,
                    35,52,36,53,37,54,38,55,39,56,40,57)]

merged_raw <- merged


#Save raw structures
RawData <- list(Ionome = melted_ionome,CFU = Tab_cfu,Suberin = df_suberin,
                PrimaryRoot = df_mroot,DryWeight = df_weight,
                Propidium = df_prop,SuberinScr = melted)


################# Zscores ########################
##### Ionome #######
Dat <- readRDS("../cleandata/dat_41isolatesmonoassociation_cfu_suberin_primroot_dryweight_ionome.RDS")
df_ionome <- Dat$Ionome

#Aggregate the dataframe to display as heatmap
melted <- df_ionome %>% acast(formula = Index~Ion,fill = 0,
                              value.var = "NormValue")  %>% scale

melted_quant <- matrix(data = NA,nrow = nrow(melted),ncol = ncol(melted))
for(i in 1:ncol(melted)){
  x <- melted[,i]
  maxvalue <- quantile(x,0.99)
  minvalue <- quantile(x,0.01)
  x[which(x > maxvalue)] <- maxvalue
  x[which((x < minvalue))] <- minvalue
  melted_quant[,i] <- x
}
colnames(melted_quant) <- colnames(melted)
melted <- melted_quant 
melted <- melted %>% melt

Map <- df_ionome[,c(1,5,11)] %>% unique %>% droplevels
colnames(melted) <- c("Index","Ion","zscore")
melted <- merge(melted,Map, by = "Index")                    
melted_ionome <- melted
sum_ionome <- ddply(melted, c("Strains","Ion"), summarise,
                    N    = length(zscore),
                    mean = mean(zscore),
                    sd   = sd(zscore),
                    se   = sd / sqrt(N),
                    ci = (qt(0.95/2 + 0.5, N- 1))*se
                    
)

sum_ionome_mean <- dcast(data = sum_ionome,Strains~Ion,value.var = "mean")
sum_ionome_se <- dcast(data = sum_ionome,Strains~Ion,value.var = "ci")
colnames(sum_ionome_mean)[2:ncol(sum_ionome_mean)] <- paste0("Mean",colnames(sum_ionome_mean)[2:ncol(sum_ionome_mean)])
colnames(sum_ionome_se)[2:ncol(sum_ionome_se)] <- paste0("CI",colnames(sum_ionome_se)[2:ncol(sum_ionome_se)])

sum_ionome <- merge(sum_ionome_mean,sum_ionome_se,by ="Strains")



Tab_ionome <- acast(data = melted,formula = Strains~Ion,fun.aggregate = mean,
                    fill = 0,value.var = "zscore")
Tab_ionome <- Tab_ionome %>% as.data.frame
Tab_ionome$Strains <- rownames(Tab_ionome)

rownames(Tab_ionome) <- NULL
Tab_ionome <- Tab_ionome[,c(18,1:17)]

##### CFU #######
df_cfu <- Dat$CFU
df_temp <- df_cfu %>% subset(Fraction == "Agar") %>% droplevels
Tab_cfu_agar <- acast(df_temp,formula = Sample_Id~Fraction,value.var = "Log_Norm_cfu_weight",fill = 0) %>% scale
df_temp <- df_cfu %>% subset(Fraction == "Root") %>% droplevels
Tab_cfu_root <- acast(df_temp,formula = Sample_Id~Fraction,value.var = "Log_Norm_cfu_weight",fill = 0) %>% scale
df_temp <- df_cfu %>% subset(Fraction == "Shoot") %>% droplevels
Tab_cfu_shoot <- acast(df_temp,formula = Sample_Id~Fraction,value.var = "Log_Norm_cfu_weight",fill = 0) %>% scale

#Quantile normalizations
maxvalue <- Tab_cfu_agar %>% melt %$% value  %>% quantile(0.99)
minvalue <- Tab_cfu_agar %>% melt %$% value  %>% quantile(0.01)
Tab_cfu_agar[Tab_cfu_agar > maxvalue]  <- maxvalue
Tab_cfu_agar[Tab_cfu_agar < minvalue]  <- minvalue

maxvalue <- Tab_cfu_root %>% melt %$% value  %>% quantile(0.99)
minvalue <- Tab_cfu_root %>% melt %$% value  %>% quantile(0.01)
Tab_cfu_root[Tab_cfu_root > maxvalue]  <- maxvalue
Tab_cfu_root[Tab_cfu_root < minvalue]  <- minvalue

maxvalue <- Tab_cfu_shoot %>% melt %$% value  %>% quantile(0.99)
minvalue <- Tab_cfu_shoot %>% melt %$% value  %>% quantile(0.01)
Tab_cfu_shoot[Tab_cfu_shoot > maxvalue]  <- maxvalue
Tab_cfu_shoot[Tab_cfu_shoot < minvalue]  <- minvalue



Tab_cfu <- rbind(Tab_cfu_agar %>% melt,Tab_cfu_root %>% melt) %>%
  rbind(Tab_cfu_shoot %>% melt)

colnames(Tab_cfu) <- c("Sample_Id","Fraction","zscore")
Tab_cfu <- merge(Tab_cfu,df_cfu[,-3], by = "Sample_Id")


sum_cfu <- ddply(Tab_cfu, c("Strains","Fraction"), summarise,
                 N    = length(zscore),
                 mean = mean(zscore),
                 sd   = sd(zscore),
                 se   = sd / sqrt(N),
                 ci = (qt(0.95/2 + 0.5, N- 1))*se
                 
)
sum_cfu_mean <- dcast(data = sum_cfu,Strains~Fraction,value.var = "mean")
sum_cfu_se <- dcast(data = sum_cfu,Strains~Fraction,value.var = "ci")
sum_cfu <- merge(sum_cfu_mean,sum_cfu_se,by ="Strains")


##### Suberin #####
df_suberin <- Dat$Suberin
df_suberin$zscore <- df_suberin$Norm_Distance_tip_cz  %>% scale

maxvalue <- df_suberin$zscore %>% quantile(0.99)
minvalue <- df_suberin$zscore %>% quantile(0.01)
df_suberin$zscore[which(df_suberin$zscore > maxvalue)] <- maxvalue
df_suberin$zscore[which(df_suberin$zscore < minvalue)] <- minvalue
sum_suberin <- ddply(df_suberin, c("Strains"), summarise,
                     N    = length(zscore),
                     mean = mean(zscore),
                     sd   = sd(zscore),
                     se   = sd / sqrt(N),
                     ci = (qt(0.95/2 + 0.5, N- 1))*se
                     
)


### Primary Root ####
df_mroot <- Dat$PrimaryRoot
df_mroot$zscore <- df_mroot$Norm_Main_root_elongation  %>% scale

maxvalue <- df_mroot$zscore %>% quantile(0.99)
minvalue <- df_mroot$zscore %>% quantile(0.01)
df_mroot$zscore[which(df_mroot$zscore > maxvalue)] <- maxvalue
df_mroot$zscore[which(df_mroot$zscore < minvalue)] <- minvalue
sum_root <- ddply(df_mroot, c("Strains"), summarise,
                  N    = length(zscore),
                  mean = mean(zscore),
                  sd   = sd(zscore),
                  se   = sd / sqrt(N),
                  ci = (qt(0.95/2 + 0.5, N- 1))*se
                  
)


### Dry weight
df_weight <- Dat$DryWeight
df_weight$zscore <- df_weight$Norm_WeightbyPlant  %>% scale

maxvalue <- df_weight$zscore %>% quantile(0.99)
minvalue <- df_weight$zscore %>% quantile(0.01)
df_weight$zscore[which(df_weight$zscore > maxvalue)] <- maxvalue
df_weight$zscore[which(df_weight$zscore < minvalue)] <- minvalue
sum_weight <- ddply(df_weight, c("Strains"), summarise,
                    N    = length(zscore),
                    mean = mean(zscore),
                    sd   = sd(zscore),
                    se   = sd / sqrt(N),
                    ci = (qt(0.95/2 + 0.5, N- 1))*se
                    
)



#Load the original Pi data
Dat <- readRDS(file = "../cleandata/dat_screeningbarriers.RDS")
df_prop <- Dat$Propidium
df_prop$zscore <- df_prop$Normvalue_ra  %>% scale
maxvalue <- df_prop$zscore %>% quantile(0.99)
minvalue <- df_prop$zscore %>% quantile(0.01)
df_prop$zscore[which(df_prop$zscore > maxvalue)] <- maxvalue
df_prop$zscore[which(df_prop$zscore < minvalue)] <- minvalue

sum_prop <- ddply(df_prop, c("Strains"), summarise,
                  N    = length(zscore),
                  mean = mean(zscore),
                  sd   = sd(zscore),
                  se   = sd / sqrt(N),
                  ci = (qt(0.95/2 + 0.5, N- 1))*se
                  
)



#Original Suberin
df_sub <- Dat$Suberin
melted <- acast(data = df_sub,formula = UId ~Localization,value.var = "Normvalue_ra")
suma <- melted[,1] + melted[,2]

melted <- cbind(melted,suma)
colnames(melted)[4] <- "SumNoDiscrete"
melted <- melted %>% scale

melted <- melted %>% melt
colnames(melted) <- c("Merged","Zone","zscore")
melted$Strains <- match(melted$Merged,df_sub$Merged) %>% df_sub$Strains[.] %>% as.character
sum_sub_original <- ddply(melted, c("Strains","Zone"), summarise,
                          N    = length(zscore),
                          mean = mean(zscore),
                          sd   = sd(zscore),
                          se   = sd / sqrt(N),
                          ci = (qt(0.95/2 + 0.5, N- 1))*se
                          
)
sum_sub_original_mean <- dcast(data = sum_sub_original,Strains~Zone,value.var = "mean")
sum_sub_original_se <- dcast(data = sum_sub_original,Strains~Zone,value.var = "ci")
sum_sub_original <- merge(sum_sub_original_mean,sum_sub_original_se,by ="Strains")



#Create the structure with the summarized data
merged <- merge(sum_prop, sum_sub_original, by = "Strains", all = TRUE) %>%
  merge(sum_suberin, by = "Strains", all = TRUE) %>%
  merge(sum_cfu, by = "Strains", all = TRUE) %>%
  merge(sum_root, by = "Strains", all = TRUE) %>%
  merge(sum_weight, by = "Strains", all = TRUE)
merged <- merged[,c(1,3,6:14,16,19,20:25,27,30,32,35)]

colnames(merged) <- c("Strains","MeanPropidium","CIPropidium",
                      "MeanNoExpSub","MeanDiscExpSub","MeanContExpSub","MeanSumNoDiscExpSub",
                      "CINoExpSub","CIDiscExpSub","CIContExpSub","CISumNoDiscExpSub",
                      "MeanDistSub","CIDistSub",
                      "MeanCFUAgar","MeanCFURoot","MeanCFUShoot",
                      "CICFUAgar","CICFURoot","CIFUShoot",
                      "MeanPrimRoot","CIPrimRoot",
                      "MeanDryWeight","CIDryWeight")

merged <- merge(merged,sum_ionome, by = "Strains", all = TRUE)
colnames(merged)
#Reorder the columns so they follow the mean +se structure
merged <- merged[,c(1,2,3,4,8,5,9,6,10,7,11,12,13,14,17,15,18,16,19,20:23,
                    24,41,25,42,26,43,27,44,28,45,29,46,30,47,31,48,32,49,33,50,34,51,
                    35,52,36,53,37,54,38,55,39,56,40,57)]

merged_zscore <- merged
#Save raw structures
RawDataZscore <- list(Ionome = melted_ionome,CFU = Tab_cfu,Suberin = df_suberin,
                      PrimaryRoot = df_mroot,DryWeight = df_weight,
                      Propidium = df_prop,SuberinScr = melted)




melted <-  list(Raw = merged_raw, zscore = merged_zscore)
raw <- list(Raw = RawData, zscore = RawDataZscore)


all <- list(Merged = melted,Raw = raw)


saveRDS(object = all,
        file = "../cleandata/dat_41isolatesmonoassociation_cfu_suberin_primroot_dryweight_ionome_zscoreRaw.distance.updated.RDS")


rm(list=ls())
gc()










