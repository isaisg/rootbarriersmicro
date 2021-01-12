library(ohchibi)
library(egg)
library(emmeans)

set.seed(130816)
setwd('/home/isai/Documents/results/rootbarriersmicro/scripts/')



mDat <- readRDS(file = "../cleandata/dat_syncom_genotypes.RDS") 


###### Suberin ##########

df_sub <- mDat$Suberin

genotypes <- df_sub$Genotype %>% unique
Res <- NULL
for(geno in genotypes){
  df_temp <- df_sub %>% subset(Genotype == geno) %>% droplevels
  nb <- df_temp %>% subset(Type == "NB") %$% NormSuberin
  sc <- df_temp %>% subset(Type == "SynCom") %$% NormSuberin
  m1 <- t.test(nb,sc)
  Res <- data.frame(Genotype = geno,p.value = m1$p.value) %>%
    rbind(Res,.)
}

Res$p.adj <- Res$p.value %>% p.adjust(method = "fdr")
Res$Letter <- format.pval(Res$p.value, digits = 8)

Res_suberin <- Res
  
#### Root#######
df_sub <- mDat$Root


genotypes <- df_sub$Genotype %>% unique
Res <- NULL
for(geno in genotypes){
  df_temp <- df_sub %>% subset(Genotype == geno) %>% droplevels
  nb <- df_temp %>% subset(Type == "NB") %$% NormRoot
  sc <- df_temp %>% subset(Type == "SynCom") %$% NormRoot
  m1 <- t.test(nb,sc)
  Res <- data.frame(Genotype = geno,p.value = m1$p.value) %>%
    rbind(Res,.)
}

Res$p.adj <- Res$p.value %>% p.adjust(method = "fdr")
Res$Letter <- format.pval(Res$p.value, digits = 8)

Res_root <- Res


#### Dry  weight #######
df_sub <- mDat$DryWeight

genotypes <- df_sub$Genotype %>% unique
Res <- NULL
for(geno in genotypes){
  df_temp <- df_sub %>% subset(Genotype == geno) %>% droplevels
  nb <- df_temp %>% subset(Type == "NB") %$% Normweight_mg_plant
  sc <- df_temp %>% subset(Type == "SynCom") %$% Normweight_mg_plant
  m1 <- t.test(nb,sc)
  Res <- data.frame(Genotype = geno,p.value = m1$p.value) %>%
    rbind(Res,.)
}

Res$p.adj <- Res$p.value %>% p.adjust(method = "fdr")
Res$Letter <- format.pval(Res$p.value, digits = 8)
Res_dw <- Res

res <- list(Res_suberin = Res_suberin,
     Res_root = Res_root,
     Res_weight = Res_dw)
saveRDS(object = res,file = "../cleandata/res_syncom_mutants_suberin_root_weight.RDS")

res_original <- res
rm(res)


 ### With lavene test ###

mDat <- readRDS(file = "../cleandata/dat_syncom_genotypes.RDS") 


###### Suberin ##########

df_sub <- mDat$Suberin

genotypes <- df_sub$Genotype %>% unique
Res <- NULL
for(geno in genotypes){
  df_temp <- df_sub %>% subset(Genotype == geno) %>% droplevels
  mres <- car::leveneTest(NormSuberin~Type,df_temp)
  nb <- df_temp %>% subset(Type == "NB") %$% NormSuberin
  sc <- df_temp %>% subset(Type == "SynCom") %$% NormSuberin
  if(mres$`Pr(>F)`[1]< 0.05){
    m1 <- t.test(nb,sc,var.equal = FALSE)
    
  }else{
    m1 <- t.test(nb,sc,var.equal = TRUE)
    
  }
  Res <- data.frame(Genotype = geno,p.value = m1$p.value) %>%
    rbind(Res,.)
}

Res$p.adj <- Res$p.value %>% p.adjust(method = "fdr")
Res$Letter <- format.pval(Res$p.adj, digits = 8)

Res_suberin <- Res

#### Root#######
df_sub <- mDat$Root


genotypes <- df_sub$Genotype %>% unique
Res <- NULL
for(geno in genotypes){
  df_temp <- df_sub %>% subset(Genotype == geno) %>% droplevels
  mres <- car::leveneTest(NormRoot~Type,df_temp)
  nb <- df_temp %>% subset(Type == "NB") %$% NormRoot
  sc <- df_temp %>% subset(Type == "SynCom") %$% NormRoot
  if(mres$`Pr(>F)`[1]< 0.05){
    m1 <- t.test(nb,sc,var.equal = FALSE)
    
  }else{
    m1 <- t.test(nb,sc,var.equal = TRUE)
    
  } 
  Res <- data.frame(Genotype = geno,p.value = m1$p.value) %>%
    rbind(Res,.)
}

Res$p.adj <- Res$p.value %>% p.adjust(method = "fdr")
Res$Letter <- format.pval(Res$p.adj, digits = 8)

Res_root <- Res


#### Dry  weight #######
df_sub <- mDat$DryWeight

genotypes <- df_sub$Genotype %>% unique
Res <- NULL
for(geno in genotypes){
  df_temp <- df_sub %>% subset(Genotype == geno) %>% droplevels
  mres <- car::leveneTest(Normweight_mg_plant~Type,df_temp)
  nb <- df_temp %>% subset(Type == "NB") %$% Normweight_mg_plant
  sc <- df_temp %>% subset(Type == "SynCom") %$% Normweight_mg_plant
  if(mres$`Pr(>F)`[1]< 0.05){
    m1 <- t.test(nb,sc,var.equal = FALSE)
    
  }else{
    m1 <- t.test(nb,sc,var.equal = TRUE)
    
  }
  Res <- data.frame(Genotype = geno,p.value = m1$p.value) %>%
    rbind(Res,.)
}

Res$p.adj <- Res$p.value %>% p.adjust(method = "fdr")
Res$Letter <- format.pval(Res$p.adj, digits = 8)
Res_dw <- Res

res <- list(Res_suberin = Res_suberin,
            Res_root = Res_root,
            Res_weight = Res_dw)

plot(res_original$Res_suberin$p.adj,res$Res_suberin$p.adj)
cor.test(res_original$Res_suberin$p.adj,res$Res_suberin$p.adj)
cor.test(res_original$Res_root$p.adj,res$Res_root$p.adj)
cor.test(res_original$Res_weight$p.adj,res$Res_weight$p.adj)
saveRDS(object = res,file = "../cleandata/res_syncom_mutants_suberin_root_weight.RDS")



res$Res_root %>% subset(p.adj < 0.05)
res$Res_weight %>% subset(p.adj < 0.05)
res$Res_suberin %>% subset(p.adj < 0.05)
res$Res_suberin
res$Res_root
res$Res_weight
