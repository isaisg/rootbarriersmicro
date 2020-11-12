library(ohchibi)
library(palettesPM)
library(Rmisc)
library(ggpubr)
library(DESeq2)
library(dendextend)

set.seed(130816)
setwd('/home/isai/Documents/results/rootbarriersmicro/scripts/')
dir.create("../cleandata")
dir.create("../figures")


####### DEseq2 models testing for genotype effect  ##

Dat <- readRDS(file = "../cleandata/dat_censusmutants_16s.RDS")
Dat_raw <- Dat$RawCounts

deseq_function <- function(Dat_model = Dat,ref = "Col_0"){
  #DEseq2 model
  dds <- DESeqDataSetFromMatrix(countData = Dat_model$Tab,colData = Dat_model$Map,
                                design =~Rep + Genotype )
  dds <-tryCatch(expr = DESeq(dds),error = function(e)Genotype = NA)
  if(is.na(dds)){
    dds <- DESeqDataSetFromMatrix(countData = Dat_model$Tab,colData = Dat_model$Map,
                                  design =~Rep + Genotype )
    gm_mean = function(x, na.rm=TRUE){
      exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
    }
    geoMeans <- apply(counts(dds), 1, gm_mean)
    dds <- estimateSizeFactors(dds, geoMeans = geoMeans)
    dds <- DESeq(dds, fitType="local")
  }
  #Test everything against Col-0
  Res <- NULL
  mgenos <- Dat_model$Map$Genotype %>% as.character %>%
    unique %>% grep(pattern = ref,value = T,invert = T)
  for(geno in mgenos){
    cat("Working on ",geno,"\n")
    mcontrast <- c("Genotype",geno,ref)
    temp <- results(object = dds,contrast = mcontrast) %>% as.data.frame
    temp$Id <- rownames(temp)
    rownames(temp) <- NULL
    temp$Genotype <- geno
    Res <- rbind(Res,temp)
  }
  return(Res)
}



######### Root #############
Dat_sub <- Dat_raw %>% subset.Dataset(Fraction == "Root",drop = T,clean = T)

mDats <- list()
counter <- 0
for(i in 3:7){
  Dat_temp <- collapse_by_taxonomy.Dataset(Dat = Dat_sub,level = i)
  counter <- counter+1
  mDats[[counter]] <- Dat_temp
}
mres <- lapply(X = mDats,FUN = deseq_function)

## apply it over the otus
total <- Dat_sub$Tab %>% sum
Dat_filter <- measurable_taxa.Dataset(Dat = Dat_sub,
                                      min_samples_otu =15,min_reads_otu = 15 )
(Dat_filter$Tab %>% sum)/total

res_otu <- deseq_function(Dat_model = Dat_filter)

mres[[6]] <- res_otu
mres_root <- mres
rm(mres)

######### Shoot #############
Dat_sub <- Dat_raw %>% subset.Dataset(Fraction == "Shoot",drop = T,clean = T)

mDats <- list()
counter <- 0
for(i in 3:7){
  Dat_temp <- collapse_by_taxonomy.Dataset(Dat = Dat_sub,level = i)
  counter <- counter+1
  mDats[[counter]] <- Dat_temp
}

mres <- lapply(X = mDats,FUN = deseq_function)

## apply it over the otus
total <- Dat_sub$Tab %>% sum
Dat_filter <- measurable_taxa.Dataset(Dat = Dat_sub,
                                      min_samples_otu =15,min_reads_otu = 15 )
(Dat_filter$Tab %>% sum)/total

res_otu <- deseq_function(Dat_model = Dat_filter)


mres[[6]] <- res_otu
mres_shoot <- mres
rm(mres)

######### Soil #############
Dat_sub <- Dat_raw %>% subset.Dataset(Fraction == "Soil",drop = T,clean = T)

mDats <- list()
counter <- 0
for(i in 3:7){
  Dat_temp <- collapse_by_taxonomy.Dataset(Dat = Dat_sub,level = i)
  counter <- counter+1
  mDats[[counter]] <- Dat_temp
}

mres <- lapply(X = mDats,FUN = deseq_function)

## apply it over the otus
total <- Dat_sub$Tab %>% sum
Dat_filter <- measurable_taxa.Dataset(Dat = Dat_sub,
                                      min_samples_otu =15,min_reads_otu = 15 )
(Dat_filter$Tab %>% sum)/total

res_otu <- deseq_function(Dat_model = Dat_filter)

mres[[6]] <- res_otu
mres_soil <- mres
rm(mres)


mlista <- list(Soil = mres_soil,
     Root = mres_root,
     Shoot = mres_shoot)

saveRDS(object = mlista,file = "../cleandata/deseq2_results_censusmutants_genotypes_inside_fractions.RDS")