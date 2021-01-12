library(ohchibi)
library(DESeq2)
library(clusterProfiler)
library(org.At.tair.db)
library(scales)
library(egg)
library(ggpubr)
library(ggvenn)

setwd('/home/isai/Documents/results/rootbarriersmicro/scripts/')

set.seed(130816)
Dat <-readRDS(file = "../cleandata/dat_rnaseq_col_myb36.RDS")
Res <- readRDS(file = "../cleandata/res_rnaseq_col_myb36.RDS")


Res_aba <- readRDS(file = "../cleandata/aba_robust_genes.RDS")
chosen_up <- Res_aba$aba_genes_up
chosen_down <- Res_aba$aba_genes_down

df <- Res$df_classification
df_totals <- df$Classification %>% table %>% data.frame
colnames(df_totals)[1] <- c("Classification")

res_hyp <- NULL

########### Chosen up ################

mdf <- data.frame(Gene = chosen_up,
           Group = match(chosen_up,df$Gene) %>%df$Classification[.])
mdf$Group <- mdf$Group %>% as.character
mdf$Group[which(is.na(mdf$Group ))] <- "Unclassified"
df_table <- mdf$Group %>% table %>% as.data.frame
colnames(df_table)[1] <- "Classification"
df_prop <- merge(df_table,df_totals,by = "Classification")
df_prop$Percentage <- (df_prop$Freq.x/df_prop$Freq.y)*100

mdf$Group <- mdf$Group %>% factor(levels = c("bacteria_up_unique","bacteria_down_unique" ,"genotype_up_unique",
                                "genotype_down_unique","bacteria_up_genotype_up","bacteria_up_genotype_down",
"bacteria_down_genotype_up","bacteria_down_genotype_down","Unclassified"))

df_prop <- rbind(df_prop ,data.frame(Classification = "bacteria_up_genotype_down", Freq.x = 0,Freq.y = 0,Percentage = 0))

df_prop$Classification <-  df_prop$Classification %>%
  factor(levels = c("bacteria_up_unique","bacteria_down_unique" ,"genotype_up_unique",
                                            "genotype_down_unique","bacteria_up_genotype_up","bacteria_up_genotype_down",
                                            "bacteria_down_genotype_up","bacteria_down_genotype_down","Unclassified"))
  

pup <- ggplot(data = mdf %>% subset(Group != "Unclassified"),aes(Group)) +
  geom_bar(,fill = "black", color = "black") + 
  theme_ohchibi() +
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1))


pup_prop <- ggplot(data =df_prop,aes(Classification,y = Percentage)) +
  geom_bar(stat = "identity",fill = "black", color = "black") + 
  theme_ohchibi() +
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1))


##### Enrichment part #######
mgroups <- df_totals$Classification %>% levels
for(group in mgroups){
  sucess_sample <- mdf %>% subset(Group == group) %>% nrow
  sample_size <- df_totals %>% subset(Classification == group) %$% Freq
  sucess_background <- length(chosen_up)
  total_genes <- Dat$Dat$Tab %>% nrow
  mpval <- phyper(sucess_sample-1,sucess_background, total_genes, sample_size, lower.tail=FALSE);
  res_hyp <- data.frame(Classification = group,p.value = mpval,Type = "ABA_Up")   %>%
    rbind(res_hyp,.)
  
}


########### Chosen down ################

mdf <- data.frame(Gene = chosen_down,
                  Group = match(chosen_down,df$Gene) %>%df$Classification[.])
mdf$Group <- mdf$Group %>% as.character
mdf$Group[which(is.na(mdf$Group ))] <- "Unclassified"
df_table <- mdf$Group %>% table %>% as.data.frame
colnames(df_table)[1] <- "Classification"
df_prop <- merge(df_table,df_totals,by = "Classification")
df_prop$Percentage <- (df_prop$Freq.x/df_prop$Freq.y)*100



mdf$Group <- mdf$Group %>% factor(levels = c("bacteria_up_unique","bacteria_down_unique" ,"genotype_up_unique",
                                             "genotype_down_unique","bacteria_up_genotype_up","bacteria_up_genotype_down",
                                             "bacteria_down_genotype_up","bacteria_down_genotype_down","Unclassified"))

df_prop <- rbind(df_prop ,data.frame(Classification = "bacteria_down_genotype_up", Freq.x = 0,Freq.y = 0,Percentage = 0))
df_prop <- rbind(df_prop ,data.frame(Classification = "bacteria_down_genotype_down", Freq.x = 0,Freq.y = 0,Percentage = 0))


df_prop$Classification <-  df_prop$Classification %>%
  factor(levels = c("bacteria_up_unique","bacteria_down_unique" ,"genotype_up_unique",
                    "genotype_down_unique","bacteria_up_genotype_up","bacteria_up_genotype_down",
                    "bacteria_down_genotype_up","bacteria_down_genotype_down","Unclassified"))


pdown <- ggplot(data = mdf %>% subset(Group != "Unclassified"),aes(Group)) +
  geom_bar(,fill = "black", color = "black") + 
  theme_ohchibi() +
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1))


pdown_prop <- ggplot(data =df_prop,aes(Classification,y = Percentage)) +
  geom_bar(stat = "identity",fill = "black", color = "black") + 
  theme_ohchibi() +
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1))


##### Enrichment part #######
mgroups <- df_totals$Classification %>% levels
for(group in mgroups){
  sucess_sample <- mdf %>% subset(Group == group) %>% nrow
  sample_size <- df_totals %>% subset(Classification == group) %$% Freq
  sucess_background <- length(chosen_up)
  total_genes <- Dat$Dat$Tab %>% nrow
  mpval <- phyper(sucess_sample-1,sucess_background, total_genes, sample_size, lower.tail=FALSE);
  res_hyp <- data.frame(Classification = group,p.value = mpval,Type = "ABA_Down")   %>%
    rbind(res_hyp,.)
  
}


#### Ethylene up ###
mregs <- readRDS(file = "../../variovorax/rnaseq/cleandata/regulons_rnaseq.RDS")
mgenes <- mregs$ethylene

mdf <- data.frame(Gene = mgenes,
                  Group = match(mgenes,df$Gene) %>%df$Classification[.])
mdf$Group <- mdf$Group %>% as.character
mdf$Group[which(is.na(mdf$Group ))] <- "Unclassified"
df_table <- mdf$Group %>% table %>% as.data.frame
colnames(df_table)[1] <- "Classification"
df_prop <- merge(df_table,df_totals,by = "Classification")
df_prop$Percentage <- (df_prop$Freq.x/df_prop$Freq.y)*100


mdf$Group <- mdf$Group %>% factor(levels = c("bacteria_up_unique","bacteria_down_unique" ,"genotype_up_unique",
                                             "genotype_down_unique","bacteria_up_genotype_up","bacteria_up_genotype_down",
                                             "bacteria_down_genotype_up","bacteria_down_genotype_down","Unclassified"))


df_prop$Classification <-  df_prop$Classification %>%
  factor(levels = c("bacteria_up_unique","bacteria_down_unique" ,"genotype_up_unique",
                    "genotype_down_unique","bacteria_up_genotype_up","bacteria_up_genotype_down",
                    "bacteria_down_genotype_up","bacteria_down_genotype_down","Unclassified"))

peth<- ggplot(data = mdf %>% subset(Group != "Unclassified"),aes(Group)) +
  geom_bar(,fill = "black", color = "black") + 
  theme_ohchibi() +
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1))

peth_prop <- ggplot(data =df_prop,aes(Classification,y = Percentage)) +
  geom_bar(stat = "identity",fill = "black", color = "black") + 
  theme_ohchibi() +
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1))



##### Enrichment part #######
mgroups <- df_totals$Classification %>% levels
for(group in mgroups){
  sucess_sample <- mdf %>% subset(Group == group) %>% nrow
  sample_size <- df_totals %>% subset(Classification == group) %$% Freq
  sucess_background <- length(chosen_up)
  total_genes <- Dat$Dat$Tab %>% nrow
  mpval <- phyper(sucess_sample-1,sucess_background, total_genes, sample_size, lower.tail=FALSE);
  res_hyp <- data.frame(Classification = group,p.value = mpval,Type = "Ethylene_Up")   %>%
    rbind(res_hyp,.)
  
}

### EMBO 
df_embo <- read.table(file = "../rawdata/embo_930genes.csv",header = T,sep = "\t")
mgenes <- df_embo %>% subset(cluster2 == 2) %$% rowname %>% as.character

mdf <- data.frame(Gene = mgenes,
                  Group = match(mgenes,df$Gene) %>%df$Classification[.])
mdf$Group <- mdf$Group %>% as.character
mdf$Group[which(is.na(mdf$Group ))] <- "Unclassified"
df_table <- mdf$Group %>% table %>% as.data.frame
colnames(df_table)[1] <- "Classification"
df_prop <- merge(df_table,df_totals,by = "Classification")
df_prop$Percentage <- (df_prop$Freq.x/df_prop$Freq.y)*100


df_prop <- rbind(df_prop ,data.frame(Classification = "bacteria_down_genotype_down", Freq.x = 0,Freq.y = 0,Percentage = 0))

mdf$Group <- mdf$Group %>% factor(levels = c("bacteria_up_unique","bacteria_down_unique" ,"genotype_up_unique",
                                             "genotype_down_unique","bacteria_up_genotype_up","bacteria_up_genotype_down",
                                             "bacteria_down_genotype_up","bacteria_down_genotype_down","Unclassified"))
df_prop$Classification <-  df_prop$Classification %>%
  factor(levels = c("bacteria_up_unique","bacteria_down_unique" ,"genotype_up_unique",
                    "genotype_down_unique","bacteria_up_genotype_up","bacteria_up_genotype_down",
                    "bacteria_down_genotype_up","bacteria_down_genotype_down","Unclassified"))

pembo<- ggplot(data = mdf %>% subset(Group != "Unclassified"),aes(Group)) +
  geom_bar(,fill = "black", color = "black") + 
  theme_ohchibi() +
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1))

pembo_prop <- ggplot(data =df_prop,aes(Classification,y = Percentage)) +
  geom_bar(stat = "identity",fill = "black", color = "black") + 
  theme_ohchibi() +
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1))



##### Enrichment part #######
mgroups <- df_totals$Classification %>% levels
for(group in mgroups){
  sucess_sample <- mdf %>% subset(Group == group) %>% nrow
  sample_size <- df_totals %>% subset(Classification == group) %$% Freq
  sucess_background <- length(chosen_up)
  total_genes <- Dat$Dat$Tab %>% nrow
  mpval <- phyper(sucess_sample-1,sucess_background, total_genes, sample_size, lower.tail=FALSE);
  res_hyp <- data.frame(Classification = group,p.value = mpval,Type = "EMBO")   %>%
    rbind(res_hyp,.)
  
}

res_hyp$p.value %>% hist
res_hyp$p.adj <- res_hyp$p.value %>% p.adjust(method = "fdr") 


res_hyp %>% subset(p.adj < 0.05)

composition <- egg::ggarrange(pup,pup_prop,
               pdown,pdown_prop,
               peth,peth_prop,
               pembo,pembo_prop,nrow = 4,ncol = 2)
oh.save.pdf(p = composition,outname = "composition_bars.pdf",outdir = "../figures/",width = 30,height = 50)
