library(ohchibi)
library(scales)
library(Rmisc)
library(ggrepel)

set.seed(130816)
setwd('/Users/isaisalasgonzalez/Documents/rootbarriersmicro/scripts/')

#Read microbiome plots 
Res_microbiome <- readRDS(file = "../cleandata/plots_correlations_ionome_vs_amplicon_cap_agarsystem.RDS")

## Load the number scheme
df_num <- read.table(file = "../rawdata/df_genotypes_to_classification.csv",header = T,sep = ",")
df_num$Number <- 1:20

#Read the feature table
df <- read.table(file = "../rawdata/Gabriel_3000_0percent_Q85_reSCALED.tsv",
                 header = T,sep = "\t",comment.char = "",quote = "")
Tab <- colnames(df) %>% grep(pattern = "NEG") %>%
  df[,.] %>% as.matrix
rownames(Tab) <- df$ID
colnames(Tab) <- colnames(Tab) %>% gsub(pattern = "NEG",replacement = "")

Map <- read.table(file = "../rawdata/metadata_metabolomics.tsv",header = T,sep = "\t")
rownames(Map) <- Map$Order

Map$Genotype <- Map$Genotype %>% 
  gsub(pattern = "^wt$",replacement = "Col_0") %>%
  gsub(pattern = " ",replacement = "")  %>%
  gsub(pattern = "^myb36$",replacement = "myb36_2") %>%
  gsub(pattern = "^esb1$",replacement = "esb1_1") %>%
  gsub(pattern = "^sgn3$",replacement = "sgn3_3") %>%
  gsub(pattern = "^sgn3/myb36$",replacement = "myb36_2_sgn3_3") %>%
  gsub(pattern = "^CDEF\\(wt\\)$",replacement = "pCASP1_CDEFwt") %>%
  gsub(pattern = "^CDEF\\(esb1\\)$",replacement = "pCASP1_CDEFesb1_1") %>%
  gsub(pattern = "CDEF\\(sgn3\\/myb36\\)",replacement = "pELTP_CDEFsgn3_3_myb36_1")

Map$Number <- match(Map$Genotype,df_num$Genotype) %>%
  df_num$Number[.]
Map$Number <- Map$Number %>% factor

Tab <- match(Map$Order,colnames(Tab)) %>%
  Tab[,.]


Tab_z <- Tab %>% t %>% scale %>% t

#Create dataset
Dat <- create_dataset(Tab = Tab_z,Map = Map)

## Perform PCA ##
distfun <- function(x,method) vegan::vegdist(x = x, method = "euclidean")

mpca <- oh.pca(Tab = Dat$Tab %>% t,Map = Dat$Map,retx = T,center = F,scale = F,id_var = "Order")
# chibi.pca(list_ohpca = mpca,col_val = "Number",shape_val = "Bacteria") +
#   scale_color_brewer(palette = "Paired")
df_sum_pc1 <- mpca$Map_pca %>%
  summarySE(data = .,measurevar = "PC1",groupvars = c("Genotype","Bacteria"))

df_sum_pc2 <- mpca$Map_pca %>%
  summarySE(data = .,measurevar = "PC2",groupvars = c("Genotype","Bacteria"))

merged <- merge(df_sum_pc1[,c("Genotype","Bacteria","PC1","ci")],
      df_sum_pc2[,c("Genotype","Bacteria","PC2","ci")],
      by = c("Genotype","Bacteria"))

colnames(merged)[c(4,6)] <- c("PC1.ci","PC2.ci")

#Append the number of the genotype
merged$Number <- match(merged$Genotype,df_num$Genotype) %>%
  df_num$Number[.]

p <- ggplot(data = merged,aes(PC1,PC2)) +
  geom_vline(xintercept = 0,size = 0.5, color = "#D9D9D9", linetype = "longdash") + 
  geom_hline(yintercept = 0, size = 0.5, color = "#D9D9D9", linetype = "longdash") +
  geom_linerange(aes(xmin = PC1 - PC1.ci,xmax = PC1 + PC1.ci),color = "#D9D9D9",alpha = 0.5) +
  geom_linerange(aes(ymin = PC2 - PC2.ci,ymax = PC2 + PC2.ci),color = "#D9D9D9",alpha = 0.5)  +
  geom_point(aes(shape = Bacteria),size = 4)  +
  geom_text_repel(aes(label = Number),size = 8) +
  scale_color_brewer(palette = "Paired") +
  scale_x_continuous(limits = c(-120,120),oob = rescale_none) +
  scale_y_continuous(limits = c(-120,120),oob = rescale_none) + 
  theme_ohchibi(font_family = "Helvetica") +
  xlab(label = "PC1 15.64%") +
  ylab(label = "PC2 11.08%") +
  scale_shape_manual(values = c(17,16))

###Test the genotypes again its respective  Col-0
df_test <- merged %>% subset(Bacteria == "NB") %>% droplevels
mgenos <- df_test$Number %>% grep(pattern = "^1$",invert = T,value = T)
df_test <- mpca$Map_pca %>% subset(Bacteria == "NB") %>% droplevels
mvars <- c("PC1","PC2")
Res <- NULL
for(mv in  mvars){
  df_wt <- df_test %>% subset(Number == 1)
  x <- df_wt[,mv]
  for(geno in mgenos){
    df_temp <- df_test %>% subset(Number == geno)
    y <- df_temp[,mv]
    m1 <- t.test(x,y)
    pval <- m1$p.value
    Res <- data.frame(Number = geno,Var = mv,p.value = pval,Type = "NB") %>%
      rbind(Res,.)
  }
}

df_test <- merged %>% subset(Bacteria != "NB") %>% droplevels
mgenos <- df_test$Number %>% grep(pattern = "^1$",invert = T,value = T)
df_test <- mpca$Map_pca %>% subset(Bacteria != "NB") %>% droplevels
mvars <- c("PC1","PC2")
for(mv in  mvars){
  df_wt <- df_test %>% subset(Number == 1)
  x <- df_wt[,mv]
  for(geno in mgenos){
    df_temp <- df_test %>% subset(Number == geno)
    y <- df_temp[,mv]
    m1 <- t.test(x,y)
    pval <- m1$p.value
    Res <- data.frame(Number = geno,Var = mv,p.value = pval,Type = "SC") %>%
      rbind(Res,.)
  }
}

##PERMANOVA ##
mdist <- distfun(Dat$Tab %>% t)
adonis2(formula = mdist~Genotype+Bacteria+Genotype:Bacteria,data = Dat$Map,permutations = 9999)

### Try correlation with the microbiome ####
data_metabolome_root <- merged %>% subset(Bacteria == "SC")
Tab_metabolome <- data_metabolome_root[,c("PC1","PC2")]
rownames(Tab_metabolome) <- data_metabolome_root$Number

#Root microbiome 
data_micro <- Res_microbiome$agar_constellation_root_amplicon$data
Tab_micro <- data_micro[,c("CAP1","CAP2")] %>% as.matrix
rownames(Tab_micro) <- data_micro$Number
dist_micro <- dist(Tab_micro)
dist_metabolome <- dist(Tab_metabolome)

m_root <- mantel(xdis = dist_micro,dist_metabolome,permutations = 4999)
mtext_root <- paste0("r = ",round(m_root$statistic,3),"\np = ",format.pval(m_root$signif))

### Plot correlations ###
merged_mantel <- dist_micro %>% as.matrix %>% melt_dist() %>%
  merge(
    dist_metabolome %>% as.matrix %>% melt_dist(),
    by = c("iso1","iso2")
    
  )
colnames(merged_mantel)[3:4] <- c("Root_microbiome","Metabolome")


p1 <- ggplot(merged_mantel,aes(Metabolome,Root_microbiome)) +
  geom_smooth(method = "lm",color = "red",size = 2,fill = "#D9D9D9") +
  geom_point(size = 3) +
  theme_ohchibi(size_panel_border = 2) +
  annotate(geom = "text",x = 2,y = 0.6,label = mtext_root) +
  #scale_y_continuous(limits = c(0,3.5)) +
  #scale_x_continuous(limits = c(0,2.5),oob = rescale_none) +
  xlab(label = "Metabolome") + ylab(label = "Root microbiome")

oh.save.pdf(p = p + theme(legend.position = "none"),outname = "resub_metabolome_pca_51818.pdf",outdir = "../figures/",width = 5,height = 5)
oh.save.pdf(p = p1,outname = "resub_metabolome_rootmicrobiome_mantel.pdf",outdir = "../figures/",width = 5,height = 5)
#Shoot microbiome 
data_micro <- Res_microbiome$agar_constellation_shoot_amplicon$data
Tab_micro <- data_micro[,c("CAP1","CAP2")] %>% as.matrix
rownames(Tab_micro) <- data_micro$Number
dist_micro <- dist(Tab_micro)
dist_metabolome <- dist(Tab_metabolome)

mantel(xdis = dist_micro,dist_metabolome,permutations = 4999)

#Shoot microbiome 
data_micro <- Res_microbiome$agar_constellation_agar_amplicon$data
Tab_micro <- data_micro[,c("CAP1","CAP2")] %>% as.matrix
rownames(Tab_micro) <- data_micro$Number
dist_micro <- dist(Tab_micro)
dist_metabolome <- dist(Tab_metabolome)

mantel(xdis = dist_micro,dist_metabolome,permutations = 4999)



### Read the second feature table ###
df <- read.table(file = "../rawdata/Gabriel_1000_100percent_reSCALED.tsv",
                 header = T,sep = "\t",comment.char = "",quote = "")
Tab <- colnames(df) %>% grep(pattern = "NEG") %>%
  df[,.] %>% as.matrix
rownames(Tab) <- df$ID
colnames(Tab) <- colnames(Tab) %>% gsub(pattern = "NEG",replacement = "")

Map <- read.table(file = "../rawdata/metadata_metabolomics.tsv",header = T,sep = "\t")
rownames(Map) <- Map$Order

Map$Genotype <- Map$Genotype %>% 
  gsub(pattern = "^wt$",replacement = "Col_0") %>%
  gsub(pattern = " ",replacement = "")  %>%
  gsub(pattern = "^myb36$",replacement = "myb36_2") %>%
  gsub(pattern = "^esb1$",replacement = "esb1_1") %>%
  gsub(pattern = "^sgn3$",replacement = "sgn3_3") %>%
  gsub(pattern = "^sgn3/myb36$",replacement = "myb36_2_sgn3_3") %>%
  gsub(pattern = "^CDEF\\(wt\\)$",replacement = "pCASP1_CDEFwt") %>%
  gsub(pattern = "^CDEF\\(esb1\\)$",replacement = "pCASP1_CDEFesb1_1") %>%
  gsub(pattern = "CDEF\\(sgn3\\/myb36\\)",replacement = "pELTP_CDEFsgn3_3_myb36_1")

Map$Number <- match(Map$Genotype,df_num$Genotype) %>%
  df_num$Number[.]
Map$Number <- Map$Number %>% factor

Tab <- match(Map$Order,colnames(Tab)) %>%
  Tab[,.]


Tab_z <- Tab %>% t %>% scale %>% t

#Create dataset
Dat <- create_dataset(Tab = Tab_z,Map = Map)

## Perform PCA ##
distfun <- function(x,method) vegan::vegdist(x = x, method = "euclidean")

mpca <- oh.pca(Tab = Dat$Tab %>% t,Map = Dat$Map,retx = T,center = F,scale = F,id_var = "Order")
# chibi.pca(list_ohpca = mpca,col_val = "Number",shape_val = "Bacteria") +
#   scale_color_brewer(palette = "Paired")
df_sum_pc1 <- mpca$Map_pca %>%
  summarySE(data = .,measurevar = "PC1",groupvars = c("Genotype","Bacteria"))

df_sum_pc2 <- mpca$Map_pca %>%
  summarySE(data = .,measurevar = "PC2",groupvars = c("Genotype","Bacteria"))

merged <- merge(df_sum_pc1[,c("Genotype","Bacteria","PC1","ci")],
                df_sum_pc2[,c("Genotype","Bacteria","PC2","ci")],
                by = c("Genotype","Bacteria"))

colnames(merged)[c(4,6)] <- c("PC1.ci","PC2.ci")

#Append the number of the genotype
merged$Number <- match(merged$Genotype,df_num$Genotype) %>%
  df_num$Number[.]

ggplot(data = merged,aes(PC1,PC2)) +
  geom_vline(xintercept = 0,size = 0.5, color = "#D9D9D9", linetype = "longdash") + 
  geom_hline(yintercept = 0, size = 0.5, color = "#D9D9D9", linetype = "longdash") +
  geom_linerange(aes(xmin = PC1 - PC1.ci,xmax = PC1 + PC1.ci),color = "#D9D9D9",alpha = 0.5) +
  geom_linerange(aes(ymin = PC2 - PC2.ci,ymax = PC2 + PC2.ci),color = "#D9D9D9",alpha = 0.5)  +
  geom_point(aes(shape = Bacteria),size = 4)  +
  geom_text_repel(aes(label = Number),size = 8) +
  scale_color_brewer(palette = "Paired") +
  scale_x_continuous(limits = c(-50,50),oob = rescale_none) +
  scale_y_continuous(limits = c(-35,35),oob = rescale_none) + 
  theme_ohchibi(font_family = "Helvetica") +
  xlab(label = "PC1 33.46%") +
  ylab(label = "PC2 18.31%") +
  scale_shape_manual(values = c(17,16))

### Try correlation with the microbiome ####
data_metabolome_root <- merged %>% subset(Bacteria == "SC")
Tab_metabolome <- data_metabolome_root[,c("PC1","PC2")]
rownames(Tab_metabolome) <- data_metabolome_root$Number

#Root microbiome 
data_micro <- Res_microbiome$agar_constellation_root_amplicon$data
Tab_micro <- data_micro[,c("CAP1","CAP2")] %>% as.matrix
rownames(Tab_micro) <- data_micro$Number
dist_micro <- dist(Tab_micro)
dist_metabolome <- dist(Tab_metabolome)

mantel(xdis = dist_micro,dist_metabolome,permutations = 4999)

#Shoot microbiome 
data_micro <- Res_microbiome$agar_constellation_shoot_amplicon$data
Tab_micro <- data_micro[,c("CAP1","CAP2")] %>% as.matrix
rownames(Tab_micro) <- data_micro$Number
dist_micro <- dist(Tab_micro)
dist_metabolome <- dist(Tab_metabolome)

mantel(xdis = dist_micro,dist_metabolome,permutations = 4999)

#Shoot microbiome 
data_micro <- Res_microbiome$agar_constellation_agar_amplicon$data
Tab_micro <- data_micro[,c("CAP1","CAP2")] %>% as.matrix
rownames(Tab_micro) <- data_micro$Number
dist_micro <- dist(Tab_micro)
dist_metabolome <- dist(Tab_metabolome)

mantel(xdis = dist_micro,dist_metabolome,permutations = 4999)
