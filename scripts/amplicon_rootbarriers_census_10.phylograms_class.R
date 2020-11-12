library(ohchibi)
library(palettesPM)
library(Rmisc)
library(ggpubr)
library(DESeq2)
library(dendextend)
library(multcomp)
library(limma)

set.seed(130816)
setwd('/home/isai/Documents/results/rootbarriersmicro/scripts/')
dir.create("../cleandata")
dir.create("../figures")


Dat <- readRDS(file = "../cleandata/dat_censusmutants_16s.RDS")

Dat_rab <- Dat$RelativeAbundance

Tab <- Dat_rab$Tab
Tab_batch <- removeBatchEffect(x = Tab,Dat_rab$Map$Rep)

Dat_batch <- create_dataset(Tab = Tab_batch,Map = Dat_rab$Map,Tax = Dat_rab$Tax)



Dat_class <- Dat_batch %>% collapse_by_taxonomy.Dataset(Dat = .,level = 4)

df_freq <- Dat_class$Tab %>% rowMeans %>% data.frame(Id = names(.),RA = .,row.names = NULL)
df_freq <- with(df_freq,order(-RA)) %>% df_freq[.,]


mchosen <- df_freq$Id[1:12] %>% as.character
Tab_chosen <- which(rownames(Dat_class$Tab) %in% mchosen) %>% Dat_class$Tab[.,] %>% as.matrix

Tab_other <- which(!(rownames(Dat_class$Tab) %in% mchosen)) %>% Dat_class$Tab[.,] %>% as.matrix
Other <- colSums(Tab_other)
Tab_chosen <- rbind(Tab_chosen,Other) 

melted <- Tab_chosen %>% melt
melted$Var1 <- melted$Var1 %>% gsub(pattern = "Root; k__Bacteria; ",replacement = "") 
colnames(melted) <- c("Id","DADA2_Header","RA")

melted <- merge(melted,Dat_class$Map, by = "DADA2_Header")
mchosen <- mchosen %>% gsub(pattern = "Root; k__Bacteria; ",replacement = "")
melted$Id <- melted$Id %>% factor(levels = c("p__Acidobacteria; c__Acidobacteriia",
                                "p__Actinobacteria; c__Thermoleophilia",
                                "p__Actinobacteria; c__Actinobacteria",
                                "p__Bacteroidetes; c__Bacteroidia",
                                "p__Chloroflexi; c__Chloroflexia",
                                "p__Cyanobacteria; c__Melainabacteria",
                                "p__Firmicutes; c__Bacilli",
                                "p__Firmicutes; c__Negativicutes",
                                "p__Proteobacteria; c__Alphaproteobacteria",
                                "p__Proteobacteria; c__Gammaproteobacteria",
                                "p__Proteobacteria; c__Deltaproteobacteria",
                                "p__Verrucomicrobia; c__Verrucomicrobiae",
                                "Other"))
paleta <- c(RColorBrewer::brewer.pal(n = 12,name = "Paired"),"grey")
names(paleta) <- melted$Id %>% levels

## Read the  deseq 2 results ############
df_deseq <- read.table(file = "../cleandata/deseq2_results_censusmutants_genotypes_inside_fractions_all.tsv",header = T,sep = "\t")
df_deseq <- df_deseq %>% subset(Level == "Class") %>% droplevels
df_deseq$Id <- df_deseq$Id %>% gsub(pattern = "Root; k__Bacteria; ",replacement = "")
df_deseq <- df_deseq[which(df_deseq$Id %in% mchosen),] %>% droplevels
df_deseq$padjClass <- df_deseq$pvalue %>% p.adjust(method = "fdr")
df_deseq$SignificanceClass <- NA
df_deseq$SignificanceClass[which(df_deseq$padjClass < 0.1)] <- "S"



df_sum <- aggregate(RA ~Id+Genotype+Fraction,melted,mean)

df_sum <- merge(df_sum,df_deseq[,c(5,7,8,10,11,12,13,14)], by = c("Id","Fraction","Genotype"),all.x = TRUE) 

df_sum$RAAbs <- abs(df_sum$RA)


df_sum %>% subset(SignificanceClass == "S")

### Adjust the panels
df_groups <- read.table(file = "../rawdata/df_genotypes_to_classification.csv",header = T,sep = ",")
df_groups$Num <- 1:20

df_sum <- merge(df_sum,df_groups, by = "Genotype") 
df_sum$Num <- df_sum$Num %>% factor
p <- ggplot(data = df_sum,aes(Num,RAAbs)) +
  geom_bar(stat = "identity",aes(fill = Id),position = "fill") +
  facet_grid(.~Fraction) +
  scale_fill_manual(values = paleta,na.value = "#D9D9D9") +
  theme_ohchibi() +
  scale_y_continuous(breaks = seq(0,1,by = 0.1),expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme(
    axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1),
    strip.background.x = element_blank(),
    strip.text.x = element_text(family = "Arial",face = "bold",size = 20)
  ) +
  ylab(label = "Relative abundance") 


oh.save.pdf(p = p,outname = "phylogram_class_level.pdf",outdir = "../figures/",width = 24,height = 12)


df_sum %>% subset(SignificanceClass == "S")
