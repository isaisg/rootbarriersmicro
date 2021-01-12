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

######## aba
df_aba <- read.table(file = "../cleandata/aag1550_Table_S1.txt",skip = 3,header = F,sep = "\t")
colnames(df_aba) <- c("Gene","h1","h4","h8","h12","h24","h36","h60")


rownames(df_aba) <- df_aba$Gene
mat_aba <- df_aba[,-1] %>% as.matrix


df_freq_up <-apply(mat_aba,MARGIN = 1,FUN = function(x)which((x > 2)) %>% length) %>% 
  data.frame(Gene = names(.),Freq = .,row.names = NULL)
df_freq_down<- apply(mat_aba,MARGIN = 1,FUN = function(x)which((x < -2)) %>% length) %>% 
  data.frame(Gene = names(.),Freq = .,row.names = NULL)



df_freq_up$Freq %>% table
df_freq_down$Freq %>% table


### Compute number of maps using different levels
mgenes_up <- df_freq_up %>% subset(Freq >= 6) %$% Gene %>% as.character
mgenes_down <- df_freq_down %>% subset(Freq >= 6) %$% Gene %>% as.character


Dat <-readRDS(file = "../cleandata/dat_rnaseq_col_myb36.RDS")
melted_av <- Dat$melted_av

which(melted_av$Gene %in% mgenes_down) %>%
  melted_av[.,] %>%
  chibi.boxplot(x_val = "Genotype",y_val = "value",style = "open"
                ,facet_formula = "Bacteria")

res <- readRDS(file = "../cleandata/res_rnaseq_col_myb36.RDS")
df_class <- res$df_classification

a <- length(mgenes_up )
b <- length(mgenes_down)

which(df_class$Gene %in% mgenes_up) %>% df_class$Classification[.] %>% table
x <- which(df_class$Gene %in% mgenes_up) %>% df_class$Classification[.] %>% table %>% sum
which(df_class$Gene %in% mgenes_down) %>% df_class$Classification[.] %>% table
y <- which(df_class$Gene %in% mgenes_down) %>% df_class$Classification[.] %>% table %>% sum

x/a
y/b

(x+y)/(a+b)

(x+y)/nrow(df_class)



## Save results
res <- list(
  aba_genes_up = mgenes_up,
  aba_genes_down = mgenes_down
)
saveRDS(object = res,file = "../cleandata/aba_robust_genes.RDS")