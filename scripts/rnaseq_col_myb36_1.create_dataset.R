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
df <- read.table(file = "../rawdata/rnaseq_pipeline_mRNA_counts.txt",
                 header = T,row.names = 1,sep = "\t",skip = 1,
                 comment.char = "",quote = "")

#Create Tab
Tab <- df[,6:ncol(df)] %>% as.matrix
colnames(Tab) <- colnames(Tab) %>%
  gsub(pattern = "...mapped.",replacement = "") %>%
  gsub(pattern ="_.*",replacement = "" )
Tab <- Tab %>% t

Tab <- Tab[,which(colSums(Tab) != 0)]


rownames(Tab) <- paste0("S_",make.unique(rownames(Tab)))

Map <-data.frame(Sample_Id = rownames(Tab))
Map$Id <- Map$Sample_Id %>% 
  gsub(pattern = "^S_",replacement = "") %>%
  gsub(pattern = "\\..*",replacement = "")

Meta <- read.table("../rawdata/metadata_rnaseq_microrootbarriers.csv",header = T,sep = "\t")
Map <- merge(Map,Meta, by = "Id")
Map$Treatment <- Map$Treatment %>% gsub(pattern = " $",replacement = "")

Map$Genotype <- Map$Treatment %>% gsub(pattern = " .*",replacement = "") %>%
  gsub(pattern = "-",replacement = "_") %>%
  factor
Map$Bacteria <- Map$Treatment %>% gsub(pattern = ".* ",replacement = "") %>%
  factor


rownames(Map) <- Map$Sample_Id

Map <- Map %>% subset(Genotype == "myb36" | Genotype == "Col_0") %>% droplevels
toremove <- Map %>% subset(Rep == "Rep2" & Genotype == "myb36" & Bacteria == "NB" ) %$% Id %>%
  unique %>% c(96,.)
Map <- which(!(Map$Id %in% toremove)) %>% Map[.,] %>% droplevels

Tab <- Tab[match(Map$Sample_Id,rownames(Tab)),]

Tab <- Tab[,which(colSums(Tab) != 0)]


Map$group <- paste0(Map$Genotype,"_",Map$Bacteria)  %>%
  factor %>% relevel(ref = "Col_0_NB")

Dat <- create_dataset(Tab = Tab %>% t ,Map = Map)

#saveRDS(object = Dat,file = "../cleandata/dat_rnaseq_col_myb36.RDS")
Dat_sub <- Dat

dds <- DESeqDataSetFromMatrix(
  countData = Dat_sub$Tab,
  colData = Dat_sub$Map,
  design = ~  Rep + group
)



dds <- DESeq2::collapseReplicates(object = dds, dds$Id)



dds <- DESeq(object = dds,betaPrior = TRUE)


#Create object to plot
mvst <- vst(object = dds,blind = F)

#Remove the batch effect from the vst matrix
mat <- assay(mvst)

design0 <- model.matrix(~ 0 + group, colData(dds))
mat <- limma::removeBatchEffect(mat, mvst$Rep,design = design0)

Tab_z <- mat %>% t %>% scale (center = T,scale = F)

Mapa <- Dat_sub$Map[,c(1,3,4,5,6,7)] %>% unique
rownames(Mapa) <- Mapa$Id

Dat_z <- create_dataset(Tab = Tab_z %>% t,Map = Mapa)


#Prepare structure for plotting
melted <- Dat_z$Tab %>% melt
colnames(melted) <- c("Gene","Id","value")
melted <- match(melted$Id,Map$Id) %>%
  Map[.,-1] %>% cbind(melted,.)

melted_av <- dcast(melted,formula = Gene ~ group,
                   fun.aggregate = mean,value.var = "value") %>% melt
colnames(melted_av)[2] <- "group"


melted_av_rep <- dcast(melted,formula = Gene ~ group +Rep,
                       fun.aggregate = mean,value.var = "value") %>% melt
colnames(melted_av_rep)[2] <- "group"

melted_av <- match(melted_av$group,Map$group) %>%
  Map[.,5:6] %>% cbind(melted_av,.)


Map$UId <- paste0(Map$group,"_",Map$Rep)
melted_av_rep <- match(melted_av_rep$group,Map$UId) %>%
  Map[.,c(4,5,6)] %>% cbind(melted_av_rep,.)



### Save structures
res <- list(Dat = Dat,
     Dat_z = Dat_z,
     dds = dds,
     melted = melted,
     melted_av = melted_av,
     melted_av_rep = melted_av_rep
     )
saveRDS(object = res,file = "../cleandata/dat_rnaseq_col_myb36.RDS")
