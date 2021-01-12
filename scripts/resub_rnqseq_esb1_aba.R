library(ohchibi)
library(DESeq2)
library(clusterProfiler)
library(org.At.tair.db)
library(scales)
library(egg)
library(ggpubr)

setwd('/Users/isaisalasgonzalez/Documents/rootbarriersmicro/scripts/')

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

toremove <- Map %>% subset(Rep == "Rep2" & Genotype == "myb36" & Bacteria == "NB" ) %$% Id %>%
  unique %>% c(96,.)

toremove <- Map %>% subset(Genotype == "esb1" & Rep == "Rep2") %$% Id %>% unique %>%
  c(toremove,.)

Map <- which(!(Map$Id %in% toremove)) %>% Map[.,] %>% droplevels


#Keep esb1 first rep only


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
a <- plotPCA(mvst,c("Genotype"))
b <- plotPCA(mvst,c("Bacteria")) + scale_color_manual(values = c("black","grey"))
c <- plotPCA(mvst,c("Rep"))
egg::ggarrange(a,b,c,nrow = 1)

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


## read aba genes ##
Res_aba <- readRDS(file = "../cleandata/aba_robust_genes.RDS")
chosen_up <- Res_aba$aba_genes_up
chosen_down <- Res_aba$aba_genes_down

melted_av$Type <- "ABA upregulated"

a <- which(melted_av$Gene %in% chosen_up) %>%
  melted_av[.,] %>%
  subset(value < 1  & value > -1) %>%
  chibi.boxplot(Map = .,x_val = "Genotype",
                y_val = "value",col_val = "Bacteria",facet_formula = "Type",
                strip_text_size = 20,mpalette = c("#525252","#b48c36"),
                median_colored_as_points = F,size_median = 2,alpha_point = 1,size_boxplot = 1) 


melted_av$Type <- "ABA downregulated"

b <- which(melted_av$Gene %in% chosen_down) %>%
  melted_av[.,] %>%
  subset(value < 1  & value > -1) %>%
  chibi.boxplot(Map = .,x_val = "Genotype",
                y_val = "value",col_val = "Bacteria",facet_formula = "Type",
                strip_text_size = 20,mpalette = c("#525252","#b48c36"),
                median_colored_as_points = F,size_median = 2,alpha_point = 1,size_boxplot = 1) 


reyt_genes <- readRDS(file = "../cleandata/reyt_up_genes.RDS")
melted_av$Type <- "Reyt esb1 myb36 Up"

c <- which(melted_av$Gene %in% reyt_genes) %>%
  melted_av[.,] %>%
  subset(value < 1  & value > -1) %>%
  chibi.boxplot(Map = .,x_val = "Genotype",
                y_val = "value",col_val = "Bacteria",facet_formula = "Type",
                strip_text_size = 20,mpalette = c("#525252","#b48c36"),
                median_colored_as_points = F,size_median = 2,alpha_point = 1,size_boxplot = 1) 

