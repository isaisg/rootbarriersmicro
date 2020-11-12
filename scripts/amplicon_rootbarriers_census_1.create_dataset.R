library(ohchibi)

set.seed(130816)
setwd('/home/isai/Documents/results/rootbarriersmicro/scripts/')
dir.create("../cleandata")
dir.create("../figures")

#Read Tab
Tab <- read.table(file = "../rawdata/dada2_res_structure_merged_mutants.tsv",
                  header = T,sep = "\t",check.names = F,quote = "",row.names = 1)

#Read metadata
Map <- read.table(file = "../rawdata/amplicon_rootbarriers_metadata_mutants.tsv",
                  header = T,sep = "\t",quote = "",comment.char = "")
Map$DADA2_Header <- Map$DADA2_Header %>% gsub(pattern = "_",replacement = "")
rownames(Map) <- Map$DADA2_Header

#Create dataset
usable_samples <- intersect(Map$DADA2_Header,rownames(Tab)) 
Map <- match(usable_samples,Map$DADA2_Header) %>%
  Map[.,] %>% droplevels
Tab <- match(usable_samples,rownames(Tab)) %>%
  Tab[.,]
Tab <- Tab[,-which(colSums(Tab) ==0)]

Map$Genotype <- Map$Genotype %>%
  gsub(pattern = "pELTP:CDEFmyb36-1_sgn3-3_",replacement = "pELTP:CDEFsgn3-3_myb36-1") %>%
  factor %>% relevel(ref = "Col-0")

#Read the taxonomy
Tax <- read.table(file = "../rawdata/seqsformothur.full_v132.wang.mutants.taxonomy",header = F,sep = "\t")

colnames(Tax) <- c("ID","Taxonomy")
rownames(Tax) <- Tax$ID
Tax$Taxonomy <- Tax$Taxonomy %>% as.character %>% 
  gsub(pattern = "\\([0-9]+\\)",replacement = "",perl = T)
ntax <- NULL
for(i in 1:nrow(Tax)){
  ttax <- as.character(Tax[i,2])
  noms <- unlist(strsplit(x = ttax,split = ";"))
  fnom <- paste0("Root; k__",noms[1],"; p__",noms[2],"; c__",noms[3],"; o__",noms[4],"; f__",noms[5],"; g__",noms[6])
  ntax <- c(ntax,fnom)
}
Tax$Taxonomy <- factor(ntax)

Tax <- Tax[match(colnames(Tab),Tax$ID),]

Dat <- create_dataset(Tab = Tab %>% t,Map = Map,Tax = Tax)

#Create contaminants to remove later
contam_otus <- c(grep(pattern = "chloroplast", as.character(Dat$Tax$Taxonomy),ignore.case = TRUE),
                 grep(pattern = "mitochondri", as.character(Dat$Tax$Taxonomy),ignore.case = TRUE),
                 grep(pattern = "oomycete", as.character(Dat$Tax$Taxonomy),ignore.case = TRUE),
                 grep(pattern = "k__unknown", as.character(Dat$Tax$Taxonomy),ignore.case = TRUE),
                 grep(pattern = "p__Bacteria_unclassified", as.character(Dat$Tax$Taxonomy),ignore.case = TRUE),
                 grep(pattern = "k__Eukaryota", as.character(Dat$Tax$Taxonomy),ignore.case = TRUE),
                 grep(pattern = "k__Archaea", as.character(Dat$Tax$Taxonomy),ignore.case = TRUE)
)
contam_otus <- row.names(Dat$Tax)[contam_otus]


# Filter Dataset
Dat_filter <- remove_taxons(Dat = Dat, taxons = contam_otus)
Dat_filter <- clean(Dat = Dat_filter,verbose = TRUE)
contam_otus <- Dat$Tax[ contam_otus, ]

Dat <- Dat_filter
Dat$Map$Depth <- colSums(Dat$Tab)

Dat <- Dat %>% subset.Dataset(Genotype != "Blank",drop = T,clean = T) %>%
  subset.Dataset(Genotype != "BulkSoil",drop = T,clean = T)

#Remove clear fraction mislabels
toremove <- c("Plate4RootBarriersRep2E2","Plate4RootBarriersRep2C4",
  "Plate3RootBarriersRep1D11","Plate3RootBarriersRep1B9")
Dat <- remove_samples(Dat = Dat,samples = toremove,droplevels = T)
Dat <- clean(Dat) 
# 1000 reads
Dat_raw <- remove_samples(Dat,row.names(Dat$Map)[ Dat$Map$Depth < 1000 ])
Dat_raw <- clean(Dat_raw)

#Clean the Genotypes
Dat_raw$Map$Genotype <- Dat_raw$Map$Genotype %>% gsub(pattern = "\\:|\\-|\\.",replacement = "_") %>%
  factor %>% relevel(ref = "Col_0")

#Add the classification
df_class <- read.table(file = "../rawdata/df_genotypes_to_classification.csv",
                       header = T,sep = ",")

Dat_raw$Map$Group <- match(Dat_raw$Map$Genotype,df_class$Genotype) %>%
  df_class$Group[.]

#Create taxonomy data.frame
mdf <- Dat_raw$Tax
df_tax <- mdf$Taxonomy %>% as.character %>%
  strsplit(split  = "\\;") %>%
  unlist %>%  gsub(pattern = "[a-z]__",replacement = "",) %>%
  gsub(pattern = " ",replacement = "") %>%
  matrix(data = .,ncol = 7,byrow = T) %>% as.data.frame
colnames(df_tax) <-c("Root","Kingdom","Phylum",
                     "Class","Order","Family","Genus")
df_tax <- cbind(mdf$ID,df_tax)
colnames(df_tax)[1] <- "Id"


### Rarefaction
set.seed(130816)
Dat_rar <- rarefaction(x = Dat_raw,sample =1000)
Dat_rar <- clean(Dat_rar)


###Relative abundance 
Dat_rab <- Dat_rar
#Here we are scaling each column (sample) by the total number of reads in that sample (250 in this case)
Dat_rab$Tab <- scale(x = Dat_rab$Tab,center = F,scale = colSums(Dat_rab$Tab))



Dat_amplicon <- list(RawCounts=Dat_raw,
                     Rarefied=Dat_rar,
                     RelativeAbundance=Dat_rab,
                     df_tax = df_tax
                     )



filename <- "../cleandata/dat_censusmutants_16s.RDS"
saveRDS(object = Dat_amplicon,file = filename)
dev.off()
rm(list = ls())
gc()

