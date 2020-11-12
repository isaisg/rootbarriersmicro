library(ohchibi)
library(paletteer)
set.seed(130816)
setwd('/home/isai/Documents/results/rootbarriersmicro/scripts//')
dir.create("../cleandata")
dir.create("../figures")

#Read Tab 
Tab <- read.table(file = "../rawdata/dada2_res_structure_syncom_mutants.tsv",
                  header = T,sep = "\t",check.names = F,quote = "",row.names = 1)

#Read metadata
Map <- read.table(file = "../rawdata/metadata_syncom_mutants.csv",
                  header = T,sep = "\t",quote = "",comment.char = "")

Map$Plate <- Map$Plate %>% paste0("Plate",.)


Map$Replica <- Map$Replica %>% 
  gsub(pattern = "inoculum",replacement = "Inoculum") %>%
  gsub(pattern = "blanck",replacement = "Blank") %>%
  gsub(pattern = "blank",replacement = "Blank") %>%
  gsub(pattern = "empty",replacement ="Empty") %>% 
  gsub(pattern = "1",replacement = "Rep1") %>%
  gsub(pattern = "2",replacement = "Rep2") %>%
  gsub(pattern = "3",replacement = "Rep3") %>%
  gsub(pattern = "4",replacement = "Rep4") %>%
  
  factor

Map$Fraction <- Map$Fraction %>%
  gsub(pattern = "inoculum",replacement = "Inoculum") %>%
  gsub(pattern = "agar",replacement = "Agar") %>%
  gsub(pattern = "blanck",replacement = "Blank") %>%
  gsub(pattern = "blank",replacement = "Blank") %>%
  gsub(pattern = "empty",replacement ="Empty") %>% 
  gsub(pattern = "root",replacement ="Root") %>% 
  gsub(pattern = "shoot",replacement ="Shoot") %>%
  factor 



Map$PosWell <- paste0(Map$PosWell,Map$PosNum)
Map <- with(Map,order(Order)) %>%
  Map[.,]


#Read equivalencies
il <- read.table(file = "../rawdata/illuminaindex_to_samplewell.tsv.csv",header = T,sep = "\t")
fin <- read.table(file = "../rawdata/equivalencies.txt",header = F,sep = "\t")
colnames(fin) <- c("Num","Label")
il <- merge(il,fin, by = "Num")
eq_gab <- Map[1:96,c(1,3)] %>% droplevels
colnames(eq_gab)[1] <- c("Num")
merged <- merge(il,eq_gab,by = "Num")

Map$PosWell <- match(Map$PosWell,merged$PosWell) %>%
  merged$Sample_Well[.]

Map$DADA2_Header <- paste0(Map$Plate,"RootBarriersSCMut",Map$PosWell)


rownames(Map) <- Map$DADA2_Header

Map$Genotype <- Map$Genotype %>% gsub(pattern = " ",replacement = "") %>%
  factor

#Create dataset
usable_samples <- intersect(Map$DADA2_Header,rownames(Tab)) 
Map <- match(usable_samples,Map$DADA2_Header) %>%
  Map[.,] %>% droplevels
Tab <- match(usable_samples,rownames(Tab)) %>%
  Tab[.,]


#Read the taxonomy
Tax<-read.table(file = "../rawdata/seqsformothur.full_v132.wang.mutants.taxonomy",header = F,sep = "\t")

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



#Create dataset
Dat <- create_dataset(Tab = Tab %>% t,Map = Map,Tax = Tax)

Dat <- Dat %>% subset.Dataset(Genotype != "myb41",drop = T,clean = T)

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

#Create taxonomy data.frame
mdf <- Dat$Tax
df_tax <- mdf$Taxonomy %>% as.character %>%
  strsplit(split  = "\\;") %>%
  unlist %>%  gsub(pattern = "[a-z]__",replacement = "",) %>%
  gsub(pattern = " ",replacement = "") %>%
  matrix(data = .,ncol = 7,byrow = T) %>% as.data.frame
colnames(df_tax) <-c("Root","Kingdom","Phylum",
                     "Class","Order","Family","Genus")
df_tax <- cbind(mdf$ID,df_tax)
colnames(df_tax)[1] <- "Id"


#Check the proportion of ASVs captured by the usable asvs
df_counts <- Dat$Tab %>% melt %>%
  aggregate(data = .,value ~Var1,FUN = "sum")
df_counts$Usable <- "No"


#Merge df_coutns with df_tax
colnames(df_counts) <- c("Id","value","Usable")
df_counts <- merge(df_counts,df_tax, all.x = TRUE)
df_counts <- df_counts[,c(1,2,3,10)]


meta <- read.table(file = "../rawdata/metadata_41_isolates_syncom_rootbarriers.tsv",
                   header = T,sep = "\t",quote = "",comment.char = "",check.names = F)

## read the mapping tabe

df_uc <- read.table(file = "../rawdata/dada2_res_structure.asv.mutants.uc",header = F,sep = "\t")
df_uc <- df_uc[,c(1,9,10)] %>% 
  subset(V1 == "H") %>% 
  droplevels
hit_twice <- which(df_uc$V10 %>% table == 2) %>% names
df_uc[df_uc$V10 %in% hit_twice,]


hits_asvs <- df_uc$V9 %>% as.character



df_counts$Usable[which(df_counts$Id %in% hits_asvs)] <- "Yes"

total_reads <- df_counts$value %>% sum
sum_asvs <- df_counts %>% subset(Usable == "Yes") %$% value %>%
  sum

cat("Percentage of reads mapped to isolates: ",(sum_asvs/total_reads)*100,"% \n")

######### Part of transformation ##############

#Here we need to transform the matrix into a new matrix in which we summarize it
chosen_asvs <- df_counts %>% subset(Usable == "Yes") %$% Id %>% as.character %>%
  unique

Tab <- Dat$Tab
Tab_chosen <- match(chosen_asvs,rownames(Tab)) %>%
  Tab[.,]
sum(Tab_chosen)/sum(Tab)

Tab_chosen <- Tab_chosen[,which(colSums(Tab_chosen) !=0)]

#Put them in the sequence context
#read stres meta
Dat_st <-readRDS(file = "../cleandata/dat_rootbarriers_amplicon_bacteria_syncom_stresses.RDS")
meta_good <- Dat_st$meta

#here map to useq using uc
rownames(Tab_chosen) <- match(rownames(Tab_chosen),df_uc$V9) %>%
  df_uc$V10[.] %>% as.character



#Check the ones with one strain covering more than one pattern
tocollapse <- names(which(meta_good$Id %>% table >1))

tosub <- which(meta_good$Id %in% tocollapse) %>% meta_good$USeq[.] %>%
  as.character %>% unique


Tab_other <- which(!(rownames(Tab_chosen) %in% tosub)) %>%
  Tab_chosen[.,]


Tab_to_add <- NULL
df_tax <- NULL

for(strain in tocollapse){
  mtc <- meta_good %>% subset(Id == strain) %$% USeq %>% as.character %>% unique
  mTab <- which(rownames(Tab_chosen) %in% mtc) %>% Tab_chosen[.,]
  if(is.null(ncol(mTab))){
    mv <- mTab %>% as.matrix %>% t
    rownames(mv) <- strain
    Tab_to_add <- rbind(Tab_to_add,mv)
    mtax <- data.frame(ID = strain,Taxonomy = meta_good %>% subset(Id == strain) %$% Taxonomy %>% as.character %>% unique)
    df_tax <- rbind(df_tax,mtax)
  }else{
    mv <- colSums(mTab)  %>% matrix(data = .,nrow = 1,byrow = T)
    rownames(mv) <- strain
    colnames(mv) <- colnames(mTab)
    Tab_to_add <- rbind(Tab_to_add,mv)
    mtax <- data.frame(ID = strain,Taxonomy = meta_good %>% subset(Id == strain) %$% Taxonomy %>% as.character %>% unique)
    df_tax <- rbind(df_tax,mtax)
  }

  
}

#Now process the Tab_other to add the correct name
new_names <- NULL
for(i in 1:nrow(Tab_other)){
  mnom <- rownames(Tab_other)[i]
  nnom <- meta_good %>% subset(USeq == mnom) %$% Id %>% 
    as.character %>% paste(collapse = "_")
  new_names <- c(new_names,nnom)
  mtax <- data.frame(ID = nnom,Taxonomy = meta_good %>% subset(USeq == mnom) %$% Taxonomy %>% 
                       as.character %>% unique)
  df_tax <- rbind(df_tax,mtax)
}


rownames(Tab_other) <- new_names


Tab <- rbind(Tab_other,Tab_to_add)
Tax <- df_tax
Map <-Dat$Map

#Put everything in the same context and create final dataset
Tax <- match(rownames(Tab),Tax$ID) %>%
  Tax[.,]
rownames(Tax) <- Tax$ID

Map <- match(colnames(Tab),Map$DADA2_Header) %>%
  Map[.,]

Dat <- create_dataset(Tab = Tab,Map = Map,Tax = Tax)

#Count depth
Dat$Map$Depth <- colSums(Dat$Tab)
Dat$Map$Depth %>% sort
# 1000 reads
Dat_raw <- remove_samples(Dat,row.names(Dat$Map)[ Dat$Map$Depth < 1000 ])
Dat_raw <- clean(Dat_raw)


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
                     meta = meta_good)

filename <- "../cleandata/dat_rootbarriers_amplicon_bacteria_syncom_mutants.RDS"
saveRDS(object = Dat_amplicon,file = filename)
rm(list=ls())
