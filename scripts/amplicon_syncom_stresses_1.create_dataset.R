library(ohchibi)
library(paletteer)
set.seed(130816)
setwd('/home/isai/Documents/results/rootbarriersmicro/scripts/')
dir.create("../cleandata")
dir.create("../figures")

#Read Tab
Tab <- read.table(file = "../rawdata/dada2_res_structure_syncom_stresses.tsv",
                  header = T,sep = "\t",check.names = F,quote = "",row.names = 1)

#Read metadata
Map <- read.table(file = "../rawdata/MetadataStresses.csv",
                  header = T,sep = "\t",quote = "",comment.char = "")

Map$Plate <- Map$Plate %>% paste0("Plate",.)

#Clean Map
Map$Treatment <- Map$Treatment %>% 
  gsub(pattern = "inoculum",replacement = "Inoculum") %>%
  gsub(pattern = "blanck",replacement = "Blank") %>%
  gsub(pattern = "blank",replacement = "Blank") %>%
  gsub(pattern = "empty",replacement ="Empty")  %>%
  gsub(pattern = "Full ",replacement = "Full") %>%
  factor

Map$Bacteria <- Map$Bacteria %>%
  gsub(pattern = "inoculum",replacement = "Inoculum") %>%
  gsub(pattern = "blanck",replacement = "Blank") %>%
  gsub(pattern = "blank",replacement = "Blank") %>%
  gsub(pattern = "empty",replacement ="Empty") %>% 
  factor

Map$Plant <- Map$Plant %>%
  gsub(pattern = "inoculum",replacement = "Inoculum") %>%
  gsub(pattern = "blanck",replacement = "Blank") %>%
  gsub(pattern = "blank",replacement = "Blank") %>%
  gsub(pattern = "empty",replacement ="Empty") %>% 
  gsub(pattern = "P",replacement ="Plant") %>% 
  gsub(pattern = "NPlant",replacement ="NoPlant") %>% 
  factor

Map$Replica <- Map$Replica %>% 
  gsub(pattern = "inoculum",replacement = "Inoculum") %>%
  gsub(pattern = "blanck",replacement = "Blank") %>%
  gsub(pattern = "blank",replacement = "Blank") %>%
  gsub(pattern = "empty",replacement ="Empty") %>% 
  gsub(pattern = "1",replacement = "Rep1") %>%
  gsub(pattern = "2",replacement = "Rep2") %>% factor

Map$Fraction <- Map$Fraction %>%
  gsub(pattern = "inoculum",replacement = "Inoculum") %>%
  gsub(pattern = "agar",replacement = "Agar") %>%
  gsub(pattern = "blanck",replacement = "Blank") %>%
  gsub(pattern = "blank",replacement = "Blank") %>%
  gsub(pattern = "empty",replacement ="Empty") %>% 
  gsub(pattern = "root",replacement ="Root") %>% 
  gsub(pattern = "shoot",replacement ="Shoot") %>%
  factor 

Map$Experiment <- Map$Experiment %>% paste0("Experiment",.)

Map$Treatment <- Map$Treatment %>% gsub(pattern = " ",replacement = "") %>%
  factor %>% relevel( x = .,ref = "Full")

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

Map$DADA2_Header <- paste0(Map$Plate,"RootBarriersSC",Map$PosWell)


rownames(Map) <- Map$DADA2_Header



#Create dataset
usable_samples <- intersect(Map$DADA2_Header,rownames(Tab)) 
Map <- match(usable_samples,Map$DADA2_Header) %>%
  Map[.,] %>% droplevels
Tab <- match(usable_samples,rownames(Tab)) %>%
  Tab[.,]


#Read the taxonomy
Tax<-read.table(file = "../rawdata/seqsformothur.full_v132.wang_syncom_stresses.taxonomy",header = F,sep = "\t")

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

#Change name of treatment column to stress
colnames(Map)[6] <- "Stress"

#Create dataset
Dat <- create_dataset(Tab = Tab %>% t,Map = Map,Tax = Tax)


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

df_uc <- read.table(file = "../rawdata/dada2_res_structure.asv.uc",header = F,sep = "\t")
df_uc <- df_uc[,c(1,9,10)] %>% 
  subset(V1 == "H") %>% 
  droplevels
hit_twice <- which(df_uc$V10 %>% table == 2) %>% names
df_uc[df_uc$V10 %in% hit_twice,]
#If they hit twice is asvs that are clear artifacts  because they are smaller asvs of bigger asvs

toremove <- c("ASV2511","ASV3724","ASV3822")

hits_asvs <- df_uc$V9 %>% as.character

hits_asvs <- which(!(hits_asvs %in% toremove)) %>%
  hits_asvs[.]


df_counts$Usable[which(df_counts$Id %in% hits_asvs)] <- "Yes"

total_reads <- df_counts$value %>% sum
sum_asvs <- df_counts %>% subset(Usable == "Yes") %$% value %>%
  sum

cat("Percentage of reads mapped to isolates: ",(sum_asvs/total_reads)*100,"% \n")

#Evaluate for which isolates we dont have sequences yet
df_missing <- meta %>% subset(HitStress == "No") %>% droplevels
df_missing <- df_missing[,c(2,3,8,16)]

df_counts <- with(df_counts,order(Id)) %>% df_counts[.,]
df_counts_sub <- df_counts %>% subset(Usable == "No") %>% droplevels

meta %>% subset(HitStress == "No")

#Manually check df_missing and df_counts_sub with global alignment in emboss to assemble the following informaiton

#We are missing two pseudomonas
#ASV1 and ASV6 are two pseudomonas that  we may be missing

#Using closest relative to Pseudoomonassp 397MF that is
#2643221502_Pseudomonas395MF we have a perfect match of ASV1 and not ASV6
#therefore ASV1  corresponds to 2643221503  RMF397

#ASV6 does not match to  2556921015  RCL58 Pseudomonas umsongensis UNC430CL58Col
#and also none of the other pseudomonas match to it
#because we dont have reference of rcl58 i used the closest relative 2643221484 Pseudomonas 113MF
#so without too much confidence we can say asv6 should be rcl58 given the fact we only inputted two pseudomonas



#L61 is ASV26 or ASV38 perfect hit except for the N in the rnammer assembly
df_counts_sub %>% subset(Genus == "Duganella")
meta %>% subset(Genus == "Duganella")

#we only have 1 duganella in dataset and those asvs, 26 and 38 are the most abundant ones
#also both asvs appear again in the mutant dataset as asv35(asv36) and asv80 (asv38)
#i could not resolve using phylogenetics due to a lot of ns in the other references
#i assume the n exist cause probably duganella has more than one copy
temp <- which(rownames(Dat$Tab) %in% c("ASV26","ASV38")) %>% Dat$Tab[.,] %>% t
cor.test(temp[,1],temp[,2])

#Both of the  asvs have a correlation of 0.99 so they are covaring across the dataset
#confidently we can say they belong to l61



#L50 There are 3 Ns in the reference. When aligned in emboss all as vs have
#identical score meaning they are different between each other by the Ns
#With confidence we have the N problem with ASV 10, ASV 18 and ASV20.
#ASv 19 has more mismatches that could be due to another 16s copy hidden in the genome

df_counts_sub %>% subset(Genus == "Rahnella")
#which of these ones appear in the mutant dataset
#all asv10,18,20 and 19 appear in the stress dataset asasv 6,15,18 and 19
#lets check their correlation pattern
temp <- which(rownames(Dat$Tab) %in% c("ASV10","ASV18","ASV19","ASV20")) %>% Dat$Tab[.,] %>% t
cor.test(temp[,1],temp[,2])
cor.test(temp[,1],temp[,3])
cor.test(temp[,1],temp[,4])
cor.test(temp[,2],temp[,4])
cor.test(temp[,2],temp[,3])
cor.test(temp[,3],temp[,4])
#All are correlating extremely well so all of them will be confidently assigned as L50



#ASV792 is the only brevundimonas in the dataset could be L363 but a lot of mismatches with rnammer
meta %>% subset(Genus == "Brevundimonas")
#two bevundimonas but the other is a right hit
#amazingly asv792 could not be found in the stress dataset
#in the stress dataset we could find a hit for l363 it is 
#ASV376
#now lets do exploration with phylogeny l363 is an outgroup and the closest relative
#being 2643221663 brevundimonas root 1279 does not have a hit neither withasv 792 nor with asv376
#therefore we di dnot detect that brevundimonas in this dataset




#Sphingomonas L226 is ASV29 and ASV47 given the same N problem 
meta %>% subset(Genus == "Sphingomonas")
df_counts_sub %>% subset(Genus == "Sphingomonas")
#both asv29 and asv47 appear also in the stress dataset

df_counts %>% subset(Genus == "Sphingomonas")
#using phylogenetics and the closert relative to l226 that is 2643221763 sphingomonas leaf30
# we do have a perfect hit with asv29 this would make sense given the fact
#it is more abundant than asv47
#also in the stress dataset asv29 matches to asv51 that is more abundant than 
#asv47 homolog in that adataset that is asv97
#just for curiosity check correlation
temp <- which(rownames(Dat$Tab) %in% c("ASV29","ASV47")) %>% Dat$Tab[.,] %>% t
cor.test(temp[,1],temp[,2])
#sort of correlated but not as perfectly
#probably a subpopulation of sphingomonas nevertheless we done have evidence for that
#with this we conclude that asv29 corresponds to l226




#ASV34 and ASV 77 are paenibacillus and should correspond to RCL52 CL130 for which we dont have reference
meta %>% subset(HitStress == "No")
df_counts_sub %>% subset(Genus == "Paenibacillus")
#to start both asv34 and asv77 have hit in the stress dayaset
#asv35 corresponds to asv36 and asv 77 to asv101
#we dont have reference for both paenibacillus so we will use closest relative

#Using closest relative to Paeni cl52 that was paeni cl91 2643221535
#we got a hit with asv34
# and not a hit with asv 77

#now the problem is that cl130 is closeset relative to cl52 so basically
#this is a circular probelm and asv34 will correspond to both
#now lets see the correlation of asv 34 and asv 77
temp <- which(rownames(Dat$Tab) %in% c("ASV34","ASV77")) %>% Dat$Tab[.,] %>% t
cor.test(temp[,1],temp[,2])

#97 not 99 so also could be some sub population intra
#both of them assigned to asv34



#ASv16 is Paenarthrobacter RMF26 for wihch there was not reference
meta %>% subset(HitStress == "No")
df_counts_sub [,1:3]
# we found asv 16 also in stress dataset as highly abundant
#also using phylogenetics closest isolate 2523533508  Paenarthrobacter nicotinovorans 231Sha2.1M6
#does match perfectly with asv16 so confidently we can
#assign asv16 to rm26



#RMF279 is supposed to be a Caulobacter but it has high level of contamination an hetogeneity
#The dada2 analysis shows that there is ASV9 is an alphaproteobacteria like caulobacter and it is quite abundant
#It can be that ASV9 is actually that isolate for which the genome denotes it is a caulobacter but given the high
#level of contamination it could be other thing

df_counts_sub [3,]
#asv9 is found also in mutant dataset as asv8
meta %>% subset(HitStress == "No")
#The only explanation is that that caulobavter is something completely different than what we think

#check correlation with asv 30 that is rhizobium and mapped to the only rhizobium in database
temp <- which(rownames(Dat$Tab) %in% c("ASV9","ASV30")) %>% Dat$Tab[.,] %>% t
cor.test(temp[,1],temp[,2])
#completely independent so most parsimonious is that caulobacter is something else



certain <- c("ASV1","ASV26","ASV38","ASV10","ASV18","ASV19","ASV20","ASV29","ASV34","ASV16")
uncertain <- c("ASV9","ASV6")

df_counts$Usable[which(df_counts$Id %in% certain)] <- "Certain"
df_counts$Usable[which(df_counts$Id %in% uncertain)] <- "Uncertain"


#Create final structure to plot
good_seqs <- df_uc$V10 %>% as.character %>% unique
meta_good <- which(meta$USeq %in% good_seqs) %>% meta[.,]
meta_good <- meta_good[,c(1,2,3,16)] %>% unique %>% droplevels

#RMf 160 is repeated
meta_good %>% subset(Id == "RMF160")
meta_good %>% subset(USeq == "Sequence_9")
meta_good %>% subset(USeq == "Sequence_6")
#It makes sense given the fact we have two copies

#Append missing ifnromation of the new strains
#add equivalent in stress syncom
meta_good$ASVStress <- match(meta_good$USeq,df_uc$V10) %>% 
  df_uc$V9[.] %>% as.character


#Pseudoomonas
meta_good <- data.frame(match("2643221503",meta$taxon_oid) %>% meta[.,c(1,2,3)],
           USeq = "Sequence_33",ASVStress = "ASV1") %>%
  rbind(meta_good,.)




#Duganella
#L61 is ASV26 or ASV38

a <- data.frame(match("L61",meta$Id) %>% meta[.,c(1,2,3)],
                USeq= "Sequence_34",ASVStress = "ASV26")

b <-  data.frame(match("L61",meta$Id) %>% meta[.,c(1,2,3)],
                 USeq= "Sequence_35",ASVStress = "ASV38")        


meta_good <- rbind(meta_good,a,b)


#serratia L50 is c("ASV10","ASV18","ASV19","ASV20")
a <- data.frame(match("L50",meta$Id) %>% meta[.,c(1,2,3)],
                USeq= "Sequence_36",ASVStress = "ASV10")
b <- data.frame(match("L50",meta$Id) %>% meta[.,c(1,2,3)],
                USeq= "Sequence_37",ASVStress = "ASV18")
c <- data.frame(match("L50",meta$Id) %>% meta[.,c(1,2,3)],
                USeq= "Sequence_38",ASVStress = "ASV19")
d <- data.frame(match("L50",meta$Id) %>% meta[.,c(1,2,3)],
                USeq= "Sequence_39",ASVStress = "ASV20")


meta_good <- rbind(meta_good,a,b,c,d)

#Brevundimonas L363
a <- data.frame(match("L363",meta$Id) %>% meta[.,c(1,2,3,16)],ASVStress = NA)
meta_good <- rbind(meta_good,a)


#Sphingomonas asv29 corresponds to l226
a <- data.frame(match("L226",meta$Id) %>% meta[.,c(1,2,3)],
                USeq= "Sequence_40",ASVStress = "ASV29")

meta_good <- rbind(meta_good,a)

#Paeniabcillus asv34 RCL52 CL130
a <- data.frame(match("RCL52",meta$Id) %>% meta[.,c(1,2,3)],
                USeq= "Sequence_41",ASVStress = "ASV34")
b <- data.frame(match("RCL130",meta$Id) %>% meta[.,c(1,2,3)],
                USeq= "Sequence_41",ASVStress = "ASV34")

meta_good <- rbind(meta_good,a,b)

#Arthrobacter
#assign asv16 to rm26
a <- data.frame(match("RMF26",meta$Id) %>% meta[.,c(1,2,3)],
                USeq= "Sequence_42",ASVStress = "ASV16")

meta_good <- rbind(meta_good,a)

#Caulobacter RMF279 uncertain asv9
a <- data.frame(match("RMF279",meta$Id) %>% meta[.,c(1,2,3)],
                USeq= "Sequence_43",ASVStress = "ASV9")

meta_good <- rbind(meta_good,a)


#Pseudomonas RCL58 weird too but matched to asv6
a <- data.frame(match("RCL58",meta$Id) %>% meta[.,c(1,2,3)],
                USeq= "Sequence_44",ASVStress = "ASV6")

meta_good <- rbind(meta_good,a)

#Add certainty
meta_good$Certainty <- "Certain"
meta_good$Certainty[which(meta_good$Id %in% c("RCL58","RMF279"))] <- "Uncertain"


sum_asvs <- df_counts %>% subset(Usable == "Yes" | Usable == "Certain" | Usable == "Uncertain" ) %$% value %>%
  sum
cat("Percentage of reads mapped to isolates: ",(sum_asvs/total_reads)*100,"% \n")

usable_asvs <- df_counts %>% subset(Usable != "No") %>% droplevels %$% Id %>% as.character

#Readjust files for mapping future syncom experiments
meta_good$USeq %>% as.character %>% gsub(pattern = "Sequence_",replacement = "") %>%
  as.numeric %>% unique %>% sort

meta_good$USeq %>% as.character %>% gsub(pattern = "Sequence_",replacement = "") %>%
  as.numeric %>% unique %>% length

df_tax <- match(meta_good$taxon_oid,meta$taxon_oid) %>%
  meta[.,c(4,5,6,7,8)]

meta_good <- cbind(meta_good,df_tax)
rownames(meta_good) <- NULL

meta_good$Taxonomy <- paste0("Root; k__Bacteria; p__",meta_good$Phylum,"; c__",meta_good$Class,
       "; o__",meta_good$Order,"; f__",meta_good$Family,"; g__",meta_good$Genus)
  

meta_good <- meta_good %>% droplevels

#Here we need to transform the matrix into a new matrix in which we summarize it
chosen_asvs <- na.omit(meta_good$ASVStress) %>% as.character %>% unique

Tab <- Dat$Tab
Tab_chosen <- match(chosen_asvs,rownames(Tab)) %>%
  Tab[.,]
sum(Tab_chosen)/sum(Tab)

Tab_chosen <- Tab_chosen[,which(colSums(Tab_chosen) !=0)]

#Put them in the sequence context

rownames(Tab_chosen) <- match(rownames(Tab_chosen),meta_good$ASVStress) %>%
  meta_good$USeq[.] %>% as.character


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
  mv <- colSums(mTab)  %>% matrix(data = .,nrow = 1,byrow = T)
  rownames(mv) <- strain
  colnames(mv) <- colnames(mTab)
  Tab_to_add <- rbind(Tab_to_add,mv)
  mtax <- data.frame(ID = strain,Taxonomy = meta_good %>% subset(Id == strain) %$% Taxonomy %>% as.character %>% unique)
  df_tax <- rbind(df_tax,mtax)
  
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

paleta_stress <- paletteer_d(package = "ochRe",palette = "namatjira_qual" )[1:8]
names(paleta_stress) <- levels(Dat_raw$Map$Stress)[1:8]


Dat_amplicon <- list(RawCounts=Dat_raw,
                     Rarefied=Dat_rar,
                     RelativeAbundance=Dat_rab,
                     meta = meta_good,
                     paleta_stress = paleta_stress)

filename <- "../cleandata/dat_rootbarriers_amplicon_bacteria_syncom_stresses.RDS"
saveRDS(object = Dat_amplicon,file = filename)

