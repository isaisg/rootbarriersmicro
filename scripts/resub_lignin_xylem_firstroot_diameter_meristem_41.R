library(ohchibi)
library(multcomp)
library(Rmisc)
library(ggrepel)
library(scales)
library(ggtree)


#Analysis of development
set.seed(130816)
setwd('/Users/isaisalasgonzalez/Documents/rootbarriersmicro/scripts/')

df <- read.table(file = "../rawdata/growthparametersmono.csv",header = T,sep = ",")
colnames(df) <- c("Date","Strain","Replicate","MeristemSize",
                  "DistanceCSLignin","DistanceXylem","DistanceFirstRootHair",
                  "RootDiameter","NumCellsCSLigning","NumCellsXylem","NumCellsFirstRootHair")

df$Date <- df$Date %>% gsub(pattern = '\\"',replacement = "") %>%
  gsub(pattern = " ",replacement = "") %>% factor
df$Strain <- df$Strain %>% factor %>% relevel(ref = "NB")

#Check consitency of names
df_names <- read.table(file = "../rawdata/41_strains_selected.txt") %$% V1
df_names <- data.frame(Name = df_names,Id = df_names %>% 
                         gsub(pattern = "^L|^R",replacement = ""))
df_names <- df_names %>%
  rbind(data.frame(Name = "NB", Id = "NB"))

mnames <- df$Strain %>% as.character %>% unique

#Check which names we cover
which(!(df_names$Id %in% mnames)) %>%
  df_names[.,]


which(!(mnames %in% df_names$Id)) %>%
  mnames[.]

df$Strain <- df$Strain %>%
  gsub(pattern = "CL50$",replacement = "50") %>%
  gsub(pattern = "^333$",replacement = "335") %>%
  gsub(pattern = "^323$",replacement = "363")

df$Strain <- match(df$Strain,df_names$Id) %>%
  df_names$Name[.]

df$Strain <- df$Strain %>% factor %>% relevel(ref = "NB")


df$Id <- 1:nrow(df)
Map <- df[,c("Date","Strain","Replicate","Id")]
rownames(Map) <- df$Id

Tab <- df[,c("MeristemSize","DistanceCSLignin","DistanceXylem","DistanceFirstRootHair","RootDiameter","NumCellsCSLigning","NumCellsXylem","NumCellsFirstRootHair")]
rownames(Tab) <- df$Id


Dat <- create_dataset(Tab = Tab %>% t,Map = Map)

Tab_z <- Dat$Tab %>% t %>% scale %>% t
Dat_z <- create_dataset(Tab = Tab_z,Map = Map)

#Save both structures
mlist <- list(Dat = Dat,
      Dat_z = Dat_z)
saveRDS(object = mlist,file = "../cleandata/dat_lignin_xylem_firstroot_diameter_meristem_41.RDS")
rm(list=ls())
gc()