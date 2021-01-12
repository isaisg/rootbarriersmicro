library(ohchibi)
library(ggpubr)
library(Rmisc)
library(ggrepel)
library(scales)
#
set.seed(130816)
setwd('/Users/isaisalasgonzalez/Documents/rootbarriersmicro/scripts/')

### Read the table
df <- read.table(file = "../rawdata/revision_root_length_vs_number_of_cells.csv",header = T,sep = ",")

## Split the dataset by the variable measured
df_root <- df[,c("Strain","Root_length","Batch")]
df_pi <- df[,c("Strain","PI_permeability","Batch")]
df_tc <- df[,c("Strain","Total_number_of_cells","Batch")]

### Calculate normalization factors using the NB present in both batches

### Root ####
overall_mean <- df_root %>% subset(Strain == "NB") %$% Root_length %>% mean
df_batch <- df_root %>% subset(Strain == "NB") %>%
  aggregate(Root_length~Batch,.,mean)
df_batch$Nfactor <- overall_mean/df_batch$Root_length

df_root <- merge(df_root,df_batch[,c("Batch","Nfactor")], by = "Batch")
df_root$Norm_Root_length <- df_root$Root_length * df_root$Nfactor

### Pi ####
df_pi <- df_pi %>% na.omit
overall_mean <- df_pi %>% subset(Strain == "NB") %$% PI_permeability %>% mean
df_batch <- df_pi %>% subset(Strain == "NB") %>%
  aggregate(PI_permeability~Batch,.,mean)
df_batch$Nfactor <- overall_mean/df_batch$PI_permeability

df_pi <- merge(df_pi,df_batch[,c("Batch","Nfactor")], by = "Batch")
df_pi$Norm_PI_permeability <- df_pi$PI_permeability * df_pi$Nfactor

### Cells  ####
df_tc <- df_tc %>% na.omit
overall_mean <- df_tc %>% subset(Strain == "NB") %$% Total_number_of_cells %>% mean
df_batch <- df_tc %>% subset(Strain == "NB") %>%
  aggregate(Total_number_of_cells~Batch,.,mean)
df_batch$Nfactor <- overall_mean/df_batch$Total_number_of_cells

df_tc <- merge(df_tc,df_batch[,c("Batch","Nfactor")], by = "Batch")
df_tc$Norm_Total_number_of_cells <- df_tc$Total_number_of_cells * df_tc$Nfactor


#Calculate number of samples per strain per batch
df_temp <- df_tc
df_temp$Count <- 1

df_num <- dcast(data = df_temp,formula = Strain~Batch,value.var = "Count",fun.aggregate = sum)

#### Loop over the primary root elongation and pica number of roots equal to the number picked
mstrains <- df_num$Strain %>% unique
df_rl <- NULL
for(st in mstrains){
  df_numtemp <- df_num %>% subset(Strain == st) %>% droplevels
  df_temp <- df_root %>% subset(Strain == st) %>% droplevels
  #Determine how many batches we have
  mbatches <- df_temp$Batch %>% unique
  for(mb in mbatches){
    df_in <- df_temp %>% subset(Batch == mb)
    ntochoose <- df_numtemp[,mb]
    df_in <- with(df_in,order(-Norm_Root_length)) %>%
      df_in[.,]
    ml <- df_in$Norm_Root_length[1:ntochoose]
    df_rl <- data.frame(Strain = st,Norm_Root_length = ml) %>%
      rbind(df_rl ,.)
  }
}



#Create summary data frames 
df_sum_rl <- summarySE(data = df_rl,measurevar = "Norm_Root_length",groupvars = "Strain")
df_sum_pi <- summarySE(data = df_pi,measurevar = "Norm_PI_permeability",groupvars = "Strain")
df_sum_tc <- summarySE(data = df_tc,measurevar = "Norm_Total_number_of_cells",groupvars = "Strain")

##Create merge structure for plotting
merged <- merge(df_sum_rl[,c("Strain","Norm_Root_length","ci")],df_sum_tc[,c("Strain","Norm_Total_number_of_cells","ci")],
                by  = "Strain") %>%
  merge(.,df_sum_pi[,c("Strain","Norm_PI_permeability","ci")], by = "Strain")


##Create a summary using all primary root measured
df_raw_rl <- summarySE(data = df_root,measurevar = "Norm_Root_length",groupvars = "Strain")
colnames(df_raw_rl)[3] <- "Norm_Root_length_All"
colnames(df_raw_rl)[6] <- "ci.Norm_Root_length_All"

merged <- merge(merged,df_raw_rl[,c("Strain","Norm_Root_length_All","ci.Norm_Root_length_All")], by = "Strain")

##### Rename columns in the structure for plotting 
colnames(merged) <- c("Strain","Norm_Root_length","ci.Norm_Root_length",
                      "Norm_Total_number_of_cells","ci.Norm_Total_number_of_cells",
                      "Norm_PI_permeability","ci.Norm_PI_permeability","Norm_Root_length_All","ci.Norm_Root_length_All")


#Adjust the names of the 41 
df_41 <- read.table(file = "../rawdata/41_isolates.tsv")
noms <- c(df_41$V2,"NB")
merged$Names <- noms
merged$Names[26] <- "RCL130"
merged$Names[27] <- "RCL4"
merged[,c("Strain","Names")]

#Append the names
df_rl$Names <- match(df_rl$Strain,merged$Strain) %>%
  merged$Names[.]
df_pi$Names <- match(df_pi$Strain,merged$Strain) %>%
  merged$Names[.]
df_tc$Names <- match(df_tc$Strain,merged$Strain) %>%
  merged$Names[.]
df_raw_rl$Names <- match(df_raw_rl$Strain,merged$Strain) %>%
  merged$Names[.]

#Save structure
mlist <- list(
  df_rl = df_rl,
  df_pi = df_pi,
  df_tc = df_tc,
  df_raw_rl = df_raw_rl,
  df_sum = merged
)
saveRDS(object = mlist,file = "../cleandata/dat_root_length_vs_number_of_cells_41.RDS")
rm(list=ls())
gc()
