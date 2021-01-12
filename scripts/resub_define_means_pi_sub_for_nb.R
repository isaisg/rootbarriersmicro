library(ohchibi)

set.seed(130816)
setwd('/Users/isaisalasgonzalez/Documents/rootbarriersmicro/scripts/')

#Read the screnning barrier dataset
Dat <- readRDS(file = "../cleandata/dat_screeningbarriers.RDS")

#Suberin
df <- Dat$Suberin
df <- dcast(data = df,formula = Strains + Plant + Batch ~Localization,value.var = "Normvalue_ra")
df$ZoneSum <- df$No_expression + df$Discrete_expression

df_sum_sub <- aggregate(ZoneSum~Strains,df,mean)
df_sum_sub %>% subset(Strains == "NB")

#Totla number of cells
df <- Dat$Suberin
df %>% subset(Localization == "Continous_expression") %>%
  aggregate(Total_cells_Normvalue~Strains,.,mean) %>%
  subset(Strains == "NB")

#Propidium

df <- Dat$Propidium
df %>% 
  aggregate(Normvalue~Strains,.,mean) %>%
  subset(Strains == "NB")

