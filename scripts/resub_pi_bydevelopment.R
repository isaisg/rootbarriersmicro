library(ohchibi)
library(Rmisc)
library(ggpubr)
library(scales)

setwd('/Users/isaisalasgonzalez/Documents/rootbarriersmicro/scripts/')

set.seed(130816)
Dat <- readRDS(file = "../cleandata/dat_screeningbarriers.RDS")

# Propidium 
Tab_prop <- Dat$Propidium 

### Do the plotting
df_pi <- Tab_prop%>%
  summarySE(data = .,measurevar = "Normvalue_ra",groupvars = "Strains")
df_numcells <- Tab_prop%>%
  summarySE(data = .,measurevar = "Total_cells_Normvalue",groupvars = "Strains")

merged <- merge(df_pi[,c("Strains","Normvalue_ra","ci")],df_numcells[,c("Strains","Total_cells_Normvalue","ci")], by = "Strains")

#Read the  information about significance
res <- read.table(file = "../../rootbarriersmicro/cleandata/res_merged_propidium_suberin.tsv",header = T,sep = "\t")
res$PISignificance <- "No"
res$PISignificance[which(res$Propidium_p.adj < 0.1)] <- "Yes"

merged$PISignificance <- match(merged$Strains,res$Strains) %>%
  res$PISignificance[.]


y_axis_nb <- merged %>% subset(Strains == "NB") %$% Normvalue_ra
x_axis_nb <- merged %>% subset(Strains == "NB") %$% Total_cells_Normvalue

### All points 
p <- ggplot(data = merged %>% subset(Strains != "NB"),aes(Total_cells_Normvalue,Normvalue_ra)) +
  geom_vline(xintercept = x_axis_nb,size =0.5,linetype = "solid",color = "black") +
  geom_hline(yintercept = y_axis_nb,size =0.5,linetype = "solid",color = "black") +
  #geom_smooth(method = "lm",se = F,size = 1,color = "black") + 
  #geom_linerange(aes(ymin = Normvalue - ci.x ,ymax = Normvalue + ci.x),color = "#D9D9D9",size = 0.5) +
  #geom_linerange(aes(xmin = Total_cells_Normvalue -ci.y,xmax = Total_cells_Normvalue + ci.y),color = "#D9D9D9",size = 0.5) +
  geom_point(aes(fill = PISignificance),shape =21,size = 4) + 
  theme_ohchibi(font_family = "Helvetica")  +
  xlab(label = "Normalized total number of cells") +
  ylab(label = "Propidium iodide") +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_line(color = "#D9D9D9",size = 0.1),
    panel.grid.major.x = element_line(color = "#D9D9D9",size = 0.1)
    
  ) +
  #scale_y_continuous(limits = c(0,11),oob = rescale_none,breaks = seq(0,11,1)) +
  scale_x_continuous(breaks = seq(20,100,10))  +
  scale_fill_manual(values = c("#D9D9D9","red")) +
  stat_cor(label.sep = "\n",label.x.npc = 0.8,label.y.npc = 0.8,size =5)


###Only chosen points
p1 <- ggplot(data = merged %>% subset(Strains != "NB") %>%
               subset(PISignificance == "Yes"),aes(Total_cells_Normvalue,Normvalue_ra)) +
  geom_vline(xintercept = x_axis_nb,size =0.5,linetype = "solid",color = "black") +
  geom_hline(yintercept = y_axis_nb,size =0.5,linetype = "solid",color = "black") +
  #geom_smooth(method = "lm",se = F,size = 1,color = "black") + 
  #geom_linerange(aes(ymin = Normvalue - ci.x ,ymax = Normvalue + ci.x),color = "#D9D9D9",size = 0.5) +
  #geom_linerange(aes(xmin = Total_cells_Normvalue -ci.y,xmax = Total_cells_Normvalue + ci.y),color = "#D9D9D9",size = 0.5) +
  geom_point(aes(fill = PISignificance),shape =21,size = 4) + 
  theme_ohchibi(font_family = "Helvetica")  +
  xlab(label = "Normalized total number of cells") +
  ylab(label = "Propidium iodide") +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_line(color = "#D9D9D9",size = 0.1),
    panel.grid.major.x = element_line(color = "#D9D9D9",size = 0.1)
    
  ) +
 # scale_y_continuous(limits = c(0,11),oob = rescale_none,breaks = seq(0,11,1)) +
  scale_x_continuous(breaks = seq(20,100,10))  +
  scale_fill_manual(values = c("red")) +
  stat_cor(label.sep = "\n",label.x.npc = 0.8,label.y.npc = 0.8,size =5)

merged_pi <- merged

### Suberin
Tab_norm_sub <- Dat$Suberin

df_sub <- Tab_norm_sub %>% subset(Localization == "Continous_expression") %>% droplevels
df_numcell <- summarySE(data = df_sub,measurevar = "Total_cells_Normvalue",groupvars = "Strains")
df_cont <- summarySE(data = df_sub,measurevar = "Normvalue_ra",groupvars = "Strains")

merged <- merge(df_numcell[,c("Strains","Total_cells_Normvalue","ci")],
                df_cont[,c("Strains","Normvalue_ra","ci")], by = "Strains")


merged$SumZones <- 1-merged$Normvalue_ra
res$SubSignificance <- "No"
res$SubSignificance[which(res$Suberin_p.adj  < 0.01)] <- "Yes"

merged$SubSignificance <- match(merged$Strains,res$Strains) %>%
  res$SubSignificance[.]

y_axis_nb <- merged %>% subset(Strains == "NB") %$% SumZones
x_axis_nb <- merged %>% subset(Strains == "NB") %$% Total_cells_Normvalue


p2 <- ggplot(data = merged %>% subset(Strains != "NB"),aes(Total_cells_Normvalue,SumZones)) +
  geom_vline(xintercept = x_axis_nb,size =0.5,linetype = "solid",color = "black") +
  geom_hline(yintercept = y_axis_nb,size =0.5,linetype = "solid",color = "black") +
  #geom_smooth(method = "lm",se = F,size = 1,color = "black") + 
  #geom_linerange(aes(ymin = Normvalue - ci.x ,ymax = Normvalue + ci.x),color = "#D9D9D9",size = 0.5) +
  #geom_linerange(aes(xmin = Total_cells_Normvalue -ci.y,xmax = Total_cells_Normvalue + ci.y),color = "#D9D9D9",size = 0.5) +
  geom_point(aes(fill = SubSignificance),shape =21,size = 4) + 
  theme_ohchibi(font_family = "Helvetica")  +
  xlab(label = "Normalized total number of cells") +
  ylab(label = "Sum Zones RA") +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_line(color = "#D9D9D9",size = 0.1),
    panel.grid.major.x = element_line(color = "#D9D9D9",size = 0.1)
    
  ) +
  scale_y_continuous(limits = c(0,1),oob = rescale_none,breaks = seq(0,1,0.1)) +
  scale_x_continuous(breaks = seq(20,100,10))  +
  scale_fill_manual(values = c("#D9D9D9","red")) +
  stat_cor(label.sep = "\n",label.x.npc = 0.8,label.y.npc = 0.8,size =5)


p3 <- ggplot(data = merged %>% subset(Strains != "NB") %>%
               subset(SubSignificance == "Yes"),aes(Total_cells_Normvalue,SumZones)) +
  geom_vline(xintercept = x_axis_nb,size =0.5,linetype = "solid",color = "black") +
  geom_hline(yintercept = y_axis_nb,size =0.5,linetype = "solid",color = "black") +
  geom_smooth(method = "lm",se = F,size = 1,color = "black") + 
  #geom_linerange(aes(ymin = Normvalue - ci.x ,ymax = Normvalue + ci.x),color = "#D9D9D9",size = 0.5) +
  #geom_linerange(aes(xmin = Total_cells_Normvalue -ci.y,xmax = Total_cells_Normvalue + ci.y),color = "#D9D9D9",size = 0.5) +
  geom_point(aes(fill = SubSignificance),shape =21,size = 4) + 
  theme_ohchibi(font_family = "Helvetica")  +
  xlab(label = "Normalized total number of cells") +
  ylab(label = "Sum Zones RA") +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_line(color = "#D9D9D9",size = 0.1),
    panel.grid.major.x = element_line(color = "#D9D9D9",size = 0.1)
    
  ) +
  scale_y_continuous(limits = c(0,1),oob = rescale_none,breaks = seq(0,1,0.1)) +
  scale_x_continuous(breaks = seq(20,100,10))  +
  scale_fill_manual(values = c("red")) +
  stat_cor(label.sep = "\n",label.x.npc = 0.8,label.y.npc = 0.8,size =5)

merged_sub <- merged
rm(merged)


##Do visualization by position 
merged_pi$Total_cells_Normvalue %>% hist
merged_pi$RangeCels <- "20-29"
merged_pi$RangeCels[(merged_pi$Total_cells_Normvalue >= 30) & (merged_pi$Total_cells_Normvalue <40)] <- "30-39"
merged_pi$RangeCels[(merged_pi$Total_cells_Normvalue >= 40) & (merged_pi$Total_cells_Normvalue <50)] <- "40-49"
merged_pi$RangeCels[(merged_pi$Total_cells_Normvalue >= 50) & (merged_pi$Total_cells_Normvalue <60)] <- "50-59"
merged_pi$RangeCels[(merged_pi$Total_cells_Normvalue >= 60) & (merged_pi$Total_cells_Normvalue <70)] <- "60-69"
merged_pi$RangeCels[(merged_pi$Total_cells_Normvalue >= 70) & (merged_pi$Total_cells_Normvalue <80)] <- "70-79"
merged_pi$RangeCels[(merged_pi$Total_cells_Normvalue >= 80) & (merged_pi$Total_cells_Normvalue <90)] <- "80-89"
merged_pi$RangeCels[(merged_pi$Total_cells_Normvalue >= 90)]  <- ">=90"

merged_pi$RangeCels <- merged_pi$RangeCels %>% factor(levels = c("20-29","30-39","40-49","50-59","60-69","70-79","80-89",">=90"))

p_pi <- ggplot(data = merged_pi,aes(RangeCels,Normvalue_ra)) +
  geom_point(size = 1,alpha = 0.5) +
  stat_summary(fun.data =  mean_cl_normal,geom = "linerange",size = 0.5,color = "Red")+
  stat_summary(fun =  mean,geom = "point",size = 2,color = "Red")+
  scale_y_continuous(limits = c(0,0.25),breaks = seq(0,0.3,0.05)) +
  theme_ohchibi(font_family = "Helvetica")  +
  theme(
    axis.text.x = element_text(family = "Helvetica",angle =45,vjust = 1,hjust = 1)
  ) +
  xlab(label = "Range total number of cells") + ylab(label = "Propidium Iodide")


##Repeat for suberin

merged_sub$Total_cells_Normvalue %>% hist
merged_sub$RangeCels <- "20-29"
merged_sub$RangeCels[(merged_sub$Total_cells_Normvalue >= 30) & (merged_sub$Total_cells_Normvalue <40)] <- "30-39"
merged_sub$RangeCels[(merged_sub$Total_cells_Normvalue >= 40) & (merged_sub$Total_cells_Normvalue <50)] <- "40-49"
merged_sub$RangeCels[(merged_sub$Total_cells_Normvalue >= 50) & (merged_sub$Total_cells_Normvalue <60)] <- "50-59"
merged_sub$RangeCels[(merged_sub$Total_cells_Normvalue >= 60) & (merged_sub$Total_cells_Normvalue <70)] <- "60-69"
merged_sub$RangeCels[(merged_sub$Total_cells_Normvalue >= 70) & (merged_sub$Total_cells_Normvalue <80)] <- "70-79"
merged_sub$RangeCels[(merged_sub$Total_cells_Normvalue >= 80) & (merged_sub$Total_cells_Normvalue <90)] <- "80-89"
merged_sub$RangeCels[(merged_sub$Total_cells_Normvalue >= 90)]  <- ">=90"

merged_sub$RangeCels <- merged_sub$RangeCels %>% factor(levels = c("20-29","30-39","40-49","50-59","60-69","70-79","80-89",">=90"))

p_sub <- ggplot(data = merged_sub,aes(RangeCels,SumZones)) +
  geom_point(size = 1,alpha = 0.5) +
  stat_summary(fun.data =  mean_cl_normal,geom = "linerange",size = 0.5,color = "Red")+
  stat_summary(fun =  mean,geom = "point",size = 2,color = "Red")+
  scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.1)) +
  theme_ohchibi(font_family = "Helvetica")  +
  theme(
    axis.text.x = element_text(family = "Helvetica",angle =45,vjust = 1,hjust = 1)
  ) +
  xlab(label = "Range total number of cells") + ylab(label = "Sum NoExp + Discr Zones RA")


p_sub_sig <- merged_sub %>% subset(SubSignificance == "Yes") %>%
  ggplot(data = .,aes(RangeCels,SumZones))  +
  geom_point(size = 1,alpha = 0.5) +
  stat_summary(fun.data =  mean_cl_normal,geom = "linerange",size = 0.5,color = "Red")+
  stat_summary(fun =  mean,geom = "point",size = 2,color = "Red")+
  scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.1)) +
  theme_ohchibi(font_family = "Helvetica")  +
  theme(
    axis.text.x = element_text(family = "Helvetica",angle =45,vjust = 1,hjust = 1)
  ) +
  xlab(label = "Range total number of cells") + ylab(label = "Sum NoExp + Discr Zones RA")

### Save graphs
oh.save.pdf(p = p_pi,outname = "resub_screeningbarriers_pi_by_tonumcells.pdf",
            outdir = "../figures/",width = 6,height = 6)

### Save graphs
oh.save.pdf(p = p_sub,outname = "resub_screeningbarriers_sub_by_tonumcells.pdf",
            outdir = "../figures/",width = 6,height = 6)


###Check the chosen ones
meta_41 <- read.table(file = "../rawdata/41_isolates.tsv")
st <- meta_41$V2
merged_sub$Chosen <- "No"
merged_sub$Chosen[which(merged_sub$Strains %in% st)] <- "Yes"

p_sub_chosen <- merged_sub %>% subset(Chosen == "Yes") %>%
  ggplot(data = .,aes(Total_cells_Normvalue,SumZones))  +
  geom_point(size = 1,alpha = 0.5) +
  stat_summary(fun.data =  mean_cl_normal,geom = "linerange",size = 0.5,color = "Red")+
  stat_summary(fun =  mean,geom = "point",size = 2,color = "Red")+
  scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.1)) +
  theme_ohchibi(font_family = "Helvetica")  +
  theme(
    axis.text.x = element_text(family = "Helvetica",angle =45,vjust = 1,hjust = 1)
  ) +
  xlab(label = "Range total number of cells") + ylab(label = "Sum NoExp + Discr Zones RA") +
  stat_cor(method = "pearson")


### Do the PI vs new measurements 
Dat_a <- readRDS(file = "../cleandata/dat_root_length_vs_number_of_cells_41.RDS")
Dat_b <- readRDS(file = "../cleandata/dat_lignin_xylem_firstroot_diameter_meristem_41.RDS")

#PReare the structures
df_pi <- Dat_a$df_pi
df_tc <- Dat_a$df_tc
df_rl <- Dat_a$df_rl

Tab <- Dat_b$Dat$Tab %>% t
Mapa <- Dat_b$Dat$Map
x <- Tab[,c("DistanceCSLignin")]
Mapa$DistanceCSLignin <- x
df_cs <- Mapa

#Compute a normalization factor
overall_me <- mean(subset(df_cs,Strain=="NB")$DistanceCSLignin)
df_norm <- NULL
for(batch in unique(df_cs$Date)){
  local_me <- subset(df_cs, Date == batch & Strain == "NB") %$% DistanceCSLignin %>% mean 
  nfactor <- overall_me/local_me
  temp <- subset(df_cs,Date==batch)
  temp$NormDistanceCSLignin <- temp$DistanceCSLignin*nfactor
  df_norm <- rbind(df_norm,temp)
}

#Summarize each data frame
df_pi <- df_pi%>%
  summarySE(data = .,measurevar = "Norm_PI_permeability",groupvars = "Names")
df_tc <- df_tc%>%
  summarySE(data = .,measurevar = "Norm_Total_number_of_cells",groupvars = "Names")
df_rl <- df_rl%>%
  summarySE(data = .,measurevar = "Norm_Root_length",groupvars = "Names")
df_cs <- df_norm%>%
  summarySE(data = .,measurevar = "NormDistanceCSLignin",groupvars = "Strain")
colnames(df_cs)[1] <- "Names"

#Merge structure in one single plot
merged <- merge(df_pi[c("Names","Norm_PI_permeability","ci")],
                df_tc[,c("Names","Norm_Total_number_of_cells","ci")], by = "Names") %>%
  merge(df_rl[,c("Names","Norm_Root_length","ci")], by = "Names") %>%
  merge(df_cs[,c("Names","NormDistanceCSLignin","ci")], by = "Names")

merged <- merged %>% subset(Names != "NB") %>%
  subset(Names != "L61") %>%
  droplevels

#Perform lineal model and test
m1 <- lm(formula = Norm_PI_permeability~Norm_Total_number_of_cells,data = merged)
summary(m1)

colnames(merged)[c(3,5,7,9)] <- c("ci.Norm_PI_permeability","ci.Norm_Total_number_of_cells",
                                  "ci.Norm_Root_length","ci.NormDistanceCSLignin")
#Plot
ggplot(data = merged,aes(Norm_Total_number_of_cells,Norm_PI_permeability)) +
  geom_smooth(method = "lm",se =T,color = "red",fill = "#D9D9D9") + 
  geom_linerange(aes(xmin = Norm_Total_number_of_cells - ci.Norm_Total_number_of_cells,
                     xmax = Norm_Total_number_of_cells + ci.Norm_Total_number_of_cells),
                 color = "#D9D9D9",size = 0.1) +
  geom_linerange(aes(ymin = Norm_PI_permeability - ci.Norm_PI_permeability,
                     ymax = Norm_PI_permeability + ci.Norm_PI_permeability),
                 color = "#D9D9D9",size = 0.1) +
  geom_point(size= 3) +
  theme_ohchibi(font_family = "Helvetica") +
  stat_cor() +
  xlab(label = "Total number of cells") +
  ylab(label = "PI permeability)") +
  scale_y_continuous(limits = c(0,20),oob = rescale_none)


m1 <- lm(formula = Norm_PI_permeability~Norm_Root_length,data = merged)
summary(m1)
ggplot(data = merged,aes(Norm_Root_length,Norm_PI_permeability)) +
  geom_smooth(method = "lm",se =T,color = "red",fill = "#D9D9D9") + 
  geom_linerange(aes(xmin = Norm_Root_length - ci.Norm_Root_length,
                     xmax = Norm_Root_length + ci.Norm_Root_length),
                 color = "#D9D9D9",size = 0.1) +
  geom_linerange(aes(ymin = Norm_PI_permeability - ci.Norm_PI_permeability,
                     ymax = Norm_PI_permeability + ci.Norm_PI_permeability),
                 color = "#D9D9D9",size = 0.1) +
  geom_point(size= 3) +
  theme_ohchibi(font_family = "Helvetica") +
  stat_cor() +
  xlab(label = "Primary root elongation cm") +
  ylab(label = "PI permeability)") +
  scale_y_continuous(limits = c(0,20),oob = rescale_none)



m1 <- lm(formula = Norm_PI_permeability~NormDistanceCSLignin,data = merged)
summary(m1)

ggplot(data = merged,aes(NormDistanceCSLignin,Norm_PI_permeability)) +
  geom_smooth(method = "lm",se =T,color = "red",fill = "#D9D9D9") + 
  geom_linerange(aes(xmin = NormDistanceCSLignin - ci.NormDistanceCSLignin,
                     xmax = NormDistanceCSLignin + ci.NormDistanceCSLignin),
                 color = "#D9D9D9",size = 0.1) +
  geom_linerange(aes(ymin = Norm_PI_permeability - ci.Norm_PI_permeability,
                     ymax = Norm_PI_permeability + ci.Norm_PI_permeability),
                 color = "#D9D9D9",size = 0.1) +
  geom_point(size= 3) +
  theme_ohchibi(font_family = "Helvetica") +
  stat_cor() +
  xlab(label = "Distance CS lignin") +
  ylab(label = "PI permeability)") +
  scale_y_continuous(limits = c(0,20),oob = rescale_none)
