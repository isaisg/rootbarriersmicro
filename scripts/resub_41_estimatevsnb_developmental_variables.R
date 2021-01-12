library(ohchibi)
library(ggpubr)
library(Rmisc)
library(ggrepel)
library(scales)
library(multcomp)

#
set.seed(130816)
setwd('/Users/isaisalasgonzalez/Documents/rootbarriersmicro/scripts/')

Dat_a <- readRDS(file = "../cleandata/dat_root_length_vs_number_of_cells_41.RDS")
Dat_b <- readRDS(file = "../cleandata/dat_lignin_xylem_firstroot_diameter_meristem_41.RDS")


### Correlation 
merged <- Dat_a$df_sum

p1 <- ggplot(data = merged,aes(Norm_Total_number_of_cells,Norm_Root_length)) +
  geom_smooth(method = "lm",se =T,color = "red",fill = "#D9D9D9") + 
  geom_linerange(aes(xmin = Norm_Total_number_of_cells - ci.Norm_Total_number_of_cells,
                     xmax = Norm_Total_number_of_cells + ci.Norm_Total_number_of_cells),
                 color = "#D9D9D9",size = 0.1) +
  geom_linerange(aes(ymin = Norm_Root_length - ci.Norm_Root_length,
                     ymax = Norm_Root_length + ci.Norm_Root_length),
                 color = "#D9D9D9",size = 0.1) +
  geom_point(size= 3) +
  theme_ohchibi(font_family = "Helvetica") +
  stat_cor() +
  xlab(label = "Total number of cells") +
  ylab(label = "Primary root elongation (cm)")


#### Violin plots ####

#Prepare the primary root information
df_temp <- Dat_a$df_rl
df_temp$Norm_Root_length <- df_temp$Norm_Root_length %>% scale

df_temp$Names <- df_temp$Names %>% factor %>%
  relevel(ref = "NB")

df_root <- aov(formula = Norm_Root_length~Names,data = df_temp) %>%
  glht(linfct = mcp(Names = "Dunnett")) %>%
  summary

df_res <- data.frame(Estimate = df_root$test$coefficients,
                     padj =df_root$test$pvalues)
df_res$Variable <- "Primaryroot"  
df_res$Strain <- rownames(df_res) %>%gsub(pattern = " .*",replacement = "")

df_res_root <- df_res

#Prepare the total cells
df_temp <- Dat_a$df_tc
df_temp$Norm_Total_number_of_cells <- df_temp$Norm_Total_number_of_cells %>% scale

df_temp$Names <- df_temp$Names %>% factor %>%
  relevel(ref = "NB")

df_root <- aov(formula = Norm_Total_number_of_cells~Names,data = df_temp) %>%
  glht(linfct = mcp(Names = "Dunnett")) %>%
  summary

df_res <- data.frame(Estimate = df_root$test$coefficients,
                     padj =df_root$test$pvalues)
df_res$Variable <- "TotalCells"  
df_res$Strain <- rownames(df_res) %>%gsub(pattern = " .*",replacement = "")

df_res_cells <- df_res


## Perform testing of the measured variables
Tab <- Dat_b$Dat_z$Tab %>% t
df <- Dat_b$Dat_z$Map

### Perform testing
df_z <- Tab %>% melt
tobind <- match(df_z$Var1,df$Id) %>%
  df[.,c("Date","Strain")]
df_z <- cbind(df_z,tobind)
mvars <- df_z$Var2 %>% as.character %>% unique

mvars <- c("DistanceCSLignin","MeristemSize","DistanceFirstRootHair","RootDiameter",
           "DistanceXylem",
           "NumCellsCSLigning","NumCellsXylem","NumCellsFirstRootHair")

Res <- NULL
for(v in mvars){
  df_sub <- df_z %>% subset(Var2 == v) %>% droplevels
  m1 <- df_sub %>%
    aov(formula = value~Strain + Date,data = .)
  m1_res <- multcomp::glht(m1,linfct=mcp(Strain="Dunnett")) %>% summary
  df_res <- data.frame(Estimate = m1_res$test$coefficients,padj =m1_res$test$pvalues)
  df_res$Variable <- v  
  Res <- rbind(Res,df_res)
}

Res$Strain <- rownames(Res) %>% gsub(pattern = " .*",replacement = "")


Res <- rbind(df_res_root,df_res_cells,Res)

Res_sub <- Res %>%
  subset (Variable != "DistanceCSLignin") %>%
  subset(Variable != "NumCellsCSLigning") %>% droplevels

Res_lignin <- Res %>%
  subset (Variable == "DistanceCSLignin" | Variable == "NumCellsCSLigning") %>%
  droplevels

Res_sub$padj <- Res_sub$padj %>% p.adjust(method = "fdr")

Res <- rbind(Res_lignin,Res_sub)
saveRDS(object = Res,file = "../cleandata/res_41_developmentvariables.RDS")

### For the plotting only use 
Res <- Res %>%
  subset (Variable != "DistanceCSLignin") %>%
  subset(Variable != "NumCellsCSLigning") %>% droplevels


### Adjust levels for plotting
Res$Type <- "Distance"
Res$Type[Res$Variable %>% grep(pattern = "Cell")] <- "Cells"
Res$Type <- Res$Type %>% factor(levels = c("Distance","Cells"))
Res$Variable <- Res$Variable %>% 
  factor(levels = c("DistanceFirstRootHair","DistanceXylem","Primaryroot","MeristemSize","RootDiameter",
                    "NumCellsFirstRootHair","NumCellsXylem","TotalCells"))


p2 <- Res %>% 
  ggplot(data = .,aes(Variable,Estimate)) +
  geom_hline(yintercept = 0,color = "#D9D9D9",size = 0.5) +
  geom_boxplot(fill = NA,aes(color = Variable),outlier.size = 0,outlier.stroke = 0) +
  facet_grid(.~Type,scales = "free",space = "free") +
  geom_sina(aes(color = Variable)) +
  theme_ohchibi(font_family = "Helvetica") +
  ylab(label = "Estimate vs NB") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(family = "Helvetica",angle = 45,hjust = 1,vjust = 1,size = 15),
    strip.text.x = element_text(family = "Helvetica",size = 15),
    strip.background.x = element_blank()
    
  ) +
  scale_color_paletteer_d("ggthemes::colorblind")

#Visualize coloring based on significance
Res$Significance <- "NS"
Res$Significance[which(Res$padj < 0.1)] <- "q < 0.1"

p3 <- Res %>% 
  ggplot(data = .,aes(Variable,Estimate)) +
  geom_hline(yintercept = 0,color = "#D9D9D9",size = 0.5) +
  geom_boxplot(fill = NA,fill = NA,outlier.size = 0,outlier.stroke = 0) +
  facet_grid(.~Type,scales = "free",space = "free") +
  geom_jitter(aes(color = Significance),width = c(0.2)) +
  theme_ohchibi(font_family = "Helvetica") +
  ylab(label = "Estimate vs NB") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(family = "Helvetica",angle = 45,hjust = 1,vjust = 1,size = 15),
    strip.text.x = element_text(family = "Helvetica",size = 15),
    strip.background.x = element_blank()
    
  ) +
  scale_color_manual(values = c("black","red")) +
  scale_y_continuous(limits = c(-4.6,1.5),breaks = c(-4,-3,-2,-1,0,1))

#Bars
df_count <- Res %>% subset(Significance != "NS") %$% Variable %>%
  table %>% data.frame()
colnames(df_count)[1] <- "Variable"
df_count$Perc <- (df_count$Freq/41)*100



#Save graphs ###
oh.save.pdf(p = p1,outname = "resub_cor_totcells_primaryroot.pdf",
            outdir = "../figures/",width =4,height = 4)


#Save graphs ##
oh.save.pdf(p = p3,outname = "resub_estimatevsnb_development_boxplots.pdf",
            outdir = "../figures/",width =12,height = 6)


rm(list=ls())
