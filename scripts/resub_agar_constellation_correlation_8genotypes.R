library(ohchibi)
library(Rmisc)
library(scales)
library(ggrepel)

set.seed(130816)
setwd('/Users/isaisalasgonzalez/Documents/rootbarriersmicro/scripts/')


### Plotting parameters 
size_legend_text <- 45
size_axis_title <- 35
size_legend_title <- 55
legend_proportion_size <- 4

mix.colors <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", 
                "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", 
                "#6A3D9A", "#FFFF99", "#B15928")

## Select only genotypes that are comparable with the agar dataset
chosen_genotypes_soil <- c("Col_0","myb36_2_sgn3_3","pELTP_CDEFsgn3_3_myb36_1","sgn3_3","esb1_1",
                           "pCASP1_CDEFesb1_1","myb36_2","pCASP1_CDEFwt")

#### Microbiome ########
distfun <- function(x,method) vegan::vegdist(x = x, method = "bray")

####### Agar system ############
Dat_agar <- readRDS(file = "../cleandata/dat_rootbarriers_amplicon_bacteria_syncom_mutants.RDS")

Dat_rar_agar <- Dat_agar$RelativeAbundance

################################
######### Root #################
################################
fraction <- "Root"
Dat_sub_agar <- Dat_rar_agar %>% 
  subset.Dataset(subset = Genotype != "blank",drop = T,clean = T) %>%
  subset.Dataset(subset = Genotype != "noPlant",drop = T,clean = T) %>%
  subset.Dataset(subset = Treatment != "NB",drop = T,clean = T) %>%
  subset.Dataset(subset = Fraction == fraction,drop = T,clean = T) 


mcap <- oh.cap(Tab = Dat_sub_agar$Tab ,Map = Dat_sub_agar$Map ,
               formula = "Genotype + Condition(Replica)",
               distfun = distfun,perms = 9999)

bray_dist <- distfun(x = Dat_sub_agar$Tab %>% t)
mperm <- adonis(formula = bray_dist ~Genotype +Replica,data = Dat_sub_agar$Map,permutations = 9999)

p_all_micro_agar_root <- chibi.cap(list_ohpco = mcap,col_val = "Genotype",
                                   mypch = 21,size = 10,alpha=1,comp_a = "CAP1",comp_b = "CAP2",
                                   size_legend_text = 10,
                                   size_axis_title = size_axis_title,size_axis_line = 2,
                                   size_title_text = 15,
                                   font_family = "Helvetica",legend_proportion_size = 1,lines_zero = T) +
  #theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_color_manual(values = mix.colors) 

### Perform PERMANOVA ####
Tab_bray <- distfun(t(Dat_sub_agar$Tab))
mypermanova <- adonis(Tab_bray ~   Genotype,
                      data = Dat_sub_agar$Map,
                      strata = Dat_sub_agar$Map$Replica,
                      permutations = 9999)
mypermanova




# Summarize
Mapa_cap <- mcap$Map_cap
Mapa_cap$CAP1 <- Mapa_cap$CAP1*-1

df_cap1 <- summarySE(data = Mapa_cap,measurevar = c("CAP1"),groupvars = "Genotype")
df_cap2 <- summarySE(data = Mapa_cap,measurevar = c("CAP2"),groupvars = "Genotype")
df_cap3 <- summarySE(data = Mapa_cap,measurevar = c("CAP3"),groupvars = "Genotype")
df_cap4 <- summarySE(data = Mapa_cap,measurevar = c("CAP4"),groupvars = "Genotype")
df_cap5 <- summarySE(data = Mapa_cap,measurevar = c("CAP5"),groupvars = "Genotype")
df_cap6 <- summarySE(data = Mapa_cap,measurevar = c("CAP6"),groupvars = "Genotype")


merged_cap <- merge(df_cap1,df_cap2, by = "Genotype") %>%
  merge(.,df_cap3, by = "Genotype") %>%
  merge(.,df_cap4, by = "Genotype") %>%
  merge(.,df_cap5, by = "Genotype") %>%
  merge(.,df_cap6, by = "Genotype")
merged_cap <- merged_cap[,c(1,3,6,8,11,13,16,18,21,23,26,28,31)]

c1 <- paste0("CAP1 ",mcap$variance_explained_axis[1],"%")
c2 <- paste0("CAP2 ",mcap$variance_explained_axis[2],"%")


#Append the number information 
## Load the number scheme
df_num <- read.table(file = "../rawdata/df_genotypes_to_classification.csv",header = T,sep = ",")
df_num$Number <- 1:20

merged_cap$Genotype
df_num$Genotype %>% as.character

merged_cap$Genotype <- merged_cap$Genotype %>% 
  gsub(pattern = "^Col-0$",replacement = "Col_0") %>%
  gsub(pattern = "^esb1-1$",replacement = "esb1_1") %>%
  gsub(pattern = "^myb36-2$",replacement = "myb36_2") %>%
  gsub(pattern = "^myb36-2sgn3-3$",replacement = "myb36_2_sgn3_3") %>%
  gsub(pattern = "^pCASP1\\:CDEF\\(esb1-1\\)$",replacement = "pCASP1_CDEFesb1_1") %>%
  gsub(pattern = "^pCASP1\\:CDEF\\(wt\\)$",replacement = "pCASP1_CDEFwt") %>%
  gsub(pattern = "^pELTP\\:CDEF\\(sgn3-3myb36-1\\)$",replacement = "pELTP_CDEFsgn3_3_myb36_1") %>%
  gsub(pattern = "^sgn3-3$",replacement = "sgn3_3")

merged_cap$Number <- match(merged_cap$Genotype,df_num$Genotype) %>%
  df_num$Number[.]

### Test the first two axis against Col-0

Mapa_cap$Genotype <- Mapa_cap$Genotype %>%
  gsub(pattern = "^Col-0$",replacement = "Col_0") %>%
  gsub(pattern = "^esb1-1$",replacement = "esb1_1") %>%
  gsub(pattern = "^myb36-2$",replacement = "myb36_2") %>%
  gsub(pattern = "^myb36-2sgn3-3$",replacement = "myb36_2_sgn3_3") %>%
  gsub(pattern = "^pCASP1\\:CDEF\\(esb1-1\\)$",replacement = "pCASP1_CDEFesb1_1") %>%
  gsub(pattern = "^pCASP1\\:CDEF\\(wt\\)$",replacement = "pCASP1_CDEFwt") %>%
  gsub(pattern = "^pELTP\\:CDEF\\(sgn3-3myb36-1\\)$",replacement = "pELTP_CDEFsgn3_3_myb36_1") %>%
  gsub(pattern = "^sgn3-3$",replacement = "sgn3_3")

#Perform difference ves col 0 in the first two axis
maxis <- c("CAP1","CAP2")
df <- NULL
mgenos <- merged_cap$Genotype %>% grep(pattern = "Col_0",invert = T,value = T)
Mapa_cap$Genotype <- Mapa_cap$Genotype %>% factor %>% relevel(ref = "Col_0")
for(ax in maxis){
  for(geno in mgenos){
    temp <- Mapa_cap %>% subset(Genotype == geno )
    y <- temp[,colnames(temp) %in% ax]
    if(length(y) ==0){
      next
    }
    temp <- Mapa_cap %>% subset(Genotype == "Col_0" )
    x <- temp[,colnames(temp) %in% ax]
    estimate <- mean(y)-mean(x)
    mlevene <- rbind(data.frame(value = y,geno = "other"),data.frame(value =x,geno = "col")) %>%
      car::leveneTest(value ~geno, data = .) 
    if(mlevene$`Pr(>F)`[1] < 0.05){
      pval <- tryCatch(expr =t.test(x,y,var.equal = F) %$% p.value,error = function(e)p.value = NA )
      df <- data.frame(Fraction = fraction,Genotype = geno,Axis = ax,Estimate = estimate,p.value = pval) %>%
        rbind(df,.)
    }else{
      pval <- tryCatch(expr =t.test(x,y,var.equal = T) %$% p.value,error = function(e)p.value = NA )
      df <- data.frame(Fraction = fraction,Genotype = geno,Axis = ax,Estimate = estimate,p.value = pval) %>%
        rbind(df,.)
    }
    
  }
}

df_root <- df
merged_cap_micro_agar_root <- merged_cap

################################
######### Shoot #################
################################
fraction <- "Shoot"
Dat_sub_agar <- Dat_rar_agar %>% 
  subset.Dataset(subset = Genotype != "blank",drop = T,clean = T) %>%
  subset.Dataset(subset = Genotype != "noPlant",drop = T,clean = T) %>%
  subset.Dataset(subset = Treatment != "NB",drop = T,clean = T) %>%
  subset.Dataset(subset = Fraction == fraction,drop = T,clean = T) 


mcap <- oh.cap(Tab = Dat_sub_agar$Tab ,Map = Dat_sub_agar$Map ,
               formula = "Genotype + Condition(Replica)",
               distfun = distfun,perms = 9999)

bray_dist <- distfun(x = Dat_sub_agar$Tab %>% t)
mperm <- adonis(formula = bray_dist ~Genotype +Replica,data = Dat_sub_agar$Map,permutations = 9999)

p_all_micro_agar_shoot <- chibi.cap(list_ohpco = mcap,col_val = "Genotype",
                                   mypch = 21,size = 10,alpha=1,comp_a = "CAP1",comp_b = "CAP2",
                                   size_legend_text = 10,
                                   size_axis_title = size_axis_title,size_axis_line = 2,
                                   size_title_text = 15,
                                   font_family = "Helvetica",legend_proportion_size = 1,lines_zero = T) +
  #theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_color_manual(values = mix.colors) 

### Perform PERMANOVA ####
Tab_bray <- distfun(t(Dat_sub_agar$Tab))
mypermanova <- adonis(Tab_bray ~   Genotype,
                      data = Dat_sub_agar$Map,
                      strata = Dat_sub_agar$Map$Replica,
                      permutations = 9999)
mypermanova



# Summarize
Mapa_cap <- mcap$Map_cap

df_cap1 <- summarySE(data = Mapa_cap,measurevar = c("CAP1"),groupvars = "Genotype")
df_cap2 <- summarySE(data = Mapa_cap,measurevar = c("CAP2"),groupvars = "Genotype")
df_cap3 <- summarySE(data = Mapa_cap,measurevar = c("CAP3"),groupvars = "Genotype")
df_cap4 <- summarySE(data = Mapa_cap,measurevar = c("CAP4"),groupvars = "Genotype")
df_cap5 <- summarySE(data = Mapa_cap,measurevar = c("CAP5"),groupvars = "Genotype")
df_cap6 <- summarySE(data = Mapa_cap,measurevar = c("CAP6"),groupvars = "Genotype")


merged_cap <- merge(df_cap1,df_cap2, by = "Genotype") %>%
  merge(.,df_cap3, by = "Genotype") %>%
  merge(.,df_cap4, by = "Genotype") %>%
  merge(.,df_cap5, by = "Genotype") %>%
  merge(.,df_cap6, by = "Genotype")
merged_cap <- merged_cap[,c(1,3,6,8,11,13,16,18,21,23,26,28,31)]

c1 <- paste0("CAP1 ",mcap$variance_explained_axis[1],"%")
c2 <- paste0("CAP2 ",mcap$variance_explained_axis[2],"%")


#Append the number information 
## Load the number scheme
df_num <- read.table(file = "../rawdata/df_genotypes_to_classification.csv",header = T,sep = ",")
df_num$Number <- 1:20

merged_cap$Genotype
df_num$Genotype %>% as.character

merged_cap$Genotype <- merged_cap$Genotype %>% 
  gsub(pattern = "^Col-0$",replacement = "Col_0") %>%
  gsub(pattern = "^esb1-1$",replacement = "esb1_1") %>%
  gsub(pattern = "^myb36-2$",replacement = "myb36_2") %>%
  gsub(pattern = "^myb36-2sgn3-3$",replacement = "myb36_2_sgn3_3") %>%
  gsub(pattern = "^pCASP1\\:CDEF\\(esb1-1\\)$",replacement = "pCASP1_CDEFesb1_1") %>%
  gsub(pattern = "^pCASP1\\:CDEF\\(wt\\)$",replacement = "pCASP1_CDEFwt") %>%
  gsub(pattern = "^pELTP\\:CDEF\\(sgn3-3myb36-1\\)$",replacement = "pELTP_CDEFsgn3_3_myb36_1") %>%
  gsub(pattern = "^sgn3-3$",replacement = "sgn3_3")

merged_cap$Number <- match(merged_cap$Genotype,df_num$Genotype) %>%
  df_num$Number[.]

### Test the first two axis against Col-0

Mapa_cap$Genotype <- Mapa_cap$Genotype %>%
  gsub(pattern = "^Col-0$",replacement = "Col_0") %>%
  gsub(pattern = "^esb1-1$",replacement = "esb1_1") %>%
  gsub(pattern = "^myb36-2$",replacement = "myb36_2") %>%
  gsub(pattern = "^myb36-2sgn3-3$",replacement = "myb36_2_sgn3_3") %>%
  gsub(pattern = "^pCASP1\\:CDEF\\(esb1-1\\)$",replacement = "pCASP1_CDEFesb1_1") %>%
  gsub(pattern = "^pCASP1\\:CDEF\\(wt\\)$",replacement = "pCASP1_CDEFwt") %>%
  gsub(pattern = "^pELTP\\:CDEF\\(sgn3-3myb36-1\\)$",replacement = "pELTP_CDEFsgn3_3_myb36_1") %>%
  gsub(pattern = "^sgn3-3$",replacement = "sgn3_3")

#Perform difference ves col 0 in the first two axis
maxis <- c("CAP1","CAP2")
df <- NULL
mgenos <- merged_cap$Genotype %>% grep(pattern = "Col_0",invert = T,value = T)
Mapa_cap$Genotype <- Mapa_cap$Genotype %>% factor %>% relevel(ref = "Col_0")
for(ax in maxis){
  for(geno in mgenos){
    temp <- Mapa_cap %>% subset(Genotype == geno )
    y <- temp[,colnames(temp) %in% ax]
    if(length(y) ==0){
      next
    }
    temp <- Mapa_cap %>% subset(Genotype == "Col_0" )
    x <- temp[,colnames(temp) %in% ax]
    estimate <- mean(y)-mean(x)
    mlevene <- rbind(data.frame(value = y,geno = "other"),data.frame(value =x,geno = "col")) %>%
      car::leveneTest(value ~geno, data = .) 
    if(mlevene$`Pr(>F)`[1] < 0.05){
      pval <- tryCatch(expr =t.test(x,y,var.equal = F) %$% p.value,error = function(e)p.value = NA )
      df <- data.frame(Fraction = fraction,Genotype = geno,Axis = ax,Estimate = estimate,p.value = pval) %>%
        rbind(df,.)
    }else{
      pval <- tryCatch(expr =t.test(x,y,var.equal = T) %$% p.value,error = function(e)p.value = NA )
      df <- data.frame(Fraction = fraction,Genotype = geno,Axis = ax,Estimate = estimate,p.value = pval) %>%
        rbind(df,.)
    }
    
  }
}

df_shoot <- df
merged_cap_micro_agar_shoot <- merged_cap


################################
######### Agar #################
################################
fraction <- "Agar"
Dat_sub_agar <- Dat_rar_agar %>% 
  subset.Dataset(subset = Genotype != "blank",drop = T,clean = T) %>%
  subset.Dataset(subset = Genotype != "noPlant",drop = T,clean = T) %>%
  subset.Dataset(subset = Treatment != "NB",drop = T,clean = T) %>%
  subset.Dataset(subset = Fraction == fraction,drop = T,clean = T) 


mcap <- oh.cap(Tab = Dat_sub_agar$Tab ,Map = Dat_sub_agar$Map ,
               formula = "Genotype + Condition(Replica)",
               distfun = distfun,perms = 9999)

bray_dist <- distfun(x = Dat_sub_agar$Tab %>% t)
mperm <- adonis(formula = bray_dist ~Genotype +Replica,data = Dat_sub_agar$Map,permutations = 9999)

p_all_micro_agar_agar <- chibi.cap(list_ohpco = mcap,col_val = "Genotype",
                                    mypch = 21,size = 10,alpha=1,comp_a = "CAP1",comp_b = "CAP2",
                                    size_legend_text = 10,
                                    size_axis_title = size_axis_title,size_axis_line = 2,
                                    size_title_text = 15,
                                    font_family = "Helvetica",legend_proportion_size = 1,lines_zero = T) +
  #theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_color_manual(values = mix.colors) 

### Perform PERMANOVA ####
Tab_bray <- distfun(t(Dat_sub_agar$Tab))
mypermanova <- adonis(Tab_bray ~   Genotype,
                      data = Dat_sub_agar$Map,
                      strata = Dat_sub_agar$Map$Replica,
                      permutations = 9999)
mypermanova


# Summarize
Mapa_cap <- mcap$Map_cap

df_cap1 <- summarySE(data = Mapa_cap,measurevar = c("CAP1"),groupvars = "Genotype")
df_cap2 <- summarySE(data = Mapa_cap,measurevar = c("CAP2"),groupvars = "Genotype")
df_cap3 <- summarySE(data = Mapa_cap,measurevar = c("CAP3"),groupvars = "Genotype")
df_cap4 <- summarySE(data = Mapa_cap,measurevar = c("CAP4"),groupvars = "Genotype")
df_cap5 <- summarySE(data = Mapa_cap,measurevar = c("CAP5"),groupvars = "Genotype")
df_cap6 <- summarySE(data = Mapa_cap,measurevar = c("CAP6"),groupvars = "Genotype")


merged_cap <- merge(df_cap1,df_cap2, by = "Genotype") %>%
  merge(.,df_cap3, by = "Genotype") %>%
  merge(.,df_cap4, by = "Genotype") %>%
  merge(.,df_cap5, by = "Genotype") %>%
  merge(.,df_cap6, by = "Genotype")
merged_cap <- merged_cap[,c(1,3,6,8,11,13,16,18,21,23,26,28,31)]

c1 <- paste0("CAP1 ",mcap$variance_explained_axis[1],"%")
c2 <- paste0("CAP2 ",mcap$variance_explained_axis[2],"%")


#Append the number information 
## Load the number scheme
df_num <- read.table(file = "../rawdata/df_genotypes_to_classification.csv",header = T,sep = ",")
df_num$Number <- 1:20

merged_cap$Genotype
df_num$Genotype %>% as.character

merged_cap$Genotype <- merged_cap$Genotype %>% 
  gsub(pattern = "^Col-0$",replacement = "Col_0") %>%
  gsub(pattern = "^esb1-1$",replacement = "esb1_1") %>%
  gsub(pattern = "^myb36-2$",replacement = "myb36_2") %>%
  gsub(pattern = "^myb36-2sgn3-3$",replacement = "myb36_2_sgn3_3") %>%
  gsub(pattern = "^pCASP1\\:CDEF\\(esb1-1\\)$",replacement = "pCASP1_CDEFesb1_1") %>%
  gsub(pattern = "^pCASP1\\:CDEF\\(wt\\)$",replacement = "pCASP1_CDEFwt") %>%
  gsub(pattern = "^pELTP\\:CDEF\\(sgn3-3myb36-1\\)$",replacement = "pELTP_CDEFsgn3_3_myb36_1") %>%
  gsub(pattern = "^sgn3-3$",replacement = "sgn3_3")

merged_cap$Number <- match(merged_cap$Genotype,df_num$Genotype) %>%
  df_num$Number[.]

### Test the first two axis against Col-0

Mapa_cap$Genotype <- Mapa_cap$Genotype %>%
  gsub(pattern = "^Col-0$",replacement = "Col_0") %>%
  gsub(pattern = "^esb1-1$",replacement = "esb1_1") %>%
  gsub(pattern = "^myb36-2$",replacement = "myb36_2") %>%
  gsub(pattern = "^myb36-2sgn3-3$",replacement = "myb36_2_sgn3_3") %>%
  gsub(pattern = "^pCASP1\\:CDEF\\(esb1-1\\)$",replacement = "pCASP1_CDEFesb1_1") %>%
  gsub(pattern = "^pCASP1\\:CDEF\\(wt\\)$",replacement = "pCASP1_CDEFwt") %>%
  gsub(pattern = "^pELTP\\:CDEF\\(sgn3-3myb36-1\\)$",replacement = "pELTP_CDEFsgn3_3_myb36_1") %>%
  gsub(pattern = "^sgn3-3$",replacement = "sgn3_3")

#Perform difference ves col 0 in the first two axis
maxis <- c("CAP1","CAP2")
df <- NULL
mgenos <- merged_cap$Genotype %>% grep(pattern = "Col_0",invert = T,value = T)
Mapa_cap$Genotype <- Mapa_cap$Genotype %>% factor %>% relevel(ref = "Col_0")
for(ax in maxis){
  for(geno in mgenos){
    temp <- Mapa_cap %>% subset(Genotype == geno )
    y <- temp[,colnames(temp) %in% ax]
    if(length(y) ==0){
      next
    }
    temp <- Mapa_cap %>% subset(Genotype == "Col_0" )
    x <- temp[,colnames(temp) %in% ax]
    estimate <- mean(y)-mean(x)
    mlevene <- rbind(data.frame(value = y,geno = "other"),data.frame(value =x,geno = "col")) %>%
      car::leveneTest(value ~geno, data = .) 
    if(mlevene$`Pr(>F)`[1] < 0.05){
      pval <- tryCatch(expr =t.test(x,y,var.equal = F) %$% p.value,error = function(e)p.value = NA )
      df <- data.frame(Fraction = fraction,Genotype = geno,Axis = ax,Estimate = estimate,p.value = pval) %>%
        rbind(df,.)
    }else{
      pval <- tryCatch(expr =t.test(x,y,var.equal = T) %$% p.value,error = function(e)p.value = NA )
      df <- data.frame(Fraction = fraction,Genotype = geno,Axis = ax,Estimate = estimate,p.value = pval) %>%
        rbind(df,.)
    }
    
  }
}

df_agar <- df
merged_cap_micro_agar_agar <- merged_cap


df <- rbind(df_root,df_shoot,df_agar)
df$p.adj <- df$p.value %>% p.adjust(method = "bonferroni")

res_lm <- df
res_lm$Significance <- "NS"
res_lm$Significance[which(res_lm$p.adj < 0.1)] <- "Significant"
res_lm$Count <- 1

df_resum <- dcast(data = res_lm %>% subset(Significance == "Significant"),
      formula = Genotype~Fraction,value.var = "Count",fun.aggregate = sum,fill = 0) %>%
  melt
df_resum$value[which(df_resum$value>1)]  <- 1

df_resum$Sig <- df_resum$value %>% gsub(pattern = "0",replacement = "NS") %>%
  gsub(pattern = "1",replacement = "Significant")

merged_cap_micro_agar_root$Fraction <- "Root"
merged_cap_micro_agar_shoot$Fraction <- "Shoot"
merged_cap_micro_agar_agar$Fraction <- "Agar"

merged_cap <- rbind(merged_cap_micro_agar_root,merged_cap_micro_agar_shoot,merged_cap_micro_agar_agar)


## Cosntellation plots
colnames(df_resum)[2] <- "Fraction"
merged_cap <- merge(merged_cap,df_resum[,c("Genotype","Fraction","Sig")],by = c("Genotype","Fraction"),all = T)
merged_cap$Sig <- merged_cap$Sig %>% factor

paleta <- c("black","red")
names(paleta) <- c("NS","Significant")

p_microbiome_root <- ggplot(data = merged_cap %>% subset(Fraction == "Root"),aes(CAP1,CAP2)) + 
  geom_vline(xintercept = 0,linetype = "longdash",size =2,color = "#D9D9D9") +
  geom_hline(yintercept = 0,linetype = "longdash",size =2,color = "#D9D9D9") +
  geom_linerange(aes(ymin = CAP2 - ci.y,ymax = CAP2 + ci.y),size = 0.1,alpha = 0.2) +
  geom_linerange(aes(xmin = CAP1 - ci.x,xmax = CAP1 + ci.x),size = 0.1,alpha = 0.2) +
  geom_text_repel(aes(label = Number,color = Sig),size = 8) +
  theme_ohchibi(size_panel_border = 2) + 
  theme(
    legend.position = "none",
    strip.background.x = element_blank(),
    strip.text.x = element_text(family = "Arial",face = "bold",size = 25) 
  )   +
  scale_color_manual(values = paleta,na.value = "black") +
  scale_y_continuous(limits = c(-2.3,1)) +
  scale_x_continuous(limits = c(-1,1.5),oob = rescale_none)

p_microbiome_shoot <- ggplot(data = merged_cap %>% subset(Fraction == "Shoot"),aes(CAP1,CAP2)) + 
  geom_vline(xintercept = 0,linetype = "longdash",size =2,color = "#D9D9D9") +
  geom_hline(yintercept = 0,linetype = "longdash",size =2,color = "#D9D9D9") +
  geom_linerange(aes(ymin = CAP2 - ci.y,ymax = CAP2 + ci.y),size = 0.1,alpha = 0.2) +
  geom_linerange(aes(xmin = CAP1 - ci.x,xmax = CAP1 + ci.x),size = 0.1,alpha = 0.2) +
  geom_text_repel(aes(label = Number,color = Sig),size = 8) +
  theme_ohchibi(size_panel_border = 2) + 
  theme(
    legend.position = "none",
    strip.background.x = element_blank(),
    strip.text.x = element_text(family = "Arial",face = "bold",size = 25) 
  )   +
  scale_color_manual(values = paleta,na.value = "black")


p_microbiome_agar <- ggplot(data = merged_cap %>% subset(Fraction == "Agar"),aes(CAP1,CAP2)) + 
  geom_vline(xintercept = 0,linetype = "longdash",size =2,color = "#D9D9D9") +
  geom_hline(yintercept = 0,linetype = "longdash",size =2,color = "#D9D9D9") +
  geom_linerange(aes(ymin = CAP2 - ci.y,ymax = CAP2 + ci.y),size = 0.1,alpha = 0.2) +
  geom_linerange(aes(xmin = CAP1 - ci.x,xmax = CAP1 + ci.x),size = 0.1,alpha = 0.2) +
  geom_text_repel(aes(label = Number,color = Sig),size = 8) +
  theme_ohchibi(size_panel_border = 2) + 
  theme(
    legend.position = "none",
    strip.background.x = element_blank(),
    strip.text.x = element_text(family = "Arial",face = "bold",size = 25) 
  )   +
  scale_color_manual(values = paleta,na.value = "black")

merged_cap_microbiome <- merged_cap

##### Ionome #############
#### Agar  system ######
Dat_ion_agar <- readRDS(file = "../cleandata/dat_syncom_genotypes.RDS") %$%
  Ionome
Tab_z <- Dat_ion_agar$Tab %>% t %>% scale %>% t 

Dat_ion_agar <- create_dataset(Tab = Tab_z,Map = Dat_ion_agar$Map)


Dat_ion_agar <- Dat_ion_agar %>% subset.Dataset(Type == "SynCom",drop = T,clean = T)

#Change the name of the genotypes to match previous consensus
Dat_ion_agar$Map$Genotype <- Dat_ion_agar$Map$Genotype %>% 
  gsub(pattern = "Col-0",replacement = "Col_0") %>%
  gsub(pattern = "^sgn3myb36$",replacement = "myb36_2_sgn3_3") %>%
  gsub(pattern = "sgn3.3myb36pELTP::CEDF",replacement = "pELTP_CDEFsgn3_3_myb36_1") %>%  
  gsub(pattern = "^sgn3$",replacement = "sgn3_3") %>%
  gsub(pattern = "^esb1$",replacement = "esb1_1") %>%
  gsub(pattern = "esb1.1pCASP1::CDEF",replacement = "pCASP1_CDEFesb1_1") %>%
  gsub(pattern = "^myb36$" ,replacement = "myb36_2") %>%
  gsub(pattern = "^pCASP1::CDEF$",replacement = "pCASP1_CDEFwt")

samples_to_remove <- which(!(Dat_ion_agar$Map$Genotype %in% chosen_genotypes_soil)) %>%
  Dat_ion_agar$Map$Index[.]

Dat_ion_agar <- remove_samples(Dat = Dat_ion_agar,samples = samples_to_remove,droplevels = T)

#Define Bray Curtis as the dissimilarity
distfun <- function(x,method) vegan::vegdist(x = x, method = "euclidean")


#Control for the plot effect observed
mcap <-oh.cap(Tab = Dat_ion_agar$Tab,Map = Dat_ion_agar$Map,
              formula = "Genotype + Condition(Replicates)",
              distfun = distfun,perms = 9999,sqrt = F)


p_all_ion_agar_shoot <- chibi.cap(list_ohpco = mcap,col_val = "Genotype",
                                  mypch = 21,size = 10,alpha=1,comp_a = "CAP1",comp_b = "CAP2",
                                  size_legend_text = 10,
                                  size_axis_title = size_axis_title,size_axis_line = 2,
                                  size_title_text = 15,
                                  font_family = "Helvetica",legend_proportion_size = 1,lines_zero = T) +
  #theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_color_manual(values = mix.colors) 

mcap_ion_agar_shoot <- mcap$Map_cap

# Summarize
Mapa_cap <- mcap$Map_cap

df_cap1 <- summarySE(data = Mapa_cap,measurevar = c("CAP1"),groupvars = "Genotype")
df_cap2 <- summarySE(data = Mapa_cap,measurevar = c("CAP2"),groupvars = "Genotype")
df_cap3 <- summarySE(data = Mapa_cap,measurevar = c("CAP3"),groupvars = "Genotype")
df_cap4 <- summarySE(data = Mapa_cap,measurevar = c("CAP4"),groupvars = "Genotype")
df_cap5 <- summarySE(data = Mapa_cap,measurevar = c("CAP5"),groupvars = "Genotype")
df_cap6 <- summarySE(data = Mapa_cap,measurevar = c("CAP6"),groupvars = "Genotype")


merged_cap <- merge(df_cap1,df_cap2, by = "Genotype") %>%
  merge(.,df_cap3, by = "Genotype") %>%
  merge(.,df_cap4, by = "Genotype") %>%
  merge(.,df_cap5, by = "Genotype") %>%
  merge(.,df_cap6, by = "Genotype")
merged_cap <- merged_cap[,c(1,3,6,8,11,13,16,18,21,23,26,28,31)]

c1 <- paste0("CAP1 ",mcap$variance_explained_axis[1],"%")
c2 <- paste0("CAP2 ",mcap$variance_explained_axis[2],"%")

merged_cap_ion_agar_shoot <- merged_cap


#Append the number information to the genotype
merged_cap$Number <- match(merged_cap$Genotype,df_num$Genotype) %>%
  df_num$Number[.]

#Perform difference ves col 0 in the first two axis
maxis <- c("CAP1","CAP2")
df <- NULL
mgenos <- merged_cap$Genotype %>% grep(pattern = "Col_0",invert = T,value = T)
Mapa_cap$Genotype <- Mapa_cap$Genotype %>% factor %>% relevel(ref = "Col_0")
for(ax in maxis){
  for(geno in mgenos){
    temp <- Mapa_cap %>% subset(Genotype == geno )
    y <- temp[,colnames(temp) %in% ax]
    if(length(y) ==0){
      next
    }
    temp <- Mapa_cap %>% subset(Genotype == "Col_0" )
    x <- temp[,colnames(temp) %in% ax]
    estimate <- mean(y)-mean(x)
    mlevene <- rbind(data.frame(value = y,geno = "other"),data.frame(value =x,geno = "col")) %>%
      car::leveneTest(value ~geno, data = .) 
    if(mlevene$`Pr(>F)`[1] < 0.05){
      pval <- tryCatch(expr =t.test(x,y,var.equal = F) %$% p.value,error = function(e)p.value = NA )
      df <- data.frame(Fraction = fraction,Genotype = geno,Axis = ax,Estimate = estimate,p.value = pval) %>%
        rbind(df,.)
    }else{
      pval <- tryCatch(expr =t.test(x,y,var.equal = T) %$% p.value,error = function(e)p.value = NA )
      df <- data.frame(Fraction = fraction,Genotype = geno,Axis = ax,Estimate = estimate,p.value = pval) %>%
        rbind(df,.)
    }
    
  }
}

df$p.adj <- df$p.value %>% p.adjust(method = "bonferroni")
df$Significance <- "NS"
df$Significance[which(df$p.adj < 0.1)] <- "Significant"

df %>% subset(Significance == "Significant") %$% Genotype %>% table %>% length

p_ionome_agar <- ggplot(data = merged_cap,aes(CAP1,CAP2)) + 
  geom_vline(xintercept = 0,linetype = "longdash",size =2,color = "#D9D9D9") +
  geom_hline(yintercept = 0,linetype = "longdash",size =2,color = "#D9D9D9") +
  geom_linerange(aes(ymin = CAP2 - ci.y,ymax = CAP2 + ci.y),size = 0.1,alpha = 0.2) +
  geom_linerange(aes(xmin = CAP1 - ci.x,xmax = CAP1 + ci.x),size = 0.1,alpha = 0.2) +
  geom_text_repel(aes(label = Number),size = 8,color = "red") +
  theme_ohchibi(size_panel_border = 2) + 
  theme(
    legend.position = "none",
    strip.background.x = element_blank(),
    strip.text.x = element_text(family = "Arial",face = "bold",size = 25) 
  )   +
  scale_y_continuous(limits = c(-2.3,1)) +
  scale_x_continuous(limits = c(-1,1.5),oob = rescale_none)


df_resum <- merged_cap_microbiome[,c("Genotype","CAP1","CAP2","Fraction")]

df_sub <- df_resum %>% subset(Fraction == "Root")
Tab_sub <- df_sub[,c("CAP1","CAP2")]
rownames(Tab_sub) <- df_sub$Genotype
mdist_root <- dist(Tab_sub)

df_sub <- df_resum %>% subset(Fraction == "Shoot")
Tab_sub <- df_sub[,c("CAP1","CAP2")]
rownames(Tab_sub) <- df_sub$Genotype
mdist_shoot <- dist(Tab_sub)

df_sub <- df_resum %>% subset(Fraction == "Agar")
Tab_sub <- df_sub[,c("CAP1","CAP2")]
rownames(Tab_sub) <- df_sub$Genotype
mdist_agar <- dist(Tab_sub)

df_ip <- merged_cap[,c("Genotype","CAP1","CAP2")]
Tab_ip <- df_ip[,c("CAP1","CAP2")] %>% as.matrix
rownames(Tab_ip) <- df_ip$Genotype 

mdist_ip <- dist(Tab_ip)


m_root <- mantel(mdist_ip,mdist_root,permutations = 9999)
m_shoot <- mantel(mdist_ip,mdist_shoot,permutations = 9999)
m_soil <- mantel(mdist_ip,mdist_agar,permutations = 9999)

mtext_root <- paste0("r = ",round(m_root$statistic,3),"\np = ",format.pval(m_root$signif))
mtext_shoot <- paste0("r = ",round(m_shoot$statistic,3),"\np = ",format.pval(m_shoot$signif))
mtext_soil <- paste0("r = ",round(m_soil$statistic,3),"\np = ",format.pval(m_soil$signif))


### Plot correlations ###
merged_mantel <- mdist_root %>% as.matrix %>% melt_dist() %>%
  merge(
    mdist_shoot %>% as.matrix %>% melt_dist(),
    by = c("iso1","iso2")
    
  )%>%
  merge(
    mdist_agar %>% as.matrix %>% melt_dist(),
    by = c("iso1","iso2")
    
  )
colnames(merged_mantel)[3:5] <- c("Root","Shoot","Agar")

melted <- merged_mantel %>% melt

melted <- mdist_ip %>% as.matrix %>% melt_dist() %>%
  merge(melted,., by = c("iso1","iso2"))


p1 <- ggplot(melted %>% subset(variable == "Agar"),aes(dist,value)) +
  geom_smooth(method = "lm",color = "red",size = 2,fill = "#D9D9D9") +
  geom_point(size = 3) +
  theme_ohchibi(size_panel_border = 2) + 
  annotate(geom = "text",x = 2,y = 0.6,label = mtext_soil) +
  #scale_y_continuous(limits = c(0,3.5)) +
  #scale_x_continuous(limits = c(0,2.5),oob = rescale_none) +
  xlab(label = "ShootIonome") + ylab(label = "Agar")

p2 <- ggplot(melted %>% subset(variable == "Root"),aes(dist,value)) +
  geom_smooth(method = "lm",color = "red",size = 2,fill = "#D9D9D9") +
  geom_point(size = 3) +
  theme_ohchibi(size_panel_border = 2) + 
  annotate(geom = "text",x = 2,y = 0.6,label = mtext_root) +
  #scale_y_continuous(limits = c(0,3.5)) +
  scale_x_continuous(limits = c(0,2.5),oob = rescale_none) +
  xlab(label = "ShootIonome") + ylab(label = "Root")


p3 <- ggplot(melted %>% subset(variable == "Shoot"),aes(dist,value)) +
  geom_smooth(method = "lm",color = "red",size = 2,fill = "#D9D9D9") +
  geom_point(size = 3) +
  theme_ohchibi(size_panel_border = 2) + 
  annotate(geom = "text",x = 2,y = 0.6,label = mtext_shoot) +
  #scale_y_continuous(limits = c(0,3.5)) +
  #scale_x_continuous(limits = c(0,2.5),oob = rescale_none) +
  xlab(label = "ShootIonome") + ylab(label = "Shoot")



## Save the plots in a structue to use later 
###Save plots ####
mplots <- list(
  agar_correlation_ionome_shoot_amplicon_agar = p1,
  agar_correlation_ionome_shoot_amplicon_root = p2,
  agar_correlation_ionome_shoot_amplicon_shoot = p3,
  agar_constellation_root_amplicon = p_microbiome_root,
  agar_constellation_shoot_amplicon = p_microbiome_shoot,
  agar_constellation_agar_amplicon = p_microbiome_agar,
  agar_constellation_shoot_ionome = p_ionome_agar
  
)

saveRDS(object = mplots,
        file = "../cleandata/plots_correlations_ionome_vs_amplicon_cap_agarsystem.RDS")
rm(list=ls())
gc()
