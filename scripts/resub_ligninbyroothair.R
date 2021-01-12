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

df <- Dat_b$Dat_z$Tab %>% t
Map <- match(rownames(df),Dat_b$Dat_z$Map$Id) %>%
  Dat_b$Dat_z$Map[.,]
df <- cbind(df,Map[,-4])

#Test distance first root hair vs lignin
mids <- df$Strain %>% levels
mids <- mids %>% grep(pattern = "NB",invert = T,value =  T)
Res <- NULL
for(id in mids){
  df_sub <- df %>% subset(Strain == id) %>% droplevels
  m1 <- lm(formula = DistanceCSLignin~DistanceFirstRootHair,data = df_sub)
  m1_res <- summary(m1)
  df_coef <- m1_res$coefficients %>% as.data.frame
  mest <- df_coef[2,1]
  mpval <- df_coef[2,4]
  mar <- m1_res$adj.r.squared
  Res <- data.frame(Strain = id,Estimate = mest, p.value = mpval,R2 = mar) %>%
    rbind(Res,.)
}

#Test for NB
df_sub <- df %>% subset(Strain == "NB") %>% droplevels
m1 <- lm(formula = DistanceCSLignin~DistanceFirstRootHair +Date,data = df_sub)
m1_res <- summary(m1)
df_coef <- m1_res$coefficients %>% as.data.frame
mest <- df_coef[2,1]
mpval <- df_coef[2,4]
mar <- m1_res$adj.r.squared
Res <- data.frame(Strain = "NB",Estimate = mest, p.value = mpval,R2 = mar) %>%
  rbind(Res,.)

Res$p.adj <- Res$p.value %>% p.adjust(method = "fdr")
Res$Significance <- NA
Res$Significance[which(Res$p.value < 0.05)] <- "p < 0.05"

Res$R2[which(Res$R2 < 0)] <- 0

## Append the PI mean information ##
df_pi <- Dat_a$df_pi
df_pi$Norm_PI_permeability <- df_pi$Norm_PI_permeability %>% scale

#Scale the area so the values can be compared to agar
mgenos <- df_pi$Names %>% unique %>%
  grep(pattern = "NB",invert = T,value = T)
df_col0 <- df_pi %>% subset(Names == "NB")

mcol <- df_col0$Norm_PI_permeability %>% mean
Res_11 <- NULL
for(mg in mgenos){
  mtemp_11 <- df_pi %>% subset(Names == mg) %>% droplevels
  #Levene test
  df_temp <- rbind(df_col0,mtemp_11)
  df_temp$Names <- df_temp$Names %>% factor %>% relevel(ref = "NB")
  m <- car::leveneTest(as.numeric(Norm_PI_permeability) ~ Names,data = df_temp)
  pval <- m$`Pr(>F)`[1]
  x <- df_col0$Norm_PI_permeability
  y <- mtemp_11$Norm_PI_permeability 
  mmean <- y %>% mean
  mest <- mmean - mcol
  if(pval < 0.05){
    m0 <- t.test(x = x,y = y,exact = F,var.equal = F)
    
  }else{
    m0 <- t.test(x = x,y = y,exact = F,var.equal = T)
    
  }
  m0 <- t.test(x = df_col0$Norm_PI_permeability,y = mtemp_11$Norm_PI_permeability)
  mpv <- m0$p.value
  Res_11 <- data.frame(Strain = mg,PIEstimate = mest,PIp.value = mpv) %>%
    rbind(Res_11,.)
}
Res_11$PIp.adj <- Res_11$PIp.value %>% p.adjust(method = "fdr")

Res <- merge(Res,Res_11, by = "Strain")

## Append the CS ligning coefficients 
Dat_vars <- readRDS(file = "../cleandata/res_41_developmentvariables.RDS")
df_tobind <- dcast(data = Dat_vars,formula = Strain~Variable,value.var = "Estimate")


#Create a variable denoting the trends
Res$Group <- "p > 0.05 & r2 < 0.2"
Res$Group[which((is.na(Res$Significance)) & (Res$R2 > 0.2))] <- "p > 0.05 & r2 < 0.6"
Res$Group[which(Res$Significance == "p < 0.05")] <- "p < 0.05 & r2 > 0.6"

Res$Group <- Res$Group %>% factor(levels = (c("p > 0.05 & r2 < 0.2","p > 0.05 & r2 < 0.6","p < 0.05 & r2 > 0.6")))

colnames(df_tobind)[2:11] <- paste0("Estimate_",colnames(df_tobind)[2:11])

Res <- merge(Res,df_tobind, by = "Strain") 

##Append the p.value
df_tobind <- dcast(data = Dat_vars,formula = Strain~Variable,value.var = "Estimate")
colnames(df_tobind)[2:11] <- paste0("padj_",colnames(df_tobind)[2:11])

Res <- merge(Res,df_tobind, by = "Strain") 

#Recolor the p value
p2 <- ggplot(data = Res %>% subset(Strain != "NB"),aes(p.value,R2)) +
  geom_vline(xintercept = 0.05,color = "#D9D9D9",size = 1) +
  geom_hline(yintercept = 0.6,color = "#D9D9D9",size = 1) +
  geom_point(shape = 21,aes(fill = Group),size = 3) +
  geom_text_repel(aes(label = Strain,color = Group)) +
  theme_ohchibi(font_family = "Helvetica") +
  scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.1))  +
  scale_fill_manual(values = c("#E68A65","#E6C54E","#59B0E6") %>% rev) +
  scale_color_manual(values = c("#E68A65","#E6C54E","#59B0E6") %>% rev) +
  ylab(label = "Variance explained (R2)") 

oh.save.pdf(p = p2,outname = "resub_lm_ligninbyroothair.pdf",outdir = "../figures/",
            width = 10,height = 8)

### Order based on teh coefficient
order_strains <- with(Res,order(-R2)) %>%
  Res[.,] %$% Strain %>% unique
Res$Strain <- Res$Strain %>% factor(levels = order_strains)
Res$Group <- Res$Group %>% factor(levels = Res$Group %>% levels %>% rev)

Res_ori <- Res
#Prepare structure for heatmap
colnames(Res)
melt_est <- melt(data = Res,measure.vars = c("Estimate_DistanceCSLignin","PIEstimate"),
     id.vars = c("Strain","Group"))

melt_padj  <- melt(data = Res,measure.vars = c("padj_DistanceCSLignin","PIp.adj"),
     id.vars = c("Strain"))


#Put both in the same context
melt_est$variable <- melt_est$variable %>% gsub(pattern = "Estimate_",replacement = "") %>%
  gsub(pattern = "Estimate",replacement = "")
melt_padj$variable <- melt_padj$variable %>% gsub(pattern = "padj_",replacement = "") %>%
  gsub(pattern = "p.adj",replacement = "")

merged <- merge(melt_est,melt_padj,c("Strain","variable"))
colnames(merged)[4:5] <- c("Estimate","p.adj")

merged$Significance <- NA
merged$Significance[which(merged$p.adj < 0.1)] <- "q < 0.1"


merged$Group <- match(merged$Strain,Res_ori$Strain) %>%
  Res_ori$Group[.]

merged <- merged %>% subset(Strain != "CDEF") %>% droplevels
merged$GroupMin <- "Development independent"
merged$GroupMin[merged$Group %>% grep(pattern = "p < 0.05 & r2 > 0.6")] <- "Development dependent"
merged$GroupMin <- merged$GroupMin %>% factor(levels = c("Development dependent","Development independent"))

merged$Strain <- merged$Strain %>% factor(levels = merged$Strain %>% levels %>% rev)

merged_chosen <- merged %>% 
  subset(GroupMin == "Development independent" & variable != "Suberin") %>% droplevels
merged_chosen$Estimate %>% sort %>% plot

##Append the original PI information
temp <- read.table(file = "../cleandata/res_merged_propidium_suberin_zscore.tsv",header = T)
temp <- temp[,c("Strains","Propidium_est","Propidium_p.adj")]

df_ori_est <- melt(data = temp,measure.vars = "Propidium_est",id.vars = "Strains")
df_ori_padj <- melt(data = temp,measure.vars = "Propidium_p.adj",id.vars = "Strains")

#Append to create the heatmap
mchosen <- merged_chosen$Strain %>% levels

df_ori_est <- which(df_ori_est$Strains %in% mchosen) %>%
  df_ori_est[.,] %>% droplevels

df_ori_padj <- which(df_ori_padj$Strains %in% mchosen) %>%
  df_ori_padj[.,] %>% droplevels


mtemp <- merged_chosen %>% subset(variable == "PI") %>% droplevels

df_a <- match(df_ori_est$Strains,mtemp$Strain) %>%
  mtemp[.,c("Group","Significance","GroupMin")]
df_ori_est <- cbind(df_ori_est,df_a)

df_b <- match(df_ori_padj$Strains,mtemp$Strain) %>%
  mtemp[.,c("Group","Significance","GroupMin")]
df_ori_padj <- cbind(df_ori_padj,df_b)


#Put all of them in the same context
df_ori_est$p.adj <- match(df_ori_est$Strains,df_ori_padj$Strains) %>%
  df_ori_padj$value[.]

colnames(df_ori_est)[3] <- "Estimate"
df_ori_est <- df_ori_est[,c("Strains","variable","Group","Estimate","p.adj","Significance","GroupMin")]
colnames(df_ori_est)[1] <- "Strain"
df_ori_est$Significance <- NA
df_ori_est$Significance[which(df_ori_est$p.adj < 0.1)] <- "q < 0.1"


merged_chosen <- rbind(merged_chosen,df_ori_est)


p3 <- ggplot(data = merged_chosen %>% subset(variable != "PI") %>% droplevels,aes(Strain,variable)) +
  geom_raster(aes(fill = Estimate)) +
  geom_tile(aes(color = Significance),fill = '#00000000', size = 0.5,width = 0.85,height = 0.85) + 
  scale_fill_paletteer_c("pals::kovesi.diverging_bwr_55_98_c37",
                         limits = c(-2,2),oob = squish) +
  geom_text(aes(label = round(Estimate,2))) +
  facet_grid(GroupMin~variable,space = "free",scales = "free") +
  theme_ohchibi(font_family = "Helvetica") +
  scale_color_manual(values = "black") + coord_flip() +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(family = "Helvetica",face = "plain",size = 10),
    axis.text.y = element_text(family = "Helvetica",face = "plain",size = 15),
  )
oh.save.pdf(p = p3,outname = "resub_heatmap_ligninbyroothair.pdf",outdir = "../figures/",width = 6,height = 10)
