library(ohchibi)
library(dplyr)
library(scales)

#Merge datasets into a single structure
set.seed(130816)
setwd('/home/isai/Documents/results/rootbarriersmicro/scripts/')


### Do the suberin plot ###
Dat <- readRDS(file = "../cleandata/dat_41isolatesmonoassociation_cfu_suberin_primroot_dryweight_ionome.RDS")
df <- Dat$Suberin



### Do the comparison
nb_mean <- df %>% subset(Strains == "NB") %$% Norm_Distance_tip_cz %>% mean
Res <- NULL
mstrains <- df$Strains %>% unique %>% grep(pattern = "NB",invert = T,value = T)
for(st in mstrains){
  df_sub <- df %>% subset(Strains == st) %>% droplevels
  mbatches <- df_sub$Batch %>% as.character %>% unique
  df_nb <- which(df$Batch %in% mbatches) %>% df[.,]  %>%
    subset(Strains == "NB") %>% droplevels
  y <- df_sub$Norm_Distance_tip_cz
  y_mean <- y %>% mean
  y_df <- df_sub
  x <- df_nb$Norm_Distance_tip_cz
  mdf <- rbind(df_nb,y_df)
  #Test for homogeneity of variances
  m1 <- car::leveneTest(Norm_Distance_tip_cz ~ Strains,mdf) 
  pvalue <- m1$`Pr(>F)`[1]
  if(pvalue < 0.05){
    m1 <- t.test(x,y,exact = F,var.equal = F)
    m_w <- wilcox.test(x,y,exact = F,var.equal = F)
  }else{
    m1 <- t.test(x,y,exact = F,var.equal = T)
    m_w <- wilcox.test(x,y,exact = F,var.equal = T)
    
  }
  pval <- m1$p.value
  temp <- data.frame(Strains = st,lavene = pvalue,
                     p.value = pval,
                     p.valueW = m_w$p.value,Mean = y_mean,fc=y_mean-nb_mean
  )
  Res <- rbind(Res,temp)
  
}

plot(Res$p.value,Res$p.valueW)
cor.test(Res$p.value,Res$p.valueW)
cor.test(Res$p.value,Res$p.valueW,method = "spearman",exact = F)



ggplot(data = Res,aes(p.value,p.valueW)) + 
  geom_point() + 
  scale_x_continuous(limits = c(0,0.1)) +
  scale_y_continuous(limits = c(0,0.1)) 


#Stay with t.test p.value
Res$p.adj <- Res$p.value %>% p.adjust(method = "fdr")

Res <- rbind(Res,
             data.frame(Strains = "NB",lavene = 1,p.value = 1,p.valueW = 1,Mean=nb_mean,fc = 0,p.adj=1))

#Order by the means
Res <- with(Res,order(Mean)) %>%
  Res[.,]

order_strains <- Res$Strains %>% as.character %>% unique

#Add significance
Res$Significance <- "NotSignificant"
Res$Significance[which(Res$p.adj < 0.1)] <- "Significant"

msig <- Res %>% subset(Significance == "Significant") %$% Strains %>% as.character

df$Strains <- df$Strains %>% factor(levels = order_strains)
df$Significance <- "NotSignificant"
df$Significance[which(df$Strains %in% msig)] <- "Significant"
df$Significance <- df$Significance %>% factor
chibi.boxplot(Map = df,x_val = "Strains",y_val = "Norm_Distance_tip_cz",col_val = "Significance") +
  theme(axis.text.x = element_text(angle = 90))


#Determine interquantile of nb and cdef and plto them along the plot
#also need to change the name so they denotoe the genus name

##Try manual implementation for figure
ggplot(data =df,aes(Strains,Norm_Distance_tip_cz) ) +
  geom_sina(alpha = 0.1,size = 2) +
  stat_summary(fun.data = mean_cl_normal, geom = "pointrange",size =1) +
  theme_ohchibi() + 
  theme(
    axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1)
  ) +
  ylab(label = "Suberin")

#Check the qunatile classification
df_mean <- aggregate(Norm_Distance_tip_cz ~ Strains,df,FUN = mean)

df_mean <- with(df_mean,order(Norm_Distance_tip_cz)) %>% df_mean[.,]

#Determine the ones that are less than nb
val_nb <- df_mean %>% subset(Strains == "NB") %$% Norm_Distance_tip_cz
df_down <- df_mean %>% subset(Norm_Distance_tip_cz < val_nb)

df_up <- df_mean %>% subset(Norm_Distance_tip_cz >= val_nb)

df_up$Norm_Distance_tip_cz %>% hist

groups <- ntile(df_up$Norm_Distance_tip_cz,n = 4)
df_up$Group <- paste0("Group",groups+1)

df_down$Group <- "Group1"

df_sum <- rbind(df_up,df_down)

df$Group <- match(df$Strains,df_sum$Strains) %>% df_sum$Group[.]


##Add the genus to the figure so we can 
df <- merge(df,Dat$Map, by = "Strains") 
df$Nom <- paste0(df$Genus," ",df$Strains)
df$Nom <- df$Nom %>% gsub(pattern = "NA ",replacement = "")

order_nom <- with(df,order(Strains)) %>% df[.,] %$% Nom %>% as.character %>% unique
df$Nom <- df$Nom %>% factor(levels = order_nom)

p <- ggplot(data =df,aes(Nom,Norm_Distance_tip_cz) ) +
  geom_sina(alpha = 0.1,size = 4) +
  stat_summary(fun.data = mean_cl_normal, geom = "pointrange",size =1,aes(color = Group)) +
  stat_summary(data = df %>% subset(Significance != "Significant"),
               fun.y = mean, geom = "point",shape = 21,size =7,stroke = 0.75,aes(fill = Group)) +
  stat_summary(data = df %>% subset(Significance == "Significant"),
               fun.y = mean, geom = "point",shape = 21,size =7,stroke = 1.5,aes(fill = Group)) +
  theme_ohchibi(size_panel_border = 1) + 
  ylab(label = "Suberin") +
  facet_grid(.~Group,space = "free",scales = "free") +
  scale_color_paletteer_d(package ="LaCroixColoR",palette = "PassionFruit" ) +
  scale_fill_paletteer_d(package ="LaCroixColoR",palette = "PassionFruit" )  +
  theme(
    panel.spacing.x=unit(0.1, "lines"),
    panel.grid.major.y = element_line(size = 0.5,color = "#D9D9D9"),
    axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1),
    panel.background = element_rect(fill = "white"),
    strip.text.x = element_text(family = "Arial",face = "bold",size = 20),
    axis.title.x = element_blank()
    
  )  + 
  scale_y_continuous(breaks = seq(0,5,1),limits = c(0,5),oob = rescale_none)


oh.save.pdf(p = p,outname = "suberin_41_groups.distance.pdf",
            outdir = "../figures/",width = 17,height =10)



#save the group information
df_mean <- aggregate(Norm_Distance_tip_cz ~ Strains,df,mean)
df_mean$Group <- match(df_mean$Strains,df$Strains) %>%
  df$Group[.]


write.table(x = df_mean,file = "../cleandata/df_41_suberin_results.distance.tsv",
            append = F,quote = F,sep = "\t",row.names = F,col.names = T)


## Second aesthetic
p <- ggplot(data =df,aes(Strains,Norm_Distance_tip_cz) ) +
  geom_rect(xmin = 0,xmax =6.5,ymin =-Inf , ymax = Inf,alpha = 0.1,fill = "#E6E6E6") +
  geom_rect(xmin = 16.5,xmax =25.5,ymin =-Inf , ymax = Inf,alpha = 0.1,fill = "#E6E6E6") +
  geom_rect(xmin = 34.5,xmax =44,ymin =-Inf , ymax = Inf,alpha = 0.1,fill = "#E6E6E6") +
  geom_sina(alpha = 0.1,size = 4) +
  stat_summary(fun.data = mean_cl_normal, geom = "pointrange",size =1,aes(color = Group)) +
  stat_summary(data = df %>% subset(Significance != "Significant"),
               fun.y = mean, geom = "point",shape = 21,size =7,stroke = 0,aes(fill = Group)) +
  stat_summary(data = df %>% subset(Significance == "Significant"),
               fun.y = mean, geom = "point",shape = 21,size =7,stroke = 1.5,aes(fill = Group)) +
  theme_ohchibi(size_panel_border = 1) + 
  ylab(label = "Suberin") +
  scale_color_paletteer_d(package ="LaCroixColoR",palette = "PassionFruit" ) +
  scale_fill_paletteer_d(package ="LaCroixColoR",palette = "PassionFruit" )  +
  theme(
    panel.spacing.x=unit(0.1, "lines"),
    panel.grid.major.y = element_line(size = 0.5,color = "#D9D9D9"),
    axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1),
    panel.background = element_rect(fill = "white"),
    strip.text.x = element_text(family = "Arial",face = "bold",size = 20),
    axis.title.x = element_blank()
    
  )  + 
  scale_y_continuous(breaks = seq(0,4,1),limits = c(0,4),oob = rescale_none)

oh.save.pdf(p = p,outname = "suberin_41_groups_continous.distance.pdf",
            outdir = "../figures/",width = 17,height =10)

