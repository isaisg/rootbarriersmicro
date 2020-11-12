library(ohchibi)
library(palettesPM)

setwd('/home/isai/Documents/results/rootbarriersmicro/scripts/')

distfun <- function(x,method) vegan::vegdist(x = x, method = "bray")
size_legend_text <- 45
size_axis_title <- 35
size_legend_title <- 55
legend_proportion_size <- 4

palette_variance <- paletteer_d(package = "dutchmasters",palette = "pearl_earring",11)[c(1:4,7)] %>%
  c(.,"white")
names(palette_variance) <- c("Fraction","Genotype","Stress",
                             "Fraction:Stress","typebyTissue:Genotype","Residual")

mix.colors <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", 
                "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", 
                "#6A3D9A", "#FFFF99", "#B15928")

Dat <- readRDS(file = "../cleandata/dat_rootbarriers_amplicon_bacteria_syncom_mutants.RDS")

Dat_rar <- Dat$RelativeAbundance



#General test using all the samples
Dat_sub <- Dat_rar %>% 
  subset.Dataset(subset = Genotype != "blank",drop = T,clean = T) %>%
  subset.Dataset(subset = Genotype != "noPlant",drop = T,clean = T) %>%
  subset.Dataset(subset = Treatment != "NB",drop = T,clean = T) 
  
  

#PCO
mpco <- oh.pco(Tab = Dat_sub$Tab %>% t,
               Map = Dat_sub$Map,ndim = 3,eig = T,
               distfun = distfun,id_var = "DADA2_Header")

chibi.pco (list_ohpco  = mpco,col_val = "Fraction",
           mypch = 21,size = 20,alpha=1,
           size_legend_text = size_legend_text,
           size_axis_title = size_axis_title,size_axis_line = 2,
           size_title_text = size_legend_title,
           font_family = "Arial",legend_proportion_size = 4,lines_zero = T) +
  scale_fill_fraction() +
  #theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())



#CAP
mcap <- oh.cap(Tab = Dat_sub$Tab,Map = Dat_sub$Map,
               formula = "Fraction +  Genotype  + Condition(Replica)",
               distfun = distfun,perms = 9999,sqrt = T)



p_fraction <- chibi.cap(list_ohpco = mcap,col_val = "Fraction",
                        mypch = 21,size = 20,alpha=1,comp_a = "CAP1",comp_b = "CAP2",
                        size_legend_text = size_legend_text,
                        size_axis_title = size_axis_title,size_axis_line = 2,
                        size_title_text = size_legend_title,
                        font_family = "Arial",legend_proportion_size = 4,lines_zero = T) +
  scale_fill_fraction() +
  #theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())



#Permanova
#Remove inoculum
Dat_sub <- Dat_rar %>% 
  subset.Dataset(subset = Genotype != "blank",drop = T,clean = T) %>%
  subset.Dataset(subset = Genotype != "noPlant",drop = T,clean = T) %>%
  subset.Dataset(subset = Treatment != "NB",drop = T,clean = T) %>%
  subset.Dataset(subset = Fraction != "Inoculum",drop = T,clean = T) 




Tab_bray <- distfun(t(Dat_sub$Tab))
mypermanova <- adonis(Tab_bray ~  Fraction + Genotype + Fraction:Genotype  ,
                      data = Dat_sub$Map,strata = Dat_sub$Map$Replica,permutations = 9999)
mypermanova

#Create chibi permanova
p_perm <- chibi.permanova(mypermanova = mypermanova,  
                          size_legend_text = size_legend_text,size_axis_title = 15,
                          size_axis_text = 30,size_title_text = size_legend_text,legend_proportion_size = 4)
p_perm <- p_perm$p 


mcap <- oh.cap(Tab = Dat_sub$Tab,Map = Dat_sub$Map,
               formula = "Fraction +  Genotype  + Condition(Replica)",
               distfun = distfun,perms = 9999,sqrt = T)


p_fraction <- chibi.cap(list_ohpco = mcap,col_val = "Fraction",
                        mypch = 21,size = 20,alpha=1,comp_a = "CAP1",comp_b = "CAP2",
                        size_legend_text = size_legend_text,
                        size_axis_title = size_axis_title,size_axis_line = 2,
                        size_title_text = size_legend_title,
                        font_family = "Arial",legend_proportion_size = 4,lines_zero = T) +
  scale_fill_fraction() +
  #theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

composition <- egg::ggarrange(p_perm,p_fraction, widths = c(0.05,1))


#### Root ####
Dat_sub <- Dat_rar %>% 
  subset.Dataset(subset = Genotype != "blank",drop = T,clean = T) %>%
  subset.Dataset(subset = Genotype != "noPlant",drop = T,clean = T) %>%
  subset.Dataset(subset = Treatment != "NB",drop = T,clean = T) %>%
  subset.Dataset(subset = Fraction == "Root",drop = T,clean = T) 

#Permanova
Tab_bray <- distfun(t(Dat_sub$Tab))
mypermanova <- adonis(Tab_bray ~   Genotype   ,
                      data = Dat_sub$Map,strata = Dat_sub$Map$Replica,permutations = 9999)
mypermanova

p_perm <- chibi.permanova(mypermanova = mypermanova,  
                          size_legend_text = size_legend_text,size_axis_title = 15,
                          size_axis_text = 30,size_title_text = size_legend_text,legend_proportion_size = 4)
p_perm <- p_perm$p + 
  scale_fill_manual(values = palette_variance) + xlab(label = "Term Model")   #theme(legend.position = "none")



#PCO
mpco <- oh.pco(Tab = Dat_sub$Tab %>% t,
               Map = Dat_sub$Map,ndim = 3,eig = T,
               distfun = distfun,id_var = "DADA2_Header")

chibi.pco (list_ohpco  = mpco,col_val = "Genotype",
           mypch = 21,size = 20,alpha=1,
           size_legend_text = size_legend_text,
           size_axis_title = size_axis_title,size_axis_line = 2,
           size_title_text = size_legend_title,
           font_family = "Arial",legend_proportion_size = 4,lines_zero = T) +
  #theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())



mcap <- oh.cap(Tab = Dat_sub$Tab,Map = Dat_sub$Map,
               formula = "Genotype  + Condition(Replica)",
               distfun = distfun,perms = 9999,sqrt = T)


p_stress <- chibi.cap(list_ohpco = mcap,col_val = "Genotype",
                      mypch = 21,size = 10,alpha=1,comp_a = "CAP1",comp_b = "CAP2",
                      size_legend_text = size_legend_text,
                      size_axis_title = size_axis_title,size_axis_line = 2,
                      size_title_text = size_legend_title,
                      font_family = "Arial",legend_proportion_size = 4,lines_zero = T) +
  #theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  stat_ellipse(aes(color = Genotype),size = 2)+
  scale_color_manual(values = mix.colors) 

composition <- egg::ggarrange(p_perm,p_stress, widths = c(0.05,1))
oh.save.pdf(p = composition,outname = "cap_cap1_cap2_root_colored_genotype_permanova.pdf",
            outdir = "../figures/",width = 25,height = 15)


#### Agar ####
Dat_sub <- Dat_rar %>% 
  subset.Dataset(subset = Genotype != "blank",drop = T,clean = T) %>%
  subset.Dataset(subset = Genotype != "noPlant",drop = T,clean = T) %>%
  subset.Dataset(subset = Treatment != "NB",drop = T,clean = T) %>%
  subset.Dataset(subset = Fraction == "Agar",drop = T,clean = T)

#Permanova
Tab_bray <- distfun(t(Dat_sub$Tab))
mypermanova <- adonis(Tab_bray ~   Genotype   ,
                      data = Dat_sub$Map,strata = Dat_sub$Map$Replica,permutations = 9999)
mypermanova

p_perm <- chibi.permanova(mypermanova = mypermanova,  
                          size_legend_text = size_legend_text,size_axis_title = 15,
                          size_axis_text = 30,size_title_text = size_legend_text,legend_proportion_size = 4)
p_perm <- p_perm$p + 
  scale_fill_manual(values = palette_variance) + xlab(label = "Term Model")   #theme(legend.position = "none")



mcap <- oh.cap(Tab = Dat_sub$Tab,Map = Dat_sub$Map,
               formula = "Genotype  + Condition(Replica)",
               distfun = distfun,perms = 9999,sqrt = T)


p_stress <- chibi.cap(list_ohpco = mcap,col_val = "Genotype",
                      mypch = 21,size = 10,alpha=1,comp_a = "CAP1",comp_b = "CAP2",
                      size_legend_text = size_legend_text,
                      size_axis_title = size_axis_title,size_axis_line = 2,
                      size_title_text = size_legend_title,
                      font_family = "Arial",legend_proportion_size = 4,lines_zero = T) +
  #theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  stat_ellipse(aes(color = Genotype),size = 2)+
  scale_color_manual(values = mix.colors) 

composition <- egg::ggarrange(p_perm,p_stress, widths = c(0.05,1))
oh.save.pdf(p = composition,outname = "cap_cap1_cap2_agar_colored_genotype_permanova.pdf",
            outdir = "../figures/",width = 25,height = 15)



#### Shoot ####
Dat_sub <- Dat_rar %>% 
  subset.Dataset(subset = Genotype != "blank",drop = T,clean = T) %>%
  subset.Dataset(subset = Genotype != "noPlant",drop = T,clean = T) %>%
  subset.Dataset(subset = Treatment != "NB",drop = T,clean = T) %>%
  subset.Dataset(subset = Fraction == "Shoot",drop = T,clean = T)

#Permanova
Tab_bray <- distfun(t(Dat_sub$Tab))
mypermanova <- adonis(Tab_bray ~   Genotype   ,
                      data = Dat_sub$Map,strata = Dat_sub$Map$Replica,permutations = 9999)
mypermanova

p_perm <- chibi.permanova(mypermanova = mypermanova,  
                          size_legend_text = size_legend_text,size_axis_title = 15,
                          size_axis_text = 30,size_title_text = size_legend_text,legend_proportion_size = 4)
p_perm <- p_perm$p + 
  scale_fill_manual(values = palette_variance) + xlab(label = "Term Model")   #theme(legend.position = "none")



mcap <- oh.cap(Tab = Dat_sub$Tab,Map = Dat_sub$Map,
               formula = "Genotype  + Condition(Replica)",
               distfun = distfun,perms = 9999,sqrt = T)


p_stress <- chibi.cap(list_ohpco = mcap,col_val = "Genotype",
                      mypch = 21,size = 10,alpha=1,comp_a = "CAP1",comp_b = "CAP2",
                      size_legend_text = size_legend_text,
                      size_axis_title = size_axis_title,size_axis_line = 2,
                      size_title_text = size_legend_title,
                      font_family = "Arial",legend_proportion_size = 4,lines_zero = T) +
  #theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  stat_ellipse(aes(color = Genotype),size = 2)+
  scale_color_manual(values = mix.colors) 


composition <- egg::ggarrange(p_perm,p_stress, widths = c(0.05,1))
oh.save.pdf(p = composition,outname = "cap_cap1_cap2_shoot_colored_genotype_permanova.pdf",
            outdir = "../figures/",width = 25,height = 15)

rm(list=ls())
dev.off()

