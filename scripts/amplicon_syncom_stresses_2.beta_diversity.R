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



Dat <- readRDS(file = "../cleandata/dat_rootbarriers_amplicon_bacteria_syncom_stresses.RDS")

Dat_rar <- Dat$RelativeAbundance



#General test using all the samples
Dat_sub <- Dat_rar %>% 
  subset.Dataset(subset = Fraction != "Blank",drop = T,clean = T) %>%
  subset.Dataset(subset = Bacteria != "NB",drop = T,clean = T) 

Dat_sub$Map$typebyTissue <- paste0(Dat_sub$Map$Fraction,Dat_sub$Map$Plant)
Dat_sub$Map$typebyTissue <- Dat_sub$Map$typebyTissue %>%
  gsub(pattern = "RootPlant",replacement = "Root") %>%
  gsub(pattern = "ShootPlant",replacement = "Shoot") %>%
  gsub(pattern = "InoculumInoculum",replacement = "Inoculum") %>%
  factor(levels = c("Inoculum","AgarNoPlant","AgarPlant","Root","Shoot"))


#PCO
mpco <- oh.pco(Tab = Dat_sub$Tab %>% t,
               Map = Dat_sub$Map,ndim = 3,eig = T,
               distfun = distfun,id_var = "DADA2_Header")

chibi.pco (list_ohpco  = mpco,col_val = "typebyTissue",
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
               formula = "typebyTissue +  Stress  + Condition(Experiment)",
               distfun = distfun,perms = 9999,sqrt = T)



p_fraction <- chibi.cap(list_ohpco = mcap,col_val = "typebyTissue",
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
  subset.Dataset(subset = Fraction != "Blank",drop = T,clean = T) %>%
  subset.Dataset(subset = Bacteria != "NB",drop = T,clean = T)  %>%
  subset.Dataset(subset = Fraction != "Inoculum",drop = T,clean = T) %>%
  subset.Dataset(subset = Plant != "NoPlant",drop = T,clean = T)

Dat_sub$Map$Fraction <- Dat_sub$Map$Fraction %>%
  gsub(pattern = "Agar",replacement = "AgarPlant") %>%
  factor


Tab_bray <- distfun(t(Dat_sub$Tab))
mypermanova <- adonis(Tab_bray ~  Fraction + Stress + Fraction:Stress  ,
                      data = Dat_sub$Map,strata = Dat_sub$Map$Experiment,permutations = 9999)
mypermanova

#Create chibi permanova
p_perm <- chibi.permanova(mypermanova = mypermanova,  
                          size_legend_text = size_legend_text,size_axis_title = 15,
                          size_axis_text = 30,size_title_text = size_legend_text,legend_proportion_size = 4)
p_perm <- p_perm$p + 
  scale_fill_manual(values = palette_variance) + xlab(label = "Term Model")   #theme(legend.position = "none")



mcap <- oh.cap(Tab = Dat_sub$Tab,Map = Dat_sub$Map,
               formula = "Fraction +  Stress  + Condition(Experiment)",
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
oh.save.pdf(p = composition,outname = "cap_cap1_cap2_alfractions_permanova.pdf",
            outdir = "../figures/",width = 25,height = 15)


#### Root ####
Dat_sub <- Dat_rar %>% 
  subset.Dataset(subset = Fraction == "Root",drop = T,clean = T) %>%
  subset.Dataset(subset = Bacteria != "NB",drop = T,clean = T)  %>%
  subset.Dataset(subset = Plant != "NoPlant",drop = T,clean = T)

#Permanova
Tab_bray <- distfun(t(Dat_sub$Tab))
mypermanova <- adonis(Tab_bray ~   Stress   ,
                      data = Dat_sub$Map,strata = Dat_sub$Map$Experiment,permutations = 9999)
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

chibi.pco (list_ohpco  = mpco,col_val = "Stress",
           mypch = 21,size = 20,alpha=1,
           size_legend_text = size_legend_text,
           size_axis_title = size_axis_title,size_axis_line = 2,
           size_title_text = size_legend_title,
           font_family = "Arial",legend_proportion_size = 4,lines_zero = T) +
  #theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


chibi.pco (list_ohpco  = mpco,col_val = "Experiment",
           mypch = 21,size = 20,alpha=1,
           size_legend_text = size_legend_text,
           size_axis_title = size_axis_title,size_axis_line = 2,
           size_title_text = size_legend_title,
           font_family = "Arial",legend_proportion_size = 4,lines_zero = T) +
  #theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())



mcap <- oh.cap(Tab = Dat_sub$Tab,Map = Dat_sub$Map,
               formula = "Stress  + Condition(Experiment)",
               distfun = distfun,perms = 9999,sqrt = T)


p_stress <- chibi.cap(list_ohpco = mcap,col_val = "Stress",
                      mypch = 21,size = 20,alpha=1,comp_a = "CAP1",comp_b = "CAP2",
                      size_legend_text = size_legend_text,
                      size_axis_title = size_axis_title,size_axis_line = 2,
                      size_title_text = size_legend_title,
                      font_family = "Arial",legend_proportion_size = 4,lines_zero = T) +
  #theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_fill_manual(values = Dat$paleta_stress) +
  stat_ellipse(aes(color = Stress),size = 2)+
scale_color_manual(values = Dat$paleta_stress) 
  

composition <- egg::ggarrange(p_perm,p_stress, widths = c(0.05,1))
oh.save.pdf(p = composition,outname = "cap_cap1_cap2_root_colored_stress_permanova.pdf",
            outdir = "../figures/",width = 25,height = 15)


#### Agar ####
Dat_sub <- Dat_rar %>% 
  subset.Dataset(subset = Fraction == "Agar",drop = T,clean = T) %>%
  subset.Dataset(subset = Bacteria != "NB",drop = T,clean = T)  %>%
  subset.Dataset(subset = Plant != "NoPlant",drop = T,clean = T)

#Permanova
Tab_bray <- distfun(t(Dat_sub$Tab))
mypermanova <- adonis(Tab_bray ~   Stress   ,
                      data = Dat_sub$Map,strata = Dat_sub$Map$Experiment,permutations = 9999)
mypermanova

p_perm <- chibi.permanova(mypermanova = mypermanova,  
                          size_legend_text = size_legend_text,size_axis_title = 15,
                          size_axis_text = 30,size_title_text = size_legend_text,legend_proportion_size = 4)
p_perm <- p_perm$p + 
  scale_fill_manual(values = palette_variance) + xlab(label = "Term Model")   #theme(legend.position = "none")



mcap <- oh.cap(Tab = Dat_sub$Tab,Map = Dat_sub$Map,
               formula = "Stress  + Condition(Experiment)",
               distfun = distfun,perms = 9999,sqrt = T)


p_stress <- chibi.cap(list_ohpco = mcap,col_val = "Stress",
                      mypch = 21,size = 20,alpha=1,comp_a = "CAP1",comp_b = "CAP2",
                      size_legend_text = size_legend_text,
                      size_axis_title = size_axis_title,size_axis_line = 2,
                      size_title_text = size_legend_title,
                      font_family = "Arial",legend_proportion_size = 4,lines_zero = T) +
  #theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_fill_manual(values = Dat$paleta_stress)+
  stat_ellipse(aes(color = Stress),size = 2)+
  scale_color_manual(values = Dat$paleta_stress) 

composition <- egg::ggarrange(p_perm,p_stress, widths = c(0.05,1))
oh.save.pdf(p = composition,outname = "cap_cap1_cap2_agar_colored_stress_permanova.pdf",
            outdir = "../figures/",width = 25,height = 15)



#### Shoot ####
Dat_sub <- Dat_rar %>% 
  subset.Dataset(subset = Fraction == "Shoot",drop = T,clean = T) %>%
  subset.Dataset(subset = Bacteria != "NB",drop = T,clean = T)  %>%
  subset.Dataset(subset = Plant != "NoPlant",drop = T,clean = T)

#Permanova
Tab_bray <- distfun(t(Dat_sub$Tab))
mypermanova <- adonis(Tab_bray ~   Stress   ,
                      data = Dat_sub$Map,strata = Dat_sub$Map$Experiment,permutations = 9999)
mypermanova

p_perm <- chibi.permanova(mypermanova = mypermanova,  
                          size_legend_text = size_legend_text,size_axis_title = 15,
                          size_axis_text = 30,size_title_text = size_legend_text,legend_proportion_size = 4)
p_perm <- p_perm$p + 
  scale_fill_manual(values = palette_variance) + xlab(label = "Term Model")   #theme(legend.position = "none")



mcap <- oh.cap(Tab = Dat_sub$Tab,Map = Dat_sub$Map,
               formula = "Stress  + Condition(Experiment)",
               distfun = distfun,perms = 9999,sqrt = T)


p_stress <- chibi.cap(list_ohpco = mcap,col_val = "Stress",
                      mypch = 21,size = 20,alpha=1,comp_a = "CAP1",comp_b = "CAP2",
                      size_legend_text = size_legend_text,
                      size_axis_title = size_axis_title,size_axis_line = 2,
                      size_title_text = size_legend_title,
                      font_family = "Arial",legend_proportion_size = 4,lines_zero = T) +
  #theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_fill_manual(values = Dat$paleta_stress)+
  stat_ellipse(aes(color = Stress),size = 2)+
  scale_color_manual(values = Dat$paleta_stress) 


composition <- egg::ggarrange(p_perm,p_stress, widths = c(0.05,1))
oh.save.pdf(p = composition,outname = "cap_cap1_cap2_shoot_colored_stress_permanova.pdf",
            outdir = "../figures/",width = 25,height = 15)


