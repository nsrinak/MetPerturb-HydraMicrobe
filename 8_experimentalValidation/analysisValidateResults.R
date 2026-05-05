library(tidyverse)
library(ggplot2)
library(ggh4x)
library(cowplot)

comb_book <- read.csv("ammoniaValidation.csv")

comb_book <- comb_book %>%
  filter(Condition != "Standard") %>% 
  mutate(Microbe = factor(
    Microbe,
    levels = c(
      "Control",
      "Curvibacter",
      "Limnobacter",
      "Pseudomonas",
      "Curvi+Pseudo",
      "Limno+Pseudo"
    ))) %>% 
  mutate(Condition = factor(
    Condition,
    levels = c("H2O", "S-medium", "L-arginine", "NH3")
  ))


p_vitro <- comb_book %>% 
  filter(Host == "none" & Condition != "S-medium") %>%
  ggplot(aes(x = Microbe, y = Concentration))+
  geom_hline(yintercept = 0, colour = "grey", linetype = "dashed")+
  geom_boxplot(fill = "#4CC9F0")+
  facet_wrap(~ Condition, scales = "free_y")+
  facetted_pos_scales(
    y = list(
      "H2O" = scale_y_continuous(limits = c(-0.85, 2.1)),
      "L-arginine" = scale_y_continuous(limits = c(-0.1, 19)),
      "NH3" = scale_y_continuous(limits = c(-0.1, 4.5))
    ))+
  theme_bw()+
  labs(y = "Relative NH3 (μg/ml)")+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())

p_gf <- comb_book %>% 
  filter(Host == "GF") %>%
  ggplot(aes(x = Microbe, y = Concentration))+
  geom_hline(yintercept = 0, colour = "grey", linetype = "dashed")+
  geom_boxplot(fill = "#4CC9F0")+
  facet_wrap(~ Condition, scales = "free_y")+
  facetted_pos_scales(
    y = list(
      "S-medium" = scale_y_continuous(limits = c(-1.3, 1)),
      "L-arginine" = scale_y_continuous(limits = c(-0.7, 1.6)),
      "NH3" = scale_y_continuous(limits = c(-0.1, 12.2))
    ))+
  theme_bw()+
  labs(y = "Relative NH3 (μg/ml)")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank())

p <- plot_grid(p_vitro, p_gf, 
               nrow = 2, align = "v", 
               rel_heights = c(1.35,2))

#ggsave(filename = "NH3conc_p_vitro_gf_recalculate_Tim.png",
#       plot = p,
#       width = 6.25, 
#       height = 4.5, 
#       units = "in", 
#       dpi = 300)

# Statistic pairwise t-test ----

dat <- comb_book %>% filter(Host == "none" & Condition == "H2O")
pairwise.t.test(x= dat$Concentration,g = dat$Microbe, p.adjust.method = "BH" )

#data:  dat$Concentration and dat$Microbe 
#
#             Control Curvibacter Limnobacter Pseudomonas Curvi+Pseudo
#  Curvibacter  0.62    -           -           -           -           
#  Limnobacter  0.71    0.71        -           -           -           
#  Pseudomonas  0.93    0.62        0.71        -           -           
#  Curvi+Pseudo 0.71    0.71        0.93        0.71        -           
#  Limno+Pseudo 0.47    0.77        0.71        0.47        0.71        
#
#P value adjustment method: BH

dat <- comb_book %>% filter(Host == "none" & Condition == "L-arginine")
pairwise.t.test(x= dat$Concentration,g = dat$Microbe, p.adjust.method = "BH" )

#data:  dat$Concentration and dat$Microbe 
#
#               Control Curvibacter Limnobacter Pseudomonas Curvi+Pseudo
#  Curvibacter  0.80365 -           -           -           -           
#  Limnobacter  9.1e-06 0.00013     -           -           -           
#  Pseudomonas  < 2e-16 < 2e-16     < 2e-16     -           -           
#  Curvi+Pseudo < 2e-16 < 2e-16     < 2e-16     0.00463     -           
#  Limno+Pseudo < 2e-16 < 2e-16     < 2e-16     0.00014     0.29834     
#
#P value adjustment method: BH 

#Control	a
#Curvibacter	a
#Limnobacter	b
#Pseudomonas	c
#Curvi+Pseudo	d
#Limno+Pseudo	d

dat <- comb_book %>% filter(Host == "none" & Condition == "NH3")
pairwise.t.test(x= dat$Concentration,g = dat$Microbe, p.adjust.method = "BH" )

#data:  dat$Concentration and dat$Microbe 
#
#               Control Curvibacter Limnobacter Pseudomonas Curvi+Pseudo
#  Curvibacter  0.7061  -           -           -           -           
#  Limnobacter  0.0738  0.1357      -           -           -           
#  Pseudomonas  0.6366  0.9224      0.1357      -           -           
#  Curvi+Pseudo 0.1357  0.1357      0.0079      0.0738      -           
#  Limno+Pseudo 0.1357  0.3279      0.6366      0.3130      0.0206      
#
#P value adjustment method: BH 

#Control	ab
#Curvibacter	a
#Limnobacter	abc
#Pseudomonas	a
#Curvi+Pseudo	b
#Limno+Pseudo	ac

dat <- comb_book %>% filter(Host == "GF" & Condition == "S-medium")
pairwise.t.test(x= dat$Concentration,g = dat$Microbe, p.adjust.method = "BH" )

#data:  dat$Concentration and dat$Microbe 
#
#             Control Curvibacter Limnobacter Pseudomonas Curvi+Pseudo
#  Curvibacter  0.79    -           -           -           -           
#  Limnobacter  0.77    0.77        -           -           -           
#  Pseudomonas  0.77    0.88        0.80        -           -           
#  Curvi+Pseudo 0.77    0.88        0.79        0.88        -           
#  Limno+Pseudo 0.88    0.77        0.77        0.77        0.77        
#
#P value adjustment method: BH 

dat <- comb_book %>% filter(Host == "GF" & Condition == "L-arginine")
pairwise.t.test(x= dat$Concentration,g = dat$Microbe, p.adjust.method = "BH" )

#data:  dat$concentration and dat$microbe 
#
#             Control Curvibacter Limnobacter Pseudomonas Curvi+Pseudo
#  Curvibacter  0.75    -           -           -           -           
#  Limnobacter  0.68    0.75        -           -           -           
#  Pseudomonas  0.68    0.60        0.60        -           -           
#  Curvi+Pseudo 0.71    0.65        0.60        0.78        -           
#  Limno+Pseudo 0.75    0.68        0.60        0.75        0.78        
#
#P value adjustment method: BH 

dat <- comb_book %>% filter(Host == "GF" & Condition == "NH3")
pairwise.t.test(x= dat$Concentration,g = dat$Microbe, p.adjust.method = "BH" )

#data:  dat$concentration and dat$microbe 
#
#               Control Curvibacter Limnobacter Pseudomonas Curvi+Pseudo
#  Curvibacter  0.2764  -           -           -           -           
#  Limnobacter  0.0789  0.0067      -           -           -           
#  Pseudomonas  0.1539  0.7531      0.0051      -           -           
#  Curvi+Pseudo 0.3446  0.0565      0.3446      0.0287      -           
#  Limno+Pseudo 0.0565  0.0052      0.8287      0.0051      0.2791      
#
#P value adjustment method: BH

#Control	ab
#Curvibacter	a
#Limnobacter	c
#Pseudomonas	a
#Curvi+Pseudo	abc
#Limno+Pseudo	bc

#-----------
# To test host production of NH3

NH3_host <- comb_book %>% filter(Microbe == "Control" &
                                   Host %in% c("GF", "none") &
                                   Condition %in% c("H2O", "S-medium"))
t.test(x = NH3_host$Concentration, g = NH3_host$Host)

#One Sample t-test
#
#data:  NH3_host$Concentration
#t = 0.21108, df = 14, p-value = 0.8359
#alternative hypothesis: true mean is not equal to 0
#95 percent confidence interval:
#  -0.3470140  0.4227744
#sample estimates:
#  mean of x 
#0.03788017 

NH3_host <- comb_book %>% filter(Microbe == "Control" &
                                          Host %in% c("GF", "none") &
                                          Condition %in% c("L-arginine"))
t.test(x = NH3_host$Concentration, g = NH3_host$Host)

#One Sample t-test
#
#data:  NH3_host$Concentration
#t = 2.9479, df = 14, p-value = 0.01059
#alternative hypothesis: true mean is not equal to 0
#95 percent confidence interval:
#  0.1086528 0.6890222
#sample estimates:
#  mean of x 
#0.3988375 

NH3_host <- comb_book %>% filter(Microbe == "Control" &
                                   Host %in% c("GF", "none") &
                                   Condition %in% c("NH3"))

t.test(x = NH3_host$Concentration, g = NH3_host$Host)

#One Sample t-test
#
#data:  NH3_host$Concentration
#t = 5.5688, df = 14, p-value = 6.919e-05
#alternative hypothesis: true mean is not equal to 0
#95 percent confidence interval:
#  2.696173 6.073870
#sample estimates:
#  mean of x 
#4.385021 

NH3_host_forp <- comb_book %>% filter(Microbe == "Control" &
                                   Host %in% c("GF", "none") &
                                   Condition %in% c("H2O", "S-medium", "L-arginine", "NH3")) %>% 
  mutate(Host = if_else(Host=="none", Host, "Hydra")) %>% 
  mutate(condition_title = if_else(Condition %in% c("H2O", "S-medium"), "H2O&S-medium", Condition),
         Host = factor(Host,levels = c("none", "Hydra")))

p_host_nh3 <- ggplot(data = NH3_host_forp, aes(x = Host, y = Concentration))+
  geom_hline(yintercept = 0, colour = "grey", linetype = "dashed")+
  geom_boxplot(fill = "#4CC9F0", colour = "black")+
  facet_wrap(~ condition_title)+
  theme_bw()+
  labs(y = "Relative NH3 (μg/ml)", x = "Host")

#ggsave(filename = "Host_NH3conc_recalculate_Tim.png",
#       plot = p_host_nh3,
#       width = 4.167, 
#       height = 2.5, 
#       units = "in", 
#       dpi = 300)