library("tidyverse")
library("ggplot2")
library("stringr")
library("cowplot")
library("vegan")
library("RColorBrewer")
library("ggvenn")
        
lm_metabolite_wFDR <- read.csv("../5_functionalCalculation/LM_sepPlate_metabolite_withT_FDR_newlmer.csv")

taxa_table <- read.csv("../6_analysisFunctionalCoverage/Filter_taxa_link16s_cluster.csv")

p_dominant_all <- taxa_table %>% 
  group_by(metabolite, taxa_group) %>%
  summarise(group_ra = sum(avg_ra), .groups = "drop_last") %>% 
  mutate(g_max = max(group_ra)) %>% 
  ungroup() %>% 
  mutate(group = case_when(
    taxa_group == 1 ~ "Curvibacter",
    taxa_group == 2 ~ "Legionella",
    taxa_group == 3 ~ "Rheinheimera",
    taxa_group == 4 ~ "Fluviicola",
    taxa_group == 5 ~ "Pseudomonas",
    TRUE            ~ "Flavobacterium"  
  )) %>% 
  filter(g_max == group_ra) %>%
  ggplot(., aes(x = group)) +
  geom_bar(fill = "#4CC9F0", color = "black") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0, size = 10)) +
  labs(y = "Dominance count", x = "Microbial Group")

#ggsave(filename = "p_dominant_all.png",
#       plot = p_dominant_all,
#       width = 1.5, 
#       height = 3, 
#      units = "in", 
#       dpi = 300)

## Venn diagram --

all <- lm_metabolite_wFDR %>% 
  filter(FDR< 0.05) %>% 
  group_by(metabolite) %>% 
  summarise(metabolite_n = n()) %>%
  arrange(metabolite_n) %>% 
  ungroup() %>% pull(metabolite)

pseu_dom <- taxa_table %>% 
  group_by(metabolite, taxa_group) %>%
  summarise(group_ra = sum(avg_ra)) %>% 
  mutate(g_max = max(group_ra)) %>% ungroup() %>% 
  filter(g_max == group_ra & taxa_group == 5) %>% 
  pull(metabolite)

curvi_dom <- taxa_table %>% 
  group_by(metabolite, taxa_group) %>%
  summarise(group_ra = sum(avg_ra)) %>% 
  mutate(g_max = max(group_ra)) %>% ungroup() %>% 
  filter(g_max == group_ra & taxa_group == 1) %>% 
  pull(metabolite)

legio_dom <- taxa_table %>% 
  group_by(metabolite, taxa_group) %>%
  summarise(group_ra = sum(avg_ra)) %>% 
  mutate(g_max = max(group_ra)) %>% ungroup() %>% 
  filter(g_max == group_ra & taxa_group == 2) %>% 
  pull(metabolite)


dom_df <- list("All significant" = all,
              "Curvibacter" = curvi_dom,
              "Legionella" = legio_dom,
              "Pseudomonas" = pseu_dom)

p_venn_dom <- ggvenn(dom_df, 
                        show_percentage = FALSE, 
                        fill_color = c("#529DCB", "#FF0054","#FF5400","#5F0F40"),
                        set_name_size = 4.3)

#ggsave(filename = "p_venn_dom.png",
#       plot = p_venn_dom,
#       width = 3.5, 
#       height = 3, 
#       units = "in", 
#       dpi = 300)

## Sig per each dom ----

dom_group <- taxa_table %>% 
  group_by(metabolite, taxa_group) %>%
  summarise(group_ra = sum(avg_ra)) %>% 
  mutate(g_max = max(group_ra)) %>% ungroup() %>% 
  filter(g_max == group_ra)

p_numSigMod_eachDomi <- lm_metabolite_wFDR %>% 
  filter(FDR< 0.05) %>% 
  group_by(metabolite) %>% 
  summarise(metabolite_n = n()) %>% 
  left_join(., dom_group, by = "metabolite") %>% 
  mutate(group = if_else(taxa_group == 1, "Curvibacter", 
                         if_else(taxa_group == 2, "Legionella",
                                 if_else(taxa_group == 3, "Rheinheimera",
                                         if_else(taxa_group == 4, "Fluviicola",
                                                 if_else(taxa_group == 5, "Pseudomonas", "Flavobacterium")))))) %>% 
  ggplot(., aes(x = group, y = metabolite_n))+
  geom_boxplot(fill = c("#FF0054","#FF5400","#5F0F40"), color = "black")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0, size = 10))+
  labs(y = "Number of significnat module", x = "Microbial Group")


ggsave(filename = "C:/Project/1_Hydra-microbe_interaction/new_calculation_figure_1/29102025_p_numSigMod_eachDomi.png",
       plot = p_numSigMod_eachDomi,
       width = 1.5, 
       height = 3.2, 
       units = "in", 
       dpi = 300)



## Trigger point ----

dom_group_all <- taxa_table %>% 
  group_by(metabolite, taxa_group) %>%
  summarise(group_ra = sum(avg_ra)) %>% 
  mutate(g_max = max(group_ra)) %>% ungroup() 

dom_group <- taxa_table %>% 
  group_by(metabolite, taxa_group) %>%
  summarise(group_ra = sum(avg_ra)) %>% 
  mutate(g_max = max(group_ra)) %>% ungroup() %>% 
  mutate(group = case_when(
    taxa_group == 1 ~ "Curvibacter",
    taxa_group == 2 ~ "Legionella",
    taxa_group == 3 ~ "Rheinheimera",
    taxa_group == 4 ~ "Fluviicola",
    taxa_group == 5 ~ "Pseudomonas",
    TRUE            ~ "Flavobacterium"  
  )) %>% 
  filter(g_max == group_ra) %>% 
  mutate(dom_group = group) %>% 
  dplyr::select(metabolite, dom_group)

no_sig_df <- lm_metabolite_wFDR %>%
  group_by(metabolite) %>%
  summarise(metabolite_n = sum(FDR < 0.05, na.rm = TRUE)) %>% 
  left_join(., dom_group_all, by = "metabolite") %>% 
  mutate(group = case_when(
    taxa_group == 1 ~ "Curvibacter",
    taxa_group == 2 ~ "Legionella",
    taxa_group == 3 ~ "Rheinheimera",
    taxa_group == 4 ~ "Fluviicola",
    taxa_group == 5 ~ "Pseudomonas",
    TRUE            ~ "Flavobacterium"  
  )) %>% 
  mutate(sig = if_else(metabolite_n > 0, 1, 0)) %>% 
  left_join(., dom_group, by = "metabolite")
  
  
legio_sig_ra <- no_sig_df %>% filter(taxa_group == 2 & dom_group == group) %>% 
  ggplot(.,aes(x = group_ra, y = metabolite_n))+
  geom_point(color = "#4CC9F0", alpha = 0.8, stroke = 0, size = 2)+
  geom_smooth(color = "black")+
  theme_bw()+
  labs(x = "Group relative abundance", y = "Number of significant module")

pseudo_sig_ra <- no_sig_df %>% filter(taxa_group == 5 & dom_group == group) %>% 
  ggplot(.,aes(x = group_ra, y = metabolite_n))+
  geom_point(color = "#4CC9F0", alpha = 0.8, stroke = 0, size = 2)+
  geom_smooth(color = "black")+
  theme_bw()+
  labs(x = "Group relative abundance", y = "Number of significant module")

pseu_sig_df <- no_sig_df %>% filter(taxa_group == 5 & dom_group == group)

glm_pseu <- glm(sig ~ group_ra, family = binomial(link = "logit"), data = pseu_sig_df)

summary(glm_pseu)

#Call:
#  glm(formula = sig ~ group_ra, family = binomial(link = "logit"), 
#      data = pseu_sig_df)
#
#Coefficients:
#  Estimate Std. Error z value Pr(>|z|)  
#(Intercept)   -8.232      3.217  -2.559   0.0105 *
#  group_ra      17.692      7.304   2.422   0.0154 *
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for binomial family taken to be 1)
#
#Null deviance: 33.104  on 23  degrees of freedom
#Residual deviance: 18.279  on 22  degrees of freedom
#AIC: 22.279
#
#Number of Fisher Scoring iterations: 6

curv_sig_df <- no_sig_df %>% filter(taxa_group == 1 & dom_group == group)

glm_curv <- glm(sig ~ group_ra, family = binomial(link = "logit"), data = curv_sig_df)

summary(glm_curv)

#Call:
#  glm(formula = sig ~ group_ra, family = binomial(link = "logit"), 
#      data = curv_sig_df)
#
#Coefficients:
#  Estimate Std. Error z value Pr(>|z|)
#(Intercept)   -1.170      1.059  -1.105    0.269
#group_ra      -2.323      2.069  -1.123    0.261
#
#(Dispersion parameter for binomial family taken to be 1)
#
#Null deviance: 156.56  on 265  degrees of freedom
#Residual deviance: 155.29  on 264  degrees of freedom
#AIC: 159.29
#
#Number of Fisher Scoring iterations: 5

legio_sig_df <- no_sig_df %>% filter(taxa_group == 2 & dom_group == group)

glm_legio <- glm(sig ~ group_ra, family = binomial(link = "logit"), data = legio_sig_df)

summary(glm_legio)

#Call:
#  glm(formula = sig ~ group_ra, family = binomial(link = "logit"), 
#      data = legio_sig_df)
#
#Coefficients:
#  Estimate Std. Error z value Pr(>|z|)
#(Intercept)    1.563      3.349   0.467    0.641
#group_ra      -9.259      8.773  -1.055    0.291
#
#(Dispersion parameter for binomial family taken to be 1)
#
#Null deviance: 24.877  on 34  degrees of freedom
#Residual deviance: 23.706  on 33  degrees of freedom
#AIC: 27.706
#
#Number of Fisher Scoring iterations: 5

p_logis_pseu <- ggplot(pseu_sig_df, aes(x = group_ra, y = sig)) +
  geom_point(color = "#4CC9F0", alpha = 0.4, stroke = 0, size = 2) +  # show raw data
  stat_smooth(
    method = "glm",
    method.args = list(family = binomial),
    se = TRUE,           # show confidence interval
    color = "black",
    size = 1
  ) +
  labs(x = "Group relative abundance", y = "Significant module") +
  theme_bw()+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())

p_logis_curv <- ggplot(curv_sig_df, aes(x = group_ra, y = sig)) +
  geom_point(color = "#4CC9F0", alpha = 0.4, stroke = 0, size = 2) +  # show raw data
  stat_smooth(
    method = "glm",
    method.args = list(family = binomial),
    se = TRUE,           # show confidence interval
    color = "black",
    size = 1
  ) +
  labs(x = "Group relative abundance", y = "Significance") +
  theme_bw()+
  theme(axis.title.x = element_blank())

p_logis_legio <- ggplot(legio_sig_df, aes(x = group_ra, y = sig)) +
  geom_point(color = "#4CC9F0", alpha = 0.4, stroke = 0, size = 2) +  # show raw data
  stat_smooth(
    method = "glm",
    method.args = list(family = binomial),
    se = TRUE,           # show confidence interval
    color = "black",
    size = 1
  ) +
  labs(x = "Group relative abundance", y = "Significant module") +
  theme_bw()+
  theme(axis.title.y = element_blank())

p_comb_sig_dom_logis <- plot_grid(p_logis_curv,p_logis_legio, p_logis_pseu,  
                                  ncol = 3, align = "hv",
                                  rel_widths = c(1,1.1,1.1))

#ggsave(filename = "p_comb_sigdom_logis.png",
#       plot = p_comb_sig_dom_logis,
#       width = 6, 
#       height = 2.5, 
#       units = "in", 
#       dpi = 300)