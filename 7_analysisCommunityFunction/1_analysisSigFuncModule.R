library("tidyverse")
library("vegan")
library("ggplot2")

# Data ----
lm_metabolite_wFDR <- read.csv("../5_functionalCalculation/LM_sepPlate_metabolite_withT_FDR_newlmer.csv")

lm_metabolite_2 <- lm_metabolite_wFDR %>% 
  group_by(metabolite) %>% 
  mutate(significant = FDR < 0.05,
         count = sum(significant, na.rm = TRUE)) %>% 
  dplyr::select(metabolite, count) %>% unique()

taxa_table <- read.csv("../1_ampliconIdentification/Filtered_taxa_withPlate_chemTaxa.csv")

taxa_table <- taxa_table %>% 
  group_by(sample) %>% 
  mutate(shannon = diversity(count, index = "shannon")) %>% 
  mutate(simpson = diversity(count, index = "simpson")) %>% 
  mutate(richness = sum(ra > 0)) %>% 
  mutate(evenness = shannon/log(richness)) %>% 
  mutate(rank = rank(-ra)) %>% 
  ungroup()

chemTaxa1 <- taxa_table %>% 
  dplyr::select(metabolite, direct_parent, sub_class, class, super_class, kingdom, sample, plate.x) %>% 
  unique() %>% mutate(sample2 = sample) %>% 
  column_to_rownames(var = "sample2")

# Plot significant modules ----

clas_ls <- taxa_table %>% dplyr::select(metabolite, class, sub_class) %>% unique()

msig_byclass <- lm_metabolite_wFDR %>% 
  filter(FDR <= 0.05) %>% 
  group_by(metabolite) %>% 
  summarise(metabolite_n = n()) %>%
  left_join(.,clas_ls, by="metabolite")

p_box_msig_byclass <-ggplot(msig_byclass, aes(y=reorder(sub_class, metabolite_n, FUN = median), x=metabolite_n))+
  geom_boxplot(fill = "#4CC9F0")+
  #geom_text(aes(x = 13, y = 300, label = "*"),size = 4)+
  theme_bw()+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 6),
        axis.title.x = element_text(size = 6))+
  labs(x= "Number of\nsignificant module")

p_bar_samplesig_byclass<-ggplot(msig_byclass, aes(y=reorder(sub_class, metabolite_n, FUN = median)))+
  geom_bar(fill = "#4CC9F0", color = "black")+
  theme_bw()+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 6), 
        axis.title.x = element_text(size = 6))+
  labs(x= "Number of\nsample")

p_comb_sigmod <- plot_grid(p_box_msig_byclass, p_bar_samplesig_byclass, ncol = 2, rel_widths = c(2,0.55))

#ggsave(filename = "p_comb_sigmod.png",
#       plot = p_comb_sigmod,
#       width = 3.1, 
#       height = 3.5, 
#       units = "in", 
#       dpi = 300)

# Correlation between significant modules and changes of diversity index ----
# Shannon, richness, evenness prepare difference from control ----

cont_alpha <- taxa_table %>% 
  filter(metabolite =="Control")


alpha_div <- taxa_table %>% dplyr::select(sample, plate.x, metabolite, shannon, richness, evenness) %>% unique()

diff_alpha_df <- data.frame("metabolite" = character(),
                            "Diff_Shannon" = numeric(),
                            "Diff_Evenness" = numeric(),
                            "Diff_Richness" = numeric()
)

for (i in unique(alpha_div$metabolite)) {
  int = alpha_div %>% filter(metabolite==i)
  test_shannon = median(int$shannon)
  cont_shannon = cont_alpha %>% filter(plate.x %in% unique(int$plate.x)) %>% pull(shannon) %>% median()
  diff_shannon = test_shannon - cont_shannon
  
  test_even= median(int$evenness)
  cont_even = cont_alpha %>% filter(plate.x %in% unique(int$plate.x)) %>% pull(evenness) %>% median()
  diff_even = test_even - cont_even
  
  test_rich = median(int$richness)
  cont_rich = cont_alpha %>% filter(plate.x %in% unique(int$plate.x)) %>% pull(richness) %>% median()
  diff_rich = test_rich - cont_rich
  
  
  int_diff_df <- data.frame("metabolite" = i,
                            "Diff_Shannon" = diff_shannon,
                            "Diff_Evenness" = diff_even,
                            "Diff_Richness" = diff_rich
                            
  )
  
  diff_alpha_df <- rbind(diff_alpha_df, int_diff_df)
}

glimpse(diff_alpha_df)


# Bray curtis distance from control to significant module prepare difference from control ----
set.seed(12091994)

matrix <- taxa_table %>%
  dplyr::select(sample, asv, ra) %>%
  pivot_wider(names_from = asv, values_from = ra) %>% 
  as.data.frame()

rownames(matrix) <- matrix$sample
matrix <- matrix[,-1]
matrix <- as.matrix(matrix)

bray <- vegdist(matrix, method = "bray")
bray <- as.matrix(bray)[chemTaxa1[["sample"]],chemTaxa1[["sample"]]]

plate1 <- taxa_table %>% dplyr::select(metabolite, plate.x, sample) %>% unique()
plate2 <- taxa_table %>% dplyr::select(metabolite, plate.x, sample) %>% unique() %>% 
  dplyr::rename(plate2=plate.x) %>% 
  dplyr::rename(sample2=sample) %>% 
  dplyr::rename(metabolite2=metabolite)

bray_df <- as.data.frame(bray)
bray_df<- bray_df %>% 
  rownames_to_column(var="sample") %>% 
  pivot_longer(cols = -sample, names_to = "sample2", values_to = "bray_index") %>% 
  left_join(., plate1, by="sample") %>% 
  left_join(., plate2, by="sample2") %>% 
  glimpse()

diff_bray <- bray_df %>% 
  filter(metabolite2 == "Control") %>% 
  group_by(metabolite) %>% 
  summarise(Median_Bray_Curtis = median(bray_index)) %>% 
  drop_na()

diff_alpha_beta <- diff_alpha_df %>% 
  left_join(., lm_metabolite_2, by = "metabolite") %>% 
  left_join(., diff_bray, by="metabolite") %>% 
  pivot_longer(cols = c(Median_Bray_Curtis,Diff_Shannon,Diff_Evenness,Diff_Richness),names_to = "index", values_to = "diff_value") %>% 
  mutate(abs_diff = abs(diff_value))

diff_alpha_beta %>%
  group_by(index) %>%
  summarise(
    test = list(cor.test(count, abs_diff, method = "spearman"))
  ) %>%
  mutate(
    rho = sapply(test, function(x) x$estimate),
    p_value = sapply(test, function(x) x$p.value)
  ) %>%
  dplyr::select(index, rho, p_value)

# A tibble: 4 × 3
#index                  rho     p_value
#<chr>                <dbl>       <dbl>
#  1 Diff_Evenness       0.0955 0.0851     
#2 Diff_Richness      -0.0173 0.755      
#3 Diff_Shannon        0.0427 0.443      
#4 Median_Bray_Curtis  0.288  0.000000126

p_reg_diff_nsigmod <- ggplot(diff_alpha_beta, aes(x = abs(diff_value), y = count)) +
  geom_point(color = "#4CC9F0", alpha = 0.5, stroke = 0, size = 2) +
  geom_smooth(method = "lm", color = "black") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 8), 
        axis.title.x = element_text(size = 8),
        strip.text = element_text(size = 8))+
  labs(x = "Absolute distance", y = "Module count")+
  facet_wrap(~ index, scales = "free")


#ggsave(filename = "reg_diff_nsigmod_abs.png",
#       plot = p_reg_diff_nsigmod,
#       width = 3, 
#       height = 3, 
#       units = "in", 
#       dpi = 300)

# lm test for each variable ----

# Initialize an empty list to store results
lm_results_list <- list()

# Loop over each unique index value
for (idx in unique(diff_alpha_beta$index)) {
  # Filter data for the current index
  subset_data <- diff_alpha_beta %>% filter(index == idx)
  
  # Fit the linear model
  model <- lm(count ~ diff_value, data = subset_data)
  
  # Extract summary
  summary_model <- summary(model)
  
  # Store relevant results
  lm_results_list[[idx]] <- tibble(
    index = idx,
    estimate = summary_model$coefficients["diff_value", "Estimate"],
    std_error = summary_model$coefficients["diff_value", "Std. Error"],
    p_value = summary_model$coefficients["diff_value", "Pr(>|t|)"],
    t_value = summary_model$coefficients["diff_value", "t value"],
    r_squared = summary_model$r.squared
  )
}

# Combine all results into one dataframe

lm_results_df <- bind_rows(lm_results_list)
lm_results_df$fdr <- p.adjust(p = lm_results_df$p_value, method = "BH")

#  A tibble: 4 × 7
#  index         estimate std_error  p_value t_value r_squared      fdr
#  <chr>            <dbl>     <dbl>    <dbl>   <dbl>     <dbl>    <dbl>
#  1 med_bray       146.       14.2   1.03e-21  10.3     0.247   4.14e-21
#  2 diff_shannon   -12.7       3.76  8.29e- 4  -3.37    0.0340  1.10e- 3
#  3 diff_evenness  -77.6      15.4   8.30e- 7  -5.03    0.0723  1.66e- 6
#  4 diff_richness    0.262     0.293 3.73e- 1   0.892   0.00245 3.73e- 1