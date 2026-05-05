# Analysis of alpha and beta diversity

library("tidyverse")
library("ggplot2")
library("stringr")
library("cowplot")
library("vegan")
library("RColorBrewer")
library("lmerTest")
library("pheatmap")
library("performance")

## Sub class colors ----
metabolite_colors <- c(
  "Amino acids, peptides, and analogues" = "#9B2226",   # Deep Red
  "Carbohydrates and carbohydrate conjugates" = "#F4A261",  # Warm Orange
  "Purine ribonucleotides" = "#2A9D8F",  # Teal Green
  "Carboxylic acid derivatives" = "#264653",  # Deep Navy Blue
  "Fatty acids and conjugates" = "#E9C46A",  # Golden Yellow
  "Phenethylamines" = "#A8DADC",  # Soft Cyan
  "Alcohols and polyols" = "#457B9D",  # Steel Blue
  "Non-metal sulfates" = "#E76F51",  # Burnt Orange
  "Non-metal thiosulfates" = "#8ECAE6",  # Sky Blue
  "Fatty acyl glycosides" = "#D62828",  # Rich Red
  "Benzylisoquinolines" = "#F77F00",  # Bright Orange
  "Other" = "darkgreen",
  "Amines" = "grey",
  # Other classes (randomized pastel but distinct shades)
  "Pyrimidine ribonucleotides" = "#A98467",
  "Quaternary ammonium salts" = "#7EBDC2",
  "Pyrimidines and pyrimidine derivatives" = "#A3B18A",
  "Control" = "#6D6875",
  "Non-metal tetrathionates" = "#D9ACF5",
  "Alkylthiols" = "#735D78",
  "Lipoamides" = "#FFB4A2",
  "Bile acids, alcohols and derivatives" = "#E29578",
  "Organosulfonic acids and derivatives" = "#E63946",
  "Sulfinic acids" = "#FFCDB2",
  "Benzenesulfonic acids and derivatives" = "#003049",
  "Gamma butyrolactones" = "#D00000",
  "Dicarboxylic acids and derivatives" = "#14213D",
  "1-benzopyrans" = "#9D4EDD",
  "Glycerophosphates" = "#5A189A",
  "Alpha hydroxy acids and derivatives" = "#7B2CBF",
  "Carboxylic acids" = "#3C096C",
  "Beta hydroxy acids and derivatives" = "#7209B7",
  "Pyrimidine 2'-deoxyribonucleosides" = "#560BAD",
  "Short-chain keto acids and derivatives" = "#B5179E",
  "Purine 2'-deoxyribonucleosides" = "#F72585",
  "Tricarboxylic acids and derivatives" = "#4CC9F0",
  "Monoterpenoids" = "#4895EF",
  "Alpha-keto acids and derivatives" = "#4361EE",
  "1-hydroxy-2-unsubstituted benzenoids" = "#3F37C9",
  "Carbonyl compounds" = "#3A0CA3",
  "Fatty alcohols" = "#480CA8",
  "Benzoic acids and derivatives" = "#560BAD",
  "Medium-chain hydroxy acids and derivatives" = "#7209B7",
  "Carboximidic acids" = "#B5179E",
  "Non-metal nitrites" = "#F72585",
  "Ureas" = "#F4A261",
  "Indolyl carboxylic acids and derivatives" = "#FF006E",
  "Guanidines" = "#8338EC",
  "Flavones" = "#3A86FF",
  "Purines and purine derivatives" = "#FB5607",
  "Imidazoles" = "#FFBE0B",
  "Imidazolines" = "#FF006E",
  "Hybrid peptides" = "#06D6A0",
  "Non-metal phosphates" = "#118AB2",
  "Non-metal pyrophosphates" = "#073B4C",
  "Phosphate esters" = "#FFD166"
)


# Alpha diversity - Read data ----

# Metabolite association to Shannon diversity

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
  select(metabolite, direct_parent, sub_class, class, super_class, kingdom, sample, plate.x) %>% 
  unique() %>% 
  mutate(sample2 = sample) %>% 
  column_to_rownames(var = "sample2")

sub_class <- chemTaxa1 %>% 
  select(sub_class, metabolite) %>% 
  unique()

uq_taxa <- taxa_table %>% select(asv, taxa) %>% unique() %>% mutate(unique_taxa = make.unique(taxa)) %>% select(-taxa)

# First, find the total number of unique samples in your study
n_samples <- n_distinct(taxa_table$sample)

taxaabun_summary <- taxa_table %>%
  left_join(uq_taxa, by = "asv") %>% 
  mutate(taxa = case_when(
    taxa == "s__Pseudomonas anguilliseptica" ~ "g__Pseudomonas",
    taxa == "s__Pseudomonas asplenii"        ~ "g__Pseudomonas",
    taxa == "s__Pseudomonas synxantha"       ~ "g__Pseudomonas",
    TRUE ~ taxa 
  )) %>% 
  group_by(taxa) %>% 
  summarise(
    med_abun = median(ra, na.rm = TRUE),
    max_abun = max(ra, na.rm=TRUE),
    # Count how many UNIQUE samples this taxon appears in, then divide by total study samples
    prev = n_distinct(sample[count > 0]) / n_samples
  )

#write.csv(taxaabun_summary, "taxaAbundance_summary.csv")

sum_mainTaxaAbun <- taxa_table %>%
  filter(grepl("Curvi|Pseu|Acido|Acine|Legio|Rhodo|Limno", taxa)) %>% 
  group_by(sample) %>% 
  summarise(core_ra = sum(ra))

hist(sum_mainTaxaAbun[[2]])

median(sum_mainTaxaAbun$core_ra)
#[1] 0.8163075

## Boxplot - Shannon diversity - metabolite subclass
shanon_control <- taxa_table %>% 
  filter(sub_class != "NA") %>% 
  group_by(sub_class) %>% 
  mutate(median_shannon_subclass = median(shannon)) %>% 
  ungroup() %>% 
  filter(sub_class == "Control") %>% 
  pull() %>% unique()

p_shannon_subclass <- taxa_table %>% 
  filter(sub_class != "NA") %>% 
  group_by(sub_class) %>% 
  mutate(median_shannon_subclass = median(shannon)) %>% 
  ungroup() %>% 
  ggplot(aes(x = shannon, 
             y= reorder(sub_class, median_shannon_subclass, FUN="max")))+
  geom_boxplot(fill = "#4CC9F0")+
  geom_vline(xintercept = shanon_control, colour = "firebrick2", linetype = "73", linewidth = 0.65)+
  labs(y= NULL, x = "Shannon index")+
  theme_bw()+
  theme(axis.text.y = element_blank(),
        axis.title.x = element_text(size = 9),
        axis.title.y = element_text(size = 9))

## Heatmap - Abundance - metabolite subclass

uq_taxa <- taxa_table %>% select(asv, taxa) %>% unique() %>% mutate(unique_taxa = make.unique(taxa)) %>% select(-taxa)

library(ggtext)

p_ra_subclass <- taxa_table %>% 
  filter(sub_class != "NA") %>% 
  group_by(sub_class) %>% 
  mutate(median_shannon_subclass = median(shannon)) %>% ungroup() %>% 
  left_join(.,uq_taxa, by="asv") %>% 
  group_by(unique_taxa, sub_class) %>%
  mutate(represent_ra_subclass = mean(ra)) %>% ungroup() %>% 
  group_by(unique_taxa) %>% 
  mutate(represent_ra = mean(ra)) %>% ungroup() %>% 
  ggplot(aes(x = reorder(unique_taxa, -represent_ra, FUN="min"), 
             y=reorder(sub_class, median_shannon_subclass, FUN="max"), fill = log10(100+represent_ra_subclass)))+
  geom_tile()+
  scale_fill_gradient(low = "white", high = "#548B54",
                      guide = guide_colorbar(
                        frame.colour = "black",   # 🔥 this adds the border
                        frame.linewidth = 0.5
                      ))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0, size = 6),
        legend.position = "bottom", 
        axis.text.y = element_text(size = 7),
        axis.title.x = element_text(size = 9),
        axis.title.y = element_text(size = 9))+
  labs(y= "Metabolite subclass", x = "Microbe", fill = "Mean relative/nabundance")

p_ra_subclass <- p_ra_subclass+
  scale_y_discrete(
    labels = function(x) {
      ifelse(x == "Control", "<b>Control</b>", x)
    }
  ) +
  theme(
    axis.text.y = ggtext::element_markdown(size = 7)
  )

## Heatmap - diff Abundance - metabolite subclass
taxa_table_sum <- taxa_table %>% 
  filter(sub_class != "NA") %>% 
  group_by(sub_class) %>% 
  mutate(median_shannon_subclass = median(shannon)) %>% 
  ungroup() %>% 
  left_join(.,uq_taxa, by="asv") %>%
  group_by(sub_class, unique_taxa) %>% 
  mutate(median_ra = median(ra)) %>% 
  ungroup() %>% 
  select(sub_class, median_ra, unique_taxa, median_shannon_subclass) %>% 
  unique()

control_ra <- taxa_table_sum %>% 
  filter(sub_class == "Control") %>%
  unique() %>% 
  rename(control_ra = median_ra) %>% 
  select(-median_shannon_subclass,-sub_class)

res <- taxa_table_sum %>% 
  left_join(., control_ra, by = "unique_taxa") %>% 
  mutate(adjusted_test_ra = ifelse(median_ra == 0, 7.406859e-05/10, median_ra)) %>% 
  mutate(adjusted_control_ra = ifelse(control_ra == 0, 7.406859e-05/10, control_ra)) %>% 
  mutate(
    ra_diff = median_ra - control_ra,
    fold_change = median_ra / control_ra,
    log2_fc = log2(fold_change)) %>% 
  mutate(adjusted_log2_fc = if_else(control_ra + median_ra == 0, 0, log2(adjusted_test_ra/adjusted_control_ra))) %>% 
  group_by(unique_taxa) %>% mutate(agjusted_sum_fraction = sum(adjusted_test_ra/adjusted_control_ra)) %>% ungroup()

#> min(res[res$median_ra > 0, ]$median_ra)
#[1] 7.406859e-05
#> min(res[res$control_ra > 0, ]$control_ra)
#[1] 0.0002696387

max_val <- unique(sort(abs((res$log2_fc)), decreasing = TRUE))[2]

p_diff_ra_subclass <- ggplot(res, aes(y = reorder(sub_class, median_shannon_subclass, FUN="max"), 
                                      x = reorder(unique_taxa, -agjusted_sum_fraction, FUN="max"), fill = adjusted_log2_fc)) +
  geom_tile(color = "white", size = 0.2) +
  scale_fill_gradient2(
    low = "#00BFFF",   # blue (depleted)
    mid = "white",
    high = "#FF3030",  # red (enriched)
    midpoint = 0,
    #limits = c(-max_val, max_val),
    oob = scales::squish,
    name = "log 2 FC",
    guide = guide_colorbar(
      frame.colour = "black",   # 🔥 this adds the border
      frame.linewidth = 0.5
    )
  ) +
  labs(x = "Microbe")+
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0, size = 6),
    legend.position = "none", 
    #legend.position = "bottom",
    axis.title.x = element_text(size = 9),
    axis.title.y = element_blank(),
    axis.text.y = element_blank())

p_comb_shannon_diff_ra_subclass <- plot_grid(p_ra_subclass, p_diff_ra_subclass, p_shannon_subclass, align = "h", ncol = 3, rel_widths = c(1,0.5,0.25))

#ggsave(filename = "p_comb_shannon_add_diff_ra_subclass.png",
#       plot = p_comb_shannon_diff_ra_subclass,
#       width = 7.5, 
#       height = 6, 
#       units = "in", 
#       dpi = 600)

# Plot only Curvibacter and Pseudomonas

p_diff_curvi <- res %>% filter(grepl("Curvi", unique_taxa)) %>% ggplot(aes(y = reorder(sub_class, median_shannon_subclass, FUN="max"), 
                x = reorder(unique_taxa, -agjusted_sum_fraction, FUN="max"), fill = adjusted_log2_fc)) +
  geom_tile(color = "white", size = 0.2) +
  scale_fill_gradient2(
    low = "#00BFFF",   # blue (depleted)
    mid = "white",
    high = "#FF3030",  # red (enriched)
    midpoint = 0,
    #limits = c(-max_val, max_val),
    oob = scales::squish,
    name = "log 2 FC",
    guide = guide_colorbar(
      frame.colour = "black",   # 🔥 this adds the border
      frame.linewidth = 0.5
    )
  ) +
  labs(x = NULL,
       y = "Metabolite subclass")+
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0, size = 8),
    #legend.position = "none", 
    legend.position = "bottom",
    legend.title = element_text(size = 8),
    legend.title.position = "top",
    legend.key.size = unit(0.35, "cm"),
    legend.text = element_text(size = 6),
    #axis.title.x = element_text(size = 9),
    #axis.title.y = element_blank(),
    #axis.text.y = element_blank()
    )

p_diff_pseudo <- res %>% filter(grepl("Pseu", unique_taxa)) %>% ggplot(aes(y = reorder(sub_class, median_shannon_subclass, FUN="max"), 
                                                                           x = reorder(unique_taxa, -agjusted_sum_fraction, FUN="max"), fill = adjusted_log2_fc)) +
  geom_tile(color = "white", size = 0.2) +
  scale_fill_gradient2(
    low = "#00BFFF",   # blue (depleted)
    mid = "white",
    high = "#FF3030",  # red (enriched)
    midpoint = 0,
    #limits = c(-max_val, max_val),
    oob = scales::squish,
    name = "log 2 FC",
    guide = guide_colorbar(
      frame.colour = "black",   # 🔥 this adds the border
      frame.linewidth = 0.5
    )
  ) +
  labs(x = NULL,
       y = "Metabolite subclass")+
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0, size = 8),
    #legend.position = "none", 
    legend.position = "bottom",
    legend.title = element_text(size = 8),
    legend.title.position = "top",
    legend.key.size = unit(0.35, "cm"),
    legend.text = element_text(size = 6),
    #axis.title.x = element_text(size = 9),
    axis.title.y = element_blank(),
    axis.text.y = element_blank()
  )

p_diff_other <- res %>% filter(grepl("Limno|Undi|Acido|Rho", unique_taxa)) %>% ggplot(aes(y = reorder(sub_class, median_shannon_subclass, FUN="max"), 
                                                                           x = reorder(unique_taxa, -agjusted_sum_fraction, FUN="max"), fill = adjusted_log2_fc)) +
  geom_tile(color = "white", size = 0.2) +
  scale_fill_gradient2(
    low = "#00BFFF",   # blue (depleted)
    mid = "white",
    high = "#FF3030",  # red (enriched)
    midpoint = 0,
    #limits = c(-max_val, max_val),
    oob = scales::squish,
    name = "log 2 FC",
    guide = guide_colorbar(
      frame.colour = "black",   # 🔥 this adds the border
      frame.linewidth = 0.5
    )
  ) +
  labs(x = NULL,
       y = "Metabolite subclass")+
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0, size = 8),
    #legend.position = "none", 
    legend.position = "bottom",
    legend.title = element_text(size = 8),
    legend.title.position = "top",
    legend.key.size = unit(0.35, "cm"),
    legend.text = element_text(size = 6),
    #axis.title.x = element_text(size = 9),
    axis.title.y = element_blank(),
    axis.text.y = element_blank()
  )

aligned <- align_plots(
  p_diff_curvi,
  p_diff_other,
  p_diff_pseudo,
  align = "h",   # horizontal + vertical alignment
  axis = "tblr"   # align all margins
)

p_diff_all <- plot_grid(plotlist = aligned, nrow = 1)

#ggsave(filename = "p_diff_all.png",
#       plot = p_diff_all,
#       width = 11.5, 
#       height = 8, 
#       units = "in", 
#       dpi = 600)

# Plot only Curvibacter with boxplot

taxa_table_uniquetaxa <- taxa_table %>% 
  filter(sub_class != "NA") %>% 
  group_by(sub_class) %>% 
  mutate(median_shannon_subclass = median(shannon)) %>% 
  ungroup() %>% 
  left_join(.,uq_taxa, by="asv")

control_ra_individual <- taxa_table_uniquetaxa %>% 
  filter(metabolite == "Control") %>%
  group_by(unique_taxa) %>% 
  summarise(control_ra = median(ra)) %>% 
  select(control_ra, unique_taxa)

res_uniquetaxa <- taxa_table_uniquetaxa %>% 
  left_join(., control_ra_individual, by = "unique_taxa") %>% 
  mutate(adjusted_test_ra = ifelse(ra == 0, 7.406859e-05/10, ra)) %>% 
  mutate(adjusted_control_ra = ifelse(control_ra == 0, 7.406859e-05/10, control_ra)) %>% 
  mutate(
    ra_diff = ra - control_ra,
    fold_change = ra / control_ra,
    log2_fc = log2(fold_change)) %>% 
  mutate(adjusted_log2_fc = if_else(control_ra + ra == 0, 0, log2(adjusted_test_ra/adjusted_control_ra))) %>% 
  group_by(unique_taxa) %>% mutate(agjusted_sum_fraction = sum(adjusted_test_ra/adjusted_control_ra)) %>% ungroup()

p_res_uniquetaxa <- res_uniquetaxa %>% 
  filter(grepl("Curvi|Pseudo", unique_taxa)) %>% 
  ggplot(aes(y = reorder(sub_class, median_shannon_subclass, FUN="max"), 
             x = adjusted_log2_fc)) +
  geom_boxplot(fill = "#4CC9F0")+
  geom_vline(xintercept = 0, colour = "firebrick2", linetype = "73", linewidth = 0.65)+
  labs(y= "Metabolite subclass", x = "Log2FC")+
  theme_bw()+
  theme(axis.title.x = element_text(size = 9),
        axis.title.y = element_text(size = 9),
        strip.text = element_text(size = 8)) +
  facet_wrap(~ unique_taxa, nrow = 1, scales = "free_x")

p_res_uniquetaxa <- p_res_uniquetaxa +
  scale_y_discrete(
    labels = function(x) {
      ifelse(x == "Control", "<b>Control</b>", x)
    }
  ) +
  theme(
    axis.text.y = ggtext::element_markdown(size = 8)
  )

#ggsave(filename = "p_res_uniquetaxa_cur_pseu.png",
#       plot = p_res_uniquetaxa,
#       width = 11.5, 
#       height = 8, 
#       units = "in", 
#       dpi = 600)

p_res_uniquetaxa_other <- res_uniquetaxa %>% 
  filter(grepl("Limno|Undi|Acido|Rho", unique_taxa)) %>% 
  ggplot(aes(y = reorder(sub_class, median_shannon_subclass, FUN="max"), 
             x = adjusted_log2_fc)) +
  geom_boxplot(fill = "#4CC9F0")+
  geom_vline(xintercept = 0, colour = "firebrick2", linetype = "73", linewidth = 0.65)+
  labs(y= "Metabolite subclass", x = "Log2FC")+
  theme_bw()+
  theme(axis.title.x = element_text(size = 9),
        axis.title.y = element_text(size = 9),
        strip.text = element_text(size = 8)) +
  facet_wrap(~ unique_taxa, nrow = 1, scales = "free_x")

p_res_uniquetaxa_other <- p_res_uniquetaxa_other +
  scale_y_discrete(
    labels = function(x) {
      ifelse(x == "Control", "<b>Control</b>", x)
    }
  ) +
  theme(
    axis.text.y = ggtext::element_markdown(size = 8)
  )

#ggsave(filename = "p_res_uniquetaxa_other.png",
#       plot = p_res_uniquetaxa_other,
#       width = 11.5, 
#       height = 8, 
#       units = "in", 
#       dpi = 600)

## Combine Boxplot and Heatmap above

p_comb_shannon_ra_subclass <- plot_grid(p_ra_subclass, p_shannon_subclass, align = "h", ncol = 2, rel_widths = c(1,0.25))

#ggsave(filename = "p_comb_shannon_ra_subclass.png",
#       plot = p_comb_shannon_ra_subclass,
#       width = 6.5, 
#       height = 6, 
#       units = "in", 
#       dpi = 300)

# Test alpha diversity at metabolite sub class level

df_lm <- taxa_table %>% 
  mutate(sub_class = factor(sub_class)) %>% 
  select(sub_class,plate.x, sample, shannon, metabolite, richness, simpson, evenness) %>% 
  mutate(metabolite = as.factor(metabolite)) %>% 
  unique()


df_lm[["sub_class"]] <- relevel(df_lm[["sub_class"]], "Control")

lmm_sub_class <- lmerTest::lmer(shannon ~ sub_class + (1 | plate.x), data = df_lm)

summary(lmm_sub_class)
r2(lmm_sub_class)

# R2 for Mixed Models
#
# Conditional R2: 0.087
# Marginal    R2: 0.051

un_sub_class <- df_lm %>% filter(sub_class != "Control") %>% pull(sub_class) %>% unique() %>% as.character()

res_df_lm_sub_class <- data.frame("sub_class" = as.character(),
                                  "Estimate" = as.numeric(),
                                  "Std.Error" = as.numeric(),
                                  "tvalue" = as.numeric(),
                                  "pvalue" = as.numeric(), 
                                  "var_rand" = as.numeric())

for (i in un_sub_class) {
  plates_presence <- df_lm %>% filter(sub_class == i) %>% pull(plate.x) %>% unique()
  sub_df_lm <- df_lm %>% filter(sub_class %in% c(i, "Control") & plate.x %in% plates_presence)
  if (length(plates_presence)==1) {
    mix.lm <- lm(shannon ~ sub_class, data = sub_df_lm)
    variance_rand <- -9
  } else {
    rand.eff <- shannon ~ sub_class + (1 | plate.x)
    mix.lm <- lmerTest::lmer(rand.eff, data = sub_df_lm)
    random_effects <- VarCorr(mix.lm)
    std_dev_rand <- random_effects$plate.x[1]
    variance_rand <- std_dev_rand^2
  }
  sum_mod.lm <- summary(mix.lm)
  int <- data.frame("sub_class" = i,
                    "Estimate" = sum_mod.lm$coefficients[, "Estimate"][[2]],
                    "Std.Error" = sum_mod.lm$coefficients[, "Std. Error"][[2]],
                    "tvalue" = sum_mod.lm$coefficients[, "t value"][[2]],
                    "pvalue" = sum_mod.lm$coefficients[, "Pr(>|t|)"][[2]],
                    "var_rand" = variance_rand)
  res_df_lm_sub_class <- rbind(res_df_lm_sub_class, int)
}
view(res_df_lm_sub_class)

# Multiple tests correction 

# The correction was conducted for individual plate, except a metabolite is found in several plates, 
# the metabolic condition for the plates would be together taken into account.

qvalue_list <- c()

for (j in res_df_lm_sub_class$sub_class) {
  plates_presence <- df_lm %>% filter(sub_class == j) %>% pull(plate.x) %>% unique()
  number_sub_class <- df_lm %>% filter(plate.x %in% plates_presence) %>% pull(sub_class) %>% unique()
  sub_df_q <- res_df_lm_sub_class %>% filter(sub_class %in% number_sub_class)
  bh <- p.adjust(sub_df_q$pvalue, method = "BH")
  # the rule can be improved by keeping result in datafram, if found j in metabolite column then skip.
  # but I am now lazy and just keep it as it is for now, the result does not going to change
  df_all <- data.frame("fdr"=bh)
  sub_df_q <- cbind(sub_df_q, df_all)
  bh_j <- sub_df_q %>% filter(sub_class == j) %>% pull(fdr)
  qvalue_list <- c(qvalue_list, bh_j)
}

bh_df <- data.frame("FDR"=qvalue_list)
res_df_lm_sub_class_q <- cbind(res_df_lm_sub_class, bh_df)


view(res_df_lm_sub_class_q %>% filter(pvalue <= 0.05))

# Test alpha diversity at metabolite level

df_lm[["metabolite"]] <- relevel(df_lm[["metabolite"]], "Control")

lmm_metabolite <- lmerTest::lmer(shannon ~ metabolite + (1 | plate.x), data = df_lm)

summary(lmm_metabolite)

r2(lmm_metabolite)

# R2 for Mixed Models

# Conditional R2: 0.375
# Marginal    R2: 0.347

un_met <- df_lm %>% filter(metabolite != "Control") %>% pull(metabolite) %>% unique() %>% as.character()

res_df_lm_metabolite <- data.frame("metabolite" = as.character(),
                                   "Estimate" = as.numeric(),
                                   "Std.Error" = as.numeric(),
                                   "tvalue" = as.numeric(),
                                   "pvalue" = as.numeric(), 
                                   "var_rand" = as.numeric())

for (i in un_met) {
  plates_presence <- df_lm %>% filter(metabolite == i) %>% pull(plate.x) %>% unique()
  sub_df_lm <- df_lm %>% filter(metabolite %in% c(i, "Control") & plate.x %in% plates_presence)
  if (length(plates_presence)==1) {
    mix.lm <- lm(shannon ~ metabolite, data = sub_df_lm)
    variance_rand <- -9
  } else {
    rand.eff <- shannon ~ metabolite + (1 | plate.x)
    mix.lm <- lmerTest::lmer(rand.eff, data = sub_df_lm)
    random_effects <- VarCorr(mix.lm)
    std_dev_rand <- random_effects$plate.x[1]
    variance_rand <- std_dev_rand^2
  }
  sum_mod.lm <- summary(mix.lm)
  int <- data.frame("metabolite" = i,
                    "Estimate" = sum_mod.lm$coefficients[, "Estimate"][[2]],
                    "Std.Error" = sum_mod.lm$coefficients[, "Std. Error"][[2]],
                    "tvalue" = sum_mod.lm$coefficients[, "t value"][[2]],
                    "pvalue" = sum_mod.lm$coefficients[, "Pr(>|t|)"][[2]],
                    "var_rand" = variance_rand)
  res_df_lm_metabolite <- rbind(res_df_lm_metabolite, int)
}
view(res_df_lm_metabolite)

# Multiple tests correction 

# The correction was conducted for individual plate, except a metabolite is found in several plates, 
# the metabolic condition for the plates would be together taken into account.

qvalue_list <- c()

for (j in res_df_lm_metabolite$metabolite) {
  plates_presence <- df_lm %>% filter(metabolite == j) %>% pull(plate.x) %>% unique()
  number_metabolite <- df_lm %>% filter(plate.x %in% plates_presence) %>% pull(metabolite) %>% unique()
  sub_df_q <- res_df_lm_metabolite %>% filter(metabolite %in% number_metabolite)
  bh <- p.adjust(sub_df_q$pvalue, method = "BH")
  # the rule can be improved by keeping result in datafram, if found j in metabolite column then skip.
  # but I am now lazy and just keep it as it is for now, the result does not going to change
  df_all <- data.frame("fdr"=bh)
  sub_df_q <- cbind(sub_df_q, df_all)
  bh_j <- sub_df_q %>% filter(metabolite == j) %>% pull(fdr)
  qvalue_list <- c(qvalue_list, bh_j)
}

bh_df <- data.frame("FDR"=qvalue_list)
res_df_lm_metabolite_q <- cbind(res_df_lm_metabolite, bh_df)

# Plot significant subclass and metabolite

p_bar_subclass <- res_df_lm_sub_class %>% 
  filter(pvalue < 0.05) %>% 
  ggplot(aes(x = Estimate, y = reorder(sub_class, Estimate, FUN="max"), fill = sub_class))+
  geom_bar(stat = "identity", colour = "black")+
  xlim(-1.1,1.1)+
  scale_fill_manual(values = metabolite_colors)+
  geom_vline(xintercept = 0)+
  labs(y= "Metabolite\nsubclass\n ")+
  theme_bw()+
  theme(legend.position = "none", 
        axis.text.y = element_text(size = 7),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 9))

p_bar_metabolite <-res_df_lm_metabolite %>% filter(pvalue< 0.05) %>% 
  left_join(sub_class, by = "metabolite") %>% 
  ggplot(aes(x = Estimate, y = reorder(metabolite, Estimate, FUN="max"), fill = sub_class))+
  geom_bar(stat = "identity", colour = "black")+
  labs(y= "Metabolite")+
  xlim(-1.1,1.1)+
  geom_vline(xintercept = 0)+
  scale_fill_manual(values = metabolite_colors)+
  theme_bw()+
  theme(legend.position = "none", 
        axis.text.y = element_text(size = 7),
        axis.title.x = element_text(size = 9),
        axis.title.y = element_text(size = 9))

p_comb_bar_sig_lm_shannon <- plot_grid(p_bar_subclass, p_bar_metabolite, nrow = 2, rel_heights = c(6, 25), align = "v", axis = "tblr")

#ggsave(filename = "p_comb_bar_sig_lm_shannon.png",
#       plot = p_comb_bar_sig_lm_shannon ,
#       width = 3, 
#       height = 3.5, 
#       units = "in", 
#       dpi = 300)

# Beta diversity
set.seed(12091994)

matrix <- taxa_table %>%
  select(sample, asv, ra) %>%
  pivot_wider(names_from = asv, values_from = ra) %>% 
  as.data.frame()

rownames(matrix) <- matrix$sample
matrix <- matrix[,-1]
matrix <- as.matrix(matrix)

bray <- vegdist(matrix, method = "bray")
bray <- as.matrix(bray)[chemTaxa1[["sample"]],chemTaxa1[["sample"]]]

# Calculate PERMANOVA pairwise comparison control to non-control at the metabolite level

bray <- as.matrix(bray)[chemTaxa1[["sample"]],chemTaxa1[["sample"]]]

groups_metabolite <- chemTaxa1 %>% filter(metabolite != "Control") %>% pull(metabolite) %>% unique()

pairwise_metabolite <- data.frame("metabolite" = character(),
                                  "sumOfSqure" = numeric(),
                                  "r2" = numeric(),
                                  "F" = numeric(),
                                  "p-value" = numeric())

# Loop through all pairs of groups
for (i in groups_metabolite) {
  # Subset data for the two groups
  subset_data <- chemTaxa1[chemTaxa1$metabolite %in% c(i, "Control"), ]
  plate <- subset_data %>% filter(metabolite == i) %>% pull(plate.x) %>% unique()  
  subset_data <- subset_data %>% filter(plate.x %in% plate) %>% mutate(plate_met = str_c(metabolite, plate.x, sep = "."))
  
  # Get the corresponding indices for the subset
  subset_indices <- chemTaxa1$sample %in% subset_data$sample
  
  # Subset the distance matrix
  bray_subset <- as.dist(as.matrix(bray)[subset_indices, subset_indices])
  
  # Perform PERMANOVA
  result <- adonis2(bray_subset ~ metabolite, data = subset_data, permutations = 50000)
  int_df <- data.frame("metabolite" = i,
                       "sumOfSqure" = result$SumOfSqs[1],
                       "r2" = result$R2[1],
                       "F" = result$F[1],
                       "p-value" = result$`Pr(>F)`[1])
  
  # Store the result
  pairwise_metabolite <- rbind(pairwise_metabolite, int_df)
  
}

pairwise_metabolite_allPlate <- pairwise_metabolite

pairwise_metabolite_allPlate$qvalue <- NA

for (j in unique(pairwise_metabolite_allPlate$metabolite)) {
  plates_presence <- df_lm %>%
    filter(metabolite == j) %>%
    pull(plate.x) %>%
    unique()
  
  related_metabolites <- df_lm %>%
    filter(plate.x %in% plates_presence) %>%
    pull(metabolite) %>%
    unique()
  
  sub_df_q <- pairwise_metabolite_allPlate %>%
    filter(metabolite %in% related_metabolites)
  
  bh <- p.adjust(sub_df_q$p.value, method = "BH")
  sub_df_q$qvalue <- bh
  
  # Update qvalue column in original data frame
  for (i in 1:nrow(sub_df_q)) {
    row_idx <- which(pairwise_metabolite_allPlate$metabolite == sub_df_q$metabolite[i] &
                       pairwise_metabolite_allPlate$p.value == sub_df_q$p.value[i])
    pairwise_metabolite_allPlate$qvalue[row_idx] <- sub_df_q$qvalue[i]
  }
}

#write.csv(pairwise_metabolite_allPlate, "pairwisePERMANOVA_metabolite_allPlate_qvalue.csv")

pairwise_metabolite_allPlate <- read.csv("pairwisePERMANOVA_metabolite_allPlate_qvalue.csv")

p_permanova_pairwise_qvalue_xy <- pairwise_metabolite_allPlate %>% 
  left_join(sub_class, by = "metabolite") %>% 
  filter(qvalue < 0.05) %>% 
  mutate(sub_class = ifelse(is.na(sub_class), "Other", sub_class)) %>% 
  ggplot(aes(y = F, x = reorder(metabolite, F, FUN = "min"), fill = sub_class))+
  geom_bar(stat = "identity", colour = "black")+
  labs(x="Metabolite")+
  scale_fill_manual(values = metabolite_colors)+
  theme_bw()+
  theme(legend.position = "none", 
        axis.text.x = element_text(size =9, angle = 270, hjust = 0, vjust = 0.3),
        axis.title.y = element_text(size = 9),
        axis.title.x = element_text(size = 9))

#ggsave(filename = "p_permanova_pairwise_qvaluexy.png",
#       plot = p_permanova_pairwise_qvalue_xy,
#       width = 2.8, 
#       height = 2.81, 
#       units = "in", 
#       dpi = 300)

## END ----
