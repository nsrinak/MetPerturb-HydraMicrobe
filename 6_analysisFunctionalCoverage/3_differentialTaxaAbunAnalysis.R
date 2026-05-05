library(compositions)
library(tidyverse)
library(ggplot2)
library(zCompositions)
library(lmerTest)

func_cluster <- read.csv("Filter_taxa_link16s_cluster.csv")
func_cluster <- func_cluster %>% dplyr::select(asv, taxa_group) %>% unique()

taxa_table <- read.csv("../1_ampliconIdentification/Filtered_taxa_withPlate_chemTaxa.csv")
taxa_table <- taxa_table %>% 
  left_join(., func_cluster, by = "asv") %>% 
  filter(!is.na(taxa_group))

taxa_table_group <- taxa_table %>% 
  group_by(sample, metabolite, plate.x, taxa_group, sub_class) %>% 
  summarise(group_ra = sum(ra)) %>% 
  ungroup()

control_ra <- taxa_table_group %>% 
  filter(metabolite == "Control") %>% 
  group_by(taxa_group, plate.x) %>% 
  summarise(control_ra = median(group_ra))

df_log2FC <- taxa_table_group %>% 
  left_join(., control_ra, by = c("taxa_group", "plate.x")) %>% 
  mutate(adjusted_test_ra = ifelse(group_ra == 0, 7.406859e-05/10, group_ra)) %>% 
  mutate(adjusted_control_ra = ifelse(control_ra == 0, 7.406859e-05/10, control_ra)) %>% 
  group_by(metabolite, taxa_group) %>% 
  mutate(
    ra_diff = group_ra - control_ra,
    fold_change = group_ra / control_ra,
    log2_fc = log2(fold_change)) %>%
  mutate(adjusted_log2_fc = if_else(control_ra + group_ra == 0, 0, log2(adjusted_test_ra/adjusted_control_ra))) %>% 
  mutate(
    taxa_group = case_when(
      taxa_group == "1" ~ "Curvibacter",
      taxa_group == "2" ~ "Legionella",
      taxa_group == "3" ~ "Rheinheimera",
      taxa_group == "4" ~ "Fluviicola",
      taxa_group == "5" ~ "Pseudomonas",
      TRUE              ~ "Flavobacterium"))

p_12 <- df_log2FC %>% filter(plate.x %in% c("PM1", "PM2A")) %>% ggplot(aes(x = taxa_group, y = metabolite, fill = adjusted_log2_fc))+
  geom_tile(color = "white", size = 0.2) +
  scale_fill_gradient2(
    low = "#00BFFF",   # blue (depleted)
    mid = "white",
    high = "#FF3030",  # red (enriched)
    midpoint = 0,
    #limits = c(-max_val, max_val),
    oob = scales::squish,
    name = "Log2FC",
    guide = guide_colorbar(
      frame.colour = "black",   # 🔥 this adds the border
      frame.linewidth = 0.5
    )
  ) +
  labs(x = "Microbial group", y = "Metabolite", fill = "Log2FC")+
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 310, vjust = -0.5, hjust = 0.3, size = 9),
    #legend.position = "none", 
    legend.position = "right",
    axis.text.y = element_text(size = 7))+
  facet_wrap(~ plate.x, nrow = 1, scales = "free_y")

p_34 <- df_log2FC %>% filter(plate.x %in% c("PM3B", "PM4A")) %>% ggplot(aes(x = taxa_group, y = metabolite, fill = adjusted_log2_fc))+
  geom_tile(color = "white", size = 0.2) +
  scale_fill_gradient2(
    low = "#00BFFF",   # blue (depleted)
    mid = "white",
    high = "#FF3030",  # red (enriched)
    midpoint = 0,
    #limits = c(-max_val, max_val),
    oob = scales::squish,
    name = "Log2FC",
    guide = guide_colorbar(
      frame.colour = "black",   # 🔥 this adds the border
      frame.linewidth = 0.5
    )
  ) +
  labs(x = "Microbial group", y = "Metabolite")+
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 310, vjust = -0.5, hjust = 0.3, size = 9),
    #legend.position = "none", 
    legend.position = "right",
    axis.text.y = element_text(size = 7))+
  facet_wrap(~ plate.x, nrow = 1, scales = "free_y")

#ggsave(filename = "group_diffAbun_Plate12.png",
#       plot = p_12,
#       width = 11, 
#       height = 11, 
#       units = "in", 
#       dpi = 600)

#ggsave(filename = "group_diffAbun_Plate34.png",
#       plot = p_34,
#       width = 11, 
#       height = 11, 
#       units = "in", 
#      dpi = 600)

## use ANCOMBc separate plate ----

library(ANCOMBC)

# Collapse ASVs into taxa_groups
collapsed_data <- taxa_table %>%
  group_by(sample, taxa_group) %>%
  summarise(total_count = sum(count, na.rm = TRUE), .groups = "drop")

# Pivot to wide format (Rows = taxa_group, Columns = sample)
count_matrix_grouped <- collapsed_data %>%
  pivot_wider(names_from = sample, values_from = total_count, values_fill = 0) %>%
  column_to_rownames("taxa_group") %>%
  as.matrix()

# Prepare metadata (ensure one row per sample)
metadata <- taxa_table %>%
  dplyr::select(sample, metabolite, plate.x) %>% # add other relevant columns
  distinct(sample, .keep_all = TRUE) %>%
  column_to_rownames("sample")

# Match ordering
metadata <- metadata[colnames(count_matrix_grouped), , drop = FALSE]

# Create the object

library(TreeSummarizedExperiment)

tse_grouped <- TreeSummarizedExperiment(
  assays = list(counts = count_matrix_grouped),
  colData = metadata
)

# Ensure metadata types are correct before TSE creation
metadata$metabolite <- as.factor(metadata$metabolite)
metadata$metabolite <- relevel(metadata$metabolite, ref = "Control")
metadata$plate.x    <- as.factor(metadata$plate.x) # Crucial for random effects

# Re-create/Update TSE to ensure colData is synced
colData(tse_grouped) <- DataFrame(metadata)

# Run ANCOM-BC2
output_grouped = ancombc2(
  data = tse_grouped, 
  assay_name = "counts", 
  tax_level = NULL, 
  fix_formula = "metabolite", 
  group = "metabolite", 
  rand_formula = "(1 | plate.x)", 
  p_adj_method = "BH", 
  prv_cut = 0.1,     
  lib_cut = 0, 
  global = TRUE,  
  struc_zero = TRUE,  
  neg_lb = TRUE,      
  verbose = TRUE
)

# Convert the results to a long format for easier filtering
res_filtered <- output_grouped$res %>%
  dplyr::select(taxon, contains("metabolite")) %>% 
  pivot_longer(
    cols = -taxon,
    names_to = c(".value", "metabolite"),
    names_sep = "_metabolite"
  )

allPlate_ANCOMBc <- res_filtered %>%
  arrange(q) %>% 
  mutate(
    taxon = case_when(
      taxon == "1" ~ "Curvibacter",
      taxon == "2" ~ "Legionella",
      taxon == "3" ~ "Rheinheimera",
      taxon == "4" ~ "Fluviicola",
      taxon == "5" ~ "Pseudomonas",
      TRUE         ~ "Flavobacterium"
    ))

#write.csv(allPlate_ANCOMBc, file = "allPlate_groupAbun_ANCOMBc.csv")

allPlate_ANCOMBc <- read.csv("allPlate_groupAbun_ANCOMBc.csv")

ggplot(allPlate_ANCOMBc, aes(x = lfc, y = -log10(q))) +
  geom_point(aes(color = diff), alpha = 0.6, size = 3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  # Match the aesthetic here: 'color' instead of 'fill'
  scale_color_manual(values = c("TRUE" = "#4CC9F0", "FALSE" = "grey")) +
  theme_bw() +
  labs(x = "Log fold change",
    y = "-Log10(q)",
    color = "Significant" # Match the aesthetic here too
  )


# Filter for only significant results
sig_data <- allPlate_ANCOMBc %>%
  filter(q<0.05) %>%
  # Optional: Sort by lfc so the plot is easy to read
  arrange(lfc) %>%
  mutate(label = paste(metabolite,taxon,  sep = " - "),
         label = factor(label, levels = label))

p_group_abund_ANCOMBC <-ggplot(sig_data, aes(y = lfc, x = label, fill = lfc > 0)) +
  geom_col(colour = "black") +
  scale_fill_manual(values = c("TRUE" = "#FF3030", "FALSE" = "#4CC9F0"), 
                    labels = c("TRUE" = "Increased", "FALSE" = "Decreased")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 55, hjust = 1),
        legend.position = "none")+
  labs(y = "Log fold change",
       x = "Test condition - Microbial group")

#ggsave(filename = "p_group_abund_ANCOMBC.png",
#       plot = p_group_abund_ANCOMBC,
#       width = 5.5, 
#       height = 4, 
#       units = "in", 
#       dpi = 300)
