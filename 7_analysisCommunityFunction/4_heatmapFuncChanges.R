library("tidyverse")
library("ggplot2")
library("stringr")
library("cowplot")
library("vegan")
library("pheatmap")
library("factoextra")
library("RColorBrewer")
library(ComplexHeatmap)
library(circlize)

# Data ----

taxa_table <- read.csv("../6_analysisFunctionalCoverage/Filter_taxa_link16s_cluster.csv")

lm_metabolite_wFDR <- read.csv("../5_functionalCalculation/LM_sepPlate_metabolite_withT_FDR_newlmer.csv")

list_module <- read.csv("../5_functionalCalculation/KEGG_moduleDB/KEGGDataBase_List_Module_reformat.csv")

# Heatmap of overall changed profiles ----

fil_sig_module <- lm_metabolite_wFDR %>% 
  mutate(tvalue = if_else(FDR < 0.05, tvalue, 0)) %>% 
  pivot_wider(id_cols = module, names_from = metabolite, values_from = tvalue) %>%
  mutate(across(everything(), ~replace_na(.x, 0))) %>% 
  column_to_rownames(var = "module")

sig_mod <- lm_metabolite_wFDR %>% 
  filter(FDR < 0.05) %>% 
  dplyr::select(module) %>% pull() %>% unique()

sig_met <- lm_metabolite_wFDR %>% 
  filter(FDR < 0.05) %>%
  dplyr::select(metabolite) %>% pull() %>% unique()

test_df_sig_module <- fil_sig_module[sig_mod, sig_met] %>%
  rownames_to_column(var = "module") %>% 
  pivot_longer(cols = -module, names_to = "metabolite", values_to = "t.value") %>% 
  left_join(list_module %>% dplyr::select(-X), by = "module") %>%
  group_by(metabolite) %>% 
  filter(sum(t.value !=0)>=5) %>% ungroup() %>% 
  group_by(module) %>% 
  filter(sum(t.value !=0)>=1) %>% ungroup()

test_mat_sig_module <- test_df_sig_module %>% 
  dplyr::select(metabolite, t.value, module_description) %>% 
  pivot_wider(id_cols = module_description, names_from = metabolite, values_from = t.value) %>% 
  column_to_rownames(var = "module_description") %>% 
  as.matrix()

g1 <- taxa_table %>%
  filter(metabolite %in% colnames(test_mat_sig_module) & taxa_group == "1") %>%
  group_by(metabolite) %>% summarise(group_ra = sum(avg_ra)) %>% column_to_rownames("metabolite")
g2 <- taxa_table %>%
  filter(metabolite %in% colnames(test_mat_sig_module) & taxa_group == "2") %>%
  group_by(metabolite) %>% summarise(group_ra = sum(avg_ra)) %>% column_to_rownames("metabolite")
g3 <- taxa_table %>%
  filter(metabolite %in% colnames(test_mat_sig_module) & taxa_group == "3") %>%
  group_by(metabolite) %>% summarise(group_ra = sum(avg_ra)) %>% column_to_rownames("metabolite")
g4 <- taxa_table %>%
  filter(metabolite %in% colnames(test_mat_sig_module) & taxa_group == "4") %>%
  group_by(metabolite) %>% summarise(group_ra = sum(avg_ra)) %>% column_to_rownames("metabolite")
g5 <- taxa_table %>%
  filter(metabolite %in% colnames(test_mat_sig_module) & taxa_group == "5") %>%
  group_by(metabolite) %>% summarise(group_ra = sum(avg_ra)) %>% column_to_rownames("metabolite")
g6 <- taxa_table %>%
  filter(metabolite %in% colnames(test_mat_sig_module) & taxa_group == "6") %>%
  group_by(metabolite) %>% summarise(group_ra = sum(avg_ra)) %>% column_to_rownames("metabolite")

g1 <- g1[colnames(test_mat_sig_module), , drop = FALSE]
g2 <- g2[colnames(test_mat_sig_module), , drop = FALSE]
g3 <- g3[colnames(test_mat_sig_module), , drop = FALSE]
g4 <- g4[colnames(test_mat_sig_module), , drop = FALSE]
g5 <- g5[colnames(test_mat_sig_module), , drop = FALSE]
g6 <- g6[colnames(test_mat_sig_module), , drop = FALSE]

number_sigModule_foeSelectedMet <- fil_sig_module[sig_mod, sig_met] %>%
  rownames_to_column(var = "module") %>% 
  pivot_longer(cols = -module, names_to = "metabolite", values_to = "t.value") %>% 
  left_join(list_module %>% dplyr::select(-X), by = "module") %>%
  group_by(module_description) %>% 
  summarise(sig_number = sum(t.value != 0)) %>% 
  ungroup() %>% 
  column_to_rownames(var = "module_description") 

number_sigModule_foeSelectedMet <- number_sigModule_foeSelectedMet[rownames(test_mat_sig_module), , drop = FALSE]

meta_data <- read.csv("../1_ampliconIdentification/Filtered_taxa_withPlate_chemTaxa.csv")

abun <- read.csv("../5_functionalCalculation/Community_module_abundance.csv", check.names = FALSE)

meta_data <- meta_data %>% 
  dplyr::select(sample, super_class, class, sub_class, direct_parent, metabolite, plate.x) %>% 
  unique() 

module_CV <- abun %>%  
  pivot_longer(-c(module, module_description),names_to = "sample", values_to = "module_abun") %>% 
  left_join(., meta_data, by = "sample") %>% 
  group_by(sample) %>% 
  mutate(rank = rank(-module_abun)) %>% 
  ungroup() %>% 
  filter(module_description %in% rownames(test_mat_sig_module)) %>% 
  group_by(module_description) %>%
  summarise(CV = sd(module_abun) / mean(module_abun) * 100, .groups = "drop") %>% 
  column_to_rownames(var = "module_description")

module_CV <- module_CV[rownames(test_mat_sig_module), , drop = FALSE]

module_median <- abun %>%  
  pivot_longer(-c(module, module_description),names_to = "sample", values_to = "module_abun") %>% 
  left_join(., meta_data, by = "sample") %>% 
  group_by(sample) %>% 
  mutate(rank = rank(-module_abun)) %>% 
  ungroup() %>% 
  filter(module_description %in% rownames(test_mat_sig_module)) %>% 
  group_by(module_description) %>%
  summarise(median = median(module_abun), .groups = "drop") %>% 
  column_to_rownames(var = "module_description")

module_median <- module_median[rownames(test_mat_sig_module), , drop = FALSE]

col_an <- HeatmapAnnotation(Curvibacter = anno_barplot(g1, gp = gpar(fill = "#FF0054"), border = FALSE),
                            Pseudomonas = anno_barplot(g5, gp = gpar(fill = "#5F0F40"), border = FALSE),
                            Legionella = anno_barplot(g2, gp = gpar(fill = "#FF5400"), border = FALSE),
                            Rheinheimera = anno_barplot(g3, gp = gpar(fill = "#FFBD00"), border = FALSE),
                            Fluviicola = anno_barplot(g4, gp = gpar(fill = "#00F5D4"), border = FALSE),
                            Flavobacterium = anno_barplot(g6, gp = gpar(fill = "#9B5DE5"), border = FALSE),
                            gap = unit(1.5, "mm"), 
                            height = unit(5.5, "cm"), 
                            annotation_name_gp = gpar(fontsize = 8),
                            annotation_name_side = "left",
                            annotation_name_rot = 0)

row_an <- rowAnnotation(#Number_Significant = anno_barplot(number_sigModule_foeSelectedMet, gp = gpar(fill = "#4CC9F0"), border = FALSE),
                        Abundance_Median = anno_barplot(module_median, gp = gpar(fill = "#4CC9F0"), border = FALSE),
                        Abundance_CV = anno_barplot(module_CV, gp = gpar(fill = "#4CC9F0"), border = FALSE),
                        annotation_name_rot = -90, 
                        annotation_name_gp = gpar(fontsize = 8),
                        gap = unit(1.5, "mm"))

col_fun <- colorRamp2(c(-12, -3, 0, 3, 12), c("#1C5679", "#529DCB", "white", "#D46934", "#5F4B3B"))

p_heatmap_funAbun <- Heatmap(test_mat_sig_module,
        top_annotation = col_an,
        left_annotation = row_an,
        show_row_names = FALSE,
        show_column_names = TRUE,
        column_title = " ",
        column_names_gp = gpar(fontsize = 8),
        column_names_side = "bottom",
        column_names_rot = -90,
        show_row_dend = FALSE, 
        heatmap_legend_param = list(
          title = "t value",
          legend_direction = "horizontal",
          border = "black"), 
        clustering_method_rows = "complete", 
        clustering_distance_rows = "spearman",
        clustering_method_columns = "complete",
        clustering_distance_columns = "spearman",
        #rect_gp = gpar(border = "black"),
        border = "black", 
        cluster_rows = TRUE, 
        #split = 25, 
        cluster_columns = TRUE,
        column_split = 2,
        col = col_fun)

png("p_heatmap_funAbun_changeColor.png", 
    width = 3.75, height = 9, units = "in", res = 600)
draw(p_heatmap_funAbun , heatmap_legend_side = "bottom")
dev.off()

## Zoom in for readability all ---- 

kegg_module <- lm_metabolite_wFDR %>% select(module)

p_main = draw(p_heatmap_funAbun)
r_order = unlist(row_order(p_main))
c_order = unlist(column_order(p_main))

## Row 1-70 ----
selected_rows = r_order[1:70]
mat_zoom = p_heatmap_funAbun@matrix[selected_rows, c_order]

row_indices = match(rownames(mat_zoom), list_module$module_description)
module_ids = list_module$module[row_indices]
short_desc = sub("[,:].*| =>.*", "", rownames(mat_zoom))
final_labels = paste0(module_ids, ": ", short_desc)
final_labels[is.na(module_ids)] = short_desc[is.na(module_ids)]
rownames(mat_zoom) = final_labels

p_row_1_70 <- Heatmap(mat_zoom, 
        col = col_fun, 
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        column_names_rot = -45, 
        column_names_gp = gpar(
          fontsize = 9),
        column_title = "Row 1-70",

        row_names_gp = gpar(fontsize = 8),
        border = "black",
        heatmap_legend_param = list(
          title = "t value",
          legend_direction = "horizontal",
          border = "black"), 
        
        show_heatmap_legend = FALSE,
       )

png("p_heatmap_funAbun_row1-70.png", 
    width = 9, height = 11, units = "in", res = 600)
draw(p_row_1_70)
dev.off()

## Row 71-140 ----
selected_rows = r_order[71:140]
mat_zoom = p_heatmap_funAbun@matrix[selected_rows, c_order]

row_indices = match(rownames(mat_zoom), list_module$module_description)
module_ids = list_module$module[row_indices]
short_desc = sub("[,:].*| =>.*", "", rownames(mat_zoom))
final_labels = paste0(module_ids, ": ", short_desc)
final_labels[is.na(module_ids)] = short_desc[is.na(module_ids)]
rownames(mat_zoom) = final_labels

p_row_71_140 <- Heatmap(mat_zoom, 
                      col = col_fun, 
                      cluster_rows = FALSE, 
                      cluster_columns = FALSE,
                      column_names_rot = -45, 
                      column_names_gp = gpar(
                        fontsize = 9),
                      column_title = "Row 71-140",
                      
                      row_names_gp = gpar(fontsize = 8),
                      border = "black",
                      heatmap_legend_param = list(
                        title = "t value",
                        legend_direction = "horizontal",
                        border = "black"), 
                      
                      show_heatmap_legend = FALSE,
)

png("p_heatmap_funAbun_row71-140.png", 
    width = 9, height = 11, units = "in", res = 600)
draw(p_row_71_140)
dev.off()

## Row 141-210 ----
selected_rows = r_order[141:210]
mat_zoom = p_heatmap_funAbun@matrix[selected_rows, c_order]

row_indices = match(rownames(mat_zoom), list_module$module_description)
module_ids = list_module$module[row_indices]
short_desc = sub("[,:].*| =>.*", "", rownames(mat_zoom))
final_labels = paste0(module_ids, ": ", short_desc)
final_labels[is.na(module_ids)] = short_desc[is.na(module_ids)]
rownames(mat_zoom) = final_labels

p_row_141_210 <- Heatmap(mat_zoom, 
                        col = col_fun, 
                        cluster_rows = FALSE, 
                        cluster_columns = FALSE,
                        column_names_rot = -45, 
                        column_names_gp = gpar(
                          fontsize = 9),
                        column_title = "Row 141-210",
                        
                        row_names_gp = gpar(fontsize = 8),
                        border = "black",
                        heatmap_legend_param = list(
                          title = "t value",
                          legend_direction = "horizontal",
                          border = "black"), 
                        
                        show_heatmap_legend = FALSE,
)

png("p_heatmap_funAbun_row141-210.png", 
    width = 9, height = 11, units = "in", res = 600)
draw(p_row_141_210)
dev.off()

## Row 211-269 ----

selected_rows = r_order[211:269]
mat_zoom = p_heatmap_funAbun@matrix[selected_rows, c_order]

row_indices = match(rownames(mat_zoom), list_module$module_description)
module_ids = list_module$module[row_indices]
short_desc = sub("[,:].*| =>.*", "", rownames(mat_zoom))
final_labels = paste0(module_ids, ": ", short_desc)
final_labels[is.na(module_ids)] = short_desc[is.na(module_ids)]
rownames(mat_zoom) = final_labels

p_row_211_269 <- Heatmap(mat_zoom, 
                         col = col_fun, 
                         cluster_rows = FALSE, 
                         cluster_columns = FALSE,
                         column_names_rot = -45, 
                         column_names_gp = gpar(
                           fontsize = 9),
                         column_title = "Row 211-269",
                         
                         row_names_gp = gpar(fontsize = 8),
                         border = "black",
                         heatmap_legend_param = list(
                           title = "t value",
                           legend_direction = "horizontal",
                           border = "black"), 
                         
                         show_heatmap_legend = FALSE,
)

png("p_heatmap_funAbun_row211-269.png", 
    width = 9, height = 11, units = "in", res = 600)
draw(p_row_211_269)
dev.off()

## Zoom in to some functional modules ----

col_order <- column_order(draw(p_heatmap_funAbun))
col_order_flat <- unlist(col_order)

# Virulence module
vir_mat <- test_df_sig_module %>%
  filter(grepl("resis|aerobactin|iron", module_description, ignore.case = TRUE)) %>%
  select(metabolite, t.value, module_description) %>% 
  pivot_wider(id_cols = module_description, names_from = metabolite, values_from = t.value) %>% 
  column_to_rownames(var = "module_description") %>% 
  as.matrix()

original_rowname <- row.names(vir_mat)

rename_map <- c(
  "Multidrug resistance, repression of porin OmpF" = "MDR OmpF repression",
  "Multidrug resistance, efflux pump AcrEF-TolC" = "MDR AcrEF-TolC",
  "Multidrug resistance, efflux pump MexEF-OprN" = "MDR MexEF-OprN",
  "Multidrug resistance, efflux pump MexJK-OprM" = "MDR MexJK-OprM",
  "Multidrug resistance, efflux pump MexXY-OprM" = "MDR MexXY-OprM",
  "Multidrug resistance, efflux pump MdtEF-TolC" = "MDR MdtEF-TolC",
  "Multidrug resistance, efflux pump MexAB-OprM" = "MDR MexAB-OprM",
  "Cationic antimicrobial peptide (CAMP) resistance, dltABCD operon" = "CAMP Dlt pathway",
  "Cationic antimicrobial peptide (CAMP) resistance, lysyl-phosphatidylglycerol (L-PG) synthase MprF" = "CAMP MprF pathway",
  "Cationic antimicrobial peptide (CAMP) resistance, protease PgtE" = "CAMP PgtE resistance",
  "Multidrug resistance, efflux pump MexPQ-OpmE" = "MDR MexPQ-OpmE",
  "Multidrug resistance, efflux pump BpeEF-OprC" = "MDR BpeEF-OprC"
)

clean_names <- recode(original_rowname, !!!rename_map)
clean_names <- gsub("^M\\d+_|,.*$", "", clean_names)
row.names(vir_mat) <- clean_names

p_heatmap_funAbun_vir <- Heatmap(vir_mat, 
                                        cluster_columns = FALSE,
                                        cluster_rows = TRUE,
                                        clustering_method_rows = "complete", 
                                        clustering_distance_rows = "spearman",
                                        row_names_gp = gpar(fontsize = 8),
                                        column_order = col_order_flat,
                                        show_column_names = FALSE,
                                        show_heatmap_legend = FALSE,
                                        show_row_dend = FALSE,
                                        col = col_fun,
                                        border = "black",
                                        width = unit(5.9, "cm"))

png("p_heatmap_vir_changeColor_shortenName.png", 
    width = 4, height = 3, units = "in", res = 600)
draw(p_heatmap_funAbun_vir)
dev.off()

# Glycocalyx
gly_mat <- test_df_sig_module %>%
  filter(grepl("keratan|Glucuronate|Heparan|glycan", module_description, ignore.case = TRUE)) %>%
  select(metabolite, t.value, module_description) %>% 
  pivot_wider(id_cols = module_description, names_from = metabolite, values_from = t.value) %>% 
  column_to_rownames(var = "module_description") %>% 
  as.matrix()

original_rowname <- row.names(gly_mat)

rename_map <- c(
  "Glucuronate pathway (uronate pathway)" = "Glucuronate pathway"
)

clean_names <- recode(original_rowname, !!!rename_map)
clean_names <- gsub("^M\\d+_|,.*$", "", clean_names)
row.names(gly_mat) <- clean_names

p_heatmap_funAbun_glycocalyx <- Heatmap(gly_mat, 
                                        cluster_columns = FALSE,
                                        cluster_rows = TRUE,
                                        clustering_method_rows = "complete", 
                                        clustering_distance_rows = "spearman",
                                        row_names_gp = gpar(fontsize = 8),
                                        column_names_gp = gpar(fontsize = 8),
                                        column_names_side = "bottom",
                                        column_names_rot = -90,
                                        column_order = col_order_flat,
                                        show_column_names = TRUE,
                                        show_heatmap_legend = FALSE,
                                        show_row_dend = FALSE,
                                        col = col_fun,
                                        border = "black",
                                        width = unit(5.9, "cm"))

png("p_heatmap_glycocalyx_changeColor_shortenName.png", 
    width = 4, height = 3, units = "in", res = 600)
draw(p_heatmap_funAbun_glycocalyx , heatmap_legend_side = "right")
dev.off()

# Nitrogen 

nitro_mat <- test_df_sig_module %>%
  filter(grepl("nitrate|nitrite|NO3", module_description, ignore.case = TRUE)) %>%
  select(metabolite, t.value, module_description) %>% 
  pivot_wider(id_cols = module_description, names_from = metabolite, values_from = t.value) %>% 
  column_to_rownames(var = "module_description") %>% 
  as.matrix()

original_rowname <- row.names(nitro_mat)

clean_names <- gsub("^M\\d+_|,.*$", "", original_rowname)
row.names(nitro_mat) <- clean_names

p_heatmap_funAbun_nitro <- Heatmap(nitro_mat, 
                                 cluster_columns = FALSE,
                                 cluster_rows = TRUE,
                                 clustering_method_rows = "complete", 
                                 clustering_distance_rows = "spearman",
                                 row_names_gp = gpar(fontsize = 8),
                                 column_order = col_order_flat,
                                 show_column_names = FALSE,
                                 show_heatmap_legend = FALSE,
                                 show_row_dend = FALSE,
                                 col = col_fun,
                                 border = "black",
                                 width = unit(5.9, "cm"))

png("p_heatmap_nitro_changeColor_shortenName.png", 
    width = 4, height = 1.5, units = "in", res = 600)
draw(p_heatmap_funAbun_nitro)
dev.off()
