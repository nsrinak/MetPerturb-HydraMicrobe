library("tidyverse")
library("ggplot2")
library("factoextra")
library("ggrepel")
library("pheatmap")
library("DescTools")
library("vegan")
library('dendextend')
library("ComplexHeatmap")

# Link functions to microbial profile - load data ----
allMAGs_MoCov <- read.csv("../5_functionalCalculation/allMAGs_MoCov.csv", check.names = FALSE)

lin_biolog_mag <- read.csv("../4_amplicaonMAGsLinking/16s_mags_linked.csv")

taxa_table <- read.csv("../1_ampliconIdentification/Filtered_taxa_withPlate_chemTaxa.csv")

asv <- taxa_table %>% 
  select(asv,taxa) %>%  
  unique()

lin_biolog_mag <- lin_biolog_mag %>% 
  left_join(., asv, by="asv")%>% 
  dplyr::rename(genome = GenomicSeq) %>% 
  mutate(unique_taxa = make.unique(taxa))

df_cov <- allMAGs_MoCov%>% 
  select(c(lin_biolog_mag$genome, module, module_description)) %>%
  pivot_longer(-c(module, module_description), names_to = "genome", values_to = "coverage") %>% 
  left_join(., lin_biolog_mag, by = "genome") %>% glimpse()

len_module <- df_cov %>% 
  group_by(module) %>% 
  summarise(sum_cov = sum(coverage)) %>% 
  filter(sum_cov > 0) %>% 
  pull() %>%  
  length()

non_zero_module <- df_cov %>% 
  group_by(module) %>% 
  summarise(sum_cov = sum(coverage)) %>% 
  filter(sum_cov > 0) %>% 
  pull(module)

df_cov <- df_cov %>% 
  group_by(module) %>% 
  mutate(sum_cov = sum(coverage)) %>% 
  filter(sum_cov > 0) %>% 
  ungroup()

## Gene count per genome ----

dir = "../3_reconstructionMAGs/eggNOG_results/csv_file/"
listfiles <- list.files(path = dir, pattern = ".csv", full.names = TRUE)
magnames <- gsub(pattern = ".emapper.annotations.csv", replacement = "", x = basename(listfiles))

df_geneCount <- data.frame("mag" = character(),
                           "count" = numeric())

for (i in 1:length(listfiles)) {
  df_x <- read.csv(file = listfiles[i])
  
  count_gene <- df_x %>% 
    filter(!grepl("K", COG_category)) %>% 
    filter(!grepl("S", COG_category)) %>% 
    rownames() %>% length()
  df_y <- data.frame("mag" = magnames[i],
                     "count" = count_gene)
  
  df_geneCount <- rbind(df_geneCount, df_y)
  
  print(paste(magnames[i], ": ", count_gene))
}

linked <- read.csv(file = "../4_amplicaonMAGsLinking/16s_mags_linked_withName.csv")
linked$mag <- linked$GenomicSeq

p_genecount <- df_geneCount %>% 
  left_join(.,linked, by ="mag") %>% 
  filter(!is.na(GenomicSeq)) %>% 
  ggplot(aes(x = reorder(unique_taxa, count, FUN = "max"), y = count))+
  geom_bar(stat = "identity", fill = "#4CC9F0", colour = "black"  )+
  theme_bw()+
  labs(y = "Gene count")+
  theme(axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0.2, size = 7.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 9))

#ggsave(filename = "p_genecount.png",
#       plot = p_genecount,
#       width = 3, 
#       height = 3, 
#       units = "in", 
#       dpi = 300)

## Community relative abundance coverage ----

p_cummucov <- taxa_table %>% 
  filter(asv %in% lin_biolog_mag$asv) %>% 
  group_by(sample, plate.x) %>% 
  summarise(sum_ra = sum(ra)) %>% 
  ggplot(aes(x = plate.x, y = sum_ra))+
  geom_boxplot(fill = "#4CC9F0")+
  theme_bw()+
  ylim(0,1)+
  labs(y = "Relative abundance")+
  theme(axis.title.x = element_blank())

#ggsave(filename = "p_cummucov.png",
#       plot = p_cummucov,
#       width = 3, 
#       height = 2.5, 
#       units = "in", 
#       dpi = 300)

taxa_table %>% 
  filter(asv %in% lin_biolog_mag$asv) %>% 
  group_by(sample, plate.x) %>% 
  summarise(sum_ra = sum(ra), .groups = "drop") %>% 
  group_by(plate.x) %>% summarise(med = median(sum_ra))

#plate.x   med
#<chr>   <dbl>
#1 PM1     0.925
#2 PM2A    0.920
#3 PM3B    0.825
#4 PM4A    0.853

# Grouping the genome similarity by their KEGG module coverage ----
# (I) Jaccard ----
# (II) Find the optimal number of cluster, so assigning the group correctly

mat_cov <- df_cov %>% 
  select(unique_taxa,module,coverage) %>% 
  pivot_wider(id_cols = unique_taxa, names_from = module, values_from = coverage) %>% 
  rowwise() %>%
  mutate(row_sum = sum(c_across(where(is.numeric)), na.rm = TRUE)) %>%
  ungroup() %>% 
  filter(row_sum > 0) %>% 
  select(-row_sum) %>%
  column_to_rownames(var = "unique_taxa") %>% as.matrix()

jacc.mat <- vegdist(mat_cov, method = "jaccard") %>% as.matrix()
optimalCluster.wss <- fviz_nbclust(1-jacc.mat, hcut, method = "wss", met = "average", k.max = 16)

jacc.dist <- vegdist(mat_cov, method = "jaccard") %>% as.dist()
hc <- hclust(jacc.dist, method = "average")  
group_assignments <- cutree(hc, k = 6)
group_list <- split(names(group_assignments), group_assignments)

group_assignments_df <- data.frame("taxa" = names(group_assignments),
                                   "taxa_group" = unname(group_assignments))

#write.csv(group_assignments_df, file = "group_assignments_df.csv", row.names = FALSE)

dist_mat <- vegdist(mat_cov, method = "jaccard")

# NMDS with 2 dimensions
nmds <- metaMDS(jacc.mat, k = 2, trymax = 100)

# NMDS site scores (samples)
nmds_points <- as.data.frame(nmds$points)
nmds_points$SampleID <- rownames(nmds_points)
nmds_points$Group <- factor(group_assignments[rownames(nmds_points)])

my_colors <- c(
  "#FF0054",  # hot pink
  "#FF5400",  # orange flame
  "#FFBD00",  # golden yellow
  "#00F5D4",  # aqua cyan
  "#5F0F40",  # deep plum
  "#9B5DE5"   # violet
)

p_NMDS_cov_text <- ggplot(nmds_points, aes(x = MDS1, y = MDS2, color = Group)) +
  geom_point(size = 1.5) +
  #geom_text(aes(label = SampleID), vjust = -0.5, size = 3) +
  labs(title = "",
       x = "NMDS1", y = "NMDS2") +
  theme_bw() +
  geom_text_repel(aes(label = SampleID), 
                  size = 2, 
                  max.overlaps = Inf, 
                  min.segment.length = 0, 
                  force = 1.6, 
                  segment.size = 0.2,
                  segment.alpha = 0.6) +
  scale_color_manual(values = my_colors)+
  theme(legend.position = "none")

#ggsave(filename = "p_NMDS_cov_text_reduceLineThickness.png",
#       plot = p_NMDS_cov_text,
#       width = 3, 
#       height = 3, 
#       units = "in", 
#       dpi = 300)

# Statistical analysis of module coverage across group ----
# (I) Kruskal Walis test-nonparametric multiple group comparison ----

annotation <- data.frame(Cluster = factor(group_assignments))
annotation <- annotation %>% rownames_to_column(var = "unique_taxa")

df_cov <- df_cov %>% 
  select(unique_taxa, module, module_description, coverage, asv) %>% 
  mutate(Module_Des = str_c(module, module_description, sep = "_"))%>% 
  left_join(., annotation, by="unique_taxa")

module_unique <- df_cov %>% pull(Module_Des) %>% unique()

print(paste("Total number of module: ", length(module_unique), sep = ""))

# Initialize an empty data frame with appropriate column names and types

df_kw <- data.frame(
  Module_des = character(),
  chi_square = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each unique module

for (i in module_unique) {
  a <- df_cov %>%
    filter(Module_Des == i) %>%
    select(coverage, Cluster)
  
  # Perform Kruskal-Wallis test
  kw <- kruskal.test(a$coverage ~ factor(a$Cluster))
  
  # Store results as a data frame row
  ls_kw <- data.frame(
    Module_des = i,
    chi_square = as.numeric(kw$statistic),
    p_value = as.numeric(kw$p.value)
  )
  
  # Append the results to the main data frame
  df_kw <- bind_rows(df_kw, ls_kw)
}

#write.csv(x = df_kw, file = "KW_moduleCOV.csv")

df_kw <- read.csv("KW_moduleCOV.csv")

bh <- p.adjust(as.vector(df_kw$p_value), method = "BH", n=len_module)
bh <- data.frame("FDR"=bh)
df_kw <- cbind(df_kw, bh)

slct_mol <- df_kw %>% filter(FDR < 0.05) %>% pull(Module_des)
# In total, there are 202 modules, significant different with KW

# (II) Gini coefficient across the group, so calculate the representative of module coverage per group first ----

df_cov <- df_cov %>% 
  filter(Module_Des %in% slct_mol) %>% 
  group_by(Module_Des) %>%
  mutate(gini = Gini(coverage, unbiased = FALSE))

filtered_cov <- df_cov %>% 
  filter(gini > 0.5) %>%
  group_by(module) %>%
  mutate(sd = sd(coverage)) %>% ungroup() %>% 
  filter(sd > 0.1) %>% 
  select(unique_taxa, Module_Des,coverage, asv) %>% 
  filter(Module_Des %in% slct_mol)

#write.csv(x = filtered_cov, file = "filtered_moduleCOV.csv")

hc_dendro <- as.dendrogram(hc)
hc_color <- color_branches(hc_dendro, k = 6, col = my_colors[c(3, 5, 1, 2, 6, 4)])

matrix_cov <- df_cov %>% 
  filter(gini > 0.5) %>%
  group_by(module) %>%
  mutate(sd = sd(coverage)) %>% ungroup() %>% 
  filter(sd > 0.1) %>%
  select(unique_taxa, Module_Des,coverage) %>% 
  filter(Module_Des %in% slct_mol) %>% 
  pivot_wider(id_cols = Module_Des, names_from = unique_taxa, values_from = coverage) %>% 
  column_to_rownames(var="Module_Des") %>% 
  as.matrix()

original_rowname <- row.names(matrix_cov)

rename_map <- c(
  "M00879_Arginine succinyltransferase pathway, arginine => glutamate" = "AST pathway",
  "M00365_C10-C20 isoprenoid biosynthesis, archaea" = "Archaeal MVA pathway",
  "M00366_C10-C20 isoprenoid biosynthesis, plants" = "Plant MVA pathway",
  "M00367_C10-C20 isoprenoid biosynthesis, non-plant eukaryotes" = "Eukaryote MVA pathway",
  "M00579_Phosphate acetyltransferase-acetate kinase pathway, acetyl-CoA => acetate" = "Pta-AckA Pathway",
  "M00641_Multidrug resistance, efflux pump MexEF-OprN" = "MDR MexEF-OprN",
  "M00642_Multidrug resistance, efflux pump MexJK-OprM" = "MDR MexJK-OprM",
  "M00725_Cationic antimicrobial peptide (CAMP) resistance, dltABCD operon" = "CAMP dltABCD operon",
  "M00726_Cationic antimicrobial peptide (CAMP) resistance, lysyl-phosphatidylglycerol (L-PG) synthase MprF" = "CAMP L-PG synthase MprF",
  "M00744_Cationic antimicrobial peptide (CAMP) resistance, protease PgtE" = "CAMP protease PgtE"
)

clean_names <- recode(original_rowname, !!!rename_map)

clean_names <- gsub("^M\\d+_|,.*$", "", clean_names)

row.names(matrix_cov) <- clean_names

p_cheat_cov <- Heatmap(matrix_cov, 
        cluster_columns = hc_color,
        column_split = 6,
        column_title = " ",
        column_names_gp = gpar(fontsize = 8),
        column_names_side = "bottom",
        column_names_rot = -90,
        row_names_gp = gpar(fontsize = 7),
        show_row_dend = FALSE,
        col = colorRampPalette(c("white", "#A3B18A", "darkgreen", "#003300"))(30),
        border_gp = gpar(col = "black", lwd = 0.3),
        border = "black",
        rect_gp = gpar(border = "black"),
        heatmap_legend_param = list(
          title = "Module coverage",
          legend_direction = "horizontal",
          legend_width = unit(2.8, "cm"),
          border = "black"
        ))

png("p_cheat_cov_re-ylabel.png", 
    width = 4.5, height = 9, units = "in", res = 600)
draw(p_cheat_cov, heatmap_legend_side = "bottom")
dev.off()

