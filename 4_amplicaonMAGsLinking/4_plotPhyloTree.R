library(ggtree)
library(tidyverse)

tree <- read.tree("2_phylogeneticTree/tree_16SBiolog_16SMatam_16SBarrnap_16SGTDB")

tree$tip.label
asv_taxa <- read.csv("../1_ampliconIdentification/Filtered_taxa_withPlate_chemTaxa.csv")
asv_taxa <- asv_taxa %>% select(asv, taxa) %>% unique() %>% mutate(unique_taxa = make.unique(taxa)) %>% select(-taxa)
fasta_file <- "16s_after_filtering.fasta"
lines <- readLines(fasta_file)

titles <- lines[grep("^>", lines)]
sequences <- lines[!grepl("^>", lines)]

# Remove ">" from titles
titles <- sub("^>", "", titles)

df_seq <- data.frame(
  title = titles,
  asv = sequences,
  stringsAsFactors = FALSE
)

df_seq <- df_seq %>% 
  left_join(.,asv_taxa, by = "asv") %>% 
  rename("source_16s" = title)

df_tip_1 <- as.data.frame(tree$tip.label)

df_tip_1 <- df_tip_1 %>%
  mutate(source_16s = str_split(`tree$tip.label`, "--", n = 2) %>% 
           sapply(`[`, 1)) %>% 
  mutate(link_genome = str_split(`tree$tip.label`, "--", n = 2) %>% 
           sapply(`[`, 2)) %>% 
  mutate(new_link_genome = if_else(
    !is.na(link_genome) & str_detect(link_genome, "^(MT|MG)"),
    str_replace(link_genome, "^.*?-([^ -]+-[^ -]+)$", "\\1"),
    link_genome
  )) %>% 
  mutate(new_link_genome = str_remove(new_link_genome, ".fa_16S_\\d+$"))%>% 
  left_join(., df_seq, by = "source_16s") %>% 
  mutate(unique_taxa = if_else(is.na(unique_taxa), source_16s, unique_taxa)) %>% 
  mutate(short_tip = if_else(is.na(new_link_genome), unique_taxa, paste(unique_taxa, new_link_genome, sep = "->")))

#-----------------------------

tree$tip.label = df_tip_1$short_tip

p_tree <- ggtree(tree, layout = "circular", branch.length = "none") +
  geom_tiplab2(size = 2, align = TRUE, linesize = 0.5) +
  theme(plot.margin = margin(60, 60, 60, 60))

ggsave(filename = "p_tree.png",
       plot = p_tree,
       width = 12, 
       height = 12, 
       units = "in", 
       dpi = 600)
