library("ape")
library("ggplot2")
library("tidyverse")
library("stringr")
library("Biostrings")

#------------------------------------------------------------------------------#
#                                  1st step                                    #
#------------------------------------------------------------------------------#

linked_file <- "Linkage_biolog_16s_stats_combined.txt"
# Read resulting file from MarkerMAG (biolog 16s data to MAGs)
linked_data <- read.csv(linked_file, stringsAsFactors = FALSE)

# Extract and replace column of MarkerGene with seq_ID (biolog_16s_***)
linked_data$MarkerGene <- str_extract(linked_data$MarkerGene, "biolog_16s_\\d+")
linked_data$GenomicSeq <- sub("GenomicSeq__", "", linked_data$GenomicSeq)

# More than 1 16s could be linked to a MAGs
# The one with highest number of linkage is kept#, to remove duplicate
uq_16s <- unique(linked_data$MarkerGene)
ft_linked_data <- data.frame()
for (i in uq_16s) {
  a <- linked_data[linked_data$MarkerGene==i,]
  a <- a[a$Number==max(a$Number),]
  ft_linked_data = rbind(ft_linked_data,a, stringsAsFactors=FALSE)
}

ft_linked_data <- ft_linked_data %>% select(-Number)


#------------------------------------------------------------------------------#
#                                  2nd step                                    #
#------------------------------------------------------------------------------#

# Read the FastTree Newick tree
tree <- read.tree("2_phylogeneticTree/tree_16SBiolog_16SMatam_16SBarrnap_16SGTDB")

# Calculate the pairwise distance matrix
dist_matrix <- cophenetic(tree)

fil_dist_matrix <- dist_matrix[unlinked_16s, !colnames(dist_matrix) %in% unlinked_16s]

# Apply function to find the lowest value in each row and its column
result <- apply(fil_dist_matrix, 1, function(row) {
  min_value <- min(row)                   # Get minimum value in the row
  col_name <- names(which.min(row))       # Get the column name of the minimum value
  row_name <- rownames(fil_dist_matrix)[which(fil_dist_matrix == row[1], arr.ind = TRUE)[1]]  # Get the row name
  
  # Return a named vector with the min value, row name, and column name
  return(c(min_value = min_value, row_name = row_name, col_name = col_name))
})

# Transpose and convert the result to a data frame
result_df <- as.data.frame(t(result), stringsAsFactors = FALSE)

# Print the result
result_df <- result_df %>% 
  mutate(min_value = as.numeric(min_value), row_name = fct_reorder(row_name, min_value)) %>%
  filter(min_value <= 0.05) %>% 
  mutate(GenomicSeq = sub(".*--", "", col_name)) %>% 
  dplyr::rename(MarkerGene = row_name) %>% 
  select(-min_value, -col_name) %>% glimpse()

#------------------------------------------------------------------------------#

linked_16s_mag <- rbind(ft_linked_data, result_df)
rownames(linked_16s_mag) <- NULL

#------------------------------------------------------------------------------#

# Load the FASTA file
fasta <- readDNAStringSet("16s_after_filtering.fasta")

# Convert to dataframe
fasta_df <- data.frame(
  MarkerGene = names(fasta),   # Sequence names
  asv = as.character(fasta)  # Sequences
)

linked_16s_mag <- linked_16s_mag %>% 
  left_join(., fasta_df, by = "MarkerGene") 

write.csv(linked_16s_mag, "16s_mags_linked.csv")