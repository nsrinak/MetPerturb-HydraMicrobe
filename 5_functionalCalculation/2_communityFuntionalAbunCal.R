library("tidyverse")
library("stringr")

# Calculate community module abundance ----

# From the taxonomic profiling table, get sub-table containing only the selected 
# taxa that successfully linked to MAGs
fil_biolog <- read.csv("Filtered_taxa_withPlate_chemTaxa.csv")
lin_biolog_mag <- read.csv("16s_mags_linked.csv")
sub_taxonomic_profile <- fil_biolog %>% filter(asv %in% unique(lin_biolog_mag$asv))
allMAGs_NormMoAbun <- read.csv("allMAGs_NormMoAbun.csv")


# The functional abundance is calculated for each sample.
# The summary table contains KEGG module as row and functional abundance in each sample as column
# Variable for storing functional abundance for all samples
community_MoAbun_table <- allMAGs_NormMoAbun[3:4]
sample_unique <- unique(sub_taxonomic_profile$sample)
forIndexModuleFile<- colnames(allMAGs_NormMoAbun)

# Loop through all sample
for (i in 1:length(sample_unique)) {
  fil_sample <- sub_taxonomic_profile %>% filter(sample==sample_unique[i])
  community_MoAbun_table[[sample_unique[i]]] <- rep(0, nrow(community_MoAbun_table))
  lin_seq <- fil_sample$asv
  
  # Loop through all taxa in a sample
  # Variable for storing functional abundance for a sample
  for (k in 1:length(lin_seq)) {
    
    # Get taxa abundance
    # Adding the module sample table (whole column addition)
    bin_name <- lin_biolog_mag %>% filter(asv==lin_seq[k]) %>% select(GenomicSeq)
    
    # "-" turns to "." when it is a column name, so there need to be a slight transformation before matching
    taxa_module <- allMAGs_NormMoAbun[,which(forIndexModuleFile==gsub("-", ".", as.character(bin_name$GenomicSeq[1])))]
    taxa_abun <- fil_sample %>% filter(asv==lin_seq[k])
    a = paste(as.character(sample_unique[i]), as.character(bin_name$GenomicSeq[1]), as.character(taxa_abun$ra), sep = " ")
    print(a)
    weighted_module <- data.frame("weight_abundance"=taxa_abun$ra*taxa_module)
    community_MoAbun_table[[sample_unique[i]]] <- community_MoAbun_table[[sample_unique[i]]] + weighted_module$weight_abundance
  }
}

write.csv(community_MoAbun_table, "Community_module_abundance.csv", row.names = FALSE)
