library("tidyverse")

taxa_table <- read.csv("../1_ampliconIdentification/Filtered_taxa_withPlate_chemTaxa.csv")
taxa_group <- read.csv(file = "group_assignments_df.csv")
lin_biolog_mag <- read.csv("../4_amplicaonMAGsLinking/16s_mags_linked.csv")

taxa_table$sub_class <- ifelse(is.na(taxa_table$sub_class), "Other", taxa_table$sub_class)

uq_taxa <- taxa_table %>% select(asv, taxa) %>% unique() %>% mutate(unique_taxa = make.unique(taxa)) %>% select(-taxa)

taxa_fil<-taxa_table %>% 
  filter(asv %in% lin_biolog_mag$asv) %>%
  group_by(metabolite, asv, taxa) %>% 
  summarise(avg_ra = mean(ra), .groups = "drop") %>%  
  left_join(.,uq_taxa, by="asv") %>% 
  left_join(., taxa_group, by = "taxa")

write.csv(x = taxa_fil, file = "Filter_taxa_link16s_cluster.csv")
