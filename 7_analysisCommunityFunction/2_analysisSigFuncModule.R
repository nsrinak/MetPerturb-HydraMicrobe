library("tidyverse")
library("ggplot2")
library("vegan")
library("RColorBrewer")

# Data ----
lm_metabolite_wFDR <- read.csv("../5_functionalCalculation/LM_sepPlate_metabolite_withT_FDR_newlmer.csv")

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
  unique() %>% 
  mutate(sample2 = sample) %>% 
  column_to_rownames(var = "sample2")

lin_biolog_mag <- read.csv("../4_amplicaonMAGsLinking/16s_mags_linked.csv")

# Per sample - number of sig module, Shannon diversity, and microbial relative abundance ----

ord <- lm_metabolite_wFDR %>% 
  filter(FDR< 0.05) %>% 
  group_by(metabolite) %>% 
  summarise(metabolite_n = n()) %>%
  arrange(metabolite_n) %>% 
  ungroup() %>% pull(metabolite)

uq_taxa <- taxa_table %>% 
  dplyr::select(asv, taxa) %>% 
  unique() %>% 
  mutate(unique_taxa = make.unique(taxa)) %>% 
  dplyr::select(-taxa)

taxa_fil<-taxa_table %>% 
  filter(asv %in% lin_biolog_mag$asv) %>%
  group_by(metabolite, asv, taxa) %>% 
  summarise(avg_ra = mean(ra)) %>% ungroup() %>% 
  left_join(.,uq_taxa, by="asv") %>%
  filter(metabolite %in% ord)

taxa_fil$metabolite <- factor(taxa_fil$metabolite, levels = ord)

cute_palette <- colorRampPalette(brewer.pal(12, "Paired"))(35)

p_bar_taxa_kosigmet <- ggplot(taxa_fil, aes(x=metabolite, y = avg_ra, fill = unique_taxa))+
  geom_bar(position = "stack", stat = "identity", color = "black")+
  scale_fill_manual(values = cute_palette)+
  theme_bw()+
  labs(y = "Relative abundance", x = "Metabolite", fill = "Taxa")+
  theme(
    legend.position = "none",
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.3, "cm"),  
    axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0.5, size = 7),
    axis.text.y = element_text(size = 7),
    axis.title.y = element_text(size = 8), 
    axis.title.x = element_text(size = 8)
  )

p_bar_msig <- lm_metabolite_wFDR %>% 
  filter(FDR<= 0.05) %>% 
  group_by(metabolite) %>%
  mutate(metabolite_n = n()) %>%
  mutate(direction = ifelse(Estimate > 0, "up", ifelse(Estimate < 0, "down", "neutral"))) %>%
  arrange(desc(metabolite_n)) %>% 
  ungroup() %>% 
  ggplot(aes(x=reorder(metabolite, metabolite_n, FUN = min), fill = direction))+
  geom_bar(stat = "count", width = 0.5, color = "black")+
  scale_fill_manual(values=c("up" = "darkgreen",
                             "down" = "brown3"))+
  theme_bw()+
  labs(y = "Number of\nsignificant module", x = NULL, fill = "Direction")+
  theme(axis.text.x = element_blank(),  
        legend.position = "none", 
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 8),
        plot.margin = margin(t = 6, unit = "pt"))

shannon <- taxa_table %>% group_by(sample) %>% 
  mutate(shannon = diversity(count, index = "shannon")) %>% 
  filter(metabolite %in% ord)

shannon$metabolite<-factor(shannon$metabolite, levels = ord)


p_shannon_msig <-ggplot(shannon, aes(x=metabolite, y=shannon))+
  geom_boxplot(fill = "#4CC9F0", color = "black")+
  theme_bw()+
  labs(y = "Shannon index", x = NULL)+
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 8))


p_comb_taxa_nkosig_shannon <- plot_grid(p_bar_msig,
                                        p_shannon_msig, 
                                        p_bar_taxa_kosigmet, 
                                        align = "v", ncol = 1, rel_heights = c(13.5,13,70))

#ggsave(filename = "p_comb_taxa_nkosig_shannon_movetick.png",
#       plot = p_comb_taxa_nkosig_shannon,
#       width = 4.5, 
#       height = 6.7, 
#       units = "in", 
#       dpi = 300)