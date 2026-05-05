library(tidyverse)

df_assem <- read_tsv("../2_reconstructionMAGs/genome_summary.tsv")

df_16slinked <- read_csv("16s_mags_linked.csv")

df_genomequality <- read_tsv("../2_reconstructionMAGs/genome_properties.tsv")

# filter for only MAGs
df_assem_MAG <- df_assem %>% filter(grepl("MG_|MT_", `Bin Id`))

median(df_assem_MAG $Completeness)
# 96.96 %
median(df_assem_MAG $Contamination)
# 0.65 %

p_MAGscomcon <- df_assem_MAG %>% select(`Bin Id`, Completeness, Contamination) %>% 
  pivot_longer(cols = -`Bin Id`, names_to = "qc", values_to = "percent") %>% 
  ggplot(aes(x = percent)) +
  geom_histogram(fill = "#4CC9F0", colour = "black")+
  facet_wrap(~ qc, scales = "free")+
  theme_bw()+
  labs(y = "Frequency")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90))
  
#ggsave(filename = "p_MAGsComplandConta.png",
#       plot = p_MAGscomcon,
#       width = 3, 
#       height = 3, 
#       units = "in", 
#       dpi = 300)

p_MAGsquality <- df_genomequality %>% 
  filter(genome %in% df_assem$`Bin Id`) %>% 
  select(genome, num_seq, sum_len) %>% 
  pivot_longer(cols = -genome, names_to = "qc", values_to = "bps") %>% 
  ggplot(aes(x = bps))+
  geom_histogram(fill = "#4CC9F0", colour = "black")+
  facet_wrap(~ qc, scales = "free")+
  theme_bw()+
  labs(y = "Frequency")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90))

#ggsave(filename = "p_MAGsquality.png",
#       plot = p_MAGsquality,
#       width = 3, 
#       height = 3, 
#       units = "in", 
#       dpi = 300)



p_selectedGenomeComCon <- df_assem %>% 
  filter(`Bin Id` %in% df_16slinked$GenomicSeq) %>% 
  select(`Bin Id`, Completeness, Contamination) %>% 
  pivot_longer(cols = -`Bin Id`, names_to = "qc", values_to = "percent") %>% 
  ggplot(aes(x = percent)) +
  geom_histogram(fill = "#4CC9F0", colour = "black")+
  facet_wrap(~ qc, scales = "free")+
  theme_bw()+
  labs(y = "Frequency")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90))

#ggsave(filename = "p_selectedGenomeComCon.png",
#         plot = p_selectedGenomeComCon,
#         width = 3, 
#         height = 3, 
#         units = "in", 
#         dpi = 300)

p_selectedGenomeQuality <- df_genomequality %>% 
  filter(genome %in% df_16slinked$GenomicSeq) %>% 
  select(genome, num_seq, sum_len) %>% 
  pivot_longer(cols = -genome, names_to = "qc", values_to = "bps") %>% 
  ggplot(aes(x = bps))+
  geom_histogram(fill = "#4CC9F0", colour = "black")+
  facet_wrap(~ qc, scales = "free")+
  theme_bw()+
  labs(y = "Frequency")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90))
  

#ggsave(filename = "p_selectedGenomeQuality.png",
#       plot = p_selectedGenomeQuality ,
#       width = 3, 
#       height = 3, 
#       units = "in", 
#       dpi = 300)

