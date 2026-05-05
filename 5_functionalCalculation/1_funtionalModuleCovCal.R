# Library load ----
library("tidyverse")
library("stringr")

# 1. Prepare the database file and also get some statistic for more understanding. ----

# All databases below were downloaded on 08 January 2025
link_ko_module <- read.csv("KEGG_moduleDB/08012025_KEGGDataBase_Linked_KO_Module.csv")
list_ko <- read.csv("KEGG_moduleDB/08012025_KEGGDataBase_List_KO.csv")
list_module <- read.csv("KEGG_moduleDB/08012025_KEGGDataBase_List_Module.csv")

link_ko_module <- link_ko_module %>% separate(KOtoModule, into = c("ko", "module"), sep = "\\s+") %>%
  mutate(ko = str_remove(ko, "ko:"), module = str_remove(module, "md:"))
list_ko <- list_ko %>% separate(KODes, into = c("ko", "ko_description"), sep = "  ")
list_module <- list_module %>% separate(ModuleDes, into = c("module", "module_description"), sep = "  ")

comb_ko_module <- link_ko_module %>% 
  left_join(., list_ko, by="ko") %>% 
  left_join(., list_module, by="module")

print(paste("Number of KO:", as.character(length(unique(list_ko$ko))),sep=" "))
print(paste("Number of Module:", as.character(length(unique(list_module$module))),sep=" "))
print(paste("Number of KO (in KEGG REST link):", as.character(length(unique(link_ko_module$ko))),sep=" "))
print(paste("Number of Module (in KEGG REST link):", as.character(length(unique(link_ko_module$module))),sep=" "))

# In total KEGG database:
# - Number of KO: 27321
# - Number of Module: 496
# However in KEGG REST link, 
# - Number of KO: 2862
# - Number of module: 490

# This suggest that some KOs are not in the module (maybe not relate to the metabolism), also some modules do not have KO.
# Further analysis would exclude KOs and modules that are not in the KEGG REST link.

list_ko <- list_ko %>% filter(ko %in% link_ko_module$ko)
list_module <- list_module %>% filter(module %in% link_ko_module$module)

write.csv(comb_ko_module, "KEGGDataBase_Linked_KO_Module_reformat.csv")
write.csv(list_ko, "KEGGDataBase_List_KO_reformat.csv")
write.csv(list_module, "KEGGDataBase_List_Module_reformat.csv")

p_Hist_noKOinModule <- comb_ko_module %>% 
  group_by(module) %>% 
  summarise(count = n()) %>% 
  ggplot(aes(count))+
  geom_histogram()+
  theme_classic()+
  labs(x="Number of KO in a KEGG module",
       y="Count")

#ggsave("Hist_noKOinModule.png", 
#       width = 4, 
#       height = 4)

# 2. Count KO profile for each genome and calculate the KEGG module coverage and abundance. ----

# Extract KEGG orthology term from each bin (eggNOG results)
# Note that all necessary files are moved to new local directory.
# Use the directory for the below script

bin_file<-read_lines("/filename.txt")

ko_number_df<-data.frame("Bin_name"=character(),
                         "All_query"=integer(),
                         "Annotated_query"=integer(),
                         "Nonannotated_query"=integer(),
                         "All_KO"=integer())
allMAGs_NormMoAbun<-list_module
allMAGs_MoCov<-list_module
allMAGs_MoAbun<-list_module
for (f in bin_file) {
  bin_name<-sub("\\.emapper\\.annotations$", "",f)
  file.path<-paste("/annotation_file/", f, sep = "")
  file.anno<-readLines(file.path)
  
  n<-length(file.anno)
  
  file.anno.keep<-file.anno[5:(n-3)]
  file.anno.tsv<-read_tsv(file = I(file.anno.keep))
  file.path.csv<-paste("/annotation_file/csv_file/", f, ".csv", sep = "")
  
  write.csv(file.anno.tsv, file.path.csv)
  
  ko_non<-length(which(file.anno.tsv$KEGG_ko=="-"))
  ko_all<-length(file.anno.tsv$KEGG_ko)
  ko_anno<-ko_all-ko_non
  
  ko_duplicate_array<-file.anno.tsv %>% 
    select(KEGG_ko) %>% 
    filter(KEGG_ko!="-") %>% 
    separate_rows(KEGG_ko,sep = ",") %>% 
    mutate(KEGG_ko = gsub("^ko:", "", KEGG_ko))
  
  ko_anno_all<-length(ko_duplicate_array$KEGG_ko)
  
  ko_number_df<-rbind(ko_number_df, data.frame("Bin_name"=bin_name,
                                               "All_query"=ko_all,
                                               "Annotated_query"=ko_anno,
                                               "Nonannotated_query"=ko_non,
                                               "All_KO"=ko_anno_all))
  
  mo_abundance_df<-data.frame("KEGG_module"=character(),
                              "Description"=character(),
                              "Module_size"=numeric(),
                              "KO_unique"=numeric(),
                              "KO_abundance"=numeric(),
                              "KO_total"=numeric(),
                              "Module_coverage"=numeric(),
                              "Module_abundance"=numeric(),
                              "Normalized_module_abundance"=numeric())
  
  for (i in 1:nrow(list_module)) {
    ko_abun<-0
    ko_unique<-0
    koInMo<-comb_ko_module %>% 
      filter(module == list_module$module[i]) %>% 
      pull(ko)
    n_koPerMo<-length(koInMo)
    
    for (j in 1:n_koPerMo) {
      ko_j<-koInMo[j]
      count_ko<-length(which(ko_duplicate_array==ko_j))
      count_ko_unique<-length(which(unique(ko_duplicate_array)==ko_j))
      ko_abun<-ko_abun+count_ko
      ko_unique<-ko_unique+count_ko_unique
    }
    if (ko_abun!=0) {
      mo_cov<-ko_unique/n_koPerMo
      mo_abun<-ko_abun*mo_cov*1000
      mo_abun_norm<-mo_abun/ko_anno_all
      mo_abundance_df<-rbind(mo_abundance_df, data.frame("KEGG_module"=list_module$module[i],
                                                         "Description"=list_module$module_description[i],
                                                         "Module_size"=n_koPerMo,
                                                         "KO_unique"=ko_unique,
                                                         "KO_abundance"=ko_abun,
                                                         "KO_total"=ko_anno_all,
                                                         "Module_coverage"=mo_cov,
                                                         "Module_abundance"=mo_abun,
                                                         "Normalized_module_abundance"=mo_abun_norm))
    } else {
      mo_abundance_df<-rbind(mo_abundance_df, data.frame("KEGG_module"=list_module$module[i],
                                                         "Description"=list_module$module_description[i],
                                                         "Module_size"=n_koPerMo,
                                                         "KO_unique"=ko_abun,
                                                         "KO_abundance"=ko_abun,
                                                         "KO_total"=ko_anno_all,
                                                         "Module_coverage"=ko_abun,
                                                         "Module_abundance"=ko_abun,
                                                         "Normalized_module_abundance"=ko_abun))
    }
  }
  
  allMAGs_NormMoAbun<-cbind(allMAGs_NormMoAbun, mo_abundance_df$Normalized_module_abundance) %>% 
    rename(!!bin_name := 'mo_abundance_df$Normalized_module_abundance')
  allMAGs_MoCov<-cbind(allMAGs_MoCov, mo_abundance_df$Module_coverage) %>% 
    rename(!!bin_name := 'mo_abundance_df$Module_coverage')
  allMAGs_MoAbun<-cbind(allMAGs_MoAbun, mo_abundance_df$Module_abundance) %>% 
    rename(!!bin_name := 'mo_abundance_df$Module_abundance')
}

write.csv(ko_number_df,"ko_stat_all.csv")
write.csv(allMAGs_NormMoAbun, "allMAGs_NormMoAbun.csv")
write.csv(allMAGs_MoCov, "allMAGs_MoCov.csv")
write.csv(allMAGs_MoAbun, "allMAGs_MoAbun.csv")