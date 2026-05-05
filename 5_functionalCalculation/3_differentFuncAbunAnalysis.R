library(tidyverse)
library(lmerTest)
library(lme4)


# Linear mixed model between control and non-control ----

taxa_table <- read.csv("../1_ampliconIdentification/Filtered_taxa_withPlate_chemTaxa.csv")

abun <- read.csv("Community_module_abundance.csv", check.names = FALSE)

meta_data <- taxa_table %>% 
  select(sample, super_class, class, sub_class, direct_parent, metabolite, plate.x) %>% 
  unique() 

abun <- abun %>%  
  pivot_longer(-c(module, module_description),names_to = "sample", values_to = "module_abun") %>% 
  left_join(., meta_data, by = "sample") %>% 
  group_by(sample) %>% 
  mutate(rank = rank(-module_abun)) %>% 
  ungroup()

unique_module <- abun %>% pull(module) %>% unique()
unique_plate <- abun %>% pull(plate.x) %>% unique()
unique_metabolite <- abun %>% filter(metabolite != "Control") %>%  pull(metabolite) %>% unique()
lm_sepPlate_metabolite <- data.frame("metabolite" = character(),
                                     "module" = character(),
                                     "Estimate" = numeric(),
                                     "Std.Error" = numeric(),
                                     "tvalue" = numeric(),
                                     "pvalue" = numeric(),
                                     "var_rand" = numeric())
met <- 0


for (i in unique_metabolite){
  print(met)
  met <- met + 1
  plate_met <- abun %>% filter(metabolite == i) %>% pull(plate.x) %>% unique()
  abun_1 <- abun %>% filter(metabolite %in% c(i, "Control") & plate.x == plate_met)
  for (j in unique_module){
    df <- abun_1 %>% filter(module==j) %>% mutate(metabolite = as.factor(metabolite))
    df[["metabolite"]] <- relevel(df[["metabolite"]], "Control")
    if (var(df$module_abun)==0){
      df_1 <- data.frame("metabolite" = i,
                         "module" = j,
                         "Estimate" = -100,
                         "Std.Error" = -100,
                         "tvalue" = -100,
                         "pvalue" = -100,
                         "var_rand" = -100)
      lm_sepPlate_metabolite <- rbind(lm_sepPlate_metabolite,df_1)  
    } else {
      if (length(plate_met) == 1){
        lm_rs <- lm(module_abun ~ metabolite, data = df)
        variance_rand <- -9
      } else {
        lm_rs <- lmerTest::lmer(module_abun ~ metabolite + (1 | plate.x), data = df)
        random_effects <- VarCorr(lm_rs)
        std_dev_rand <- random_effects$plate.x[1]
        variance_rand <- std_dev_rand^2
      }
      sum_mod.lm <- summary(lm_rs)
      df_1 <- data.frame("metabolite" = i,
                         "module" = j,
                         "Estimate" = sum_mod.lm$coefficients[, "Estimate"][[2]],
                         "Std.Error" = sum_mod.lm$coefficients[, "Std. Error"][[2]],
                         "tvalue" = sum_mod.lm$coefficients[, "t value"][[2]],
                         "pvalue" = sum_mod.lm$coefficients[, "Pr(>|t|)"][[2]],
                         "var_rand" = variance_rand)
      lm_sepPlate_metabolite <- rbind(lm_sepPlate_metabolite,df_1)
    }
  }
}
view(lm_sepPlate_metabolite)

# From 490 modules, there are 312 modules that at appear in least 1 genome.
# Here, I am using 312 modules as a total size of the module in this system. This way, the result could be more explainable, but also more specific.

allMAGs_MoCov <- read.csv("allMAGs_MoCov.csv", check.names = FALSE)

lin_biolog_mag <- read.csv("../3_amplicaonMAGsLinking/16s_mags_linked.csv")

asv <- taxa_table %>% 
  select(asv,taxa) %>%  
  unique()

lin_biolog_mag <- lin_biolog_mag %>% 
  left_join(., asv, by="asv") %>% 
  rename(genome=GenomicSeq) %>% 
  mutate(unique_taxa = make.unique(taxa))

df_cov <- allMAGs_MoCov%>% 
  select(c(lin_biolog_mag$genome, module, module_description)) %>%
  pivot_longer(-c(module, module_description), names_to = "genome", values_to = "coverage") %>% 
  left_join(., lin_biolog_mag, by = "genome") %>% glimpse()

non_zero_module <- df_cov %>% group_by(module) %>% 
  summarise(sum_cov = sum(coverage)) %>% 
  filter(sum_cov > 0) %>% 
  pull(module)

lm_metabolite_wFDR <- data.frame("metabolite"=character(),
                                 "KEGG_module"=character(),
                                 "coefficient"=numeric(),
                                 "pvalue"=numeric(),
                                 "plate"=character(),
                                 "FDR"=numeric(),
                                 "t"=numeric(),
                                 "SE"=numeric(),
                                 "CIL"=numeric(),
                                 "CIH"=numeric(),
                                 "df_error"=numeric())
n<-0

for (i in unique_metabolite){
  print(n)
  n<- n+1
  test <- lm_sepPlate_metabolite %>% filter(metabolite==i & var_rand != -100 & module %in% non_zero_module)
  bh <- p.adjust(as.vector(test$pvalue), method = "BH", n=length(non_zero_module))
  bh <- data.frame("FDR"=bh)
  c <- cbind(test, bh)
  lm_metabolite_wFDR <- rbind(lm_metabolite_wFDR, c)
}

lm_metabolite_wFDR %>% filter(FDR < 0.05) %>% view()
write.csv(lm_metabolite_wFDR, "LM_sepPlate_metabolite_withT_FDR_newlmer.csv")
