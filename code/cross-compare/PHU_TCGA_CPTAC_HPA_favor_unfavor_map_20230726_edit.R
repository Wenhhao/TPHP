rm(list = ls())
setwd('Y:/members/jiangwenhao/TPHP/20220908/HPA_map')
library(tidyverse)
library(magrittr)

read_excel_allsheets <- function(filename, tibble = FALSE, col_types = NULL) { # such as col_types = 'numeric'
  # I prefer straight data.frames
  # but if you like tidyverse tibbles (the default with read_excel)
  # then just pass tibble = TRUE
  sheets <- readxl::excel_sheets(filename)
  if(is.null(col_types)){
    x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  } else{
    x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X, col_types = col_types))
  }
  
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}




# 1. TCGA, CPTAC ----------------------------------------------------------
res <- readxl::read_xlsx("//172.16.13.136/share/members/jiangwenhao/TPHP/TCGA/PUH_TCGA_CPTAC_DEA_nofilterP_chr_20230626.xlsx")

df_TCGA <- res %>%
  filter(!is.na(organ)) %>% 
  filter((abs(PUH_log2FC) > log2(1.6)) & (abs(TCGA_log2FC) > log2(1.6))) %>% 
  filter(abs(TCGA_consistency) == 1) #%>% # filter
  # select(organ, PUH_log2FC, TCGA_log2FC, TCGA_consistency)

df_CPTAC <- res %>% 
  filter(!is.na(organ)) %>% 
  filter((abs(PUH_log2FC) > log2(1.6)) & (abs(CPTAC_log2FC) > log2(1.2))) %>%
  filter(abs(CPTAC_consistency) == 1) #%>% # filter
  # select(organ, PUH_log2FC, CPTAC_log2FC, CPTAC_consistency)

# organ name to cancer name
ht_organ2cancer <- c('Glioblastoma', 'Cervical carcinoma', 'Colon carcinoma', 'Esophageal carcinoma', 'Fallopian tube carcinoma', 'Gallbladder carcinoma', 'Renal carcinoma', 'Hepatocellular carcinoma', 'Lung carcinoma', 'Diffused large B-cell carcinoma', 'Breast carcinoma', 'Muscle tumor', 'Ovarian carcinoma', 'Pancreas carcinoma', 'Male genitalia carcinoma', 'Prostate carcinoma', 'Rectum carcinoma', 'Gastrointestinal stromal tumors', 'Gastric carcinoma', 'Testis carcinoma', 'Laryngocarcinoma', 'Thymoma and thymic carcinoma', 'Thyroid carcinoma', 'Laryngocarcinoma', 'Endometrial carcinoma')
names(ht_organ2cancer) <- c('brain', 'cervix_uteri', 'colon', 'esophagus', 'fallopian tube', 'gall bladder', 'kidney', 'liver', 'lung', 'lymph node', 'mammary gland', 'muscle', 'ovary', 'pancreas', 'penis', 'prostate', 'rectum', 'small intestine', 'stomach', 'testis', 'throat', 'thymus', 'thyroid', 'tongue', 'uterus')
df_TCGA$cancer <- ht_organ2cancer[df_TCGA$organ]
df_CPTAC$cancer <- ht_organ2cancer[df_CPTAC$organ]

df_TCGA %<>% add_column(Mapped = 'TCGA', .before = 1)
df_CPTAC %<>% add_column(Mapped = 'CPTAC', .before = 1)

df_comb <- rbind(df_TCGA, df_CPTAC)
df_comb %>% rio::export('comb_TCGA_CPTAC_filtered.xlsx')








# 2. HPA ------------------------------------------------------------------
df_dea <- rio::import('//172.16.13.136/share/members/jiangwenhao/TPHP/20220908/diff_expr/20230318_TPHP_dysregulated_proteins_filter50NAByOrgan.xlsx')
df_dea[str_detect(df_dea$organ, 'cervix'), 'organ'] <- 'cervix uteri'
df_dep <- df_dea %>% filter(abs(log2FC) >= log2(1.5), pAdj_t < 0.05)
df_dep %<>% rename(Uniprot = protein)
df_dep$dysregulation <- ifelse(df_dep$log2FC > 0, 'Up', 'Down')
df_dep$organ %>% unique
# [1] "brain"           "cervix uteri"    "colon"           "gall bladder"    "liver"           "lung"        
# [7] "lymph node"      "mammary gland"   "muscle"          "ovary"           "pancreas"        "penis"       
# [13] "rectum"          "small intestine" "stomach"         "testis"          "throat"          "thyroid"    
# [19] "tongue"          "uterus"  

# map PUH organ names to HPA cancer names
ht_organ2cancer <- c('', 'cervical', 'colorectal', '', 'liver', 'lung', '', 'breast', '', 'ovarian', 'pancreatic', '', 'colorectal', '', 'stomach', 'testis', 'head_and_neck', 'thyroid', '', 'endometrial')
names(ht_organ2cancer) <- c('brain', 'cervix uteri', 'colon', 'gall bladder', 'liver', 'lung', 'lymph node', 'mammary gland', 'muscle', 'ovary', 'pancreas', 'penis', 'rectum', 'small intestine', 'stomach', 'testis', 'throat', 'thyroid', 'tongue', 'uterus')

df_dep$cancer <- ht_organ2cancer[df_dep$organ]
df_dep$cancer %<>% str_c(., ' cancer') %>% str_replace_all('_', '') %>% str_to_sentence()


# read HPA supplementary table 6
hpa_tbls <- read_excel_allsheets('Y:/members/jiangwenhao/TPHP/20220908/HPA_map/Table S6.xlsx')
df_hpa <- plyr::ldply(hpa_tbls, .id = 'cancer')
df_hpa %<>% filter(cancer %in% unique(df_dep$cancer)) # 11cancers




# read uniprot data
# df_uniprot <- rio::import('Y:/members/jiangwenhao/TPHP/20220908/drug/uniprot/uniprot-download_true_fields_accession_2Cid_2Cprotein_name_2Cgene_na-2022.09.21-14.18.23.22.xlsx')
df_uniprot <- rio::import('Y:/members/jiangwenhao/TPHP/20220908/HPA_map/uniprot-download_true_fields_accession_2Creviewed_2Cid_2Cgene_names_-2023.06.16-09.12.08.30.tsv')

# combine
df <- df_uniprot %>% select(From, Entry) %>% rename(EnsemblIDs = From, Uniprot = Entry) %>% 
  left_join(df_hpa, .) %>% 
  inner_join(df_dep, .)


# 
df$MATCHED <- (df$`Prognostic Types` == 'Unfavorable' & df$dysregulation == 'Up') | (df$`Prognostic Types` == 'Favorable' & df$dysregulation == 'Down')

# &
df %<>% select(-p_t, -p_wilcox, -pAdj_wilcox) # 
df %<>% select(cancer, organ, Uniprot, dysregulation, `Prognostic Types`, MATCHED, everything())

# log-rank p value
df %<>% filter(`log-rank P Values` < 0.05)


df %>% count(cancer)
# cancer    n
# 1       Breast cancer 2289
# 2     Cervical cancer 1882
# 3   Colorectal cancer 1254
# 4  Endometrial cancer 2813
# 5        Liver cancer  533
# 6         Lung cancer  379
# 7      Ovarian cancer 1778
# 8   Pancreatic cancer  216
# 9      Stomach cancer  169
# 10      Testis cancer  906
# 11     Thyroid cancer  645

df %>% count(dysregulation)
# dysregulation     n
# 1          Down   925
# 2            Up 11939

df %>% count(`Prognostic Types`)
# Prognostic Types    n
# 1        Favorable 6881
# 2      Unfavorable 5983

df_dea <- rio::import('//172.16.13.136/share/members/jiangwenhao/TPHP/20220908/diff_expr/20230318_TPHP_dysregulated_proteins_filter50NAByOrgan.xlsx')
df_dea[str_detect(df_dea$organ, 'cervix'), 'organ'] <- 'cervix uteri'
df_dep <- df_dea %>% filter(abs(log2FC) >= log2(2), pAdj_t < 0.05)
df_dep %<>% rename(Uniprot = protein)
df_dep$dysregulation <- ifelse(df_dep$log2FC > 0, 'Up', 'Down')
df_dep$organ %>% unique
# [1] "brain"           "cervix uteri"    "colon"           "gall bladder"    "liver"           "lung"        
# [7] "lymph node"      "mammary gland"   "muscle"          "ovary"           "pancreas"        "penis"       
# [13] "rectum"          "small intestine" "stomach"         "testis"          "throat"          "thyroid"    
# [19] "tongue"          "uterus"  

# map PUH organ names to HPA cancer names
ht_organ2cancer <- c('', 'cervical', 'colorectal', '', 'liver', 'lung', '', 'breast', '', 'ovarian', 'pancreatic', '', 'colorectal', '', 'stomach', 'testis', 'head_and_neck', 'thyroid', '', 'endometrial')
names(ht_organ2cancer) <- c('brain', 'cervix uteri', 'colon', 'gall bladder', 'liver', 'lung', 'lymph node', 'mammary gland', 'muscle', 'ovary', 'pancreas', 'penis', 'rectum', 'small intestine', 'stomach', 'testis', 'throat', 'thyroid', 'tongue', 'uterus')

df_dep$cancer <- ht_organ2cancer[df_dep$organ]
df_dep$cancer %<>% str_c(., ' cancer') %>% str_replace_all('_', '') %>% str_to_sentence()


# read HPA supplementary table 6
hpa_tbls <- read_excel_allsheets('Y:/members/jiangwenhao/TPHP/20220908/HPA_map/Table S6.xlsx')
df_hpa <- plyr::ldply(hpa_tbls, .id = 'cancer')
df_hpa %<>% filter(cancer %in% unique(df_dep$cancer)) # 11cancers




# read uniprot data
# df_uniprot <- rio::import('Y:/members/jiangwenhao/TPHP/20220908/drug/uniprot/uniprot-download_true_fields_accession_2Cid_2Cprotein_name_2Cgene_na-2022.09.21-14.18.23.22.xlsx')
df_uniprot <- rio::import('uniprot-download_true_fields_accession_2Creviewed_2Cid_2Cgene_names_-2023.06.16-09.12.08.30.tsv')

# combine
df <- df_uniprot %>% select(From, Entry) %>% rename(EnsemblIDs = From, Uniprot = Entry) %>% 
  left_join(df_hpa, .) %>% 
  inner_join(df_dep, .)


# 
df$MATCHED <- (df$`Prognostic Types` == 'Unfavorable' & df$dysregulation == 'Up') | (df$`Prognostic Types` == 'Favorable' & df$dysregulation == 'Down')

# &
df %<>% select(-p_t, -p_wilcox, -pAdj_wilcox) # 
df %<>% select(cancer, organ, Uniprot, dysregulation, `Prognostic Types`, MATCHED, everything())

# log-rank p value
df %<>% filter(`log-rank P Values` < 0.05)
rio::export(df, 'PHU_dysregulation_HPA_favor_unfavor_overlap_20230616.xlsx')





df %>% count(cancer)
# cancer    n
# 1       Breast cancer 2289
# 2     Cervical cancer 1882
# 3   Colorectal cancer 1254
# 4  Endometrial cancer 2813
# 5        Liver cancer  533
# 6         Lung cancer  379
# 7      Ovarian cancer 1778
# 8   Pancreatic cancer  216
# 9      Stomach cancer  169
# 10      Testis cancer  906
# 11     Thyroid cancer  645

df %>% count(dysregulation)
# dysregulation     n
# 1          Down   925
# 2            Up 11939

df %>% count(`Prognostic Types`)
# Prognostic Types    n
# 1        Favorable 6881
# 2      Unfavorable 5983


X <- df_dep %>% count(cancer, name = 'x')
Y <- df %>% count(cancer, name = 'y')
Z <- X %>% inner_join(Y, by = 'cancer') %>% mutate(ratio = round(100 * y / x, 2))
rio::export(Z, 'PUH_dep_HPA_FavorUnfavor_ratio_20230628.xlsx')




# # sankey plot
# colnames(df)
# # [1] "cancer"             "organ"              "Uniprot"            "dysregulation"     
# # [5] "Prognostic Types"   "MATCHED"            "log2FC"             "pAdj_t"            
# # [9] "EnsemblIDs"         "Symbols"            "Mean Expression"    "Sample Numbers"    
# # [13] "Expression Cutoffs" "log-rank P Values" 
library(ggalluvial)
# Sheet1 <- df %>% count(cancer)
Sheet1 <- Z %>% select(cancer, y, ratio) %>% mutate(n = str_glue("{y}, {ratio}%")) %>% select(cancer, n)
Sheet2 <- df %>% count(dysregulation)
Sheet3 <- df %>% count(dysregulation, `Prognostic Types`)
Sheet4 <- df %>% count(`Prognostic Types`)
df_n <- list(Sheet1, Sheet2, Sheet4) %>% plyr::ldply(function(dfsub){
  colnames(dfsub) <- c('from', 'to')
  return(dfsub)
})
df_n %<>% mutate(to = str_glue("{from} (n = {to})"))
ht_fromto <- df_n$to
names(ht_fromto) <- df_n$from

df_alluvium <- df %>% count(cancer, dysregulation, `Prognostic Types`)
df_alluvium %<>% mutate_at(vars(-n), function(x) ht_fromto[x])

# set.seed(2023)
# my_colors <- sample(colorRampPalette(c("#2266BB", "#66BB22", "#BB2266"),bias=1)(n=12))
my_colors <- c('#2B2B6B', '#9D7DDB', '#F6D151', '#D8894E', '#E58989', '#BA6CA4', '#D34A8F', '#103D47', '#8080BC', '#55BC6D', '#7C0823', '#2578A0')
p <- ggplot(df_alluvium,
            aes(y = n, axis1 = cancer, axis2 = dysregulation, axis3 = `Prognostic Types`)) +
  geom_alluvium(aes(fill = cancer), alpha=1,
                width = 4/8, knot.pos = 0, reverse = FALSE) +
  scale_fill_manual(values = my_colors) +
  guides(fill = "none") +
  geom_stratum(alpha = .5, width = 4/8, reverse = FALSE) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            reverse = FALSE) +
  scale_x_continuous(breaks = 1:3, labels = c("Caner name", "Dysrequlation (PUH dataset)", "Prognostic Types (HPA dataset)")) +
  # coord_flip() +
  ggtitle("PUH dysregulated-HPA favorable/unfavorable proteins") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        # axis.line = element_line(colour = 'black'),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title.x = element_blank(),
        # axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
  )


ggsave('Y:/members/jiangwenhao/TPHP/20220908/HPA_map/PUH_dep_HPA_FavorUnfavor_sankey_20230628.pdf', p, height = 10, width = 7)



# 3. combine --------------------------------------------------------------
dfprot <- rio::import('//172.16.13.136/TPHP/TPL/libs/20220616_fraglib/protein.tsv')
df_comb %<>% left_join(dfprot, by = c('protein' = 'Protein ID'))

df %<>% rename(HPA_cancer = cancer)
df %<>% add_column(protein = df$Uniprot, .before = 1)

intersect(colnames(df), colnames(df_comb)) # organ
unique(df$organ)
unique(df_comb$organ)

df_final_outer <- df_comb %>% full_join(df, by = c('organ', 'protein')) %>% select(cancer, organ, protein, everything())
df_final <- df_comb %>% inner_join(df, by = c('organ', 'protein')) %>% select(cancer, organ, protein, everything())
df_final1 <- df_final %>% filter(MATCHED)

list(
  `PUH-TCGA&CPTAC-HPA` = df_final1 %>% filter(TCGA_significant, CPTAC_significant) %>% arrange(cancer, protein),
  `PUH-TCGA|CPTAC-HPA` = df_final1 %>% arrange(cancer, protein),
  `PUH-TCGA-CPTAC-HPA_concat` = df_final_outer %>% arrange(cancer, protein)
) %>% rio::export('PUH_TCGA_CPTAC_HPA_mapping_20230726.xlsx')



df_final1 %>% filter(CPTAC_significant)  %>% select(cancer, organ, protein, everything())
