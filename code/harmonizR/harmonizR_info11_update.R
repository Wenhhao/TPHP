rm(list = ls())
pacman::p_unload(pacman::p_loaded(), character.only = T)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(tidyverse)
library(magrittr)
library(RColorBrewer)
library(scales)
library(pheatmap)
library(Rtsne)
library(umap)
library(ggpubr)
library(ggsci)
library(job)
source('../source/my_fun.R')




# 0.Data preparation ---------
## dat1 --------
dat1 <- read.delim('mapped_pg_matrix_1780_14062.csv', sep = ',', header = T, row.names = 1, check.names = F, stringsAsFactors = F)
info1 <- rio::import('20250725_PUH_sample_information_1900files_v10_addLabels.xlsx')
# meta1 <- info1 %>% select(FileName, file_id, sample_id, sample_type, Inner_ring_label, tissue_name, date)
pool1 <- read.delim('pool_pg_matrix_120_14063.csv', sep = ',', header = T, row.names = 1, check.names = F, stringsAsFactors = F)

dfrm1 <- rio::import('tims_QC_source_checked.xlsx', sheet = 'sample.identity')
rm1 <- dfrm1$FileName[dfrm1$Is.Lower.Ingroup]
rm1.pool <- 'N20210825yuel_nail_pool_Slot1-6_1_7467'

info1 %<>% filter(!(FileName %in% c(rm1, rm1.pool)))
dat1 <- dat1[intersect(rownames(dat1), info1$FileName), ]
pool1 <- pool1[setdiff(rownames(pool1), rm1.pool), ]

pm1 <- rbind(dat1, pool1) %>% t()
dim(pm1) # 14062  1811
pm1 <- pm1[rowSums(!is.na(pm1)) != 0, colSums(!is.na(pm1)) != 0] %>% .[, info1$FileName]
dim(pm1) # 14058  1811



# 
# # # NA imputation; min * 0.8 * e~N(1, 0.01)
# # # pm1[is.finite(pm1)] <- NA
# # set.seed(2023)
# # na_pool <- min(pm1, na.rm = T) + log2(0.8) + log2(rnorm(sum(is.na(pm1)), 1, 0.01))
# # pm1[is.na(pm1)] <- na_pool
# # dim(pm1) # 7181  975



## dat2 --------
info_merge <- rio::import('20250809_PUH_RCA_sample_and_pool_1179files_info.xlsx')
dat2 <- read.delim('mapped_dat2_pg_matrix_1152_12500.csv', sep = ',', header = T, row.names = 1, check.names = F, stringsAsFactors = F)
pool2 <- read.delim('mapped_dat2_pool_matrix_27_12501.csv', sep = ',', header = T, row.names = 1, check.names = F, stringsAsFactors = F)

pm2 <- rbind(dat2, pool2) %>% t()
dim(pm2) # 12500  1179
pm2 <- pm2[rowSums(!is.na(pm2)) != 0, colSums(!is.na(pm2)) != 0] %>% .[, info_merge$FileName]
dim(pm2) # 12500  1179







# 
## merge datasets ----
identical(colnames(pm1), info1$FileName) # TRUE
identical(colnames(pm2), info_merge$FileName) # TRUE
info_comb <- list(tphp1 = info1, tphp2 = info_merge)

prots <- Reduce(intersect, list(rownames(pm1), rownames(pm2)))
mat_comb <- cbind(pm1[prots, ], pm2[prots, ])
des <- data.frame(
  ID = c(colnames(pm1), colnames(pm2)),
  sample = 1,
  batch = c(rep(1, ncol(pm1)), rep(2, ncol(pm2)))
) %>% mutate(sample = 1:nrow(.))

# save(pm1, info1, pm2, info_merge, file = 'harmonizR_data.RData')
# load('harmonizR_data.RData')

# save(mat_comb, des, info_comb, file = 'harmonizR_data_combined.RData')
# load('harmonizR_data_combined.RData')



# 1.Datasets integration ----
load('harmonizR_data_combined.RData')

library(HarmonizR)

quantile(mat_comb, na.rm = T)
result <- harmonizR(log2(mat_comb), des, cores = 14,
                    algorithm = 'ComBat', ComBat_mode = 1,
                    output_file = 'PUH_data1_data2_log2_HarmonizR_cured', plot = 'samplemeans')
# ComBat_mode     Corresponding ComBat Arguments
# 1 (default)     par.prior = TRUE, mean.only = FALSE
# 2               par.prior = TRUE, mean.only = TRUE
# 3               par.prior = FALSE, mean.only = FALSE
# 4               par.prior = FALSE, mean.only = TRUE

# save(result, file = 'PUH_data1_data2_log2_HarmonizR_cured.RData')


# 2.Sort information ------------
# info1 <- rio::import('20250725_PUH_sample_information_1900files_v10_addLabels.xlsx')
# info_merge <- rio::import('20250809_PUH_RCA_sample_and_pool_1179files_info.xlsx')
load('harmonizR_data_combined.RData')

info_comb$tphp1[which(info_comb$tphp1$tissue_name == 'prostate cacinoma'), 'tissue_name'] <- 'prostate carcinoma' # correct typo



## uniform and clear info -----
# Chinese cancer name -> standardized cancer name
cn2cancer <- c(
  "胃腺癌" = "gastric adenocarcinoma",
  "大肠粘液腺癌" = "large intestine mucinous adenocarcinoma",
  "肝癌" = "hepatocellular carcinoma",
  "肾透明细胞癌" = "renal clear cell carcinoma",
  "结肠腺癌" = "large intestine adenocarcinoma",
  "宫颈癌" = "cervical carcinoma",
  "胰腺癌" = "pancreatic carcinoma",
  "舌癌" = "tongue carcinoma",
  "前列腺癌" = "prostate carcinoma",
  "食管癌" = "esophagus carcinoma",
  "肺腺癌" = "lung adenocarcinoma",
  "子宫内膜癌" = "endometrial carcinoma",
  "阴茎癌" = "penis carcinoma",
  "胆总管癌" = "cholangiocarcinoma",
  "直肠癌" = "rectum adenocarcinoma",
  "胸腺瘤" = "thymoma and thymic carcinoma",
  "胚胎癌" = "embryonal carcinoma",
  "精原细胞瘤" = "seminoma",
  "胆囊癌" = "gallbladder carcinoma",
  "肝内胆管癌" = "cholangiocarcinoma",
  "输卵管癌" = "fallopian tube carcinoma",
  "胃肠道间质瘤" = "gastrointestinal stromal tumors",
  "甲状腺" = "thyroid carcinoma",
  "喉癌" = "laryngocarcinoma",
  "平滑肌肉瘤" = "leiomyosarcoma",
  "胶质母细胞瘤" = "glioblastoma",
  "弥漫大B" = "diffuse large b-cell lymphoma",
  "反应性增生淋巴结" = "reactive hyperplastic lymph nodes",  # reactive, not malignant
  "卵巢高级别浆液性癌" = "high grade serous ovarian cancer", # HGSOC subtype
  "卵巢透明细胞样癌" = "clear cell ovarian cancer",  # OCCC subtype
  "乳腺癌" = "breast carcinoma",  # OR breast lobular carcinoma?
  "横纹肌肉瘤" = "rhabdomyosarcoma",  # OR pleomorphic rhabdomyosarcoma?
  "膀胱癌" = "bladder carcinoma"  # not existed before
)
# unique(str_subset(info_comb$tphp1$tissue_name[info_comb$tphp1$sample_type == 'T'], 'lymphoma'))

tn2ana <- info_comb$tphp1 %>%
  filter(sample_type %in% c('T', 'NT')) %>%
  distinct(tissue_name, anatomical_classification) %>% 
  drop_na() %>% 
  pull(anatomical_classification, tissue_name) %>% 
  append(c('reactive hyperplastic lymph nodes' = 'lymph node',
           'high grade serous ovarian cancer' = 'ovary',
           'clear cell ovarian cancer' = 'ovary',
           'breast carcinoma' = 'mammary gland',
           'bladder carcinoma' = 'bladder'))
# cn2cancer[which(is.na(tn2ana[cn2cancer]))] # should be character(0)

## artificial checking
# cnmatch <- cbind(cn2cancer, tn2ana[cn2cancer]) %>% as.data.frame() %>%
#   setNames(c('tissue_name', 'anatomical_classification')) %>%
#   rownames_to_column('癌症类型')
# # rio::export(cnmatch, '癌症类型match.xlsx')
# tmp1 <- info_comb$tphp1 %>%
#   filter(sample_type %in% c('T', 'NT')) %>%
#   distinct(tissue_name, anatomical_classification) %>% 
#   rename(anatomical_classification1 = anatomical_classification) %>% 
#   inner_join(cnmatch)



info1 <- info_comb$tphp1 %>% select(FileName, file_id, sample_id, sample_type, tissue_name, anatomical_classification, patient_ID, instrument, date, new) %>% 
  rename(Date = date) %>% 
  mutate(new = ifelse(new == 'old', 'data1', 'data2'))

info_merge <- info_comb$tphp2 %>%
  rename(pathological_ID = 病理号, 
         tissue_name = `癌症类型`, 
         patient_ID = 住院号) %>% 
  mutate(sample_id = str_remove(tmp1$file_id, '_(brep|trep|pool).*$'),
         sample_type = ifelse(!is.na(Cancer_type), Cancer_type, 'p'),
         tissue_name = cn2cancer[tissue_name],
         anatomical_classification = tn2ana[tissue_name],
         instrument = str_extract(FileName, '^[A-Z]+'),
         new = 'data3') %>% 
  select(FileName, file_id, sample_id, sample_type, tissue_name, anatomical_classification, patient_ID, instrument, Date, new)

## artificial checking
View(info1 %>% drop_na() %>% anti_join(info1, .))
View(info_merge %>% drop_na() %>% anti_join(info_merge, .))

info_merge <- rbind(info1, info_merge)
identical(colnames(mat_comb), info_merge$FileName) # TRUE

## unique file_id
tmp1 <- info_merge %>% count(file_id) %>% filter(n > 1) %>% 
  semi_join(info_merge, .) %>% 
  select(sample_id) %>% 
  semi_join(info_comb$tphp1, .) %>%
  group_by(sample_id) %>%
  arrange(DateTime) %>% 
  mutate(Rep = ifelse(is.na(Rep), DIA_ID, Rep)) %>% 
  mutate(Rep_type = ifelse(row_number() == 1, NA, 't'),
         Rep_type = ifelse(DIA_ID != Rep, 'b', Rep_type))

# edit
tmp1.1 <- tmp1 %>% 
  filter(is.na(Rep_type))
tmp1.2 <- tmp1 %>% 
  filter(!is.na(Rep_type)) %>% 
  group_by(sample_id, Rep_type) %>%
  arrange(DateTime) %>%
  mutate(
    seq_num  = row_number(),
    file_id  = str_c(sample_id, '_', Rep_type, seq_num),
    seq_num = NULL, .after = sample_id
  ) %>%
  ungroup()
tmp1.c <- rbind(tmp1.1, tmp1.2) %>% arrange(sample_id, DateTime)

# formatting
tmp1.cf <- info_merge %>% filter(FileName %in% tmp1.c$FileName) %>% 
  select(FileName, Date) %>% 
  inner_join(tmp1.c) %>% 
  select(all_of(colnames(info_merge)))

# corrected full info
info_merge.c <- info_merge %>% filter(!(FileName %in% tmp1.c$FileName)) %>%
  rbind(tmp1.cf) %>% 
  set_rownames(.$FileName) %>% 
  .[info_merge$FileName, ] %>% 
  set_rownames(NULL)

# correct `new`
new.c <- des %>% rename(FileName = ID) %>%
  full_join(info_comb$tphp1 %>% select(FileName, new)) %>% 
  mutate(new = ifelse(is.na(new), 'data3', new),
         new = ifelse(batch == 1 & new == 'old', 'data1', new),
         new = ifelse(batch == 1 & new == 'new', 'data2', new)) %>% 
  select(FileName, new)
info_merge.c %<>% select(-new) %>% inner_join(new.c)

# correct `tissue name`
info_merge.c$tissue_name <- str_remove(info_merge.c$tissue_name, '[\\r\\n]+$')


rio::export(info_merge.c, 'PUH_data1_data2_info_simple_20250812.xlsx')

# also correct info_comb
# info_comb$tphp1 <- info_comb$tphp1 %>% filter(FileName %in% tmp1.c$FileName) %>% 
#   select(FileName) %>% 
#   inner_join(tmp1.c) %>% 
#   rbind(info_comb$tphp1 %>% filter(!(FileName %in% tmp1.c$FileName))) %>% 
#   set_rownames(.$FileName) %>% 
#   .[info_comb$tphp1$FileName, ] %>% 
#   set_rownames(NULL)
info_comb$tphp1 %<>%
  select(-all_of(intersect(colnames(info_comb$tphp1), colnames(info_merge.c)) %>%
                  setdiff('FileName'))) %>% 
  inner_join(info_merge.c) %>% 
  select(all_of(colnames(info_comb$tphp1)))


save(mat_comb, des, info_comb, file = 'harmonizR_data_combined.RData')


## correct info1 version 11 -----
info11 <- rio::import('20250811_PUH_sample_information_1900files_v11.xlsx')
# info11$FileName %>% setdiff(info_merge$FileName)
info11.1 <- info11 %>% filter(!(FileName %in% info_comb$tphp1$FileName))
info11.2 <- info11 %>% filter(FileName %in% info_comb$tphp1$FileName)


which(!sapply(1:ncol(info11.2), function(j){
  identical(info_comb$tphp1[, j], info11.2[, j])
}))
View(data.frame(one = info_comb$tphp1[, 13],
                another = info11.2[, 13]) %>% 
       filter(one != another | xor(is.na(one), is.na(another))))
# remove: 13 14 15 18
# remove: 7  8 11 30

info11.c <- info11.2 %>% select(-c(7, 8, 11, 30)) %>% 
  full_join(info_comb$tphp1 %>% select(-c(13, 14, 15, 18))) %>% 
  rbind(info11.1) %>% 
  set_rownames(.$FileName) %>% 
  .[info11$FileName, ] %>% 
  set_rownames(NULL)
identical(info11.c$FileName, info11$FileName) # TRUE

info11.c$tissue_name <- str_remove(info11.c$tissue_name, '[\\r\\n]+$') # remove \r\n
info11.c[which(info11.c$tissue_name == 'prostate cacinoma'), 'tissue_name'] <- 'prostate carcinoma' # correct typo
info11.c %<>%
  mutate(new = ifelse(new == 'old', 'data1', new),
         new = ifelse(new == 'new', 'data2', new))
rio::export(info11.c, '20250811_PUH_sample_information_1900files_v11.xlsx')


## correct info2 -----
info_merge <- rio::import('PUH_data1_data2_info_simple_20250812.xlsx')
info2 <- rio::import('20250809_PUH_RCA_sample_and_pool_1179files_info.xlsx')

info2$FileName %>% setdiff(info_merge$FileName) # character(0)
colnames(info_merge) %>% setdiff(colnames(info2))
colnames(info_merge) %>% setdiff(colnames(info2), .)
colnames(info_merge) %>% intersect(colnames(info2)) # "FileName" "file_id"  "Date"

info2.c <- info_merge %>% filter(FileName %in% info2$FileName) %>%
  inner_join(info2 %>% select(-all_of(colnames(info_merge) %>% intersect(colnames(info2)) %>% setdiff('FileName')))) %>% 
  set_rownames(.$FileName) %>% 
  .[info2$FileName, ] %>% 
  set_rownames(NULL) %>% 
  select(all_of(colnames(info_merge)), everything())

rio::export(info2.c, '20250812_PUH_RCA_sample_and_pool_1179files_info.xlsx')



# 3.Estimation ------------
info_merge <- rio::import('PUH_data1_data2_info_simple_20250812.xlsx')
dat_merge <- rio::import('PUH_data1_data2_log2_HarmonizR_cured.tsv') %>%
  rename(protein = V1)

pm <- dat_merge %>% column_to_rownames('protein') %>% as.matrix()
dim(pm) # 12019  2990
sum(is.na(pm)) / nrow(pm) / ncol(pm) # 0.4132438
quantile(pm, na.rm = T)
# table(pm < 0) # 7.351111% < 0
# # FALSE     TRUE 
# # 19536081  1550066 
# hist(pm[pm < 0])
quantile(colSums(!is.na(pm)))
quantile(rowSums(!is.na(pm)))

df_est <- pm %>% t() %>% as.data.frame() %>%
  rownames_to_column('FileName') %>% 
  inner_join(info_merge, .)


## hclust ----------
tmp_heat <- pm %>% set_colnames(info_merge$file_id)
tmp_heat[is.na(tmp_heat)] <- 1
tmp_ann <- info_merge %>%
  column_to_rownames('file_id') %>%
  select(tissue_name, anatomical_classification,
         sample_type, instrument, Date, new) %>% 
  mutate(Date = str_extract(Date, '^\\d{4}'))
tmp_ann_clr <- list(
  sample_type = mycolors[1:5] %>% setNames(unique(tmp_ann$sample_type)),
  instrument = mycolors[6:10] %>% setNames(unique(tmp_ann$instrument)),
  new = mycolors[11:13] %>% setNames(unique(tmp_ann$new)),
  Date = mycolors[14:17] %>% setNames(unique(tmp_ann$Date))
)

p <- pheatmap::pheatmap(
  tmp_heat, scale = 'none',
  annotation_col = tmp_ann,
  annotation_colors = tmp_ann_clr,
  clustering_method = 'ward.D2',
  clustering_distance_rows = 'euclidean',
  clustering_distance_cols = 'euclidean',
  show_colnames = F, show_rownames = F,
  filename = 'data1_data2_heatmap.pdf',
  width = 10, height = 5)
pheatmap::pheatmap(
  tmp_heat[1:2, ], scale = 'none',
  annotation_col = tmp_ann,
  annotation_colors = tmp_ann_clr,
  cluster_rows = F, cluster_cols = p$tree_col,
  # clustering_method = 'ward.D2',
  # clustering_distance_rows = 'euclidean',
  # clustering_distance_cols = 'euclidean',
  show_colnames = F, show_rownames = F,
  filename = 'data1_data2_heatmap_legend.pdf',
  width = 10, height = 20)


## dimention reduction ----------
meta_df <- meta2 %>% 
  rename(Gender = 性别, Age = 年龄) %>%
  select(file_id, Batch_m, Date, Gender, Cancer_type) %>% 
  mutate(Date = as.character(Date),
         Batch_m = factor(Batch_m))
id_col <- 'file_id'
# date_col <- 'Age'
var_col <- colnames(meta_df) %>% setdiff(c(id_col, date_col))
seed <- 10
res.dr <- beca.DR(tmp_heat, meta_df, id_col, var_col, seed)












