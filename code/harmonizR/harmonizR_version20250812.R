setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list = ls())
pacman::p_unload(pacman::p_loaded(), character.only = T)
source('../source/my_fun.R')



# 0.Data preparation ---------
## dat1 --------
dat1 <- read.delim('mapped_pg_matrix_1780_14062.csv', sep = ',', header = T, row.names = 1, check.names = F, stringsAsFactors = F)
info1 <- rio::import('20250811_PUH_sample_information_1900files_v11.xlsx')
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
info2 <- rio::import('20250812_PUH_RCA_sample_and_pool_1179files_info.xlsx')
dat2 <- read.delim('mapped_dat2_pg_matrix_1152_12500.csv', sep = ',', header = T, row.names = 1, check.names = F, stringsAsFactors = F)
pool2 <- read.delim('mapped_dat2_pool_matrix_27_12501.csv', sep = ',', header = T, row.names = 1, check.names = F, stringsAsFactors = F)

pm2 <- rbind(dat2, pool2) %>% t()
dim(pm2) # 12500  1179
pm2 <- pm2[rowSums(!is.na(pm2)) != 0, colSums(!is.na(pm2)) != 0] %>% .[, info2$FileName]
dim(pm2) # 12500  1179







# 
## merge datasets ----
identical(colnames(pm1), info1$FileName) # TRUE
identical(colnames(pm2), info2$FileName) # TRUE
info_comb <- list(tphp1 = info1, tphp2 = info2)

prots <- Reduce(intersect, list(rownames(pm1), rownames(pm2)))
mat_comb <- cbind(pm1[prots, ], pm2[prots, ])
des <- data.frame(
  ID = c(colnames(pm1), colnames(pm2)),
  sample = 1,
  batch = c(rep(1, ncol(pm1)), rep(2, ncol(pm2)))
) %>% mutate(sample = 1:nrow(.))

# save(pm1, info1, pm2, info2, file = 'harmonizR_data.RData')
# load('harmonizR_data.RData')

# save(mat_comb, des, info_comb, file = 'harmonizR_data_combined.RData')
# load('harmonizR_data_combined.RData')


## simple data info extraction ------
col_names <- c('FileName', 'file_id', 'sample_id', 'sample_type', 'tissue_name', 'anatomical_classification', 'patient_ID', 'instrument', 'Date', 'new')
info_merge <- info_comb$tphp1 %>%
  mutate(Date = date) %>%
  select(all_of(col_names)) %>% 
  rbind(info_comb$tphp2 %>% select(all_of(col_names)))
rio::export(info_merge, 'PUH_simple_info_2990files_10labels.xlsx')



# 1.Datasets integration ----
load('harmonizR_data_combined.RData')

library(HarmonizR)

dim(mat_comb) # 12019  2990
quantile(colSums(!is.na(mat_comb))) # 497-9347
quantile(rowSums(!is.na(mat_comb))) # 7-2991
quantile(mat_comb, na.rm = T) # 5.8 - 1e12
result <- harmonizR(log2(mat_comb), des, cores = 14,
                    algorithm = 'ComBat', ComBat_mode = 1,
                    output_file = 'PUH_pg_12019_2990_HarmonizR_cured', plot = 'samplemeans')
# ComBat_mode     Corresponding ComBat Arguments
# 1 (default)     par.prior = TRUE, mean.only = FALSE
# 2               par.prior = TRUE, mean.only = TRUE
# 3               par.prior = FALSE, mean.only = FALSE
# 4               par.prior = FALSE, mean.only = TRUE

# save(result, file = 'PUH_pg_12019_2990_HarmonizR_cured.RData')



# 2.Estimation ------------
# rm(list = ls())
# pacman::p_unload(pacman::p_loaded(), character.only = T)
# source('../source/my_fun.R')


info <- rio::import('PUH_simple_info_2990files_10labels.xlsx')
dat <- rio::import('PUH_pg_12019_2990_HarmonizR_cured.tsv') %>%
  rename(protein = V1)

pm <- dat %>% column_to_rownames('protein') %>%
  as.matrix() %>% set_colnames(info$file_id)
dim(pm) # 12019  2990
sum(is.na(pm)) / nrow(pm) / ncol(pm) # 0.4132438
quantile(pm, na.rm = T) # [4, 33]

quantile(colSums(!is.na(pm)))
quantile(rowSums(!is.na(pm)))

## NA impute -------
pm_fill <- pm
pm_fill[is.na(pm_fill)] <- 0


## metadata ----
ann <- info %>%
  column_to_rownames('file_id') %>%
  select(tissue_name, anatomical_classification,
         sample_type, instrument, Date, new) %>% 
  mutate(Date = str_extract(Date, '^\\d{4}'))

ann_clr <- list(
  sample_type = mycolors[1:5] %>% setNames(unique(ann$sample_type)),
  instrument = mycolors[6:10] %>% setNames(unique(ann$instrument)),
  new = mycolors[11:13] %>% setNames(unique(ann$new)),
  Date = mycolors[14:17] %>% setNames(unique(ann$Date))
)

# df_est <- pm %>% t() %>% as.data.frame() %>%
#   rownames_to_column('FileName') %>% 
#   inner_join(info, .)


## hclust ----------


p <- pheatmap::pheatmap(
  pm_fill, scale = 'none',
  annotation_col = ann,
  annotation_colors = ann_clr,
  clustering_method = 'ward.D2',
  clustering_distance_rows = 'euclidean',
  clustering_distance_cols = 'euclidean',
  show_colnames = F, show_rownames = F,
  filename = 'data1_data2_heatmap.pdf',
  width = 10, height = 5)
pheatmap::pheatmap(
  pm_fill[1:2, ], scale = 'none',
  annotation_col = ann,
  annotation_colors = ann_clr,
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
res.dr <- beca.DR(pm_fill, meta_df, id_col, var_col, seed)












