rstudioapi::getActiveDocumentContext()$path
setwd('~/GitHub/TPHP/code/BECA/')
source('../../source/my_fun.R')
library(tidyverse)
library(magrittr)
library(RColorBrewer)
library(scales)
library(pheatmap)
library(Rtsne)
library(umap)
# library(broom)

# 0.Input ----------
## data matrix and sample annotation ------
df <- rio::import('//172.16.13.136/TPHP/results/rlt_combine_entraplib_G3/mapped_pg_matrix_1780_13964.csv')
df_filter <- rio::import('//172.16.13.136/TPHP/results/rlt_combine_entraplib_G3/data_matrix_filtered_05NA_1899_13516.csv') %>% mutate(V1 = str_replace_all(V1, '\\.', '-'))
setdiff(df_filter$V1, df$V1) # pooling samples

# df_q_log2 <- rio::import('//172.16.13.136/TPHP/results/rlt_combine_entraplib_G3/data_matrix_q_log.csv') %>% mutate(V1 = str_replace_all(V1, '\\.', '-'))
# info <- rio::import('//172.16.13.136/TPHP/results/rlt_combine_entraplib_G3/20250714_PUH_sample_information_1781files_info_edited_v9.xlsx')
info <- rio::import('//172.16.13.136/TPHP/results/rlt_combine_entraplib_G3/20250725_PUH_sample_information_1781files_info_edited_v10.xlsx')

info_datetime <- rio::import('~/GitHub/TPHP/code/datetime_extract/rawdata_datetime.xlsx')

identical(df$V1, info$FileName) # TRUE

info_full <- info %>% filter(!is.na(anatomical_classification)) %>% 
  full_join(info_datetime)
info %<>% filter(!is.na(anatomical_classification)) %>% 
  inner_join(info_datetime)

df %<>% filter(V1 %in% info$FileName)



## filter out decoys -----
matpg <- df_filter %>% column_to_rownames('V1') %>% 
  select(matches('^[A-Z0-9]+$')) %>% t()
dim(matpg) # 13516  1899
min(colSums(!is.na(matpg))) # 551
min(rowSums(!is.na(matpg))) # 1


info_fit <- info_of_pm(matpg) %>% rename(FileName = Filename) %>% 
  full_join(info_full, .) %>%
  mutate(
    sample_type = c('Normal' = 'N', 'adjacent' = 'NT', 'carcinoma' = 'T', 'Fetal' = 'F')[sample_type],
    date.in.filename = date,
    date = str_remove_all(str_extract(DateTime, '[^_]+'), '-'),
    year = str_sub(date, 1, 4),
    month = str_sub(date, 5, 6)
  ) %>% 
  mutate(
    Is.Rep = !is.na(Rep_type),
    Rep_type = ifelse(Rep_type == 'biological replicates', 'b', ifelse(Rep_type == 'technical replicates', 't', NA)),
    sample_id = ifelse(!Is.Rep, DIA_ID, Rep),
    .after = Rep_type
  )

# append missing replicates
rep_sup <- info_fit %>%
  filter(!Is.Rep, !is.na(sample_id)) %>%
  count(sample_id) %>% filter(n > 1) %>%
  semi_join(info_fit, .) %>% 
  filter(!Is.Rep) %>% pull(FileName)
info_fit %<>% mutate(
  Rep = ifelse(FileName %in% rep_sup, sample_id, Rep),
  Is.Rep = ifelse(FileName %in% rep_sup, T, Is.Rep),
  Rep_type = ifelse(FileName %in% rep_sup, 't', Rep_type)
)

# remove wrong replicates
rep_smp <- info_fit %>% filter(Is.Rep) %>% pull(sample_id) %>% unique()
rep_smp_wrong <- info_fit %>% filter(sample_id %in% rep_smp) %>% count(sample_id) %>% filter(n == 1) %>% pull(sample_id)
info_fit %<>% mutate(
  Rep = ifelse(sample_id %in% rep_smp_wrong, NA, Rep),
  Is.Rep = ifelse(sample_id %in% rep_smp_wrong, F, Is.Rep),
  Rep_type = ifelse(sample_id %in% rep_smp_wrong, NA, Rep_type)
)

# info_fit %>% filter(Is.Rep) %>%
#   # mutate(suffix = as.numeric(str_extract(FileName, '\\d+$'))) %>% 
#   group_by(sample_id) %>%
#   # arrange(suffix) %>% 
#   arrange(DateTime) %>% 
#   mutate(file_id = str_c(sample_id, 1:nrow(.)), ifelse(Rep_type == 'biological replicates', 'b', 't'))

info_fit1 <- info_fit %>% filter(Is.Rep) %>%
  group_by(sample_id, Rep_type) %>%
  arrange(DateTime) %>%
  mutate(
    seq_num  = row_number(),
    file_id  = str_c(sample_id, '_', Rep_type, seq_num),
    seq_num = NULL, .after = sample_id
  ) %>%
  ungroup()
info_fit2 <- info_fit %>% filter(!Is.Rep) %>% mutate(file_id = sample_id)

info_fitc <- info_fit1 %>% rbind(info_fit2) %>% as.data.frame() %>%
  set_rownames(.$FileName) %>% .[colnames(matpg), ] %>% 
  remove_rownames() %>%
  mutate(file_id = ifelse(is.na(file_id), str_c('pool', str_extract(FileName, '^[A-Z]+'), str_extract(FileName, '\\d+$')), file_id),
         sample_type = ifelse(is.na(sample_type), 'p', sample_type)) # for pooling samples


rio::export(info_fitc, '~/GitHub/TPHP/input/20250729_1899files_add_fields.xlsx')

identical(colnames(matpg), info_fitc$FileName) # TRUE
info_fitc %>% count(file_id) %>% count(n) %>% nrow() # 1

colnames(matpg) <- info_fitc$file_id
metadata <- info_fitc %>%
  select(file_id, sample_id, Rep, Rep_type, sample_type:patient_ID, DateTime:last_col())

# info_fitc %>% count(sample_type)
# info_fitc %>% count(tissue_name)
# info_fitc %>% filter(sample_type %in% c('T', 'NT')) %>% count(tissue_name)
# info_fitc %>% count(Rep_type)


## filter out low-quanlity samples -----
rm.fid <- c()

# setdiff(colnames(df), colnames(df_filter))
# colSums(!is.na(df[, -1]))[setdiff(colnames(df), colnames(df_filter))]

pm_filter <- df_filter %>% column_to_rownames('V1') %>% t()
pm_log2 <- log2(pm_filter)
pm_q <- preprocessCore::normalize.quantiles(pm_filter, copy = T) %>%
  set_colnames(colnames(pm_filter)) %>%
  set_rownames(rownames(pm_filter))
pm_q_log2 <- log2(pm_q)




# df_est <- matpg %>% t() %>% as.data.frame() %>% 
#   rownames_to_column('file_id') %>% 
#   inner_join(metadata %>% select(file_id:anatomical_classification), .)
df_est <- pm_filter %>% t() %>% as.data.frame() %>%
  rownames_to_column('FileName') %>% 
  inner_join(info_fitc %>% select(FileName, file_id), .) %>% 
  mutate(file_id = ifelse(is.na(file_id), str_c('pool', str_extract(FileName, '^[A-Z]+'), str_extract(FileName, '\\d+$')), file_id)) %>% 
  select(-FileName) %>% 
  full_join(metadata %>% select(file_id:anatomical_classification), .) %>% 
  mutate(sample_id = ifelse(is.na(sample_id), 'pool', sample_id))

df_est_qlog <- pm_q_log2 %>% t() %>% as.data.frame() %>%
  rownames_to_column('FileName') %>% 
  inner_join(info_fitc %>% select(FileName, file_id), .) %>% 
  select(-FileName) %>% 
  inner_join(metadata %>% select(file_id:anatomical_classification), .)


### replicates check -----
### replicates identity
rep.identity <- df_est %>% filter(!is.na(Rep_type)) %>%
  select(sample_id) %>% 
  semi_join(df_est, .) %>% 
  plyr::ddply('sample_id', function(dfsub){
    cat(dfsub$sample_id[1], '...\r')
    X <- dfsub %>%
      column_to_rownames('file_id') %>%
      select(-(sample_id:anatomical_classification)) %>% t()
    ret <- data.frame(file_id = dfsub$file_id,
                      `# proteins` = colSums(!is.na(X)),
                      check.names = F)
    ret$outlier.lower.ingroup <- get_outliers(ret$`# proteins`)[1]
    return(ret)
  })
rep.identity %<>%
  group_by(sample_id) %>%
  summarise(identity.mean = mean(`# proteins`)) %>%
  inner_join(rep.identity, .) %>%
  arrange(identity.mean) %>%
  mutate(sample_id = factor(sample_id, levels = unique(sample_id)),
         Is.Lower.Ingroup = `# proteins` < outlier.lower.ingroup) %>% 
  inner_join(info_fitc %>% select(file_id, Rep_type))

### correlation
rep.pearson <- df_est %>% filter(!is.na(Rep_type)) %>%
  select(sample_id) %>% 
  semi_join(df_est, .) %>% 
  plyr::ddply('sample_id', function(dfsub){
    cat(dfsub$sample_id[1], '...\r')
    X <- dfsub %>%
      column_to_rownames('file_id') %>%
      select(-(sample_id:anatomical_classification)) %>% t()
    corX <- cor(X, use = 'pairwise.complete.obs', method = 'pearson')
    corX_pair <- corX %>% as.data.frame() %>%
      rownames_to_column('ID1') %>%
      pivot_longer(cols = -ID1, names_to = 'ID2') %>%
      unite('pair', ID1, ID2, sep = '-') %>%
      pull(pair) %>%
      matrix(nrow = nrow(corX), ncol = ncol(corX),
             byrow = T, dimnames = dimnames(corX))
    ret <- data.frame(ID.pair = corX_pair %>% .[upper.tri(., diag = F)],
                      pearson.r = corX %>% .[upper.tri(., diag = F)],
                      check.names = F)
    ret$outlier.lower.ingroup <- get_outliers(ret$pearson.r)[1]
    return(ret)
  })
rep.pearson %<>%
  group_by(sample_id) %>%
  summarise(pearson.r.mean = mean(pearson.r)) %>%
  inner_join(rep.pearson, .) %>%
  arrange(pearson.r.mean) %>%
  mutate(sample_id = factor(sample_id, levels = unique(sample_id)))

### intensity
rep.quantity <- df_est %>% filter(!is.na(Rep_type)) %>%
  select(sample_id) %>% 
  semi_join(df_est, .) %>% 
  plyr::ddply('sample_id', function(dfsub){
    cat(dfsub$sample_id[1], '...\r')
    X <- dfsub %>%
      column_to_rownames('file_id') %>%
      select(-(sample_id:anatomical_classification)) %>% t()
    ret <- data.frame(file_id = dfsub$file_id,
                      Total.Quantity = log2(colSums(X, na.rm = T)),
                      check.names = F)
    ret$outlier.lower.ingroup <- get_outliers(ret$Total.Quantity)[1]
    return(ret)
  })
rep.quantity %<>%
  group_by(sample_id) %>%
  summarise(Total.Quantity.mean = mean(Total.Quantity)) %>%
  inner_join(rep.quantity, .) %>%
  arrange(Total.Quantity.mean) %>%
  mutate(sample_id = factor(sample_id, levels = unique(sample_id))) %>% 
  inner_join(info_fitc %>% select(file_id, Rep_type))


### figures
sid_order <- rep.identity %>% arrange(`# proteins`) %>% distinct(sample_id) %>% pull()
plot_identity <- rep.identity %>% 
  mutate(sample_id = factor(sample_id, sid_order)) %>% 
  ggplot() +
  aes(x = `# proteins`, y = sample_id) +
  geom_boxplot(color = '#000000', outlier.color = '#c23190', outlier.size = 3) +
  geom_point(aes(color = Rep_type), size = 1) +
  labs(x = "# proteins", y = "Sample ID", subtitle = "Protein identification") +
  ggsci::scale_color_npg() +
  theme_bw() +
  theme(text = element_text(size = 10))
plot_quantity <- rep.quantity %>% 
  mutate(sample_id = factor(sample_id, sid_order)) %>% 
  ggplot() +
  aes(x = Total.Quantity, y = sample_id) +
  geom_boxplot(color = '#000000', outlier.color = '#c23190', outlier.size = 3) +
  geom_point(aes(color = Rep_type), size = 1) +
  labs(x = "Log2(total quantity)", y = "Sample ID", subtitle = "Protein quantification") +
  ggsci::scale_color_npg() +
  theme_bw() +
  theme(text = element_text(size = 10))
plot_pearson <- ggplot(rep.pearson) +
  aes(x = pearson.r, y = sample_id) +
  geom_boxplot(color = '#000000', outlier.color = '#c23190', outlier.size = 3) +
  geom_point(color = '#000000', size = 1) +
  labs(x = "Pearson's r", y = "Sample ID", subtitle = "Correlation - same.sample.id") +
  theme_bw() +
  theme(text = element_text(size = 10))

ggsave(filename = 'replicates_quality_check.pdf',
       ggpubr::ggarrange(plot_identity, plot_quantity, plot_pearson,
                         nrow = 1, ncol = 3, common.legend = T),
       width = 12, height = 12)


rm.fid.cor <- c('DIA_594_b1')
rm.fid %<>% union(rm.fid.cor)


### pooling check ---------
pool.identity <- df_est %>% filter(sample_id == 'pool') %>%
  plyr::ddply('sample_id', function(dfsub){
    cat(dfsub$sample_id[1], '...\r')
    X <- dfsub %>%
      column_to_rownames('file_id') %>%
      select(-(sample_id:anatomical_classification)) %>% t()
    ret <- data.frame(file_id = dfsub$file_id,
                      `# proteins` = colSums(!is.na(X)),
                      check.names = F)
    ret$outlier.lower.ingroup <- get_outliers(ret$`# proteins`)[1]
    return(ret)
  })
pool.identity %<>%
  group_by(sample_id) %>%
  summarise(identity.mean = mean(`# proteins`)) %>%
  inner_join(pool.identity, .) %>%
  arrange(identity.mean) %>%
  mutate(sample_id = factor(sample_id, levels = unique(sample_id)),
         Is.Lower.Ingroup = `# proteins` < outlier.lower.ingroup) %>% 
  mutate(Instrument = str_extract(file_id, '[A-Z]+'))


p.pool_identity <- pool.identity %>% 
  ggplot() +
  aes(x = `# proteins`, y = sample_id, color = Instrument) +
  geom_boxplot(width = 0.2, outlier.size = 0.5, position = position_dodge(width = 0.5)) +
  geom_boxplot(color = 'black', width = 0.05, outlier.size = 0.5) +
  geom_point(size = 1) +
  labs(x = "# proteins", y = "", subtitle = "Protein identification") +
  ggsci::scale_color_npg() +
  theme_bw() +
  theme(text = element_text(size = 10))
ggsave('dia_pool_quality_check.pdf', p.pool_identity, width = 5.5, height = 4.5)



### samples check -----
### identity
est.identity <- df_est %>% filter(sample_type != 'p') %>% 
  plyr::ddply(c('sample_type', 'tissue_name'), function(dfsub){
  cat(dfsub$sample_type[1], dfsub$tissue_name[1], '...\r')
  X <- dfsub %>%
    column_to_rownames('file_id') %>%
    select(-(sample_id:Rep_type), -(tissue_name:anatomical_classification)) %>% t()
  ret <- data.frame(file_id = dfsub$file_id,
                    `# proteins` = colSums(!is.na(X)),
                    check.names = F)
  ret$outlier.lower.ingroup <- get_outliers(ret$`# proteins`)[1]
  return(ret)
})
est.identity %<>%
  group_by(sample_type) %>%
  summarise(identity.mean = mean(`# proteins`)) %>%
  inner_join(est.identity, .) %>%
  arrange(identity.mean) %>%
  mutate(sample_type = factor(sample_type, levels = unique(sample_type)),
         Is.Lower.Ingroup = `# proteins` < outlier.lower.ingroup)
est.identity <- rbind(
  est.identity %>%
    mutate(sample_type = 'total', tissue_name = 'total', identity.mean = mean(`# proteins`)),
  est.identity
) %>%
  mutate(y = str_c(sample_type, '-', tissue_name),
         y = factor(y, unique(y)),
         outlier.lower.total = get_outliers(`# proteins`)[1],
         Is.Lower.Ingroup = `# proteins` < outlier.lower.ingroup,
         Is.Lower.Total = `# proteins` < outlier.lower.total)

est.identity.low <- est.identity %>%
  filter(sample_type != 'total', tissue_name != 'total') %>% 
  filter(Is.Lower.Ingroup) %>% 
  distinct(sample_type, tissue_name) %>% 
  semi_join(est.identity, .)
rm.fid.ident <- est.identity.low %>% filter(Is.Lower.Ingroup, `# proteins` < 4000) %>% pull(file_id)
rm.fid %<>% union(rm.fid.ident)


### correlation
est.pearson <- df_est %>%
  filter(sample_type != 'p',
         !(file_id %in% rm.fid)) %>%
  count(sample_type, tissue_name) %>% filter(n > 1) %>%
  semi_join(df_est, .) %>% filter(!(file_id %in% rm.fid)) %>% 
  plyr::ddply(c('sample_type', 'tissue_name'), function(dfsub){
  cat(dfsub$sample_type[1], dfsub$tissue_name[1], '...\r')
  X <- dfsub %>%
    column_to_rownames('file_id') %>%
    select(-(sample_id:Rep_type), -(sample_type:anatomical_classification)) %>% t()
  corX <- cor(X, use = 'pairwise.complete.obs', method = 'pearson')
  corX_pair <- corX %>% as.data.frame() %>%
    rownames_to_column('ID1') %>%
    pivot_longer(cols = -ID1, names_to = 'ID2') %>%
    unite('pair', ID1, ID2, sep = '-') %>%
    pull(pair) %>%
    matrix(nrow = nrow(corX), ncol = ncol(corX),
           byrow = T, dimnames = dimnames(corX))
  ret <- data.frame(ID.pair = corX_pair %>% .[upper.tri(., diag = F)],
                    pearson.r = corX %>% .[upper.tri(., diag = F)],
                    check.names = F)
  ret$outlier.lower.ingroup <- get_outliers(ret$pearson.r)[1]
  return(ret)
})
est.pearson %<>%
  group_by(sample_type, tissue_name) %>%
  summarise(pearson.r.mean = mean(pearson.r)) %>%
  inner_join(est.pearson, .) %>%
  arrange(pearson.r.mean) %>%
  mutate(y = str_c(sample_type, '-', tissue_name),
         y = factor(y, unique(y)))

est.pearson_log <- df_est %>%
  filter(sample_type != 'p',
         !(file_id %in% rm.fid)) %>%
  count(sample_type, tissue_name) %>% filter(n > 1) %>%
  semi_join(df_est, .) %>% filter(!(file_id %in% rm.fid)) %>% 
  mutate_if(is.numeric, log2) %>% 
  plyr::ddply(c('sample_type', 'tissue_name'), function(dfsub){
    cat(dfsub$sample_type[1], dfsub$tissue_name[1], '...\r')
    X <- dfsub %>%
      column_to_rownames('file_id') %>%
      select(-(sample_id:Rep_type), -(sample_type:anatomical_classification)) %>% t()
    corX <- cor(X, use = 'pairwise.complete.obs', method = 'pearson')
    corX_pair <- corX %>% as.data.frame() %>%
      rownames_to_column('ID1') %>%
      pivot_longer(cols = -ID1, names_to = 'ID2') %>%
      unite('pair', ID1, ID2, sep = '-') %>%
      pull(pair) %>%
      matrix(nrow = nrow(corX), ncol = ncol(corX),
             byrow = T, dimnames = dimnames(corX))
    ret <- data.frame(ID.pair = corX_pair %>% .[upper.tri(., diag = F)],
                      pearson.r = corX %>% .[upper.tri(., diag = F)],
                      check.names = F)
    ret$outlier.lower.ingroup <- get_outliers(ret$pearson.r)[1]
    return(ret)
  })
est.pearson_log %<>%
  group_by(sample_type, tissue_name) %>%
  summarise(pearson.r.mean = mean(pearson.r)) %>%
  inner_join(est.pearson_log, .) %>%
  arrange(pearson.r.mean) %>%
  mutate(y = str_c(sample_type, '-', tissue_name),
         y = factor(y, unique(y)))


est.pearson_qlog <- df_est_qlog %>%
  filter(sample_type != 'p',
         !(file_id %in% rm.fid)) %>%
  count(sample_type, tissue_name) %>% filter(n > 1) %>%
  semi_join(df_est_qlog, .) %>% filter(!(file_id %in% rm.fid)) %>% 
  plyr::ddply(c('sample_type', 'tissue_name'), function(dfsub){
    cat(dfsub$sample_type[1], dfsub$tissue_name[1], '...\r')
    X <- dfsub %>%
      column_to_rownames('file_id') %>%
      select(-(sample_id:Rep_type), -(sample_type:anatomical_classification)) %>% t()
    corX <- cor(X, use = 'pairwise.complete.obs', method = 'pearson')
    corX_pair <- corX %>% as.data.frame() %>%
      rownames_to_column('ID1') %>%
      pivot_longer(cols = -ID1, names_to = 'ID2') %>%
      unite('pair', ID1, ID2, sep = '-') %>%
      pull(pair) %>%
      matrix(nrow = nrow(corX), ncol = ncol(corX),
             byrow = T, dimnames = dimnames(corX))
    ret <- data.frame(ID.pair = corX_pair %>% .[upper.tri(., diag = F)],
                      pearson.r = corX %>% .[upper.tri(., diag = F)],
                      check.names = F)
    ret$outlier.lower.ingroup <- get_outliers(ret$pearson.r)[1]
    return(ret)
  })
est.pearson_qlog %<>%
  group_by(sample_type, tissue_name) %>%
  summarise(pearson.r.mean = mean(pearson.r)) %>%
  inner_join(est.pearson_qlog, .) %>%
  arrange(pearson.r.mean) %>%
  mutate(y = str_c(sample_type, '-', tissue_name),
         y = factor(y, unique(y)))

# 
# ### intensity
# rep.quantity <- df_est %>% filter(!is.na(Rep_type)) %>%
#   select(sample_type) %>%
#   semi_join(df_est, .) %>%
#   plyr::ddply(c('sample_type', 'tissue_name'), function(dfsub){
#     cat(dfsub$sample_type[1], '...\r')
#     X <- dfsub %>%
#       column_to_rownames('file_id') %>%
#       select(-(sample_type:anatomical_classification)) %>% t()
#     ret <- data.frame(file_id = dfsub$file_id,
#                       Total.Quantity = log2(colSums(X, na.rm = T)),
#                       check.names = F)
#     ret$outlier.lower.ingroup <- get_outliers(ret$Total.Quantity)[1]
#     return(ret)
#   })
# rep.quantity %<>%
#   group_by(sample_type) %>%
#   summarise(Total.Quantity.mean = mean(Total.Quantity)) %>%
#   inner_join(rep.quantity, .) %>%
#   arrange(Total.Quantity.mean) %>%
#   mutate(sample_type = factor(sample_type, levels = unique(sample_type))) %>%
#   inner_join(info_fitc %>% select(file_id, Rep_type))


### figures
plot_identity <- ggplot(est.identity) +
  aes(x = `# proteins`, y = y) +
  geom_boxplot(color = '#000000', outlier.color = '#c23190', outlier.size = 3) +
  geom_point(color = '#000000', size = 1) +
  labs(x = "# proteins", y = "Group", subtitle = "Protein identification") +
  theme_bw() +
  theme(text = element_text(size = 10))
ggsave('dia_quality_check.pdf', plot_identity, width = 10, height = 30)


plot_identity_low <- ggplot(est.identity.low) +
  aes(x = `# proteins`, y = y, group = y) +
  geom_boxplot(color = '#000000', outlier.shape = NA, position = position_dodge(width = 0.6)) +
  geom_jitter(data = est.identity.low %>% filter(!Is.Lower.Ingroup),
              color = '#000000', size = 1, position = position_dodge(width = 0.6)) +
  geom_point(data = est.identity.low %>% filter(Is.Lower.Ingroup, `# proteins` < 4000),
             color = '#c23190', size = 1) +
  annotate('text', x = 2500, y = 33, size = 3.5,
           label = str_glue('# outliers = {sum(est.identity.low$Is.Lower.Ingroup & est.identity.low$`# proteins` < 4000)}')) +
  labs(x = "# proteins", y = "Group", subtitle = "Protein identification") +
  theme_bw() +
  theme(text = element_text(size = 10))
ggsave('dia_quality_check_low.pdf', plot_identity_low, width = 6, height = 5)




# plot_quantity <- rep.quantity %>% 
#   mutate(sample_type = factor(sample_type, sid_order)) %>% 
#   ggplot() +
#   aes(x = Total.Quantity, y = sample_type) +
#   geom_boxplot(color = '#000000', outlier.color = '#c23190', outlier.size = 3) +
#   geom_point(aes(color = Rep_type), size = 1) +
#   labs(x = "Log2(total quantity)", y = "Group", subtitle = "Protein quantification") +
#   ggsci::scale_color_npg() +
#   theme_bw() +
#   theme(text = element_text(size = 10))
plot_pearson <- ggplot(est.pearson) +
  aes(x = pearson.r, y = y) +
  geom_boxplot(color = '#000000', outlier.color = '#c23190', outlier.size = 3) +
  geom_point(color = '#000000', size = 1) +
  labs(x = "Pearson's r", y = "Group", subtitle = "Correlation") +
  theme_bw() +
  theme(text = element_text(size = 10))
ggsave('dia_quality_pearson_check.pdf', plot_pearson, width = 10, height = 30)
plot_pearson_log <- ggplot(est.pearson_log) +
  aes(x = pearson.r, y = y) +
  geom_boxplot(color = '#000000', outlier.color = '#c23190', outlier.size = 3) +
  geom_point(color = '#000000', size = 1) +
  labs(x = "Pearson's r", y = "Group", subtitle = "Correlation") +
  theme_bw() +
  theme(text = element_text(size = 10))
ggsave('dia_quality_pearson_check_log2.pdf', plot_pearson_log, width = 10, height = 30)
plot_pearson_qlog <- ggplot(est.pearson_qlog) +
  aes(x = pearson.r, y = y) +
  geom_boxplot(color = '#000000', outlier.color = '#c23190', outlier.size = 3) +
  geom_point(color = '#000000', size = 1) +
  labs(x = "Pearson's r", y = "Group", subtitle = "Correlation") +
  theme_bw() +
  theme(text = element_text(size = 10))
ggsave('dia_quality_pearson_check_q_log2.pdf', plot_pearson_qlog, width = 10, height = 30)
est.pearson_compare <- est.pearson %>% select(y, ID.pair, pearson.r) %>%
  inner_join(est.pearson_log %>%
               select(y, ID.pair, pearson.r) %>%
               rename(pearson.r_log = pearson.r)) %>% 
  inner_join(est.pearson_qlog %>%
               select(y, ID.pair, pearson.r) %>%
               rename(pearson.r_qlog = pearson.r))
tmp_p1 <- est.pearson_compare %>%
  pivot_longer(cols = -c('ID.pair', 'y', 'pearson.r'),
               names_to = 'Normalization', values_to = 'pearson.r_log') %>% 
  mutate(Normalization = c(pearson.r_log = 'None', pearson.r_qlog = 'Quantile-norm')[Normalization]) %>% 
  ggplot() +
  aes(x = pearson.r, y = pearson.r_log, color = Normalization) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = 'black') +
  labs(x = 'Pearson.r (RawIntensity)', y = 'Pearson.r (Log2Intensity)') +
  # scale_x_continuous(limits = c(-1, 1)) +
  # scale_y_continuous(limits = c(-1, 1)) +
  ggsci::scale_color_startrek() +
  theme_bw() +
  theme(text = element_text(size = 10), legend.position = 'bottom')
tmp_p2 <- ggExtra::ggMarginal(tmp_p1, type = 'density', margins = 'both',
                              groupColour = T, groupFill = T, alpha = 0.3)
ggsave('dia_quality_pearson_check_compare_normalize.pdf', tmp_p2, width = 5, height = 5.5)


# 
# ggsave(filename = 'sample_quality_check.pdf',
#        ggpubr::ggarrange(plot_identity, plot_quantity, plot_pearson,
#                          nrow = 1, ncol = 3, common.legend = T),
#        width = 12, height = 12)



## filter out decoys -----
info_fitc %>%
  filter(sample_type %in% c('T', 'NT')) %>% 
  count(patient_ID) %>% 
  count(n) %>%
  filter(n == 1) %>% pull(nn) # 12 patients

info_fitc %>%
  filter(sample_type %in% c('T', 'NT'),
         !(file_id %in% rm.fid)) %>%
  count(patient_ID) %>% 
  count(n) %>%
  filter(n == 1) %>% pull(nn) # 18 patients


# 1.Diagnosis of batch effect --------
fn2fid <- info_fitc$file_id %>% setNames(info_fitc$FileName)
pm_list <- list()
pm_list$original <- pm_q_log2 %>% set_colnames(fn2fid[colnames(pm_q_log2)])

# mds <- naimpute.methods[c(1:4, 12:13, 5:6)]
mds <- naimpute.methods[c(1:4)]
for(i in seq_along(mds)){
  md <- mds[i]
  cat(i, md, '...\n')
  pm_list[[i+1]] <- as.matrix(naimpute(pm_list$original, md))
}
names(pm_list)[seq_along(pm_list)[-1]] <- mds

# save(pm_list, file = 'pm_list.RData')



## Dimension reduction ------

beca.DR <- function(X, meta_df, id_col, var_col, date_col = NULL, seed = NA){
  # Transpose to prepare data for PCA/t-SNE/UMAP (samples as rows, features as columns)
  # each row = one sample, each column = a feature
  t_expr <- t(X)
  
  # PCA computation
  print('PCA computation...')
  pca_res <- prcomp(t_expr, center = TRUE, scale. = FALSE)
  pca_df <- data.frame(PC1 = pca_res$x[,1], PC2 = pca_res$x[,2]) %>% 
    rownames_to_column(id_col) %>% 
    inner_join(meta_df)
    
  # Plot PCA
  plots_pca <- lapply(seq_along(var_col), function(i){
    ggplot(pca_df) +
      aes(x = PC1, y = PC2) +
      geom_point(aes_string(color = var_col[i]), size = 2, alpha = 0.8) +
      labs(title = str_c("PCA: Samples colored by ", var_col[i]),
           x = str_c("PC1 (", round(summary(pca_res)$importance[2,1]*100,1), "%)"),
           y = str_c("PC2 (", round(summary(pca_res)$importance[2,2]*100,1), "%)")) +
      scale_color_manual(values = mycolors) +
      theme_bw() +
      theme(text = element_text(size = 10))
  }) %>% setNames(var_col)
  if(!is.null(date_col)){
    plots_pca$DateTime <- ggplot(pca_df) +
      aes(x = PC1, y = PC2) +
      geom_point(aes_string(color = date_col), size = 2, alpha = 0.8) +
      labs(title = str_c("PCA: Samples colored by ", date_col),
           x = str_c("PC1 (", round(summary(pca_res)$importance[2,1]*100,1), "%)"),
           y = str_c("PC2 (", round(summary(pca_res)$importance[2,2]*100,1), "%)")) +
      scale_color_viridis_c(
        option = 'G', begin = 0.05, end = 0.95,
        name   = "Sample date",
        breaks  = as.numeric(date_breaks("2 month")(range(pca_df$DateTime))),
        labels  = date_format("%Y-%m-%d")(date_breaks("2 month")(range(pca_df$DateTime)))
      ) +
      theme_bw() +
      theme(text = element_text(size = 10))
  }
  # ggpubr::ggarrange(plotlist = plots_pca)
  
  # t-SNE computation (perplexity adjusted for dataset size)
  print('t-SNE computation...')
  if(is.integer(seed)) set.seed(seed)  # for reproducibility
  tsne_res <- Rtsne(t_expr, dims = 2, perplexity = 30, verbose = FALSE)
  tsne_df <- data.frame(Dim1 = tsne_res$Y[,1], Dim2 = tsne_res$Y[,2], row.names = rownames(t_expr)) %>% 
    rownames_to_column(id_col) %>% 
    inner_join(meta_df)
  # Plot t-SNE
  plots_tsne <- lapply(seq_along(var_col), function(i){
    ggplot(tsne_df) +
      aes(x = Dim1, y = Dim2) +
      geom_point(aes_string(color = var_col[i]), size = 2, alpha = 0.8) +
      labs(title = str_c("t-SNE: Samples colored by ", date_col), x = "t-SNE Dim1", y = "t-SNE Dim2") +
      scale_color_manual(values = mycolors) +
      theme_bw() +
      theme(text = element_text(size = 10))
  }) %>% setNames(var_col)
  if(!is.null(date_col)){
    plots_tsne$DateTime <- ggplot(tsne_df) +
      aes(x = Dim1, y = Dim2) +
      geom_point(aes_string(color = date_col), size = 2, alpha = 0.8) +
      labs(title = str_c("t-SNE: Samples colored by ", date_col), x = "t-SNE Dim1", y = "t-SNE Dim2") +
      scale_color_viridis_c(
        option = 'G', begin = 0.05, end = 0.95,
        name   = "Sample date",
        breaks  = as.numeric(date_breaks("2 month")(range(pca_df$DateTime))),
        labels  = date_format("%Y-%m-%d")(date_breaks("2 month")(range(pca_df$DateTime)))
      ) +
      theme_bw() +
      theme(text = element_text(size = 10))
  }
  # ggpubr::ggarrange(plotlist = plots_tsne)
  
  # UMAP computation
  print('UMAP computation...')
  if(is.integer(seed)) set.seed(seed)  # for reproducibility
  umap_res <- umap(t_expr)             # default 2D UMAP
  umap_df <- data.frame(UM1 = umap_res$layout[,1], UM2 = umap_res$layout[,2]) %>% 
    rownames_to_column(id_col) %>% 
    inner_join(meta_df)
  # Plot UMAP
  plots_umap <- lapply(seq_along(var_col), function(i){
    ggplot(umap_df) +
      aes(x = UM1, y = UM2) +
      geom_point(aes_string(color = var_col[i]), size = 2, alpha = 0.8) +
      labs(title = str_c("UMAP: Samples colored by ", date_col), x = "UMAP1", y = "UMAP2") +
      scale_color_manual(values = mycolors) +
      theme_bw() +
      theme(text = element_text(size = 10))
  }) %>% setNames(var_col)
  if(!is.null(date_col)){
    plots_umap$DateTime <- ggplot(umap_df) +
      aes(x = UM1, y = UM2) +
      geom_point(aes_string(color = date_col), size = 2, alpha = 0.8) +
      labs(title = str_c("UMAP: Samples colored by ", date_col), x = "UMAP1", y = "UMAP2") +
      scale_color_viridis_c(
        option = 'G', begin = 0.05, end = 0.95,
        name   = "Sample date",
        breaks  = as.numeric(date_breaks("2 month")(range(umap_df$DateTime))),
        labels  = date_format("%Y-%m-%d")(date_breaks("2 month")(range(umap_df$DateTime)))
      ) +
      theme_bw() +
      theme(text = element_text(size = 10))
  }
  # ggpubr::ggarrange(plotlist = plots_umap)
  
  ret <- list(pca_res = pca_res, pca_df = pca_df, plots_pca = plots_pca,
              tsne_res = tsne_res, tsne_df = tsne_df, plots_tsne = plots_tsne,
              umap_res = umap_res, umap_df = umap_df, plots_umap = plots_umap)
  return(ret)
}


meta_df <- metadata %>%
  mutate(yearmonth = str_c(year, month),
         DateTime = with_tz(ymd_hms(DateTime, tz = 'Asia/Shanghai'), 'UTC')) %>% 
  select(file_id, instrument, trans, batch, year, new, sample_type, DateTime)
id_col <- 'file_id'
date_col <- 'DateTime'
var_col <- colnames(meta_df) %>% setdiff(c(id_col, date_col))
seed <- 10

res_dr_ls <- list()
for(i in seq_along(names(pm_list))[-1]){
  cat(i, '...\n')
  nm <- names(pm_list)[i]
  X <- pm_list[[nm]]
  
  res_dr_ls[[i-1]] <- beca.DR(X, meta_df, id_col, var_col, date_col, seed)
}
names(res_dr_ls) <- names(pm_list)[-1]

graphics.off()
pdf('BECA_diag_dr.pdf', width = 6*3, height = 4*3)
for(nm in names(res_dr_ls)){
  pout <- ggpubr::ggarrange(plotlist = res_dr_ls[[nm]]$plots_pca)
  pout2 <- ggpubr::annotate_figure(pout, top = ggpubr::text_grob(str_c('NA impute method: ', nm), color = "black", face = "bold", size = 14))
  print(pout2)
  
  pout <- ggpubr::ggarrange(plotlist = res_dr_ls[[nm]]$plots_tsne)
  pout2 <- ggpubr::annotate_figure(pout, top = ggpubr::text_grob(str_c('NA impute method: ', nm), color = "black", face = "bold", size = 14))
  print(pout2)
  
  pout <- ggpubr::ggarrange(plotlist = res_dr_ls[[nm]]$plots_tsne)
  pout2 <- ggpubr::annotate_figure(pout, top = ggpubr::text_grob(str_c('NA impute method: ', nm), color = "black", face = "bold", size = 14))
  print(pout2)
}
graphics.off()

tbl_dr <- plyr::ldply(res_dr_ls, function(res){
  res$pca_df %>% inner_join(res$tsne_df) %>% 
    inner_join(res$umap_df)
}, .id = 'naimpute.method')

save(res_dr_ls, file = 'res_dr_ls.RData')



## Hierarchical Clustering ------------
# Create a annotation data frame for samples (to show batch colors on heatmap)
ann_col <- meta_df %>% column_to_rownames('file_id') %>% select(-DateTime)
ann_clrs <- list(
  batch = viridis::viridis_pal(option = 'D', begin = 0.05, end = 0.95)(length(unique(ann_col$batch))) %>% setNames(unique(ann_col$batch)),
  year = viridis::viridis_pal(option = 'A', begin = 0.05, end = 0.95)(length(unique(ann_col$year))) %>% setNames(unique(ann_col$year)),
  sample_type = mycolors[1:length(unique(ann_col$sample_type))] %>% setNames(unique(ann_col$sample_type)),
  instrument = mycolors[(1:length(unique(ann_col$instrument))) + 5] %>% setNames(unique(ann_col$instrument)),
  trans = mycolors[(1:length(unique(ann_col$trans))) + 9] %>% setNames(unique(ann_col$trans)),
  new = c(new = 'orange3', old = 'purple3')
)

# Plot heatmap with hierarchical clustering
res_hc_ls <- list()
for(i in seq_along(names(pm_list))){
  cat(i, '...\n')
  nm <- names(pm_list)[i]
  X <- pm_list[[nm]]
  
  vars <- apply(X, 1, var)
  features <- names(sort(vars, decreasing = TRUE))
  top_features <- features %>%
    head(ifelse(length(features) < 100, length(features), length(features) * 0.5))
  mat_top <- X[top_features, ]
  
  res_hc_ls[[i]] <- pheatmap(
    mat_top, scale = "row",
    clustering_distance_cols = "euclidean",
    clustering_distance_rows = "euclidean",
    clustering_method = "ward.D2",
    annotation_col = ann_col,
    annotation_colors = ann_clrs,
    show_rownames = FALSE, show_colnames = FALSE,
    main = str_glue("Hierarchical clustering (NA impute method: {nm})"),
    filename = str_glue('BECA_diag_hc_{nm}.pdf'), width = 12, height = 5
  )
}
names(res_hc_ls) <- names(pm_list)
save(res_hc_ls, file = 'res_hc_ls.RData')



## Relative Log Expression -----
beca.RLE <- function(X, meta_df, id_col){
  # Calculate Relative Log Expression (RLE) values
  # For each gene j, calculate its median expression across the m samples, i.e. Med(y*j), then calculate the deviations from this median, i.e. calculate yij − Med(y*j), across the is
  # For each sample, generate a boxplot of all the deviations for that sample.
  feature_med <- apply(X, 1, median, na.rm = T)
  rle_matrix <- X - feature_med
  
  df_rle <- rle_matrix %>% t() %>% as.data.frame() %>% 
    rownames_to_column(id_col) %>% 
    inner_join(meta_df, .)
  return(df_rle)
}

res_rle_ls <- list()
for(i in seq_along(names(pm_list))){
  cat(i, '...\n')
  nm <- names(pm_list)[i]
  X <- pm_list[[nm]]
  
  res_rle_ls[[i]] <- beca.RLE(X,
                              meta_df %>%
                                mutate(sample_type = factor(sample_type, c('F', 'T', 'NT', 'N', 'p'))) %>%
                                arrange(sample_type, DateTime),
                              'file_id')
}
names(res_rle_ls) <- names(pm_list)
save(res_rle_ls, file = 'res_rle_ls.RData')


graphics.off()
pdf('BECA_diag_rle.pdf', width = 50, height = 10)
for(nm in names(res_rle_ls)){
  cat(nm, '...\n')
  plot_rle <- res_rle_ls[[nm]] %>%
    select(-(instrument:new), -DateTime) %>%
    pivot_longer(cols = -(file_id:sample_type), names_to = 'Protein', values_to = 'RLE') %>% 
    mutate(file_id = factor(file_id, df_rle$file_id)) %>% 
    ggplot() +
    aes(x = file_id, y = RLE, color = sample_type) +
    geom_boxplot(outlier.size = 0.1, width = 0.1) +
    labs(x = "File ID", y = "RLE", subtitle = str_glue('RLE boxplot: {nm}')) +
    theme_classic() +
    theme(text = element_text(size = 10))
  print(plot_rle)
}
graphics.off()
# ggsave('rle.pdf', plot_rle, width = 50, height = 10, limitsize = F)









ann_col <- df_rle %>% column_to_rownames('file_id') %>% select(1:DateTime)
mat_rle <- df_rle %>% column_to_rownames('file_id') %>%
  select(-(1:DateTime)) %>% 
  t()

# res_rle_ls <- list()
# res_rle_ls[[i]] <- pheatmap(
#   mat_rle, scale = "none",
#   cluster_cols = F, cluster_rows = T,
#   # clustering_distance_cols = "euclidean",
#   clustering_distance_rows = "euclidean",
#   clustering_method = "ward.D2",
#   annotation_col = ann_col %>% select(-DateTime),
#   annotation_colors = ann_clrs,
#   show_rownames = FALSE, show_colnames = FALSE,
#   main = str_glue("Relative Log Expression (NA impute method: {nm})"),
#   filename = str_glue('BECA_diag_rle_{nm}.pdf'), width = 12, height = 5
# )



# plot_rle <- ggplot(rle_df, aes(x = Sample, y = RLE, fill = Batch)) +
#   geom_boxplot(outlier.shape = NA) +
#   labs(title = "Relative Log Expression (RLE) Plot by Sample",
#        y = "Relative Log Expression (log2 scale)",
#        x = "Sample") +
#   theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
# ggsave('rle.pdf', plot_rle, width = 10, height = 10)



## Principal Variance Component Analysis --------
my_pvcaBatchAssess <- function (theDataMatrix, expInfo, threshold) {
  # Cite: Bushel P (2024). pvca: Principal Variance Component Analysis (PVCA). R package version 1.44.0.
  # https://doi.org/10.1002/9780470685983.ch12
  # theDataMatrix, row as probs, col as samples;
  dataRowN <- nrow(theDataMatrix)
  dataColN <- ncol(theDataMatrix)
  theDataMatrixCentered_transposed = apply(theDataMatrix, 1, 
                                           scale, center = TRUE, scale = FALSE)
  theDataMatrixCentered = t(theDataMatrixCentered_transposed)
  theDataCor <- cor(theDataMatrixCentered)
  eigenData <- eigen(theDataCor)
  eigenValues = eigenData$values
  ev_n <- length(eigenValues)
  eigenVectorsMatrix = eigenData$vectors
  eigenValuesSum = sum(eigenValues)
  percents_PCs = eigenValues/eigenValuesSum
  
  exp_design <- as.data.frame(expInfo)
  expDesignRowN <- nrow(exp_design)
  expDesignColN <- ncol(exp_design)
  my_counter_2 = 0
  my_sum_2 = 1
  for (i in ev_n:1) {
    my_sum_2 = my_sum_2 - percents_PCs[i]
    if ((my_sum_2) <= threshold) {
      my_counter_2 = my_counter_2 + 1
    }
  }
  pc_n <- ifelse(my_counter_2 < 3, 3, my_counter_2)
  pc_data_matrix <- matrix(data = 0, nrow = (expDesignRowN * 
                                               pc_n), ncol = 1)
  mycounter = 0
  for (i in 1:pc_n) {
    for (j in 1:expDesignRowN) {
      mycounter <- mycounter + 1
      pc_data_matrix[mycounter, 1] = eigenVectorsMatrix[j, 
                                                        i]
    }
  }
  AAA <- exp_design[rep(1:expDesignRowN, pc_n), ]
  Data <- cbind(AAA, pc_data_matrix)
  variables <- c(colnames(exp_design))
  for (i in 1:length(variables)) {
    Data$variables[i] <- as.factor(Data$variables[i])
  }
  op <- options(warn = (-1))
  on.exit(options(op))
  effects_n = expDesignColN + choose(expDesignColN, 2) + 1
  randomEffectsMatrix <- matrix(data = 0, nrow = pc_n, ncol = effects_n)
  model.func <- c()
  index <- 1
  for (i in 1:length(variables)) {
    mod = paste("(1|", variables[i], ")", sep = "")
    model.func[index] = mod
    index = index + 1
  }
  for (i in 1:(length(variables) - 1)) {
    for (j in (i + 1):length(variables)) {
      mod = paste("(1|", variables[i], ":", variables[j], 
                  ")", sep = "")
      model.func[index] = mod
      index = index + 1
    }
  }
  function.mods <- paste(model.func, collapse = " + ")
  for (i in 1:pc_n) {
    y = (((i - 1) * expDesignRowN) + 1)
    funct <- paste("pc_data_matrix", function.mods, sep = " ~ ")
    Rm1ML <- lme4::lmer(funct, Data[y:(((i - 1) * expDesignRowN) + 
                                         expDesignRowN), ], REML = TRUE, verbose = FALSE, 
                        na.action = na.omit)
    randomEffects <- Rm1ML
    randomEffectsMatrix[i, ] <- c(unlist(lme4::VarCorr(Rm1ML)), 
                                  resid = sigma(Rm1ML)^2)
  }
  effectsNames <- c(names(lme4::getME(Rm1ML, "cnms")), "resid")
  randomEffectsMatrixStdze <- matrix(data = 0, nrow = pc_n, 
                                     ncol = effects_n)
  for (i in 1:pc_n) {
    mySum = sum(randomEffectsMatrix[i, ])
    for (j in 1:effects_n) {
      randomEffectsMatrixStdze[i, j] = randomEffectsMatrix[i, j]/mySum
    }
  }
  randomEffectsMatrixWtProp <- matrix(data = 0, nrow = pc_n, 
                                      ncol = effects_n)
  for (i in 1:pc_n) {
    weight = eigenValues[i]/eigenValuesSum
    for (j in 1:effects_n) {
      randomEffectsMatrixWtProp[i, j] = randomEffectsMatrixStdze[i, j] * weight
    }
  }
  randomEffectsSums <- matrix(data = 0, nrow = 1, ncol = effects_n)
  randomEffectsSums <- colSums(randomEffectsMatrixWtProp)
  totalSum = sum(randomEffectsSums)
  randomEffectsMatrixWtAveProp <- matrix(data = 0, nrow = 1, 
                                         ncol = effects_n)
  for (j in 1:effects_n) {
    randomEffectsMatrixWtAveProp[j] = randomEffectsSums[j]/totalSum
  }
  return(list(dat = randomEffectsMatrixWtAveProp, label = effectsNames, matrixWtProp = set_colnames(randomEffectsMatrixWtProp, effectsNames), eigenData = eigenData))
}
# pca_colors <- c('#D7DA86', '#330A5FFF', '#FCB519FF', '#ED6925FF', '#781C6DFF', '#BB3754FF', '#000004FF') %>% setNames(c('resid', 'instrument', 'trans', 'Patient', 'instrument:trans'))

res_pvca_ls <- list()
for(i in seq_along(names(pm_list))[-1]){
  cat(i, '...\n')
  nm <- names(pm_list)[i]
  X <- pm_list[[nm]]
  theDataMatrix <- t(X)
  expInfo <- metadata %>%
    column_to_rownames('file_id') %>%
    select(instrument, trans, batch, sample_type, new, anatomical_classification)
  threshold <- 0.9
  pvca <- my_pvcaBatchAssess(theDataMatrix, expInfo, threshold)

  res_pvca_ls[[i-1]] <- pvca
}
names(res_pvca_ls) <- names(pm_list)[-1]

save(res_pvca_ls, file = 'res_pvca_ls.RData')



pvca <- res_pvca_ls$minimum


# Average proportion
df_pvca <- data.frame(t(pvca$dat)) %>% 
  cbind(pvca$label) %>% 
  setNames(c('RandomEffectWtAveProp', 'EffectName')) %>% 
  arrange(RandomEffectWtAveProp) %>% 
  mutate(EffectName = factor(EffectName, levels = EffectName)) %>% 
  arrange(desc(RandomEffectWtAveProp)) %>% 
  mutate(y_label = cumsum(RandomEffectWtAveProp) - 0.5 * RandomEffectWtAveProp)
plot_pvca <- ggplot(df_pvca, aes(x = 2, y = RandomEffectWtAveProp, fill = EffectName)) +
  geom_bar(stat = 'identity', color = 'white') +
  coord_polar(theta = 'y', start = 0) +
  theme(legend.position = 'none') +
  geom_text(aes(y = y_label, label = round(RandomEffectWtAveProp*100, 2)), color = 'white', size = 6) +
  scale_fill_manual(values = pca_colors) +
  theme_void() +
  xlim(0.5, 2.5)

# every PC
df_prop <- data.frame(pvca$matrixWtProp)
df_prop_percent <- data.frame(t(apply(df_prop, 1, function(x) x / sum(x))))
df_prop_percent$PC <- str_glue('PC{1:nrow(df_prop_percent)}: ({round(100 * pvca$eigenData$values[1:nrow(df_prop_percent)] / sum(pvca$eigenData$values[1:nrow(df_prop_percent)]), 2)} %)')
df_prop_percent$PC %<>% factor(., levels = .)

plot_pvca_individual <- plyr::dlply(df_prop_percent, 'PC', function(dfsub){
  tbl <- dfsub %>% select(-PC) %>% 
    setNames(., str_replace(names(.), '\\.', ':')) %>% 
    t() %>% data.frame() %>% 
    setNames('RandomEffectWtAveProp') %>% 
    rownames_to_column('EffectName') %>% 
    arrange(RandomEffectWtAveProp) %>% 
    mutate(EffectName = factor(EffectName, levels = EffectName)) %>% 
    arrange(desc(RandomEffectWtAveProp)) %>% 
    mutate(y_label = cumsum(RandomEffectWtAveProp) - 0.5 * RandomEffectWtAveProp)
  
  set.seed(1000)
  ggplot(tbl, aes(x = 2, y = RandomEffectWtAveProp, fill = EffectName)) +
    geom_bar(stat = 'identity', color = 'white') +
    coord_polar(theta = 'y', start = 0) +
    geom_text(aes(y = y_label, label = round(RandomEffectWtAveProp*100, 2)), color = 'white', size = 6) +
    scale_fill_manual(values = pca_colors) +
    labs(subtitle = dfsub$PC) +
    theme_void() +
    theme(legend.position = 'none') +
    xlim(0.5, 2.5)
})
length(plot_pvca_individual) # 771
legend_pies <- ggpubr::get_legend(plot_pvca)
p <- ggpubr::ggarrange(plotlist = plot_pvca_individual, nrow = 20, ncol = 40) %>%
  ggpubr::ggarrange(legend_pies, widths = c(0.95, 0.05))
plot(cumsum(pvca$eigenData$values / sum(pvca$eigenData$values)), xlab = 'PC_n', ylab = 'Summed eigenvalues proportion')
ggsave('BECA_diag_pvca_individual_minimum.pdf', p, width = 4*40, height = 4*20, limitsize = F)


# top3 PCs
pvca$eigenData$values[1:3]
dfsub <- apply(df_prop_percent[1:3, 1:22], 2, function(y) {
  y * pvca$eigenData$values[1:3]
}) %>% as.data.frame() %>%
  mutate_all(sum) %>% slice(1)
dfsub <- dfsub / sum(dfsub)
PVCA_top3 <- dfsub %>%
  setNames(., str_replace(names(.), '\\.', ':')) %>% 
  t() %>% data.frame() %>% 
  setNames('RandomEffectWtAveProp') %>% 
  rownames_to_column('EffectName') %>% 
  arrange(RandomEffectWtAveProp) %>% 
  mutate(EffectName = factor(EffectName, levels = EffectName)) %>% 
  arrange(desc(RandomEffectWtAveProp)) %>% 
  mutate(y_label = cumsum(RandomEffectWtAveProp) - 0.5 * RandomEffectWtAveProp)


plot_PVCA_top3 <- ggplot(PVCA_top3, aes(x = 2, y = RandomEffectWtAveProp, fill = EffectName)) +
  geom_bar(stat = 'identity', color = 'white') +
  coord_polar(theta = 'y', start = 0) +
  geom_text(aes(y = y_label, label = round(RandomEffectWtAveProp*100, 2)), color = 'white', size = 6) +
  # geom_text(aes(x = 1, y = y_label, label = round(RandomEffectWtAveProp*100, 2)), color = 'black', size = 6) +
  scale_fill_manual(values = pca_colors) +
  labs(subtitle = dfsub$PC) +
  theme_void() +
  xlim(0.5, 2.5)

ggsave('BECA_diag_pvca_top3PCs_minimum.pdf',plot_PVCA_top3, width = 7, height = 4)




# 2.Batch effect correction ---------








# 3.Assessment of BECA --------



# Try proBatch ------------
if (!requireNamespace("proBatch", quietly = TRUE)) {
  devtools::install_github("symbioticMe/proBatch", build_vignettes = TRUE)
}
library(proBatch)





# SPLIT -------
# Prepare the data and metadata ------------
data <- df %>% column_to_rownames('V1')

metadata <- info %>%
  column_to_rownames('FileName') %>% 
  select(sample_type, tissue_name, anatomical_location, anatomical_classification, patient_ID)

valid_ids <- intersect(rownames(data), rownames(metadata)) %>% sort()
data_fit <- data[valid_ids, , drop = F]
metadata_fit <- metadata[valid_ids, , drop = F]

dim(metadata_fit)  # 1780 files × 5 fields
dim(data_fit)      # 1780 files × 13794 proteins
all(rownames(data_fit) == rownames(metadata_fit)) # row names match exactly


# library(imputeLCMD)     # for impute.MinProb()
# library(preprocessCore) # for normalize.quantiles()

# Log2‑transform
matpg <- t(data_fit)
matpg[matpg < 1] <- 1
matpg_log2 <- log2(matpg)

# keep proteins quantified in ≥70% of samples
min_samples <- ceiling(0.5 * ncol(matpg_log2))
keep_prot   <- rowSums(!is.na(matpg_log2)) >= min_samples
matpg_filt  <- matpg_log2[keep_prot, , drop = FALSE]
cat("Kept", nrow(matpg_filt), "of", nrow(matpg_log2), "proteins.\n")



# test naimpute
set.seed(1)
x <- rnorm(10000, 20, 3)
set.seed(1)
x[sample(1:length(x), 1000, replace = F)] <- NA
x <- matrix(x, nrow = 50, ncol = 200)
xx <- lapply(seq_along(namethods) %>% setdiff(c(1:21)), function(i){
  print(i)
  naimpute(x, namethods[i])
} )

matpg_imp <- impute.MinProb(as.matrix(matpg_filt), q = 0.01) # left‐shifted Gaussian (q=0.01: bottom 1%)



# ─── 5. Quantile‑normalize across samples ────────────────────────────────────
# forces each sample (column) to have identical distribution
matpg_norm <- normalize.quantiles(matpg_imp)
rownames(matpg_norm) <- rownames(matpg_imp)
colnames(matpg_norm) <- colnames(matpg_imp)

# ─── 6. Result ───────────────────────────────────────────────────────────────
# 'matpg_norm' is your preprocessed intensity matrix (proteins × samples)
matpg_preprocessed <- matpg_norm

# Sanity checks
stopifnot(!any(is.na(matpg_preprocessed)))
summary(as.vector(matpg_preprocessed))

# # quick sanity checks
# summary(data_preprocessed)
# any(is.na(data_preprocessed))  # should be FALSE






# now proceed with your downstream analysis, e.g. t‑SNE:
# set.seed(42)
# tsne_out <- Rtsne(as.matrix(data_preprocessed), perplexity = 30)




# 8. Run t‑SNE on the protein abundance matrix
set.seed(42)  # for reproducibility
tsne_out <- Rtsne::Rtsne(
  as.matrix(data_fit),
  dims = 2, initial_dims = 50, perplexity = 30, theta = 0.5,
  check_duplicates = TRUE, pca = TRUE, partial_pca = FALSE,
  max_iter = 1000, is_distance = FALSE, pca_center = TRUE, pca_scale = FALSE,
  normalize = TRUE, momentum = 0.5, final_momentum = 0.8, eta = 200,
  exaggeration_factor = 12, num_threads = 1
)

# 9. Combine the t‑SNE coordinates with your metadata
tsne_df <- as.data.frame(tsne_out$Y)
colnames(tsne_df) <- c("TSNE1", "TSNE2")
tsne_df <- cbind(tsne_df, metadata_fit)

# 10. Plot in ggplot2, colouring points by sample_type
ggplot(tsne_df, aes(x = TSNE1, y = TSNE2, color = sample_type)) +
  geom_point(size = 2, alpha = 0.8) +
  labs(title = "t-SNE of Protein Matrix",
       x     = "Dimension 1",
       y     = "Dimension 2",
       color = "Sample Type") +
  theme_minimal()





# 1. INITIAL QUALITY ASSESSMENT ---------------

## 1.1 Dimensions of data & metadata
cat(sprintf("Data matrix: %d runs × %d proteins\n", nrow(data_fit), ncol(data_fit)))
cat(sprintf("Metadata table: %d runs × %d fields\n\n", nrow(metadata_fit), ncol(metadata_fit)))

## 1.2 Missing‑value diagnostics
missing_by_run     <- rowSums(is.na(data_fit))
missing_by_protein <- colSums(is.na(data_fit))
cat("Missing values per run:\n")
cat(sprintf("  • min = %d, max = %d, total = %d\n",
            min(missing_by_run), max(missing_by_run), sum(missing_by_run)))
cat("Missing values per protein:\n")
cat(sprintf("  • min = %d, max = %d, total = %d\n\n",
            min(missing_by_protein), max(missing_by_protein), sum(missing_by_protein)))

# Heatmap of missing values (1 = missing, 0 = observed)
missing_mat <- is.na(data_fit) * 1
pheatmap(missing_mat,
         cluster_rows = FALSE, cluster_cols = FALSE,
         show_rownames = FALSE, show_colnames = FALSE,
         main = "Missing‑Value Heatmap")

## 1.3 Outlier detection via Z‑scores (|Z| > 3 per run)
z_scores <- t(apply(data_fit, 1, function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}))
n_outliers <- sum(abs(z_scores) > 3, na.rm = TRUE)
cat(sprintf("Total outlier measurements (|Z| > 3): %d\n\n", n_outliers))

## 1.4 Descriptive statistics per protein
desc_stats <- data.frame(
  Protein = colnames(data_fit),
  Min     = apply(data_fit, 2, min, na.rm = TRUE),
  `1st_Qu`= apply(data_fit, 2, quantile, probs = 0.25, na.rm = TRUE),
  Median  = apply(data_fit, 2, median, na.rm = TRUE),
  `3rd_Qu`= apply(data_fit, 2, quantile, probs = 0.75, na.rm = TRUE),
  Max     = apply(data_fit, 2, max, na.rm = TRUE)
)
print(head(desc_stats))  # show first few proteins’ stats

## 1.5 Overall intensity distribution
all_vals <- as.numeric(as.matrix(data_fit))
ggplot(data.frame(Intensity = all_vals), aes(x = Intensity)) +
  geom_histogram(bins = 100, aes(y = ..density..), color = "black", fill = "grey80") +
  geom_density() +
  ggtitle("Distribution of All Protein Intensities") +
  theme_minimal()

## 1.6 Sample–sample correlation & clustering
corr_mat <- cor(t(data_fit), use = "pairwise.complete.obs")
pheatmap(corr_mat,
         main = "Sample–Sample Correlation Clustermap",
         show_rownames = FALSE, show_colnames = FALSE)

# 2. NORMALIZATION PIPELINE ---------------

## 2.1 Log₂(x + 1) transformation
data_log2 <- log2(data_fit + 1)
cat("Applied log2(x + 1) transformation.\n")

## 2.2 Quantile normalization to a Gaussian distribution
map_to_normal <- function(x) {
  rk      <- rank(x, na.last = "keep", ties.method = "average")
  uniform <- (rk - 0.5) / sum(!is.na(x))
  qnorm(uniform)
}
data_qn <- apply(data_log2, 2, map_to_normal)
cat("Applied quantile normalization to a normal distribution.\n")

## 2.3 Standard (Z‑score) scaling across proteins
data_scaled <- scale(data_qn, center = TRUE, scale = TRUE)
cat("Applied standard (Z‑score) scaling across proteins.\n")

## 2.4 Unit‑norm normalization per run
row_norms <- apply(data_scaled, 1, function(x) sqrt(sum(x^2, na.rm = TRUE)))
data_norm <- data_scaled / row_norms
cat("Applied unit‑norm normalization per run.\n\n")

## 2.5 Save normalized data
write.csv(data_norm,
          file = "normalized_proteomics_data.csv",
          quote = FALSE,
          row.names = TRUE)
cat("Normalization complete. Saved normalized data to 'normalized_proteomics_data.csv'.\n")








# OUTPUT ----------
list(rep.identity = rep.identity,
     rep.quantity = rep.quantity,
     rep.pearson = rep.pearson,
     est.identity = est.identity,
     est.pearson = est.pearson,
     BECA.DR = tbl_dr) %>%
  rio::export('source_data.xlsx')

# save.image('BEAC_v20250729.RData')

