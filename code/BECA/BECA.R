rstudioapi::getActiveDocumentContext()$path
setwd('~/GitHub/TPHP/code/BECA/')
source('../../source/my_fun.R')
library(tidyverse)
library(magrittr)
library(RColorBrewer)
library(pheatmap)
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
  full_join(info, .) %>%
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
rep_sup <- info_fit %>% filter(!Is.Rep) %>%
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
  remove_rownames()

rio::export(info_fitc, '~/GitHub/TPHP/input/20250725_1780files_add_fields.xlsx')

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




df_est <- matpg %>% t() %>% as.data.frame() %>% 
  rownames_to_column('file_id') %>% 
  inner_join(metadata %>% select(file_id:anatomical_classification), .)

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
pool.identity <- df_est %>% filter(!is.na(Rep_type)) %>%
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
pool.identity %<>%
  group_by(sample_id) %>%
  summarise(identity.mean = mean(`# proteins`)) %>%
  inner_join(rep.identity, .) %>%
  arrange(identity.mean) %>%
  mutate(sample_id = factor(sample_id, levels = unique(sample_id)),
         Is.Lower.Ingroup = `# proteins` < outlier.lower.ingroup) %>% 
  inner_join(info_fitc %>% select(file_id, Rep_type))





### samples check -----
### identity
est.identity <- plyr::ddply(df_est, c('sample_type', 'tissue_name'), function(dfsub){
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
est.pearson <- df_est %>% filter(!(file_id %in% rm.fid)) %>%
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

est.pearson_log <- df_est %>% filter(!(file_id %in% rm.fid)) %>%
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


est.pearson_qlog <- df_est_qlog %>% filter(!(file_id %in% rm.fid)) %>%
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
  count(n) # 12 patients

info_fitc %>%
  filter(sample_type %in% c('T', 'NT'),
         !(file_id %in% rm.fid)) %>%
  count(patient_ID) %>% 
  count(n) # 18 patients


# 1.Diagnosis of batch effect --------








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
     est.pearson = est.pearson) %>%
  rio::export('source_data.xlsx')

save.image('preprocess.RData')
