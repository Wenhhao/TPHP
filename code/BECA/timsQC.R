rm(list = ls())
pacman::p_unload(pacman::p_loaded(), character.only = T)

rstudioapi::getActiveDocumentContext()$path
library(tidyverse)
library(magrittr)
library(ggpubr)
library(ggsci)
source('../source/my_fun.R')


# 1.read data ----------
info_datetime <- rio::import('rawdata_datetime.xlsx')
info1 <- rio::import('20250725_PUH_sample_information_1781files_info_edited_v10.xlsx')
dat1 <- read.delim('mapped_pg_matrix_1780_14062.csv', sep = ',', header = T, row.names = 1, check.names = F, stringsAsFactors = F)
pool1 <- read.delim('pool_pg_matrix_120_14063.csv', sep = ',', header = T, row.names = 1, check.names = F, stringsAsFactors = F)

# log2 transform
dat1 <- rbind(dat1, pool1) %>%
  t() %>% log2()
dim(dat1) # 14062  1900

# append labels
info1_full <- info1 %>% filter(!is.na(anatomical_classification)) %>% 
  full_join(info_datetime)
info1 <- info1_full %>% filter(FileName %in% colnames(dat1))
setequal(colnames(dat1), info1$FileName) # TRUE
info1 <- info1 %>% set_rownames(info1$FileName) %>% .[colnames(dat1), ]
identical(colnames(dat1), info1$FileName) # TRUE

info1_fit <- info_of_pm(dat1) %>% rename(FileName = Filename) %>% 
  inner_join(info1_full, .) %>%
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
rep_smp <- info1_fit %>% filter(Is.Rep) %>% pull(sample_id) %>% unique()
rep_smp_wrong <- info1_fit %>% filter(sample_id %in% rep_smp) %>% count(sample_id) %>% filter(n == 1) %>% pull(sample_id)
info1_fit %<>% mutate(
  Rep = ifelse(sample_id %in% rep_smp_wrong, NA, Rep),
  Is.Rep = ifelse(sample_id %in% rep_smp_wrong, F, Is.Rep),
  Rep_type = ifelse(sample_id %in% rep_smp_wrong, NA, Rep_type)
)

info1_fit1 <- info1_fit %>% filter(Is.Rep) %>%
  group_by(sample_id, Rep_type) %>%
  arrange(DateTime) %>%
  mutate(
    seq_num  = row_number(),
    file_id  = str_c(sample_id, '_', Rep_type, seq_num),
    seq_num = NULL, .after = sample_id
  ) %>%
  ungroup()
info1_fit2 <- info1_fit %>% filter(!Is.Rep) %>% mutate(file_id = sample_id)

info1_fitc <- info1_fit1 %>% rbind(info1_fit2) %>% as.data.frame() %>%
  set_rownames(.$FileName) %>% .[colnames(dat1), ] %>% 
  remove_rownames() %>%
  mutate(file_id = ifelse(is.na(file_id), str_c('pool', str_extract(FileName, '^[A-Z]+'), str_extract(FileName, '\\d+$')), file_id),
         sample_type = ifelse(is.na(sample_type), 'p', sample_type),
         sample_id = ifelse(!is.na(sample_id), sample_id, sample_type)) # for pooling samples

identical(info1$FileName, info1_fitc$FileName) # TRUE

# 
# meta1 <- info1 %>%
#   select(file_id, BatchID, Batch_m, Batch_n, Date,
#          Is.trep, Is.brep, Is.pool, Rep_type,
#          Inner_ring_label, `癌症类型`, `性别`, `年龄`) %>% 
#   mutate(Rep_type = ifelse(Is.pool, 'Pooling', Rep_type)) %>% 
#   mutate(sample_id = ifelse(!is.na(BatchID), BatchID, Rep_type), .after = BatchID)
# 
# 
rio::export(info1_fitc, '20250725_PUH_sample_information_1900files_v10_addLabels.xlsx')
rm(list = str_subset(ls(), '^info1'))
info1 <- rio::import('20250725_PUH_sample_information_1900files_v10_addLabels.xlsx')

meta1 <- info1 %>% select(FileName, file_id, sample_id, sample_type, Inner_ring_label, tissue_name, date)



# 2.check -----
df_ident <- data.frame(FileName = colnames(dat1),
                       `# proteins` = colSums(!is.na(dat1)),
                       check.names = F)


tmp_check <- info1 %>% count(file_id) %>% filter(n > 1) %>% semi_join(info1, .) %>% 
  inner_join(df_ident) %>% 
  arrange(file_id, date) %>% 
  mutate(file_id = factor(file_id, unique(file_id))) %>% 
  group_by(file_id) %>% mutate(Group = c('Earlier', 'Later'))
rio::export(tmp_check, 'PUH_double_fileID_check.xlsx')
p.dfcheck <- tmp_check %>%
  ggplot() +
  aes(x = Group, y = `# proteins`, fill = Group) +
  geom_boxplot(alpha = 0.5) +
  geom_line(aes(group = file_id), 
            color = "grey70", linetype = "dashed", size = 0.5) +
  geom_point(aes(color = file_id), fill = NA, size = 3) +
  labs(x = '', y = '# proteins') +
  scale_color_manual(values = unique(c(pal_nejm()(8), pal_npg()(10), pal_lancet()(9), mycolors))) +
  scale_fill_bmj() +
  stat_compare_means(method = "t.test", paired = TRUE) +
  theme_classic2()
ggsave('PUH_double_fileID_check.pdf', p.dfcheck, width = 8, height = 6)


# 3.replicates check --------
df_est_rep <- dat1 %>% t() %>% as.data.frame() %>%
  rownames_to_column('FileName') %>% 
  inner_join(info1 %>% select(FileName, file_id), .) %>% 
  select(-FileName) %>% 
  filter(!(file_id %in% tmp_check$file_id)) %>%
  inner_join(info1 %>% filter(!(file_id %in% tmp_check$file_id)) %>%
               count(sample_id) %>% filter(n > 1) %>%
               semi_join(info1, .) %>% arrange(sample_id), .)

rm.fid <- c()


rep.identity <- plyr::ddply(df_est_rep, 'sample_id', function(dfsub){
  cat(dfsub$sample_id[1], '...\r')
  X <- dfsub %>%
    column_to_rownames('file_id') %>%
    select(-(1:date.in.filename)) %>% t()
  ret <- data.frame(file_id = dfsub$file_id,
                    `# proteins` = colSums(!is.na(X)),
                    check.names = F)
  ret$outlier.lower.ingroup <- get_outliers(ret$`# proteins`)[1]
  return(ret)
})
rep.identity %<>%
  group_by(sample_id) %>%
  summarise(identity.mean = mean(`# proteins`),
            identity.min = min(`# proteins`),
            identity.max = max(`# proteins`),
            identity.range.length = identity.max - identity.min) %>%
  inner_join(rep.identity, .) %>%
  arrange(identity.mean) %>%
  mutate(sample_id = factor(sample_id, levels = unique(sample_id)),
         Is.Lower.Ingroup = `# proteins` < outlier.lower.ingroup) %>% 
  inner_join(info1 %>% select(file_id, Rep_type, Inner_ring_label, sample_type))
rep.sid_order <- rep.identity %>% arrange(`# proteins`) %>% distinct(sample_id) %>% pull()
rep.identity %<>% 
  mutate(sample_id = factor(sample_id, rep.sid_order),
         Rep_type = ifelse(sample_type == 'p', 'p', Rep_type),
         Rep_type = ifelse(!is.na(Rep_type), Rep_type, 'Reference'),
         Rep_type = factor(Rep_type, c('Reference', 'b', 't', 'p')))

low.rep.ident <- rep.identity %>% filter(identity.range.length > 1000) %>%
  group_by(sample_id) %>% arrange(`# proteins`) %>% slice(1) %>% 
  pull(file_id)
rm.fid %<>% append(low.rep.ident)
rm.sid <- df_est_rep %>% filter(file_id %in% rm.fid) %>% distinct(sample_id) %>% pull() %>% setdiff(c('p'))


rep.pearson <- plyr::ddply(df_est_rep %>%
                             filter(!(sample_id %in% rm.sid),
                                    !(file_id %in% rm.fid)),
                           'sample_id', function(dfsub){
  cat(dfsub$sample_id[1], '...\r')
  X <- dfsub %>%
    column_to_rownames('file_id') %>%
    select(-(1:date.in.filename)) %>% t()
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


# 4.samples -----
df_est <- dat1 %>% t() %>% as.data.frame() %>%
  rownames_to_column('FileName') %>% 
  inner_join(meta1, .)


est.identity <- df_est %>%
  filter(sample_id != 'p') %>% 
  filter(!(file_id %in% rm.fid)) %>% 
  plyr::ddply(c('sample_type', 'Inner_ring_label'), function(dfsub){
    cat(dfsub$sample_type[1], dfsub$Inner_ring_label[1], '...\r')
    X <- dfsub %>%
      column_to_rownames('FileName') %>%
      select(-(1:date)) %>% t()
    ret <- data.frame(FileName = dfsub$FileName,
                      `# proteins` = colSums(!is.na(X)),
                      check.names = F)
    ret$outlier.lower.ingroup <- get_outliers(ret$`# proteins`)[1]
    return(ret)
  })
est.identity %<>%
  group_by(sample_type, Inner_ring_label) %>%
  summarise(identity.mean = mean(`# proteins`),
            identity.min = min(`# proteins`),
            identity.max = max(`# proteins`),
            identity.range.length = identity.max - identity.min,
            .groups = 'drop') %>%
  inner_join(est.identity, .) %>%
  mutate(Group = str_c(sample_type, ' - ', Inner_ring_label)) %>% 
  arrange(identity.mean) %>%
  mutate(Group = factor(Group, levels = unique(Group)),
         Is.Lower.Ingroup = `# proteins` < outlier.lower.ingroup) %>% 
  inner_join(info1 %>% select(FileName, Rep_type, Inner_ring_label, sample_type, tissue_name))
group_order <- est.identity %>% arrange(`# proteins`) %>% distinct(Group) %>% pull()
est.identity %<>% mutate(Group = factor(Group, group_order))




# Output -------
plot_ident_rep <- ggplot(rep.identity) +
  aes(x = `# proteins`, y = sample_id) +
  geom_boxplot(data = rep.identity %>% filter(identity.range.length <= 1000),
               color = 'black', outlier.color = '#c23190', outlier.size = 3) +
  geom_boxplot(data = rep.identity %>% filter(identity.range.length > 1000),
               color = 'red4', outlier.color = '#c23190', outlier.size = 3) +
  geom_point(aes(color = Rep_type), size = 1.2) +
  labs(x = "# proteins", y = "Sample ID", subtitle = "Protein identification") +
  ggsci::scale_color_aaas(name = 'Replicates type') +
  theme_bw() +
  theme(text = element_text(size = 10, color = 'black'))
plot_pearson_rep <- ggplot(rep.pearson) +
  aes(x = pearson.r, y = sample_id) +
  # geom_boxplot(color = '#000000') +
  geom_point(color = '#000000', size = 1.2) +
  labs(x = "Pearson's r", y = "Sample ID", subtitle = "Correlation - same.sample.id") +
  theme_bw() +
  theme(text = element_text(size = 10, color = 'black'))
plot_rep <- ggarrange(plot_ident_rep, plot_pearson_rep, nrow = 1, ncol = 2, widths = c(4, 3))
ggsave('tims_rep.pdf', plot_rep, width = 10, height = 8)


plot_ident <- ggplot(est.identity) +
  aes(x = `# proteins`, y = Group) +
  geom_boxplot(data = est.identity %>% filter(identity.range.length <= 1000),
               color = 'black', outlier.color = '#c23190', outlier.size = 3) +
  geom_boxplot(data = est.identity %>% filter(identity.range.length > 1000),
               color = 'red4', outlier.color = '#c23190', outlier.size = 3) +
  geom_point(aes(color = tissue_name), size = 1.2) +
  labs(x = "# proteins", y = "Sample ID", subtitle = "Protein identification") +
  # ggsci::scale_color_aaas(name = 'Sample type') +
  theme_bw() +
  theme(text = element_text(size = 10, color = 'black'))
ggsave('tims_identity.pdf', plot_ident, width = 40, height = 15)


list(rep.identity = rep.identity,
     rep.pearson = rep.pearson,
     sample.identity = est.identity) %>% 
  rio::export('tims_QC_source.xlsx')






