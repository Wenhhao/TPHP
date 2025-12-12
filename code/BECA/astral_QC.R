rstudioapi::getActiveDocumentContext()$path
library(tidyverse)
library(magrittr)
library(ggpubr)
library(ggsci)
source('../source/my_fun.R')

# 1.read data ----------
info2 <- rio::import('20250807_PUH_RCA_sample_and_pool_1227files_info.xlsx')
# dat2 <- rio::import('mapped_dat2_pg_matrix_1200_12501.csv')
# pool2 <- rio::import('mapped_dat2_pool_matrix_27_12501.csv')
dat2 <- read.delim('mapped_dat2_pg_matrix_1200_12501.csv', sep = ',', header = T, row.names = 1, check.names = F, stringsAsFactors = F)
pool2 <- read.delim('mapped_dat2_pool_matrix_27_12501.csv', sep = ',', header = T, row.names = 1, check.names = F, stringsAsFactors = F)

cn_cn <- c('胃腺癌','大肠粘液腺癌','乳腺癌','肝癌','肾透明细胞癌','结肠腺癌','宫颈癌','胰腺癌','舌癌','前列腺癌','食管癌','肺腺癌','子宫内膜癌','阴茎癌','胆总管癌','卵巢高级别浆液性癌','卵巢透明细胞样癌','直肠癌','胸腺瘤','胚胎癌','精原细胞瘤','肝内胆管癌','输卵管癌','胃肠道间质瘤','甲状腺','喉癌','平滑肌肉瘤','胶质母细胞瘤','弥漫大B','反应性增生淋巴结','胆囊癌','横纹肌肉瘤','膀胱癌')
cancer_cn2en <- c("胃腺癌" = "gastric adenocarcinoma", "大肠粘液腺癌" = "colonic mucinous adenocarcinoma", "乳腺癌" = "breast cancer", "肝癌" = "liver cancer", "肾透明细胞癌" = "clear cell renal cell carcinoma", "结肠腺癌" = "colon adenocarcinoma", "宫颈癌" = "cervical cancer", "胰腺癌" = "pancreatic cancer", "舌癌" = "tongue cancer", "前列腺癌" = "prostate cancer", "食管癌" = "esophageal cancer", "肺腺癌" = "lung adenocarcinoma", "子宫内膜癌" = "endometrial cancer", "阴茎癌" = "penile cancer", "胆总管癌" = "common bile duct carcinoma", "卵巢高级别浆液性癌" = "high-grade serous ovarian carcinoma", "卵巢透明细胞样癌" = "ovarian clear cell carcinoma", "直肠癌" = "rectal cancer", "胸腺瘤" = "thymoma", "胚胎癌" = "embryonal carcinoma", "精原细胞瘤" = "seminoma", "肝内胆管癌" = "intrahepatic cholangiocarcinoma", "输卵管癌" = "fallopian tube carcinoma", "胃肠道间质瘤" = "gastrointestinal stromal tumor", "甲状腺" = "thyroid cancer", "喉癌" = "laryngeal cancer", "平滑肌肉瘤" = "leiomyosarcoma", "胶质母细胞瘤" = "glioblastoma", "弥漫大B" = "diffuse large B-cell lymphoma", "反应性增生淋巴结" = "reactive lymph node hyperplasia", "胆囊癌" = "gallbladder cancer", "横纹肌肉瘤" = "rhabdomyosarcoma", "膀胱癌" = "bladder cancer")

# append labels
meta2 <- info2 %>%
  select(file_id, BatchID, Batch_m, Batch_n, Date,
         Is.trep, Is.brep, Is.pool, Rep_type,
         Cancer_type, `癌症类型`, `性别`, `年龄`) %>% 
  mutate(Rep_type = ifelse(Is.pool, 'Pooling', Rep_type)) %>% 
  mutate(SampleID = ifelse(!is.na(BatchID), BatchID, Rep_type), .after = BatchID) %>% 
  mutate(Cancer_name = cancer_cn2en[癌症类型])

# log2 transform
df2 <- rbind(dat2, pool2) %>%
  t() %>% log2()


# 2.replicates ----
df_est_rep <- df2 %>% t() %>% as.data.frame() %>%
  rownames_to_column('FileName') %>% 
  inner_join(info2 %>% select(FileName, file_id), .) %>% 
  select(-FileName) %>% 
  inner_join(meta2 %>% count(SampleID) %>% filter(n > 1) %>%
              semi_join(meta2, .) %>% arrange(SampleID), .)
rm.fid <- c()


rep.identity <- plyr::ddply(df_est_rep, 'SampleID', function(dfsub){
    cat(dfsub$SampleID[1], '...\r')
    X <- dfsub %>%
      column_to_rownames('file_id') %>%
      select(-(BatchID:年龄)) %>% t()
    ret <- data.frame(file_id = dfsub$file_id,
                      `# proteins` = colSums(!is.na(X)),
                      check.names = F)
    ret$outlier.lower.ingroup <- get_outliers(ret$`# proteins`)[1]
    return(ret)
  })
rep.identity %<>%
  group_by(SampleID) %>%
  summarise(identity.mean = mean(`# proteins`),
            identity.min = min(`# proteins`),
            identity.max = max(`# proteins`),
            identity.range.length = identity.max - identity.min) %>%
  inner_join(rep.identity, .) %>%
  arrange(identity.mean) %>%
  mutate(SampleID = factor(SampleID, levels = unique(SampleID)),
         Is.Lower.Ingroup = `# proteins` < outlier.lower.ingroup) %>% 
  inner_join(meta2 %>% select(file_id, Rep_type))
sid_order <- rep.identity %>% arrange(`# proteins`) %>% distinct(SampleID) %>% pull()
rep.identity %<>% 
  mutate(SampleID = factor(SampleID, sid_order),
         Rep_type = ifelse(!is.na(Rep_type), Rep_type, 'Reference'),
         Rep_type = factor(Rep_type, c('Reference', 'biological replicates', 'technical replicates', 'Pooling')))

low.rep.ident <- rep.identity %>% filter(identity.range.length > 1000) %>%
  group_by(SampleID) %>% arrange(`# proteins`) %>% slice(1) %>% 
  pull(file_id)
rm.fid %<>% append(low.rep.ident)
rm.sid <- df_est_rep %>% filter(file_id %in% rm.fid) %>% distinct(SampleID) %>% pull()


rep.pearson <- plyr::ddply(df_est_rep %>% filter(!(SampleID %in% rm.sid)), 'SampleID', function(dfsub){
  cat(dfsub$SampleID[1], '...\r')
  X <- dfsub %>%
    column_to_rownames('file_id') %>%
    select(-(BatchID:Cancer_name)) %>% t()
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
  group_by(SampleID) %>%
  summarise(pearson.r.mean = mean(pearson.r)) %>%
  inner_join(rep.pearson, .) %>%
  arrange(pearson.r.mean) %>%
  mutate(SampleID = factor(SampleID, levels = unique(SampleID)))


# 3.samples -----
df_est <- df2 %>% t() %>% as.data.frame() %>%
  rownames_to_column('FileName') %>% 
  inner_join(info2 %>% select(FileName, file_id), .) %>% 
  select(-FileName) %>% 
  full_join(meta2, .)


est.identity <- df_est %>% filter(SampleID != 'Pooling') %>% 
  plyr::ddply(c('Cancer_name', 'Cancer_type'), function(dfsub){
    cat(dfsub$Cancer_name[1], dfsub$Cancer_type[1], '...\r')
    X <- dfsub %>%
      column_to_rownames('file_id') %>%
      select(-(1:年龄)) %>% t()
    ret <- data.frame(file_id = dfsub$file_id,
                      `# proteins` = colSums(!is.na(X)),
                      check.names = F)
    ret$outlier.lower.ingroup <- get_outliers(ret$`# proteins`)[1]
    return(ret)
  })
est.identity %<>%
  group_by(Cancer_name, Cancer_type) %>%
  summarise(identity.mean = mean(`# proteins`),
            identity.min = min(`# proteins`),
            identity.max = max(`# proteins`),
            identity.range.length = identity.max - identity.min) %>%
  inner_join(est.identity, .) %>%
  arrange(identity.mean) %>%
  mutate(Group = str_c(Cancer_type, ' - ', Cancer_name)) %>% 
  mutate(Group = factor(Group, levels = unique(Group)),
         Is.Lower.Ingroup = `# proteins` < outlier.lower.ingroup) %>% 
  inner_join(meta2 %>% select(file_id, Rep_type))
group_order <- est.identity %>% arrange(`# proteins`) %>% distinct(Group) %>% pull()
est.identity %<>% 
  mutate(Group = factor(Group, group_order))
lower_outliers <- est.identity %>%
  filter(!(file_id %in% rm.fid)) %>%
  group_by(Group) %>%
  filter(Is.Lower.Ingroup) %>%
  ungroup() %>% 
  arrange(desc(identity.range.length <= 1000), desc(Group))

df_est_rmp <- df_est %>% filter(SampleID != 'Pooling', !(file_id %in% rm.fid))



# pearson's correlation by Cancer_name
cor_by_cancer <- plyr::ddply(df_est_rmp, 'Cancer_name', function(tmp_df){
  X <- tmp_df %>% column_to_rownames('file_id') %>%
    select(-(1:Cancer_name)) %>% t()
  X <- X[rowSums(!is.na(X)) > 0, ]
  
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
cor_by_cancer %<>%
  group_by(Cancer_name) %>%
  summarise(pearson.r.mean = mean(pearson.r)) %>%
  inner_join(cor_by_cancer, .)


# hclust by Cancer_name
heat_by_cancer <- plyr::dlply(df_est_rmp, 'Cancer_name', function(tmp_df){
  tmp_heat <- tmp_df %>% column_to_rownames('file_id') %>%
    select(-(1:Cancer_name)) %>% t()
  tmp_heat <- tmp_heat[rowSums(!is.na(tmp_heat)) > 0, ]
  tmp_heat[is.na(tmp_heat)] <- 1
  tmp_ann <- tmp_df %>% column_to_rownames('file_id') %>%
    select(1:Cancer_name)
  
  p <- pheatmap::pheatmap(
    tmp_heat, scale = 'none',
    annotation_col = tmp_ann %>%
      rename(Gender = `性别`, Age = `年龄`) %>% 
      select(Cancer_name, Cancer_type, Batch_m, Date, Age, Gender) %>% 
      mutate(Batch_m = factor(Batch_m)),
    clustering_method = 'ward.D2',
    clustering_distance_rows = 'euclidean',
    clustering_distance_cols = 'euclidean',
    show_colnames = T, show_rownames = F,
    filename = str_glue('check_rca_{tmp_df$Cancer_name[1]}_heat.pdf'), width = 10, height = 5
  )
  return(p)
})

membership_phmap <- plyr::ldply(heat_by_cancer, function(phmap){
  hc <- phmap$tree_col
  grp <- cutree(hc, k = 2)
  membership_df <- data.frame(
    file_id  = hc$labels,
    cluster = as.integer(grp[hc$labels]),
    stringsAsFactors = FALSE
  )
  membership_df$order <- match(membership_df$file_id, hc$labels[hc$order])
  membership_df <- membership_df[order(membership_df$order), ]
  return(membership_df)
}, .id = 'Cancer_name')
membership_phmap %<>% inner_join(est.identity)
# 1) Compute T/NT counts and proportions within each Cancer_name × cluster,
#    and determine the dominant type (>50%).
dominant_map <- membership_phmap %>%
  count(Cancer_name, cluster, Cancer_type, name = "n") %>%
  group_by(Cancer_name, cluster) %>%
  mutate(total = sum(n), prop = n / total) %>%
  arrange(Cancer_name, cluster, desc(n)) %>%
  slice_max(n, n = 1, with_ties = TRUE) %>%       # allow ties to identify 50–50 cases
  mutate(
    status = case_when(
      prop > 0.5 ~ "dominant",
      prop == 0.5 ~ "tie"       # exact 50–50
    ),
    dominant_type = ifelse(status == "dominant", Cancer_type, NA_character_)
  ) %>%
  ungroup() %>%
  # For a cluster with a tie, keep one row (dominant_type = NA);
  # if there is a dominant type, keep only that type.
  distinct(Cancer_name, cluster, dominant_type, status, total, .keep_all = FALSE)

# 2) Map the dominant type back to samples and find file_id values that do not match it.
mismatched_file_id <- membership_phmap %>%
  left_join(
    dominant_map %>% select(Cancer_name, cluster, dominant_type, status),
    by = c("Cancer_name", "cluster")
  ) %>%
  # Check consistency only for groups with a dominant type; do not evaluate ties.
  filter(!is.na(dominant_type)) %>%
  mutate(is_match = Cancer_type == dominant_type) %>%
  filter(!is_match) %>%
  select(Cancer_name, file_id, cluster, Cancer_type,
         expected_Cancer_type = dominant_type)

# 3) Provide a detailed count table per group (for review).
by_group_counts <- membership_phmap %>%
  count(Cancer_name, cluster, Cancer_type, name = "n") %>%
  group_by(Cancer_name, cluster) %>%
  mutate(total = sum(n), prop = n / total) %>%
  arrange(Cancer_name, cluster, desc(n)) %>%
  ungroup()

# Inspect the results:
# dominant_map        # Dominant type (or tie) for each Cancer_name × cluster
# mismatched_file_id  # List of file_id values that do not match the dominant type
# by_group_counts     # Detailed counts and proportions

list(`hclust&identity` = membership_phmap,
     Pearson.r = cor_by_cancer,
     Dominant.map = dominant_map,
     mismatched.fid = mismatched_file_id,
     mismatched.fid.stat = by_group_counts,
     mismatched.fid.cor = cor_by_cancer %>% filter(str_detect(ID.pair, str_c(mismatched_file_id$file_id, collapse = '|')))) %>% 
  rio::export('check_rca_heat.xlsx')




# check breast cancer
tmp_df <- df_est %>% filter(SampleID != 'Pooling') %>% 
  filter(Cancer_name == 'breast cancer')
tmp_heat <- tmp_df %>% column_to_rownames('file_id') %>%
  select(-(1:Cancer_name)) %>% t()
tmp_heat <- tmp_heat[rowSums(!is.na(tmp_heat)) > 0, ]
tmp_heat[is.na(tmp_heat)] <- 1
tmp_ann <- tmp_df %>% column_to_rownames('file_id') %>%
  select(1:Cancer_name)

pheatmap::pheatmap(
  tmp_heat, scale = 'none',
  annotation_col = tmp_ann %>%
    rename(Gender = `性别`, Age = `年龄`) %>% 
    select(Cancer_name, Cancer_type, Batch_m, Date, Age, Gender) %>% 
    mutate(Batch_m = factor(Batch_m)),
  clustering_method = 'ward.D2',
  clustering_distance_rows = 'euclidean',
  clustering_distance_cols = 'euclidean',
  show_colnames = F, show_rownames = F,
  filename = 'check_rca_breast_cancer_heat.pdf', width = 10, height = 5
)

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

pdf('atral_qc_dr_breast_cancer.pdf', width = 4*2, height = 3.5*2)
for(nm in names(res.dr)){
  pout <- ggpubr::ggarrange(plotlist = res.dr$plots_pca, nrow = 2, ncol = 2)
  pout2 <- ggpubr::annotate_figure(pout, top = ggpubr::text_grob(str_c('breast cancer'), color = "black", face = "bold", size = 14))
  print(pout2)
  
  pout <- ggpubr::ggarrange(plotlist = res.dr$plots_tsne, nrow = 2, ncol = 2)
  pout2 <- ggpubr::annotate_figure(pout, top = ggpubr::text_grob(str_c('breast cancer'), color = "black", face = "bold", size = 14))
  print(pout2)
  
  pout <- ggpubr::ggarrange(plotlist = res.dr$plots_tsne, nrow = 2, ncol = 2)
  pout2 <- ggpubr::annotate_figure(pout, top = ggpubr::text_grob(str_c('NA impute method: ', nm), color = "black", face = "bold", size = 14))
  print(pout2)
}
graphics.off()



## Output -------
plot_ident_rep <- ggplot(rep.identity) +
  aes(x = `# proteins`, y = SampleID) +
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
  aes(x = pearson.r, y = SampleID) +
  # geom_boxplot(color = '#000000') +
  geom_point(color = '#000000', size = 1.2) +
  labs(x = "Pearson's r", y = "Sample ID", subtitle = "Correlation - same.sample.id") +
  theme_bw() +
  theme(text = element_text(size = 10, color = 'black'))
plot_rep <- ggarrange(plot_ident_rep, plot_pearson_rep, nrow = 1, ncol = 2, widths = c(4, 3))
ggsave('astral_rep.pdf', plot_rep, width = 10, height = 6)

plot_ident <- ggplot(est.identity %>% filter(!(file_id %in% rm.fid))) +
  aes(x = `# proteins`, y = Group) +
  geom_boxplot(data = est.identity %>% filter(identity.range.length <= 1000),
               color = 'black', outlier.color = '#c23190', outlier.size = 2) +
  geom_boxplot(data = est.identity %>% filter(identity.range.length > 1000),
               color = 'black', outlier.color = '#c23190', outlier.size = 2) +
  geom_point(aes(color = Cancer_type), size = 1, alpha = 0.8) +
  geom_text(data = lower_outliers, aes(label = file_id), size = 3) +
  labs(x = "# proteins", y = "Sample ID", subtitle = "Protein identification") +
  ggsci::scale_color_nejm(name = 'Sample type') +
  theme_bw() +
  theme(text = element_text(size = 10, color = 'black'))
ggsave('astral_identity.pdf', plot_ident, width = 16, height = 8)

list(rep.identity = rep.identity,
     rep.pearson = rep.pearson,
     est.identity = est.identity %>% filter(!(file_id %in% rm.fid)),
     est.identity.lower = lower_outliers) %>% 
  rio::export('astral_QC_source.xlsx')


# 4.remove low-identify files -------
rm.fid %<>% append(lower_outliers$file_id %>% setdiff(c('b2_65', 'b3_1', 'b3_19')))

identical(info2$FileName, rownames(dat2))
info2 %<>% filter(!(file_id %in% rm.fid))
dat2 <- dat2[info2$FileName[!info2$Is.pool], ]
# quantile(colSums(!is.na(dat2)))
dim(dat2)

rio::export(info2, '20250809_PUH_RCA_sample_and_pool_1179files_info.xlsx')
write.csv(dat2, 'mapped_dat2_pg_matrix_1152_12500.csv', row.names = T)     














