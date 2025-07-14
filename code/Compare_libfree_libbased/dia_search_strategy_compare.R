rstudioapi::getActiveDocumentContext()$path
setwd('~/GitHub/TPHP/code/Compare_libfree_libbased/')
library(magrittr)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(broom)


TPHP_HOME <- '//172.16.13.136/TPHP/'


cor_color <- colorRampPalette(brewer.pal(11, 'PuOr')[c(11:8, 4:1)])(50)

my_cor_smooth <- function(df, col1, col2, cor_method = "pearson", ggpmisc_label = c("eq", "R2", "p"), plot_title = NULL){
  data_points <- df %>% select(all_of(col1), all_of(col2)) %>% drop_na() %>% nrow()
  df %<>% select(all_of(col1), all_of(col2)) %>% setNames(c('x', 'y'))
  p <- ggplot(df, aes(x = x, y = y))+
    geom_point(na.rm = T)+
    ggpmisc::stat_poly_line(method = 'lm', se = TRUE, level = 0.95,
                            linewidth = 1, color = '#1D439E') +
    ggpmisc::stat_poly_eq(ggpmisc::use_label(ggpmisc_label),
                          method = 'lm', coef.digits = 3, rr.digits = 3, p.digits = 3) +
    # geom_smooth(formula = y ~ x, method = 'lm', color = "#00AFBB", na.rm = T)+
    # ggpubr::stat_cor(method = cor_method, na.rm = T, cor.coef.name = 'rho', vjust = 3)+
    labs(subtitle = str_glue("{col1}: {col2} (x: y, n = {data_points})"))+
    theme_bw()+
    theme(legend.text = element_text(size = 12, color = "black"),legend.position = 'top',
          legend.title = element_text(size = 12, color="black") ,
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.title=element_text(size = 12, hjust = 0.5, color="black"),
          axis.text = element_text(size = 10,color = "black"),
          plot.title = element_text(size = 12, face = "bold"),
          plot.subtitle=element_text(size = 12, hjust = 0, color="black")
    )
  if (!is.na(plot_title)){
    p <- p + labs(title = plot_title)
  }
  return(p)
}

# 1.Input -------
df_info <- readxl::read_excel('../../input/20220708TPHP_1781file_info_edited_v3_2.xlsx') %>% as.data.frame()
df_info %<>% mutate(ID = sapply(FileName, function(x){
  str_extract(x, '^(\\w+?)\\d+.+_(\\d+)\\.d$', group = c(1, 2)) %>%
    str_c(collapse = '_')
}), .before = 1)
df_info %<>% mutate(name = str_remove(FileName, '\\.d$'))


libbased <- arrow::read_parquet(file.path(TPHP_HOME, 'results/swisslib_entraplib_compare/tphp_swiss_lib/report.parquet'))
libfree <- arrow::read_parquet(file.path(TPHP_HOME, 'results/library_free_test_comapre/report.parquet'))

## data filter
pr <- libbased %>% filter(Proteotypic == 1, Global.Q.Value < 0.01, Q.Value < 0.01,
                    Precursor.Quantity != 0, Lib.Q.Value < 0.005)
pg <- pr %>% filter(Lib.PG.Q.Value < 0.005, Global.PG.Q.Value < 0.01)

prfree <- libfree %>% filter(Proteotypic == 1, Global.Q.Value < 0.01, Q.Value < 0.01,
                              Precursor.Quantity != 0, Lib.Q.Value < 0.005)
pgfree <- prfree %>% filter(Lib.PG.Q.Value < 0.005, Global.PG.Q.Value < 0.01)

# 2.Compare -------
length(setdiff(pr$Run, prfree$Run)) # 0
length(setdiff(prfree$Run, pr$Run)) # 1
setdiff(prfree$Run, pr$Run) # K20200825yuel_TPHP_DIA_pool_Slot2-48_1_730

prfree %<>% filter(Run %in% intersect(prfree$Run, pr$Run))
pgfree %<>% filter(Run %in% intersect(pgfree$Run, pg$Run))

pr_mat <- pr %>% select(Precursor.Id, Run, Precursor.Normalised) %>% 
  pivot_wider(id_cols = Precursor.Id, names_from = 'Run', values_from = 'Precursor.Normalised')
prfree_mat <- prfree %>% select(Precursor.Id, Run, Precursor.Normalised) %>% 
  pivot_wider(id_cols = Precursor.Id, names_from = 'Run', values_from = 'Precursor.Normalised')
pg_mat <- pg %>% distinct(Protein.Ids, Run, PG.MaxLFQ) %>% 
  pivot_wider(id_cols = Protein.Ids, names_from = 'Run', values_from = 'PG.MaxLFQ')
pgfree_mat <- pgfree %>% distinct(Protein.Ids, Run, PG.MaxLFQ) %>% 
  pivot_wider(id_cols = Protein.Ids, names_from = 'Run', values_from = 'PG.MaxLFQ')

# rio::export(pg, 'pg.tsv')
# rio::export(pgfree, 'pgfree.tsv')
# iq::process_long_format(
#   input_filename = "pg.tsv",
#   output_filename = "pg_maxlfq.tsv",
#   sample_id = 'Run',
#   primary_id = 'Protein.Ids',
#   annotation_col = c("Genes"),
#   filter_double_less = c("Global.Q.Value" = "0.01", "Global.PG.Q.Value" = "0.01"),
#   pdf_out = 'pg_maxlfq-median-normalize'
# )


## 2.1 pr and pg VennDiagram -----
pg_ls <- list(lib.based = pg_mat$Protein.Ids,
              lib.free = pgfree_mat$Protein.Ids)
pr_ls <- list(lib.based = pr_mat$Precursor.Id,
              lib.free = prfree_mat$Precursor.Id)
venn_pr <- VennDiagram::venn.diagram(x = pr_ls,
                                      resolution = 300,
                                      alpha=rep(0.95, length(pr_ls)),
                                      fill = 'white',
                                      main=stringr::str_glue("Precursors ({length(unique(unlist(pr_ls)))} in total)"),
                                      #sub = rep,
                                      main.cex = 4,
                                      sub.cex = 3,
                                      cex = 4,
                                      cex.lab=4,
                                      cat.cex=4,
                                      imagetype = "tiff",
                                      filename = NULL, disable.logging = T
)
venn_pg <- VennDiagram::venn.diagram(x = pg_ls,
                                       resolution = 300,
                                       alpha=rep(0.95, length(pg_ls)),
                                       # fill=allFills[c(1, 4, 5)],
                                       fill = 'white',
                                       main=stringr::str_glue("Proteins ({length(unique(unlist(pg_ls)))} in total)"),
                                       #sub = rep,
                                       main.cex = 4,
                                       sub.cex = 3,
                                       cex = 4,
                                       cex.lab=4,
                                       cat.cex=4,
                                       imagetype = "tiff",
                                       filename = NULL, disable.logging = T
)
pdf('TPHP_search_method_compare_venndiagram.pdf', width = 20, height = 20)
grid::grid.newpage(); grid::grid.draw(venn_pr)
grid::grid.newpage(); grid::grid.draw(venn_pg)
graphics.off()


## 2.2 pr and pg correlation -----
### just load .RData --
# # pg linear
# pglong1 <- pg_mat %>% 
#   pivot_longer(cols = -Protein.Ids, values_drop_na = T) %>%
#   mutate(Search.method = 'lib.based')
# pglong2 <- pgfree_mat %>% 
#   pivot_longer(cols = -Protein.Ids, values_drop_na = T) %>%
#   mutate(Search.method = 'lib.free')
# pglong <- rbind(pglong1, pglong2) %>%
#   pivot_wider(id_cols = c('Protein.Ids', 'name'), names_from = 'Search.method') %>%
#   drop_na()
# 
# ##option1
# # plot_cor_pg <- plyr::dlply(pglong, 'name', function(dfsub){
# #   cat(dfsub$name[1], '\r')
# #   my_cor_smooth(dfsub, 'TPHP', 'Entrapment', ggpmisc_label = c('eq', 'R2'),
# #                 plot_title = str_extract(dfsub$name[1], '\\\\(\\w+?)\\d+.+_(\\d+)\\.d$', group = c(1, 2)) %>% str_c(collapse = '_'))
# # })
# # pdf('TPHP_entrap_proteinGroup_Linear_correlation.pdf', width = 8, height = 6)
# # for(i in seq_along(plot_cor_pg)){
# #   cat(i, '\r')
# #   print(plot_cor_pg[[i]])
# # }
# # graphics.off()
# 
# ##option2
# lm_pg <- pglong %>%
#   group_by(name) %>%
#   do({
#     model <- lm(lib.free ~ lib.based, data = .)
#     tidy_model <- tidy(model)  # Coefficients (intercept, slope)
#     glance_model <- glance(model)  # Model fit statistics
#     data.frame(
#       Intercept = tidy_model$estimate[1],
#       Slope = tidy_model$estimate[2],
#       p_value = tidy_model$p.value[2],
#       R_squared = glance_model$r.squared,
#       adj_R_squared = glance_model$adj.r.squared
#     )
#   }) %>%
#   ungroup()
# lm_pg %<>%
#   mutate(Response = 'lib.free', Predictor = 'lib.based', .after = name) %>%
#   mutate(adj_p_value = p.adjust(p_value, method = 'BH')) %>%
#   as.data.frame()
# 
# # pr linear
# prlong1 <- pr_mat %>% 
#   pivot_longer(cols = -Precursor.Id, values_drop_na = T) %>%
#   mutate(Search.method = 'lib.based')
# prlong2 <- prfree_mat %>% 
#   pivot_longer(cols = -Precursor.Id, values_drop_na = T) %>%
#   mutate(Search.method = 'lib.free')
# prlong <- rbind(prlong1, prlong2) %>%
#   pivot_wider(id_cols = c('Precursor.Id', 'name'), names_from = 'Search.method') %>%
#   drop_na()
# 
# ##option1
# # plot_cor_pr <- plyr::dlply(prlong, 'name', function(dfsub){
# #   my_cor_smooth(dfsub, 'TPHP', 'Entrapment', ggpmisc_label = c('eq', 'R2'),
# #                 plot_title = str_extract(dfsub$name[1], '\\\\(\\w+?)\\d+.+_(\\d+)\\.d$', group = c(1, 2)) %>% str_c(collapse = '_'))
# # })
# # ggsave('TPHP_entrap_precursor_Linear_correlation.pdf',
# #        ggpubr::ggarrange(plotlist = plot_cor_pr, nrow = 6, ncol = 8),
# #        width = 4*8, height = 4*6)
# 
# ##option2
# lm_pr <- prlong %>%
#   group_by(name) %>%
#   do({
#     model <- lm(lib.free ~ lib.based, data = .)
#     tidy_model <- tidy(model)  # Coefficients (intercept, slope)
#     glance_model <- glance(model)  # Model fit statistics
#     data.frame(
#       Intercept = tidy_model$estimate[1],
#       Slope = tidy_model$estimate[2],
#       p_value = tidy_model$p.value[2],
#       R_squared = glance_model$r.squared,
#       adj_R_squared = glance_model$adj.r.squared
#     )
#   }) %>%
#   ungroup()
# lm_pr %<>%
#   mutate(Response = 'lib.free', Predictor = 'lib.based', .after = name) %>%
#   mutate(adj_p_value = p.adjust(p_value, method = 'BH')) %>%
#   as.data.frame()
# 
# lm_pgpr <- rbind(lm_pg %>% mutate(report = 'pg'), lm_pr %>% mutate(report = 'pr'))
# save(lm_pg, lm_pr, lm_pgpr, file = 'lm_pgpr.RData')
load('lm_pgpr.RData')


## pg plots
pbox_lmpgpr <- ggplot(lm_pgpr) +
  aes(x = report, y = adj_R_squared, fill = report) +
  geom_boxplot(width = 0.6, position = position_dodge()) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.6)) +
  labs(x = 'Report level', y = 'adj.R2') +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0, 0.05, 0)) +
  ggsci::scale_fill_npg() +
  theme_classic() +
  theme(text = element_text(size = 12))
ggsave('TPHP_search_method_compare_box_linear_correlation.pdf', pbox_lmpgpr, width = 6, height = 5)

# arrange
tmp <- lm_pg %>% arrange(adj_R_squared) %>% pull(name)
pglong %<>% mutate(name = factor(name, tmp))
tmp <- lm_pr %>% arrange(adj_R_squared) %>% pull(name)
prlong %<>% mutate(name = factor(name, tmp))

plot_cor_pg <- pglong %>%
  plyr::dlply('name', function(dfsub){
    cat(dfsub$name[1], '\r')
    my_cor_smooth(dfsub, 'lib.based', 'lib.free', ggpmisc_label = c('eq', 'R2'),
                  plot_title = str_extract(dfsub$name[1], '\\\\(\\w+?)\\d+.+_(\\d+)\\.d$', group = c(1, 2)) %>% str_c(collapse = '_'))
  })
pdf('TPHP_search_method_compare_box_linear_correlation_pgfit.pdf', width = 8, height = 6)
for(i in seq_along(plot_cor_pg)){
  cat(i, '\r')
  print(plot_cor_pg[[i]])
}
graphics.off()

# after log2-transformed
plot_cor_pg <- pglong %>%
  mutate(lib.based = log2(lib.based), lib.free = log2(lib.free)) %>% 
  plyr::dlply('name', function(dfsub){
    cat(dfsub$name[1], '\r')
    my_cor_smooth(dfsub, 'lib.based', 'lib.free', ggpmisc_label = c('eq', 'R2'),
                  plot_title = str_extract(dfsub$name[1], '\\\\(\\w+?)\\d+.+_(\\d+)\\.d$', group = c(1, 2)) %>% str_c(collapse = '_'))
  })
pdf('TPHP_search_method_compare_box_linear_correlation_pgfit_afterLog2.pdf', width = 8, height = 6)
for(i in seq_along(plot_cor_pg)){
  cat(i, '\r')
  print(plot_cor_pg[[i]])
}
graphics.off()



## pr plots
plot_cor_pr <- prlong %>%
  plyr::dlply('name', function(dfsub){
    cat(dfsub$name[1], '\r')
    my_cor_smooth(dfsub, 'lib.based', 'lib.free', ggpmisc_label = c('eq', 'R2'),
                  plot_title = str_extract(dfsub$name[1], '\\\\(\\w+?)\\d+.+_(\\d+)\\.d$', group = c(1, 2)) %>% str_c(collapse = '_'))
  })
pdf('TPHP_search_method_compare_box_linear_correlation_prfit.pdf', width = 8, height = 6)
for(i in seq_along(plot_cor_pr)){
  cat(i, '\r')
  print(plot_cor_pr[[i]])
}
graphics.off()

# after log2-transformed
plot_cor_pr <- prlong %>%
  mutate(lib.based = log2(lib.based), lib.free = log2(lib.free)) %>% 
  plyr::dlply('name', function(dfsub){
    cat(dfsub$name[1], '\r')
    my_cor_smooth(dfsub, 'lib.based', 'lib.free', ggpmisc_label = c('eq', 'R2'),
                  plot_title = str_extract(dfsub$name[1], '\\\\(\\w+?)\\d+.+_(\\d+)\\.d$', group = c(1, 2)) %>% str_c(collapse = '_'))
  })
pdf('TPHP_search_method_compare_box_linear_correlation_prfit_afterLog2.pdf', width = 8, height = 6)
for(i in seq_along(plot_cor_pr)){
  cat(i, '\r')
  print(plot_cor_pr[[i]])
}
graphics.off()



## 2.3 pr and pg Q.Value -----
# pg linear
pglong1 <- pg %>% distinct(Protein.Ids, Run, PG.Q.Value) %>% 
  pivot_wider(id_cols = Protein.Ids, names_from = 'Run', values_from = 'PG.Q.Value') %>%
  pivot_longer(cols = -Protein.Ids, values_drop_na = T) %>%
  mutate(Search.method = 'lib.based')
pglong2 <- pgfree %>% distinct(Protein.Ids, Run, PG.Q.Value) %>% 
  pivot_wider(id_cols = Protein.Ids, names_from = 'Run', values_from = 'PG.Q.Value') %>%
  pivot_longer(cols = -Protein.Ids, values_drop_na = T) %>%
  mutate(Search.method = 'lib.free')
pglong <- rbind(pglong1, pglong2) %>%
  pivot_wider(id_cols = c('Protein.Ids', 'name'), names_from = 'Search.method') %>%
  drop_na()

##option1
# plot_cor_pg <- plyr::dlply(pglong, 'name', function(dfsub){
#   cat(dfsub$name[1], '\r')
#   my_cor_smooth(dfsub, 'TPHP', 'Entrapment', ggpmisc_label = c('eq', 'R2'),
#                 plot_title = str_extract(dfsub$name[1], '\\\\(\\w+?)\\d+.+_(\\d+)\\.d$', group = c(1, 2)) %>% str_c(collapse = '_'))
# })
# pdf('TPHP_entrap_proteinGroup_Linear_correlation.pdf', width = 8, height = 6)
# for(i in seq_along(plot_cor_pg)){
#   cat(i, '\r')
#   print(plot_cor_pg[[i]])
# }
# graphics.off()

##option2
lm_pg <- pglong %>%
  group_by(name) %>%
  do({
    model <- lm(lib.free ~ lib.based, data = .)
    tidy_model <- tidy(model)  # Coefficients (intercept, slope)
    glance_model <- glance(model)  # Model fit statistics
    data.frame(
      Intercept = tidy_model$estimate[1],
      Slope = tidy_model$estimate[2],
      p_value = tidy_model$p.value[2],
      R_squared = glance_model$r.squared,
      adj_R_squared = glance_model$adj.r.squared
    )
  }) %>%
  ungroup()
lm_pg %<>%
  mutate(Response = 'lib.free', Predictor = 'lib.based', .after = name) %>%
  mutate(adj_p_value = p.adjust(p_value, method = 'BH')) %>%
  as.data.frame()

# pr linear
prlong1 <- pr %>% select(Precursor.Id, Run, Q.Value) %>% 
  pivot_wider(id_cols = Precursor.Id, names_from = 'Run', values_from = 'Q.Value') %>% 
  pivot_longer(cols = -Precursor.Id, values_drop_na = T) %>%
  mutate(Search.method = 'lib.based')
prlong2 <- prfree %>% select(Precursor.Id, Run, Q.Value) %>% 
  pivot_wider(id_cols = Precursor.Id, names_from = 'Run', values_from = 'Q.Value') %>% 
  pivot_longer(cols = -Precursor.Id, values_drop_na = T) %>%
  mutate(Search.method = 'lib.free')
prlong <- rbind(prlong1, prlong2) %>%
  pivot_wider(id_cols = c('Precursor.Id', 'name'), names_from = 'Search.method') %>%
  drop_na()

##option1
# plot_cor_pr <- plyr::dlply(prlong, 'name', function(dfsub){
#   my_cor_smooth(dfsub, 'TPHP', 'Entrapment', ggpmisc_label = c('eq', 'R2'),
#                 plot_title = str_extract(dfsub$name[1], '\\\\(\\w+?)\\d+.+_(\\d+)\\.d$', group = c(1, 2)) %>% str_c(collapse = '_'))
# })
# ggsave('TPHP_entrap_precursor_Linear_correlation.pdf',
#        ggpubr::ggarrange(plotlist = plot_cor_pr, nrow = 6, ncol = 8),
#        width = 4*8, height = 4*6)

##option2
lm_pr <- prlong %>%
  group_by(name) %>%
  do({
    model <- lm(lib.free ~ lib.based, data = .)
    tidy_model <- tidy(model)  # Coefficients (intercept, slope)
    glance_model <- glance(model)  # Model fit statistics
    data.frame(
      Intercept = tidy_model$estimate[1],
      Slope = tidy_model$estimate[2],
      p_value = tidy_model$p.value[2],
      R_squared = glance_model$r.squared,
      adj_R_squared = glance_model$adj.r.squared
    )
  }) %>%
  ungroup()
lm_pr %<>%
  mutate(Response = 'lib.free', Predictor = 'lib.based', .after = name) %>%
  mutate(adj_p_value = p.adjust(p_value, method = 'BH')) %>%
  as.data.frame()

long_q <- rbind(pglong %>% mutate(report = 'pg') %>% rename(Identifier = Protein.Ids),
                prlong %>% mutate(report = 'pr') %>% rename(Identifier = Precursor.Id))
# save(long_q, file = 'long_q.RData')
load('long_q.RData')


## pg plots
pbox_q <- ggplot(long_q) +
  facet_wrap(~report, scales = 'free') +
  aes(x = -log10(lib.based), y = -log10(lib.free), color = report) +
  geom_point(size = 0.5) +
  # geom_boxplot(width = 0.6, position = position_dodge()) +
  # geom_jitter(position = position_jitterdodge(jitter.width = 0.6)) +
  labs(x = '-Log10.Q.Value-lib.based', y = '-Log10.Q.Value-lib.free') +
  # scale_y_continuous(limits = c(0, 1), expand = c(0, 0, 0.05, 0)) +
  ggsci::scale_color_npg() +
  theme_classic() +
  theme(text = element_text(size = 12))
ggsave('TPHP_search_method_compare_qvalue.pdf', pbox_q, width = 10, height = 5)



## OLD 2.3 pr and pg Q.Value -----
pr_heat <- pg %>% select(Precursor.Id, Run, matches('Q\\.Value')) %>% 
  select(-Peptidoform.Q.Value, -Translated.Q.Value, -Channel.Q.Value) %>% 
  pivot_longer(cols = -c('Precursor.Id', 'Run'), values_drop_na = T) %>%
  unite(remove = F, col = 'rowid', Precursor.Id, name, sep = '.')
ann_row <- pr_heat %>% distinct(rowid, name) %>%
  column_to_rownames('rowid') %>% setNames('Q.Value')
pr_heat %<>%
  pivot_wider(id_cols = rowid, names_from = 'Run', values_from = 'value') %>% 
  mutate(Search.method = 'lib.based') %>% 
  column_to_rownames('rowid') %>% as.matrix()
ann_col <- df_info %>% distinct(name, cancer_type, tissue_name) %>% 
  filter(name %in% colnames(pr_heat)) %>% 
  column_to_rownames('name')
# save(pr_heat, ann_row, ann_col, file = 'pr_heat_with_ann.RData')

pheatmap(
  pr_heat, scale = 'none',
  cluster_cols = T, cluster_rows = T,
  clustering_distance_cols = 'euclidean',
  clustering_distance_rows = 'euclidean',
  clustering_method = 'ward.D2',
  annotation_col = ann_col,
  annotation_row = ann_row,
  # annotation_colors = ann_clr,
  color = cor_color, #breaks = cor_bk,
  show_colnames = F, show_rownames = F,
  fontsize = 8,
  main = 'Q values',
  width = 8, height = 6,
  filename = 'TPHP_search_method_compare_qvalues_heatmap.pdf'
)



# OUTPUT ----------
tbl1 <- lm_pgpr %>% arrange(adj_R_squared) %>% inner_join(df_info)
list(linear.regression = tbl1) %>% 
  rio::export('TPHP_search_method_compare.xlsx')


