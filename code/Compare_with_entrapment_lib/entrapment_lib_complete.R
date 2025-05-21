rstudioapi::getActiveDocumentContext()$path
setwd('~/GitHub/TPHP/code/Compare_with_entrapment_lib/')
source('../../source/my_fun.R')
library(tidyverse)
library(magrittr)
library(RColorBrewer)
library(pheatmap)
library(broom)


TPHP_HOME <- '//172.16.13.136/TPHP/'

target_decoy_color <- ggsci::pal_d3()(2) %>% setNames(c('Target', 'Decoy'))
library_color <- ggsci::pal_d3()(2) %>% setNames(c('TPHP', 'Entrapment'))
cancer_color <- ggsci::pal_nejm()(5) %>% setNames(c('C_pool', 'N_pool', 'Normal', 'adjacent', 'carcinoma'))
df_color <- rio::import('//172.16.13.136/share/members/jiangwenhao/TPHP/input/PUH_tissue_colorset_20230210.xlsx')
tissue_color <- str_c('#', df_color$color[1:51]) %>% setNames(df_color$tissue[1:51])
tissue_color %<>% append(c('FO' = '#23087C'))

format_custom <- function(x, cutoff = 0.01) {
  format_custom_individual <- function(x, cutoff = 0.01){
    if (x >= cutoff) {
      sprintf('%.3f', x)
    } else {
      as.character(signif(x, 3))
      # format(x, scientific = F)
    }
  }
  sapply(x, format_custom_individual, cutoff = cutoff)
  
}
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
# # 0.Entrapment library prepare ----------
# lib_tphp <- rio::import('//172.16.13.136/TPHP/library/TPHPlib_frag1025_swissprot_final.tsv')
# lib_entrap <- rio::import('//172.16.13.136/share/members/jiangwenhao/20250401_GNHSF_lib_tims90minIGC_88/IGC_humanswiss_irt_contam_88_ddafile_NEW_rmone_20220102.tsv')
# setdiff(names(lib_tphp), names(lib_entrap)) # AverageExperimentalRetentionTime
# 
# 
# length(unique(lib_tphp$PeptideSequence)) # 484,391
# length(unique(lib_entrap$PeptideSequence)) # 506,470
# length(intersect(lib_tphp$PeptideSequence, lib_entrap$PeptideSequence)) # 9184
# 
# 
# 
# # bind
# lib_entrap$AverageExperimentalRetentionTime <- NA
# lib_bind <- lib_entrap %>% 
#   filter(!(PeptideSequence %in% lib_tphp$PeptideSequence)) %>% 
#   rbind(lib_tphp, .)
# rio::export(lib_bind, 'entrapment_lib.tsv')

# 1.Data readin ------
libpep <- rio::import(file.path(TPHP_HOME, 'TPL/libs/rlt_all/TPHP_swissprot/peptide.tsv'))
tphp_pg <- rio::import(file.path(TPHP_HOME, 'results/rlt_combine_swissprot1025/TPHP_swissprot1025_report.pg_matrix.tsv')) %>% setNames(., str_remove(names(.), '.+\\\\'))
entra_pg <- rio::import(file.path(TPHP_HOME, 'results/rlt_combine_entraplib_G3/report.pg_matrix.tsv')) %>% setNames(., str_remove(names(.), '.+\\\\'))
tphp_pr <- rio::import(file.path(TPHP_HOME, 'results/rlt_combine_swissprot1025/TPHP_swissprot1025_report.pr_matrix.tsv')) %>% setNames(., str_remove(names(.), '.+\\\\'))
entra_pr <- rio::import(file.path(TPHP_HOME, 'results/rlt_combine_entraplib_G3/report.pr_matrix.tsv')) %>% setNames(., str_remove(names(.), '.+\\\\'))

df_info <- readxl::read_excel('../../input/20220708TPHP_1781file_info_edited_v3_2.xlsx') %>% as.data.frame()
# df_info %<>% mutate(ID = sapply(FileName, function(x){
#   str_extract(x, '^(\\w+?)\\d+.+_(\\d+)\\.d$', group = c(1, 2)) %>%
#     str_c(collapse = '_')
# }), .before = 1)


# 2.For entrapmant pg report -------
# a posteriori FDR
sum(!(entra_pg$Protein.Group %in% libpep$`Protein ID`)) / nrow(entra_pg) * 100 # 5.467587 (%)
sum(!(entra_pr$Stripped.Sequence %in% libpep$Peptide)) / nrow(entra_pr) * 100 # 0.2709605 (%)

# df_td <- data.frame(Identifier = ifelse(entra_pg$Protein.Group %in% libpep$`Protein ID`, 'Target', 'Decoy'))
# ggplot(df_td) +
#   aes(x = Identifier) +
#   geom_bar(fill = '#000000') +
#   labs(x = 'Identifier', y = 'Count') +
#   scale_y_continuous(limits = c(0, NA)) +
#   theme_classic() +
#   theme(text = element_text(size = 12))

df_td <- entra_pg %>% # target-decoy
  select(-(Protein.Names:First.Protein.Description)) %>% 
  mutate(Type = ifelse(Protein.Group %in% libpep$`Protein ID`, 'Target', 'Decoy'),
         .before = 2) %>% 
  pivot_longer(cols = -c('Protein.Group', 'Type'),
               names_to = 'FileName', values_drop_na = T) %>% 
  count(FileName, Type) %>% 
  arrange(desc(n)) %>% 
  mutate(FileName = factor(FileName, levels = unique(FileName)),
         Type = factor(Type, levels = c('Target', 'Decoy'))) %>%
  group_by(FileName) %>%
  mutate(ratio = n / sum(n)) %>% 
  as.data.frame()
plot_target_decoy <- ggplot(df_td) +
  facet_wrap(~Type, nrow = 2, ncol = 1, scales = 'free') +
  aes(x = FileName, y = n, fill = Type) +
  geom_col(color = '#000000', width = 1, position = position_stack()) +
  labs(x = 'DIA runs', y = 'Count') +
  scale_y_continuous(limits = c(0, NA), expand = c(0, 0, 0.05, 0)) +
  scale_fill_manual(values = target_decoy_color) +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  theme(text = element_text(size = 12))


df_decoy <- df_td %>% filter(Type == 'Decoy')
df_decoy %<>% ggpubr::get_summary_stats(ratio) %>%
  select(min:ci) %>%
  mutate(Type = 'Decoy', .before = 1) %>% 
  mutate(box.lower = q1 - 1.5*iqr,
         box.upper = q3 + 1.5*iqr) %>% 
  inner_join(df_td, .) %>% 
  mutate(Is.Outlier = !between(ratio, box.lower, box.upper))
# plot_decoy_box <- ggplot(df_decoy) +
#   aes(x = Type, y = ratio, fill = Type) +
#   geom_boxplot(width = 0.6, position = position_dodge()) +
#   geom_jitter(data = df_decoy %>% filter(!Is.Outlier),
#               position = position_jitterdodge(dodge.width = 0.6, seed = 100)) +
#   geom_text(data = df_decoy %>% filter(Is.Outlier),
#             aes(label = format_custom(ratio)), hjust = -1) +
#   labs(x = str_glue("DIA runs\n(n={nrow(df_decoy)})"), y = 'Posterior run-specific FDR') +
#   scale_y_continuous(limits = c(0, NA), expand = c(0, 0, 0.05, 0)) +
#   scale_fill_manual(values = target_decoy_color) +
#   theme_classic() +
#   theme(axis.text.x = element_blank()) +
#   theme(text = element_text(size = 12))

df_target <- df_td %>% filter(Type == 'Target')
df_target %<>% ggpubr::get_summary_stats(ratio) %>%
  select(min:ci) %>%
  mutate(Type = 'Target', .before = 1) %>% 
  mutate(box.lower = q1 - 1.5*iqr,
         box.upper = q3 + 1.5*iqr) %>% 
  inner_join(df_td, .) %>% 
  mutate(Is.Outlier = !between(ratio, box.lower, box.upper))
# plot_target_box <- ggplot(df_target) +
#   aes(x = Type, y = ratio, fill = Type) +
#   geom_boxplot(width = 0.6, position = position_dodge()) +
#   geom_jitter(data = df_target %>% filter(!Is.Outlier),
#               position = position_jitterdodge(dodge.width = 0.6, seed = 100)) +
#   geom_text(data = df_target %>% filter(Is.Outlier),
#             aes(label = format_custom(ratio)), hjust = -1) +
#   labs(x = str_glue("DIA runs\n(n={nrow(df_target)})"), y = 'Posterior run-specific TDR') +
#   scale_y_continuous(limits = c(0.986, NA), expand = c(0, 0, 0.05, 0)) +
#   scale_fill_manual(values = target_decoy_color) +
#   theme_classic() +
#   theme(axis.text.x = element_blank()) +
#   theme(text = element_text(size = 12))
df_td_box <- rbind(df_target, df_decoy) %>%
  mutate(Type = factor(Type, levels = c('Target', 'Decoy')))
plot_target_decoy_box <- ggplot(df_td_box) +
  facet_wrap(~Type, nrow = 2, ncol = 1, scales = 'free') +
  aes(x = Type, y = ratio, fill = Type) +
  geom_boxplot(width = 0.6, position = position_dodge()) +
  # geom_jitter(data = df_td_box %>% filter(!Is.Outlier),
  #             position = position_jitterdodge(dodge.width = 0.6, seed = 100)) +
  ggrepel::geom_text_repel(
    data = df_td_box %>% filter(Is.Outlier) %>% 
      filter(Type == 'Target' & ratio < 0.99 |
               Type == 'Decoy' & ratio > 0.01),
    aes(label = format_custom(ratio)), hjust = -1
  ) +
  labs(x = str_glue("DIA runs\n(n={nrow(df_td_box) / 2})"), y = 'Posterior run-specific discovery rate') +
  scale_fill_manual(values = target_decoy_color) +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  theme(text = element_text(size = 12))

plot_td <- ggpubr::ggarrange(plot_target_decoy, plot_target_decoy_box, nrow = 1, ncol = 2,
                             widths = c(3, 2), common.legend = T, labels = c('a', 'b'))
ggsave('entrapment_lib_result.pdf', plot_td, width = 8, height = 6)


# 3.Compare two reports ------
intersect(colnames(tphp_pg), colnames(entra_pg)) # 1961; -> 1957 files
tphp_pg <- tphp_pg[, intersect(colnames(tphp_pg), colnames(entra_pg))]
entra_pg <- entra_pg[, intersect(colnames(tphp_pg), colnames(entra_pg))]
tphp_pr <- tphp_pr[, intersect(colnames(tphp_pr), colnames(entra_pr))]
entra_pr <- entra_pr[, intersect(colnames(tphp_pr), colnames(entra_pr))]

## 3.1 pr and pg VennDiagram -----
prot_ls <- list(tphp = tphp_pg$Protein.Group,
                entrap = entra_pg$Protein.Group)
pep_ls <- list(tphp = tphp_pr$Stripped.Sequence,
               entrap = entra_pr$Stripped.Sequence)

venn_pep <- VennDiagram::venn.diagram(x = pep_ls,
                                resolution = 300,
                                alpha=rep(0.95, length(pep_ls)),
                                fill = 'white',
                                main=stringr::str_glue("Peptides ({length(unique(unlist(pep_ls)))} in total)"),
                                #sub = rep,
                                main.cex = 4,
                                sub.cex = 3,
                                cex = 4,
                                cex.lab=4,
                                cat.cex=4,
                                imagetype = "tiff",
                                filename = NULL, disable.logging = T
)
venn_prot <- VennDiagram::venn.diagram(x = prot_ls,
                                resolution = 300,
                                alpha=rep(0.95, length(prot_ls)),
                                # fill=allFills[c(1, 4, 5)],
                                fill = 'white',
                                main=stringr::str_glue("Proteins ({length(unique(unlist(prot_ls)))} in total)"),
                                #sub = rep,
                                main.cex = 4,
                                sub.cex = 3,
                                cex = 4,
                                cex.lab=4,
                                cat.cex=4,
                                imagetype = "tiff",
                                filename = NULL, disable.logging = T
)
pdf('TPHP_entrap_VennDiagram.pdf', width = 20, height = 20)
grid::grid.newpage(); grid::grid.draw(venn_pep)
grid::grid.newpage(); grid::grid.draw(venn_prot)
graphics.off()


## 3.2 pr and pg correlation -----
#### just load previous results
# # pg linear
# pglong1 <- tphp_pg %>% select(-(Protein.Names:First.Protein.Description)) %>% 
#   pivot_longer(cols = -Protein.Group, values_drop_na = T) %>% 
#   mutate(Library = 'TPHP')
# pglong2 <- entra_pg %>% select(-(Protein.Names:First.Protein.Description)) %>% 
#   pivot_longer(cols = -Protein.Group, values_drop_na = T) %>% 
#   mutate(Library = 'Entrapment')
# pglong <- rbind(pglong1, pglong2) %>% 
#   pivot_wider(id_cols = c('Protein.Group', 'name'), names_from = 'Library') %>% 
#   drop_na()
# 
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
# 
# #option2
# lm_pg <- pglong %>%
#   group_by(name) %>%
#   do({
#     model <- lm(Entrapment ~ TPHP, data = .)
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
#   mutate(Response = 'Entrapment', Predictor = 'TPHP', .after = name) %>% 
#   mutate(adj_p_value = p.adjust(p_value, method = 'BH')) %>% 
#   as.data.frame()
# 
# 
# 
# # pr linear
# ### Run on G3 --
# # prlong1 <- tphp_pr %>% select(-(Protein.Group:Precursor.Charge)) %>% 
# #   pivot_longer(cols = -Precursor.Id, values_drop_na = T) %>% 
# #   mutate(Library = 'TPHP')
# # prlong2 <- entra_pr %>% select(-(Protein.Group:Precursor.Charge)) %>% 
# #   pivot_longer(cols = -Precursor.Id, values_drop_na = T) %>% 
# #   mutate(Library = 'Entrapment')
# # prlong <- rbind(prlong1, prlong2) %>% 
# #   pivot_wider(id_cols = c('Precursor.Id', 'name'), names_from = 'Library') %>% 
# #   drop_na()
# # 
# # ##option1
# # # plot_cor_pr <- plyr::dlply(prlong, 'name', function(dfsub){
# # #   my_cor_smooth(dfsub, 'TPHP', 'Entrapment', ggpmisc_label = c('eq', 'R2'),
# # #                 plot_title = str_extract(dfsub$name[1], '\\\\(\\w+?)\\d+.+_(\\d+)\\.d$', group = c(1, 2)) %>% str_c(collapse = '_'))
# # # })
# # # ggsave('TPHP_entrap_precursor_Linear_correlation.pdf',
# # #        ggpubr::ggarrange(plotlist = plot_cor_pr, nrow = 6, ncol = 8),
# # #        width = 4*8, height = 4*6)
# # 
# # #option2
# # lm_pr <- prlong %>%
# #   group_by(name) %>%
# #   do({
# #     model <- lm(Entrapment ~ TPHP, data = .)
# #     tidy_model <- tidy(model)  # Coefficients (intercept, slope)
# #     glance_model <- glance(model)  # Model fit statistics
# #     data.frame(
# #       Intercept = tidy_model$estimate[1],
# #       Slope = tidy_model$estimate[2],
# #       p_value = tidy_model$p.value[2],
# #       R_squared = glance_model$r.squared,
# #       adj_R_squared = glance_model$adj.r.squared
# #     )
# #   }) %>%
# #   ungroup()
# # lm_pr %<>%
# #   mutate(Response = 'Entrapment', Predictor = 'TPHP', .after = name) %>% 
# #   mutate(adj_p_value = p.adjust(p_value, method = 'BH')) %>% 
# #   as.data.frame()
# load('lm_pr.RData')
# 
# lm_pgpr <- rbind(lm_pg %>% mutate(report = 'pg'), lm_pr %>% mutate(report = 'pr'))
# save(lm_pg, lm_pr, lm_pgpr, file = 'lm_pgpr.RData')
# 
# 
load('lm_pgpr.RData')

# output
volcano_lm_pg <- lm_pg %>% 
  mutate(`-Log10.adj.p` = ifelse(adj_p_value == 0, 1e-300, adj_p_value),
         `-Log10.adj.p` = -log10(`-Log10.adj.p`)) %>% 
  ggplot()+
  aes(x = R_squared, `-Log10.adj.p`) +
  geom_point() +
  # geom_density_2d_filled(alpha = 0.8, contour_var = "density") + 
  labs(x = 'R squared', y = '-Log10(adjusted p value)', subtitle = 'PG.report') +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 300), breaks = c(0, 100, 200, 300), labels = c(0, 100, 200, 'Inf')) +
  # scale_color_viridis_d(option = "D") +
  theme_bw() + theme(text = element_text(size = 10))
volcano_lm_pg <- ggExtra::ggMarginal(volcano_lm_pg, type = "density", margins = "both", size = 5, fill = 'orange3')

volcano_lm_pr <- lm_pr %>% 
  mutate(`-Log10.adj.p` = ifelse(adj_p_value == 0, 1e-300, adj_p_value),
         `-Log10.adj.p` = -log10(`-Log10.adj.p`)) %>% 
  ggplot()+
  aes(x = R_squared, `-Log10.adj.p`) +
  geom_point() +
  # geom_density_2d_filled(alpha = 0.8, contour_var = "density") + 
  labs(x = 'R squared', y = '-Log10(adjusted p value)', subtitle = 'PR.report') +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 300), breaks = c(0, 100, 200, 300), labels = c(0, 100, 200, 'Inf')) +
  # scale_color_viridis_d(option = "D") +
  theme_bw() + theme(text = element_text(size = 10))
volcano_lm_pr <- ggExtra::ggMarginal(volcano_lm_pr, type = "density", margins = "both", size = 5, fill = 'orange3')

volcano_lm <- ggpubr::ggarrange(volcano_lm_pg, volcano_lm_pr)
ggsave('Entrapment_TPHP_lm_volcano.pdf', volcano_lm, width = 8, height = 4)
list(entrapment_tphp_lm_pg = lm_pg, entrapment_tphp_lm_pr = lm_pr) %>% 
  rio::export('Entrapment_TPHP_lm_volcano.xlsx')



### correlation matrix ----------
X_pg1 <- tphp_pg %>%
  rename(`N20210825yuel_TPHP_nail_pool_Slot1-6_1_7467.d` = `N20210825yuel_nail_pool_Slot1-6_1_7467.d`)
colnames(X_pg1)[-(1:4)] %<>% sapply(function(x){
  str_extract(x, '^(\\w+?).+TPHP_(.+)_Slot.+_(\\d+)\\.d$', group = c(1, 2, 3)) %>%
    str_c(collapse = '_') %>% str_c('.tphp')
})
X_pg2 <- entra_pg %>%
  rename(`N20210825yuel_TPHP_nail_pool_Slot1-6_1_7467.d` = `N20210825yuel_nail_pool_Slot1-6_1_7467.d`)
colnames(X_pg2)[-(1:4)] %<>% sapply(function(x){
  str_extract(x, '^(\\w+?).+TPHP_(.+)_Slot.+_(\\d+)\\.d$', group = c(1, 2, 3)) %>%
    str_c(collapse = '_') %>% str_c('.entrap')
})
X_pg <- X_pg1 %>% inner_join(X_pg2) %>% 
  column_to_rownames('Protein.Group') %>% 
  select(-(Protein.Names:First.Protein.Description)) %>% 
  as.matrix()
# cor_pg <- X_pg %>% cor(use = 'pairwise.complete.obs', method = 'spearman')
# save(cor_pg, file = 'cor_pg.RData')
load('cor_pg.RData')
ann <- data.frame(FileName = c(colnames(tphp_pg)[-(1:4)], colnames(entra_pg)[-(1:4)]),
                  ID = c(colnames(X_pg1)[-(1:4)], colnames(X_pg2)[-(1:4)]),
                  Library = c(rep('TPHP', ncol(X_pg1)-4),
                              rep('Entrapment', ncol(X_pg2)-4))) %>% 
  inner_join(df_info %>% select(FileName, cancer_type, tissue_type, DDA_lib_type), .) %>% 
  column_to_rownames('ID') %>% 
  select(-FileName)
df_abbr <- read.delim('//172.16.13.136/share/members/jiangwenhao/TPHP/input/sample_types_abbr_20230113.txt', stringsAsFactors = F, check.names = F, na.strings = '')
ann_abbr <- get_abbr(ann, 'DDA_lib_type', df_abbr = df_abbr) %>% 
  select(DDA_lib_type, cancer_type, Library)

ann_clr <- list(Library = library_color, cancer_type = cancer_color, DDA_lib_type = tissue_color)

cor_color <- colorRampPalette(brewer.pal(11, 'PuOr')[c(11:8, 4:1)])(50)
cor_bk <- unique(c(seq(0.4, 1, length=50)))
pheatmap(
  cor_pg, scale = 'none',
  cluster_cols = T, cluster_rows = T,
  clustering_distance_cols = 'euclidean',
  clustering_distance_rows = 'euclidean',
  clustering_method = 'ward.D2',
  annotation_col = ann_abbr,
  annotation_row = ann_abbr,
  annotation_colors = ann_clr,
  color = cor_color, breaks = cor_bk,
  show_colnames = F, show_rownames = F,
  fontsize = 8,
  main = 'Spearman correlation - proteins',
  width = 8, height = 6,
  filename = 'TPHP_entrap_proteinGroup_spearman_matrix.pdf'
)






X_pr1 <- tphp_pr %>%
  rename(`N20210825yuel_TPHP_nail_pool_Slot1-6_1_7467.d` = `N20210825yuel_nail_pool_Slot1-6_1_7467.d`)
colnames(X_pr1)[-(1:10)] %<>% sapply(function(x){
  str_extract(x, '^(\\w+?).+TPHP_(.+)_Slot.+_(\\d+)\\.d$', group = c(1, 2, 3)) %>%
    str_c(collapse = '_') %>% str_c('.tphp')
})
X_pr2 <- entra_pr %>%
  rename(`N20210825yuel_TPHP_nail_pool_Slot1-6_1_7467.d` = `N20210825yuel_nail_pool_Slot1-6_1_7467.d`)
colnames(X_pr2)[-(1:10)] %<>% sapply(function(x){
  str_extract(x, '^(\\w+?).+TPHP_(.+)_Slot.+_(\\d+)\\.d$', group = c(1, 2, 3)) %>%
    str_c(collapse = '_') %>% str_c('.entrap')
})
X_pr <- X_pr1 %>% inner_join(X_pr2) %>% 
  column_to_rownames('Precursor.Id') %>% 
  select(-(Protein.Group:Precursor.Charge)) %>% 
  as.matrix()
# cor_pr <- X_pr %>% cor(use = 'pairwise.complete.obs', method = 'spearman')
# save(cor_pr, file = 'cor_pr.RData')
load('cor_pr.RData')
ann <- data.frame(FileName = c(colnames(tphp_pr)[-(1:10)], colnames(entra_pr)[-(1:10)]),
                  ID = c(colnames(X_pr1)[-(1:10)], colnames(X_pr2)[-(1:10)]),
                  Library = c(rep('TPHP', ncol(X_pr1)-10),
                              rep('Entrapment', ncol(X_pr2)-10))) %>% 
  inner_join(df_info %>% select(FileName, cancer_type, tissue_type, DDA_lib_type), .) %>% 
  column_to_rownames('ID') %>% 
  select(-FileName)
df_abbr <- read.delim('//172.16.13.136/share/members/jiangwenhao/TPHP/input/sample_types_abbr_20230113.txt', stringsAsFactors = F, check.names = F, na.strings = '')
ann_abbr <- get_abbr(ann, 'DDA_lib_type', df_abbr = df_abbr) %>% 
  select(DDA_lib_type, cancer_type, Library)

ann_clr <- list(Library = library_color, cancer_type = cancer_color, DDA_lib_type = tissue_color)

cor_color <- colorRampPalette(brewer.pal(11, 'PuOr')[c(11:8, 4:1)])(50)
cor_bk <- unique(c(seq(0.3, 1, length=50)))
pheatmap(
  cor_pr, scale = 'none',
  cluster_cols = T, cluster_rows = T,
  clustering_distance_cols = 'euclidean',
  clustering_distance_rows = 'euclidean',
  clustering_method = 'ward.D2',
  annotation_col = ann_abbr,
  annotation_row = ann_abbr,
  annotation_colors = ann_clr,
  color = cor_color, breaks = cor_bk,
  show_colnames = F, show_rownames = F,
  fontsize = 8,
  main = 'Spearman correlation - precursors',
  width = 8, height = 6,
  filename = 'TPHP_entrap_precursor_spearman_matrix.pdf'
)

