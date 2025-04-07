rstudioapi::getActiveDocumentContext()$path
setwd('~/GitHub/TPHP/code/Compare_with_entrapment_lib/')
library(tidyverse)
library(magrittr)
library(RColorBrewer)
library(pheatmap)


TPHP_HOME <- '//172.16.13.136/TPHP/'

target_decoy_color <- ggsci::pal_d3()(2) %>% setNames(c('Target', 'Decoy'))

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
tphp_pg <- rio::import('//172.16.13.136/TPHP/results/swisslib_entraplib_compare/tphp_swiss_lib/report.pg_matrix.tsv')
entra_pg <- rio::import('//172.16.13.136/TPHP/results/swisslib_entraplib_compare/entrapment_lib/report.pg_matrix.tsv')
tphp_pr <- rio::import('//172.16.13.136/TPHP/results/swisslib_entraplib_compare/tphp_swiss_lib/report.pr_matrix.tsv')
entra_pr <- rio::import('//172.16.13.136/TPHP/results/swisslib_entraplib_compare/entrapment_lib/report.pr_matrix.tsv')

df_info <- readxl::read_excel('../../input/20220706TPHP_1781file_117pool_info_onlyfilename-id_pm_swissprot1025.xlsx', col_types = c(rep('guess', 21), rep('numeric', 13479)))
df_info %<>% select(1:21)

# 2.For entrapmant pg report -------
# a posteriori FDR
sum(!(entra_pg$Protein.Group %in% libpep$`Protein ID`)) / nrow(entra_pg) * 100 # 1.14387 (%)
sum(!(entra_pr$Stripped.Sequence %in% libpep$Peptide)) / nrow(entra_pr) * 100 # 0.06751949 (%)

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
  geom_jitter(data = df_td_box %>% filter(!Is.Outlier),
              position = position_jitterdodge(dodge.width = 0.6, seed = 100)) +
  geom_text(data = df_td_box %>% filter(Is.Outlier),
            aes(label = format_custom(ratio)), hjust = -1) +
  labs(x = str_glue("DIA runs\n(n={nrow(df_td_box) / 2})"), y = 'Posterior run-specific discovery rate') +
  scale_fill_manual(values = target_decoy_color) +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  theme(text = element_text(size = 12))

plot_td <- ggpubr::ggarrange(plot_target_decoy, plot_target_decoy_box, nrow = 1, ncol = 2,
                             widths = c(3, 2), common.legend = T, labels = c('a', 'b'))
ggsave('entrapment_lib_result.pdf', plot_td, width = 8, height = 6)


# 3.Compare two reports ------
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
# pg linear
pglong1 <- tphp_pg %>% select(-(Protein.Names:First.Protein.Description)) %>% 
  pivot_longer(cols = -Protein.Group, values_drop_na = T) %>% 
  mutate(Library = 'TPHP')
pglong2 <- entra_pg %>% select(-(Protein.Names:First.Protein.Description)) %>% 
  pivot_longer(cols = -Protein.Group, values_drop_na = T) %>% 
  mutate(Library = 'Entrapment')
pglong <- rbind(pglong1, pglong2) %>% 
  pivot_wider(id_cols = c('Protein.Group', 'name'), names_from = 'Library') %>% 
  drop_na()

plot_cor_pg <- plyr::dlply(pglong, 'name', function(dfsub){
  my_cor_smooth(dfsub, 'TPHP', 'Entrapment', ggpmisc_label = c('eq', 'R2'),
                plot_title = str_extract(dfsub$name[1], '\\\\(\\w+?)\\d+.+_(\\d+)\\.d$', group = c(1, 2)) %>% str_c(collapse = '_'))
})
ggsave('TPHP_entrap_proteinGroup_Linear_correlation.pdf',
       ggpubr::ggarrange(plotlist = plot_cor_pg, nrow = 6, ncol = 8),
       width = 4*8, height = 4*6)

# pr linear
prlong1 <- tphp_pr %>% select(-(Protein.Group:Precursor.Charge)) %>% 
  pivot_longer(cols = -Precursor.Id, values_drop_na = T) %>% 
  mutate(Library = 'TPHP')
prlong2 <- entra_pr %>% select(-(Protein.Group:Precursor.Charge)) %>% 
  pivot_longer(cols = -Precursor.Id, values_drop_na = T) %>% 
  mutate(Library = 'Entrapment')
prlong <- rbind(prlong1, prlong2) %>% 
  pivot_wider(id_cols = c('Precursor.Id', 'name'), names_from = 'Library') %>% 
  drop_na()

plot_cor_pr <- plyr::dlply(prlong, 'name', function(dfsub){
  my_cor_smooth(dfsub, 'TPHP', 'Entrapment', ggpmisc_label = c('eq', 'R2'),
                plot_title = str_extract(dfsub$name[1], '\\\\(\\w+?)\\d+.+_(\\d+)\\.d$', group = c(1, 2)) %>% str_c(collapse = '_'))
})
ggsave('TPHP_entrap_precursor_Linear_correlation.pdf',
       ggpubr::ggarrange(plotlist = plot_cor_pr, nrow = 6, ncol = 8),
       width = 4*8, height = 4*6)


# correlation matrix
X_pg1 <- tphp_pg
colnames(X_pg1)[-(1:4)] %<>% sapply(function(x){
  str_extract(x, '\\\\(\\w+?)\\d+.+_(\\d+)\\.d$', group = c(1, 2)) %>%
    str_c(collapse = '_') %>% str_c('.tphp')
})
X_pg2 <- entra_pg
colnames(X_pg2)[-(1:4)] %<>% sapply(function(x){
  str_extract(x, '\\\\(\\w+?)\\d+.+_(\\d+)\\.d$', group = c(1, 2)) %>%
    str_c(collapse = '_') %>% str_c('.entrap')
})
X_pg <- X_pg1 %>% inner_join(X_pg2) %>% 
  column_to_rownames('Protein.Group') %>% 
  select(-(Protein.Names:First.Protein.Description)) %>% 
  as.matrix()
cor_pg <- X_pg %>% cor(use = 'pairwise.complete.obs', method = 'spearman')
ann <- data.frame(row.names = c(colnames(X_pg1)[-(1:4)], colnames(X_pg2)[-(1:4)]),
                  Library = c(rep('TPHP', ncol(X_pg1)-4),
                              rep('Entrapment', ncol(X_pg2)-4)))
cor_color <- colorRampPalette(brewer.pal(11, 'PuOr')[c(11:8, 4:1)])(50)
cor_bk <- unique(c(seq(0.4, 1, length=50)))
pheatmap(
  cor_pg, scale = 'none',
  cluster_cols = T, cluster_rows = T,
  clustering_distance_cols = 'euclidean',
  clustering_distance_rows = 'euclidean',
  clustering_method = 'ward.D2',
  annotation_col = ann,
  annotation_row = ann,
  color = cor_color, breaks = cor_bk,
  show_colnames = F, show_rownames = F,
  main = 'Spearman correlation - proteins',
  width = 10, height = 10,
  filename = 'TPHP_entrap_proteinGroup_spearman_matrix.pdf'
)








