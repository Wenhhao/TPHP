setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(magrittr)
library(tidyverse)



df <- rio::import('X:/results/rlt_pool/report.pg_matrix.tsv')
pg_mat <- df %>% select(-(Protein.Ids:First.Protein.Description)) %>% 
  column_to_rownames('Protein.Group') %>% as.matrix()
colnames(pg_mat) %<>% str_remove_all('^.+\\\\|\\.d$')
df_pr <- rio::import('X:/results/rlt_pool/report.pr_matrix.tsv')
pr_mat <- df_pr %>% select(-(1:Precursor.Charge)) %>% 
  column_to_rownames('Precursor.Id') %>% as.matrix()
colnames(pr_mat) %<>% str_remove_all('^.+\\\\|\\.d$')



info <- data.frame(Run = colnames(pg_mat))
info %<>% mutate(Instrument = str_extract(Run, '^[A-Z]+'),
                 Date = str_extract(Run, '^[A-Z]+(\\d+)', group = 1) %>% as.numeric(),
                 id = str_extract(Run, '\\d+$'),
                 UniqueID = str_c(Instrument, Date, '_', id)) 
# colnames(pg_mat) %>% identical(info$Run) # TRUE
colnames(pg_mat) <- info$UniqueID
# colnames(pr_mat) %>% identical(info$Run) # TRUE
colnames(pr_mat) <- info$UniqueID



ann <- info %>% select(-Run, -id) %>% column_to_rownames('UniqueID')

cor_pg <- cor(pg_mat, use = 'pairwise.complete.obs', method = 'spearman')
cor_pr <- cor(pr_mat, use = 'pairwise.complete.obs', method = 'spearman')

p_pgcor <- pheatmap::pheatmap(cor_pg, scale = 'none',
                              annotation_row = ann, annotation_col = ann,
                              display_numbers = T, number_color = 'black',
                              fontsize = 8,
                              main = "Spearman's correlation -pool pg") %>% 
  ggplotify::as.ggplot()
p_prcor <- pheatmap::pheatmap(cor_pr, scale = 'none',
                              annotation_row = ann, annotation_col = ann,
                              display_numbers = T, number_color = 'black',
                              fontsize = 8,
                              main = "Spearman's correlation -pool pr") %>% 
  ggplotify::as.ggplot()
p_cor <- ggpubr::ggarrange(p_pgcor, p_prcor, nrow = 1, ncol = 2)
ggsave('QC_35pool_spearman.pdf', p_cor, width = 22, height = 10)

df_ident <- data.frame(pg = colSums(!is.na(pg_mat), na.rm = T),
                       pr = colSums(!is.na(pr_mat), na.rm = T)) %>% 
  rownames_to_column('UniqueID') %>% 
  inner_join(info, .) %>% 
  arrange(Instrument, Date) %>% 
  mutate(UniqueID = factor(UniqueID, levels = UniqueID))

tbl1 <- df_ident %>%
  pivot_longer(cols = pg:pr, names_to = 'Report.Type', values_to = 'Identity') %>% 
  mutate(label = as.integer(Identity),
         label = ifelse(str_count(label, '\\d') < 5, label, format(label, big.mark = ',')))
p_ident <- ggplot(tbl1) +
  facet_wrap(~Report.Type, scales = 'free_x') +
  aes(x = UniqueID, y = Identity) +
  geom_col(aes(fill = Instrument), color = NA) +
  geom_text(aes(label = label, color = Instrument, hjust = 0)) +
  labs(x = '# proteins', y = 'UniqueID') +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12), add = c(0, 0))) +
  ggsci::scale_color_cosmic() +
  ggsci::scale_fill_cosmic() +
  coord_flip() +
  theme_classic() +
  theme(text = element_text(size = 10))
ggsave('QC_35pool_identity.pdf', p_ident, width = 12, height = 6)


