rstudioapi::getActiveDocumentContext()$path
setwd('~/GitHub/TPHP/code/Compare_with_entrapment_lib/')
# source('../../source/my_fun.R')
library(tidyverse)
library(magrittr)
library(RColorBrewer)
library(pheatmap)
library(broom)


TPHP_HOME <- '//172.16.13.136/TPHP/'



fdp_est <- function(N.τ, N.ε, r) {
  FDP.combined <- N.ε * (1 + 1/r) / (N.τ + N.ε)
  return(FDP.combined)
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
    labs(x = col1, y = col2, subtitle = str_glue("{col1}: {col2} (x: y, n = {data_points})"))+
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
  if (!is.null(plot_title)){
    p <- p + labs(title = plot_title)
  }
  return(p)
}


# 1.Entrapment library estimate ----------
lib_tphp <- rio::import('//172.16.13.136/TPHP/library/TPHPlib_frag1025_swissprot_final.tsv')
lib_entrap1 <- rio::import('//172.16.13.136/share/members/jiangwenhao/20250401_GNHSF_lib_tims90minIGC_88/IGC_humanswiss_irt_contam_88_ddafile_NEW_rmone_20220102.tsv') %>%
  mutate(AverageExperimentalRetentionTime = NA) #%>% 
  #filter(!(PeptideSequence %in% lib_tphp$PeptideSequence))

lib_entrap2 <- rio::import('//172.16.13.136/share/members/jiangwenhao/20250327_GNHSF_lib_tims60minIGC/20250718_MPW_GNHS_IGC/library.tsv')
# setdiff(names(lib_tphp), names(lib_entrap)) # AverageExperimentalRetentionTime
length(unique(lib_tphp$PeptideSequence)) # 484,391
length(unique(lib_entrap1$PeptideSequence)) # 506,470
length(unique(lib_entrap2$PeptideSequence)) # 273,779
table(str_detect(unique(lib_entrap1$ProteinId), 'Human_swiss_sp')) # 1841 / (224427+1841)
table(str_detect(unique(lib_entrap2$ProteinId), 'Human_swiss_sp')) # 1202 / (66955+1202)
# human_swiss <- str_subset(unique(lib_entrap$ProteinId), '86ggHuman_swiss_sp') %>% 
#   str_extract_all('Human_swiss_sp\\|[A-Z0-9]+')
# table(sapply(human_swiss, length)) # 1-10
# human_swiss_pro <- str_remove(unlist(human_swiss), '^Human_swiss_sp\\|')
# length(human_swiss_pro) # 418
# setdiff(human_swiss_pro, lib_tphp$ProteinId) # 14 proteins

# add PrecursorId
lib_tphp %<>%
  mutate(PrecursorId = str_c(ModifiedPeptideSequence, PrecursorCharge))
lib_entrap1 %<>%
  mutate(PrecursorId = str_c(ModifiedPeptideSequence, PrecursorCharge))
lib_entrap2 %<>%
  mutate(PrecursorId = str_c(ModifiedPeptideSequence, PrecursorCharge))

# merge entrap1 and entrap2
lib_entrap <- lib_entrap2 %>%
  filter(!(PeptideSequence %in% lib_entrap1$PeptideSequence),
         !(ProteinId %in% lib_entrap1$ProteinId)) %>% 
  rbind(lib_entrap1, .)

# check spiked-iRT
irt_tphp <- lib_tphp %>%
  distinct(ProteinId, PrecursorId, ModifiedPeptideSequence, PeptideSequence) %>% 
  filter(str_detect(ProteinId, 'iRT'))

irt_entrap <- lib_entrap %>%
  distinct(ProteinId, PrecursorId, ModifiedPeptideSequence, PeptideSequence) %>% 
  filter(PeptideSequence %in% irt_tphp$PeptideSequence) # from entrap2

lib_entrap %<>% filter(!(PeptideSequence %in% irt_tphp$PeptideSequence))


intersect(lib_tphp$PrecursorId,
          lib_entrap$PrecursorId) %>% length() # 12681
intersect(lib_tphp$ModifiedPeptideSequence,
          lib_entrap$ModifiedPeptideSequence) %>% length() # 10386
intersect(lib_tphp$PeptideSequence,
          lib_entrap$PeptideSequence) %>% length() # 9506
pep_overlap <- intersect(lib_tphp$PeptideSequence, lib_entrap$PeptideSequence)
lib_overlap <- list(tphp = lib_tphp %>% filter(PeptideSequence %in% pep_overlap),
     entrap = lib_entrap %>% filter(PeptideSequence %in% pep_overlap)) %>% 
  plyr::ldply(.id = 'Lib')

rio::export(lib_overlap, 'TPHP_entrapment_lib_overlapped_peptides.xlsx')



# regression of iRT values
irt_overlap <- lib_overlap %>%
  distinct(Lib, ModifiedPeptideSequence, PrecursorCharge, NormalizedRetentionTime) %>% 
  pivot_wider(id_cols = c('ModifiedPeptideSequence', 'PrecursorCharge'),
              names_from = Lib, values_from = NormalizedRetentionTime) %>% 
  rename(iRT.Universal = tphp, iRT.Entrapment = entrap)

p_irt <- my_cor_smooth(irt_overlap, col1 = 'iRT.Universal', col2 = 'iRT.Entrapment')
p_irt2 <- p_irt + aes(color = factor(irt_overlap$PrecursorCharge)) +
  ggsci::scale_color_lancet(name = 'PrecursorCharge', alpha = 0.6)
ggsave('overlapped_peptide_irt_v1_20250718.pdf', p_irt, width = 5, height = 5)
ggsave('overlapped_peptide_irt_v2_20250718.pdf', p_irt2, width = 8, height = 8)


# 2.Target-decoy library prepare ----------
lib_d %>% distinct(ProteinId) %>% rio::export('tphpSwiss1025_mpwIGC_IGC88_library_proteinid_check.tsv')
lib_d %>% distinct(ProteinId) %>% filter(str_detect(ProteinId, '')) %>% rio::export('tphpSwiss1025_mpwIGC_IGC88_library_proteinid_check_remove.tsv')

# entrapment library filter
peps_t <- unique(lib_tphp$PeptideSequence) # 484,391 targets
peps_d_neq <- union(peps_t, stringi::stri_reverse(peps_t)) # decoys cannot be in these 968,767 entries
lib_d <- lib_entrap %>%
  filter(!(PeptideSequence %in% peps_d_neq)) %>% 
  filter(!str_detect(ProteinId, 'Human_swiss_sp|SWISS-PROT|MAX_'))


# stat
df_stat <- data.frame(
  pr.e = length(unique(lib_d$PrecursorId)),
  pr.t = length(unique(lib_tphp$PrecursorId)),
  pep.e = length(unique(lib_d$PeptideSequence)),
  pep.t = length(unique(lib_tphp$PeptideSequence)),
  pg.e = length(unique(lib_d$ProteinId)),
  pg.t = length(unique(lib_tphp$ProteinId))
)
df_stat %<>% mutate(
  r_pr = pr.e / pr.t, # 1.072
  r_pep = pep.e / pep.t, # 1.152
  r_pg = pg.e / pg.t # 16.714
)

# append entrapment to universal library
lib_final <- plyr::ldply(list(t = lib_tphp, e = lib_d), .id = 'Tag')
rio::export(lib_final, 'tphpSwiss1025_mpwIGC_IGC88_library.tsv')


# EOF ------------
# # 
# length(unique(lib_tphp$PeptideSequence)) # 484,391
# length(unique(lib_entrap$PeptideSequence)) # 497,286
list(iRT.overlap = irt_overlap,
     r.values = df_stat) %>%
  rio::export('source_data_mixed_lib_20250718.xlsx')

