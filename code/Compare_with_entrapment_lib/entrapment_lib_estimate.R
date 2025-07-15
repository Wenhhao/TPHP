# Reference: https://doi.org/10.1038/s41592-025-02719-x
##  FDP.combined <- N.ε * (1 + 1/r) / (N.τ + N.ε)
rstudioapi::getActiveDocumentContext()$path
setwd('~/GitHub/TPHP/code/Compare_with_entrapment_lib/')
source('../../source/my_fun.R')
library(tidyverse)
library(magrittr)
library(RColorBrewer)
library(pheatmap)
library(broom)


TPHP_HOME <- '//172.16.13.136/TPHP/'
# 
# target_decoy_color <- ggsci::pal_d3()(2) %>% setNames(c('Target', 'Decoy'))
# library_color <- ggsci::pal_d3()(2) %>% setNames(c('TPHP', 'Entrapment'))
# cancer_color <- ggsci::pal_nejm()(5) %>% setNames(c('C_pool', 'N_pool', 'Normal', 'adjacent', 'carcinoma'))
# df_color <- rio::import('//172.16.13.136/share/members/jiangwenhao/TPHP/input/PUH_tissue_colorset_20230210.xlsx')
# tissue_color <- str_c('#', df_color$color[1:51]) %>% setNames(df_color$tissue[1:51])
# tissue_color %<>% append(c('FO' = '#23087C'))
# 
# format_custom <- function(x, cutoff = 0.01) {
#   format_custom_individual <- function(x, cutoff = 0.01){
#     if (x >= cutoff) {
#       sprintf('%.3f', x)
#     } else {
#       as.character(signif(x, 3))
#       # format(x, scientific = F)
#     }
#   }
#   sapply(x, format_custom_individual, cutoff = cutoff)
#   
# }
# my_cor_smooth <- function(df, col1, col2, cor_method = "pearson", ggpmisc_label = c("eq", "R2", "p"), plot_title = NULL){
#   data_points <- df %>% select(all_of(col1), all_of(col2)) %>% drop_na() %>% nrow()
#   df %<>% select(all_of(col1), all_of(col2)) %>% setNames(c('x', 'y'))
#   p <- ggplot(df, aes(x = x, y = y))+
#     geom_point(na.rm = T)+
#     ggpmisc::stat_poly_line(method = 'lm', se = TRUE, level = 0.95,
#                             linewidth = 1, color = '#1D439E') +
#     ggpmisc::stat_poly_eq(ggpmisc::use_label(ggpmisc_label),
#                           method = 'lm', coef.digits = 3, rr.digits = 3, p.digits = 3) +
#     # geom_smooth(formula = y ~ x, method = 'lm', color = "#00AFBB", na.rm = T)+
#     # ggpubr::stat_cor(method = cor_method, na.rm = T, cor.coef.name = 'rho', vjust = 3)+
#     labs(subtitle = str_glue("{col1}: {col2} (x: y, n = {data_points})"))+
#     theme_bw()+
#     theme(legend.text = element_text(size = 12, color = "black"),legend.position = 'top',
#           legend.title = element_text(size = 12, color="black") ,
#           panel.grid = element_blank(),
#           panel.background = element_blank(),
#           axis.line = element_line(colour = "black"),
#           axis.title=element_text(size = 12, hjust = 0.5, color="black"),
#           axis.text = element_text(size = 10,color = "black"),
#           plot.title = element_text(size = 12, face = "bold"),
#           plot.subtitle=element_text(size = 12, hjust = 0, color="black")
#     )
#   if (!is.na(plot_title)){
#     p <- p + labs(title = plot_title)
#   }
#   return(p)
# }

fdp_est <- function(N.τ, N.ε, r) {
  FDP.combined <- N.ε * (1 + 1/r) / (N.τ + N.ε)
  return(FDP.combined)
}


# # 0.Entrapment library prepare ----------
lib_tphp <- rio::import('//172.16.13.136/TPHP/library/TPHPlib_frag1025_swissprot_final.tsv')
lib_entrap <- rio::import('//172.16.13.136/share/members/jiangwenhao/20250401_GNHSF_lib_tims90minIGC_88/IGC_humanswiss_irt_contam_88_ddafile_NEW_rmone_20220102.tsv') %>%
  filter(!(PeptideSequence %in% lib_tphp$PeptideSequence)) %>% 
  mutate(AverageExperimentalRetentionTime = NA)
# setdiff(names(lib_tphp), names(lib_entrap)) # AverageExperimentalRetentionTime
# 
# 
length(unique(lib_tphp$PeptideSequence)) # 484,391
length(unique(lib_entrap$PeptideSequence)) # 497,286

X1 <- lib_tphp %>% distinct(ModifiedPeptideSequence, PrecursorCharge) %>% 
  unite(col = 'Precursor.Id', ModifiedPeptideSequence, PrecursorCharge, sep = '', remove = F)
nrow(X1) # 689,568
X2 <- lib_entrap %>% distinct(ModifiedPeptideSequence, PrecursorCharge) %>% 
  unite(col = 'Precursor.Id', ModifiedPeptideSequence, PrecursorCharge, sep = '', remove = F)
nrow(X2) # 669,502

length(unique(lib_tphp$ProteinId)) # 15,332
length(unique(lib_entrap$ProteinId)) # 224,670


r_pr <- length(unique(lib_entrap$PeptideSequence)) / length(unique(lib_tphp$PeptideSequence))
r_pep <- nrow(X2) / nrow(X1)
r_pg <- length(unique(lib_entrap$ProteinId)) / length(unique(lib_tphp$ProteinId))

tbl1 <- data.frame(name = c('r_pr', 'r_pep', 'r_pg'),
                   value = c(r_pr, r_pep, r_pg))

# 1.Data readin ------
tphp <- rio::import(file.path(TPHP_HOME, 'results/rlt_combine_entraplib_G3/report.pro_matrix_full_0005LibQval.tsv'), fill = TRUE)

df_info <- readxl::read_excel('../../input/20220708TPHP_1781file_info_edited_v3_2.xlsx') %>% as.data.frame()
df_info %<>% mutate(ID = sapply(FileName, function(x){
  str_extract(x, '^(\\w+?)\\d+.+_(\\d+)\\.d$', group = c(1, 2)) %>%
    str_c(collapse = '_')
}), .before = 1)

tphp %<>%
  mutate(Is.pr.Target = Precursor.Id %in% X1$Precursor.Id,
         Is.pg.Target = Protein.Group %in% unique(lib_tphp$ProteinId))

# 2.FDP estimate -------
## global pr FDP
tmp <- tphp %>% distinct(Precursor.Id, Is.pr.Target) %>% count(Is.pr.Target) %>% pull(n)
global.pr.FDP.combined <- fdp_est(N.τ = tmp[2], N.ε = tmp[1], r = r_pr)

## global pg FDP
tmp <- tphp %>% distinct(Protein.Group, Is.pg.Target) %>% count(Is.pg.Target) %>% pull(n)
global.pg.FDP.combined <- fdp_est(N.τ = tmp[2], N.ε = tmp[1], r = r_pg)

## specific pr FDP
runSpecific.pr.FDP.combined <- tphp %>% plyr::ddply('Run', function(dfsub){
  tmp <- dfsub %>% distinct(Protein.Group, Is.pr.Target) %>% count(Is.pr.Target) %>% pull(n)
  if(length(tmp) == 1) tmp <- c(0, tmp)
  data.frame(
    Run = dfsub$Run[1],
    runSpecific.pr.FDP.combined = fdp_est(N.τ = tmp[2], N.ε = tmp[1], r = r_pr)
  )
}) %>% arrange(desc(runSpecific.pr.FDP.combined))

## specific pg FDP
runSpecific.pg.FDP.combined <- tphp %>% plyr::ddply('Run', function(dfsub){
  tmp <- dfsub %>% distinct(Protein.Group, Is.pg.Target) %>% count(Is.pg.Target) %>% pull(n)
  if(length(tmp) == 1) tmp <- c(0, tmp)
  data.frame(
    Run = dfsub$Run[1],
    runSpecific.pg.FDP.combined = fdp_est(N.τ = tmp[2], N.ε = tmp[1], r = r_pg)
  )
}) %>% arrange(desc(runSpecific.pg.FDP.combined))

# identity
df_identity <- tphp %>% plyr::ddply('Run', function(dfsub){
    npg <- dfsub %>%
      distinct(Protein.Group, Is.pg.Target) %>%
      nrow()
    npr <- dfsub %>%
      distinct(Precursor.Id, Is.pr.Target) %>%
      nrow()
    data.frame(`# precursors` = npr, `# protein groups` = npg, check.names = F)
})

tbl2 <- data.frame(
  global.pr.FDP.combined = global.pr.FDP.combined,
  global.pg.FDP.combined = global.pg.FDP.combined
)

dfbox <- rbind(runSpecific.pg.FDP.combined %>% setNames(c('Run', 'FDP')) %>% mutate(Type = 'pg'),
               runSpecific.pr.FDP.combined %>% setNames(c('Run', 'FDP')) %>% mutate(Type = 'pr'))
p_runFDP <- ggplot(dfbox) +
  aes(x = Type, y = FDP, fill = Type) +
  geom_boxplot(width = 0.6, position = position_dodge()) +
  # geom_text(data = dfbox %>% filter(FDP > 0.01),
  #           aes(label = round(FDP, 4)), hjust = -0.5) +
  labs(x = str_glue("DIA runs\n(n={length(unique(dfbox$Run))})"), y = 'Run-specific FDP') +
  scale_y_continuous(limits = c(0, NA), expand = c(0, 0, 0.05, 0)) +
  # scale_fill_manual(values = target_decoy_color) +
  ggsci::scale_fill_npg() +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  theme(text = element_text(size = 12))
ggsave('entrapment_estimate_combined_runspecific.pdf', p_runFDP, width = 5, height = 4)


bad_runs <- runSpecific.pr.FDP.combined %>%
  filter(runSpecific.pr.FDP.combined > 0.01) %>% 
  pull(Run)
bad_runs_info <- df_info %>%
  filter(str_detect(FileName, str_c(bad_runs, collapse = '|')))

tbl3 <- df_identity %>%
  inner_join(runSpecific.pr.FDP.combined) %>% 
  inner_join(runSpecific.pg.FDP.combined)

# OUTPUT -------
list(r.values = tbl1,
     global.FDP = tbl2,
     # run.specific.pr.FDP = runSpecific.pr.FDP.combined,
     # run.specific.pg.FDP = runSpecific.pg.FDP.combined,
     identity = tbl3,
     bad.runs = bad_runs_info) %>% 
  rio::export('source_data.xlsx')
