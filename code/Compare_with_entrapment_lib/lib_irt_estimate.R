rstudioapi::getActiveDocumentContext()$path
setwd('~/GitHub/TPHP/code/Compare_with_entrapment_lib/')
source('../../source/my_fun.R')
library(tidyverse)
library(magrittr)
library(RColorBrewer)
library(pheatmap)
library(broom)


TPHP_HOME <- '//172.16.13.136/TPHP/'


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


lib_tphp %<>%
  mutate(PrecursorId = str_c(ModifiedPeptideSequence, PrecursorCharge))
lib_entrap  %<>%
  mutate(PrecursorId = str_c(ModifiedPeptideSequence, PrecursorCharge))

lib_tphp$ProteinId %>% intersect(lib_entrap$precurid)




