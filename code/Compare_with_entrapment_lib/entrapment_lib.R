rstudioapi::getActiveDocumentContext()$path
library(tidyverse)
library(magrittr)

lib_tphp <- rio::import('//172.16.13.136/TPHP/library/TPHPlib_frag1025_swissprot_final.tsv')
lib_entrap <- rio::import('//172.16.13.136/share/members/jiangwenhao/20250401_GNHSF_lib_tims90minIGC_88/IGC_humanswiss_irt_contam_88_ddafile_NEW_rmone_20220102.tsv')
setdiff(names(lib_tphp), names(lib_entrap)) # AverageExperimentalRetentionTime


length(unique(lib_tphp$PeptideSequence)) # 484,391
length(unique(lib_entrap$PeptideSequence)) # 506,470
length(intersect(lib_tphp$PeptideSequence, lib_entrap$PeptideSequence)) # 9184



# bind
lib_entrap$AverageExperimentalRetentionTime <- NA
lib_bind <- lib_entrap %>% 
  filter(!(PeptideSequence %in% lib_tphp$PeptideSequence)) %>% 
  rbind(lib_tphp, .)
rio::export(lib_bind, 'entrapment_lib.tsv')


