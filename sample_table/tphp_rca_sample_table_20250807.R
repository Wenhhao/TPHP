rstudioapi::getActiveDocumentContext()$path
setwd('~/GitHub/TPHP/sample_table')

library(magrittr)
library(tidyverse)

read_excel_allsheets <- function(filename, tibble = FALSE) {
  # I prefer straight data.frames
  # but if you like tidyverse tibbles (the default with read_excel)
  # then just pass tibble = TRUE
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

read_excel_allsheets_from_list <- function(flist, tibble = FALSE) {
  lapply(flist, read_excel_allsheets) %>% do.call(c, .) %>% return
}

date_rescue <- function(x){
  if (length(x) != 1)
    print('Input should be only one element.')
  flag <- tryCatch(as.numeric(x), warning = function(w) 'coercion')
  if (flag != 'coercion'){
    # if coercion occurred, do not convert the date
    x <- flag
    # To rescue the wrong Excel date
    if(str_length(x) >= 9) { # covering [1973-03-03 09:46:40 UTC, +∞)
      # Likely a Unix timestamp
      x <- as.POSIXct((as.numeric(x)), origin = "1970-01-01", tz = "UTC")
    } else if (str_length(x) == 5){ # covering [1927-05-18, 2173-10-13]
      # Likely an Excel date serial number
      x <- as.Date(as.numeric(x), origin = "1899-12-30")
    }
  }
  return(as.character(x))
}

# 1.patient table -------
tbl_ls <- read_excel_allsheets('2025.06.01 西湖样本全部列表 共计33个类型_称量表_edit_append_2nd.xlsx')

# lapply(tbl_ls, function(X) colnames(X))
# table(sapply(tbl_ls, function(X) ncol(X)))
# tbl_ls[sapply(tbl_ls, function(X) ncol(X)) == 8]
# tbl_ls[sapply(tbl_ls, function(X) ncol(X)) == 9]
# tbl_ls[sapply(tbl_ls, function(X) ncol(X)) == 10]
# tbl_ls[sapply(tbl_ls, function(X) ncol(X)) == 10][[1]] %>%
#   rename(`病理诊断` = `...7`, `病理诊断2` = `...8`) %>% 
#   unite(col = '病理诊断', 病理诊断, 病理诊断2, sep = '', na.rm = T, remove = T)


tbl_ls_tidy <- lapply(tbl_ls, function(tbl){
  # revise
  if ('Column1' %in% colnames(tbl))
    tbl %<>% rename(
      '病理号' = Column1,
      '住院号' = Column2,
      '姓名' = Column3,
      '性别' = Column4,
      '年龄' = Column5,
      '取材日期' = Column6,
      '病理诊断' = Column7,
      '肿物' = Column8,
      '癌旁' = Column9
    )
  if ('均无癌旁' %in% colnames(tbl))
    tbl %<>% rename(`癌旁` = `均无癌旁`)
  if ('报告日期' %in% colnames(tbl))
    tbl %<>% rename(`取材日期` = `报告日期`)
  if ('日期' %in% colnames(tbl))
    tbl %<>% rename(`取材日期` = `日期`)
  if ('...7' %in% colnames(tbl))
    tbl %<>% rename(`病理诊断` = `...7`)
  if ('...8' %in% colnames(tbl))
    tbl %<>% rename(`病理诊断2` = `...8`) %>% 
    unite(col = '病理诊断', 病理诊断, 病理诊断2, sep = '', na.rm = T, remove = T)
  # if ('...10' %in% colnames(tbl))
  #   tbl %<>% rename(Note = `...10`)
  if (sum(c('男', '女') %in% tbl$年龄) > 0)
    tbl %<>% rename(a = 性别, b = 年龄) %>% rename(性别 = b, 年龄 = a)
  
  # append
  if (!('病理诊断' %in% colnames(tbl)))
    tbl %<>% mutate(`病理诊断` = NA)
  # if (!('Note' %in% colnames(tbl)))
  #   tbl %<>% mutate(`Note` = NA)
  
  return(tbl)
})

Reduce(union, lapply(tbl_ls_tidy, colnames))
Reduce(intersect, lapply(tbl_ls_tidy, colnames))
unique(sapply(tbl_ls_tidy, ncol))

info <- plyr::ldply(tbl_ls_tidy, .id = '癌症类型')
info_tidy <- info %>%
  filter(!is.na(`病理号`)) %>% 
  relocate(病理诊断, .after = 取材日期) %>% 
  mutate(`年龄` = as.numeric(str_remove(`年龄`, '岁$'))) %>% 
  mutate(病理诊断 = ifelse(str_remove(取材日期, '^[\\d\\-]+') != '',
                           str_c(str_remove(取材日期, '^[\\d\\-]+'), 病理诊断),
                           病理诊断)) %>% 
  mutate(取材日期 = str_extract(取材日期, '^[\\d\\-]+'))

# rio::export(info_tidy, '2025.06.01patient_info.xlsx')

sum(sapply(tbl_ls, nrow)) # 33 cancer types
length(unique(info_tidy$姓名)) # 974 patients
nrow(info_tidy) # 991 patients
12*88+82 == 1138 # samples


# 2.sample table -------
info_smp <- read_excel_allsheets('sample_info_20250714_corrected.xlsx') %>% 
  plyr::ldply(.id = 'Batch')

tmp1 <- info_tidy %>% filter(病理号 %in% info_smp$病理号) %>%
  count(病理号) %>% filter(n > 1) %>% semi_join(info_tidy, .)
tmp2 <- info_tidy %>% filter(病理号 %in% info_smp$病理号) %>%
  count(病理号) %>% filter(n > 1) %>% semi_join(info_smp, .)
# rio::export(list(tmp1, tmp2), 'multiple_info_check.xlsx')

# after checking the FFPE samples
tmp1_filter <- tibble(
  `癌症类型` = c("横纹肌肉瘤", "输卵管癌", "输卵管癌", "大肠粘液腺癌", "横纹肌肉瘤", "输卵管癌"),
  `病理号` = as.character(c(2001948, 2004084, 2029513, 2210851, 2223877, 2233444)),
  `肿物` = c("6-11", "1-2", "1", "5-8", "4-10", "4-5")
)
tmp1_filtered <- tmp1 %>% semi_join(tmp1_filter)

info_tidy %<>% anti_join(tmp1) %>% rbind(tmp1_filtered)
info_smp <- info_tidy %>% inner_join(info_smp)

info_smp %<>% mutate(
  Cancer_type = setNames(c('T', 'NT'), c('C', 'N'))[Cancer_type],
  BatchID = str_c(Batch, '-', Batch_n), .before = Batch
)
info_smp %<>%
  mutate(Date.Repaired = sapply(info_smp$`取材日期`, date_rescue),
         .after = 取材日期)

rio::export(info_smp, 'TPHP_RCA_sample_info_20250807.xlsx')



