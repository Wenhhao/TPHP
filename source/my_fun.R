
# ------- transform color name from 0-255 decimal-format to hexadecimal-format -------
hex_col <- function(ls){
  for(e in ls){
    if((FALSE %in% unique(e %in% 0:255)) | (length(e) != 3)){
      stop('wrong decimal-format (0-255)')
    }
  }
  hex <- lapply(ls, function(vec){
    stringr::str_c('#', stringr::str_c(as.hexmode(vec), collapse = ''))
  })
  return(unlist(hex))
  
}

# ------ get abbreviation of sample types, return input name which not matched -------
# df_abbr$Entrie should be total lowercase !!!!
get_abbr <- function(df, sample_type, df_abbr = NULL){
  if(is.null(df_abbr)){
    df_abbr <- read.delim('//172.16.13.114/share/members/jiangwenhao/TPHP/input/sample_types_abbr_20220722.txt', stringsAsFactors = F, check.names = F, na.strings = '')
    # df_abbr[is.na(df_abbr)] <- 'NA'
  }
  h <- hash::hash(keys = df_abbr$'Entire',
            values = df_abbr$'Abbr')
  df[, sample_type] <- unlist(lapply(df[, sample_type], function(e){
    if(is.null(h[[tolower(e)]])){
      return(e)
    }else{
      return(h[[tolower(e)]])
    }
  }))
  return(df)
}

# ------------- read all sheets from excel table(s) ----------------------
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


# peptide matrix to protein -----------------------------------------------
easyp2p <- function(pep_data, roundn = NULL){
  # peptide matrix (row: peptides; column: files)
  #   - peptide sequence as the first column
  #   - protein id as the second column
  #   - data should be in raw scale
  
  # Step 1
  t1 <- proc.time()
  cat("1) Data preparation: \n"); print(proc.time() - t1)
  
  pep_data[is.infinite(as.matrix(pep_data))] <- NA # assign infinity as NA
  pep_data[is.na(pep_data)] <- NA
  pep_data <- pep_data[complete.cases(pep_data[, 1]), ] # remove NA peptides
  pep_data <- pep_data[!grepl("^1/CON", pep_data[, 2], fixed = F), ] # remove contaminants
  pep_data[pep_data == 0] <- NA # assign zeros as NA
  
  
  # Step 2
  cat("2) log2 transformation: \n"); print(proc.time() - t1)
  pep_data_log2 <- log2(pep_data[, 3:ncol(pep_data), drop = F]) %>% tibble::add_column(., prot = pep_data$prot, .before = 1)
  rownames(pep_data_log2) <- pep_data[, 1]
  
  
  # Step 3
  cat("3) Quantile normalization: \n"); print(proc.time() - t1)
  pep_data_log2_qn <- preprocessCore::normalize.quantiles(as.matrix(pep_data_log2))
  colnames(pep_data_log2_qn) <- colnames(pep_data_log2)
  rownames(pep_data_log2_qn) <- rownames(pep_data_log2)
  
  data_tech_rep <- cbind(pep_data[, 1:2], pep_data_log2_qn)
  #is.null(tech_rep_f)
  #is.null(batchf)
  
  rm(list = ls()[grep('^pep_data', ls())])
  
  
  # Step 4
  cat("4) Arrangement: \n"); print(proc.time() - t1)
  data <- data_tech_rep
  colnames(data)[1:2] <- c("tg", "prot")
  
  n <- ncol(data)
  pep2 <- apply(data[, -c(1, 2), drop = F], 1, function(x) { # log2 then mean of all files
    NAs <- length(which(is.na(x)))
    meanexpr1 <- sum(as.numeric(x), na.rm = TRUE) / (n - NAs)
    meanexpr2 <- sum(as.numeric(x), na.rm = TRUE) / n
    d <- c(NAs, meanexpr1, meanexpr2)
    return(d)
  })
  pep2 <- t(pep2)
  colnames(pep2) <- c("NAs", "meanexpr1", "meanexpr2")
  pep_expr <- cbind(data[, 1], pep2, data[, c(-1)])
  
  #order by pg ,#NA,intesity
  pep_order <- pep_expr[order(pep_expr[, 5], pep_expr[, 2], -pep_expr[, 3]), ]
  colnames(pep_order)[1] <- "tg"
  pep_order2 <- pep_order[, c(-2, -3, -4)]
  
  rm(list = ls()[grep('^data', ls())])
  rm(list = c('pep2', 'pep_expr', 'pep_order'))
  gc() # Garbage Collection
  
  # Step 5
  cat("5) Top 3 peptides selection: \n"); print(proc.time() - t1)
  pep_order2_dup1 <- pep_order2[duplicated(pep_order2$prot), ]
  pep_order2_dup2 <- pep_order2_dup1[duplicated(pep_order2_dup1$prot), ]
  pep_order2_dup3 <- pep_order2_dup2[duplicated(pep_order2_dup2$prot), ]
  pep_order2_top3 <- pep_order2[!(rownames(pep_order2) %in% rownames(pep_order2_dup3)), ]
  
  rm(list = ls()[grep('^pep_order2_dup', ls())])
  
  
  # Step 6
  cat("6) Protein matrix generation: \n"); print(proc.time() - t1)
  pep_order2_top3 <- pep_order2_top3[c("prot", "tg", colnames(pep_order2_top3)[3:ncol(pep_order2_top3)])]
  pep_order2_top3[pep_order2_top3 == 0] <- NA
  
  lr_top3 <- "top3"
  if(lr_top3 == "top3"){ # mean of top 3
    top3_mean <- plyr::ddply(pep_order2_top3,
                             .variables = "prot",
                             .fun = function(df_sub){
                               mean_ls <- colMeans(df_sub[, -c(1, 2), drop = F], na.rm = T)
                               if(!is.null(roundn)){
                                 mean_ls <- round(mean_ls, roundn)
                               }
                               return(mean_ls)
                             })
    df_prot <- top3_mean
    #readr::write_csv(top3_mean, 'prot_matrix_top3.csv', na = '')
    
  }else{ # LR
    prot_matrix <- pep2prot(pep_order2_top3)
    prot_matrix <- prot_matrix[, -2]
    
    df_prot <- prot_matrix
    #readr::write_csv(prot_matrix, 'prot_matrix_lr.csv', na = '')
  }
  
  return(df_prot)
  
}





# crawl protein name from uniprot.org -------------------------------------
my_uniprot2prot <- function(vec_uni){
  vec_prot <- lapply(vec_uni, function(uniprotid){
    # my_html <- rvest::read_html(stringr::str_glue('https://www.uniprot.org/uniprotkb/{uniprotid}/entry'))
    my_html <- rvest::read_html(stringr::str_glue('https://rest.uniprot.org/genecentric/{uniprotid}')) # 2022-08-08 update; get json
    my_json <- my_html %>%
      rvest::html_text() %>%
      jsonlite::fromJSON()
    
    prot <- my_json$canonicalProtein$proteinName
    return(prot)
  }) %>% unlist
  return(stringr::str_c(vec_uni, vec_prot, sep = '_'))
}


# scatter plot of protein abundance ---------------------------------------
my_plot <- function(df){
  df %<>%
    add_column(label = stringr::str_c(df$organ, ' (', df$cancer_type, ')'), .before = 1) %>%
    dplyr::mutate(cancer_type = factor(cancer_type, levels = c('N', 'Adj', 'C'))) %>%
    arrange(organ, cancer_type) %>%
    dplyr::mutate(cancer_type = as.character(cancer_type))
  
  #remove organ with all subtypes >= 90% NA ratio
  na_ratio <- function(x) sum(is.na(x)) / length(x)
  organ_selected <- df %>% dplyr::group_by(organ, cancer_type) %>%
    summarise_at(vars(2), list(na_ratio = 'na_ratio')) %>%
    ungroup() %>%
    filter(na_ratio < 0.9) %>%
    pull(organ) %>%
    unique()
  
  df %<>% filter(organ %in% organ_selected)
  
  rlt <- list()
  #df[is.na(df)] <- 6.759341
  df_fillrate <- df
  my_uni <- colnames(df)[ncol(df)]
  my_title <- my_uniprot2prot(my_uni)
  colnames(df_fillrate)[ncol(df_fillrate)] <- colnames(df)[ncol(df)] <- c('intensity')
  
  
  df <- na.omit(df)
  df_fillrate <- df_fillrate[df_fillrate$organ %in% unique(df$organ), c(2, 3, 4)]
  
  # color setting
  df$intensity <- log10(2 ^ df$intensity) # log2 -> log10
  # pseudo_ht_color <- c('#F01025', '#EAD41A', '#1E72BA')
  pseudo_ht_color <- c('#F01025', '#EAD41A', '#1E72BA')
  names(pseudo_ht_color) <- c('C', 'Adj', 'N')
  my_fills <- pseudo_ht_color[df$cancer_type]
  names(my_fills) <- df$cancer_type
  
  # text plotting
  my_texts_pos <- 1:length(unique(df$organ)) - 0.5
  my_vlines <- 1:length(unique(df$organ))
  x <- c()
  for(i in 1:length(unique(df$organ))){
    df_tmp <- df[df$organ == unique(df$organ)[i],]
    x_tmp <- seq(from = i-1, to = i, length.out = nrow(df_tmp)+2)[2:(nrow(df_tmp)+1)]
    x <- append(x, x_tmp)
  }
  df$x <- x
  
  # barplot preparation
  df_bar <- lapply(unique(df_fillrate$organ), function(e){
    df_tmp <- df_fillrate[df_fillrate$organ == e, c('cancer_type', 'intensity')]
    df_ttmp <- lapply(unique(df_tmp$cancer_type), function(ee){
      vec_ttmp <- df_tmp$intensity[df_tmp$cancer_type == ee]
      return(sum(!is.na(vec_ttmp)) / length(vec_ttmp))
    }) %>% as.data.frame(stringsAsFactor = F)
    colnames(df_ttmp) <- unique(df_tmp$cancer_type)
    rownames(df_ttmp) <- e
    return(df_ttmp)
  }) %>% do.call(plyr::rbind.fill, .)
  df_bar$organ <- unique(df_fillrate$organ)
  df_bar %<>% reshape2::melt()
  colnames(df_bar)[2:3] <- c('cancer_type', 'fillingvalue')
  df_bar$cancer_type <- factor(df_bar$cancer_type, levels = c('C', 'Adj', 'N'), ordered = T)
  rlt$table <- df_bar %<>% dplyr::arrange(organ, cancer_type)
  df_bar %<>% na.omit
  df_bar$label <- stringr::str_c(df_bar$organ, df_bar$cancer_type, sep = '_')
  pseudo_ht_bar <- seq(from = 0, to = length(unique(df_bar$organ)), by = 1/4) %>% .[. %% 1 != 0]
  names(pseudo_ht_bar) <- stringr::str_c(lapply(unique(df_bar$organ), function(e) rep(e, 3)) %>% unlist, c('N', 'Adj', 'C'), sep = '_')
  
  df_bar$x <- pseudo_ht_bar[df_bar$label]
  df_bar <- df_bar[df_bar$fillingvalue != 0, ]
  nmax <- max(df$intensity); nmin <- min(df$intensity)
  #df_bar$fillingvalue_lt <- my_linear_trans(c(df_bar$fillingvalue, 0.5), nmax, nmin) %>% .[1:(length(.) - 1)]
  df_bar$fillingvalue_lt <- df_bar$fillingvalue * (nmax - nmin) + nmin
  # double coordinates
  fillingvalue_x <- max(df$x) + 0.5
  fillingvalue_scale <- seq(nmin, nmax, length.out = 6)
  fillingvalue_str <- stringr::str_c(seq(0, 100, length.out = 6), '%', sep = '')
  fillingtitle <- '1-missing rate'
  fillingtitle_pos <- mean(c(nmin, nmax))
  #hline0.5_lt <- my_linear_trans(c(df_bar$fillingvalue, 0.5), nmax, nmin) %>% .[length(.)]
  
  # figure
  tmp <- stringr::str_c("+annotate('text', x = ", my_texts_pos, ", y = max(df$intensity) * 1.02, label = '", unique(sort(df$organ)), "', color = '#000000', size = 4)", collapse = ' ') # labels of sample type
  p <- ggplot()+
    geom_bar(data = df_bar, mapping = aes(x = x, y = fillingvalue_lt, fill = cancer_type), fill = pseudo_ht_color[df_bar$cancer_type], stat = 'identity', alpha = 0.2, width = 0.25)+
    geom_bar(data = df_bar, mapping = aes(x = x, y = nmin, fill = cancer_type), fill = '#FFFFFF', stat = 'identity', width = 0.25)+# block non-missing rate lower than 0
    geom_point(data = df, mapping = aes(x = x, y = intensity), fill = my_fills, color = '#000000', shape = 21, size = 5)+
    geom_vline(xintercept = my_vlines, color = '#000000', linetype = 'dashed', size = 0.5)+
    #geom_hline(yintercept = hline0.5_lt, color = pseudo_ht_color['C'], linetype = 'dashed', size = 0.5, alpha = 0.5)+
    #annotate('text', x = 0.2, y = hline0.5_lt * 1.01, label = 'Missing value < 50%', color = pseudo_ht_color['C'], size = 4, alpha = 0.5, hjust = 0)+
    annotate('text', x = fillingvalue_x * 1.06, y = fillingtitle_pos, label = fillingtitle, size = 4.5, hjust = 0.5, vjust = 0.5, angle = 270)+# mock y-axis title
    annotate('text', x = fillingvalue_x * 1.01, y = fillingvalue_scale, label = fillingvalue_str, size = 4, hjust = 0, vjust = 0.5)+# mock y-axis text
    labs(
      y = 'log10 Intensity', title = my_title
    )+
    geom_vline(xintercept = c(placeholder_x = fillingvalue_x * 1.1), color = '#FFFFFF')+# placeholder for x-axis
    coord_cartesian(ylim = c(nmin, nmax * 1.02))+
    scale_x_continuous(expand = c(0, 0))+
    theme(panel.grid = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.text = element_text(size = 12,color = 'black'),
          axis.line = element_line(color = 'black'),
          axis.line.x = element_blank(),
          plot.subtitle = element_text(size = 30, hjust = 0, color = 'black')
    )+
    theme(legend.text = element_text(size = 12, color = 'black'), legend.position = 'top',
          legend.title = element_text(size = 15, color = 'black'))+
    theme(axis.title.x = element_blank(),
          #axis.text.x = element_text(size = 12, vjust = 0.5 ,color = 'black', angle = 90),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()
    )+
    theme(axis.title.y = element_text(size = 12, hjust = 0.5, color = 'black', angle = 90),
          axis.text.y = element_text(size = 10,color = 'black', angle = 0),
          axis.ticks.y = element_blank()
    )
  rlt$p <- eval(parse(text = stringr::str_c('p', tmp)))
  return(rlt)
}

