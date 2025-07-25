# General functions ----------
get_outliers <- function(vec, coef = 1.5){
  # outliers based on Q1 and Q3
  stats <- quantile(vec, na.rm = T)
  iqr <- diff(stats[c(2, 4)])
  ret <- c(stats[2] - coef * iqr, stats[4] + coef * iqr)
  return(ret)
}

hex_col <- function(ls){
  # transform color name from 0-255 decimal-format to hexadecimal-format
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

get_random_colors <- function(n, from_color_palette = brewer.pal(8, "Set2")){
  # get random colors
  my_colors <- sample(from_color_palette, 1)
  i <- 1
  while(length(my_colors) < n){
    color <- sample(from_color_palette, 1)
    if(color != my_colors[i]){
      my_colors[i+1] <- color
      i <- i + 1
    }
  }
  return(my_colors)
}

info_of_pm <- function(pm){
  info <- data.frame(Filename = colnames(pm))
  info$instrument <- as.character(stringr::str_extract(colnames(pm), '^[A-Z]+'))
  info$date <- as.character(stringr::str_extract(colnames(pm), '[0-9]+')) %>% factor(., levels = unique(sort(.)), ordered = T)
  info$year <- stringr::str_sub(info$date, end = 4)
  info$month <- stringr::str_sub(info$date, end = -3)
  
  info$batch <- NA
  info$batch[grep("M202006|K202006|K202007|M202007|K202008|M202008|N202008", info$Filename)] <- 'b1'
  info$batch[grep("N202011", info$Filename)] <- 'b2'
  info$batch[grep("M202101", info$Filename)] <- 'b3'
  info$batch[grep("N202104", info$Filename)] <- 'b4'
  info$batch[grep("N202105|N202106", info$Filename)] <- 'b5'
  info$batch[is.na(info$batch)] <- 'b6'
  
  info$trans <- NA
  info$trans[grep("b1|b2", info$batch)] <- "T1"
  info$trans[grep("b3|b4", info$batch)] <- "T2"
  info$trans[grep("b5", info$batch)] <- "T3"
  info$trans[grep("b6", info$batch)] <- "T4"
  info$trans[grep("^N202203", info$Filename)] <- "T5"
  
  info$new <- ifelse(grepl('^2022', info$month), 'new', 'old')
  return(info)
}

read_excel_allsheets <- function(filename, tibble = FALSE) {
  # Read all sheets from excel table(s) --
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


easyp2p <- function(pep_data, roundn = NULL){
  # Summarise peptide matrix to protein --
  # @peptide matrix (row: peptides; column: files)
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





my_uniprot2prot <- function(vec_uni){
  # crawl protein name from uniprot.org
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


# Project-specific functions -------------
get_abbr <- function(df, sample_type, df_abbr = NULL){
  # get abbreviation of sample types, return input name which not matched
  # df_abbr$Entrie should be total lowercase !!!!
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

my_plot <- function(df){
  # Scatter plot of protein abundance
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

# NAguideR methods -------------------
nafunctions<-function(x,method="zero"){
  # cite: Shisheng Wang, Wenxue Li, Liqiang Hu, Jingqiu Cheng, Hao Yang, Yansheng Liu, NAguideR: performing and prioritizing missing value imputations for consistent bottom-up proteomic analyses, Nucleic Acids Research, gkaa498, https://doi.org/10.1093/nar/gkaa498.
  # Code from GitHub: Commit 15ec862 wangshisheng authored on Aug 20, 2021
  df<-df1<-as.data.frame(x)
  method<-tolower(method)
  if(method=="zero"){
    df[is.na(df)]<-0
  }
  else if(method=="minimum"){
    df[is.na(df)]<-min(df1,na.rm = TRUE)
  }
  else if(method=="colmedian"){
    library(e1071)
    df<-impute(df1,what ="median")
  }
  else if(method=="rowmedian"){
    library(e1071)
    dfx<-impute(t(df1),what ="median")
    df<-t(dfx)
  }
  else if(method=="knnmethod"){
    library(impute)
    data_zero1<-impute.knn(as.matrix(df1),k = 10, rowmax = 1, colmax = 1)#rowmax = 0.9, colmax = 0.9
    df<-data_zero1$data
  }
  else if(method=="seqknn"){
    library(SeqKnn)
    df <- SeqKNN(df1,k = 10)
  }
  else if(method=="bpca"){
    library(pcaMethods)
    data_zero1<-pcaMethods::pca(as.matrix(df1), nPcs = ncol(df1)-1, method = "bpca", maxSteps =100)
    df<-completeObs(data_zero1)
  }
  else if(method=="svdmethod"){
    library(pcaMethods)
    data_zero1<-pcaMethods::pca(as.matrix(df1), nPcs = ncol(df1)-1, method = "svdImpute")
    df<-completeObs(data_zero1)
  }
  else if(method=="lls"){
    library(pcaMethods)
    data_zero1<-llsImpute(t(df1), k = 10)
    df<-t(completeObs(data_zero1))
  }
  else if(method=="mle"){
    library(norm)
    xxm<-as.matrix(df1)
    ss <- norm::prelim.norm(xxm)
    thx <- norm::em.norm(ss)
    norm::rngseed(123)
    df <- norm::imp.norm(ss, thx, xxm)
  }
  else if(method=="qrilc"){
    library(imputeLCMD)
    xxm<-t(df1)
    data_zero1 <- imputeLCMD::impute.QRILC(xxm, tune.sigma = 1)[[1]]
    df<-t(data_zero1)
  }
  else if(method=="mindet"){
    library(imputeLCMD)
    xxm<-as.matrix(df1)
    df <- imputeLCMD::impute.MinDet(xxm, q = 0.01)
  }
  else if(method=="minprob"){
    library(imputeLCMD)
    xxm<-as.matrix(df1)
    df <- imputeLCMD::impute.MinProb(xxm, q = 0.01, tune.sigma = 1)
  }
  else if(method=="irm"){
    library(VIM)
    df <- irmi(df1, trace = TRUE,imp_var=FALSE)
    rownames(df)<-rownames(df1)
  }
  else if(method=="impseq"){
    library(rrcovNA)
    df <- impSeq(df1)
  }
  else if(method=="impseqrob"){
    library(rrcovNA)
    data_zero1 <- impSeqRob(df1, alpha=0.9)
    df<-data_zero1$x
  }
  else if(method=="mice-norm"){
    library(mice)
    minum<-5
    datareadmi<-mice(df1,m=minum,seed = 1234, method ="norm")
    newdatareadmi<-0
    for (i in 1:minum) {
      newdatareadmi<-complete(datareadmi,action = i)+newdatareadmi
    }
    df<-newdatareadmi/minum
    rownames(df)<-rownames(df1)
  }
  else if(method=="mice-cart"){
    library(mice)
    minum<-5
    datareadmi<-mice(df1,m=minum,seed = 1234, method ="cart")
    newdatareadmi<-0
    for (i in 1:minum) {
      newdatareadmi<-complete(datareadmi,action = i)+newdatareadmi
    }
    df<-newdatareadmi/minum
    rownames(df)<-rownames(df1)
  }
  else if(method=="trknn"){
    # source('Trunc_KNN/Imput_funcs.r') # only call the minimal requirements
    {
      ##################################################################################
      #### MLE for the Truncated Normal
      #### Creating a Function that Returns the Log Likelihood, Gradient and
      #### Hessian Functions
      ##################################################################################
      
      ## data = numeric vector
      ## t    = truncation limits
      mklhood <- function(data, t, ...) {
        
        data <- na.omit(data)
        n <- length(data)
        t <- sort(t)
        
        psi<-function(y, mu, sigma){
          exp(-(y-mu)^2/(2*sigma^2))/(sigma*sqrt(2*pi))
        }
        
        psi.mu<-function(y,mu,sigma){
          exp(-(y-mu)^2/(2*sigma^2)) * ((y-mu)/(sigma^3*sqrt(2*pi)))
        }
        
        psi.sigma<-function(y,mu,sigma){
          exp(-(y-mu)^2/(2*sigma^2)) *
            (((y-mu)^2)/(sigma^4*sqrt(2*pi)) - 1/(sigma^2*sqrt(2*pi)))
        }
        
        psi2.mu<-function(y,mu,sigma){
          exp(-(y - mu)^2/(2*sigma^2)) *
            (((y - mu)^2)/(sigma^5*sqrt(2*pi))-1/(sigma^3*sqrt(2*pi)))
        }
        
        psi2.sigma<-function(y,mu,sigma){
          exp(-(y-mu)^2/(2*sigma^2)) *
            ((2)/(sigma^3*sqrt(2*pi)) - (5*(y-mu))/(sigma^5*sqrt(2*pi)) +
               ((y-mu)^4)/(sigma^7*sqrt(2*pi)))
        }
        
        psi12.musig<-function(y,mu,sigma){
          exp(-(y-mu)^2/(2*sigma^2)) *
            (((y-mu)^3)/(sigma^6*sqrt(2*pi)) - (3*(y-mu))/(sigma^4*sqrt(2*pi)))
        }
        
        ll.tnorm2<-function(p){
          out <- (-n*log(pnorm(t[2],p[1],p[2])-pnorm(t[1],p[1],p[2]))) -
            (n*log(sqrt(2*pi*p[2]^2))) - (sum((data-p[1])^2)/(2*p[2]^2))
          -1*out
        }
        
        grad.tnorm<-function(p){
          g1 <- (-n*(integrate(psi.mu,t[1],t[2],mu=p[1],sigma=p[2], stop.on.error = FALSE)$value) /
                   (pnorm(max(t),p[1],p[2])-pnorm(min(t),p[1],p[2]))) - ((n*p[1]-sum(data))/p[2]^2)
          g2 <- (-n*(integrate(psi.sigma,t[1],t[2],mu=p[1],sigma=p[2], stop.on.error = FALSE)$value) /
                   (pnorm(max(t),p[1],p[2])-pnorm(min(t),p[1],p[2]))) - ((n)/(p[2])) + ((sum((data-p[1])^2))/(p[2]^3))
          out <- c(g1,g2)
          return(out)
        }
        
        hessian.tnorm<-function(p){
          
          h1<- -n*(integrate(psi,t[1],t[2],mu=p[1],sigma=p[2], stop.on.error = FALSE)$value *
                     integrate(psi2.mu,t[1],t[2],mu=p[1],sigma=p[2], stop.on.error = FALSE)$value -
                     integrate(psi.mu,t[1],t[2],mu=p[1],sigma=p[2], stop.on.error = FALSE)$value^2) /
            (integrate(psi,t[1],t[2],mu=p[1],sigma=p[2], stop.on.error = FALSE)$value^2) -
            n/(p[2]^2)
          
          h3<- -n*(integrate(psi,t[1],t[2],mu=p[1],sigma=p[2], stop.on.error = FALSE)$value *
                     integrate(psi12.musig,t[1],t[2],mu=p[1],sigma=p[2], stop.on.error = FALSE)$value -
                     integrate(psi.mu,t[1],t[2],mu=p[1],sigma=p[2], stop.on.error = FALSE)$value *
                     integrate(psi.sigma,t[1],t[2],mu=p[1],sigma=p[2], stop.on.error = FALSE)$value) /
            (integrate(psi,t[1],t[2],mu=p[1],sigma=p[2], stop.on.error = FALSE)$value^2) +
            (2*(n*p[1]-sum(data)))/(p[2]^3)
          
          h2<- -n*(integrate(psi,t[1],t[2],mu=p[1],sigma=p[2], stop.on.error = FALSE)$value *
                     integrate(psi2.sigma,t[1],t[2],mu=p[1],sigma=p[2], stop.on.error = FALSE)$value -
                     integrate(psi.sigma,t[1],t[2],mu=p[1],sigma=p[2], stop.on.error = FALSE)$value^2) /
            (integrate(psi,t[1],t[2],mu=p[1],sigma=p[2], stop.on.error = FALSE)$value^2) +
            (n)/(p[2]^2)-(3*sum((data-p[1])^2))/(p[2]^4)
          
          H<-matrix(0,nrow=2,ncol=2)
          H[1,1]<-h1
          H[2,2]<-h2
          H[1,2]<-H[2,1]<-h3
          return(H)
        }
        
        
        return(list(ll.tnorm2 = ll.tnorm2, grad.tnorm = grad.tnorm, hessian.tnorm = hessian.tnorm))
      }
      ##################################################################################
      ###### Newton Raphson Function
      ###### This takes in the Objects Returned from mklhood Function above
      ##################################################################################
      
      NewtonRaphsonLike <- function(lhood, p, tol = 1e-07, maxit = 100) {
        
        cscore <- lhood$grad.tnorm(p)
        if(sum(abs(cscore)) < tol)
          return(list(estimate = p, value = lhood$ll.tnorm2(p), iter = 0))
        cur <- p
        for(i in 1:maxit) {
          inverseHess <- solve(lhood$hessian.tnorm(cur))
          cscore <- lhood$grad.tnorm(cur)
          new <- cur - cscore %*% inverseHess
          if (new[2] <= 0) stop("Sigma < 0")
          cscore <- lhood$grad.tnorm(new)
          
          if(((abs(lhood$ll.tnorm2(cur)- lhood$ll.tnorm2(new))/(lhood$ll.tnorm2(cur))) < tol))
            return(list(estimate = new, value= lhood$ll.tnorm2(new), iter = i))
          cur <- new
        }
        
        return(list(estimate = new, value= lhood$ll.tnorm2(new), iter = i))
      }
      
      ##################################################################################
      ###### Based on the MLE Functions (mklhood) and NewtonRaphson Function
      ###### (NewtonRaphsonLike), This function estimates the MEAN and SD from the
      ###### Truncated using Newton Raphson.
      ##################################################################################
      
      ## missingdata = matrix where rows = features, columns = samples
      ## perc = if %MVs > perc then just sample mean / SD
      ## iter = # iterations in NR algorithm
      
      EstimatesComputation <- function(missingdata, perc, iter=50) {
        
        ## 2 column matrix where column 1 = means, column 2 = SD
        ParamEstim <- matrix(NA, nrow = nrow(missingdata), ncol = 2)
        nsamp <- ncol(missingdata)
        
        ## sample means / SDs
        ParamEstim[,1] <- rowMeans(missingdata, na.rm = TRUE)
        ParamEstim[,2] <- apply(missingdata, 1, function(x) sd(x, na.rm = TRUE))
        
        ## Case 1: missing % > perc => use sample mean / SD
        na.sum <- apply(missingdata, 1, function(x) sum(is.na(x)))
        idx1 <- which(na.sum/nsamp >= perc)
        
        ## Case 2: sample mean > 3 SD away from LOD => use sample mean / SD
        lod <- min(missingdata, na.rm=TRUE) ## why use the min of whole data set??????
        idx2 <- which(ParamEstim[,1] > 3*ParamEstim[,2] + lod)
        
        ## Case 3: for all others, use NR method to obtain truncated mean / SD estimate
        idx.nr <- setdiff(1:nrow(missingdata), c(idx1, idx2))
        ## t = limits of integration (LOD and upper)
        upplim <- max(missingdata, na.rm=TRUE) + 2*max(ParamEstim[,2])
        for (i in idx.nr) {
          Likelihood <- mklhood(missingdata[i,], t=c(lod, upplim))
          res <- tryCatch(NewtonRaphsonLike(Likelihood, p = ParamEstim[i,]),
                          error = function(e) 1000)
          
          if (length(res) == 1) {
            next
          } else if (res$iter >= iter) {
            next
          } else {
            ParamEstim[i,] <- as.numeric(res$estimate)
          }
        }
        return(ParamEstim)
      }
      
      
      
      ####################################################################
      #### This Function imputes the data BASED on KNN-EUCLIDEAN
      ####################################################################
      
      ## data = data set to be imputed, where rows = features, columns = samples
      ## k    = number of neighbors for imputing values
      ## rm.na, rm.nan, rm.inf = whether NA, NaN, and Inf values should be imputed
      
      
      KNNEuc <- function (data, k, rm.na = TRUE, rm.nan = TRUE, rm.inf = TRUE) {
        
        nr <- dim(data)[1]
        
        imp.knn <- data
        imp.knn[is.finite(data) == FALSE] <- NA
        t.data<-t(data)
        
        mv.ind <- which(is.na(imp.knn), arr.ind = TRUE)
        arrays <- unique(mv.ind[, 2])
        array.ind <- match(arrays, mv.ind[, 2])
        nfeatures <- 1:nr
        
        for (i in 1:length(arrays)) {
          set <- array.ind[i]:min((array.ind[(i + 1)] - 1), dim(mv.ind)[1], na.rm = TRUE)
          cand.features <- nfeatures[-unique(mv.ind[set, 1])]
          cand.vectors <- t.data[,cand.features]
          exp.num <- arrays[i]
          
          for (j in set) {
            feature.num <- mv.ind[j, 1]
            tar.vector <- data[feature.num,]
            
            dist <- sqrt(colMeans((tar.vector-cand.vectors)^2, na.rm = TRUE))
            dist[is.nan(dist) | is.na(dist)] <- Inf
            dist[dist==0] <- ifelse(is.finite(min(dist[dist>0])), min(dist[dist>0])/2, 1)
            
            if (sum(is.finite(dist)) < k) {
              stop(message = "Fewer than K finite distances found")
            }
            k.features.ind <- order(dist)[1:k]
            k.features <- cand.features[k.features.ind]
            wghts <- 1/dist[k.features.ind]/sum(1/dist[k.features.ind])
            imp.knn[feature.num, exp.num] <- wghts %*% data[k.features, exp.num]
          }
        }
        
        if (!rm.na) {
          imp.knn[is.na(data) == TRUE & is.nan(data) == FALSE] <- NA
        }
        if (!rm.inf) {
          index <- is.finite(data) == FALSE & is.na(data) == FALSE &
            is.nan(data) == FALSE
          imp.knn[index] <- data[index]
        }
        if (!rm.nan) {
          imp.knn[is.nan(data) == TRUE] <- NaN
        }
        return(imp.knn)
      }
      
      ####################################################################
      #### This Function imputes the data based on KNN-CORRELATION or
      #### KNN-TRUNCATION. The Parameter Estimates based on the Truncated
      #### Normal from EstimateComputation function is run on this function
      ####################################################################
      
      imputeKNN <- function (data, k , distance = "correlation",
                             rm.na = TRUE, rm.nan = TRUE, rm.inf = TRUE, perc=1,...) {
        
        if (!(is.matrix(data))) {
          stop(message = paste(deparse(substitute(data)),
                               " is not a matrix.", sep = ""))
        }
        
        distance <- match.arg(distance, c("correlation","truncation"))
        
        nr <- dim(data)[1]
        if (k < 1 | k > nr) {
          stop(message = "k should be between 1 and the number of rows")
        }
        
        if (distance=="correlation"){
          genemeans<-rowMeans(data,na.rm=TRUE)
          genesd<-apply(data, 1, function(x) sd(x, na.rm = TRUE))
          data<-(data-genemeans)/genesd
        }
        
        if (distance=="truncation"){
          
          ParamMat <- EstimatesComputation(data, perc = perc)
          
          genemeans<-ParamMat[,1]
          genesd<-ParamMat[,2]
          data<-(data-genemeans)/genesd
        }
        
        imp.knn <- data
        imp.knn[is.finite(data) == FALSE] <- NA
        t.data<-t(data)
        
        mv.ind <- which(is.na(imp.knn), arr.ind = TRUE)
        arrays <- unique(mv.ind[, 2])
        array.ind <- match(arrays, mv.ind[, 2])
        ngenes <- 1:nr
        
        for (i in 1:length(arrays)) {
          set <- array.ind[i]:min((array.ind[(i + 1)] - 1), dim(mv.ind)[1],
                                  na.rm = TRUE)
          cand.genes <- ngenes[-unique(mv.ind[set, 1])]
          cand.vectors <- t.data[,cand.genes]
          exp.num<- arrays[i]
          for (j in set) {
            
            gene.num <- mv.ind[j, 1]
            tar.vector <- data[gene.num,]
            
            r <- (cor(cand.vectors,tar.vector, use = "pairwise.complete.obs"))
            dist <- switch(distance,
                           correlation = (1 - abs(r)),
                           truncation = (1 - abs(r)))
            dist[is.nan(dist) | is.na(dist)] <- Inf
            dist[dist==0]<-ifelse(is.finite(min(dist[dist>0])), min(dist[dist>0])/2, 1)
            dist[abs(r) == 1] <- Inf
            
            if (sum(is.finite(dist)) < k) {
              stop(message = "Fewer than K finite distances found")
            }
            k.genes.ind <- order(dist)[1:k]
            k.genes <- cand.genes[k.genes.ind]
            
            wghts <- (1/dist[k.genes.ind]/sum(1/dist[k.genes.ind])) * sign(r[k.genes.ind])
            imp.knn[gene.num, exp.num] <- wghts %*% data[k.genes, exp.num]
          }
        }
        
        if (distance=="correlation") {
          imp.knn <- (imp.knn * genesd) + genemeans
        }
        
        if(distance=="truncation") {
          imp.knn <- (imp.knn * genesd) + genemeans
        }
        
        if (!rm.na) {
          imp.knn[is.na(data) == TRUE & is.nan(data) == FALSE] <- NA
        }
        if (!rm.inf) {
          index <- is.finite(data) == FALSE & is.na(data) == FALSE &
            is.nan(data) == FALSE
          imp.knn[index] <- data[index]
        }
        if (!rm.nan) {
          imp.knn[is.nan(data) == TRUE] <- NaN
        }
        return(imp.knn)
      }
    }
    sim_trKNN_wrapper <- function(data) {
      result <- data %>% as.matrix %>% t %>% imputeKNN(., k=10, distance='truncation', perc=0) %>% t
      return(result)
    }
    df1x <- sim_trKNN_wrapper(t(df1))
    df<-as.data.frame(t(df1x))
  }
  else if(method=="rf"){
    library(missForest)
    data_zero1 <- missForest(t(df1), maxiter =10,ntree = input$rfntrees,mtry=floor(nrow(df1)^(1/3)),verbose = TRUE)
    df<-t(data_zero1$ximp)
  }
  else if(method=="pi"){
    width <- input$piwidth
    downshift <- input$pidownshift
    for(i in 1:ncol(df1)){
      temp <- df1[[i]]
      if(sum(is.na(temp))>0){
        temp.sd <- width * sd(temp[!is.na(temp)], na.rm = TRUE)
        temp.mean <- mean(temp[!is.na(temp)], na.rm = TRUE) - downshift * sd(temp[!is.na(temp)], na.rm = TRUE)
        n.missing <- sum(is.na(temp))
        temp[is.na(temp)] <- rnorm(n.missing, mean = temp.mean, sd = temp.sd)
        df[[i]]<-temp
      }
    }
    df
  }
  else if(method=="grr"){
    library(DreamAI)
    df<-impute.RegImpute(data=as.matrix(df1), fillmethod = "row_mean", maxiter_RegImpute = 10,conv_nrmse = 1e-03)
  }
  else if(method=="gms"){
    library(GMSimpute)
    df<-GMS.Lasso(df1,nfolds=3,log.scale=FALSE,TS.Lasso=TRUE)
  }
  else{
    stop("Unspported methods so far~~")
  }
  df<-as.data.frame(df)
  df
}



#' Impute missing values in a data frame using various methods
#'
#' @param x A data frame or matrix containing missing values.
#' @param method Character string specifying the imputation method.
#'   One of: "zero", "minimum", "colmedian", "rowmedian", "knn",
#'   "seqknn", "bpca", "svd", "lls", "mle", "qrilc", "mindet", "minprob",
#'   "irm", "impseq", "impseqrob", "mice_norm", "mice_cart", "trknn",
#'   "rf", "pi", "grr", "gms".
#' @return A data frame of the same dimension as `x`, with NAs imputed.
#' @references
#' Shisheng Wang et al., *NAguideR: performing and prioritizing missing value imputations*,
#' Nucleic Acids Research, 2021; doi:10.1093/nar/gkaa498
#' @examples
#' df_filled <- nafunctions(my_data, method = "knn")
#' @export
naimpute <- function(x,
                     method = c("zero", "minimum", "colmedian", "rowmedian",
                                "knn", "seqknn", "bpca", "svd", "lls",
                                "mle", "qrilc", "mindet", "minprob",
                                "irm", "impseq", "impseqrob",
                                "mice_norm", "mice_cart", "trknn",
                                "rf", "pi", "grr", "gms")) {
  method <- match.arg(tolower(method), method)
  df_orig <- as.data.frame(x)
  df_work <- df_orig
  
  imputed <- switch(
    method,
    
    # ---- Simple Replacements ----
    zero = {
      df_work[is.na(df_work)] <- 0
      df_work
    },
    minimum = {
      df_work[is.na(df_work)] <- min(df_orig, na.rm = TRUE)
      df_work
    },
    
    # ---- Median Imputation ----
    colmedian = {
      if (!requireNamespace("e1071", quietly = TRUE)) {
        stop("Package 'e1071' required for colmedian", call. = FALSE)
      }
      e1071::impute(df_orig, what = "median")
    },
    rowmedian = {
      if (!requireNamespace("e1071", quietly = TRUE)) {
        stop("Package 'e1071' required for rowmedian", call. = FALSE)
      }
      t(e1071::impute(t(df_orig), what = "median"))
    },
    
    # ---- KNN-based Methods ----
    knn = {
      if (!requireNamespace("impute", quietly = TRUE)) {
        stop("Package 'impute' required for knn", call. = FALSE)
      }
      impute::impute.knn(as.matrix(df_orig), k = 10)$data
    },
    seqknn = {
      if (!requireNamespace("SeqKnn", quietly = TRUE)) {
        stop("Package 'SeqKnn' required for seqknn", call. = FALSE)
      }
      SeqKnn::SeqKNN(df_orig, k = 10)
    },
    
    # ---- PCA / SVD Methods ----
    bpca = {
      if (!requireNamespace("pcaMethods", quietly = TRUE)) {
        stop("Package 'pcaMethods' required for bpca", call. = FALSE)
      }
      pca <- pcaMethods::pca(as.matrix(df_orig),
                             nPcs = ncol(df_orig) - 1,
                             method = "bpca", maxSteps = 100)
      pcaMethods::completeObs(pca)
    },
    svd = {
      if (!requireNamespace("pcaMethods", quietly = TRUE)) {
        stop("Package 'pcaMethods' required for svd", call. = FALSE)
      }
      pca <- pcaMethods::pca(as.matrix(df_orig),
                             nPcs = ncol(df_orig) - 1,
                             method = "svdImpute")
      pcaMethods::completeObs(pca)
    },
    lls = {
      if (!requireNamespace("pcaMethods", quietly = TRUE)) {
        stop("Package 'pcaMethods' required for lls", call. = FALSE)
      }
      # Local least squares on transposed data
      pca <- pcaMethods::llsImpute(t(df_orig), k = 10)
      t(pcaMethods::completeObs(pca))
    },
    
    # ---- Model-based Methods ----
    mle = {
      if (!requireNamespace("norm", quietly = TRUE)) {
        stop("Package 'norm' required for mle", call. = FALSE)
      }
      mat <- as.matrix(df_orig)
      pre <- norm::prelim.norm(mat)
      em  <- norm::em.norm(pre)
      norm::rngseed(123)
      norm::imp.norm(pre, em, mat)
    },
    qrilc = {
      if (!requireNamespace("imputeLCMD", quietly = TRUE)) {
        stop("Package 'imputeLCMD' required for qrilc", call. = FALSE)
      }
      tmp   <- t(df_orig)
      filled <- imputeLCMD::impute.QRILC(tmp, tune.sigma = 1)[[1]]
      t(filled)
    },
    mindet = {
      if (!requireNamespace("imputeLCMD", quietly = TRUE)) {
        stop("Package 'imputeLCMD' required for mindet", call. = FALSE)
      }
      imputeLCMD::impute.MinDet(as.matrix(df_orig), q = 0.01)
    },
    minprob = {
      if (!requireNamespace("imputeLCMD", quietly = TRUE)) {
        stop("Package 'imputeLCMD' required for minprob", call. = FALSE)
      }
      imputeLCMD::impute.MinProb(as.matrix(df_orig), q = 0.01, tune.sigma = 1)
    },
    
    # ---- Iterative / Sequential Methods ----
    irm = {
      if (!requireNamespace("VIM", quietly = TRUE)) {
        stop("Package 'VIM' required for irm", call. = FALSE)
      }
      out <- VIM::irmi(df_orig, trace = TRUE, imp_var = FALSE)
      `rownames<-`(out, rownames(df_orig))
    },
    impseq = {
      if (!requireNamespace("rrcovNA", quietly = TRUE)) {
        stop("Package 'rrcovNA' required for impseq", call. = FALSE)
      }
      rrcovNA::impSeq(df_orig)
    },
    impseqrob = {
      if (!requireNamespace("rrcovNA", quietly = TRUE)) {
        stop("Package 'rrcovNA' required for impseqrob", call. = FALSE)
      }
      out <- rrcovNA::impSeqRob(df_orig, alpha = 0.9)
      out$x
    },
    
    # ---- MICE Methods ----
    mice_norm = ,
    mice_cart = {
      if (!requireNamespace("mice", quietly = TRUE)) {
        stop("Package 'mice' required for mice methods", call. = FALSE)
      }
      m <- 5
      meth <- if (method == "mice_norm") "norm" else "cart"
      imp <- mice::mice(df_orig, m = m, seed = 1234, method = meth)
      avg <- Reduce(`+`, lapply(1:m, function(i) mice::complete(imp, action = i)))
      avg / m
    },
    
    # ---- Specialty Methods ----
    trknn = {
      # Assumes that 'Trunc_KNN/Imput_funcs.r' defines imputeKNN()
      # source('Trunc_KNN/Imput_funcs.r') # only call the minimal requirements
      {
        ##################################################################################
        #### MLE for the Truncated Normal
        #### Creating a Function that Returns the Log Likelihood, Gradient and
        #### Hessian Functions
        ##################################################################################
        
        ## data = numeric vector
        ## t    = truncation limits
        mklhood <- function(data, t, ...) {
          
          data <- na.omit(data)
          n <- length(data)
          t <- sort(t)
          
          psi<-function(y, mu, sigma){
            exp(-(y-mu)^2/(2*sigma^2))/(sigma*sqrt(2*pi))
          }
          
          psi.mu<-function(y,mu,sigma){
            exp(-(y-mu)^2/(2*sigma^2)) * ((y-mu)/(sigma^3*sqrt(2*pi)))
          }
          
          psi.sigma<-function(y,mu,sigma){
            exp(-(y-mu)^2/(2*sigma^2)) *
              (((y-mu)^2)/(sigma^4*sqrt(2*pi)) - 1/(sigma^2*sqrt(2*pi)))
          }
          
          psi2.mu<-function(y,mu,sigma){
            exp(-(y - mu)^2/(2*sigma^2)) *
              (((y - mu)^2)/(sigma^5*sqrt(2*pi))-1/(sigma^3*sqrt(2*pi)))
          }
          
          psi2.sigma<-function(y,mu,sigma){
            exp(-(y-mu)^2/(2*sigma^2)) *
              ((2)/(sigma^3*sqrt(2*pi)) - (5*(y-mu))/(sigma^5*sqrt(2*pi)) +
                 ((y-mu)^4)/(sigma^7*sqrt(2*pi)))
          }
          
          psi12.musig<-function(y,mu,sigma){
            exp(-(y-mu)^2/(2*sigma^2)) *
              (((y-mu)^3)/(sigma^6*sqrt(2*pi)) - (3*(y-mu))/(sigma^4*sqrt(2*pi)))
          }
          
          ll.tnorm2<-function(p){
            out <- (-n*log(pnorm(t[2],p[1],p[2])-pnorm(t[1],p[1],p[2]))) -
              (n*log(sqrt(2*pi*p[2]^2))) - (sum((data-p[1])^2)/(2*p[2]^2))
            -1*out
          }
          
          grad.tnorm<-function(p){
            g1 <- (-n*(integrate(psi.mu,t[1],t[2],mu=p[1],sigma=p[2], stop.on.error = FALSE)$value) /
                     (pnorm(max(t),p[1],p[2])-pnorm(min(t),p[1],p[2]))) - ((n*p[1]-sum(data))/p[2]^2)
            g2 <- (-n*(integrate(psi.sigma,t[1],t[2],mu=p[1],sigma=p[2], stop.on.error = FALSE)$value) /
                     (pnorm(max(t),p[1],p[2])-pnorm(min(t),p[1],p[2]))) - ((n)/(p[2])) + ((sum((data-p[1])^2))/(p[2]^3))
            out <- c(g1,g2)
            return(out)
          }
          
          hessian.tnorm<-function(p){
            
            h1<- -n*(integrate(psi,t[1],t[2],mu=p[1],sigma=p[2], stop.on.error = FALSE)$value *
                       integrate(psi2.mu,t[1],t[2],mu=p[1],sigma=p[2], stop.on.error = FALSE)$value -
                       integrate(psi.mu,t[1],t[2],mu=p[1],sigma=p[2], stop.on.error = FALSE)$value^2) /
              (integrate(psi,t[1],t[2],mu=p[1],sigma=p[2], stop.on.error = FALSE)$value^2) -
              n/(p[2]^2)
            
            h3<- -n*(integrate(psi,t[1],t[2],mu=p[1],sigma=p[2], stop.on.error = FALSE)$value *
                       integrate(psi12.musig,t[1],t[2],mu=p[1],sigma=p[2], stop.on.error = FALSE)$value -
                       integrate(psi.mu,t[1],t[2],mu=p[1],sigma=p[2], stop.on.error = FALSE)$value *
                       integrate(psi.sigma,t[1],t[2],mu=p[1],sigma=p[2], stop.on.error = FALSE)$value) /
              (integrate(psi,t[1],t[2],mu=p[1],sigma=p[2], stop.on.error = FALSE)$value^2) +
              (2*(n*p[1]-sum(data)))/(p[2]^3)
            
            h2<- -n*(integrate(psi,t[1],t[2],mu=p[1],sigma=p[2], stop.on.error = FALSE)$value *
                       integrate(psi2.sigma,t[1],t[2],mu=p[1],sigma=p[2], stop.on.error = FALSE)$value -
                       integrate(psi.sigma,t[1],t[2],mu=p[1],sigma=p[2], stop.on.error = FALSE)$value^2) /
              (integrate(psi,t[1],t[2],mu=p[1],sigma=p[2], stop.on.error = FALSE)$value^2) +
              (n)/(p[2]^2)-(3*sum((data-p[1])^2))/(p[2]^4)
            
            H<-matrix(0,nrow=2,ncol=2)
            H[1,1]<-h1
            H[2,2]<-h2
            H[1,2]<-H[2,1]<-h3
            return(H)
          }
          
          
          return(list(ll.tnorm2 = ll.tnorm2, grad.tnorm = grad.tnorm, hessian.tnorm = hessian.tnorm))
        }
        ##################################################################################
        ###### Newton Raphson Function
        ###### This takes in the Objects Returned from mklhood Function above
        ##################################################################################
        
        NewtonRaphsonLike <- function(lhood, p, tol = 1e-07, maxit = 100) {
          
          cscore <- lhood$grad.tnorm(p)
          if(sum(abs(cscore)) < tol)
            return(list(estimate = p, value = lhood$ll.tnorm2(p), iter = 0))
          cur <- p
          for(i in 1:maxit) {
            inverseHess <- solve(lhood$hessian.tnorm(cur))
            cscore <- lhood$grad.tnorm(cur)
            new <- cur - cscore %*% inverseHess
            if (new[2] <= 0) stop("Sigma < 0")
            cscore <- lhood$grad.tnorm(new)
            
            if(((abs(lhood$ll.tnorm2(cur)- lhood$ll.tnorm2(new))/(lhood$ll.tnorm2(cur))) < tol))
              return(list(estimate = new, value= lhood$ll.tnorm2(new), iter = i))
            cur <- new
          }
          
          return(list(estimate = new, value= lhood$ll.tnorm2(new), iter = i))
        }
        
        ##################################################################################
        ###### Based on the MLE Functions (mklhood) and NewtonRaphson Function
        ###### (NewtonRaphsonLike), This function estimates the MEAN and SD from the
        ###### Truncated using Newton Raphson.
        ##################################################################################
        
        ## missingdata = matrix where rows = features, columns = samples
        ## perc = if %MVs > perc then just sample mean / SD
        ## iter = # iterations in NR algorithm
        
        EstimatesComputation <- function(missingdata, perc, iter=50) {
          
          ## 2 column matrix where column 1 = means, column 2 = SD
          ParamEstim <- matrix(NA, nrow = nrow(missingdata), ncol = 2)
          nsamp <- ncol(missingdata)
          
          ## sample means / SDs
          ParamEstim[,1] <- rowMeans(missingdata, na.rm = TRUE)
          ParamEstim[,2] <- apply(missingdata, 1, function(x) sd(x, na.rm = TRUE))
          
          ## Case 1: missing % > perc => use sample mean / SD
          na.sum <- apply(missingdata, 1, function(x) sum(is.na(x)))
          idx1 <- which(na.sum/nsamp >= perc)
          
          ## Case 2: sample mean > 3 SD away from LOD => use sample mean / SD
          lod <- min(missingdata, na.rm=TRUE) ## why use the min of whole data set??????
          idx2 <- which(ParamEstim[,1] > 3*ParamEstim[,2] + lod)
          
          ## Case 3: for all others, use NR method to obtain truncated mean / SD estimate
          idx.nr <- setdiff(1:nrow(missingdata), c(idx1, idx2))
          ## t = limits of integration (LOD and upper)
          upplim <- max(missingdata, na.rm=TRUE) + 2*max(ParamEstim[,2])
          for (i in idx.nr) {
            Likelihood <- mklhood(missingdata[i,], t=c(lod, upplim))
            res <- tryCatch(NewtonRaphsonLike(Likelihood, p = ParamEstim[i,]),
                            error = function(e) 1000)
            
            if (length(res) == 1) {
              next
            } else if (res$iter >= iter) {
              next
            } else {
              ParamEstim[i,] <- as.numeric(res$estimate)
            }
          }
          return(ParamEstim)
        }
        
        
        
        ####################################################################
        #### This Function imputes the data BASED on KNN-EUCLIDEAN
        ####################################################################
        
        ## data = data set to be imputed, where rows = features, columns = samples
        ## k    = number of neighbors for imputing values
        ## rm.na, rm.nan, rm.inf = whether NA, NaN, and Inf values should be imputed
        
        
        KNNEuc <- function (data, k, rm.na = TRUE, rm.nan = TRUE, rm.inf = TRUE) {
          
          nr <- dim(data)[1]
          
          imp.knn <- data
          imp.knn[is.finite(data) == FALSE] <- NA
          t.data<-t(data)
          
          mv.ind <- which(is.na(imp.knn), arr.ind = TRUE)
          arrays <- unique(mv.ind[, 2])
          array.ind <- match(arrays, mv.ind[, 2])
          nfeatures <- 1:nr
          
          for (i in 1:length(arrays)) {
            set <- array.ind[i]:min((array.ind[(i + 1)] - 1), dim(mv.ind)[1], na.rm = TRUE)
            cand.features <- nfeatures[-unique(mv.ind[set, 1])]
            cand.vectors <- t.data[,cand.features]
            exp.num <- arrays[i]
            
            for (j in set) {
              feature.num <- mv.ind[j, 1]
              tar.vector <- data[feature.num,]
              
              dist <- sqrt(colMeans((tar.vector-cand.vectors)^2, na.rm = TRUE))
              dist[is.nan(dist) | is.na(dist)] <- Inf
              dist[dist==0] <- ifelse(is.finite(min(dist[dist>0])), min(dist[dist>0])/2, 1)
              
              if (sum(is.finite(dist)) < k) {
                stop(message = "Fewer than K finite distances found")
              }
              k.features.ind <- order(dist)[1:k]
              k.features <- cand.features[k.features.ind]
              wghts <- 1/dist[k.features.ind]/sum(1/dist[k.features.ind])
              imp.knn[feature.num, exp.num] <- wghts %*% data[k.features, exp.num]
            }
          }
          
          if (!rm.na) {
            imp.knn[is.na(data) == TRUE & is.nan(data) == FALSE] <- NA
          }
          if (!rm.inf) {
            index <- is.finite(data) == FALSE & is.na(data) == FALSE &
              is.nan(data) == FALSE
            imp.knn[index] <- data[index]
          }
          if (!rm.nan) {
            imp.knn[is.nan(data) == TRUE] <- NaN
          }
          return(imp.knn)
        }
        
        ####################################################################
        #### This Function imputes the data based on KNN-CORRELATION or
        #### KNN-TRUNCATION. The Parameter Estimates based on the Truncated
        #### Normal from EstimateComputation function is run on this function
        ####################################################################
        
        imputeKNN <- function (data, k , distance = "correlation",
                               rm.na = TRUE, rm.nan = TRUE, rm.inf = TRUE, perc=1,...) {
          
          if (!(is.matrix(data))) {
            stop(message = paste(deparse(substitute(data)),
                                 " is not a matrix.", sep = ""))
          }
          
          distance <- match.arg(distance, c("correlation","truncation"))
          
          nr <- dim(data)[1]
          if (k < 1 | k > nr) {
            stop(message = "k should be between 1 and the number of rows")
          }
          
          if (distance=="correlation"){
            genemeans<-rowMeans(data,na.rm=TRUE)
            genesd<-apply(data, 1, function(x) sd(x, na.rm = TRUE))
            data<-(data-genemeans)/genesd
          }
          
          if (distance=="truncation"){
            
            ParamMat <- EstimatesComputation(data, perc = perc)
            
            genemeans<-ParamMat[,1]
            genesd<-ParamMat[,2]
            data<-(data-genemeans)/genesd
          }
          
          imp.knn <- data
          imp.knn[is.finite(data) == FALSE] <- NA
          t.data<-t(data)
          
          mv.ind <- which(is.na(imp.knn), arr.ind = TRUE)
          arrays <- unique(mv.ind[, 2])
          array.ind <- match(arrays, mv.ind[, 2])
          ngenes <- 1:nr
          
          for (i in 1:length(arrays)) {
            set <- array.ind[i]:min((array.ind[(i + 1)] - 1), dim(mv.ind)[1],
                                    na.rm = TRUE)
            cand.genes <- ngenes[-unique(mv.ind[set, 1])]
            cand.vectors <- t.data[,cand.genes]
            exp.num<- arrays[i]
            for (j in set) {
              
              gene.num <- mv.ind[j, 1]
              tar.vector <- data[gene.num,]
              
              r <- (cor(cand.vectors,tar.vector, use = "pairwise.complete.obs"))
              dist <- switch(distance,
                             correlation = (1 - abs(r)),
                             truncation = (1 - abs(r)))
              dist[is.nan(dist) | is.na(dist)] <- Inf
              dist[dist==0]<-ifelse(is.finite(min(dist[dist>0])), min(dist[dist>0])/2, 1)
              dist[abs(r) == 1] <- Inf
              
              if (sum(is.finite(dist)) < k) {
                stop(message = "Fewer than K finite distances found")
              }
              k.genes.ind <- order(dist)[1:k]
              k.genes <- cand.genes[k.genes.ind]
              
              wghts <- (1/dist[k.genes.ind]/sum(1/dist[k.genes.ind])) * sign(r[k.genes.ind])
              imp.knn[gene.num, exp.num] <- wghts %*% data[k.genes, exp.num]
            }
          }
          
          if (distance=="correlation") {
            imp.knn <- (imp.knn * genesd) + genemeans
          }
          
          if(distance=="truncation") {
            imp.knn <- (imp.knn * genesd) + genemeans
          }
          
          if (!rm.na) {
            imp.knn[is.na(data) == TRUE & is.nan(data) == FALSE] <- NA
          }
          if (!rm.inf) {
            index <- is.finite(data) == FALSE & is.na(data) == FALSE &
              is.nan(data) == FALSE
            imp.knn[index] <- data[index]
          }
          if (!rm.nan) {
            imp.knn[is.nan(data) == TRUE] <- NaN
          }
          return(imp.knn)
        }
      }
      result <- df_orig %>%
        as.matrix() %>%
        t() %>%
        imputeKNN(k = 10, distance = "truncation", perc = 0) %>%
        t()
      as.data.frame(result)
    },
    rf = {
      if (!requireNamespace("missForest", quietly = TRUE)) {
        stop("Package 'missForest' required for rf", call. = FALSE)
      }
      out <- missForest::missForest(t(df_orig), maxiter = 10,
                                    mtry = floor(nrow(df_orig)^(1/3)),
                                    verbose = TRUE)
      t(out$ximp)
    },
    pi = {
      # Paramters `width` and `downshift` should be passed in via `...` or environment
      width     <- getOption("nafuncs.pi.width", default = 1)
      downshift <- getOption("nafuncs.pi.downshift", default = 1)
      df_work[] <- lapply(df_orig, function(col) {
        if (anyNA(col)) {
          non_na   <- col[!is.na(col)]
          sd_val   <- width * sd(non_na, na.rm = TRUE)
          mean_val <- mean(non_na, na.rm = TRUE) - downshift * sd(non_na, na.rm = TRUE)
          col[is.na(col)] <- rnorm(sum(is.na(col)), mean = mean_val, sd = sd_val)
        }
        col
      })
      df_work
    },
    grr = {
      if (!requireNamespace("DreamAI", quietly = TRUE)) {
        stop("Package 'DreamAI' required for grr", call. = FALSE)
      }
      if (!requireNamespace("glmnet", quietly = TRUE)) {
        stop("Package 'glmnet' required for grr", call. = FALSE)
      }
      if (!requireNamespace("import", quietly=TRUE)) install.packages("import")
      import::into(.into = "DreamAI", .from = "glmnet", cv.glmnet) # I really don’t want to attach glmnet, just inject the function into DreamAI’s namespace
      DreamAI::impute.RegImpute(as.matrix(df_orig),
                                fillmethod = "row_mean",
                                maxiter_RegImpute = 10,
                                conv_nrmse = 1e-03)
    },
    gms = {
      if (!requireNamespace("GMSimpute", quietly = TRUE)) {
        stop("Package 'GMSimpute' required for gms", call. = FALSE)
      }
      GMSimpute::GMS.Lasso(df_orig, nfolds = 3, log.scale = FALSE, TS.Lasso = TRUE)
    },
    
    stop(sprintf("Unsupported method '%s'", method), call. = FALSE)
  )
  
  as.data.frame(imputed)
}

