# Reference: DIALib-QC_RPlot.pl

# --------- reset environment --------------------
rm(list = ls())
libs <- c('RColorBrewer', 'ggplot2', 'scales', 'ggpubr')
table(sapply(libs, require, character.only=TRUE))
rm(libs)
setwd('//172.16.13.114/share/members/jiangwenhao/TPHP/202200720/LibQC')

# ------- function: plotting ------------------------
plotting <- function(df, fasPath = '//172.16.13.114/share/members/jiangwenhao/fasta/swissprot-homo+sapiens-20200507_20377_iRT_11.fasta', pdfName, rlt_dir = './', allColor = c('#A2565B', '#6E8E84', '#1A476F', '#E37E00', '#90353A'),
                     allFill = c('#C79A9D', '#B7C7C2', '#8DA3B7', '#F1BE80', '#C89A9D')){
  # DIA-LibQC (For FragPipe-Easypqp library format)
  df_mz <- df %>% dplyr::distinct(ModifiedPeptideSequence, PrecursorCharge, PrecursorMz)
  # Plot A: precursor m/z distribution
  p1 <- ggplot(df_mz, aes(x = PrecursorMz))+
    geom_histogram(aes(y = (..count..)/sum(..count..)),
                   color = allColor[1],
                   fill = allFill[1])+
    #geom_text(aes(y = ((..count..)/sum(..count..)), label = scales::percent((..count..)/sum(..count..))), stat = 'count', vjust = -10)+
    #xlim(350, 1300)+
    scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0.01, 0.1), add = c(0, 0)))+
    labs(title = 'Precursor m/z', x = "Precursor m/z", y = "Frequency in percent (%)")+
    theme_bw()+
    theme(legend.text = element_text(size = 18, color = "black"),legend.position = 'top',
          legend.title = element_text(size = 18, color="black") ,
          panel.grid.major =element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+
    theme(panel.grid =element_blank())+
    theme(axis.text = element_text(size = 15,color = "black"))+
    theme(axis.text.x = element_text(size = 15,color = "black", angle = 0))+
    theme(plot.title = element_text(size = 20, face = "bold"))+
    theme(plot.subtitle=element_text(size = 20, hjust = 0, color="black"))+
    theme(axis.title.x=element_text(size = 18, hjust = 0.5, color="black"))+
    theme(axis.title.y=element_text(size = 18, hjust = 0.5, color="black"))
  # Plot B: precursor charge distribution
  p2 <- ggplot(df_mz, aes(x = PrecursorCharge))+
    geom_bar(aes(y = (..count..)/sum(..count..)),
             width = 0.6,
             color = allColor[2],
             fill = allFill[2],
    )+
    geom_text(aes(y = ((..count..)/sum(..count..)), label = scales::percent((..count..)/sum(..count..))), stat = 'count', size = 6, hjust = 0.5, vjust = -0.2, position = "stack", color = 'black')+
    scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0.01, 0.1), add = c(0, 0)))+
    labs(title = 'Precursor charge', x = "Precursor charge", y = "Precursor in percent (%)")+
    theme_bw()+
    theme(legend.text = element_text(size = 18, color = "black"),legend.position = 'top',
          legend.title = element_text(size = 18, color="black") ,
          panel.grid.major =element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+
    theme(panel.grid =element_blank())+
    theme(axis.text = element_text(size = 15,color = "black"))+
    theme(axis.text.x = element_text(size = 15,color = "black", angle = 0))+
    theme(plot.title = element_text(size = 20, face = "bold"))+
    theme(plot.subtitle=element_text(size = 20, hjust = 0, color="black"))+
    theme(axis.title.x=element_text(size = 18, hjust = 0.5, color="black"))+
    theme(axis.title.y=element_text(size = 18, hjust = 0.5, color="black"))
  
  
  df_rt <- df %>% dplyr::distinct(ModifiedPeptideSequence, PrecursorCharge, NormalizedRetentionTime)
  df_rt <- df_rt[df_rt$PrecursorCharge %in% c(2, 3), ]
  tb_tmp <- table(df_rt$ModifiedPeptideSequence)
  pairs <- names(tb_tmp[tb_tmp == 2])
  df_rt <- df_rt[df_rt$ModifiedPeptideSequence %in% pairs, ]
  df_rt <- reshape2::dcast(df_rt, ModifiedPeptideSequence~PrecursorCharge, value.var = 'NormalizedRetentionTime')
  colnames(df_rt) <- c('modSeq', 'RT2value', 'RT3value')
  #df_rt %<>% dplyr::filter(RT2value <= 250) %>% dplyr::filter(RT3value <= 250)
  # Plot C: +2/+3 RT linear regression
  correlation <- sprintf("%1f", cor(df_rt$RT2value, df_rt$RT3value, method = "pearson"))
  N <- nrow(df_rt)
  lbl = paste0('n=', N, ', R2=', correlation)
  # scatter plot
  p3 <- ggplot(df_rt, aes(x = RT2value, y = RT3value))+
    stat_density2d(aes(fill = ..density..), geom = "tile", contour = FALSE, n = 300) +
    scale_fill_continuous(low = "white", high = allFill[3])+
    scale_y_continuous(expand = expansion(mult = c(0.01, 0.1), add = c(0, 0)))+
    #geom_point(shape = 20, size = 1, color = "dodgerblue4", alpha = 0.1)+
    labs(title = 'Library +2/+3 pair iRT correlation', x = "+2 iRT", y = "+3 iRT")+
    annotate("text", x = min(df_rt$RT2value), y = max(df_rt$RT3value) * 1.05, label = lbl, hjust = 0, vjust = 0, size = 5)+
    theme_bw()+
    theme(#legend.text = element_text(size = 18, color = "black"),legend.position = 'top',
      #legend.title = element_text(size = 18, color="black") ,
      legend.position = 'none',
      panel.grid.major =element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"))+
    theme(panel.grid =element_blank())+
    theme(axis.text = element_text(size = 15,color = "black"))+
    theme(axis.text.x = element_text(size = 15,color = "black", angle = 0))+
    theme(plot.title = element_text(size = 20, face = "bold"))+
    theme(plot.subtitle=element_text(size = 20, hjust = 0, color="black"))+
    theme(axis.title.x=element_text(size = 18, hjust = 0.5, color="black"))+
    theme(axis.title.y=element_text(size = 18, hjust = 0.5, color="black"))
  
  df_len <- df %>% dplyr::distinct(PeptideSequence)
  df_len$PeptideLength <- nchar(df_len$PeptideSequence)
  # Plot D: Peptide length distribution
  p4 <- ggplot(df_len, aes(x = PeptideLength))+
    geom_bar(aes(y = (..count..)/sum(..count..)),
             width = 1,
             color = allColor[4],
             fill = allFill[4])+
    scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0.01, 0.1), add = c(0, 0)))+
    labs(title = 'Peptide length', x = 'Peptide length', y = 'Peptide in percent (%)')+
    theme_bw()+
    theme(legend.text = element_text(size = 18, color = "black"),legend.position = 'top',
          legend.title = element_text(size = 18, color="black") ,
          panel.grid.major =element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+
    theme(panel.grid =element_blank())+
    theme(axis.text = element_text(size = 15,color = "black"))+
    theme(axis.text.x = element_text(size = 15,color = "black", angle = 0))+
    theme(plot.title = element_text(size = 20, face = "bold"))+
    theme(plot.subtitle=element_text(size = 20, hjust = 0, color="black"))+
    theme(axis.title.x=element_text(size = 18, hjust = 0.5, color="black"))+
    theme(axis.title.y=element_text(size = 18, hjust = 0.5, color="black"))
  
  
  df_mod <- df %>% dplyr::distinct(ModifiedPeptideSequence)
  df_mod['[+42]'] <- unlist(lapply(df_mod$ModifiedPeptideSequence, function(e){ return(grepl('(UniMod:1)', e, fixed = T)) }))
  df_mod['[+57]'] <- unlist(lapply(df_mod$ModifiedPeptideSequence, function(e){ return(grepl('(UniMod:4)', e, fixed = T)) }))
  df_mod['[+16]'] <- unlist(lapply(df_mod$ModifiedPeptideSequence, function(e){ return(grepl('(UniMod:35)', e, fixed = T)) }))
  vec_mod <- unlist(lapply(df_mod[, -1], sum))
  df_mod <- data.frame(mod_type = names(vec_mod), freq = vec_mod)
  # Plot E: Modification counts
  p5 <- ggplot(df_mod, aes(mod_type, freq))+
    geom_col(color = allColor[5],
             fill = allFill[5],
             width = 0.5)+
    geom_text(aes(label = freq),
              size = 6, hjust = 0.5, vjust = -0.2, position = "stack", color = 'black')+
    scale_y_continuous(expand = expansion(mult = c(0.01, 0.1), add = c(0, 0)))+
    labs(title = 'Modifications', x = 'Modification type', y = 'Number of modifications')+
    theme_bw()+
    theme(legend.text = element_text(size = 18, color = "black"),legend.position = 'top',
          legend.title = element_text(size = 18, color="black") ,
          panel.grid.major =element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(color = "black"))+
    theme(panel.grid =element_blank())+
    theme(axis.text = element_text(size = 15,color = "black"))+
    theme(axis.text.x = element_text(size = 15,color = "black", angle = 0))+
    theme(plot.title = element_text(size = 20, face = "bold"))+
    theme(plot.subtitle=element_text(size = 20, hjust = 0, color="black"))+
    theme(axis.title.x=element_text(size = 18, hjust = 0.5, color="black"))+
    theme(axis.title.y=element_text(size = 18, hjust = 0.5, color="black"))
  
  
  df_peppro <- df %>% dplyr::distinct(ProteinId, PeptideSequence)
  df_peppro <- data.frame(table(table(df_peppro$ProteinId))); colnames(df_peppro) <- c('pep_num', 'count')
  if (sum(as.numeric(as.character(df_peppro$pep_num)) >= 8)){
    df_peppro$pep_num <- as.numeric(as.character(df_peppro$pep_num))
    over8 <- sum(df_peppro[df_peppro$pep_num >= 8, 'count'])
    df_peppro <- df_peppro[df_peppro$pep_num < 8, ]
    df_peppro <- rbind(df_peppro, c('>=8', over8))
    df_peppro$pep_num <- factor(df_peppro$pep_num, levels = c('1', '2', '3', '4', '5', '6', '7', '>=8'), ordered = T)
    df_peppro$count <- as.numeric(df_peppro$count)
    #df_peppro[1, 'sum'] <- df_peppro$count[1]
    #for (i in 2:nrow(df_peppro)){
    #  df_peppro[i, 'sum'] <- df_peppro[i-1, 'sum'] + df_peppro$count[i]
    #}
    #df_peppro <- df_peppro[which(df_peppro$sum < (0.95 * sum(df_peppro$count))), ]
  }
  # Plot F: Peptides per protein
  p6 <- ggplot(df_peppro, aes(pep_num, count))+
    geom_col(color = allColor[6], fill = allFill[6])+
    geom_text(aes(label = count), size = 6, hjust = 0.5, vjust = -0.2, position = "stack", color = 'black')+
    labs(title = 'Peptides per protein', x = 'Peptides per protein', y = 'Number of protiens')+
    scale_x_discrete("Peptides per protein",
                     #limits = 
    )+
    scale_y_continuous(expand = expansion(mult = c(0.01, 0.1), add = c(0, 0)))+
    theme_bw()+
    theme(legend.text = element_text(size = 18, color = "black"),legend.position = 'top',
          legend.title = element_text(size = 18, color="black") ,
          panel.grid.major =element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+
    theme(panel.grid =element_blank())+
    theme(axis.text = element_text(size = 15,color = "black"))+
    theme(axis.text.x = element_text(size = 15,color = "black", angle = 0))+
    theme(plot.title = element_text(size = 20, face = "bold"))+
    theme(plot.subtitle=element_text(size = 20, hjust = 0, color="black"))+
    theme(axis.title.x=element_text(size = 18, hjust = 0.5, color="black"))+
    theme(axis.title.y=element_text(size = 18, hjust = 0.5, color="black"))
  
  
  df$Precursor <- paste0(df$ModifiedPeptideSequence, '_', df$'PrecursorCharge')
  df_frpr <- df %>% dplyr::distinct(Precursor, Annotation)
  df_frpr <- data.frame(table(table(df_frpr$Precursor))); colnames(df_frpr) <- c('fr_num', 'count')
  if(sum(as.numeric(as.character(df_frpr$fr_num)) >= 6)){
    df_frpr$fr_num <- as.numeric(as.character(df_frpr$fr_num))
    over6 <- sum(df_frpr[df_frpr$fr_num >= 6, 'count'])
    df_frpr <- df_frpr[df_frpr$fr_num < 6, ]
    df_frpr <-  rbind(df_frpr, c('>=6', over6))
    df_frpr$fr_num <- factor(df_frpr$fr_num, levels = c('1', '2', '3', '4', '5', '>=6'), ordered = T)
    df_frpr$count <- as.numeric(df_frpr$count)
  }
  df_frpr$freq <- df_frpr$count / sum(df_frpr$count)
  # Plot G: Fragments per precursor
  p7 <- ggplot(df_frpr, aes(fr_num, freq))+
    geom_col(width = 0.5, color = allColor[7], fill = allFill[7])+
    geom_text(aes(label = paste0(sprintf('%.3f', freq * 100), '%')), size = 6, hjust = 0.5, vjust = -0.2, position = "stack", color = 'black')+
    labs(title = 'Fragments per precursor', x = 'Fragments per precursor', y = 'Frequency in percent (%)')+
    scale_x_discrete("Fragments per precursor",
                     #limits = 
    )+
    scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0.01, 0.1), add = c(0, 0)))+
    theme_bw()+
    theme(legend.text = element_text(size = 18, color = "black"),legend.position = 'top',
          legend.title = element_text(size = 18, color="black") ,
          panel.grid.major =element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+
    theme(panel.grid =element_blank())+
    theme(axis.text = element_text(size = 15,color = "black"))+
    theme(axis.text.x = element_text(size = 15,color = "black", angle = 0))+
    theme(plot.title = element_text(size = 20, face = "bold"))+
    theme(plot.subtitle=element_text(size = 20, hjust = 0, color="black"))+
    theme(axis.title.x=element_text(size = 18, hjust = 0.5, color="black"))+
    theme(axis.title.y=element_text(size = 18, hjust = 0.5, color="black"))
  
  
  df_frtype <- df %>% dplyr::distinct(Precursor, Annotation, FragmentType)
  df_frtype <- data.frame(table(df_frtype$FragmentType)); colnames(df_frtype) <- c('fr_type', 'count')
  df_frtype$freq <- df_frtype$count / sum(df_frtype$count)
  # Plot H: Fragment ion type distribution
  p8 <- ggplot(df_frtype, aes(fr_type, freq))+ theme_bw()+
    geom_col(color = allColor[8], fill = allFill[8], width = 0.3)+
    geom_text(aes(label = paste0(sprintf('%.2f', freq * 100), '%')),
              size = 6, hjust = 0.5, vjust = -0.2, position = "stack", color = 'black')+
    labs(title = 'Fragment ion', x = 'Fragment ion type', y = 'Frequency in percent (%)')+
    scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0.01, 0.1), add = c(0, 0)))+
    theme_bw()+
    theme(legend.text = element_text(size = 18, color = "black"),legend.position = 'top',
          legend.title = element_text(size = 18, color="black") ,
          panel.grid.major =element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+
    theme(panel.grid =element_blank())+
    theme(axis.text = element_text(size = 15,color = "black"))+
    theme(axis.text.x = element_text(size = 15,color = "black", angle = 0))+
    theme(plot.title = element_text(size = 20, face = "bold"))+
    theme(plot.subtitle=element_text(size = 20, hjust = 0, color="black"))+
    theme(axis.title.x=element_text(size = 18, hjust = 0.5, color="black"))+
    theme(axis.title.y=element_text(size = 18, hjust = 0.5, color="black"))
  
  
  df_frz <- df %>% dplyr::distinct(Precursor, Annotation, FragmentCharge)
  df_frz <- data.frame(table(df_frz$FragmentCharge)); colnames(df_frz) <- c('fr_z', 'count')
  df_frz$freq <- df_frz$count / sum(df_frz$count)
  # Plot I: Fragment ion charge distribution
  p9 <- ggplot(df_frz, aes(fr_z, freq))+ theme_bw()+
    geom_col(color = allColor[9], fill = allFill[9], width = 0.5)+
    geom_text(aes(label = paste0(sprintf('%.2f', freq * 100), '%')),
              size = 6, hjust = 0.5, vjust = -0.2, position = "stack", color = 'black')+
    labs(title = 'Fragment ion charge', x = 'Fragment ion charge', y = 'Frequency in percent (%)')+
    scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0.01, 0.1), add = c(0, 0)))+
    theme_bw()+
    theme(legend.text = element_text(size = 18, color = "black"),legend.position = 'top',
          legend.title = element_text(size = 18, color="black") ,
          panel.grid.major =element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+
    theme(panel.grid =element_blank())+
    theme(axis.text = element_text(size = 15,color = "black"))+
    theme(axis.text.x = element_text(size = 15,color = "black", angle = 0))+
    theme(plot.title = element_text(size = 20, face = "bold"))+
    theme(plot.subtitle=element_text(size = 20, hjust = 0, color="black"))+
    theme(axis.title.x=element_text(size = 18, hjust = 0.5, color="black"))+
    theme(axis.title.y=element_text(size = 18, hjust = 0.5, color="black"))
  
  
  print('Calculating protein sequence coverage which will cost several minutes')
  fastafile <- seqinr::read.fasta(file = fasPath, seqtype = 'AA', as.string = TRUE)
  df_fas <- data.frame(protein = names(fastafile) %>% stringr::str_split(pattern = '\\|') %>% sapply(., function(e){e[2]}), sequence = unlist(fastafile))
  prots <- sort(unique(df$ProteinId))
  
  # 20,000 proteins cost 13 minutes
  coverage <- rep(NaN, length(prots))
  for(i in seq_along(prots)){
    prot <- prots[i]
    peptides <- unique(df$PeptideSequence[df$ProteinId == prot])
    full_seq <- df_fas$sequence[df_fas$protein == prot]
    covered <- rep(0, nchar(full_seq))# recording matched peptide sequence with number 1
    loc <- na.omit(stringr::str_locate(full_seq, peptides))
    if(nrow(loc) == 0) next;
    loc <- lapply(seq_len(nrow(loc)), function(j) {loc[j, 1]:loc[j, 2]} )
    loc <- Reduce(base::union, loc)
    covered[loc] <- 1
    
    coverage[i] <- sum(covered) / length(covered)
  }
  
  # Plot J: protein sequence coverage
  p10 <- ggplot(data.frame(coverage = coverage), aes(x = coverage)) + 
    geom_histogram(aes(y = (..count..)/sum(..count..)),
                   color = allColor[10],
                   fill = allFill[10])+
    #geom_density(alpha = .2, fill = "#FF6666")+
    labs(title = 'Protein coverage', x = "Protein coverage", y = "Frequency in percent (%)")+
    scale_y_continuous(expand = expansion(mult = c(0.01, 0.1), add = c(0, 0)))+
    theme_bw()+
    theme(legend.text = element_text(size = 18, color = "black"),legend.position = 'top',
          legend.title = element_text(size = 18, color="black") ,
          panel.grid.major =element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+
    theme(panel.grid =element_blank())+
    theme(axis.text = element_text(size = 15,color = "black"))+
    theme(axis.text.x = element_text(size = 15,color = "black", angle = 0))+
    theme(plot.title = element_text(size = 20, face = "bold"))+
    theme(plot.subtitle=element_text(size = 20, hjust = 0, color="black"))+
    theme(axis.title.x=element_text(size = 18, hjust = 0.5, color="black"))+
    theme(axis.title.y=element_text(size = 18, hjust = 0.5, color="black"))
  
  
  print('Calculating missed cleavage')
  df_peps <- data.frame(PeptideSequence = unique(df$PeptideSequence))
  df_peps$MissedCleavage <- 0
  df_peps$MissedCleavage <- apply(df_peps, 1, function(row){
    Seq <- row['PeptideSequence']
    # exclude the last cleavage site
    if (stringr::str_sub(Seq, -1, -1) == 'P' & stringr::str_sub(Seq, -2, -2) %in% c('K', 'R')){
      Seq <- stringr::str_sub(Seq, 1, -3)
    }else{
      Seq <- stringr::str_sub(Seq, 1, -2)
    }
    #
    Seqs <- stringr::str_sub(Seq, 1:nchar(Seq), 1:nchar(Seq))
    is_KR <- Seqs %in% c('K', 'R')
    isnt_P <- Seqs != 'P'
    flag <- c(isnt_P[2:length(isnt_P)], T)
    mis_cleav <- is_KR & flag
    return(sum(mis_cleav))
  })
  # Plot K: missed cleavage
  df_mis <- as.data.frame(table(df_peps$MissedCleavage)) %>% apply(., c(1, 2), as.character)
  vec_mis <- c()
  for(i in 1:nrow(df_mis)){
    vec_mis <- append(vec_mis, rep(df_mis[i, 1], df_mis[i, 2]))
  }
  p11 <- ggplot(data.frame(V1 = vec_mis), aes(x=`V1`)) +
    geom_bar(aes(y = (..count..)/sum(..count..)),
             width = 0.6,
             position="stack",
             color = allColor[11],
             fill = allFill[11],
    )+
    geom_text(aes(y = ((..count..)/sum(..count..)), label = scales::percent((..count..)/sum(..count..))), stat = 'count', size = 6, hjust = 0.5, vjust = -0.2, position = "stack", color = 'black')+
    scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0.01, 0.1), add = c(0, 0)))+
    labs(title = 'Missed cleavage', x = "Missed cleavage", y = "Frequency in percent (%)")+
    theme_bw()+
    theme(legend.text = element_text(size = 18, color = "black"),legend.position = 'top',
          legend.title = element_text(size = 18, color="black") ,
          panel.grid.major =element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+
    theme(panel.grid =element_blank())+
    theme(axis.text = element_text(size = 15,color = "black"))+
    theme(axis.text.x = element_text(size = 15,color = "black", angle = 0))+
    theme(plot.title = element_text(size = 20, face = "bold"))+
    theme(plot.subtitle=element_text(size = 20, hjust = 0, color="black"))+
    theme(axis.title.x=element_text(size = 18, hjust = 0.5, color="black"))+
    theme(axis.title.y=element_text(size = 18, hjust = 0.5, color="black"))
  
  pdf(file = paste0(rlt_dir, '/', pdfName), width = 20, height = 16)
  p <- ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11 , labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K"), font.label = list(size = 14, face = "bold"), ncol = 3, nrow = 4)
  print(p)
  while (length(dev.list()) > 0){
    dev.off()
  }
}


# --------- main DIA-LibQC -------------
# read library
#libPath <- 'frag_755_final.tsv'
libPath <- '//172.16.13.136/TPHP/library/TPHPlib_frag1025_swissprot_final.tsv'
df <- data.table::fread(libPath, sep = '\t', stringsAsFactors = F, check.names = F, data.table = F)
# read fasta
fasPath <- '//172.16.13.114/share/members/jiangwenhao/fasta/swissprot-homo+sapiens-20200507_20377_iRT_11.fasta'
# set color style
allColors <- c('#639d98', '#1e5d64', '#203d26', '#b49259', '#7d2a16')
allFills <- c('#75bcb6', '#2b828c', '#2f5b38', '#e0b66e', '#a6381d')
image(x = 1:5, y = 1, z = as.matrix(1:5), col = colorRampPalette(allFills)(5))

plotting(df, fasPath = fasPath, pdfName = 'frag1025_swissprot_libqc_20220722.pdf',
         allColor = allColors[c(1, 2, 3, 4, 5, 2, 3, 1, 4, 2, 4)],
         allFill = allFills[c(1, 2, 3, 4, 5, 2, 3, 1, 4, 2, 4)])

# # 200K isoform library
# libPath <- '//172.16.13.136/TPHP/TPL/libs/frag755_200k_isoform_v4.2.tsv'
# df <- data.table::fread(libPath, sep = '\t', stringsAsFactors = F, check.names = F, data.table = F)
# # read fasta
# fasPath <- '//172.16.13.136/TPHP/TPL/required_files/uniprot_human_200k_iRT_20220119.fasta'
# plotting(df, fasPath = fasPath,
#          rlt_dir = '//172.16.13.136/share/members/jiangwenhao/TPHP_git/supplementary/LibQC/',
#          pdfName = 'frag755_200K_isoform_libqc_20220615.pdf',
#          allColor = allColors[c(1, 2, 3, 4, 5, 2, 3, 1, 4, 2, 4)],
#          allFill = allFills[c(1, 2, 3, 4, 5, 2, 3, 1, 4, 2, 4)])

# --------- comparison with some database -------------
atlas_pep <- read.delim('//172.16.13.136/share/members/jiangwenhao/TPHP/DIALibQC_plot/input/peptideAtlas/atlas_tables_502-peptide.tsv', stringsAsFactors = F, check.names = F)$peptide_sequence %>% unique
atlas_prot <- read.csv('//172.16.13.136/share/members/jiangwenhao/TPHP/DIALibQC_plot/input/peptideAtlas/query_guest_20210705-051831.csv', stringsAsFactors = F, check.names = F)$nextprot_accession %>% unique

df_protdb <- readxl::read_excel('//172.16.13.136/share/members/jiangwenhao/TPHP/DIALibQC_plot/input/proteomicsDB/41586_2014_BFnature13319_MOESM94_ESM.xlsx', sheet = 2)[, c('Unique Identifier', 'Sequence')]

# TPL - PeptideAtlas
pep_list <- list(TPL = unique(df$PeptideSequence), PeptideAtlas = atlas_pep)
prot_list <- list(TPL = unique(df$ProteinId), PeptideAtlas = atlas_prot)
prot_app <- setdiff(prot_list$TPL, prot_list$PeptideAtlas)
writexl::write_xlsx(data.frame(protein = prot_app), '//172.16.13.136/share/members/jiangwenhao/TPHP/DIALibQC_plot/output/TPL_over_peptideatlas.xlsx')

image(x = 1:8, y = 1, z = as.matrix(1:8), col = brewer.pal(8,"Pastel2"))
allColors <- brewer.pal(8,"Pastel2")
p1 <- VennDiagram::venn.diagram(x = pep_list,
                                resolution = 300,
                                alpha=rep(0.5, length(pep_list)),
                                fill=allColors[c(2, 5)], 
                                main="Peptides",
                                #sub = rep,
                                main.cex = 4, 
                                sub.cex = 3,
                                cex = 4,
                                cex.lab=4,
                                cat.cex=4,
                                imagetype = "tiff", 
                                filename = NULL
)
p2 <- VennDiagram::venn.diagram(x = prot_list,
                                resolution = 300,
                                alpha=rep(0.5, length(prot_list)),
                                fill=allColors[c(2, 5)], 
                                main="Proteins",
                                #sub = rep,
                                main.cex = 4, 
                                sub.cex = 3,
                                cex = 4,
                                cex.lab=4,
                                cat.cex=4,
                                imagetype = "tiff", 
                                filename = NULL
)
pdf('//172.16.13.136/share/members/jiangwenhao/TPHP/DIALibQC_plot/output/TPL_peptideAtlas_VennDiagram_20220722.pdf', width = 20, height = 20)
grid::grid.newpage(); grid::grid.draw(p1)
grid::grid.newpage(); grid::grid.draw(p2)
dev.off()

# TPL - ProteomicsDB
pep_list <- list(TPL = unique(df$PeptideSequence), ProteomicsDB = df_protdb$Sequence)
prot_list <- list(TPL = unique(df$ProteinId), ProteomicsDB = df_protdb$`Unique Identifier`)
prot_app <- setdiff(prot_list$TPL, prot_list$ProteomicsDB)
writexl::write_xlsx(data.frame(protein = prot_app), '//172.16.13.136/share/members/jiangwenhao/TPHP/DIALibQC_plot/output/TPL_over_proteomicsdb_20220722.xlsx')
p3 <- VennDiagram::venn.diagram(x = pep_list,
                                resolution = 300,
                                alpha=rep(0.5, length(pep_list)),
                                fill=allColors[c(1, 5)], 
                                main="Peptides",
                                #sub = rep,
                                main.cex = 4, 
                                sub.cex = 3,
                                cex = 4,
                                cex.lab=4,
                                cat.cex=4,
                                imagetype = "tiff", 
                                filename = NULL
)
p4 <- VennDiagram::venn.diagram(x = prot_list,
                                resolution = 300,
                                alpha=rep(0.5, length(prot_list)),
                                fill=allColors[c(1, 5)], 
                                main="Proteins",
                                #sub = rep,
                                main.cex = 4, 
                                sub.cex = 3,
                                cex = 4,
                                cex.lab=4,
                                cat.cex=4,
                                imagetype = "tiff", 
                                filename = NULL
)
pdf('//172.16.13.136/share/members/jiangwenhao/TPHP/DIALibQC_plot/output/TPL_proteomicsDB_VennDiagram_20220722.pdf', width = 12, height = 12)
grid::grid.newpage(); grid::grid.draw(p3)
grid::grid.newpage(); grid::grid.draw(p4)
dev.off()


# TPL - PeptideAtlas - ProteomicsDB
pep_list <- list(PUH = unique(df$PeptideSequence),
                 PeptideAtlas = atlas_pep,
                 ProteomicsDB = df_protdb$Sequence)
prot_list <- list(PUH = unique(df$ProteinId),
                  PeptideAtlas = atlas_prot,
                  ProteomicsDB = df_protdb$`Unique Identifier`)

p1 <- VennDiagram::venn.diagram(x = pep_list,
                                resolution = 300,
                                #alpha=rep(0.5, length(pep_list)),
                                fill=c("#3C65D6", "#2DCC62FD", "#EB495C"), 
                                main="Peptides",
                                #sub = rep,
                                main.cex = 4, 
                                sub.cex = 3,
                                cex = 4,
                                cex.lab=4,
                                cat.cex=4,
                                imagetype = "tiff", 
                                filename = NULL,
                                width = 5000, height = 5000
)
p2 <- VennDiagram::venn.diagram(x = prot_list,
                                resolution = 300,
                                #alpha=rep(0.5, length(pep_list)),
                                fill=c("#3C65D6", "#2DCC62FD", "#EB495C"), 
                                main="Proteins",
                                #sub = rep,
                                main.cex = 4, 
                                sub.cex = 3,
                                cex = 4,
                                cex.lab=4,
                                cat.cex=4,
                                imagetype = "tiff", 
                                filename = NULL,
                                width = 5000, height = 5000
)


pdf('//172.16.13.136/share/members/jiangwenhao/TPHP/DIALibQC_plot/output/PUH_peptideAtlas_ProteomicsDB_VennDiagram_20220722.pdf', width = 20, height = 20)
grid::grid.newpage(); grid::grid.draw(p1)
grid::grid.newpage(); grid::grid.draw(p2)
dev.off()
