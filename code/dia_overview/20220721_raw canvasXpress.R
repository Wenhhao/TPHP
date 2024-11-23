# 20210623 prepare date for canvasXpress
# order data and prepare data for jiangwenhao
# first order: DDA_lib_type; second order: cancer_type
# library(xlsx)

setwd("//172.16.13.136/share/members/yuel")
source("//172.16.13.136/share/members/jiangwenhao/code/myQC.R")

library(readxl)
library(writexl)
library(magrittr)
library(tidyverse)
# use all files types
df1<-readxl::read_excel('//172.16.13.114/share/members/jiangwenhao/TPHP/preparation/20220701TPHP_1781file_117pool_13478prot_info.xlsx',
                       col_types = c(rep('guess', 16), rep('numeric', 13477)))


dim(df1) # 1898 13493

df1 %<>% filter(!str_detect(cancer_type, "pool")) # remove pool
df1[, -(1:16)] %>% dim() # 1781 13477
df1[, -(1:16)] %>%
  removeColsAllNa() %>%
  dim() # 1781 13477
dim(df1) # 1781 13493
df2 <- df1[, 1:16]


# refine the information table
# library(dplyr)
df3 <- df2
df3$cancer_type[which(df2$cancer_type == "Normal")] <- -10000
df3$cancer_type[which(df2$cancer_type == "adjacent")] <- 30000
df3$cancer_type[which(df2$cancer_type == "carcinoma")] <- 70000
df3 <- arrange(df3, DDA_lib_type, cancer_type)
df3 %<>% filter(!is.na(DDA_lib_type))
write_xlsx(df3, "sample_info/figure/F1B/20220721TPHP_canvas_plot.xlsx")

datay <- as.data.frame(t(df3[, c("Peptides", "Proteins", "cancer_type", "DDA_lib_pro")]))
colnames(datay) <- df3$FileName
write.csv(datay, "sample_info/figure/F1B/datay_v2_20220721.csv")

datax <- as.data.frame(df3[, "DDA_lib_type"])
names(datax) <- c("tissue_type")
rownames(datax) <- df3$FileName
write.csv(datax, "sample_info/figure/F1B/datax_v2_20220721.csv")

dataz <- as.data.frame(matrix(nrow = 4, ncol = 1))
rownames(dataz) <- c("Peptides", "Proteins", "cancer_type", "DDA_lib_pro")
colnames(dataz) <- "ring"
dataz$ring <- c(1, 2, 3, 4)

if (!require("jsonlite")) install.packages("jsonlite", update = F, ask = F)
library(canvasXpress)

setwd("//172.16.13.136/share/members/yuel/sample_info/figure/F1B")

canvasXpress(
  data = datay,
  smpAnnot = datax,
  varAnnot = dataz,
  circularArc = 340, # remain blank sector of 20 degree
  circularType = "radar",
  graphType = "Circular",
  legendPosition = "top",
  ringGraphType = list("bar", "bar", "heatmap", "bar"),
  smpDendrogramPosition = "ins",
  arcSegmentsSeparation = 3,
  showTransition = FALSE,
  showYaxis = FALSE,
  ringGraphWeight = list(22, 22, 3, 53), # rings from outside to inside
  ringsOrder = list("labels", "overlays", "dendrogram", "data"),
  segregateVariablesBy = list("ring"),
  segregateSamplesBy = list("tissue_type"),
  title = "Overview of TPHP proteome data",
  smpOverlays = list("tissue_type", "-", "cancer_type"),
  colorScheme = "Tableau",
  circularAnchors2Align = "inside",
  showSampleNames = FALSE
)

datax %>%
  count(tissue_type) %>%
  arrange(n) %>%
  writexl::write_xlsx('TPHP_canvas_DDALibType_count.xlsx')


