library(tidyverse)
library(magrittr)
setwd('Y:/members/jiangwenhao/TPHP/20220908/F2_sample_compare')

df <- rio::import('20230710_tissue_comparison_edited.xlsx')
df_in <- df %>% select(Inner_ring_label, Inner_class) %>% distinct()
df_out <- df %>% select(Detailed_tissue_type, Outer_class) %>% distinct()


df_in %>% count(Inner_class)
in_PUH <- df_in$Inner_ring_label
in_Jiang <- df_in %>% filter(Inner_class %in% c('Included in all', 'Overlapped with Jiang')) %>% pull(Inner_ring_label)
in_HPA <- df_in %>% filter(Inner_class %in% c('Included in all', 'Overlapped with HPA')) %>% pull(Inner_ring_label)

df_out %>% count(Outer_class)
out_PUH <- df_out$Detailed_tissue_type
out_Jiang <- df_out %>% filter(Outer_class %in% c('Included in all', 'Overlapped with Jiang')) %>% pull(Detailed_tissue_type)
out_HPA <- df_out %>% filter(Outer_class %in% c('Included in all', 'Overlapped with HPA')) %>% pull(Detailed_tissue_type)

# # devtools::install_github("vqf/nVennR")
# library(nVennR)
# inVenn <- plotVenn(list(PUH = in_PUH, Jiang = in_Jiang, HPA = in_HPA))
# showSVG(inVenn, opacity = 0.9, setColors=c('#AD1E26', '#4679AC', '#CBA4CB'), fontScale = 2, borderWidth = 0, outFile = 'Tissue_comparison_inner_PUH_Jiang_HPA_20230721.svg')
# 
# getVennRegion(inVenn, c('PUH', 'HPA'))
# getVennRegion(inVenn, c('PUH', 'Jiang'))
# getVennRegion(inVenn, c('PUH', 'HPA', 'HPA'))

p1 <- VennDiagram::venn.diagram(x = list(PUH = in_PUH, Jiang = in_Jiang, HPA = in_HPA),
                                resolution = 300,
                                col=c('#AD1E26', '#4679AC', '#CBA4CB'), 
                                fill=c('#AD1E26', '#4679AC', '#CBA4CB'), 
                                alpha=rep(0.1, 3),
                                main="Inner ring comparison",
                                #sub = rep,
                                main.cex = 4, 
                                sub.cex = 3,
                                cex = 4,
                                cex.lab=4,
                                cat.cex=4,
                                imagetype = "tiff", 
                                filename = NULL,disable.logging = T
)

p2 <- VennDiagram::venn.diagram(x = list(PUH = out_PUH, Jiang = out_Jiang, HPA = out_HPA),
                                resolution = 300,
                                col=c('#AD1E26', '#4679AC', '#CBA4CB'), 
                                fill=c('#AD1E26', '#4679AC', '#CBA4CB'), 
                                alpha=rep(0.1, 3),
                                main="Outer ring comparison",
                                #sub = rep,
                                main.cex = 4, 
                                sub.cex = 3,
                                cex = 4,
                                cex.lab=4,
                                cat.cex=4,
                                imagetype = "tiff", 
                                filename = NULL,disable.logging = T
)


pdf('Tissue_comparison_PUH_Jiang_HPA_20230721.pdf', width = 10, height = 10)
grid::grid.newpage(); grid::grid.draw(p1)
grid::grid.newpage(); grid::grid.draw(p2)
graphics.off()


