# set working directory
setwd("~/Work/PhD/Projects/Comput MM Pain/Computational Analysis/RRHO")

# load packages
library(RRHO)
library(RRHO2)
library(lattice)
library(tidyverse)
library(RColorBrewer)
library(viridisLite)
library(circlize)
library(cowplot)
library(gridExtra)

# color gradient

#col <- colorRampPalette(rev(brewer.pal(8, "RdYlBu")))(25)
#col <- viridis(100)

#col2 <- colorRampPalette(c("#ffffff","#0066ff","#0000ff","#6ae26e","#28CC2D","#1d9521","#FFF44F","#FFE135","#ff9900","#ff6666", "#ff0000"), bias > 1)(100)

#col3 <- colorRampPalette(c("#ee82ee","#0000ff","#00ffff","#008000","#ffff00","#ffa500","#ff0000")) (370)

#col4 <- colorRampPalette(c('white','#D6F8F7',"#BEDAFF",'#5DA4FF',"#0000FF","#D4F9E2","#00FF7F",
#                           "#008000","#FFFF00","#FFD27F", "#FFB732" ,"#EE7600", "#D53E4F","#FF6A6A", "gray")) (400)

col_fun <- colorRamp2(c(0, 25, 50, 75,100), c("#440154", "#3b528b", "#21918c", "#5ec962", "#fde725"))

#lattice.options(
#  layout.heights=list(bottom.padding=list(x=0), top.padding=list(x=0)),
#  layout.widths=list(left.padding=list(x=0), right.padding=list(x=0)))
#)

# load data and prepare them for rrho package
# For up-regulated genes, the input score should be calculated using-log10(pvalue) * 1;
# For down-regulated genes, the input score should be calculated using-log10(pvalue) * (-1);

#_________________________________________________________________________________________________
# load data

## plasma cells
pc_data <- read.csv("data/pc/PC_MA.csv")
pc_data <- pc_data[,c(2,3,12,15)]
pc_data <- transform(pc_data, fisherPvalMin = pmin(fisherPvalUp, fisherPvalDown))
pc_data$score <- with(pc_data, ifelse(effectSize < 0, -log10(fisherPvalMin)*-1, -log10(fisherPvalMin)*1))
pc_data <- pc_data[,c(1,6)]
colnames(pc_data)[1] <- "gene_symbol"
pc_data <- pc_data[complete.cases(pc_data), ]

## bone marrow stromal cells
bmsc_data <- read.csv("data/bmsc/bmsc_MA.csv")
bmsc_data <- bmsc_data[,c(2,3,12,15)]
bmsc_data <- transform(bmsc_data, fisherPvalMin = pmin(fisherPvalUp, fisherPvalDown))
bmsc_data$score <- with(bmsc_data, ifelse(effectSize < 0, -log10(fisherPvalMin)*-1, -log10(fisherPvalMin)*1))
bmsc_data <- bmsc_data[,c(1,6)]
colnames(bmsc_data)[1] <- "gene_symbol"
bmsc_data <- bmsc_data[complete.cases(bmsc_data), ]

## endothelial cells
ec_data <- read.csv("data/ec/GSE14230 DEG.ID.csv")
ec_data <- ec_data[,c(3,6,9)]
ec_data$score <- with(ec_data, ifelse(logFC < 0, -log10(P.Value)*-1, -log10(P.Value)*1))
ec_data <- ec_data[,c(1,4)]
colnames(ec_data)[1] <- "gene_symbol"
ec_data <- ec_data[complete.cases(ec_data), ]

## osteocyte
ocy_data <- read.csv("data/ocy/GSE27372 DEG.ID.csv")
ocy_data <- ocy_data[,c(3,6,9)]
ocy_data$score <- with(ocy_data, ifelse(logFC < 0, -log10(P.Value)*-1, -log10(P.Value)*1))
ocy_data <- ocy_data[,c(1,4)]
colnames(ocy_data)[1] <- "gene_symbol"
ocy_data <- ocy_data[complete.cases(ocy_data), ]

## osteogenic precursor cell
opc_data <- read.csv("data/opc/GSE87073 DEG.ID.csv")
opc_data <- opc_data[,c(3,6,9)]
opc_data$score <- with(opc_data, ifelse(logFC < 0, -log10(P.Value)*-1, -log10(P.Value)*1))
opc_data <- opc_data[,c(1,4)]
colnames(opc_data)[1] <- "gene_symbol"
opc_data <- opc_data[complete.cases(opc_data), ]

## adipocyte
adcy_data <- read.csv("data/adcy/GSE143269_DEG.ID.HS.Mapped.csv")
adcy_data <- adcy_data[,c(4,7,10)]
adcy_data$score <- with(adcy_data, ifelse(logFC < 0, -log10(P.Value)*-1, -log10(P.Value)*1))
adcy_data <- adcy_data[,c(3,4)]
colnames(adcy_data)[1] <- "gene_symbol"
adcy_data <- adcy_data[complete.cases(adcy_data), ]
sum(duplicated(adcy_data$gene_symbol))
adcy_data <- data.frame(adcy_data %>% 
                          group_by(gene_symbol) %>% 
                          summarise_all(funs(max)))

## natural killer cell
nkc_data <- read.csv("data/nkc/GSE27838 DEG.ID.csv")
nkc_data <- nkc_data[,c(3,6,9)]
nkc_data$score <- with(nkc_data, ifelse(logFC < 0, -log10(P.Value)*-1, -log10(P.Value)*1))
nkc_data <- nkc_data[,c(1,4)]
colnames(nkc_data)[1] <- "gene_symbol"
nkc_data <- nkc_data[complete.cases(nkc_data), ]

## neutrophil
neut_data <- read.csv("data/neut/GSE150021_DEG_limma.csv")
neut_data <- neut_data[,c(2,3,6)]
neut_data$score <- with(neut_data, ifelse(logFC < 0, -log10(P.Value)*-1, -log10(P.Value)*1))
neut_data <- neut_data[,c(1,4)]
colnames(neut_data)[1] <- "gene_symbol"
neut_data <- neut_data[complete.cases(neut_data), ]
colnames(neut_data)[1] <- "gene_symbol"
sum(duplicated(neut_data$gene_symbol))
neut_data <- data.frame(neut_data %>% 
                          group_by(gene_symbol) %>% 
                          summarise_all(funs(max)))

## macrophages
mp_data <- read.csv("data/mp/MP_MA.csv")
mp_data <- mp_data[,c(2,4,13,16)]
mp_data <- transform(mp_data, fisherPvalMin = pmin(fisherPvalUp, fisherPvalDown))
mp_data$score <- with(mp_data, ifelse(effectSize < 0, -log10(fisherPvalMin)*-1, -log10(fisherPvalMin)*1))
mp_data <- mp_data[,c(1,6)]
colnames(mp_data)[1] <- "gene_symbol"
mp_data <- mp_data[complete.cases(mp_data), ]
sum(duplicated(mp_data$gene_symbol))
mp_data <- data.frame(mp_data %>% 
                          group_by(gene_symbol) %>% 
                          summarise_all(funs(max)))

## regulatory t cell
treg_data <- read.csv("data/treg/GSE109533_DEG_limma.csv")
treg_data <- treg_data[,c(4,7,11)]
treg_data$score <- with(treg_data, ifelse(logFC < 0, -log10(P.Value)*-1, -log10(P.Value)*1))
treg_data <- treg_data[,c(3,4)]
colnames(treg_data)[1] <- "gene_symbol"
treg_data <- treg_data[complete.cases(treg_data), ]
sum(duplicated(treg_data$gene_symbol))
treg_data <- data.frame(treg_data %>% 
                        group_by(gene_symbol) %>% 
                        summarise_all(funs(max)))

## hematopoietic stem cell
hsc_data <- read.csv("data/hsc/GSE24870_HSC DEG.ID.csv")
hsc_data <- hsc_data[,c(3,6,9)]
hsc_data$score <- with(hsc_data, ifelse(logFC < 0, -log10(P.Value)*-1, -log10(P.Value)*1))
hsc_data <- hsc_data[,c(1,4)]
colnames(hsc_data)[1] <- "gene_symbol"
hsc_data <- hsc_data[complete.cases(hsc_data), ]

## whole bone marrow
wbm_data <- read.csv("data/wbm/GSE118985 DEG.ID.SYMBOL.csv")
wbm_data <- wbm_data[,c(4,21,24)]
wbm_data$score <- with(wbm_data, ifelse(logFC < 0, -log10(P.Value)*-1, -log10(P.Value)*1))
wbm_data <- wbm_data[,c(1,4)]
colnames(wbm_data)[1] <- "gene_symbol"
wbm_data <- wbm_data[complete.cases(wbm_data), ]
sum(duplicated(wbm_data$gene_symbol))
wbm_data <- data.frame(wbm_data %>% 
                        group_by(gene_symbol) %>% 
                        summarise_all(funs(max)))

#_________________________________________________________________________________________________

# take intersect

## take plasma cell intersect
pc_bmsc_intersect <- intersect(pc_data$gene_symbol,bmsc_data$gene_symbol)
pc_bmsc <- pc_data[pc_data$gene_symbol %in% pc_bmsc_intersect, ]
bmsc_pc <- bmsc_data[bmsc_data$gene_symbol %in% pc_bmsc_intersect, ]

pc_ec_intersect <- intersect(pc_data$gene_symbol,ec_data$gene_symbol)
pc_ec <- pc_data[pc_data$gene_symbol %in% pc_ec_intersect, ]
ec_pc <- ec_data[ec_data$gene_symbol %in% pc_ec_intersect, ]

pc_ocy_intersect <- intersect(pc_data$gene_symbol,ocy_data$gene_symbol)
pc_ocy <- pc_data[pc_data$gene_symbol %in% pc_ocy_intersect, ]
ocy_pc <- ocy_data[ocy_data$gene_symbol %in% pc_ocy_intersect, ]

pc_opc_intersect <- intersect(pc_data$gene_symbol,opc_data$gene_symbol)
pc_opc <- pc_data[pc_data$gene_symbol %in% pc_opc_intersect, ]
opc_pc <- opc_data[opc_data$gene_symbol %in% pc_opc_intersect, ]

pc_adcy_intersect <- intersect(pc_data$gene_symbol,adcy_data$gene_symbol)
pc_adcy <- pc_data[pc_data$gene_symbol %in% pc_adcy_intersect, ]
adcy_pc <- adcy_data[adcy_data$gene_symbol %in% pc_adcy_intersect, ]

pc_nkc_intersect <- intersect(pc_data$gene_symbol,nkc_data$gene_symbol)
pc_nkc <- pc_data[pc_data$gene_symbol %in% pc_nkc_intersect, ]
nkc_pc <- nkc_data[nkc_data$gene_symbol %in% pc_nkc_intersect, ]

pc_neut_intersect <- intersect(pc_data$gene_symbol,neut_data$gene_symbol)
pc_neut <- pc_data[pc_data$gene_symbol %in% pc_neut_intersect, ]
neut_pc <- neut_data[neut_data$gene_symbol %in% pc_neut_intersect, ]

pc_mp_intersect <- intersect(pc_data$gene_symbol,mp_data$gene_symbol)
pc_mp <- pc_data[pc_data$gene_symbol %in% pc_mp_intersect, ]
mp_pc <- mp_data[mp_data$gene_symbol %in% pc_mp_intersect, ]

pc_treg_intersect <- intersect(pc_data$gene_symbol,treg_data$gene_symbol)
pc_treg <- pc_data[pc_data$gene_symbol %in% pc_treg_intersect, ]
treg_pc <- treg_data[treg_data$gene_symbol %in% pc_treg_intersect, ]

pc_hsc_intersect <- intersect(pc_data$gene_symbol,hsc_data$gene_symbol)
pc_hsc <- pc_data[pc_data$gene_symbol %in% pc_hsc_intersect, ]
hsc_pc <- hsc_data[hsc_data$gene_symbol %in% pc_hsc_intersect, ]

pc_wbm_intersect <- intersect(pc_data$gene_symbol,wbm_data$gene_symbol)
pc_wbm <- pc_data[pc_data$gene_symbol %in% pc_wbm_intersect, ]
wbm_pc <- wbm_data[wbm_data$gene_symbol %in% pc_wbm_intersect, ]


#_________________________________________________________________________________________________

## compute plasma cell overlap  significance
rrho_obj_pc_bmsc <- RRHO(pc_bmsc, bmsc_pc, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_pc_bmsc <-  RRHO2_initialize(pc_bmsc, bmsc_pc, log10.ind=TRUE)
save(rrho_obj_pc_bmsc, file = "data/rrho_obj_pc_bmsc.RData")
load("data/pc/rrho_obj_pc_bmsc.RData")

rrho_obj_pc_ec <- RRHO(pc_ec, ec_pc, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_pc_ec <-  RRHO2_initialize(pc_ec, ec_pc, log10.ind=TRUE)
save(rrho_obj_pc_ec, file = "data/rrho_obj_pc_ec.RData")
load("data/pc/rrho_obj_pc_ec.RData")

rrho_obj_pc_ocy <- RRHO(pc_ocy, ocy_pc, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_pc_ocy <-  RRHO2_initialize(pc_ocy, ocy_pc, log10.ind=TRUE)
save(rrho_obj_pc_ocy, file = "data/rrho_obj_pc_ocy.RData")
load("data/pc/rrho_obj_pc_ocy.RData")

rrho_obj_pc_opc <- RRHO(pc_opc, opc_pc, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_pc_opc <-  RRHO2_initialize(pc_opc, opc_pc, log10.ind=TRUE)
save(rrho_obj_pc_opc, file = "data/rrho_obj_pc_opc.RData")
load("data/pc/rrho_obj_pc_opc.RData")

rrho_obj_pc_adcy <- RRHO(pc_adcy, adcy_pc, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_pc_adcy <-  RRHO2_initialize(pc_adcy, adcy_pc, log10.ind=TRUE)
save(rrho_obj_pc_adcy, file = "data/rrho_obj_pc_adcy.RData")
load("data/pc/rrho_obj_pc_adcy.RData")

rrho_obj_pc_nkc <- RRHO(pc_nkc, nkc_pc, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_pc_nkc <-  RRHO2_initialize(pc_nkc, nkc_pc, log10.ind=TRUE)
save(rrho_obj_pc_nkc, file = "data/rrho_obj_pc_nkc.RData")
load("data/pc/rrho_obj_pc_nkc.RData")

rrho_obj_pc_neut <- RRHO(pc_neut, neut_pc, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_pc_neut <-  RRHO2_initialize(pc_neut, neut_pc, log10.ind=TRUE)
save(rrho_obj_pc_neut, file = "data/rrho_obj_pc_neut.RData")
load("data/pc/rrho_obj_pc_neut.RData")

rrho_obj_pc_mp <- RRHO(pc_mp, mp_pc, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_pc_mp <-  RRHO2_initialize(pc_mp, mp_pc, log10.ind=TRUE)
save(rrho_obj_pc_mp, file = "data/rrho_obj_pc_mp.RData")
load("data/pc/rrho_obj_pc_mp.RData")

rrho_obj_pc_treg <- RRHO(pc_treg, treg_pc, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_pc_treg <-  RRHO2_initialize(pc_treg, treg_pc, log10.ind=TRUE)
save(rrho_obj_pc_treg, file = "data/rrho_obj_pc_treg.RData")
load("data/pc/rrho_obj_pc_treg.RData")

rrho_obj_pc_hsc <- RRHO(pc_hsc, hsc_pc, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_pc_hsc <-  RRHO2_initialize(pc_hsc, hsc_pc, log10.ind=TRUE)
save(rrho_obj_pc_hsc, file = "data/rrho_obj_pc_hsc.RData")
load("data/pc/rrho_obj_pc_hsc.RData")

rrho_obj_pc_wbm <- RRHO(pc_wbm, wbm_pc, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_pc_wbm <-  RRHO2_initialize(pc_wbm, wbm_pc, log10.ind=TRUE)
save(rrho_obj_pc_wbm, file = "data/rrho_obj_pc_wbm.RData")
load("data/pc/rrho_obj_pc_wbm.RData")

#_________________________________________________________________________________________________

## plot plasma cells heatmaps

pc_bmsc_plot <- levelplot(rrho_obj_pc_bmsc$hypermat, col.regions = col_fun(min(rrho_obj_pc_bmsc$hypermat):max(rrho_obj_pc_bmsc$hypermat)), 
                          scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                          xlab=textGrob("Plasma Cells", vjust = 0.5, gp = gpar(cex = 2)), 
                          ylab= textGrob("Bone Marrow Stromal Cells", vjust = 0.5, gp = gpar(cex = 1.5), rot = 90), 
                          pretty=TRUE)
#RRHO2_heatmap(rrho2_obj_pc_bmsc, colorGradient=col)

pc_ec_plot <- levelplot(rrho_obj_pc_ec$hypermat, col.regions = col_fun(min(rrho_obj_pc_ec$hypermat):max(rrho_obj_pc_ec$hypermat)), 
          scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
          xlab="", 
          ylab= textGrob("Endothelial Cells", vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
          pretty=TRUE)  
#RRHO2_heatmap(rrho2_obj_pc_ec, colorGradient=col)

pc_ocy_plot <- levelplot(rrho_obj_pc_ocy$hypermat, col.regions = col_fun(min(rrho_obj_pc_ocy$hypermat):max(rrho_obj_pc_ocy$hypermat)), 
          scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
          xlab="", 
          ylab= textGrob("Osteocyte", vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
          pretty=TRUE)  
#RRHO2_heatmap(rrho2_obj_pc_ocy, colorGradient=col)

pc_opc_plot <- levelplot(rrho_obj_pc_opc$hypermat, col.regions = col_fun(min(rrho_obj_pc_opc$hypermat):max(rrho_obj_pc_opc$hypermat)), 
          scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
          xlab="", 
          ylab= textGrob("Osteogenic Precursor Cells", vjust = 0.5, gp = gpar(cex = 2), rot = 90),  
          pretty=TRUE)
#RRHO2_heatmap(rrho2_obj_pc_opc, colorGradient=col)

pc_adcy_plot <-levelplot(rrho_obj_pc_adcy$hypermat, col.regions = col_fun(min(rrho_obj_pc_adcy$hypermat):max(rrho_obj_pc_adcy$hypermat)), 
          scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
          xlab="", 
          ylab= textGrob("Adipocyte" , vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
          pretty=TRUE)  
#RRHO2_heatmap(rrho2_obj_pc_adcy, colorGradient=col)

pc_nkc_plot <-levelplot(rrho_obj_pc_nkc$hypermat, col.regions = col_fun(min(rrho_obj_pc_nkc$hypermat):max(rrho_obj_pc_nkc$hypermat)), 
          scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
          xlab="",
          ylab=  textGrob("Natural Killer Cells", vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
          pretty=TRUE)
#RRHO2_heatmap(rrho2_obj_pc_nkc, colorGradient=col)

pc_neut_plot <-levelplot(rrho_obj_pc_neut$hypermat, col.regions = col_fun(min(rrho_obj_pc_neut$hypermat):max(rrho_obj_pc_neut$hypermat)), 
          scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
          xlab="", 
          ylab= textGrob("Neutrophil", vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
          pretty=TRUE) 
#RRHO2_heatmap(rrho2_obj_pc_neut, colorGradient=col)

pc_mp_plot <-levelplot(rrho_obj_pc_mp$hypermat, col.regions = col_fun(min(rrho_obj_pc_mp$hypermat):max(rrho_obj_pc_mp$hypermat)), 
          scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
          xlab="", 
          ylab= textGrob("Macrophage", vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
          pretty=TRUE)
#RRHO2_heatmap(rrho2_obj_pc_mp, colorGradient=col)

pc_treg_plot <-levelplot(rrho_obj_pc_treg$hypermat, col.regions = col_fun(min(rrho_obj_pc_treg$hypermat):max(rrho_obj_pc_treg$hypermat)), 
          scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
          xlab="", 
          ylab=  textGrob("Regulatory T Cells", vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
          pretty=TRUE)  
#RRHO2_heatmap(rrho2_obj_pc_treg, colorGradient=col)

pc_hsc_plot <-levelplot(rrho_obj_pc_hsc$hypermat, col.regions = col_fun(min(rrho_obj_pc_hsc$hypermat):max(rrho_obj_pc_hsc$hypermat)), 
          scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
          xlab="", 
          ylab= textGrob("Hematopoietic Stem Cells" , vjust = 0.5, gp = gpar(cex = 1.5), rot = 90), 
          pretty=TRUE) 
#RRHO2_heatmap(rrho2_obj_pc_hsc, colorGradient=col)


pc_wbm_plot <-levelplot(rrho_obj_pc_wbm$hypermat, col.regions = col_fun(min(rrho_obj_pc_wbm$hypermat):max(rrho_obj_pc_wbm$hypermat)), 
          scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
          xlab="", 
          ylab=  textGrob("Whole Bone Marrow", vjust = 0.5, gp = gpar(cex = 1.5), rot = 90), 
          pretty=TRUE)  
#RRHO2_heatmap(rrho2_obj_pc_wbm, colorGradient=col)

#_________________________________________________________________________________________________
#_________________________________________________________________________________________________

## take bone marrow stromal cell intersect


bmsc_ec_intersect <- intersect(bmsc_data$gene_symbol,ec_data$gene_symbol)
bmsc_ec <- bmsc_data[bmsc_data$gene_symbol %in% bmsc_ec_intersect, ]
ec_bmsc <- ec_data[ec_data$gene_symbol %in% bmsc_ec_intersect, ]

bmsc_ocy_intersect <- intersect(bmsc_data$gene_symbol,ocy_data$gene_symbol)
bmsc_ocy <- bmsc_data[bmsc_data$gene_symbol %in% bmsc_ocy_intersect, ]
ocy_bmsc <- ocy_data[ocy_data$gene_symbol %in% bmsc_ocy_intersect, ]

bmsc_opc_intersect <- intersect(bmsc_data$gene_symbol,opc_data$gene_symbol)
bmsc_opc <- bmsc_data[bmsc_data$gene_symbol %in% bmsc_opc_intersect, ]
opc_bmsc <- opc_data[opc_data$gene_symbol %in% bmsc_opc_intersect, ]

bmsc_adcy_intersect <- intersect(bmsc_data$gene_symbol,adcy_data$gene_symbol)
bmsc_adcy <- bmsc_data[bmsc_data$gene_symbol %in% bmsc_adcy_intersect, ]
adcy_bmsc <- adcy_data[adcy_data$gene_symbol %in% bmsc_adcy_intersect, ]

bmsc_nkc_intersect <- intersect(bmsc_data$gene_symbol,nkc_data$gene_symbol)
bmsc_nkc <- bmsc_data[bmsc_data$gene_symbol %in% bmsc_nkc_intersect, ]
nkc_bmsc <- nkc_data[nkc_data$gene_symbol %in% bmsc_nkc_intersect, ]

bmsc_neut_intersect <- intersect(bmsc_data$gene_symbol,neut_data$gene_symbol)
bmsc_neut <- bmsc_data[bmsc_data$gene_symbol %in% bmsc_neut_intersect, ]
neut_bmsc <- neut_data[neut_data$gene_symbol %in% bmsc_neut_intersect, ]

bmsc_mp_intersect <- intersect(bmsc_data$gene_symbol,mp_data$gene_symbol)
bmsc_mp <- bmsc_data[bmsc_data$gene_symbol %in% bmsc_mp_intersect, ]
mp_bmsc <- mp_data[mp_data$gene_symbol %in% bmsc_mp_intersect, ]

bmsc_treg_intersect <- intersect(bmsc_data$gene_symbol,treg_data$gene_symbol)
bmsc_treg <- bmsc_data[bmsc_data$gene_symbol %in% bmsc_treg_intersect, ]
treg_bmsc <- treg_data[treg_data$gene_symbol %in% bmsc_treg_intersect, ]


bmsc_hsc_intersect <- intersect(bmsc_data$gene_symbol,hsc_data$gene_symbol)
bmsc_hsc <- bmsc_data[bmsc_data$gene_symbol %in% bmsc_hsc_intersect, ]
hsc_bmsc <- hsc_data[hsc_data$gene_symbol %in% bmsc_hsc_intersect, ]

bmsc_wbm_intersect <- intersect(bmsc_data$gene_symbol,wbm_data$gene_symbol)
bmsc_wbm <- bmsc_data[bmsc_data$gene_symbol %in% bmsc_wbm_intersect, ]
wbm_bmsc <- wbm_data[wbm_data$gene_symbol %in% bmsc_wbm_intersect, ]


#_________________________________________________________________________________________________

## compute bone marrow stromal cell overlap and significance

rrho_obj_bmsc_ec <- RRHO(bmsc_ec, ec_bmsc, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_bmsc_ec <-  RRHO2_initialize(bmsc_ec, ec_bmsc, log10.ind=TRUE)
save(rrho_obj_bmsc_ec, file = "data/bmsc/rrho_obj_bmsc_ec.RData")
load("data/bmsc/rrho_obj_bmsc_ec.RData")

rrho_obj_bmsc_ocy <- RRHO(bmsc_ocy, ocy_bmsc, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_bmsc_ocy <-  RRHO2_initialize(bmsc_ocy, ocy_bmsc, log10.ind=TRUE)
save(rrho_obj_bmsc_ocy, file = "data/bmsc/rrho_obj_bmsc_ocy.RData")
load("data/bmsc/rrho_obj_bmsc_ocy.RData")

rrho_obj_bmsc_opc <- RRHO(bmsc_opc, opc_bmsc, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_bmsc_opc <-  RRHO2_initialize(bmsc_opc, opc_bmsc, log10.ind=TRUE)
save(rrho_obj_bmsc_opc, file = "data/bmsc/rrho_obj_bmsc_opc.RData")
load("data/bmsc/rrho_obj_bmsc_opc.RData")

rrho_obj_bmsc_adcy <- RRHO(bmsc_adcy, adcy_bmsc, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_bmsc_adcy <-  RRHO2_initialize(bmsc_adcy, adcy_bmsc, log10.ind=TRUE)
save(rrho_obj_bmsc_adcy, file = "data/bmsc/rrho_obj_bmsc_adcy.RData")
load("data/bmsc/rrho_obj_bmsc_adcy.RData")

rrho_obj_bmsc_nkc <- RRHO(bmsc_nkc, nkc_bmsc, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_bmsc_nkc <-  RRHO2_initialize(bmsc_nkc, nkc_bmsc, log10.ind=TRUE)
save(rrho_obj_bmsc_nkc, file = "data/bmsc/rrho_obj_bmsc_nkc.RData")
load("data/bmsc/rrho_obj_bmsc_nkc.RData")

rrho_obj_bmsc_neut <- RRHO(bmsc_neut, neut_bmsc, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_bmsc_neut <-  RRHO2_initialize(bmsc_neut, neut_bmsc, log10.ind=TRUE)
save(rrho_obj_bmsc_neut, file = "data/bmsc/rrho_obj_bmsc_neut.RData")
load("data/bmsc/rrho_obj_bmsc_neut.RData")

rrho_obj_bmsc_mp <- RRHO(bmsc_mp, mp_bmsc, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_bmsc_mp <-  RRHO2_initialize(bmsc_mp, mp_bmsc, log10.ind=TRUE)
save(rrho_obj_bmsc_mp, file = "data/bmsc/rrho_obj_bmsc_mp.RData")
load("data/bmsc/rrho_obj_bmsc_mp.RData")

rrho_obj_bmsc_treg <- RRHO(bmsc_treg, treg_bmsc, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_bmsc_treg <-  RRHO2_initialize(bmsc_treg, treg_bmsc, log10.ind=TRUE)
save(rrho_obj_bmsc_treg, file = "data/bmsc/rrho_obj_bmsc_treg.RData")
load("data/bmsc/rrho_obj_bmsc_treg.RData")

rrho_obj_bmsc_hsc <- RRHO(bmsc_hsc, hsc_bmsc, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_bmsc_hsc <-  RRHO2_initialize(bmsc_hsc, hsc_bmsc, log10.ind=TRUE)
save(rrho_obj_bmsc_hsc, file = "data/bmsc/rrho_obj_bmsc_hsc.RData")
load("data/bmsc/rrho_obj_bmsc_hsc.RData")

rrho_obj_bmsc_wbm <- RRHO(bmsc_wbm, wbm_bmsc, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_bmsc_wbm <-  RRHO2_initialize(bmsc_wbm, wbm_bmsc, log10.ind=TRUE)
save(rrho_obj_bmsc_wbm, file = "data/bmsc/rrho_obj_bmsc_wbm.RData")
load("data/bmsc/rrho_obj_bmsc_wbm.RData")


#_________________________________________________________________________________________________

## plot bone marrow stromal cell heatmaps

bmsc_ec_plot <- levelplot(rrho_obj_bmsc_ec$hypermat, col.regions = col_fun(min(rrho_obj_bmsc_ec$hypermat):max(rrho_obj_bmsc_ec$hypermat)), 
                          scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                          xlab=textGrob("Bone Marrow Stromal Cells", vjust = 1, gp = gpar(cex = 1.5)), 
                          ylab= "", #textGrob("Endothelial Cells", vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
                          pretty=TRUE, par.settings=list(layout.heights=list(top.padding=-3))) 
#RRHO2_heatmap(rrho2_obj_bmsc_ec, colorGradient=col)

bmsc_ocy_plot <- levelplot(rrho_obj_bmsc_ocy$hypermat, col.regions = col_fun(min(rrho_obj_bmsc_ocy$hypermat):max(rrho_obj_bmsc_ocy$hypermat)), 
                           scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                           xlab="", 
                           ylab= "", #textGrob("Osteocyte", vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
                           pretty=TRUE) 
#RRHO2_heatmap(rrho2_obj_bmsc_ocy, colorGradient=col)

bmsc_opc_plot <- levelplot(rrho_obj_bmsc_opc$hypermat, col.regions = col_fun(min(rrho_obj_bmsc_opc$hypermat):max(rrho_obj_bmsc_opc$hypermat)), 
                           scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                           xlab="", 
                           ylab= "", #textGrob("Osteogenic Precursor Cells", vjust = 0.5, gp = gpar(cex = 2), rot = 90),  
                           pretty=TRUE) 
#RRHO2_heatmap(rrho2_obj_bmsc_opc, colorGradient=col)

bmsc_adcy_plot <-levelplot(rrho_obj_bmsc_adcy$hypermat, col.regions = col_fun(min(rrho_obj_bmsc_adcy$hypermat):max(rrho_obj_bmsc_adcy$hypermat)), 
                           scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                           xlab="", 
                           ylab= "", #textGrob("Adipocyte" , vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
                           pretty=TRUE) 
#RRHO2_heatmap(rrho2_obj_bmsc_adcy, colorGradient=col)

bmsc_nkc_plot <-levelplot(rrho_obj_bmsc_nkc$hypermat, col.regions = col_fun(min(rrho_obj_bmsc_nkc$hypermat):max(rrho_obj_bmsc_nkc$hypermat)), 
                          scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                          xlab="",
                          ylab= "", #textGrob("Natural Killer Cells", vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
                          pretty=TRUE) 
#RRHO2_heatmap(rrho2_obj_bmsc_nkc, colorGradient=col)

bmsc_neut_plot <-levelplot(rrho_obj_bmsc_neut$hypermat, col.regions = col_fun(min(rrho_obj_bmsc_neut$hypermat):max(rrho_obj_bmsc_neut$hypermat)), 
                           scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                           xlab="", 
                           ylab= "", #textGrob("Neutrophil", vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
                           pretty=TRUE) 
#RRHO2_heatmap(rrho2_obj_bmsc_neut, colorGradient=col)

bmsc_mp_plot <-levelplot(rrho_obj_bmsc_mp$hypermat, col.regions = col_fun(min(rrho_obj_bmsc_mp$hypermat):max(rrho_obj_bmsc_mp$hypermat)), 
                         scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                         xlab="", 
                         ylab= "", #textGrob("Macrophage", vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
                         pretty=TRUE) 
#RRHO2_heatmap(rrho2_obj_bmsc_mp, colorGradient=col)

bmsc_treg_plot <-levelplot(rrho_obj_bmsc_treg$hypermat, col.regions = col_fun(min(rrho_obj_bmsc_treg$hypermat):max(rrho_obj_bmsc_treg$hypermat)), 
                           scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                           xlab="", 
                           ylab= "", #textGrob("Regulatory T Cells", vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
                           pretty=TRUE) 
#RRHO2_heatmap(rrho2_obj_bmsc_treg, colorGradient=col)

bmsc_hsc_plot <-levelplot(rrho_obj_bmsc_hsc$hypermat, col.regions = col_fun(min(rrho_obj_bmsc_hsc$hypermat):max(rrho_obj_bmsc_hsc$hypermat)), 
                          scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                          xlab="", 
                          ylab= "", #textGrob("Hematopoietic Stem Cells" , vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
                          pretty=TRUE) 
#RRHO2_heatmap(rrho2_obj_bmsc_hsc, colorGradient=col)


bmsc_wbm_plot <-levelplot(rrho_obj_bmsc_wbm$hypermat, col.regions = col_fun(min(rrho_obj_bmsc_wbm$hypermat):max(rrho_obj_bmsc_wbm$hypermat)), 
                          scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                          xlab="", 
                          ylab= "", #textGrob("Whole Bone Marrow", vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
                          pretty=TRUE) 
#RRHO2_heatmap(rrho2_obj_bmsc_wbm, colorGradient=col)

#_________________________________________________________________________________________________
#_________________________________________________________________________________________________

## take endothelial cell intersect

ec_ocy_intersect <- intersect(ec_data$gene_symbol,ocy_data$gene_symbol)
ec_ocy <- ec_data[ec_data$gene_symbol %in% ec_ocy_intersect, ]
ocy_ec <- ocy_data[ocy_data$gene_symbol %in% ec_ocy_intersect, ]

ec_opc_intersect <- intersect(ec_data$gene_symbol,opc_data$gene_symbol)
ec_opc <- ec_data[ec_data$gene_symbol %in% ec_opc_intersect, ]
opc_ec <- opc_data[opc_data$gene_symbol %in% ec_opc_intersect, ]

ec_adcy_intersect <- intersect(ec_data$gene_symbol,adcy_data$gene_symbol)
ec_adcy <- ec_data[ec_data$gene_symbol %in% ec_adcy_intersect, ]
adcy_ec <- adcy_data[adcy_data$gene_symbol %in% ec_adcy_intersect, ]

ec_nkc_intersect <- intersect(ec_data$gene_symbol,nkc_data$gene_symbol)
ec_nkc <- ec_data[ec_data$gene_symbol %in% ec_nkc_intersect, ]
nkc_ec <- nkc_data[nkc_data$gene_symbol %in% ec_nkc_intersect, ]

ec_neut_intersect <- intersect(ec_data$gene_symbol,neut_data$gene_symbol)
ec_neut <- ec_data[ec_data$gene_symbol %in% ec_neut_intersect, ]
neut_ec <- neut_data[neut_data$gene_symbol %in% ec_neut_intersect, ]

ec_mp_intersect <- intersect(ec_data$gene_symbol,mp_data$gene_symbol)
ec_mp <- ec_data[ec_data$gene_symbol %in% ec_mp_intersect, ]
mp_ec <- mp_data[mp_data$gene_symbol %in% ec_mp_intersect, ]

ec_treg_intersect <- intersect(ec_data$gene_symbol,treg_data$gene_symbol)
ec_treg <- ec_data[ec_data$gene_symbol %in% ec_treg_intersect, ]
treg_ec <- treg_data[treg_data$gene_symbol %in% ec_treg_intersect, ]


ec_hsc_intersect <- intersect(ec_data$gene_symbol,hsc_data$gene_symbol)
ec_hsc <- ec_data[ec_data$gene_symbol %in% ec_hsc_intersect, ]
hsc_ec <- hsc_data[hsc_data$gene_symbol %in% ec_hsc_intersect, ]

ec_wbm_intersect <- intersect(ec_data$gene_symbol,wbm_data$gene_symbol)
ec_wbm <- ec_data[ec_data$gene_symbol %in% ec_wbm_intersect, ]
wbm_ec <- wbm_data[wbm_data$gene_symbol %in% ec_wbm_intersect, ]


#_________________________________________________________________________________________________

## compute endothelial cell overlap and significance

rrho_obj_ec_ocy <- RRHO(ec_ocy, ocy_ec, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_ec_ocy <-  RRHO2_initialize(ec_ocy, ocy_ec, log10.ind=TRUE)
save(rrho_obj_ec_ocy, file = "data/ec/rrho_obj_ec_ocy.RData")
load("data/ec/rrho_obj_ec_ocy.RData")

rrho_obj_ec_opc <- RRHO(ec_opc, opc_ec, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_ec_opc <-  RRHO2_initialize(ec_opc, opc_ec, log10.ind=TRUE)
save(rrho_obj_ec_opc, file = "data/ec/rrho_obj_ec_opc.RData")
load("data/ec/rrho_obj_ec_opc.RData")

rrho_obj_ec_adcy <- RRHO(ec_adcy, adcy_ec, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_ec_adcy <-  RRHO2_initialize(ec_adcy, adcy_ec, log10.ind=TRUE)
save(rrho_obj_ec_adcy, file = "data/ec/rrho_obj_ec_adcy.RData")
load("data/ec/rrho_obj_ec_adcy.RData")

rrho_obj_ec_nkc <- RRHO(ec_nkc, nkc_ec, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_ec_nkc <-  RRHO2_initialize(ec_nkc, nkc_ec, log10.ind=TRUE)
save(rrho_obj_ec_nkc, file = "data/ec/rrho_obj_ec_nkc.RData")
load("data/ec/rrho_obj_ec_nkc.RData")

rrho_obj_ec_neut <- RRHO(ec_neut, neut_ec, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_ec_neut <-  RRHO2_initialize(ec_neut, neut_ec, log10.ind=TRUE)
save(rrho_obj_ec_neut, file = "data/ec/rrho_obj_ec_neut.RData")
load("data/ec/rrho_obj_ec_neut.RData")

rrho_obj_ec_mp <- RRHO(ec_mp, mp_ec, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_ec_mp <-  RRHO2_initialize(ec_mp, mp_ec, log10.ind=TRUE)
save(rrho_obj_ec_mp, file = "data/ec/rrho_obj_ec_mp.RData")
load("data/ec/rrho_obj_ec_mp.RData")

rrho_obj_ec_treg <- RRHO(ec_treg, treg_ec, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_ec_treg <-  RRHO2_initialize(ec_treg, treg_ec, log10.ind=TRUE)
save(rrho_obj_ec_treg, file = "data/ec/rrho_obj_ec_treg.RData")
load("data/ec/rrho_obj_ec_treg.RData")

rrho_obj_ec_hsc <- RRHO(ec_hsc, hsc_ec, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_ec_hsc <-  RRHO2_initialize(ec_hsc, hsc_ec, log10.ind=TRUE)
save(rrho_obj_ec_hsc, file = "data/ec/rrho_obj_ec_hsc.RData")
load("data/ec/rrho_obj_ec_hsc.RData")

rrho_obj_ec_wbm <- RRHO(ec_wbm, wbm_ec, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_ec_wbm <-  RRHO2_initialize(ec_wbm, wbm_ec, log10.ind=TRUE)
save(rrho_obj_ec_wbm, file = "data/ec/rrho_obj_ec_wbm.RData")
load("data/ec/rrho_obj_ec_wbm.RData")


#_________________________________________________________________________________________________

## plot endothelial cell heatmaps

ec_ocy_plot <- levelplot(rrho_obj_ec_ocy$hypermat, col.regions = col_fun(min(rrho_obj_ec_ocy$hypermat):max(rrho_obj_ec_ocy$hypermat)), 
                         scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                         xlab= textGrob("Endothelial Cells", vjust = 1, gp = gpar(cex = 2)),  
                         ylab= "", #textGrob("Osteocyte", vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
                         pretty=TRUE) 
#RRHO2_heatmap(rrho2_obj_ec_ocy, colorGradient=col)

ec_opc_plot <- levelplot(rrho_obj_ec_opc$hypermat, col.regions = col_fun(min(rrho_obj_ec_opc$hypermat):max(rrho_obj_ec_opc$hypermat)), 
                         scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                         xlab="", 
                         ylab= "", #textGrob("Osteogenic Precursor Cells", vjust = 0.5, gp = gpar(cex = 2), rot = 90),  
                         pretty=TRUE) 
#RRHO2_heatmap(rrho2_obj_ec_opc, colorGradient=col)

ec_adcy_plot <-levelplot(rrho_obj_ec_adcy$hypermat, col.regions = col_fun(min(rrho_obj_ec_adcy$hypermat):max(rrho_obj_ec_adcy$hypermat)), 
                         scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                         xlab="", 
                         ylab= "", #textGrob("Adipocyte" , vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
                         pretty=TRUE) 
#RRHO2_heatmap(rrho2_obj_ec_adcy, colorGradient=col)

ec_nkc_plot <-levelplot(rrho_obj_ec_nkc$hypermat, col.regions = col_fun(min(rrho_obj_ec_nkc$hypermat):max(rrho_obj_ec_nkc$hypermat)), 
                        scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                        xlab="",
                        ylab= "", #textGrob("Natural Killer Cells", vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
                        pretty=TRUE) 
#RRHO2_heatmap(rrho2_obj_ec_nkc, colorGradient=col)

ec_neut_plot <-levelplot(rrho_obj_ec_neut$hypermat, col.regions = col_fun(min(rrho_obj_ec_neut$hypermat):max(rrho_obj_ec_neut$hypermat)), 
                         scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                         xlab="", 
                         ylab= "", #textGrob("Neutrophil", vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
                         pretty=TRUE) 
#RRHO2_heatmap(rrho2_obj_ec_neut, colorGradient=col)

ec_mp_plot <-levelplot(rrho_obj_ec_mp$hypermat, col.regions = col_fun(min(rrho_obj_ec_mp$hypermat):max(rrho_obj_ec_mp$hypermat)), 
                       scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                       xlab="", 
                       ylab= "", #textGrob("Macrophage", vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
                       pretty=TRUE) 
#RRHO2_heatmap(rrho2_obj_ec_mp, colorGradient=col)

ec_treg_plot <-levelplot(rrho_obj_ec_treg$hypermat, col.regions = col_fun(min(rrho_obj_ec_treg$hypermat):max(rrho_obj_ec_treg$hypermat)), 
                         scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                         xlab="", 
                         ylab= "", #textGrob("Regulatory T Cells", vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
                         pretty=TRUE) 
#RRHO2_heatmap(rrho2_obj_ec_treg, colorGradient=col)

ec_hsc_plot <-levelplot(rrho_obj_ec_hsc$hypermat, col.regions = col_fun(min(rrho_obj_ec_hsc$hypermat):max(rrho_obj_ec_hsc$hypermat)), 
                        scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                        xlab="", 
                        ylab= "", #textGrob("Hematopoietic Stem Cells" , vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
                        pretty=TRUE) 
#RRHO2_heatmap(rrho2_obj_ec_hsc, colorGradient=col)


ec_wbm_plot <-levelplot(rrho_obj_ec_wbm$hypermat, col.regions = col_fun(min(rrho_obj_ec_wbm$hypermat):max(rrho_obj_ec_wbm$hypermat)), 
                        scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                        xlab="", 
                        ylab= "", #textGrob("Whole Bone Marrow", vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
                        pretty=TRUE) 
#RRHO2_heatmap(rrho2_obj_ec_wbm, colorGradient=col)
#_________________________________________________________________________________________________
#_________________________________________________________________________________________________

## take osteocyte intersect

ocy_opc_intersect <- intersect(ocy_data$gene_symbol,opc_data$gene_symbol)
ocy_opc <- ocy_data[ocy_data$gene_symbol %in% ocy_opc_intersect, ]
opc_ocy <- opc_data[opc_data$gene_symbol %in% ocy_opc_intersect, ]

ocy_adcy_intersect <- intersect(ocy_data$gene_symbol,adcy_data$gene_symbol)
ocy_adcy <- ocy_data[ocy_data$gene_symbol %in% ocy_adcy_intersect, ]
adcy_ocy <- adcy_data[adcy_data$gene_symbol %in% ocy_adcy_intersect, ]

ocy_nkc_intersect <- intersect(ocy_data$gene_symbol,nkc_data$gene_symbol)
ocy_nkc <- ocy_data[ocy_data$gene_symbol %in% ocy_nkc_intersect, ]
nkc_ocy <- nkc_data[nkc_data$gene_symbol %in% ocy_nkc_intersect, ]

ocy_neut_intersect <- intersect(ocy_data$gene_symbol,neut_data$gene_symbol)
ocy_neut <- ocy_data[ocy_data$gene_symbol %in% ocy_neut_intersect, ]
neut_ocy <- neut_data[neut_data$gene_symbol %in% ocy_neut_intersect, ]

ocy_mp_intersect <- intersect(ocy_data$gene_symbol,mp_data$gene_symbol)
ocy_mp <- ocy_data[ocy_data$gene_symbol %in% ocy_mp_intersect, ]
mp_ocy <- mp_data[mp_data$gene_symbol %in% ocy_mp_intersect, ]

ocy_treg_intersect <- intersect(ocy_data$gene_symbol,treg_data$gene_symbol)
ocy_treg <- ocy_data[ocy_data$gene_symbol %in% ocy_treg_intersect, ]
treg_ocy <- treg_data[treg_data$gene_symbol %in% ocy_treg_intersect, ]


ocy_hsc_intersect <- intersect(ocy_data$gene_symbol,hsc_data$gene_symbol)
ocy_hsc <- ocy_data[ocy_data$gene_symbol %in% ocy_hsc_intersect, ]
hsc_ocy <- hsc_data[hsc_data$gene_symbol %in% ocy_hsc_intersect, ]

ocy_wbm_intersect <- intersect(ocy_data$gene_symbol,wbm_data$gene_symbol)
ocy_wbm <- ocy_data[ocy_data$gene_symbol %in% ocy_wbm_intersect, ]
wbm_ocy <- wbm_data[wbm_data$gene_symbol %in% ocy_wbm_intersect, ]


#_________________________________________________________________________________________________

## compute osteocyte overlap and significance


rrho_obj_ocy_opc <- RRHO(ocy_opc, opc_ocy, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_ocy_opc <-  RRHO2_initialize(ocy_opc, opc_ocy, log10.ind=TRUE)
save(rrho_obj_ocy_opc, file = "data/ocy/rrho_obj_ocy_opc.RData")
load("data/ocy/rrho_obj_ocy_opc.RData")

rrho_obj_ocy_adcy <- RRHO(ocy_adcy, adcy_ocy, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_ocy_adcy <-  RRHO2_initialize(ocy_adcy, adcy_ocy, log10.ind=TRUE)
save(rrho_obj_ocy_adcy, file = "data/ocy/rrho_obj_ocy_adcy.RData")
load("data/ocy/rrho_obj_ocy_adcy.RData")

rrho_obj_ocy_nkc <- RRHO(ocy_nkc, nkc_ocy, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_ocy_nkc <-  RRHO2_initialize(ocy_nkc, nkc_ocy, log10.ind=TRUE)
save(rrho_obj_ocy_nkc, file = "data/ocy/rrho_obj_ocy_nkc.RData")
load("data/ocy/rrho_obj_ocy_nkc.RData")

rrho_obj_ocy_neut <- RRHO(ocy_neut, neut_ocy, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_ocy_neut <-  RRHO2_initialize(ocy_neut, neut_ocy, log10.ind=TRUE)
save(rrho_obj_ocy_neut, file = "data/ocy/rrho_obj_ocy_neut.RData")
load("data/ocy/rrho_obj_ocy_neut.RData")

rrho_obj_ocy_mp <- RRHO(ocy_mp, mp_ocy, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_ocy_mp <-  RRHO2_initialize(ocy_mp, mp_ocy, log10.ind=TRUE)
save(rrho_obj_ocy_mp, file = "data/ocy/rrho_obj_ocy_mp.RData")
load("data/ocy/rrho_obj_ocy_mp.RData")

rrho_obj_ocy_treg <- RRHO(ocy_treg, treg_ocy, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_ocy_treg <-  RRHO2_initialize(ocy_treg, treg_ocy, log10.ind=TRUE)
save(rrho_obj_ocy_treg, file = "data/ocy/rrho_obj_ocy_treg.RData")
load("data/ocy/rrho_obj_ocy_treg.RData")

rrho_obj_ocy_hsc <- RRHO(ocy_hsc, hsc_ocy, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_ocy_hsc <-  RRHO2_initialize(ocy_hsc, hsc_ocy, log10.ind=TRUE)
save(rrho_obj_ocy_hsc, file = "data/ocy/rrho_obj_ocy_hsc.RData")
load("data/ocy/rrho_obj_ocy_hsc.RData")

rrho_obj_ocy_wbm <- RRHO(ocy_wbm, wbm_ocy, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_ocy_wbm <-  RRHO2_initialize(ocy_wbm, wbm_ocy, log10.ind=TRUE)
save(rrho_obj_ocy_wbm, file = "data/ocy/rrho_obj_ocy_wbm.RData")
load("data/ocy/rrho_obj_ocy_wbm.RData")


#_________________________________________________________________________________________________

## plot osteocyte heatmaps


ocy_opc_plot <- levelplot(rrho_obj_ocy_opc$hypermat, col.regions = col_fun(min(rrho_obj_ocy_opc$hypermat):max(rrho_obj_ocy_opc$hypermat)), 
                          scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                          xlab= textGrob("Osteocyte", vjust = 1, gp = gpar(cex = 2)),  
                          ylab= "", #textGrob("Osteogenic Procyursor Cells", vjust = 0.5, gp = gpar(cex = 2), rot = 90),  
                          pretty=TRUE) 
#RRHO2_heatmap(rrho2_obj_ocy_opc, colorGradient=col)

ocy_adcy_plot <-levelplot(rrho_obj_ocy_adcy$hypermat, col.regions = col_fun(min(rrho_obj_ocy_adcy$hypermat):max(rrho_obj_ocy_adcy$hypermat)), 
                          scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                          xlab="", 
                          ylab= "", #textGrob("Adipocyte" , vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
                          pretty=TRUE) 
#RRHO2_heatmap(rrho2_obj_ocy_adcy, colorGradient=col)

ocy_nkc_plot <-levelplot(rrho_obj_ocy_nkc$hypermat, col.regions = col_fun(min(rrho_obj_ocy_nkc$hypermat):max(rrho_obj_ocy_nkc$hypermat)), 
                         scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                         xlab="",
                         ylab= "", #textGrob("Natural Killer Cells", vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
                         pretty=TRUE) 
#RRHO2_heatmap(rrho2_obj_ocy_nkc, colorGradient=col)

ocy_neut_plot <-levelplot(rrho_obj_ocy_neut$hypermat, col.regions = col_fun(min(rrho_obj_ocy_neut$hypermat):max(rrho_obj_ocy_neut$hypermat)), 
                          scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                          xlab="", 
                          ylab= "", #textGrob("Neutrophil", vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
                          pretty=TRUE) 
#RRHO2_heatmap(rrho2_obj_ocy_neut, colorGradient=col)

ocy_mp_plot <-levelplot(rrho_obj_ocy_mp$hypermat, col.regions = col_fun(min(rrho_obj_ocy_mp$hypermat):max(rrho_obj_ocy_mp$hypermat)), 
                        scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                        xlab="", 
                        ylab= "", #textGrob("Macrophage", vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
                        pretty=TRUE) 
#RRHO2_heatmap(rrho2_obj_ocy_mp, colorGradient=col)

ocy_treg_plot <-levelplot(rrho_obj_ocy_treg$hypermat, col.regions = col_fun(min(rrho_obj_ocy_treg$hypermat):max(rrho_obj_ocy_treg$hypermat)), 
                          scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                          xlab="", 
                          ylab= "", #textGrob("Regulatory T Cells", vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
                          pretty=TRUE) 
#RRHO2_heatmap(rrho2_obj_ocy_treg, colorGradient=col)

ocy_hsc_plot <-levelplot(rrho_obj_ocy_hsc$hypermat, col.regions = col_fun(min(rrho_obj_ocy_hsc$hypermat):max(rrho_obj_ocy_hsc$hypermat)), 
                         scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                         xlab="", 
                         ylab= "", #textGrob("Hematopoietic Stem Cells" , vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
                         pretty=TRUE) 
#RRHO2_heatmap(rrho2_obj_ocy_hsc, colorGradient=col)


ocy_wbm_plot <-levelplot(rrho_obj_ocy_wbm$hypermat, col.regions = col_fun(min(rrho_obj_ocy_wbm$hypermat):max(rrho_obj_ocy_wbm$hypermat)), 
                         scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                         xlab="", 
                         ylab= "", #textGrob("Whole Bone Marrow", vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
                         pretty=TRUE) 
#RRHO2_heatmap(rrho2_obj_ocy_wbm, colorGradient=col)


#_________________________________________________________________________________________________
#_________________________________________________________________________________________________

## take osteogenic precursor cell intersect

opc_adcy_intersect <- intersect(opc_data$gene_symbol,adcy_data$gene_symbol)
opc_adcy <- opc_data[opc_data$gene_symbol %in% opc_adcy_intersect, ]
adcy_opc <- adcy_data[adcy_data$gene_symbol %in% opc_adcy_intersect, ]

opc_nkc_intersect <- intersect(opc_data$gene_symbol,nkc_data$gene_symbol)
opc_nkc <- opc_data[opc_data$gene_symbol %in% opc_nkc_intersect, ]
nkc_opc <- nkc_data[nkc_data$gene_symbol %in% opc_nkc_intersect, ]

opc_neut_intersect <- intersect(opc_data$gene_symbol,neut_data$gene_symbol)
opc_neut <- opc_data[opc_data$gene_symbol %in% opc_neut_intersect, ]
neut_opc <- neut_data[neut_data$gene_symbol %in% opc_neut_intersect, ]

opc_mp_intersect <- intersect(opc_data$gene_symbol,mp_data$gene_symbol)
opc_mp <- opc_data[opc_data$gene_symbol %in% opc_mp_intersect, ]
mp_opc <- mp_data[mp_data$gene_symbol %in% opc_mp_intersect, ]

opc_treg_intersect <- intersect(opc_data$gene_symbol,treg_data$gene_symbol)
opc_treg <- opc_data[opc_data$gene_symbol %in% opc_treg_intersect, ]
treg_opc <- treg_data[treg_data$gene_symbol %in% opc_treg_intersect, ]


opc_hsc_intersect <- intersect(opc_data$gene_symbol,hsc_data$gene_symbol)
opc_hsc <- opc_data[opc_data$gene_symbol %in% opc_hsc_intersect, ]
hsc_opc <- hsc_data[hsc_data$gene_symbol %in% opc_hsc_intersect, ]

opc_wbm_intersect <- intersect(opc_data$gene_symbol,wbm_data$gene_symbol)
opc_wbm <- opc_data[opc_data$gene_symbol %in% opc_wbm_intersect, ]
wbm_opc <- wbm_data[wbm_data$gene_symbol %in% opc_wbm_intersect, ]


#_________________________________________________________________________________________________

## compute osteogenic precursor cell overlap and significance


rrho_obj_opc_adcy <- RRHO(opc_adcy, adcy_opc, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_opc_adcy <-  RRHO2_initialize(opc_adcy, adcy_opc, log10.ind=TRUE)
save(rrho_obj_opc_adcy, file = "data/opc/rrho_obj_opc_adcy.RData")
load("data/opc/rrho_obj_opc_adcy.RData")

rrho_obj_opc_nkc <- RRHO(opc_nkc, nkc_opc, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_opc_nkc <-  RRHO2_initialize(opc_nkc, nkc_opc, log10.ind=TRUE)
save(rrho_obj_opc_nkc, file = "data/opc/rrho_obj_opc_nkc.RData")
load("data/opc/rrho_obj_opc_nkc.RData")

rrho_obj_opc_neut <- RRHO(opc_neut, neut_opc, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_opc_neut <-  RRHO2_initialize(opc_neut, neut_opc, log10.ind=TRUE)
save(rrho_obj_opc_neut, file = "data/opc/rrho_obj_opc_neut.RData")
load("data/opc/rrho_obj_opc_neut.RData")

rrho_obj_opc_mp <- RRHO(opc_mp, mp_opc, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_opc_mp <-  RRHO2_initialize(opc_mp, mp_opc, log10.ind=TRUE)
save(rrho_obj_opc_mp, file = "data/opc/rrho_obj_opc_mp.RData")
load("data/opc/rrho_obj_opc_mp.RData")

rrho_obj_opc_treg <- RRHO(opc_treg, treg_opc, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_opc_treg <-  RRHO2_initialize(opc_treg, treg_opc, log10.ind=TRUE)
save(rrho_obj_opc_treg, file = "data/opc/rrho_obj_opc_treg.RData")
load("data/opc/rrho_obj_opc_treg.RData")

rrho_obj_opc_hsc <- RRHO(opc_hsc, hsc_opc, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_opc_hsc <-  RRHO2_initialize(opc_hsc, hsc_opc, log10.ind=TRUE)
save(rrho_obj_opc_hsc, file = "data/opc/rrho_obj_opc_hsc.RData")
load("data/opc/rrho_obj_opc_hsc.RData")

rrho_obj_opc_wbm <- RRHO(opc_wbm, wbm_opc, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_opc_wbm <-  RRHO2_initialize(opc_wbm, wbm_opc, log10.ind=TRUE)
save(rrho_obj_opc_wbm, file = "data/opc/rrho_obj_opc_wbm.RData")
load("data/opc/rrho_obj_opc_wbm.RData")


#_________________________________________________________________________________________________

## plot osteogenic precursor cell heatmaps


opc_adcy_plot <-levelplot(rrho_obj_opc_adcy$hypermat, col.regions = col_fun(min(rrho_obj_opc_adcy$hypermat):max(rrho_obj_opc_adcy$hypermat)), 
                          scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                          xlab= textGrob("Osteogenic Propcursor Cells", vjust = 1, gp = gpar(cex = 2)),  
                          ylab= "", # textGrob("Adipopcte" , vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
                          pretty=TRUE) 
#RRHO2_heatmap(rrho2_obj_opc_adcy, colorGradient=col)

opc_nkc_plot <-levelplot(rrho_obj_opc_nkc$hypermat, col.regions = col_fun(min(rrho_obj_opc_nkc$hypermat):max(rrho_obj_opc_nkc$hypermat)), 
                         scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                         xlab="",
                         ylab= "", #textGrob("Natural Killer Cells", vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
                         pretty=TRUE) 
#RRHO2_heatmap(rrho2_obj_opc_nkc, colorGradient=col)

opc_neut_plot <-levelplot(rrho_obj_opc_neut$hypermat, col.regions = col_fun(min(rrho_obj_opc_neut$hypermat):max(rrho_obj_opc_neut$hypermat)), 
                          scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                          xlab="", 
                          ylab= "", #textGrob("Neutrophil", vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
                          pretty=TRUE) 
#RRHO2_heatmap(rrho2_obj_opc_neut, colorGradient=col)

opc_mp_plot <-levelplot(rrho_obj_opc_mp$hypermat, col.regions = col_fun(min(rrho_obj_opc_mp$hypermat):max(rrho_obj_opc_mp$hypermat)), 
                        scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                        xlab="", 
                        ylab= "", # textGrob("Macrophage", vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
                        pretty=TRUE) 
#RRHO2_heatmap(rrho2_obj_opc_mp, colorGradient=col)

opc_treg_plot <-levelplot(rrho_obj_opc_treg$hypermat, col.regions = col_fun(min(rrho_obj_opc_treg$hypermat):max(rrho_obj_opc_treg$hypermat)), 
                          scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                          xlab="", 
                          ylab= "", #textGrob("Regulatory T Cells", vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
                          pretty=TRUE) 
#RRHO2_heatmap(rrho2_obj_opc_treg, colorGradient=col)

opc_hsc_plot <-levelplot(rrho_obj_opc_hsc$hypermat, col.regions = col_fun(min(rrho_obj_opc_hsc$hypermat):max(rrho_obj_opc_hsc$hypermat)), 
                         scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                         xlab="", 
                         ylab= "", #textGrob("Hematopoietic Stem Cells" , vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
                         pretty=TRUE) 
#RRHO2_heatmap(rrho2_obj_opc_hsc, colorGradient=col)


opc_wbm_plot <-levelplot(rrho_obj_opc_wbm$hypermat, col.regions = col_fun(min(rrho_obj_opc_wbm$hypermat):max(rrho_obj_opc_wbm$hypermat)), 
                         scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                         xlab="", 
                         ylab= "", #textGrob("Whole Bone Marrow", vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
                         pretty=TRUE) 
#RRHO2_heatmap(rrho2_obj_opc_wbm, colorGradient=col)

#_________________________________________________________________________________________________
#_________________________________________________________________________________________________

## take adipocyte intersect


adcy_nkc_intersect <- intersect(adcy_data$gene_symbol,nkc_data$gene_symbol)
adcy_nkc <- adcy_data[adcy_data$gene_symbol %in% adcy_nkc_intersect, ]
nkc_adcy <- nkc_data[nkc_data$gene_symbol %in% adcy_nkc_intersect, ]

adcy_neut_intersect <- intersect(adcy_data$gene_symbol,neut_data$gene_symbol)
adcy_neut <- adcy_data[adcy_data$gene_symbol %in% adcy_neut_intersect, ]
neut_adcy <- neut_data[neut_data$gene_symbol %in% adcy_neut_intersect, ]

adcy_mp_intersect <- intersect(adcy_data$gene_symbol,mp_data$gene_symbol)
adcy_mp <- adcy_data[adcy_data$gene_symbol %in% adcy_mp_intersect, ]
mp_adcy <- mp_data[mp_data$gene_symbol %in% adcy_mp_intersect, ]

adcy_treg_intersect <- intersect(adcy_data$gene_symbol,treg_data$gene_symbol)
adcy_treg <- adcy_data[adcy_data$gene_symbol %in% adcy_treg_intersect, ]
treg_adcy <- treg_data[treg_data$gene_symbol %in% adcy_treg_intersect, ]


adcy_hsc_intersect <- intersect(adcy_data$gene_symbol,hsc_data$gene_symbol)
adcy_hsc <- adcy_data[adcy_data$gene_symbol %in% adcy_hsc_intersect, ]
hsc_adcy <- hsc_data[hsc_data$gene_symbol %in% adcy_hsc_intersect, ]

adcy_wbm_intersect <- intersect(adcy_data$gene_symbol,wbm_data$gene_symbol)
adcy_wbm <- adcy_data[adcy_data$gene_symbol %in% adcy_wbm_intersect, ]
wbm_adcy <- wbm_data[wbm_data$gene_symbol %in% adcy_wbm_intersect, ]


#_________________________________________________________________________________________________

## compute adipocyte overlap and significance


rrho_obj_adcy_nkc <- RRHO(adcy_nkc, nkc_adcy, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_adcy_nkc <-  RRHO2_initialize(adcy_nkc, nkc_adcy, log10.ind=TRUE)
save(rrho_obj_adcy_nkc, file = "data/adcy/rrho_obj_adcy_nkc.RData")
load("data/adcy/rrho_obj_adcy_nkc.RData")

rrho_obj_adcy_neut <- RRHO(adcy_neut, neut_adcy, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_adcy_neut <-  RRHO2_initialize(adcy_neut, neut_adcy, log10.ind=TRUE)
save(rrho_obj_adcy_neut, file = "data/adcy/rrho_obj_adcy_neut.RData")
load("data/adcy/rrho_obj_adcy_neut.RData")

rrho_obj_adcy_mp <- RRHO(adcy_mp, mp_adcy, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_adcy_mp <-  RRHO2_initialize(adcy_mp, mp_adcy, log10.ind=TRUE)
save(rrho_obj_adcy_mp, file = "data/adcy/rrho_obj_adcy_mp.RData")
load("data/adcy/rrho_obj_adcy_mp.RData")

rrho_obj_adcy_treg <- RRHO(adcy_treg, treg_adcy, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_adcy_treg <-  RRHO2_initialize(adcy_treg, treg_adcy, log10.ind=TRUE)
save(rrho_obj_adcy_treg, file = "data/adcy/rrho_obj_adcy_treg.RData")
load("data/adcy/rrho_obj_adcy_treg.RData")

rrho_obj_adcy_hsc <- RRHO(adcy_hsc, hsc_adcy, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_adcy_hsc <-  RRHO2_initialize(adcy_hsc, hsc_adcy, log10.ind=TRUE)
save(rrho_obj_adcy_hsc, file = "data/adcy/rrho_obj_adcy_hsc.RData")
load("data/adcy/rrho_obj_adcy_hsc.RData")

rrho_obj_adcy_wbm <- RRHO(adcy_wbm, wbm_adcy, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_adcy_wbm <-  RRHO2_initialize(adcy_wbm, wbm_adcy, log10.ind=TRUE)
save(rrho_obj_adcy_wbm, file = "data/adcy/rrho_obj_adcy_wbm.RData")
load("data/adcy/rrho_obj_adcy_wbm.RData")


#_________________________________________________________________________________________________

## plot adipocyte heatmaps

adcy_nkc_plot <-levelplot(rrho_obj_adcy_nkc$hypermat, col.regions = col_fun(min(rrho_obj_adcy_nkc$hypermat):max(rrho_obj_adcy_nkc$hypermat)), 
                          scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                          xlab= textGrob("Adipocyte", vjust = 1, gp = gpar(cex = 2)),  
                          ylab= "", #  textGrob("Natural Killer Cells", vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
                          pretty=TRUE) 
#RRHO2_heatmap(rrho2_obj_adcy_nkc, colorGradient=col)

adcy_neut_plot <-levelplot(rrho_obj_adcy_neut$hypermat, col.regions = col_fun(min(rrho_obj_adcy_neut$hypermat):max(rrho_obj_adcy_neut$hypermat)), 
                           scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                           xlab="", 
                           ylab= "", # textGrob("Neutrophil", vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
                           pretty=TRUE) 
#RRHO2_heatmap(rrho2_obj_adcy_neut, colorGradient=col)

adcy_mp_plot <-levelplot(rrho_obj_adcy_mp$hypermat, col.regions = col_fun(min(rrho_obj_adcy_mp$hypermat):max(rrho_obj_adcy_mp$hypermat)), 
                         scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                         xlab="", 
                         ylab= "", # textGrob("Macrophage", vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
                         pretty=TRUE) 
#RRHO2_heatmap(rrho2_obj_adcy_mp, colorGradient=col)

adcy_treg_plot <-levelplot(rrho_obj_adcy_treg$hypermat, col.regions = col_fun(min(rrho_obj_adcy_treg$hypermat):max(rrho_obj_adcy_treg$hypermat)), 
                           scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                           xlab="", 
                           ylab="", #   textGrob("Regulatory T Cells", vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
                           pretty=TRUE) 
#RRHO2_heatmap(rrho2_obj_adcy_treg, colorGradient=col)

adcy_hsc_plot <-levelplot(rrho_obj_adcy_hsc$hypermat, col.regions = col_fun(min(rrho_obj_adcy_hsc$hypermat):max(rrho_obj_adcy_hsc$hypermat)), 
                          scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                          xlab="", 
                          ylab= "", # textGrob("Hematopoietic Stem Cells" , vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
                          pretty=TRUE) 
#RRHO2_heatmap(rrho2_obj_adcy_hsc, colorGradient=col)


adcy_wbm_plot <-levelplot(rrho_obj_adcy_wbm$hypermat, col.regions = col_fun(min(rrho_obj_adcy_wbm$hypermat):max(rrho_obj_adcy_wbm$hypermat)), 
                          scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                          xlab="", 
                          ylab= "", #textGrob("Whole Bone Marrow", vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
                          pretty=TRUE) 
#RRHO2_heatmap(rrho2_obj_adcy_wbm, colorGradient=col)

#_________________________________________________________________________________________________
#_________________________________________________________________________________________________

## take natural killer cell intersect


nkc_neut_intersect <- intersect(nkc_data$gene_symbol,neut_data$gene_symbol)
nkc_neut <- nkc_data[nkc_data$gene_symbol %in% nkc_neut_intersect, ]
neut_nkc <- neut_data[neut_data$gene_symbol %in% nkc_neut_intersect, ]

nkc_mp_intersect <- intersect(nkc_data$gene_symbol,mp_data$gene_symbol)
nkc_mp <- nkc_data[nkc_data$gene_symbol %in% nkc_mp_intersect, ]
mp_nkc <- mp_data[mp_data$gene_symbol %in% nkc_mp_intersect, ]

nkc_treg_intersect <- intersect(nkc_data$gene_symbol,treg_data$gene_symbol)
nkc_treg <- nkc_data[nkc_data$gene_symbol %in% nkc_treg_intersect, ]
treg_nkc <- treg_data[treg_data$gene_symbol %in% nkc_treg_intersect, ]


nkc_hsc_intersect <- intersect(nkc_data$gene_symbol,hsc_data$gene_symbol)
nkc_hsc <- nkc_data[nkc_data$gene_symbol %in% nkc_hsc_intersect, ]
hsc_nkc <- hsc_data[hsc_data$gene_symbol %in% nkc_hsc_intersect, ]

nkc_wbm_intersect <- intersect(nkc_data$gene_symbol,wbm_data$gene_symbol)
nkc_wbm <- nkc_data[nkc_data$gene_symbol %in% nkc_wbm_intersect, ]
wbm_nkc <- wbm_data[wbm_data$gene_symbol %in% nkc_wbm_intersect, ]


#_________________________________________________________________________________________________

## compute natural killer cell overlap and significance


rrho_obj_nkc_neut <- RRHO(nkc_neut, neut_nkc, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_nkc_neut <-  RRHO2_initialize(nkc_neut, neut_nkc, log10.ind=TRUE)
save(rrho_obj_nkc_neut, file = "data/nkc/rrho_obj_nkc_neut.RData")
load("data/nkc/rrho_obj_nkc_neut.RData")

rrho_obj_nkc_mp <- RRHO(nkc_mp, mp_nkc, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_nkc_mp <-  RRHO2_initialize(nkc_mp, mp_nkc, log10.ind=TRUE)
save(rrho_obj_nkc_mp, file = "data/nkc/rrho_obj_nkc_mp.RData")
load("data/nkc/rrho_obj_nkc_mp.RData")

rrho_obj_nkc_treg <- RRHO(nkc_treg, treg_nkc, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_nkc_treg <-  RRHO2_initialize(nkc_treg, treg_nkc, log10.ind=TRUE)
save(rrho_obj_nkc_treg, file = "data/nkc/rrho_obj_nkc_treg.RData")
load("data/nkc/rrho_obj_nkc_treg.RData")

rrho_obj_nkc_hsc <- RRHO(nkc_hsc, hsc_nkc, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_nkc_hsc <-  RRHO2_initialize(nkc_hsc, hsc_nkc, log10.ind=TRUE)
save(rrho_obj_nkc_hsc, file = "data/nkc/rrho_obj_nkc_hsc.RData")
load("data/nkc/rrho_obj_nkc_hsc.RData")

rrho_obj_nkc_wbm <- RRHO(nkc_wbm, wbm_nkc, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_nkc_wbm <-  RRHO2_initialize(nkc_wbm, wbm_nkc, log10.ind=TRUE)
save(rrho_obj_nkc_wbm, file = "data/nkc/rrho_obj_nkc_wbm.RData")
load("data/nkc/rrho_obj_nkc_wbm.RData")


#_________________________________________________________________________________________________

## plot natural killer cell heatmaps


nkc_neut_plot <-levelplot(rrho_obj_nkc_neut$hypermat, col.regions = col_fun(min(rrho_obj_nkc_neut$hypermat):max(rrho_obj_nkc_neut$hypermat)), 
                          scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                          xlab= textGrob("Natural Killer Cells", vjust = 1, gp = gpar(cex = 2)),  
                          ylab="", #  textGrob("Neutrophil", vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
                          pretty=TRUE) 
#RRHO2_heatmap(rrho2_obj_nkc_neut, colorGradient=col)

nkc_mp_plot <-levelplot(rrho_obj_nkc_mp$hypermat, col.regions = col_fun(min(rrho_obj_nkc_mp$hypermat):max(rrho_obj_nkc_mp$hypermat)), 
                        scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                        xlab="", 
                        ylab="", #  textGrob("Macrophage", vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
                        pretty=TRUE) 
#RRHO2_heatmap(rrho2_obj_nkc_mp, colorGradient=col)

nkc_treg_plot <-levelplot(rrho_obj_nkc_treg$hypermat, col.regions = col_fun(min(rrho_obj_nkc_treg$hypermat):max(rrho_obj_nkc_treg$hypermat)), 
                          scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                          xlab="", 
                          ylab="", #   textGrob("Regulatory T Cells", vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
                          pretty=TRUE) 
#RRHO2_heatmap(rrho2_obj_nkc_treg, colorGradient=col)

nkc_hsc_plot <-levelplot(rrho_obj_nkc_hsc$hypermat, col.regions = col_fun(min(rrho_obj_nkc_hsc$hypermat):max(rrho_obj_nkc_hsc$hypermat)), 
                         scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                         xlab="", 
                         ylab= "", # textGrob("Hematopoietic Stem Cells" , vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
                         pretty=TRUE) 
#RRHO2_heatmap(rrho2_obj_nkc_hsc, colorGradient=col)


nkc_wbm_plot <-levelplot(rrho_obj_nkc_wbm$hypermat, col.regions = col_fun(min(rrho_obj_nkc_wbm$hypermat):max(rrho_obj_nkc_wbm$hypermat)), 
                         scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                         xlab="", 
                         ylab="", #   textGrob("Whole Bone Marrow", vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
                         pretty=TRUE) 
#RRHO2_heatmap(rrho2_obj_nkc_wbm, colorGradient=col)


#_________________________________________________________________________________________________
#_________________________________________________________________________________________________

## take neutrophil intersect


neut_mp_intersect <- intersect(neut_data$gene_symbol,mp_data$gene_symbol)
neut_mp <- neut_data[neut_data$gene_symbol %in% neut_mp_intersect, ]
mp_neut <- mp_data[mp_data$gene_symbol %in% neut_mp_intersect, ]

neut_treg_intersect <- intersect(neut_data$gene_symbol,treg_data$gene_symbol)
neut_treg <- neut_data[neut_data$gene_symbol %in% neut_treg_intersect, ]
treg_neut <- treg_data[treg_data$gene_symbol %in% neut_treg_intersect, ]


neut_hsc_intersect <- intersect(neut_data$gene_symbol,hsc_data$gene_symbol)
neut_hsc <- neut_data[neut_data$gene_symbol %in% neut_hsc_intersect, ]
hsc_neut <- hsc_data[hsc_data$gene_symbol %in% neut_hsc_intersect, ]

neut_wbm_intersect <- intersect(neut_data$gene_symbol,wbm_data$gene_symbol)
neut_wbm <- neut_data[neut_data$gene_symbol %in% neut_wbm_intersect, ]
wbm_neut <- wbm_data[wbm_data$gene_symbol %in% neut_wbm_intersect, ]


#_________________________________________________________________________________________________

## compute neutrophil overlap and significance


rrho_obj_neut_mp <- RRHO(neut_mp, mp_neut, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_neut_mp <-  RRHO2_initialize(neut_mp, mp_neut, log10.ind=TRUE)
save(rrho_obj_neut_mp, file = "data/neut/rrho_obj_neut_mp.RData")
load("data/neut/rrho_obj_neut_mp.RData")

rrho_obj_neut_treg <- RRHO(neut_treg, treg_neut, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_neut_treg <-  RRHO2_initialize(neut_treg, treg_neut, log10.ind=TRUE)
save(rrho_obj_neut_treg, file = "data/neut/rrho_obj_neut_treg.RData")
load("data/neut/rrho_obj_neut_treg.RData")

rrho_obj_neut_hsc <- RRHO(neut_hsc, hsc_neut, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_neut_hsc <-  RRHO2_initialize(neut_hsc, hsc_neut, log10.ind=TRUE)
save(rrho_obj_neut_hsc, file = "data/neut/rrho_obj_neut_hsc.RData")
load("data/neut/rrho_obj_neut_hsc.RData")

rrho_obj_neut_wbm <- RRHO(neut_wbm, wbm_neut, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
#rrho2_obj_neut_wbm <-  RRHO2_initialize(neut_wbm, wbm_neut, log10.ind=TRUE)
save(rrho_obj_neut_wbm, file = "data/neut/rrho_obj_neut_wbm.RData")
load("data/neut/rrho_obj_neut_wbm.RData")


#_________________________________________________________________________________________________

## plot neutrophil heatmaps


neut_mp_plot <-levelplot(rrho_obj_neut_mp$hypermat, col.regions = col_fun(min(rrho_obj_neut_mp$hypermat):max(rrho_obj_neut_mp$hypermat)), 
                         scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                         xlab= textGrob("Neutrophil", vjust = 1, gp = gpar(cex = 2)),  
                         ylab="", #  textGrob("Macrophage", vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
                         pretty=TRUE) 
#RRHO2_heatmap(rrho2_obj_neut_mp, colorGradient=col)

neut_treg_plot <-levelplot(rrho_obj_neut_treg$hypermat, col.regions = col_fun(min(rrho_obj_neut_treg$hypermat):max(rrho_obj_neut_treg$hypermat)), 
                           scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                           xlab="", 
                           ylab="", #   textGrob("Regulatory T Cells", vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
                           pretty=TRUE) 
#RRHO2_heatmap(rrho2_obj_neut_treg, colorGradient=col)

neut_hsc_plot <-levelplot(rrho_obj_neut_hsc$hypermat, col.regions = col_fun(min(rrho_obj_neut_hsc$hypermat):max(rrho_obj_neut_hsc$hypermat)), 
                          scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                          xlab="", 
                          ylab="", #  textGrob("Hematopoietic Stem Cells" , vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
                          pretty=TRUE) 
#RRHO2_heatmap(rrho2_obj_neut_hsc, colorGradient=col)


neut_wbm_plot <-levelplot(rrho_obj_neut_wbm$hypermat, col.regions = col_fun(min(rrho_obj_neut_wbm$hypermat):max(rrho_obj_neut_wbm$hypermat)), 
                          scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                          xlab="", 
                          ylab="", #   textGrob("Whole Bone Marrow", vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
                          pretty=TRUE) 
#RRHO2_heatmap(rrho2_obj_neut_wbm, colorGradient=col)
#_________________________________________________________________________________________________
#_________________________________________________________________________________________________

## take macrophages intersect


mp_treg_intersect <- intersect(mp_data$gene_symbol,treg_data$gene_symbol)
mp_treg <- mp_data[mp_data$gene_symbol %in% mp_treg_intersect, ]
treg_mp <- treg_data[treg_data$gene_symbol %in% mp_treg_intersect, ]


mp_hsc_intersect <- intersect(mp_data$gene_symbol,hsc_data$gene_symbol)
mp_hsc <- mp_data[mp_data$gene_symbol %in% mp_hsc_intersect, ]
hsc_mp <- hsc_data[hsc_data$gene_symbol %in% mp_hsc_intersect, ]


mp_wbm_intersect <- intersect(mp_data$gene_symbol,wbm_data$gene_symbol)
mp_wbm <- mp_data[mp_data$gene_symbol %in% mp_wbm_intersect, ]
wbm_mp <- wbm_data[wbm_data$gene_symbol %in% mp_wbm_intersect, ]


#_________________________________________________________________________________________________

## compute macrophages overlap and significance


rrho_obj_mp_treg <- RRHO(mp_treg, treg_mp, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_mp_treg <-  RRHO2_initialize(mp_treg, treg_mp, log10.ind=TRUE)
save(rrho_obj_mp_treg, file = "data/mp/rrho_obj_mp_treg.RData")
load("data/mp/rrho_obj_mp_treg.RData")

rrho_obj_mp_hsc <- RRHO(mp_hsc, hsc_mp, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_mp_hsc <-  RRHO2_initialize(mp_hsc, hsc_mp, log10.ind=TRUE)
save(rrho_obj_mp_hsc, file = "data/mp/rrho_obj_mp_hsc.RData")
load("data/mp/rrho_obj_mp_hsc.RData")

rrho_obj_mp_wbm <- RRHO(mp_wbm, wbm_mp, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
#rrho2_obj_mp_wbm <-  RRHO2_initialize(mp_wbm, wbm_mp, log10.ind=TRUE)
save(rrho_obj_mp_wbm, file = "data/mp/rrho_obj_mp_wbm.RData")
load("data/mp/rrho_obj_mp_wbm.RData")


#_________________________________________________________________________________________________

## plot macrophages heatmaps


mp_treg_plot <-levelplot(rrho_obj_mp_treg$hypermat, col.regions = col_fun(min(rrho_obj_mp_treg$hypermat):max(rrho_obj_mp_treg$hypermat)), 
                         scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                         xlab= textGrob("Macrophage", vjust = 1, gp = gpar(cex = 2)),  
                         ylab=  "", # textGrob("Regulatory T Cells", vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
                         pretty=TRUE) 
#RRHO2_heatmap(rrho2_obj_mp_treg, colorGradient=col)

mp_hsc_plot <-levelplot(rrho_obj_mp_hsc$hypermat, col.regions = col_fun(min(rrho_obj_mp_hsc$hypermat):max(rrho_obj_mp_hsc$hypermat)), 
                        scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                        xlab="", 
                        ylab= "", # textGrob("Hematopoietic Stem Cells" , vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
                        pretty=TRUE) 
#RRHO2_heatmap(rrho2_obj_mp_hsc, colorGradient=col)


mp_wbm_plot <-levelplot(rrho_obj_mp_wbm$hypermat, col.regions = col_fun(min(rrho_obj_mp_wbm$hypermat):max(rrho_obj_mp_wbm$hypermat)), 
                        scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                        xlab="", 
                        ylab= "", #   textGrob("Whole Bone Marrow", vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
                        pretty=TRUE) 
#RRHO2_heatmap(rrho2_obj_mp_wbm, colorGradient=col)

#_________________________________________________________________________________________________
#_________________________________________________________________________________________________

## take regulatory t cell intersect


treg_hsc_intersect <- intersect(treg_data$gene_symbol,hsc_data$gene_symbol)
treg_hsc <- treg_data[treg_data$gene_symbol %in% treg_hsc_intersect, ]
hsc_treg <- hsc_data[hsc_data$gene_symbol %in% treg_hsc_intersect, ]

treg_wbm_intersect <- intersect(treg_data$gene_symbol,wbm_data$gene_symbol)
treg_wbm <- treg_data[treg_data$gene_symbol %in% treg_wbm_intersect, ]
wbm_treg <- wbm_data[wbm_data$gene_symbol %in% treg_wbm_intersect, ]


#_________________________________________________________________________________________________

## cotregute regulatory t cell overlap and significance


rrho_obj_treg_hsc <- RRHO(treg_hsc, hsc_treg, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_treg_hsc <-  RRHO2_initialize(treg_hsc, hsc_treg, log10.ind=TRUE)
save(rrho_obj_treg_hsc, file = "data/treg/rrho_obj_treg_hsc.RData")
load("data/treg/rrho_obj_treg_hsc.RData")

rrho_obj_treg_wbm <- RRHO(treg_wbm, wbm_treg, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_treg_wbm <-  RRHO2_initialize(treg_wbm, wbm_treg, log10.ind=TRUE)
save(rrho_obj_treg_wbm, file = "data/treg/rrho_obj_treg_wbm.RData")
load("data/treg/rrho_obj_treg_wbm.RData")


#_________________________________________________________________________________________________

## plot regulatory t cell heatmaps


treg_hsc_plot <-levelplot(rrho_obj_treg_hsc$hypermat, col.regions = col_fun(min(rrho_obj_treg_hsc$hypermat):max(rrho_obj_treg_hsc$hypermat)), 
                          scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                          xlab= textGrob("Regulatory T Cells", vjust = 1, gp = gpar(cex = 2)),  
                          ylab= "", # textGrob("Hematopoietic Stem Cells" , vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
                          pretty=TRUE) 
#RRHO2_heatmap(rrho2_obj_treg_hsc, colorGradient=col)


treg_wbm_plot <-levelplot(rrho_obj_treg_wbm$hypermat, col.regions = col_fun(min(rrho_obj_treg_wbm$hypermat):max(rrho_obj_treg_wbm$hypermat)), 
                          scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                          xlab="", 
                          ylab= "", #  textGrob("Whole Bone Marrow", vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
                          pretty=TRUE) 
#RRHO2_heatmap(rrho2_obj_treg_wbm, colorGradient=col)

#_________________________________________________________________________________________________
#_________________________________________________________________________________________________

## take hematopoietic stem cell intersect


hsc_wbm_intersect <- intersect(hsc_data$gene_symbol,wbm_data$gene_symbol)
hsc_wbm <- hsc_data[hsc_data$gene_symbol %in% hsc_wbm_intersect, ]
wbm_hsc <- wbm_data[wbm_data$gene_symbol %in% hsc_wbm_intersect, ]


#_________________________________________________________________________________________________

## cohscute hematopoietic stem cell overlap and significance


rrho_obj_hsc_wbm <- RRHO(hsc_wbm, wbm_hsc, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_hsc_wbm <-  RRHO2_initialize(hsc_wbm, wbm_hsc, log10.ind=TRUE)
save(rrho_obj_hsc_wbm, file = "data/hsc/rrho_obj_hsc_wbm.RData")
load("data/hsc/rrho_obj_hsc_wbm.RData")


#_________________________________________________________________________________________________

## plot hematopoietic stem cell heatmaps


hsc_wbm_plot <-levelplot(rrho_obj_hsc_wbm$hypermat, col.regions = col_fun(min(rrho_obj_hsc_wbm$hypermat):max(rrho_obj_hsc_wbm$hypermat)), 
                         scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                         xlab=textGrob("Hematopoietic Stem Cells", vjust = 1, gp = gpar(cex = 2)),  
                         ylab= "", #  textGrob("Whole Bone Marrow", vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
                         pretty=TRUE) 
#RRHO2_heatmap(rrho2_obj_hsc_wbm, colorGradient=col)



#_________________________________________________________________________________________________
#_________________________________________________________________________________________________

## arrange plots in a grid

tiff("results_rrho/pc_plot.tiff", units="in", width=10, height=50, res=300)
plot_grid(pc_wbm_plot,pc_mep_plot,pc_gmp_plot,pc_cmp_plot,pc_hsc_plot,pc_treg_plot,pc_mp_plot,pc_neut_plot,pc_nkc_plot,pc_adcy_plot,pc_opc_plot,pc_ocy_plot,pc_ec_plot,pc_bmsc_plot, align="v", ncol = 1)
dev.off()


tiff("results_rrho/bmsc_plot.tiff", units="in", width=10, height=50, res=300)
plot_grid(bmsc_wbm_plot,bmsc_mep_plot,bmsc_gmp_plot,bmsc_cmp_plot,bmsc_hsc_plot,bmsc_treg_plot,bmsc_mp_plot,bmsc_neut_plot,bmsc_nkc_plot,bmsc_adcy_plot,bmsc_opc_plot,bmsc_ocy_plot,bmsc_ec_plot,NULL, align="v", ncol = 1)
dev.off()

tiff("results_rrho/ec_plot.tiff", units="in", width=10, height=50, res=300)
plot_grid(ec_wbm_plot,ec_mep_plot,ec_gmp_plot,ec_cmp_plot,ec_hsc_plot,ec_treg_plot,ec_mp_plot,ec_neut_plot,ec_nkc_plot,ec_adcy_plot,ec_opc_plot,ec_ocy_plot,NULL,NULL, align="v", ncol = 1)
dev.off()


# combine figures
pc_plot <- plot_grid(pc_wbm_plot,pc_hsc_plot,pc_treg_plot,pc_mp_plot,pc_neut_plot,pc_nkc_plot,pc_adcy_plot,pc_opc_plot,pc_ocy_plot,pc_ec_plot,pc_bmsc_plot, align="v", ncol = 1)
save(pc_plot, file = "data/pc/pc_plot.RData")
bmsc_plot <- plot_grid(bmsc_wbm_plot,bmsc_hsc_plot,bmsc_treg_plot,bmsc_mp_plot,bmsc_neut_plot,bmsc_nkc_plot,bmsc_adcy_plot,bmsc_opc_plot,bmsc_ocy_plot,bmsc_ec_plot,NULL, align="v", ncol = 1)
ec_plot <-plot_grid(ec_wbm_plot,ec_hsc_plot,ec_treg_plot,ec_mp_plot,ec_neut_plot,ec_nkc_plot,ec_adcy_plot,ec_opc_plot,ec_ocy_plot,NULL,NULL, align="v", ncol = 1)
ocy_plot <-plot_grid(ocy_wbm_plot,ocy_hsc_plot,ocy_treg_plot,ocy_mp_plot,ocy_neut_plot,ocy_nkc_plot,ocy_adcy_plot,ocy_opc_plot,NULL,NULL,NULL, align="v", ncol = 1)
opc_plot <-plot_grid(opc_wbm_plot,opc_hsc_plot,opc_treg_plot,opc_mp_plot,opc_neut_plot,opc_nkc_plot,opc_adcy_plot,NULL,NULL,NULL,NULL, align="v", ncol = 1)
adcy_plot <-plot_grid(adcy_wbm_plot,adcy_hsc_plot,adcy_treg_plot,adcy_mp_plot,adcy_neut_plot,adcy_nkc_plot,NULL,NULL,NULL,NULL,NULL, align="v", ncol = 1)
nkc_plot <-plot_grid(nkc_wbm_plot,nkc_hsc_plot,nkc_treg_plot,nkc_mp_plot,nkc_neut_plot,NULL,NULL,NULL,NULL,NULL,NULL, align="v", ncol = 1)
neut_plot <-plot_grid(neut_wbm_plot,neut_hsc_plot,neut_treg_plot,neut_mp_plot,NULL,NULL,NULL,NULL,NULL,NULL,NULL, align="v", ncol = 1)
mp_plot <-plot_grid(mp_wbm_plot,mp_hsc_plot,mp_treg_plot,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL, align="v", ncol = 1)
treg_plot <-plot_grid(treg_wbm_plot,treg_hsc_plot,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL, align="v", ncol = 1)
hsc_plot <-plot_grid(hsc_wbm_plot,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL, align="v", ncol = 1)


tiff("results_rrho/rrho_plot.tiff", units="in", width=35, height=35, res=300)
plot_grid(pc_plot, bmsc_plot,ec_plot,ocy_plot,opc_plot, adcy_plot, nkc_plot,
          neut_plot, mp_plot, treg_plot, hsc_plot,
          align="h", ncol = 11)
dev.off()


#_________________________________________________________________________________________________
#_________________________________________________________________________________________________

## plot only significant ones 

pc_bmsc_plot <- levelplot(rrho_obj_pc_bmsc$hypermat, col.regions = col_fun(min(rrho_obj_pc_bmsc$hypermat):max(rrho_obj_pc_bmsc$hypermat)), 
                          scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                          xlab=textGrob("Bone Marrow Stromal Cells", vjust = 0.5, gp = gpar(cex = 2.2)), 
                          ylab= textGrob("Plasma Cells", vjust = 0.5, gp = gpar(cex = 2.2), rot = 90), 
                          pretty=TRUE, colorkey=list(labels=list(cex=1.5)))
pc_bmsc_plot

pc_treg_plot <-levelplot(rrho_obj_pc_treg$hypermat, col.regions = col_fun(min(rrho_obj_pc_treg$hypermat):max(rrho_obj_pc_treg$hypermat)), 
                         scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                         xlab=textGrob("Regulatory T Cells", vjust = 0.5, gp = gpar(cex = 2)), 
                         ylab="",#textGrob("Plasma Cells", vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
                         pretty=TRUE, colorkey=list(labels=list(cex=1.5)))
pc_treg_plot

pc_hsc_plot <-levelplot(rrho_obj_pc_hsc$hypermat, col.regions = col_fun(min(rrho_obj_pc_hsc$hypermat):max(rrho_obj_pc_hsc$hypermat)), 
                        scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                        xlab=textGrob("Hematopoietic Stem Cells", vjust = 0.5, gp = gpar(cex = 2)), 
                        ylab= "",#textGrob("Plasma Cells" , vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
                        pretty=TRUE, colorkey=list(labels=list(cex=1.5)))
pc_hsc_plot

pc_wbm_plot <-levelplot(rrho_obj_pc_wbm$hypermat, col.regions = col_fun(min(rrho_obj_pc_wbm$hypermat):max(rrho_obj_pc_wbm$hypermat)), 
                        scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                        xlab=textGrob("Whole Bone Marrow", vjust = 0.5, gp = gpar(cex = 2)), 
                        ylab="",#  textGrob("Plasma Cells", vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
                        pretty=TRUE, colorkey=list(labels=list(cex=1.5)))

pc_wbm_plot

bmsc_ec_plot <- levelplot(rrho_obj_bmsc_ec$hypermat, col.regions = col_fun(min(rrho_obj_bmsc_ec$hypermat):max(rrho_obj_bmsc_ec$hypermat)), 
                          scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                          xlab=textGrob("Endothelial Cells", vjust = 1, gp = gpar(cex = 2)), 
                          ylab= textGrob("Bone Marrow Stromal Cells", vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
                          pretty=TRUE, colorkey=list(labels=list(cex=1.5)))
bmsc_ec_plot

bmsc_opc_plot <- levelplot(rrho_obj_bmsc_opc$hypermat, col.regions = col_fun(min(rrho_obj_bmsc_opc$hypermat):max(rrho_obj_bmsc_opc$hypermat)), 
                           scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                           xlab=textGrob("Osteogenic Precursor Cells", vjust = 1, gp = gpar(cex = 2)), 
                           ylab= "",#textGrob("Bone Marrow Stromal Cells", vjust = 0.5, gp = gpar(cex = 2), rot = 90),  
                           pretty=TRUE, colorkey=list(labels=list(cex=1.5)))
bmsc_opc_plot

bmsc_wbm_plot <-levelplot(rrho_obj_bmsc_wbm$hypermat, col.regions = col_fun(min(rrho_obj_bmsc_wbm$hypermat):max(rrho_obj_bmsc_wbm$hypermat)), 
                          scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                          xlab=textGrob("Whole Bone Marrow", vjust = 1, gp = gpar(cex = 2)), 
                          ylab="",#textGrob("Bone Marrow Stromal Cells", vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
                          pretty=TRUE, colorkey=list(labels=list(cex=1.5)))
bmsc_wbm_plot

ec_opc_plot <- levelplot(rrho_obj_ec_opc$hypermat, col.regions = col_fun(min(rrho_obj_ec_opc$hypermat):max(rrho_obj_ec_opc$hypermat)), 
                         scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                         xlab= textGrob("Osteogenic Precursor Cells", vjust = 1, gp = gpar(cex = 2)),  
                         ylab= textGrob("Endothelial Cells", vjust = 0.5, gp = gpar(cex = 2), rot = 90),  
                         pretty=TRUE, colorkey=list(labels=list(cex=1.5)))
ec_opc_plot

ec_wbm_plot <-levelplot(rrho_obj_ec_wbm$hypermat, col.regions = col_fun(min(rrho_obj_ec_wbm$hypermat):max(rrho_obj_ec_wbm$hypermat)), 
                        scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                        xlab= textGrob("Whole Bone Marrow", vjust = 1, gp = gpar(cex = 2)),  
                        ylab= "",#textGrob("Endothelial Cells", vjust = 0.5, gp = gpar(cex = 2), rot = 90), 
                        pretty=TRUE, colorkey=list(labels=list(cex=1.5)))
ec_wbm_plot

# combine figures
pc_plot <- plot_grid(pc_bmsc_plot, pc_treg_plot, pc_hsc_plot, pc_wbm_plot, align="h", nrow = 1)
bmsc_plot <- plot_grid(bmsc_ec_plot, bmsc_opc_plot, bmsc_wbm_plot, NULL, align="h", nrow = 1)
ec_plot <-plot_grid(ec_opc_plot, ec_wbm_plot, NULL, NULL, align="h", nrow = 1)


tiff("results_rrho/rrho_plot_sign.tiff", units="in", width=20, height=15, res=300)
plot_grid(pc_plot, bmsc_plot,ec_plot,
          align="v", nrow = 3)
dev.off()

tiff("results_rrho/pc_plot_sign.tiff", units="in", width=20, height=7, res=100)
pc_plot
dev.off()

#_________________________________________________________________________________________________
#_________________________________________________________________________________________________

## assess overlap between human and mouse pcs

#_________________________________________________________________________________________________
# load data

mouse_pc_data <- read.csv("data/pc/mPC_MA.csv")
mouse_pc_data <- mouse_pc_data[,c(1,4,13,16)]
mouse_pc_data <- transform(mouse_pc_data, fisherPvalMin = pmin(fisherPvalUp, fisherPvalDown))
mouse_pc_data$score <- with(mouse_pc_data, ifelse(effectSize < 0, -log10(fisherPvalMin)*-1, -log10(fisherPvalMin)*1))
mouse_pc_data <- mouse_pc_data[,c(1,6)]
colnames(mouse_pc_data)[1] <- "gene_symbol"
mouse_pc_data <- mouse_pc_data[complete.cases(mouse_pc_data), ]


#_________________________________________________________________________________________________
# take intersect

## take plasma cell intersect
hs_mus_pc_intersect <- intersect(pc_data$gene_symbol,mouse_pc_data$gene_symbol)
hs_mus_pc <- pc_data[pc_data$gene_symbol %in% hs_mus_pc_intersect, ]
mus_hs_pc <- mouse_pc_data[mouse_pc_data$gene_symbol %in% hs_mus_pc_intersect, ]

#_________________________________________________________________________________________________

## compute overlap  significance
rrho_obj_hs_mus_pc <- RRHO(hs_mus_pc, mus_hs_pc, BY=TRUE, alternative='enrichment', log10.ind=TRUE)
# rrho2_obj_pc_bmsc <-  RRHO2_initialize(pc_bmsc, bmsc_pc, log10.ind=TRUE)
save(rrho_obj_hs_mus_pc, file = "data/rrho_obj_hs_mus_pc.RData")
load("data/pc/rrho_obj_hs_mus_pc.RData")

#_________________________________________________________________________________________________

## plot only significant ones 

hs_mus_pc_plot <- levelplot(rrho_obj_hs_mus_pc$hypermat, col.regions = col_fun(min(rrho_obj_hs_mus_pc$hypermat):max(rrho_obj_hs_mus_pc$hypermat)), 
                          scales= list(tck = c(0,0),y=(list(cex=0)), x=list(cex=0, rot=90)),
                          xlab=textGrob("Human plasma cells", vjust = 0.5, gp = gpar(cex = 2.2)), 
                          ylab= textGrob("Mouse plasma cells", vjust = 0.5, gp = gpar(cex = 2.2), rot = 90), 
                          pretty=TRUE, colorkey=list(labels=list(cex=1.5)))
hs_mus_pc_plot

#_________________________________________________________________________________________________

tiff("results_rrho/human_mouse_pc_plot.tiff", units="in", width=8, height=7, res=100)
hs_mus_pc_plot
dev.off()

#_________________________________________________________________________________________________
#_________________________________________________________________________________________________

## make summary plot of rrho plot

rrho_summary <- as.data.frame(matrix(nrow = 12, ncol = 12))
rownames(rrho_summary) <-  c("Plasma Cells","Bone Marrow Stromal Cells", "Endothelial Cells", "Osteocytes","Osteogenic Precursor Cells",
                             "Adipocytes", "Natural Killer Cells", "Neutrophils", "Macrophages",
                             "Regulatory T Cells", "Hematopoietic Stem Cells","Whole Bone Marrow")

colnames(rrho_summary) <-  c("Plasma Cells","Bone Marrow Stromal Cells", "Endothelial Cells", "Osteocytes","Osteogenic Precursor Cells",
                             "Adipocytes", "Natural Killer Cells", "Neutrophils", "Macrophages",
                             "Regulatory T Cells", "Hematopoietic Stem Cells","Whole Bone Marrow")



rrho_summary[]<-"Non-significant"
rrho_summary["Bone Marrow Stromal Cells", "Plasma Cells"] <- "Significant"
rrho_summary["Regulatory T Cells", "Plasma Cells"] <- "Significant"
rrho_summary["Hematopoietic Stem Cells", "Plasma Cells"] <- "Significant"
rrho_summary["Whole Bone Marrow", "Plasma Cells"] <- "Significant"
rrho_summary["Endothelial Cells", "Bone Marrow Stromal Cells"] <- "Significant"
rrho_summary["Osteogenic Precursor Cells", "Bone Marrow Stromal Cells"] <- "Significant"
rrho_summary["Whole Bone Marrow", "Bone Marrow Stromal Cells"] <- "Significant"
rrho_summary["Osteogenic Precursor Cells", "Endothelial Cells"] <- "Significant"
rrho_summary["Whole Bone Marrow", "Endothelial Cells"] <- "Significant"
diag(rrho_summary) <- NA
rrho_summary[upper.tri(rrho_summary)]=NA
rrho_summary <- as.matrix(rrho_summary)

colors = structure(c("#440154","#fde725"), names = c("Non-significant", "Significant")) # black, red, green, blue

rrho_summary_plot <-  ComplexHeatmap::Heatmap(rrho_summary, rect_gp = gpar(type = "none"), 
                          cluster_rows = FALSE, cluster_columns = FALSE, col = colors, row_names_side= "left",
                          height = unit(5, "cm"), width  = unit(5, "cm"),column_names_gp = gpar(fontsize = 9),
                          row_names_gp = gpar(fontsize = 9),
                          show_heatmap_legend = TRUE, 
                          heatmap_legend_param= list(title = "RRHO significance",direction = "horizontal", 
                                                     legend_width = unit(3, "cm"),
                                                     title_gp = gpar(fontsize = 9)),
                          cell_fun = function(j, i, x, y, width, height, fill) {
                            if(i >= j) {
                              grid.rect(x = x, y = y, width = width * 0.8, height = height*0.8, 
                                      gp = gpar(col = fill, fill = fill))
                            }
                              })

rrho_summary_plot
tiff("results_rrho/rrho_summary_plot.tiff", units="in", width=6, height=6, res=300)
rrho_summary_plot
dev.off()















