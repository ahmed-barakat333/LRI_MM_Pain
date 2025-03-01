# set working directory
setwd("~/Work/PhD/Projects/LRI MM Pain/Computational Analysis/Functional_Analysis")

# load packages
library(fgsea)
library(tidyverse)
library(tibble)
library(data.table)
library(ComplexHeatmap)
library(circlize)

#_________________________________________________________________________________________________
#_________________________________________________________________________________________________

# test enriched pathways in human pc

## load data
pc_data <- read.csv("data/pc/PC_MA.csv")
pc_data <- pc_data[,c(2,3,12,15)]
pc_data <- transform(pc_data, fisherPvalMin = pmin(fisherPvalUp, fisherPvalDown))
pc_data$score <- with(pc_data, ifelse(effectSize < 0, -log10(fisherPvalMin)*-1, -log10(fisherPvalMin)*1))
pc_data <- pc_data[,c(1,6)]
pc_data <- pc_data[complete.cases(pc_data), ]
pc_ranks <- deframe(pc_data)
head(pc_ranks, 20)

### test kegg enrichment
pathways_kegg <- gmtPathways("data/c2.cp.kegg.v2022.1.Hs.symbols.gmt")

pc_kegg <- fgsea(pathways=pathways_kegg, stats=pc_ranks)

pc_kegg <- pc_kegg[pc_kegg$padj < 0.05,]

pc_kegg <- pc_kegg %>%
  as_tibble() %>%
  arrange(desc(NES))

### test wikiPathways enrichment
pathways_wikiPathways <- gmtPathways("data/c2.cp.wikipathways.v2022.1.Hs.symbols.gmt")

pc_wikiPathways <- fgsea(pathways=pathways_wikiPathways, stats=pc_ranks)

pc_wikiPathways <- pc_wikiPathways[pc_wikiPathways$padj < 0.05,]

pc_wikiPathways <- pc_wikiPathways %>%
  as_tibble() %>%
  arrange(desc(NES))

### test reactome enrichment
pathways_reactome <- gmtPathways("data/c2.cp.reactome.v2022.1.Hs.symbols.gmt")

pc_reactome <- fgsea(pathways=pathways_reactome, stats=pc_ranks)

pc_reactome <- pc_reactome[pc_reactome$padj < 0.05,]

pc_reactome <- pc_reactome %>%
  as_tibble() %>%
  arrange(desc(NES))

pc_pathways <- rbind(pc_kegg, pc_reactome, pc_wikiPathways)

write_csv(pc_pathways,"results/PC_Pathways.csv")

#_________________________________________________________________________________________________
#_________________________________________________________________________________________________


# test enriched pathways in bmsc

## load data
bmsc_data <- read.csv("data/bmsc/bmsc_MA.csv")
bmsc_data <- bmsc_data[,c(2,3,12,15)]
bmsc_data <- transform(bmsc_data, fisherPvalMin = pmin(fisherPvalUp, fisherPvalDown))
bmsc_data$score <- with(bmsc_data, ifelse(effectSize < 0, -log10(fisherPvalMin)*-1, -log10(fisherPvalMin)*1))
bmsc_data <- bmsc_data[,c(1,6)]
bmsc_data <- bmsc_data[complete.cases(bmsc_data), ]

bmsc_ranks <- deframe(bmsc_data)
head(bmsc_ranks, 20)

### test kegg enrichment
bmsc_kegg <- fgsea(pathways=pathways_kegg, stats=bmsc_ranks)

bmsc_kegg <- bmsc_kegg[bmsc_kegg$padj < 0.05,]

bmsc_kegg <- bmsc_kegg %>%
  as_tibble() %>%
  arrange(desc(NES))

### test wikiPathways enrichment
bmsc_wikiPathways <- fgsea(pathways=pathways_wikiPathways, stats=bmsc_ranks)

bmsc_wikiPathways <- bmsc_wikiPathways[bmsc_wikiPathways$padj < 0.05,]

bmsc_wikiPathways <- bmsc_wikiPathways %>%
  as_tibble() %>%
  arrange(desc(NES))

### test reactome enrichment
bmsc_reactome <- fgsea(pathways=pathways_reactome, stats=bmsc_ranks)

bmsc_reactome <- bmsc_reactome[bmsc_reactome$padj < 0.05,]

bmsc_reactome <- bmsc_reactome %>%
  as_tibble() %>%
  arrange(desc(NES))

bmsc_pathways <- rbind(bmsc_kegg, bmsc_reactome, bmsc_wikiPathways)

write_csv(bmsc_pathways,"results/BMSC_Pathways.csv")
#_________________________________________________________________________________________________
#_________________________________________________________________________________________________

# test enriched pathways in ec

## load data
ec_data <- read_csv("data/ec/GSE14230 DEG.ID.csv")
ec_data <- ec_data[,c(3,6,9)]
ec_data$score <- with(ec_data, ifelse(logFC < 0, -log10(P.Value)*-1, -log10(P.Value)*1))
ec_data <- ec_data[,c(1,4)]
ec_data <- ec_data[complete.cases(ec_data), ]

sum(duplicated(ec_data$SYMBOL))

ec_ranks <- deframe(ec_data)
head(ec_ranks, 20)

### test kegg enrichment
ec_kegg <- fgsea(pathways=pathways_kegg, stats=ec_ranks)

ec_kegg <- ec_kegg[ec_kegg$padj < 0.05,]

ec_kegg <- ec_kegg %>%
  as_tibble() %>%
  arrange(desc(NES))

### test wikiPathways enrichment
ec_wikiPathways <- fgsea(pathways=pathways_wikiPathways, stats=ec_ranks)

ec_wikiPathways <- ec_wikiPathways[ec_wikiPathways$padj < 0.05,]

ec_wikiPathways <- ec_wikiPathways %>%
  as_tibble() %>%
  arrange(desc(NES))

### test reactome enrichment
ec_reactome <- fgsea(pathways=pathways_reactome, stats=ec_ranks)

ec_reactome <- ec_reactome[ec_reactome$padj < 0.05,]

ec_reactome <- ec_reactome %>%
  as_tibble() %>%
  arrange(desc(NES))

ec_pathways <- rbind(ec_kegg, ec_reactome, ec_wikiPathways)

write_csv(ec_pathways,"results/EC_Pathways.csv")
#_________________________________________________________________________________________________
#_________________________________________________________________________________________________

# test enriched pathways in ocy

## load data
ocy_data <- read_csv("data/ocy/GSE27372 DEG.ID.csv")
ocy_data <- ocy_data[,c(3,6,9)]
ocy_data$score <- with(ocy_data, ifelse(logFC < 0, -log10(P.Value)*-1, -log10(P.Value)*1))
ocy_data <- ocy_data[,c(1,4)]
ocy_data <- ocy_data[complete.cases(ocy_data), ]

sum(duplicated(ocy_data$SYMBOL))

ocy_ranks <- deframe(ocy_data)
head(ocy_ranks, 20)

### test kegg enrichment
ocy_kegg <- fgsea(pathways=pathways_kegg, stats=ocy_ranks)

ocy_kegg <- ocy_kegg[ocy_kegg$padj < 0.05,]

ocy_kegg <- ocy_kegg %>%
  as_tibble() %>%
  arrange(desc(NES))

### test wikiPathways enrichment
ocy_wikiPathways <- fgsea(pathways=pathways_wikiPathways, stats=ocy_ranks)

ocy_wikiPathways <- ocy_wikiPathways[ocy_wikiPathways$padj < 0.05,]

ocy_wikiPathways <- ocy_wikiPathways %>%
  as_tibble() %>%
  arrange(desc(NES))

### test reactome enrichment
ocy_reactome <- fgsea(pathways=pathways_reactome, stats=ocy_ranks)

ocy_reactome <- ocy_reactome[ocy_reactome$padj < 0.05,]

ocy_reactome <- ocy_reactome %>%
  as_tibble() %>%
  arrange(desc(NES))

ocy_pathways <- rbind(ocy_kegg, ocy_reactome, ocy_wikiPathways)

write_csv(ocy_pathways,"results/OCY_Pathways.csv")


#_________________________________________________________________________________________________
#_________________________________________________________________________________________________

# test enriched pathways in opc

## load data
opc_data <- read_csv("data/opc/GSE87073 DEG.ID.csv")
opc_data <- opc_data[,c(3,6,9)]
opc_data$score <- with(opc_data, ifelse(logFC < 0, -log10(P.Value)*-1, -log10(P.Value)*1))
opc_data <- opc_data[,c(1,4)]
opc_data <- opc_data[complete.cases(opc_data), ]

sum(duplicated(opc_data$SYMBOL))

opc_ranks <- deframe(opc_data)
head(opc_ranks, 20)

### test kegg enrichment
opc_kegg <- fgsea(pathways=pathways_kegg, stats=opc_ranks)

opc_kegg <- opc_kegg[opc_kegg$padj < 0.05,]

opc_kegg <- opc_kegg %>%
  as_tibble() %>%
  arrange(desc(NES))

### test wikiPathways enrichment
opc_wikiPathways <- fgsea(pathways=pathways_wikiPathways, stats=opc_ranks)

opc_wikiPathways <- opc_wikiPathways[opc_wikiPathways$padj < 0.05,]

opc_wikiPathways <- opc_wikiPathways %>%
  as_tibble() %>%
  arrange(desc(NES))

### test reactome enrichment
opc_reactome <- fgsea(pathways=pathways_reactome, stats=opc_ranks)

opc_reactome <- opc_reactome[opc_reactome$padj < 0.05,]

opc_reactome <- opc_reactome %>%
  as_tibble() %>%
  arrange(desc(NES))

opc_pathways <- rbind(opc_kegg, opc_reactome, opc_wikiPathways)

write_csv(opc_pathways,"results/OPC_Pathways.csv")

#_________________________________________________________________________________________________
#_________________________________________________________________________________________________

# test enriched pathways in adcy

## load data
adcy_data <- read_csv("data/adcy/GSE143269_DEG.ID.HS.Mapped.csv")
adcy_data <- adcy_data[,c(10,4,7)]
adcy_data$score <- with(adcy_data, ifelse(logFC < 0, -log10(P.Value)*-1, -log10(P.Value)*1))
adcy_data <- adcy_data[,c(1,4)]
adcy_data <- adcy_data[complete.cases(adcy_data), ]

sum(duplicated(adcy_data$HGNC.symbol))
adcy_data <- data.frame(adcy_data %>% 
                          group_by(HGNC.symbol) %>% 
                          summarise_all(funs(max)))
adcy_ranks <- deframe(adcy_data)
head(adcy_ranks, 20)

### test kegg enrichment
adcy_kegg <- fgsea(pathways=pathways_kegg, stats=adcy_ranks)

adcy_kegg <- adcy_kegg[adcy_kegg$padj < 0.05,]

adcy_kegg <- adcy_kegg %>%
  as_tibble() %>%
  arrange(desc(NES))

### test wikiPathways enrichment
adcy_wikiPathways <- fgsea(pathways=pathways_wikiPathways, stats=adcy_ranks)

adcy_wikiPathways <- adcy_wikiPathways[adcy_wikiPathways$padj < 0.05,]

adcy_wikiPathways <- adcy_wikiPathways %>%
  as_tibble() %>%
  arrange(desc(NES))

### test reactome enrichment
adcy_reactome <- fgsea(pathways=pathways_reactome, stats=adcy_ranks)

adcy_reactome <- adcy_reactome[adcy_reactome$padj < 0.05,]

adcy_reactome <- adcy_reactome %>%
  as_tibble() %>%
  arrange(desc(NES))

adcy_pathways <- rbind(adcy_kegg, adcy_reactome, adcy_wikiPathways)

write_csv(adcy_pathways,"results/ADCY_Pathways.csv")

#_________________________________________________________________________________________________
#_________________________________________________________________________________________________

# test enriched pathways in nkc

## load data
nkc_data <- read_csv("data/nkc/GSE27838 DEG.ID.csv")
nkc_data <- nkc_data[,c(3,6,9)]
nkc_data$score <- with(nkc_data, ifelse(logFC < 0, -log10(P.Value)*-1, -log10(P.Value)*1))
nkc_data <- nkc_data[,c(1,4)]
nkc_data <- nkc_data[complete.cases(nkc_data), ]

sum(duplicated(nkc_data$SYMBOL))

nkc_ranks <- deframe(nkc_data)
head(nkc_ranks, 20)

### test kegg enrichment
nkc_kegg <- fgsea(pathways=pathways_kegg, stats=nkc_ranks)

nkc_kegg <- nkc_kegg[nkc_kegg$padj < 0.05,]

nkc_kegg <- nkc_kegg %>%
  as_tibble() %>%
  arrange(desc(NES))

### test wikiPathways enrichment
nkc_wikiPathways <- fgsea(pathways=pathways_wikiPathways, stats=nkc_ranks)

nkc_wikiPathways <- nkc_wikiPathways[nkc_wikiPathways$padj < 0.05,]

nkc_wikiPathways <- nkc_wikiPathways %>%
  as_tibble() %>%
  arrange(desc(NES))

### test reactome enrichment
nkc_reactome <- fgsea(pathways=pathways_reactome, stats=nkc_ranks)

nkc_reactome <- nkc_reactome[nkc_reactome$padj < 0.05,]

nkc_reactome <- nkc_reactome %>%
  as_tibble() %>%
  arrange(desc(NES))

nkc_pathways <- rbind(nkc_kegg, nkc_reactome, nkc_wikiPathways)

write_csv(nkc_pathways,"results/NKC_Pathways.csv")

#_________________________________________________________________________________________________
#_________________________________________________________________________________________________
# test enriched pathways in neut

## load data
neut_data <- read_csv("data/neut/GSE150021_DEG_limma.csv")
neut_data <- neut_data[,c(2,3,6)]
neut_data$score <- with(neut_data, ifelse(logFC < 0, -log10(P.Value)*-1, -log10(P.Value)*1))
neut_data <- neut_data[,c(1,4)]
neut_data <- neut_data[complete.cases(neut_data), ]
colnames(neut_data)[1] <- "SYMBOL" 

sum(duplicated(neut_data$SYMBOL))
neut_data <- data.frame(neut_data %>% 
                          group_by(SYMBOL) %>% 
                          summarise_all(funs(max)))


neut_ranks <- deframe(neut_data)
head(neut_ranks, 20)

### test kegg enrichment
neut_kegg <- fgsea(pathways=pathways_kegg, stats=neut_ranks)

neut_kegg <- neut_kegg[neut_kegg$padj < 0.05,]

neut_kegg <- neut_kegg %>%
  as_tibble() %>%
  arrange(desc(NES))

### test wikiPathways enrichment
neut_wikiPathways <- fgsea(pathways=pathways_wikiPathways, stats=neut_ranks)

neut_wikiPathways <- neut_wikiPathways[neut_wikiPathways$padj < 0.05,]

neut_wikiPathways <- neut_wikiPathways %>%
  as_tibble() %>%
  arrange(desc(NES))

### test reactome enrichment
neut_reactome <- fgsea(pathways=pathways_reactome, stats=neut_ranks)

neut_reactome <- neut_reactome[neut_reactome$padj < 0.05,]

neut_reactome <- neut_reactome %>%
  as_tibble() %>%
  arrange(desc(NES))

neut_pathways <- rbind(neut_kegg, neut_reactome, neut_wikiPathways)

write_csv(neut_pathways,"results/NEUT_Pathways.csv")

#_________________________________________________________________________________________________
#_________________________________________________________________________________________________

# test enriched pathways in macrophages

## load data
mp_data <- read_csv("data/mp/MP_MA.csv")
mp_data <- mp_data[,c(2,4,13,16)]
mp_data <- transform(mp_data, fisherPvalMin = pmin(fisherPvalUp, fisherPvalDown))
mp_data$score <- with(mp_data, ifelse(effectSize < 0, -log10(fisherPvalMin)*-1, -log10(fisherPvalMin)*1))
mp_data <- mp_data[,c(1,6)]
mp_data <- mp_data[complete.cases(mp_data), ]
colnames(mp_data)[1] <- "SYMBOL" 

sum(duplicated(mp_data$SYMBOL))
mp_data <- data.frame(mp_data %>% 
                          group_by(SYMBOL) %>% 
                          summarise_all(funs(max)))

mp_ranks <- deframe(mp_data)
head(mp_ranks, 20)

### test kegg enrichment
mp_kegg <- fgsea(pathways=pathways_kegg, stats=mp_ranks)

mp_kegg <- mp_kegg[mp_kegg$padj < 0.05,]

mp_kegg <- mp_kegg %>%
  as_tibble() %>%
  arrange(desc(NES))

### test wikiPathways enrichment
mp_wikiPathways <- fgsea(pathways=pathways_wikiPathways, stats=mp_ranks)

mp_wikiPathways <- mp_wikiPathways[mp_wikiPathways$padj < 0.05,]

mp_wikiPathways <- mp_wikiPathways %>%
  as_tibble() %>%
  arrange(desc(NES))

### test reactome enrichment
mp_reactome <- fgsea(pathways=pathways_reactome, stats=mp_ranks)

mp_reactome <- mp_reactome[mp_reactome$padj < 0.05,]

mp_reactome <- mp_reactome %>%
  as_tibble() %>%
  arrange(desc(NES))

mp_pathways <- rbind(mp_kegg, mp_reactome, mp_wikiPathways)

write_csv(mp_pathways,"results/MP_Pathways.csv")

#_________________________________________________________________________________________________
#_________________________________________________________________________________________________

# test enriched pathways in treg

## load data
treg_data <- read_csv("data/treg/GSE109533_DEG_limma.csv")
treg_data <- treg_data[,c(11,4,7)]
treg_data$score <- with(treg_data, ifelse(logFC < 0, -log10(P.Value)*-1, -log10(P.Value)*1))
treg_data <- treg_data[,c(1,4)]
treg_data <- treg_data[complete.cases(treg_data), ]
colnames(treg_data)[1] <- "SYMBOL" 

sum(duplicated(treg_data$SYMBOL))
treg_data <- data.frame(treg_data %>% 
                        group_by(SYMBOL) %>% 
                        summarise_all(funs(max)))

treg_ranks <- deframe(treg_data)
head(treg_ranks, 20)

### test kegg enrichment
treg_kegg <- fgsea(pathways=pathways_kegg, stats=treg_ranks)

treg_kegg <- treg_kegg[treg_kegg$padj < 0.05,]

treg_kegg <- treg_kegg %>%
  as_tibble() %>%
  arrange(desc(NES))

### test wikiPathways enrichment
treg_wikiPathways <- fgsea(pathways=pathways_wikiPathways, stats=treg_ranks)

treg_wikiPathways <- treg_wikiPathways[treg_wikiPathways$padj < 0.05,]

treg_wikiPathways <- treg_wikiPathways %>%
  as_tibble() %>%
  arrange(desc(NES))

### test reactome enrichment
treg_reactome <- fgsea(pathways=pathways_reactome, stats=treg_ranks)

treg_reactome <- treg_reactome[treg_reactome$padj < 0.05,]

treg_reactome <- treg_reactome %>%
  as_tibble() %>%
  arrange(desc(NES))

treg_pathways <- rbind(treg_kegg, treg_reactome, treg_wikiPathways)

write_csv(treg_pathways,"results/Treg_Pathways.csv")

#_________________________________________________________________________________________________
#_________________________________________________________________________________________________
# test enriched pathways in hsc

## load data
hsc_data <- read_csv("data/hsc/GSE24870_HSC DEG.ID.csv")
hsc_data <- hsc_data[,c(3,6,9)]
hsc_data$score <- with(hsc_data, ifelse(logFC < 0, -log10(P.Value)*-1, -log10(P.Value)*1))
hsc_data <- hsc_data[,c(1,4)]
hsc_data <- hsc_data[complete.cases(hsc_data), ]

sum(duplicated(hsc_data$SYMBOL))

hsc_ranks <- deframe(hsc_data)
head(hsc_ranks, 20)

### test kegg enrichment
hsc_kegg <- fgsea(pathways=pathways_kegg, stats=hsc_ranks)

hsc_kegg <- hsc_kegg[hsc_kegg$padj < 0.05,]

hsc_kegg <- hsc_kegg %>%
  as_tibble() %>%
  arrange(desc(NES))

### test wikiPathways enrichment
hsc_wikiPathways <- fgsea(pathways=pathways_wikiPathways, stats=hsc_ranks)

hsc_wikiPathways <- hsc_wikiPathways[hsc_wikiPathways$padj < 0.05,]

hsc_wikiPathways <- hsc_wikiPathways %>%
  as_tibble() %>%
  arrange(desc(NES))

### test reactome enrichment
hsc_reactome <- fgsea(pathways=pathways_reactome, stats=hsc_ranks)

hsc_reactome <- hsc_reactome[hsc_reactome$padj < 0.05,]

hsc_reactome <- hsc_reactome %>%
  as_tibble() %>%
  arrange(desc(NES))

hsc_pathways <- rbind(hsc_kegg, hsc_reactome, hsc_wikiPathways)

write_csv(hsc_pathways,"results/HSC_Pathways.csv")


#_________________________________________________________________________________________________
#_________________________________________________________________________________________________
# test enriched pathways in wbm

## load data
wbm_data <- read_csv("data/wbm/GSE118985 DEG.ID.SYMBOL.csv")
wbm_data <- wbm_data[,c(4,21,24)]
wbm_data$score <- with(wbm_data, ifelse(logFC < 0, -log10(P.Value)*-1, -log10(P.Value)*1))
wbm_data <- wbm_data[,c(1,4)]
wbm_data <- wbm_data[complete.cases(wbm_data), ]

colnames(wbm_data)[1] <- "SYMBOL" 

sum(duplicated(wbm_data$SYMBOL))
wbm_data <- data.frame(wbm_data %>% 
                          group_by(SYMBOL) %>% 
                          summarise_all(funs(max)))

wbm_ranks <- deframe(wbm_data)
head(wbm_ranks, 20)

### test kegg enrichment
wbm_kegg <- fgsea(pathways=pathways_kegg, stats=wbm_ranks)

wbm_kegg <- wbm_kegg[wbm_kegg$padj < 0.05,]

wbm_kegg <- wbm_kegg %>%
  as_tibble() %>%
  arrange(desc(NES))

### test wikiPathways enrichment
wbm_wikiPathways <- fgsea(pathways=pathways_wikiPathways, stats=wbm_ranks)

wbm_wikiPathways <- wbm_wikiPathways[wbm_wikiPathways$padj < 0.05,]

wbm_wikiPathways <- wbm_wikiPathways %>%
  as_tibble() %>%
  arrange(desc(NES))

### test reactome enrichment
wbm_reactome <- fgsea(pathways=pathways_reactome, stats=wbm_ranks)

wbm_reactome <- wbm_reactome[wbm_reactome$padj < 0.05,]

wbm_reactome <- wbm_reactome %>%
  as_tibble() %>%
  arrange(desc(NES))

wbm_pathways <- rbind(wbm_kegg, wbm_reactome, wbm_wikiPathways)

write_csv(wbm_pathways,"results/WBM_Pathways.csv")

#_________________________________________________________________________________________________
#_________________________________________________________________________________________________

# heatmap of kegg pathways
pc_kegg <- pc_kegg[c(1,6)]
colnames(pc_kegg)[2] <- "Plasma Cells"
pc_kegg$pathway <- gsub('KEGG_', '', pc_kegg$pathway)
pc_kegg$pathway <- gsub('_', ' ', pc_kegg$pathway)


bmsc_kegg <- bmsc_kegg[c(1,6)]
colnames(bmsc_kegg)[2] <- "Bone Marrow Stromal Cells"
bmsc_kegg$pathway <- gsub('KEGG_', '', bmsc_kegg$pathway)
bmsc_kegg$pathway <- gsub('_', ' ', bmsc_kegg$pathway)

ec_kegg <- ec_kegg[c(1,6)]
colnames(ec_kegg)[2] <- "Endothelial Cells"
ec_kegg$pathway <- gsub('KEGG_', '', ec_kegg$pathway)
ec_kegg$pathway <- gsub('_', ' ', ec_kegg$pathway)

ocy_kegg <- ocy_kegg[c(1,6)]
colnames(ocy_kegg)[2] <- "Osteocytes"
ocy_kegg$pathway <- gsub('KEGG_', '', ocy_kegg$pathway)
ocy_kegg$pathway <- gsub('_', ' ', ocy_kegg$pathway)

opc_kegg <- opc_kegg[c(1,6)]
colnames(opc_kegg)[2] <- "Osteogenic Precursor Cells"
opc_kegg$pathway <- gsub('KEGG_', '', opc_kegg$pathway)
opc_kegg$pathway <- gsub('_', ' ', opc_kegg$pathway)

adcy_kegg <- adcy_kegg[c(1,6)]
colnames(adcy_kegg)[2] <- "Adipocytes"
adcy_kegg$pathway <- gsub('KEGG_', '', adcy_kegg$pathway)
adcy_kegg$pathway <- gsub('_', ' ', adcy_kegg$pathway)

nkc_kegg <- nkc_kegg[c(1,6)]
colnames(nkc_kegg)[2] <- "Natural Killer Cells"
nkc_kegg$pathway <- gsub('KEGG_', '', nkc_kegg$pathway)
nkc_kegg$pathway <- gsub('_', ' ', nkc_kegg$pathway)

neut_kegg <- neut_kegg[c(1,6)]
colnames(neut_kegg)[2] <- "Neutrophils"
neut_kegg$pathway <- gsub('KEGG_', '', neut_kegg$pathway)
neut_kegg$pathway <- gsub('_', ' ', neut_kegg$pathway)

mp_kegg <- mp_kegg[c(1,6)]
colnames(mp_kegg)[2] <- "Macrophages"
mp_kegg$pathway <- gsub('KEGG_', '', mp_kegg$pathway)
mp_kegg$pathway <- gsub('_', ' ', mp_kegg$pathway)

treg_kegg <- treg_kegg[c(1,6)]
colnames(treg_kegg)[2] <- "Regulatory T Cells"
treg_kegg$pathway <- gsub('KEGG_', '', treg_kegg$pathway)
treg_kegg$pathway <- gsub('_', ' ', treg_kegg$pathway)

hsc_kegg <- hsc_kegg[c(1,6)]
colnames(hsc_kegg)[2] <- "Hematopoietic Stem Cells"
hsc_kegg$pathway <- gsub('KEGG_', '', hsc_kegg$pathway)
hsc_kegg$pathway <- gsub('_', ' ', hsc_kegg$pathway)

wbm_kegg <- wbm_kegg[c(1,6)]
colnames(wbm_kegg)[2] <- "Whole Bone Marrow"
wbm_kegg$pathway <- gsub('KEGG_', '', wbm_kegg$pathway)
wbm_kegg$pathway <- gsub('_', ' ', wbm_kegg$pathway)

kegg <- list(pc_kegg, bmsc_kegg, ec_kegg, ocy_kegg, opc_kegg, adcy_kegg,nkc_kegg,
             neut_kegg, mp_kegg, treg_kegg, hsc_kegg, wbm_kegg)

kegg <- kegg %>% reduce(full_join, by='pathway')
kegg <- kegg %>% remove_rownames %>% column_to_rownames(var="pathway")

kegg$sum.na <- rowSums(is.na(kegg))
hist(kegg$sum.na)
kegg <- kegg[kegg$sum.na < 8,]
kegg$sum.na <- NULL
kegg[is.na(kegg)] <- 0



kegg_heatmap <- as.matrix(kegg)
col_fun = colorRamp2(c(-4, 0, 4), c("blue", "#e6e6e6", "red"))
ht_opt$TITLE_PADDING = unit(c(3.5, 3.5), "points")
kegg_plot <-Heatmap(kegg_heatmap, col=col_fun, row_names_side= "right",row_dend_side= "left", name="NES", 
                    row_names_gp = gpar(fontsize = 9), #column_names_max_height = max_text_height(colnames(pathways)),
                    row_title = "KEGG", row_title_gp = gpar(fontsize = 9, fill = "#33cc33", col = "white", border = NA),
                    rect_gp = gpar(col = "white", lwd = 2),
                    width = unit(3.96, "cm"), height = unit(2.97, "cm"),column_names_gp = gpar(fontsize = 9), 
                    show_column_dend = FALSE, #row_dend_width = unit(0.2, "cm") 
                    show_row_dend = FALSE,#row_title = "Multiple Myeloma Cells",
                    column_order = colnames(kegg),
                    show_heatmap_legend = FALSE)
kegg_plot

tiff("results/kegg_plot_2.tiff", units="in", width=8, height=8, res=100)
kegg_plot
dev.off()


kegg <- kegg[order(row.names(kegg)), ]
kegg_heatmap <- as.matrix(kegg)

kegg_heatmap <- t(kegg_heatmap)

col_fun = colorRamp2(c(-4, 0, 4), c("blue", "#e6e6e6", "red"))
ht_opt$TITLE_PADDING = unit(c(1, 1), "points")
kegg_plot <-Heatmap(kegg_heatmap, col=col_fun, row_names_side= "left",row_dend_side= "left", name="NES", 
        row_names_gp = gpar(fontsize = 8), #column_names_max_height = max_text_height(colnames(pathways)),
        column_title = "KEGG", column_title_gp = gpar(fill = "#33cc33", col = "white", border = NA),
        rect_gp = gpar(col = "white", lwd = 2),
        width = unit(3.96, "cm"), height = unit(4, "cm"),column_names_gp = gpar(fontsize = 8), 
        show_column_dend = FALSE, #row_dend_width = unit(0.2, "cm") 
        show_row_dend = FALSE,#row_title = "Multiple Myeloma Cells",
        row_order = colnames(kegg),
        row_title_gp = gpar(fontsize = 8),
        show_heatmap_legend = FALSE)
kegg_plot

tiff("results/kegg_plot.tiff", units="in", width=8, height=8, res=100)
kegg_plot
dev.off()

#_________________________________________________________________________________________________
#_________________________________________________________________________________________________

# heatmap of reactome pathways
pc_reactome <- pc_reactome[c(1,6)]
colnames(pc_reactome)[2] <- "Plasma Cells"
pc_reactome$pathway <- gsub('REACTOME_', '', pc_reactome$pathway)
pc_reactome$pathway <- gsub('_', ' ', pc_reactome$pathway)

bmsc_reactome <- bmsc_reactome[c(1,6)]
colnames(bmsc_reactome)[2] <- "Bone Marrow Stromal Cells"
bmsc_reactome$pathway <- gsub('REACTOME_', '', bmsc_reactome$pathway)
bmsc_reactome$pathway <- gsub('_', ' ', bmsc_reactome$pathway)

ec_reactome <- ec_reactome[c(1,6)]
colnames(ec_reactome)[2] <- "Endothelial Cells"
ec_reactome$pathway <- gsub('REACTOME_', '', ec_reactome$pathway)
ec_reactome$pathway <- gsub('_', ' ', ec_reactome$pathway)

ocy_reactome <- ocy_reactome[c(1,6)]
colnames(ocy_reactome)[2] <- "Osteocytes"
ocy_reactome$pathway <- gsub('REACTOME_', '', ocy_reactome$pathway)
ocy_reactome$pathway <- gsub('_', ' ', ocy_reactome$pathway)

opc_reactome <- opc_reactome[c(1,6)]
colnames(opc_reactome)[2] <- "Osteogenic Precursor Cells"
opc_reactome$pathway <- gsub('REACTOME_', '', opc_reactome$pathway)
opc_reactome$pathway <- gsub('_', ' ', opc_reactome$pathway)

adcy_reactome <- adcy_reactome[c(1,6)]
colnames(adcy_reactome)[2] <- "Adipocytes"
adcy_reactome$pathway <- gsub('REACTOME_', '', adcy_reactome$pathway)
adcy_reactome$pathway <- gsub('_', ' ', adcy_reactome$pathway)

nkc_reactome <- nkc_reactome[c(1,6)]
colnames(nkc_reactome)[2] <- "Natural Killer Cells"
nkc_reactome$pathway <- gsub('REACTOME_', '', nkc_reactome$pathway)
nkc_reactome$pathway <- gsub('_', ' ', nkc_reactome$pathway)

neut_reactome <- neut_reactome[c(1,6)]
colnames(neut_reactome)[2] <- "Neutrophils"
neut_reactome$pathway <- gsub('REACTOME_', '', neut_reactome$pathway)
neut_reactome$pathway <- gsub('_', ' ', neut_reactome$pathway)

mp_reactome <- mp_reactome[c(1,6)]
colnames(mp_reactome)[2] <- "Macrophages"
mp_reactome$pathway <- gsub('REACTOME_', '', mp_reactome$pathway)
mp_reactome$pathway <- gsub('_', ' ', mp_reactome$pathway)

treg_reactome <- treg_reactome[c(1,6)]
colnames(treg_reactome)[2] <- "Regulatory T Cells"
treg_reactome$pathway <- gsub('REACTOME_', '', treg_reactome$pathway)
treg_reactome$pathway <- gsub('_', ' ', treg_reactome$pathway)

hsc_reactome <- hsc_reactome[c(1,6)]
colnames(hsc_reactome)[2] <- "Hematopoietic Stem Cells"
hsc_reactome$pathway <- gsub('REACTOME_', '', hsc_reactome$pathway)
hsc_reactome$pathway <- gsub('_', ' ', hsc_reactome$pathway)

wbm_reactome <- wbm_reactome[c(1,6)]
colnames(wbm_reactome)[2] <- "Whole Bone Marrow"
wbm_reactome$pathway <- gsub('REACTOME_', '', wbm_reactome$pathway)
wbm_reactome$pathway <- gsub('_', ' ', wbm_reactome$pathway)


reactome <- list(pc_reactome, bmsc_reactome, ec_reactome, ocy_reactome, opc_reactome, adcy_reactome,nkc_reactome,
                 neut_reactome, mp_reactome, treg_reactome, hsc_reactome, wbm_reactome)
reactome <- reactome %>% reduce(full_join, by='pathway')
reactome <- reactome %>% remove_rownames %>% column_to_rownames(var="pathway")

reactome$sum.na <- rowSums(is.na(reactome))
hist(reactome$sum.na)
reactome <- reactome[reactome$sum.na < 5,]
reactome$sum.na <- NULL
reactome[is.na(reactome)] <- 0

reactome_heatmap <- as.matrix(reactome)
ht_opt$TITLE_PADDING = unit(c(3.5, 3.5), "points")
reactome_plot <-Heatmap(reactome_heatmap, col=col_fun, row_names_side= "right",row_dend_side= "left", name="NES", 
                    row_names_gp = gpar(fontsize = 9), #column_names_max_height = max_text_height(colnames(pathways)),
                    row_title = "REACTOME", row_title_gp = gpar(fontsize = 9, fill = "#ff9900", col = "white", border = NA),
                    rect_gp = gpar(col = "white", lwd = 2),
                    width = unit(3.96, "cm"), height = unit(4.95, "cm"),column_names_gp = gpar(fontsize = 9), show_column_names = FALSE, 
                    show_column_dend = FALSE, #row_dend_width = unit(0.2, "cm") 
                    show_row_dend = FALSE,#row_title = "Multiple Myeloma Cells",
                    column_order = colnames(reactome),
                    show_heatmap_legend = FALSE)
reactome_plot

tiff("results/reactome_plot_2.tiff", units="in", width=8, height=8, res=100)
reactome_plot
dev.off()


reactome <- reactome[order(row.names(reactome)), ]


reactome_heatmap <- as.matrix(reactome)
reactome_heatmap <- t(reactome_heatmap)

ht_opt$TITLE_PADDING = unit(c(3.5, 3.5), "points")
reactome_plot <- Heatmap(reactome_heatmap, col=col_fun, row_names_side= "left",row_dend_side= "left", name="NES", 
        row_names_gp = gpar(fontsize = 8), #column_names_max_height = max_text_height(colnames(pathways)),
        column_title = "REACTOME", column_title_gp = gpar(fill =  "#ff9900", col = "white", border = NA),
        rect_gp = gpar(col = "white", lwd = 2),
        width = unit(5.76, "cm"), height = unit(5, "cm"),column_names_gp = gpar(fontsize = 8), 
        show_column_dend = FALSE, #row_dend_width = unit(0.2, "cm") 
        show_row_dend = FALSE,#row_title = "Multiple Myeloma Cells",
        row_order = colnames(reactome),
        row_title_gp = gpar(fontsize = 8, fontface = "bold"), 
        show_heatmap_legend = FALSE, show_row_names=FALSE)
reactome_plot

tiff("results/reactome_plot.tiff", units="in", width=8, height=8, res=100)
reactome_plot
dev.off()

#_________________________________________________________________________________________________
#_________________________________________________________________________________________________

# heatmap of wikiPathways pathways
pc_wikiPathways <- pc_wikiPathways[c(1,6)]
colnames(pc_wikiPathways)[2] <- "Plasma Cells"
pc_wikiPathways$pathway <- gsub('WP_', '', pc_wikiPathways$pathway)
pc_wikiPathways$pathway <- gsub('_', ' ', pc_wikiPathways$pathway)

bmsc_wikiPathways <- bmsc_wikiPathways[c(1,6)]
colnames(bmsc_wikiPathways)[2] <- "Bone Marrow Stromal Cells"
bmsc_wikiPathways$pathway <- gsub('WP_', '', bmsc_wikiPathways$pathway)
bmsc_wikiPathways$pathway <- gsub('_', ' ', bmsc_wikiPathways$pathway)

ec_wikiPathways <- ec_wikiPathways[c(1,6)]
colnames(ec_wikiPathways)[2] <- "Endothelial Cells"
ec_wikiPathways$pathway <- gsub('WP_', '', ec_wikiPathways$pathway)
ec_wikiPathways$pathway <- gsub('_', ' ', ec_wikiPathways$pathway)

ocy_wikiPathways <- ocy_wikiPathways[c(1,6)]
colnames(ocy_wikiPathways)[2] <- "Osteocytes"
ocy_wikiPathways$pathway <- gsub('WP_', '', ocy_wikiPathways$pathway)
ocy_wikiPathways$pathway <- gsub('_', ' ', ocy_wikiPathways$pathway)

opc_wikiPathways <- opc_wikiPathways[c(1,6)]
colnames(opc_wikiPathways)[2] <- "Osteogenic Precursor Cells"
opc_wikiPathways$pathway <- gsub('WP_', '', opc_wikiPathways$pathway)
opc_wikiPathways$pathway <- gsub('_', ' ', opc_wikiPathways$pathway)

adcy_wikiPathways <- adcy_wikiPathways[c(1,6)]
colnames(adcy_wikiPathways)[2] <- "Adipocytes"
adcy_wikiPathways$pathway <- gsub('WP_', '', adcy_wikiPathways$pathway)
adcy_wikiPathways$pathway <- gsub('_', ' ', adcy_wikiPathways$pathway)

nkc_wikiPathways <- nkc_wikiPathways[c(1,6)]
colnames(nkc_wikiPathways)[2] <- "Natural Killer Cells"
nkc_wikiPathways$pathway <- gsub('WP_', '', nkc_wikiPathways$pathway)
nkc_wikiPathways$pathway <- gsub('_', ' ', nkc_wikiPathways$pathway)

neut_wikiPathways <- neut_wikiPathways[c(1,6)]
colnames(neut_wikiPathways)[2] <- "Neutrophils"
neut_wikiPathways$pathway <- gsub('WP_', '', neut_wikiPathways$pathway)
neut_wikiPathways$pathway <- gsub('_', ' ', neut_wikiPathways$pathway)

mp_wikiPathways <- mp_wikiPathways[c(1,6)]
colnames(mp_wikiPathways)[2] <- "Macrophages"
mp_wikiPathways$pathway <- gsub('WP_', '', mp_wikiPathways$pathway)
mp_wikiPathways$pathway <- gsub('_', ' ', mp_wikiPathways$pathway)

treg_wikiPathways <- treg_wikiPathways[c(1,6)]
colnames(treg_wikiPathways)[2] <- "Regulatory T Cells"
treg_wikiPathways$pathway <- gsub('WP_', '', treg_wikiPathways$pathway)
treg_wikiPathways$pathway <- gsub('_', ' ', treg_wikiPathways$pathway)

hsc_wikiPathways <- hsc_wikiPathways[c(1,6)]
colnames(hsc_wikiPathways)[2] <- "Hematopoietic Stem Cells"
hsc_wikiPathways$pathway <- gsub('WP_', '', hsc_wikiPathways$pathway)
hsc_wikiPathways$pathway <- gsub('_', ' ', hsc_wikiPathways$pathway)

wbm_wikiPathways <- wbm_wikiPathways[c(1,6)]
colnames(wbm_wikiPathways)[2] <- "Whole Bone Marrow"
wbm_wikiPathways$pathway <- gsub('WP_', '', wbm_wikiPathways$pathway)
wbm_wikiPathways$pathway <- gsub('_', ' ', wbm_wikiPathways$pathway)


wikiPathways <- list(pc_wikiPathways, bmsc_wikiPathways, ec_wikiPathways, ocy_wikiPathways, opc_wikiPathways, adcy_wikiPathways,nkc_wikiPathways,
                     neut_wikiPathways, mp_wikiPathways, treg_wikiPathways, hsc_wikiPathways, wbm_wikiPathways)
wikiPathways <- wikiPathways %>% reduce(full_join, by='pathway')
wikiPathways <- wikiPathways %>% remove_rownames %>% column_to_rownames(var="pathway")

wikiPathways$sum.na <- rowSums(is.na(wikiPathways))
hist(wikiPathways$sum.na)
wikiPathways <- wikiPathways[wikiPathways$sum.na < 8,]
wikiPathways$sum.na <- NULL
wikiPathways[is.na(wikiPathways)] <- 0

wikiPathways_heatmap <- as.matrix(wikiPathways)
ht_opt$TITLE_PADDING = unit(c(3.5, 3.5), "points")
wikiPathways_plot <-Heatmap(wikiPathways_heatmap, col=col_fun, row_names_side= "right",row_dend_side= "left", name="NES", 
                        row_names_gp = gpar(fontsize = 9), #column_names_max_height = max_text_height(colnames(pathways)),
                        row_title = "WikiPathways", row_title_gp = gpar(fontsize = 9, fill = "#00cccc", col = "white", border = NA),
                        rect_gp = gpar(col = "white", lwd = 2),
                        width = unit(3.96, "cm"), height = unit(5.28, "cm"),column_names_gp = gpar(fontsize = 9), show_column_names = FALSE, 
                        show_column_dend = FALSE, #row_dend_width = unit(0.2, "cm") 
                        show_row_dend = FALSE,#row_title = "Multiple Myeloma Cells",
                        column_order = colnames(wikiPathways),
                        show_heatmap_legend = TRUE, 
                        heatmap_legend_param= list(title = "NES",direction = "horizontal", 
                                                    legend_width = unit(3, "cm"), title_position = "topcenter",
                                                   title_gp = gpar(fontsize = 9)
                                                   ,at = c(-4, 0, 4)))
wikiPathways_plot

tiff("results/wikiPathways_plot_2.tiff", units="in", width=8, height=8, res=100)
wikiPathways_plot
dev.off()

wikiPathways <- wikiPathways[order(row.names(wikiPathways)), ]


wikiPathways_heatmap <- as.matrix(wikiPathways)
wikiPathways_heatmap <- t(wikiPathways_heatmap)

ht_opt$TITLE_PADDING = unit(c(3.5, 3.5), "points")
wikiPathways_plot <- Heatmap(wikiPathways_heatmap, col=col_fun,row_names_side= "left",row_dend_side= "left", name="NES", 
        row_names_gp = gpar(fontsize = 8), #column_names_max_height = max_text_height(colnames(pathways)),
        column_title = "WikiPathways", column_title_gp = gpar(fill =  "#00cccc", col = "white", border = NA),
        rect_gp = gpar(col = "white", lwd = 2),
        width = unit(5.76, "cm"), height = unit(5, "cm"),column_names_gp = gpar(fontsize = 8), 
        show_column_dend = FALSE, #row_dend_width = unit(0.2, "cm") 
        show_row_dend = FALSE,#row_title = "Multiple Myeloma Cells",
        row_order = colnames(wikiPathways),
        row_title_gp = gpar(fontsize = 8, fontface = "bold"),
        show_row_names=FALSE)
wikiPathways_plot

tiff("results/wikiPathways_plot.tiff", units="in", width=8, height=8, res=100)
wikiPathways_plot
dev.off()
#_________________________________________________________________________________________________
#_________________________________________________________________________________________________


# export figures
tiff("results/heatmap_plot_2.tiff", units="in", width=10, height=12, res=300)
draw(wikiPathways_plot %v% reactome_plot %v% kegg_plot, heatmap_legend_side = "top")
dev.off()


tiff("results/heatmap_plot.tiff", units="in", width=9, height=8, res=300)
kegg_plot+reactome_plot+wikiPathways_plot
dev.off()

#_________________________________________________________________________________________________
#_________________________________________________________________________________________________
# test enriched pathways in human pc

## load data
pc_data <- read.csv("data/pc/PC_MA.csv")
pc_data <- pc_data[,c(2,3,12,15)]
pc_data <- transform(pc_data, fisherPvalMin = pmin(fisherPvalUp, fisherPvalDown))
pc_data$score <- with(pc_data, ifelse(effectSize < 0, -log10(fisherPvalMin)*-1, -log10(fisherPvalMin)*1))
pc_data <- pc_data[,c(1,6)]
pc_data <- pc_data[complete.cases(pc_data), ]
pc_ranks <- deframe(pc_data)
head(pc_ranks, 20)


### test kegg enrichment
pathways_kegg <- gmtPathways("data/c2.cp.kegg.v2022.1.Hs.symbols.gmt")

pc_kegg <- fgsea(pathways=pathways_kegg, stats=pc_ranks)

pc_kegg <- pc_kegg[pc_kegg$padj < 0.001,]

pc_kegg <- pc_kegg %>%
  as_tibble() %>%
  arrange(desc(NES))

pc_kegg$pathway <- gsub('KEGG_', '', pc_kegg$pathway)
pc_kegg$pathway <- gsub('_', ' ', pc_kegg$pathway)

### test wikiPathways enrichment
pathways_wikiPathways <- gmtPathways("data/c2.cp.wikipathways.v2022.1.Hs.symbols.gmt")

pc_wikiPathways <- fgsea(pathways=pathways_wikiPathways, stats=pc_ranks)

pc_wikiPathways <- pc_wikiPathways[pc_wikiPathways$padj < 0.001,]

pc_wikiPathways <- pc_wikiPathways %>%
  as_tibble() %>%
  arrange(desc(NES))
pc_wikiPathways$pathway <- gsub('WP_', '', pc_wikiPathways$pathway)
pc_wikiPathways$pathway <- gsub('_', ' ', pc_wikiPathways$pathway)

### test reactome enrichment
pathways_reactome <- gmtPathways("data/c2.cp.reactome.v2022.1.Hs.symbols.gmt")

pc_reactome <- fgsea(pathways=pathways_reactome, stats=pc_ranks)

pc_reactome <- pc_reactome[pc_reactome$padj < 0.001,]

pc_reactome <- pc_reactome %>%
  as_tibble() %>%
  arrange(desc(NES))
pc_reactome$pathway <- gsub('REACTOME_', '', pc_reactome$pathway)
pc_reactome$pathway <- gsub('_', ' ', pc_reactome$pathway)

pc_pathways <- rbind(pc_kegg, pc_reactome, pc_wikiPathways)

write_csv(pc_pathways,"results/PC_Pathways_0.001.csv")

#_________________________________________________________________________________________________
#_________________________________________________________________________________________________

# test enriched pathways in mouse pc

## load data
mus_pc_data <- read.csv("data/pc/mPC_MA.csv")
mus_pc_data <- mus_pc_data[,c(1,4,13,16)]
mus_pc_data <- transform(mus_pc_data, fisherPvalMin = pmin(fisherPvalUp, fisherPvalDown))
mus_pc_data$score <- with(mus_pc_data, ifelse(effectSize < 0, -log10(fisherPvalMin)*-1, -log10(fisherPvalMin)*1))
mus_pc_data <- mus_pc_data[,c(1,6)]
mus_pc_data <- mus_pc_data[complete.cases(mus_pc_data), ]
mus_pc_ranks <- deframe(mus_pc_data)
head(mus_pc_ranks, 20)

### test kegg enrichment
pathways_kegg <- gmtPathways("data/c2.cp.kegg.v2022.1.Hs.symbols.gmt")

mus_pc_kegg <- fgsea(pathways=pathways_kegg, stats=mus_pc_ranks)

mus_pc_kegg <- mus_pc_kegg[mus_pc_kegg$padj < 0.001,]

mus_pc_kegg <- mus_pc_kegg %>%
  as_tibble() %>%
  arrange(desc(NES))
mus_pc_kegg$pathway <- gsub('KEGG_', '', mus_pc_kegg$pathway)
mus_pc_kegg$pathway <- gsub('_', ' ', mus_pc_kegg$pathway)

### test wikiPathways enrichment
pathways_wikiPathways <- gmtPathways("data/c2.cp.wikipathways.v2022.1.Hs.symbols.gmt")

mus_pc_wikiPathways <- fgsea(pathways=pathways_wikiPathways, stats=mus_pc_ranks)

mus_pc_wikiPathways <- mus_pc_wikiPathways[mus_pc_wikiPathways$padj < 0.001,]

mus_pc_wikiPathways <- mus_pc_wikiPathways %>%
  as_tibble() %>%
  arrange(desc(NES))
mus_pc_wikiPathways$pathway <- gsub('WP_', '', mus_pc_wikiPathways$pathway)
mus_pc_wikiPathways$pathway <- gsub('_', ' ', mus_pc_wikiPathways$pathway)

### test reactome enrichment
pathways_reactome <- gmtPathways("data/c2.cp.reactome.v2022.1.Hs.symbols.gmt")

mus_pc_reactome <- fgsea(pathways=pathways_reactome, stats=mus_pc_ranks)

mus_pc_reactome <- mus_pc_reactome[mus_pc_reactome$padj < 0.001,]

mus_pc_reactome <- mus_pc_reactome %>%
  as_tibble() %>%
  arrange(desc(NES))
mus_pc_reactome$pathway <- gsub('REACTOME_', '', mus_pc_reactome$pathway)
mus_pc_reactome$pathway <- gsub('_', ' ', mus_pc_reactome$pathway)

mus_pc_pathways <- rbind(mus_pc_kegg, mus_pc_reactome, mus_pc_wikiPathways)

write_csv(mus_pc_pathways,"results/MUS_PC_Pathways_0.001.csv")

#_________________________________________________________________________________________________

# heatmap of hs_mus_pc_kegg pathways
hs_pc_kegg <- pc_kegg[c(1,6)]
colnames(hs_pc_kegg)[2] <- "Human Plasma Cells"

mus_pc_kegg <- mus_pc_kegg[c(1,6)]
colnames(mus_pc_kegg)[2] <- "Mouse Plasma Cells"

hs_mus_pc_kegg <- list(hs_pc_kegg, mus_pc_kegg)

hs_mus_pc_kegg <- hs_mus_pc_kegg %>% reduce(full_join, by='pathway')
hs_mus_pc_kegg <- hs_mus_pc_kegg %>% remove_rownames %>% column_to_rownames(var="pathway")

hs_mus_pc_kegg$sum.na <- rowSums(is.na(hs_mus_pc_kegg))
hist(hs_mus_pc_kegg$sum.na)
hs_mus_pc_kegg <- hs_mus_pc_kegg[hs_mus_pc_kegg$sum.na < 1,]
hs_mus_pc_kegg$sum.na <- NULL
hs_mus_pc_kegg[is.na(hs_mus_pc_kegg)] <- 0

hs_mus_pc_kegg_heatmap <- as.matrix(hs_mus_pc_kegg)
col_fun = colorRamp2(c(-4, 0, 4), c("blue", "#e6e6e6", "red"))
ht_opt$TITLE_PADDING = unit(c(3.5, 3.5), "points")
hs_mus_pc_kegg_plot <-Heatmap(hs_mus_pc_kegg_heatmap, col=col_fun, row_names_side= "right",row_dend_side= "left", name="NES", 
                              row_names_gp = gpar(fontsize = 9), #column_names_max_height = max_text_height(colnames(pathways)),
                              row_title = "Kegg", row_title_gp = gpar(fontsize = 9, fill = "#33cc33", col = "white", border = NA),
                              rect_gp = gpar(col = "white", lwd = 2),
                              width = unit(0.66, "cm"), height = unit(6.27, "cm"),column_names_gp = gpar(fontsize = 9), 
                              show_column_dend = FALSE, #row_dend_width = unit(0.2, "cm") 
                              show_row_dend = FALSE,#row_title = "Multiple Myeloma Cells",
                              column_order = colnames(hs_mus_pc_kegg),
                              show_heatmap_legend = FALSE)
hs_mus_pc_kegg_plot

tiff("results/hs_mus_pc_kegg_plot_2.tiff", units="in", width=8, height=8, res=100)
hs_mus_pc_kegg_plot
dev.off()


hs_mus_pc_kegg <- hs_mus_pc_kegg[order(row.names(hs_mus_pc_kegg)), ]


hs_mus_pc_kegg_heatmap <- as.matrix(hs_mus_pc_kegg)
hs_mus_pc_kegg_heatmap <- t(hs_mus_pc_kegg_heatmap)

col_fun = colorRamp2(c(-4, 0, 4), c("blue", "#e6e6e6", "red"))
ht_opt$TITLE_PADDING = unit(c(1, 1), "points")
hs_mus_pc_kegg_plot <-Heatmap(hs_mus_pc_kegg_heatmap, col=col_fun, row_names_side= "left",row_dend_side= "left", name="NES", 
                              row_names_gp = gpar(fontsize = 7), #column_names_max_height = max_text_height(colnames(pathways)),
                              column_title = "KEGG", column_title_gp = gpar(fill = "#33cc33", col = "white", border = NA),
                              rect_gp = gpar(col = "white", lwd = 2),
                              width = unit(10.56, "cm"), height = unit(0.66, "cm"),column_names_gp = gpar(fontsize = 5), 
                              show_column_dend = FALSE, #row_dend_width = unit(0.2, "cm") 
                              show_row_dend = FALSE,#row_title = "Multiple Myeloma Cells",
                              row_order = colnames(hs_mus_pc_kegg),
                              row_title_gp = gpar(fontsize = 8),
                              show_heatmap_legend = FALSE)
hs_mus_pc_kegg_plot

tiff("results/hs_mus_pc_kegg_plot.tiff", units="in", width=8, height=8, res=100)
hs_mus_pc_kegg_plot
dev.off()

#_________________________________________________________________________________________________
#_________________________________________________________________________________________________

# heatmap of hs_mus_pc_reactome pathways
hs_pc_reactome <- pc_reactome[c(1,6)]
colnames(hs_pc_reactome)[2] <- "Human Plasma Cells"

mus_pc_reactome <- mus_pc_reactome[c(1,6)]
colnames(mus_pc_reactome)[2] <- "Mouse Plasma Cells"


hs_mus_pc_reactome <- list(hs_pc_reactome, mus_pc_reactome)

hs_mus_pc_reactome <- hs_mus_pc_reactome %>% reduce(full_join, by='pathway')
hs_mus_pc_reactome <- hs_mus_pc_reactome %>% remove_rownames %>% column_to_rownames(var="pathway")

hs_mus_pc_reactome$sum.na <- rowSums(is.na(hs_mus_pc_reactome))
hist(hs_mus_pc_reactome$sum.na)
hs_mus_pc_reactome <- hs_mus_pc_reactome[hs_mus_pc_reactome$sum.na < 1,]
hs_mus_pc_reactome$sum.na <- NULL
hs_mus_pc_reactome[is.na(hs_mus_pc_reactome)] <- 0

hs_mus_pc_reactome_heatmap <- as.matrix(hs_mus_pc_reactome)
ht_opt$TITLE_PADDING = unit(c(3.5, 3.5), "points")
hs_mus_pc_reactome_plot <-Heatmap(hs_mus_pc_reactome_heatmap, col=col_fun, row_names_side= "right",row_dend_side= "left", name="NES", 
                                  row_names_gp = gpar(fontsize = 9), #column_names_max_height = max_text_height(colnames(pathways)),
                                  row_title = "Reactome", row_title_gp = gpar(fontsize = 9, fill = "#ff9900", col = "white", border = NA),
                                  rect_gp = gpar(col = "white", lwd = 2),
                                  width = unit(0.66, "cm"), height = unit(33.33, "cm"),column_names_gp = gpar(fontsize = 9), show_column_names = FALSE, 
                                  show_column_dend = FALSE, #row_dend_width = unit(0.2, "cm") 
                                  show_row_dend = FALSE,#row_title = "Multiple Myeloma Cells",
                                  column_order = colnames(hs_mus_pc_reactome),
                                  show_heatmap_legend = FALSE)
hs_mus_pc_reactome_plot

tiff("results/hs_mus_pc_reactome_plot_2.tiff", units="in", width=8, height=15, res=100)
hs_mus_pc_reactome_plot
dev.off()

hs_mus_pc_reactome <- hs_mus_pc_reactome[order(row.names(hs_mus_pc_reactome)), ]


hs_mus_pc_reactome_heatmap <- as.matrix(hs_mus_pc_reactome)
hs_mus_pc_reactome_heatmap <- t(hs_mus_pc_reactome_heatmap)

ht_opt$TITLE_PADDING = unit(c(3.5, 3.5), "points")
hs_mus_pc_reactome_plot <- Heatmap(hs_mus_pc_reactome_heatmap, col=col_fun, row_names_side= "left",row_dend_side= "left", name="NES", 
                                   row_names_gp = gpar(fontsize = 7), #column_names_max_height = max_text_height(colnames(pathways)),
                                   column_title = "Reactome", column_title_gp = gpar(fill =  "#ff9900", col = "white", border = NA),
                                   rect_gp = gpar(col = "white", lwd = 2),
                                   width = unit(58.74, "cm"), height = unit(0.66, "cm"),column_names_gp = gpar(fontsize = 5), 
                                   show_column_dend = FALSE, #row_dend_width = unit(0.2, "cm") 
                                   show_row_dend = FALSE,#row_title = "Multiple Myeloma Cells",
                                   row_order = colnames(hs_mus_pc_reactome),
                                   row_title_gp = gpar(fontsize = 8, fontface = "bold"), 
                                   show_heatmap_legend = FALSE, show_row_names=FALSE)
hs_mus_pc_reactome_plot

tiff("results/hs_mus_pc_reactome_plot.tiff", units="in", width=30, height=8, res=100)
hs_mus_pc_reactome_plot
dev.off()

#_________________________________________________________________________________________________
#_________________________________________________________________________________________________

# heatmap of hs_mus_pc_wikiPathways pathways
hs_pc_wikiPathways <- pc_wikiPathways[c(1,6)]
colnames(hs_pc_wikiPathways)[2] <- "Human Plasma Cells"

mus_pc_wikiPathways <- mus_pc_wikiPathways[c(1,6)]
colnames(mus_pc_wikiPathways)[2] <- "Mouse Plasma Cells"


hs_mus_pc_wikiPathways <- list(hs_pc_wikiPathways, mus_pc_wikiPathways)
hs_mus_pc_wikiPathways <- hs_mus_pc_wikiPathways %>% reduce(full_join, by='pathway')
hs_mus_pc_wikiPathways <- hs_mus_pc_wikiPathways %>% remove_rownames %>% column_to_rownames(var="pathway")

hs_mus_pc_wikiPathways$sum.na <- rowSums(is.na(hs_mus_pc_wikiPathways))
hist(hs_mus_pc_wikiPathways$sum.na)
hs_mus_pc_wikiPathways <- hs_mus_pc_wikiPathways[hs_mus_pc_wikiPathways$sum.na < 1,]
hs_mus_pc_wikiPathways$sum.na <- NULL
hs_mus_pc_wikiPathways[is.na(hs_mus_pc_wikiPathways)] <- 0

hs_mus_pc_wikiPathways_heatmap <- as.matrix(hs_mus_pc_wikiPathways)
ht_opt$TITLE_PADDING = unit(c(3.5, 3.5), "points")
hs_mus_pc_wikiPathways_plot <-Heatmap(hs_mus_pc_wikiPathways_heatmap, col=col_fun, row_names_side= "right",row_dend_side= "left", name="NES", 
                                      row_names_gp = gpar(fontsize = 9), #column_names_max_height = max_text_height(colnames(pathways)),
                                      row_title = "WikiPathways", row_title_gp = gpar(fontsize = 9, fill = "#00cccc", col = "white", border = NA),
                                      rect_gp = gpar(col = "white", lwd = 2),
                                      width = unit(0.66, "cm"), height = unit(3.63, "cm"),column_names_gp = gpar(fontsize = 9), show_column_names = FALSE, 
                                      show_column_dend = FALSE, #row_dend_width = unit(0.2, "cm") 
                                      show_row_dend = FALSE,#row_title = "Multiple Myeloma Cells",
                                      column_order = colnames(hs_mus_pc_wikiPathways),
                                      show_heatmap_legend = TRUE, 
                                      heatmap_legend_param= list(title = "NES",direction = "horizontal", 
                                                                 legend_width = unit(3, "cm"), title_position = "topcenter",
                                                                 title_gp = gpar(fontsize = 9)
                                                                 ,at = c(-4, 0, 4)))
hs_mus_pc_wikiPathways_plot

tiff("results/hs_mus_pc_wikiPathways_plot_2.tiff", units="in", width=8, height=15, res=100)
hs_mus_pc_wikiPathways_plot
dev.off()

hs_mus_pc_wikiPathways <- hs_mus_pc_wikiPathways[order(row.names(hs_mus_pc_wikiPathways)), ]


hs_mus_pc_wikiPathways_heatmap <- as.matrix(hs_mus_pc_wikiPathways)
hs_mus_pc_wikiPathways_heatmap <- t(hs_mus_pc_wikiPathways_heatmap)

ht_opt$TITLE_PADDING = unit(c(3.5, 3.5), "points")
hs_mus_pc_wikiPathways_plot <- Heatmap(hs_mus_pc_wikiPathways_heatmap, col=col_fun,row_names_side= "left",row_dend_side= "left", name="NES", 
                                       row_names_gp = gpar(fontsize = 7), #column_names_max_height = max_text_height(colnames(pathways)),
                                       column_title = "WikiPathways", column_title_gp = gpar(fill =  "#00cccc", col = "white", border = NA),
                                       rect_gp = gpar(col = "white", lwd = 2),
                                       width = unit(11.88, "cm"), height = unit(0.66, "cm"),column_names_gp = gpar(fontsize = 5), 
                                       show_column_dend = FALSE, #row_dend_width = unit(0.2, "cm") 
                                       show_row_dend = FALSE,#row_title = "Multiple Myeloma Cells",
                                       row_order = colnames(hs_mus_pc_wikiPathways),
                                       row_title_gp = gpar(fontsize = 8, fontface = "bold"),
                                       show_row_names=FALSE)
hs_mus_pc_wikiPathways_plot

tiff("results/hs_mus_pc_wikiPathways_plot.tiff", units="in", width=8, height=8, res=100)
hs_mus_pc_wikiPathways_plot
dev.off()
#_________________________________________________________________________________________________
#_________________________________________________________________________________________________

# export figures

tiff("results/hs_mus_pc_heatmap_plot_2.tiff", units="in", width=18, height=20, res=300)
draw(hs_mus_pc_wikiPathways_plot %v% hs_mus_pc_reactome_plot %v% hs_mus_pc_kegg_plot, heatmap_legend_side = "top")
dev.off()

tiff("results/hs_mus_pc_heatmap_plot.tiff", units="in", width=40, height=10, res=300)
hs_mus_pc_kegg_plot+hs_mus_pc_reactome_plot+hs_mus_pc_wikiPathways_plot
dev.off()

#_________________________________________________________________________________________________
#_________________________________________________________________________________________________


         