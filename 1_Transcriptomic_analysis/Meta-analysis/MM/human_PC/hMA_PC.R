# set working directory
setwd("~/Work/PhD/Projects/MoD TrD LRI MM-Pain/Computational Analysis/Ligand_Receptor/Data/MM/PC/hDatasets/MA")

# load packages
library(MetaIntegrator)
library(org.Hs.eg.db)
library(annotate)
library(tidyverse)

# prepare GSE6477 
## read file 
GSE6477 <- read.csv("data/GSE6477 Normalized.Filtered(ID).Annotated.Data.csv")

## prepare datasetObject
expr <- data.matrix(GSE6477[,c(4:119)])
rownames(expr) <- GSE6477$PROBEID

keys <- as.character(GSE6477$SYMBOL)
names(keys) <- GSE6477$PROBEID

class <- c(1,1,1,1,	1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0)
names(class) <- colnames(GSE6477)[-c(1:3)]

pheno <- data.frame(c("MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM", "MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM", "MM",	"MM",	"MM",	"MM",	"MM", "MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",
                      "MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM", "MM", "MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"HD" ,"MM", "MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",	
                      "MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"HD",	"MM",	"MM",	"MM",	"MM",	"MM","MM",	"MM",	"MM",	"MM",	"MM",	"HD",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",
                      "MM",	"MM",	"HD",	"HD",	"HD",	"HD",	"HD",	"HD",	"HD",	"HD",	"HD",	"HD",	"HD"	,"HD"))
rownames(pheno) <- colnames(GSE6477)[-c(1:3)]

GSE6477_datasetObject <- list(expr = expr, 
                               keys = keys,
                               class = class, pheno = pheno, formattedName = "GSE6477")

## remove objects 
rm(GSE6477,expr,pheno,class,keys) 
 
## check for completeness
checkDataObject(GSE6477_datasetObject, "Dataset")

# prepare GSE39754 
## read file
GSE39754 <- read.csv("data/GSE39754 Normalized.Filtered(ID).Annotated.Data.csv")

## prepare datasetObject
expr <- data.matrix(GSE39754[,c(4:179)])
rownames(expr) <- GSE39754$PROBEID

keys <- as.character(GSE39754$SYMBOL)
names(keys) <- GSE39754$PROBEID

class <- c(0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
names(class) <- colnames(GSE39754)[-c(1:3)]

pheno <- data.frame(c("HD","HD","HD","HD","HD","HD","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM"))
rownames(pheno) <- colnames(GSE39754)[-c(1:3)]

GSE39754_datasetObject <- list(expr = expr, 
                               keys = keys,
                               class = class, pheno = pheno, formattedName = "GSE39754")

## remove objects 
rm(GSE39754,expr,pheno,class,keys) 

## check for completeness
checkDataObject(GSE39754_datasetObject, "Dataset")

# prepare GSE47552 
## read file 
GSE47552 <- read.csv("data/GSE47552 Normalized.Filtered(ID).Annotated.Data.csv")

## prepare datasetObject
expr <- data.matrix(GSE47552[,c(4:49)])
rownames(expr) <- GSE47552$PROBEID

keys <- as.character(GSE47552$SYMBOL)
names(keys) <- GSE47552$PROBEID

class <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0)
names(class) <- colnames(GSE47552)[-c(1:3)]

pheno <- data.frame(c("MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","HD","HD","HD","HD","HD"))
rownames(pheno) <- colnames(GSE47552)[-c(1:3)]

GSE47552_datasetObject <- list(expr = expr, 
                                         keys = keys,
                                         class = class, pheno = pheno, formattedName = "GSE47552")

## remove objects 
rm(GSE47552,expr,pheno,class,keys) 

## check for completeness
checkDataObject(GSE47552_datasetObject, "Dataset")

# prepare GSE76425 
## read file
GSE76425 <- read.csv("data/GSE76425 Filtered(Counts).Normalized.Annotated.Data_voom.csv")
GSE76425$gene_symbol <-getSYMBOL(as.character(GSE76425$id), data='org.Hs.eg')


## prepare datasetObject
expr <- data.matrix(GSE76425[,c(2:12)])
rownames(expr) <- GSE76425$id

keys <- as.character(GSE76425$gene_symbol)
names(keys) <- GSE76425$id

class <- c(0,0,0,0,1,1,1,1,1,1,1)
names(class) <- colnames(GSE76425)[-c(1,13)]

pheno <- data.frame(c("HD","HD","HD","HD","MM","MM","MM","MM","MM","MM","MM"))
rownames(pheno) <- colnames(GSE76425)[-c(1,13)]

GSE76425_datasetObject <- list(expr = expr, 
                               keys = keys,
                               class = class, pheno = pheno, formattedName = "GSE76425")

## remove objects 
rm(GSE76425,expr,pheno,class,keys) 

## check for completeness
checkDataObject(GSE76425_datasetObject, "Dataset")

# prepare GSE116294 
## read file
GSE116294  <- read.csv("data/GSE116294 Normalized.Filtered(ID).Annotated.Data.csv")

## prepare datasetObject
expr <- data.matrix(GSE116294[,c(4:57)])
rownames(expr) <- GSE116294$PROBEID

keys <- as.character(GSE116294$SYMBOL)
names(keys) <- GSE116294$PROBEID

class <- c(0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
names(class) <- colnames(GSE116294)[-c(1:3)]

pheno <- data.frame(c("HD","HD","HD","HD","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM"))
rownames(pheno) <- colnames(GSE116294)[-c(1:3)]

GSE116294_datasetObject <- list(expr = expr, 
                               keys = keys,
                               class = class, pheno = pheno, formattedName = "GSE116294")

## remove objects 
rm(GSE116294,expr,pheno,class,keys) 

## check for completeness
checkDataObject(GSE116294_datasetObject, "Dataset")

# prepare GSE125361 
## read file
GSE125361 <- read.csv("data/GSE125361 Normalized.Filtered(ID).Annotated.Data.csv")

## prepare datasetObject
expr <- data.matrix(GSE125361[,c(4:51)])
rownames(expr) <- GSE125361$LOCUSLINK_ID

keys <- as.character(GSE125361$GENE_SYMBOL)
names(keys) <- GSE125361$LOCUSLINK_ID

class <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,0)
names(class) <- colnames(GSE125361)[-c(1:3)]

pheno <- data.frame(c("MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","HD","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","HD","MM","MM","MM","MM","MM","MM","MM","HD"))
rownames(pheno) <- colnames(GSE125361)[-c(1:3)]

GSE125361_datasetObject <- list(expr = expr, 
                               keys = keys,
                               class = class, pheno = pheno, formattedName = "GSE125361")

## remove objects 
rm(GSE125361,expr,pheno,class,keys) 

## check for completeness
checkDataObject(GSE125361_datasetObject, "Dataset")


# prepare GSE153380
## read file
GSE153380 <- read.csv("data/GSE153380 Filtered(ID,Counts).Normalized.Annotated.Data_voom.csv")
GSE153380$gene_symbol <-getSYMBOL(as.character(GSE153380$entrezgene_id), data='org.Hs.eg')

## prepare datasetObject
expr <- data.matrix(GSE153380[,c(2:32)])
rownames(expr) <- GSE153380$entrezgene_id

keys <- as.character(GSE153380$gene_symbol)
names(keys) <- GSE153380$entrezgene_id

class <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0)
names(class) <- colnames(GSE153380)[-c(1,33)]

pheno <- data.frame(c("MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","HD","HD","HD"))
rownames(pheno) <- colnames(GSE153380)[-c(1,33)]

GSE153380_datasetObject <- list(expr = expr, 
                               keys = keys,
                               class = class, pheno = pheno, formattedName = "GSE153380")

## remove objects 
rm(GSE153380,expr,pheno,class,keys) 

## check for completeness
checkDataObject(GSE153380_datasetObject, "Dataset")

# prepare GSE175384
## read file
GSE175384 <- read.csv("data/GSE175384 Filtered(ID,Counts).Normalized.Annotated.Data_voom.csv")
GSE175384$gene_symbol <-getSYMBOL(as.character(GSE175384$entrezgene_id), data='org.Hs.eg')

## prepare datasetObject
expr <- data.matrix(GSE175384[,c(2:65)])
rownames(expr) <- GSE175384$entrezgene_id

keys <- as.character(GSE175384$gene_symbol)
names(keys) <- GSE175384$entrezgene_id

class <- c(1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
names(class) <- colnames(GSE175384)[-c(1,66)]

pheno <- data.frame(c("MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","HD","HD","HD","HD","HD","HD","HD","HD","HD","HD","HD","HD","HD","HD","HD","HD","HD","HD","HD","HD","HD","HD","HD","HD","HD","HD","HD","HD","HD","HD","HD","HD",
                      "MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM"))
rownames(pheno) <- colnames(GSE175384)[-c(1,66)]

GSE175384_datasetObject <- list(expr = expr, 
                                keys = keys,
                                class = class, pheno = pheno, formattedName = "GSE175384")

## remove objects 
rm(GSE175384,expr,pheno,class,keys) 

## check for completeness
checkDataObject(GSE175384_datasetObject, "Dataset")

# create metaObject
pc_datasets <- list(GSE6477_datasetObject, GSE39754_datasetObject, GSE47552_datasetObject,GSE76425_datasetObject,GSE116294_datasetObject,GSE125361_datasetObject,GSE153380_datasetObject,GSE175384_datasetObject)
names(pc_datasets) = c(GSE6477_datasetObject$formattedName, GSE39754_datasetObject$formattedName, GSE47552_datasetObject$formattedName,GSE76425_datasetObject$formattedName,GSE116294_datasetObject$formattedName,GSE125361_datasetObject$formattedName,GSE153380_datasetObject$formattedName,GSE175384_datasetObject$formattedName)
MetaObj=list() 
MetaObj$originalData <- pc_datasets

## check for completeness
checkDataObject(MetaObj, "Meta", "Pre-Analysis")


# coconutRes <- coconutMetaIntegrator( metaObject = MetaObj)


## run meta-analysis
MetaObj <- runMetaAnalysis(MetaObj, maxCores=1)

MetaObj <- filterGenes(MetaObj, isLeaveOneOut = FALSE #,FDRThresh = 0.05 
                       )

save(MetaObj, file = "data/MetaObj.RData")
load("data/MetaObj.RData")


## get resutls 
pc_ma <- MetaObj[["metaAnalysis"]][["pooledResults"]]

pc_ma <- transform(pc_ma, fisherFDRMin = pmin(fisherFDRUp, fisherFDRDown))

pc_ma$gene_symbol <- rownames(pc_ma)

sum(duplicated(pc_ma$gene_symbol))

pc_ma <- pc_ma %>%
  select(gene_symbol, everything())

write_csv(pc_ma,"results/hPC_MA.csv")

deg <- summarizeFilterResults(MetaObj, getMostRecentFilter(MetaObj))
deg_pos <- deg[["pos"]]
deg_neg <- deg[["neg"]]

## forest plots
# LGMN
forestPlot(MetaObj, "5641")
# ADM
forestPlot(MetaObj, "133")
# NMB
forestPlot(MetaObj, "4828")
# FGF7
forestPlot(MetaObj, "2252")
# VEGFB
forestPlot(MetaObj, "7423")
# CXCL12
forestPlot(MetaObj, "6387")
# NTN4
forestPlot(MetaObj, "59277")
# MIF
forestPlot(MetaObj, "4282")
# PGF
forestPlot(MetaObj, "5228")
# THBS1
forestPlot(MetaObj, "7057")
#EFEMP1
forestPlot(MetaObj, "2202")
#NLGN2
forestPlot(MetaObj, "57555")
#OLFM2
forestPlot(MetaObj, "93145")
#TFPI
forestPlot(MetaObj, "7035")
#GPC1
forestPlot(MetaObj, "2817")
#SNCA
forestPlot(MetaObj, "6622")
#GPNMB
forestPlot(MetaObj, "10457")
#CPNE3
forestPlot(MetaObj, "8895")
#SEMA6A
forestPlot(MetaObj, "57556")
#NMU
forestPlot(MetaObj, "10874")
#IGF1
forestPlot(MetaObj, "3479")
#DKK1
forestPlot(MetaObj, "22943")
#NCAM1
forestPlot(MetaObj, "4684")
#NRG3
forestPlot(MetaObj, "10718")
#WNT5A
forestPlot(MetaObj, "7474")
#HGF
forestPlot(MetaObj, "3082")
#NRG2
forestPlot(MetaObj, "9542")

#NGF
forestPlot(MetaObj, "NGF")

#TNFSF11
forestPlot(MetaObj, "TNFSF11")

#IL6
forestPlot(MetaObj, "IL6")

#IL1B
forestPlot(MetaObj, "IL1B")

#TNF
forestPlot(MetaObj, "TNF")

#CCL3
forestPlot(MetaObj, "CCL3")

#VEGFA
forestPlot(MetaObj, "VEGFA")



tiff("results/pc_ligands.tiff", units="in", width=14.5, height=8, res=300)

ggForestPlot(MetaObj, genes=c("IGF1", "MIF", "NRG2" ,"WNT5A"), facetCols=4)  +
  facet_grid(~factor(gene, levels=c("IGF1", "MIF", "NRG2" ,"WNT5A")),scales="free") +
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(color = "black", size=20)) 
  #theme(axis.text.x = element_text(color="#000000", size=20),
  #      axis.text.y = element_text(color="#000000", size=20))


dev.off()


tiff("results/pain_pc_ligands.tiff", units="in", width=15, height=8, res=300)

ggForestPlot(MetaObj, 
             genes=c("NGF", "IL6", "IL1B", "TNF", "VEGFA"), 
             facetCols=5)  +
  facet_grid(~factor(gene, 
                     levels=c("NGF", "IL6", "IL1B", "TNF", "VEGFA"))
             ,scales="free") +
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(color = "black", size=20)) 
#theme(axis.text.x = element_text(color="#000000", size=20),
#      axis.text.y = element_text(color="#000000", size=20))


dev.off()

#########################################

calc_overlap <- function(list1, list2) {
  length(intersect(list1, list2)) / min(length(list1), length(list2))
}

lists <- list('GSE116294'= GSE116294$SYMBOL, 'GSE125361'=GSE125361$GENE_SYMBOL, 
              'GSE153380'=GSE153380$gene_symbol, 'GSE175384'=GSE175384$gene_symbol, 
              'GSE39754'=GSE39754$SYMBOL, 'GSE47552'=GSE47552$SYMBOL, 
              'GSE6477'=GSE6477$SYMBOL, 'GSE76425'=GSE76425$gene_symbol)


lists_overlap <- sapply(lists, function(x) sapply(lists, function(y) calc_overlap(x,y)))

write.csv(lists_overlap, file="results/gene_lists_overlap.csv")


