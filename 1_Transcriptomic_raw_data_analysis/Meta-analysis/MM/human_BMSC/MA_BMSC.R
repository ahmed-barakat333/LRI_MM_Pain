## set working directory
setwd("~/Work/PhD/Projects/MoD TrD LRI MM-Pain/Computational Analysis/Ligand_Receptor/Data/MM/BMSC/hDatasets/MA")

## load packages
library(MetaIntegrator)
library(org.Hs.eg.db)
library(annotate)
library(tidyverse)


## prepare GSE36474 
# read file 
GSE36474 <- read.csv("data/GSE36474 Normalized.Filtered(ID).Annotated.Data.csv")

# prepare datasetObject
expr <- data.matrix(GSE36474[,c(4:10)])
rownames(expr) <- GSE36474$PROBEID

keys <- as.character(GSE36474$SYMBOL)
names(keys) <- GSE36474$PROBEID

class <- c(0,0,0,1,1,1,1)
names(class) <- colnames(GSE36474)[-c(1:3)]

pheno <- data.frame(c("HD","HD","HD","MM","MM","MM","MM"))
rownames(pheno) <- colnames(GSE36474)[-c(1:3)]

GSE36474_datasetObject <- list(expr = expr, 
                               keys = keys,
                               class = class, pheno = pheno, formattedName = "GSE36474")

# remove objects 
rm(GSE36474,expr,pheno,class,keys) 
 
# check for completeness
checkDataObject(GSE36474_datasetObject, "Dataset")

## prepare GSE46053_Coculture 
# read file
GSE46053_Coculture <- read.csv("data/GSE46053 Coculture Normalized.Filtered(ID).Annotated.Data.csv")

# prepare datasetObject
expr <- data.matrix(GSE46053_Coculture[,c(4:19)])
rownames(expr) <- GSE46053_Coculture$PROBEID

keys <- as.character(GSE46053_Coculture$SYMBOL)
names(keys) <- GSE46053_Coculture$PROBEID

class <- c(0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1)
names(class) <- colnames(GSE46053_Coculture)[-c(1:3)]

pheno <- data.frame(c("HD","HD","HD","HD","HD","HD","HD","HD","MM","MM","MM","MM","MM","MM","MM","MM"))
rownames(pheno) <- colnames(GSE46053_Coculture)[-c(1:3)]

GSE46053_Coculture_datasetObject <- list(expr = expr, 
                               keys = keys,
                               class = class, pheno = pheno, formattedName = "GSE46053_Coculture")

# remove objects 
rm(GSE46053_Coculture,expr,pheno,class,keys) 

# check for completeness
checkDataObject(GSE46053_Coculture_datasetObject, "Dataset")

## prepare GSE78235 
# read file
GSE78235 <- read.csv("data/GSE78235 Normalized.Filtered(ID).Annotated.Data.csv")

# prepare datasetObject
expr <- data.matrix(GSE78235[,c(4:9)])
rownames(expr) <- GSE78235$PROBEID

keys <- as.character(GSE78235$SYMBOL)
names(keys) <- GSE78235$PROBEID

class <- c(0,0,0,0,1,1)
names(class) <- colnames(GSE78235)[-c(1:3)]

pheno <- data.frame(c("HD","HD","HD","HD","MM","MM"))
rownames(pheno) <- colnames(GSE78235)[-c(1:3)]

GSE78235_datasetObject <- list(expr = expr, 
                               keys = keys,
                               class = class, pheno = pheno, formattedName = "GSE78235")

# remove objects 
rm(GSE78235,expr,pheno,class,keys) 

# check for completeness
checkDataObject(GSE78235_datasetObject, "Dataset")

## prepare GSE80608 
# read file
GSE80608  <- read.csv("data/GSE80608 Normalized.Filtered(ID).Annotated.Data.csv")

# prepare datasetObject
expr <- data.matrix(GSE80608[,c(4:23)])
rownames(expr) <- GSE80608$PROBEID

keys <- as.character(GSE80608$SYMBOL)
names(keys) <- GSE80608$PROBEID

class <- c(0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1)
names(class) <- colnames(GSE80608)[-c(1:3)]

pheno <- data.frame(c("HD","HD","HD","HD","HD","HD","HD","HD","HD","HD","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM"))
rownames(pheno) <- colnames(GSE80608)[-c(1:3)]

GSE80608_datasetObject <- list(expr = expr, 
                               keys = keys,
                               class = class, pheno = pheno, formattedName = "GSE80608")

# remove objects 
rm(GSE80608,expr,pheno,class,keys) 

# check for completeness
checkDataObject(GSE80608_datasetObject, "Dataset")

## prepare GSE87073 
# read file
GSE87073 <- read.csv("data/GSE87073 Normalized.Filtered(ID).Annotated.Data.csv")

# prepare datasetObject
expr <- data.matrix(GSE87073[,c(4:13)])
rownames(expr) <- GSE87073$PROBEID

keys <- as.character(GSE87073$SYMBOL)
names(keys) <- GSE87073$PROBEID

class <- c(0,0,0,0,0,1,1,1,1,1)
names(class) <- colnames(GSE87073)[-c(1:3)]

pheno <- data.frame(c("HD","HD","HD","HD","HD","MM","MM","MM","MM","MM"))
rownames(pheno) <- colnames(GSE87073)[-c(1:3)]

GSE87073_datasetObject <- list(expr = expr, 
                               keys = keys,
                               class = class, pheno = pheno, formattedName = "GSE87073")

# remove objects 
rm(GSE87073,expr,pheno,class,keys) 

# check for completeness
checkDataObject(GSE87073_datasetObject, "Dataset")


## create metaObject
bmsc_datasets <- list(GSE36474_datasetObject, GSE46053_Coculture_datasetObject, GSE78235_datasetObject,GSE80608_datasetObject,GSE87073_datasetObject)
names(bmsc_datasets) = c(GSE36474_datasetObject$formattedName, GSE46053_Coculture_datasetObject$formattedName, GSE78235_datasetObject$formattedName,GSE80608_datasetObject$formattedName,GSE87073_datasetObject$formattedName)
MetaObj=list() 
MetaObj$originalData <- bmsc_datasets

# check for completeness
checkDataObject(MetaObj, "Meta", "Pre-Analysis")

# run meta-analysis
MetaObj <- runMetaAnalysis(MetaObj, maxCores=1)

MetaObj <- filterGenes(MetaObj, isLeaveOneOut = FALSE#, FDRThresh = 0.05
                       )

save(MetaObj, file = "data/MetaObj.RData")
load("data/MetaObj.RData")

# get results
bmsc_ma <- MetaObj[["metaAnalysis"]][["pooledResults"]]

bmsc_ma <- transform(bmsc_ma, fisherFDRMin = pmin(fisherFDRUp, fisherFDRDown))

bmsc_ma$gene_symbol <- rownames(bmsc_ma)

sum(duplicated(bmsc_ma$gene_symbol))
                      
bmsc_ma <- bmsc_ma %>%
  select(gene_symbol, everything())

write_csv(bmsc_ma,"results /BMSC_MA.csv")

deg <- summarizeFilterResults(MetaObj, getMostRecentFilter(MetaObj))
deg_pos <- deg[["pos"]]
deg_neg <- deg[["neg"]]

# forest plots
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
#Gbmsc1
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
forestPlot(MetaObj, "NRG2")

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




tiff("results /bmsc_ligands.tiff", units="in", width=11.9, height=8, res=300)

ggForestPlot(MetaObj, genes=c("THBS1","FGF7" ,"SEMA6A"), facetCols=3)  +
facet_grid(~factor(gene, levels=c("THBS1","FGF7" ,"SEMA6A")),scales="free") +
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(color = "black", size=20)) 
#theme(axis.text.x = element_text(color="#000000", size=20),
#      axis.text.y = element_text(color="#000000", size=20))


dev.off()

tiff("results /pain_bmsc_ligands.tiff", units="in", width=15, height=8, res=300)

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

lists <- list('GSE36474'=GSE36474$SYMBOL, 'GSE78235'=GSE78235$SYMBOL, 'GSE80608'=GSE80608$SYMBOL, 'GSE87073'=GSE87073$SYMBOL, 'GSE46053_Coculture'=GSE46053_Coculture$SYMBOL)


lists_overlap <- sapply(lists, function(x) sapply(lists, function(y) calc_overlap(x,y)))

write.csv(lists_overlap, file="gene_lists_overlap.csv")
