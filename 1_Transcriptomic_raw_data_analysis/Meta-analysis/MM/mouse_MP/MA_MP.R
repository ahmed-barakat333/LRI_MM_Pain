# set working directory
setwd("~/Work/PhD/Projects/MM Pain/Computational Analysis/Ligand_Receptor/Data/MM/MP/mDatasets/MA")

# load packages
library(MetaIntegrator)
library(org.Mm.eg.db)
library(annotate)
library(tidyverse)

# prepare GSE143025 
## read file 
GSE143025 <- read.csv("data/GSE143025_Normalized.Filtered(ID).Annotated.Data.csv")

## prepare datasetObject
expr <- data.matrix(GSE143025[,c(4:9)])
rownames(expr) <- GSE143025$PROBEID

keys <- as.character(GSE143025$SYMBOL)
names(keys) <- GSE143025$PROBEID

class <- c(0,0,0,1,1,1)
names(class) <- colnames(GSE143025)[-c(1:3)]

pheno <- data.frame(c("HD","HD","HD","MM","MM","MM"))
rownames(pheno) <- colnames(GSE143025)[-c(1:3)]

GSE143025_datasetObject <- list(expr = expr, 
                               keys = keys,
                               class = class, pheno = pheno, formattedName = "GSE143025")

## remove objects 
rm(GSE143025,expr,pheno,class,keys) 
 
## check for completeness
checkDataObject(GSE143025_datasetObject, "Dataset")


# prepare GSE176385
## read file
GSE176385 <- read.csv("data/GSE176385_Filtered(ID,Counts).Normalized.Annotated.Data_voom.csv")
GSE176385$gene_symbol <-getSYMBOL(as.character(GSE176385$entrezgene_id), data='org.Mm.eg.db')

## prepare datasetObject
expr <- data.matrix(GSE176385[,c(2:5)])
rownames(expr) <- GSE176385$entrezgene_id

keys <- as.character(GSE176385$gene_symbol)
names(keys) <- GSE176385$entrezgene_id

class <- c(1,1,0,0)
names(class) <- colnames(GSE176385)[-c(1,6)]

pheno <- data.frame(c("MM","MM","HD","HD"))
rownames(pheno) <- colnames(GSE176385)[-c(1,6)]

GSE176385_datasetObject <- list(expr = expr, 
                               keys = keys,
                               class = class, pheno = pheno, formattedName = "GSE176385")

## remove objects 
rm(GSE176385,expr,pheno,class,keys) 

## check for completeness
checkDataObject(GSE176385_datasetObject, "Dataset")

# create metaObject
mp_datasets <- list(GSE143025_datasetObject,GSE176385_datasetObject)
names(mp_datasets) = c(GSE143025_datasetObject$formattedName, GSE176385_datasetObject$formattedName)
MetaObj=list() 
MetaObj$originalData <- mp_datasets

## check for completeness
checkDataObject(MetaObj, "Meta", "Pre-Analysis")


# coconutRes <- coconutMetaIntegrator( metaObject = MetaObj)


## run meta-analysis
MetaObj <- runMetaAnalysis(MetaObj, maxCores=1)

MetaObj <- filterGenes(MetaObj, isLeaveOneOut = FALSE, #FDRThresh = 0.1
                       )

## get resutls 
mp_ma_mus <- MetaObj[["metaAnalysis"]][["pooledResults"]]

mp_ma_mus$gene_symbol_mus <- rownames(mp_ma_mus)

sum(duplicated(mp_ma_mus$gene_symbol_mus))

mp_ma_mus <- transform(mp_ma_mus, fisherFDRMin = pmin(fisherFDRUp, fisherFDRDown))

mp_ma_mus <- mp_ma_mus %>%
            dplyr::select(gene_symbol_mus, everything())

# ortholog mapping from mouse to human
musGenes <- mp_ma_mus$gene_symbol_mus

convertMouseGeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl",host="https://dec2021.archive.ensembl.org")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl",host="https://dec2021.archive.ensembl.org")
  
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  humanx <- unique(genesV2[,c(1,2)])
  
  return(humanx)
}

gene_symbol_mapped.hs <- convertMouseGeneList(musGenes)

#mouse_human_genes = read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")

#convert_mouse_to_human <- function(gene_list){

#output = c()

#for(gene in gene_list){
#  class_key = (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name=="mouse, laboratory"))[['DB.Class.Key']]
#  if(!identical(class_key, integer(0)) ){
#    human_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="human"))[,"Symbol"]
#    for(human_gene in human_genes){
#      output = append(output,human_gene)
#    }
#  }
#}
#return (output)
#}

#gene_symbol_mapped.hs <- convert_mouse_to_human(musGenes)


colnames(gene_symbol_mapped.hs)[1] <- "gene_symbol_mus"
colnames(gene_symbol_mapped.hs)[2] <- "gene_symbol_hs.mapped"
mp_ma_mus <- merge(gene_symbol_mapped.hs,mp_ma_mus, by = "gene_symbol_mus")


# check for duplicate symbols
sum(duplicated(mp_ma_mus$gene_symbol_hs.mapped))
sum(duplicated(mp_ma_mus$gene_symbol_mus))
## There are many to many mapping

# summarize multiple MGI.symbol mapping to same HGNC.symbol 
mp_ma_mus <- mp_ma_mus %>% 
  group_by(gene_symbol_hs.mapped) %>% 
  mutate(gene_symbol_mus = paste0(gene_symbol_mus, collapse = ","), gene_entrez.id = paste0(gene_entrez.id, collapse = ",")) %>% 
  slice_min(order_by =fisherFDRMin)

write_csv(mp_ma_mus,"results/MP_MA.csv")


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
#Gmp1
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



