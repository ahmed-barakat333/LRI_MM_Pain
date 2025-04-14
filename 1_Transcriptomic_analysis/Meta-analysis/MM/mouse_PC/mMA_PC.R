# set working directory
setwd("~/Work/PhD/Projects/MM Pain/Computational Analysis/Ligand_Receptor/Data/MM/PC/mDatasets/MA")

# load packages
library(MetaIntegrator)
library(org.Mm.eg.db)
library(annotate)
library(tidyverse)
library(dplyr)


# prepare GSE111921 
## read file 
GSE111921 <- read.csv("data/GSE111921 Normalized.Filtered(ID).Annotated.Data.csv")

## prepare datasetObject
expr <- data.matrix(GSE111921[,c(4:49)])
rownames(expr) <- GSE111921$PROBEID

keys <- as.character(GSE111921$ENTREZID)
names(keys) <- GSE111921$PROBEID

class <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0)
names(class) <- colnames(GSE111921)[-c(1:3)]

pheno <- data.frame(c("MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM",	"HD",	"HD"	,"HD"))
rownames(pheno) <- colnames(GSE111921)[-c(1:3)]

GSE111921_datasetObject <- list(expr = expr, 
                               keys = keys,
                               class = class, pheno = pheno, formattedName = "GSE111921")

## remove objects 
rm(GSE111921,expr,pheno,class,keys) 
 
## check for completeness
checkDataObject(GSE111921_datasetObject, "Dataset")


# prepare GSE153626
## read file
GSE153626 <- read.csv("data/GSE153626_Filtered(ID,Counts).Normalized.Annotated.Data_voom.csv")

## prepare datasetObject
expr <- data.matrix(GSE153626[,c(2:9)])
rownames(expr) <- GSE153626$entrezgene_id

keys <- as.character(GSE153626$entrezgene_id)
names(keys) <- GSE153626$entrezgene_id

class <- c(0,0,0,1,1,1,1,1)
names(class) <- colnames(GSE153626)[-1]

pheno <- data.frame(c("HD","HD","HD","MM","MM","MM","MM","MM"))
rownames(pheno) <- colnames(GSE153626)[-1]

GSE153626_datasetObject <- list(expr = expr, 
                               keys = keys,
                               class = class, pheno = pheno, formattedName = "GSE153626")

## remove objects 
rm(GSE153626,expr,pheno,class,keys) 

## check for completeness
checkDataObject(GSE153626_datasetObject, "Dataset")

# create metaObject
pc_datasets <- list(GSE111921_datasetObject,GSE153626_datasetObject)
names(pc_datasets) = c(GSE111921_datasetObject$formattedName, GSE153626_datasetObject$formattedName)
MetaObj=list() 
MetaObj$originalData <- pc_datasets

## check for completeness
checkDataObject(MetaObj, "Meta", "Pre-Analysis")


## run meta-analysis
MetaObj <- runMetaAnalysis(MetaObj, maxCores=1)

MetaObj <- filterGenes(MetaObj, isLeaveOneOut = FALSE #, FDRThresh = 0.05
                       )

## get resutls 
pc_ma_mus <- MetaObj[["metaAnalysis"]][["pooledResults"]]

pc_ma_mus$gene_entrez.id <- rownames(pc_ma_mus)


pc_ma_mus$gene_symbol_mus <-getSYMBOL(pc_ma_mus$gene_entrez.id, data='org.Mm.eg')
sum(duplicated(pc_ma_mus$gene_symbol_mus))

pc_ma_mus <- data.frame(pc_ma_mus %>% 
                      group_by(gene_symbol_mus) %>% 
                      summarise_all(funs(max)))

pc_ma_mus <- pc_ma_mus %>%
  select(gene_entrez.id, everything())

# ortholog mapping from mouse to human
musGenes <- pc_ma_mus$gene_symbol_mus

convertMouseGeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl",host="http://dec2021.archive.ensembl.org")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl",host="http://dec2021.archive.ensembl.org")
  
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  humanx <- unique(genesV2[,c(1,2)])
  
  return(humanx)
}

gene_symbol_mapped.hs <- convertMouseGeneList(musGenes)

#mouse_human_genes = read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")

#convert_mouse_to_human <- function(gene_list){
  
  output = c()
  
  for(gene in gene_list){
    class_key = (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name=="mouse, laboratory"))[['DB.Class.Key']]
    if(!identical(class_key, integer(0)) ){
      human_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="human"))[,"Symbol"]
      for(human_gene in human_genes){
        output = append(output,human_gene)
      }
    }
  }
  
  return (output)
#}

#gene_symbol_mapped.hs <- convert_mouse_to_human(musGenes)


colnames(gene_symbol_mapped.hs)[1] <- "gene_symbol_mus"
colnames(gene_symbol_mapped.hs)[2] <- "gene_symbol_hs.mapped"
pc_ma_mus <- merge(gene_symbol_mapped.hs,pc_ma_mus, by = "gene_symbol_mus")


# check for duplicate symbols
sum(duplicated(pc_ma_mus$gene_symbol_hs.mapped))
sum(duplicated(pc_ma_mus$gene_symbol_mus))
## There are many to many mapping

# summarize multiple MGI.symbol mapping to same HGNC.symbol 
pc_ma_mus <- pc_ma_mus %>% 
  group_by(gene_symbol_hs.mapped) %>% 
  mutate(gene_symbol_mus = paste0(gene_symbol_mus, collapse = ","), gene_entrez.id = paste0(gene_entrez.id, collapse = ",")) %>% 
  summarise_all(funs(max))

write_csv(pc_ma_mus,"results/mPC_MA.csv")






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



