## Workflow protocols
## https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html
## https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html

# Needed packages
library(biomaRt)
library(edgeR)
library(limma)
library(dplyr)


# setting working directory 
setwd("~/Work/PhD/Projects/MM Pain/Computational Analysis/Ligand_Receptor/Data/MM/Treg/mDatasets/GSE109533")

# reading count files
counts_1 <- read.table("./data/GSM2945891_Treg_CT1.rsem.genes.results.txt", sep="\t", header = TRUE)
counts_1 <- counts_1[ , c("gene_id", "expected_count")]
counts_1 <- counts_1[complete.cases(counts_1), ] 
colnames(counts_1)[2] <- "GSM2945891_Treg_CT1.rsem.genes.results"

counts_2 <- read.table("./data/GSM2945892_Treg_CT3.rsem.genes.results.txt", sep="\t", header = TRUE)
counts_2 <- counts_2[ , c("gene_id", "expected_count")]
counts_2 <- counts_2[complete.cases(counts_2), ] 
colnames(counts_2)[2] <- "GSM2945892_Treg_CT3.rsem.genes.results"

counts_3 <- read.table("./data/GSM2945893_Treg_CT5.rsem.genes.results.txt", sep="\t", header = TRUE)
counts_3 <- counts_3[ , c("gene_id", "expected_count")]
counts_3 <- counts_3[complete.cases(counts_3), ] 
colnames(counts_3)[2] <- "GSM2945893_Treg_CT5.rsem.genes.results"

counts_4 <- read.table("./data/GSM2945888_Treg_MM3.rsem.genes.results.txt", sep="\t", header = TRUE)
counts_4 <- counts_4[ , c("gene_id", "expected_count")]
counts_4 <- counts_4[complete.cases(counts_4), ] 
colnames(counts_4)[2] <- "GSM2945888_Treg_MM3.rsem.genes.results"

counts_5 <- read.table("./data/GSM2945889_Treg_MM4.rsem.genes.results.txt", sep="\t", header = TRUE)
counts_5 <- counts_5[ , c("gene_id", "expected_count")]
counts_5 <- counts_5[complete.cases(counts_5), ] 
colnames(counts_5)[2] <- "GSM2945889_Treg_MM4.rsem.genes.results"

counts_6 <- read.table("./data/GSM2945890_Treg_MM5.rsem.genes.results.txt", sep="\t", header = TRUE)
counts_6 <- counts_6[ , c("gene_id", "expected_count")]
counts_6 <- counts_6[complete.cases(counts_6), ] 
colnames(counts_6)[2] <- "GSM2945890_Treg_MM5.rsem.genes.results"


counts_GSE109533 <- Reduce(function(...) merge(..., by="gene_id", all=TRUE), list(counts_1, counts_2, counts_3,counts_4,counts_5,counts_6))
remove(counts_1,counts_2,counts_3,counts_4,counts_5,counts_6)

# moving first column to rownames
rownames(counts_GSE109533)<- counts_GSE109533$gene_id
counts_GSE109533$gene_id <- NULL

# number of rows
nrow(counts_GSE109533)


# print number of reads which could be aligned to the different genes
print(colSums(counts_GSE109533))

# stripping version number from Ensembl gene ID
library(stringr)
rownames(counts_GSE109533) <- str_replace(rownames(counts_GSE109533),
                        pattern = ".[0-9]+$",
                        replacement = "")

# saving files
id_GSE109533 <- row.names(counts_GSE109533)
counts_GSE109533_save <- cbind("id" = id_GSE109533, counts_GSE109533)
write.csv(counts_GSE109533_save, file="./Unfiltered.Unnormalized.Annotated.Data.csv",row.names = FALSE)

# saving NetworkAnalyst files (cancelled)
class_GSE109533 <- data.frame("#CLASS","HD","HD","HD","MM","MM","MM",stringsAsFactors = FALSE) 
names(class_GSE109533) <- colnames(counts_GSE109533_save)
counts_GSE109533_save_NA <- rbind(class_GSE109533,counts_GSE109533_save)
colnames(counts_GSE109533_save_NA)[1] <- "#NAME"
write.table(counts_GSE109533_save_NA, file="./NA.UnFiltered.Unnormalized.Annotated.Data.txt", row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)

# getting ensembl_gene_id annotation
mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
database_gene_symbols <- getBM(
  attributes = c("ensembl_gene_id", "entrezgene_id","mgi_symbol"),
  filters = "ensembl_gene_id", values = rownames(counts_GSE109533),
  mart = mart, useCache = F)
database_gene_symbols<- na.omit(database_gene_symbols)
colnames(database_gene_symbols)[1] <- "id"

# filtering by mapping ids
ids_to_keep_GSE109533 <- database_gene_symbols$id
counts_GSE109533_filtered <- counts_GSE109533[ids_to_keep_GSE109533, ]
counts_GSE109533_filtered[is.na(counts_GSE109533_filtered)] <- 0
nrow(counts_GSE109533_filtered)
# 25745

# saving files with annotation and filtering genes (ID)
counts_GSE109533_save <- merge(counts_GSE109533_save,database_gene_symbols,by = "id")
counts_GSE109533_save <- counts_GSE109533_save[colnames(counts_GSE109533_save)[c(8,1:7)]]
counts_GSE109533_save$id <- NULL
write.csv(counts_GSE109533_save, file="./Filtered(ID).Unnormalized.Annotated.Data.csv",row.names = FALSE)

# saving NetworkAnalyst files (cancelled)
names(class_GSE109533) <- colnames(counts_GSE109533_save)
counts_GSE109533_save_NA <- rbind(class_GSE109533,counts_GSE109533_save)
colnames(counts_GSE109533_save_NA)[1] <- "#NAME"
write.table(counts_GSE109533_save_NA, file="./NA.Filtered(ID).Unnormalized.Annotated.Data.txt", row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)

# defining groups
gsms_GSE109533 <- "000111"
sml_GSE109533 <- c()
conditions_GSE109533 <- c("HD","MM")
for (i in 1:nchar(gsms_GSE109533)) { sml_GSE109533[i] <- conditions_GSE109533[as.integer(substr(gsms_GSE109533,i,i))+1] }
fl_GSE109533 <- as.factor(sml_GSE109533)

# create DGEList object from the count matrix ??
merged_dge_GSE109533 <- DGEList(counts_GSE109533_filtered)
merged_dge_GSE109533$samples$group <- fl_GSE109533
print(merged_dge_GSE109533$samples)
dim(merged_dge_GSE109533)
# 25745     6

#  filter genes using edgeR, while keeping as many genes as possible with worthwhile counts 
table(rowSums(merged_dge_GSE109533$counts==0))

keep.exprs_GSE109533 <- filterByExpr(merged_dge_GSE109533, group=fl_GSE109533, min.count = 1)
merged_filtered_GSE109533 <- merged_dge_GSE109533[keep.exprs_GSE109533, keep.lib.sizes=FALSE]
dim(merged_filtered_GSE109533)

table(rowSums(merged_filtered_GSE109533$counts==0))


# saving files with annotation and filtering genes ID and counts
id_GSE109533 <- row.names(merged_filtered_GSE109533)
counts_GSE109533_save <- cbind("id" = id_GSE109533, merged_filtered_GSE109533$counts)
counts_GSE109533_save <- merge(counts_GSE109533_save,database_gene_symbols,by = "id")
counts_GSE109533_save <- counts_GSE109533_save[colnames(counts_GSE109533_save)[c(8,1:7)]]
counts_GSE109533_save$id <- NULL
write.csv(counts_GSE109533_save, file="./Filtered(ID,Counts).Unnormalized.Annotated.Data.csv",row.names = FALSE)

# saving NetworkAnalyst files (cancelled)
counts_GSE109533_save_NA <- rbind(class_GSE109533,counts_GSE109533_save)
colnames(counts_GSE109533_save_NA)[1] <- "#NAME"
write.table(counts_GSE109533_save_NA, file="./NA.Filtered(ID,Counts).Unnormalized.Annotated.Data.txt", row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)

# normalising gene expression distribution using edgeR
merged_norm_GSE109533 <- calcNormFactors(merged_filtered_GSE109533, method = "TMM")

# saving files with annotation and filtering genes (ID), counts and normalization
## https://support.bioconductor.org/p/103747/#103769
merged_norm_GSE109533_counts <- cpm(merged_norm_GSE109533)
counts_GSE109533_save <- cbind("id" = id_GSE109533, merged_norm_GSE109533_counts)
counts_GSE109533_save <- merge(counts_GSE109533_save,database_gene_symbols,by = "id")
counts_GSE109533_save <- counts_GSE109533_save[colnames(counts_GSE109533_save)[c(8,1:7)]]
counts_GSE109533_save$id <- NULL
write.csv(counts_GSE109533_save, file="./Filtered(ID,Counts).Normalized.Annotated.Data_cpm.csv",row.names = FALSE)

# saving NetworkAnalyst files (cancelled)
counts_GSE109533_save_NA <- rbind(class_GSE109533,counts_GSE109533_save)
colnames(counts_GSE109533_save_NA)[1] <- "#NAME"
write.table(counts_GSE109533_save_NA, file="./NA.Filtered(ID,Counts).Normalized.Annotated.Data_cpm.txt", row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)

# [optional] multidimensional scaling (MDS) plot
par(mfrow=c(1,1))
lcpm_norm_GSE109533 <- cpm(merged_norm_GSE109533, log=TRUE)
plotMDS(lcpm_norm_GSE109533, labels=fl_GSE109533, col = as.numeric(fl_GSE109533))
title(main="Sample groups")


# [optional] filter based on variance 
#filtered_var_GSE175384 <- varFilter(lcpm_norm_GSE175384, var.func=IQR, var.cutoff=0.5, filterByQuantile=TRUE)
#rows_to_keep_var_GSE175384 <- rownames(filtered_var_GSE175384)
#merged_norm_filtered_GSE175384 <- merged_norm_GSE175384[rows_to_keep_var_GSE175384,]
#dim(merged_norm_filtered_GSE175384)
# 6012   64

# creating a design matrix and contrasts
designMat.id_GSE109533 <- model.matrix(~0+fl_GSE109533)
colnames(designMat.id_GSE109533) <- levels(fl_GSE109533)
rownames(designMat.id_GSE109533) <- colnames(merged_norm_GSE109533)
print(designMat.id_GSE109533)
cont.matrix.id_GSE109533 <- makeContrasts(MMvsHD = MM-HD,levels=colnames(designMat.id_GSE109533))
print(cont.matrix.id_GSE109533)

# remove heteroscedascity from count data
par(mfrow=c(1,2))
v_GSE109533 <- voom(merged_norm_GSE109533, designMat.id_GSE109533, plot=TRUE)
# counts in v_GSE109533$E
# When operating on a DGEList-object, voom converts raw counts to log-CPM values by automatically extracting library sizes and normalisation factors from x itself. 
# Additional normalisation to log-CPM values can be specified within voom using the normalize.method argument.

# saving files with annotation and filtering genes (ID), counts  and normalization using voom
## see https://support.bioconductor.org/p/62541/
merged_norm_GSE109533_voom_counts <- data.frame(2^(v_GSE109533$E)*v_GSE109533$targets$norm.factors)
counts_GSE109533_save <- cbind("id" = id_GSE109533, merged_norm_GSE109533_voom_counts)
counts_GSE109533_save <- merge(counts_GSE109533_save,database_gene_symbols,by = "id")
counts_GSE109533_save <- counts_GSE109533_save[colnames(counts_GSE109533_save)[c(8,1:7)]]
counts_GSE109533_save$id <- NULL
write.csv(counts_GSE109533_save, file="./Filtered(ID,Counts).Normalized.Annotated.Data_voom.csv",row.names = FALSE)

## saving NetworkAnalyst files (cancelled)
counts_GSE109533_save_NA <- rbind(class_GSE109533,counts_GSE109533_save)
colnames(counts_GSE109533_save_NA)[1] <- "#NAME"
write.table(counts_GSE109533_save_NA, file="./NA.Filtered(ID,Counts).Normalized.Annotated.Data_voom.txt", row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)

## fitting linear models for comparisons of interest
vfit_GSE109533 <- lmFit(v_GSE109533, designMat.id_GSE109533)
vfit_GSE109533 <- contrasts.fit(vfit_GSE109533, contrasts=cont.matrix.id_GSE109533)
efit_GSE109533 <- eBayes(vfit_GSE109533)
plotSA(efit_GSE109533, main="Final model: Mean-variance trend")

# Linear modelling in limma is carried out using the lmFit and contrasts.fit functions originally written for application to microarrays. 
# The functions can be used for both microarray and RNA-seq data and fit a separate model to the expression values for each gene. Next, 
# empirical Bayes moderation is carried out by borrowing information across all the genes to obtain more precise estimates of gene-wise 
# variability (Smyth 2004). The model’s residual variances are plotted against average expression values in the next figure. It can be 
# seen from this plot that the variance is no longer dependent on the mean expression level.

## examining the number of DE genes
dt_fdr_GSE109533 = decideTests(efit_GSE109533, adjust.method	= "fdr", p.value=0.05)
summary(dt_fdr_GSE109533)
#MMvsHD
#Down     680
#NotSig   12005
#Up       784
dt_fdr_lfc_GSE109533 = decideTests(efit_GSE109533, adjust.method	= "fdr", p.value=0.05, lfc = 1)
summary(dt_fdr_lfc_GSE109533)
#MMvsHD
#Down     67
#NotSig   13187
#Up       218
topTab_fdr_GSE109533 <- topTable (efit_GSE109533, number=nrow(efit_GSE109533), coef="MMvsHD", adjust="fdr")
#head(topTab_fdr_GSE109533)
dim(topTab_fdr_GSE109533)
# 13469         6

## add Entrez IDs and save
topTab_fdr_GSE109533_annotated <- merge(as.data.frame(topTab_fdr_GSE109533),database_gene_symbols, by.x = 0, by.y = "id")
write.csv(topTab_fdr_GSE109533_annotated, file="./GSE109533_DEG_limma.csv")

## Create table of significant ones
topTab_fdr_sig_GSE109533 <- subset(topTab_fdr_GSE109533, topTab_fdr_GSE109533$adj.P.Val < 0.05 & (topTab_fdr_GSE109533$logFC < -1 | topTab_fdr_GSE109533$logFC > 1))
dim(topTab_fdr_sig_GSE109533)
# 285       6

## Add Entrez IDs and save
topTab_fdr_sig_GSE109533_annotated <- merge(as.data.frame(topTab_fdr_sig_GSE109533),database_gene_symbols, by.x = 0, by.y = "id")
write.csv(topTab_fdr_sig_GSE109533_annotated, file="./GSE109533_DEG_limma_sig.csv")


# ortholog mapping from mouse to human of topTab_fdr_GSE109533_annotated
musGenes <- topTab_fdr_GSE109533_annotated$mgi_symbol

convertMouseGeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl",host="http://dec2021.archive.ensembl.org")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl",host="http://dec2021.archive.ensembl.org")
  
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  humanx <- unique(genesV2[,c(1,2)])
  
  return(humanx)
}

gene_symbol_mapped.hs <- convertMouseGeneList(musGenes)


colnames(topTab_fdr_GSE109533_annotated)[9] <- "MGI.symbol"
topTab_fdr_GSE109533_annotated_hs.mapped <- merge(topTab_fdr_GSE109533_annotated,gene_symbol_mapped.hs, by = "MGI.symbol")


# check for duplicate symbols
sum(duplicated(topTab_fdr_GSE109533_annotated_hs.mapped$HGNC.symbol))
sum(duplicated(topTab_fdr_GSE109533_annotated_hs.mapped$MGI.symbol))
## There are many to many mapping

# summarize multiple MGI.symbol mapping to same HGNC.symbol 
colnames(topTab_fdr_GSE109533_annotated_hs.mapped)[2] <- "ensembl_id"
topTab_fdr_GSE109533_annotated_hs.mapped <- topTab_fdr_GSE109533_annotated_hs.mapped %>% 
  group_by(HGNC.symbol) %>% 
  mutate(MGI.symbol = paste0(MGI.symbol, collapse = ","), entrezgene_id = paste0(entrezgene_id, collapse = ","), ensembl_id = paste0(ensembl_id, collapse = ",")) %>% 
  slice_min(order_by = adj.P.Val) %>% as.data.frame()

write.csv(topTab_fdr_GSE109533_annotated_hs.mapped,"GSE109533_DEG_limma.csv")

# ortholog mapping from mouse to human of topTab_fdr_sig_GSE109533_annotated
musGenes <- topTab_fdr_sig_GSE109533_annotated$mgi_symbol

gene_symbol_mapped.hs <- convertMouseGeneList(musGenes)

colnames(topTab_fdr_sig_GSE109533_annotated)[9] <- "MGI.symbol"
topTab_fdr_sig_GSE109533_annotated_hs.mapped <- merge(topTab_fdr_sig_GSE109533_annotated,gene_symbol_mapped.hs, by = "MGI.symbol")


# check for duplicate symbols
sum(duplicated(topTab_fdr_sig_GSE109533_annotated_hs.mapped$HGNC.symbol))
sum(duplicated(topTab_fdr_sig_GSE109533_annotated_hs.mapped$MGI.symbol))
## There are many to many mapping

# summarize multiple MGI.symbol mapping to same HGNC.symbol 
colnames(topTab_fdr_sig_GSE109533_annotated_hs.mapped)[2] <- "ensembl_id"
topTab_fdr_sig_GSE109533_annotated_hs.mapped <- topTab_fdr_sig_GSE109533_annotated_hs.mapped %>% 
  group_by(HGNC.symbol) %>% 
  mutate(MGI.symbol = paste0(MGI.symbol, collapse = ","), entrezgene_id = paste0(entrezgene_id, collapse = ","), ensembl_id = paste0(ensembl_id, collapse = ",")) %>% 
  slice_min(order_by = adj.P.Val) %>% as.data.frame()

write.csv(topTab_fdr_sig_GSE109533_annotated_hs.mapped,"GSE109533_DEG_limma_sig.csv")






















































