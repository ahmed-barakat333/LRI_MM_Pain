## Workflow protocols
## https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html
## https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html

# Needed packages
library(biomaRt)
library(edgeR)
library(limma)

# setting working directory 
setwd("~/Work/PhD/Projects/MoD TrD LRI MM-Pain/Computational Analysis/Ligand_Receptor/Data/MM/PC/hDatasets/GSE153380")

# reading count files
counts_GSE153380 <- read.table("GSE153380_Raw_counts_in_57992_quantified_genes.tsv", header=T)
rownames(counts_GSE153380)<- counts_GSE153380$id
counts_GSE153380$id <- NULL

# number of rows
nrow(counts_GSE153380)

# print number of reads which could be aligned to the different genes
colSums(counts_GSE153380)

# saving files
id_GSE153380 <- row.names(counts_GSE153380)
counts_GSE153380_save <- cbind("id" = id_GSE153380, counts_GSE153380)
write.csv(counts_GSE153380, file="./Unfiltered.Unnormalized.Annotated.Data.csv",row.names = FALSE)

# getting ensembl_gene_id annotation
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
database_gene_symbols <- getBM(
  attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id", values = rownames(counts_GSE153380),
  mart = mart, useCache = F)
database_gene_symbols<- na.omit(database_gene_symbols)
colnames(database_gene_symbols)[1] <- "id"

# filtering by mapping ids
ids_to_keep_GSE153380 <- database_gene_symbols$id
counts_GSE153380_filtered <- counts_GSE153380[ids_to_keep_GSE153380,]
counts_GSE153380_filtered[is.na(counts_GSE153380_filtered)] <- 0
nrow(counts_GSE153380_filtered)
# 28600

# saving files with annotation and filtering genes (ID)
counts_GSE153380_save <- merge(counts_GSE153380_save,database_gene_symbols,by = "id")
counts_GSE153380_save <- counts_GSE153380_save[colnames(counts_GSE153380_save)[c(33,34,1:32)]]
# counts_GSE153380_save$id <- NULL
write.csv(counts_GSE153380_save, file="./Filtered(ID).Unnormalized.Annotated.Data.csv",row.names = FALSE)

# defining groups
gsms_GSE153380 <- "1111111111111111111111111111000"
sml_GSE153380 <- c()
conditions_GSE153380 <- c("HD","MM")
for (i in 1:nchar(gsms_GSE153380)) { sml_GSE153380[i] <- conditions_GSE153380[as.integer(substr(gsms_GSE153380,i,i))+1] }
fl_GSE153380 <- as.factor(sml_GSE153380)

# create DGEList object from the count matrix ??
merged_dge_GSE153380 <- DGEList(counts_GSE153380_filtered)
merged_dge_GSE153380$samples$group <- fl_GSE153380
print(merged_dge_GSE153380$samples)
dim(merged_dge_GSE153380)
# 28600    31

#  filter genes using edgeR, while keeping as many genes as possible with worthwhile counts 
table(rowSums(merged_dge_GSE153380$counts==0))

keep.exprs_GSE153380 <- filterByExpr(merged_dge_GSE153380, group=fl_GSE153380, min.count = 1)
merged_filtered_GSE153380 <- merged_dge_GSE153380[keep.exprs_GSE153380, keep.lib.sizes=FALSE]
dim(merged_filtered_GSE153380)
# 25460    31

table(rowSums(merged_filtered_GSE153380$counts==0))


# saving files with annotation and filtering genes ID and counts
id_GSE153380 <- row.names(merged_filtered_GSE153380)
counts_GSE153380_save <- cbind("id" = id_GSE153380, merged_filtered_GSE153380$counts)
counts_GSE153380_save <- merge(counts_GSE153380_save,database_gene_symbols,by = "id")
counts_GSE153380_save <- counts_GSE153380_save[colnames(counts_GSE153380_save)[c(33,34,1:32)]]
#counts_GSE153380_save$id <- NULL
write.csv(counts_GSE153380_save, file="./Filtered(ID,Counts).Unnormalized.Annotated.Data.csv",row.names = FALSE)

# saving NetworkAnalyst files (cancelled)
counts_GSE153380_save_NA <- rbind(class_GSE153380,counts_GSE153380_save)
colnames(counts_GSE153380_save_NA)[1] <- "#NAME"
write.table(counts_GSE153380_save_NA, file="./NA.Filtered(ID,Counts).Unnormalized.Annotated.Data.txt", row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)

# normalising gene expression distribution using edgeR
merged_norm_GSE153380 <- calcNormFactors(merged_filtered_GSE153380, method = "TMM")

# saving files with annotation and filtering genes (ID), counts and normalization
## https://support.bioconductor.org/p/103747/#103769
merged_norm_GSE153380_counts <- cpm(merged_norm_GSE153380)
counts_GSE153380_save <- cbind("id" = id_GSE153380, merged_norm_GSE153380_counts)
counts_GSE153380_save <- merge(counts_GSE153380_save,database_gene_symbols,by = "id")
counts_GSE153380_save <- counts_GSE153380_save[colnames(counts_GSE153380_save)[c(33,34,1:32)]]
#counts_GSE153380_save$id <- NULL
write.csv(counts_GSE153380_save, file="./Filtered(ID,Counts).Normalized.Annotated.Data_cpm.csv",row.names = FALSE)

# saving NetworkAnalyst files (cancelled)
counts_GSE153380_save_NA <- rbind(class_GSE153380,counts_GSE153380_save)
colnames(counts_GSE153380_save_NA)[1] <- "#NAME"
write.table(counts_GSE153380_save_NA, file="./NA.Filtered(ID,Counts).Normalized.Annotated.Data_cpm.txt", row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)

# [optional] multidimensional scaling (MDS) plot
par(mfrow=c(1,1))
lcpm_norm_GSE153380 <- cpm(merged_norm_GSE153380, log=TRUE)
plotMDS(lcpm_norm_GSE153380, labels=fl_GSE153380, col = as.numeric(fl_GSE153380))
title(main="Sample groups")


# [optional] filter based on variance 
#filtered_var_GSE175384 <- varFilter(lcpm_norm_GSE175384, var.func=IQR, var.cutoff=0.5, filterByQuantile=TRUE)
#rows_to_keep_var_GSE175384 <- rownames(filtered_var_GSE175384)
#merged_norm_filtered_GSE175384 <- merged_norm_GSE175384[rows_to_keep_var_GSE175384,]
#dim(merged_norm_filtered_GSE175384)
# 6012   64

# creating a design matrix and contrasts
designMat.id_GSE153380 <- model.matrix(~0+fl_GSE153380)
colnames(designMat.id_GSE153380) <- levels(fl_GSE153380)
rownames(designMat.id_GSE153380) <- colnames(merged_norm_GSE153380)
print(designMat.id_GSE153380)
cont.matrix.id_GSE153380 <- makeContrasts(MMvsHD = MM-HD,levels=colnames(designMat.id_GSE153380))
print(cont.matrix.id_GSE153380)

# remove heteroscedascity from count data
par(mfrow=c(1,2))
v_GSE153380 <- voom(merged_norm_GSE153380, designMat.id_GSE153380, plot=TRUE)
# counts in v_GSE175384$E
# When operating on a DGEList-object, voom converts raw counts to log-CPM values by automatically extracting library sizes and normalisation factors from x itself. 
# Additional normalisation to log-CPM values can be specified within voom using the normalize.method argument.

# saving files with annotation and filtering genes (ID), counts  and normalization using voom
## see https://support.bioconductor.org/p/62541/
merged_norm_GSE153380_voom_counts <- data.frame(2^(v_GSE153380$E)*v_GSE153380$targets$norm.factors)
counts_GSE153380_save <- cbind("id" = id_GSE153380, merged_norm_GSE153380_voom_counts)
counts_GSE153380_save <- merge(counts_GSE153380_save,database_gene_symbols,by = "id")
counts_GSE153380_save <- counts_GSE153380_save[colnames(counts_GSE153380_save)[c(33,34,1:32)]]
#counts_GSE153380_save$id <- NULL
write.csv(counts_GSE153380_save, file="./Filtered(ID,Counts).Normalized.Annotated.Data_voom.csv",row.names = FALSE)

## saving NetworkAnalyst files (cancelled)
counts_GSE153380_save_NA <- rbind(class_GSE153380,counts_GSE153380_save)
colnames(counts_GSE153380_save_NA)[1] <- "#NAME"
write.table(counts_GSE153380_save_NA, file="./NA.Filtered(ID,Counts).Normalized.Annotated.Data_voom.txt", row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)

## fitting linear models for comparisons of interest
vfit_GSE153380 <- lmFit(v_GSE153380, designMat.id_GSE153380)
vfit_GSE153380 <- contrasts.fit(vfit_GSE153380, contrasts=cont.matrix.id_GSE153380)
efit_GSE153380 <- eBayes(vfit_GSE153380)
plotSA(efit_GSE153380, main="Final model: Mean-variance trend")

# Linear modelling in limma is carried out using the lmFit and contrasts.fit functions originally written for application to microarrays. 
# The functions can be used for both microarray and RNA-seq data and fit a separate model to the expression values for each gene. Next, 
# empirical Bayes moderation is carried out by borrowing information across all the genes to obtain more precise estimates of gene-wise 
# variability (Smyth 2004). The model’s residual variances are plotted against average expression values in the next figure. It can be 
# seen from this plot that the variance is no longer dependent on the mean expression level.

## examining the number of DE genes
dt_fdr_GSE153380 = decideTests(efit_GSE153380, adjust.method	= "fdr", p.value=0.05)
summary(dt_fdr_GSE153380)
#MMvsHD
#Down     196
#NotSig   22377
#Up       13
dt_fdr_lfc_GSE153380 = decideTests(efit_GSE153380, adjust.method	= "fdr", p.value=0.05, lfc = 1)
summary(dt_fdr_lfc_GSE153380)
#MMvsHD
#Down     195
#NotSig   22380
#Up       11
topTab_fdr_GSE153380 <- topTable (efit_GSE153380, number=nrow(efit_GSE153380), coef="MMvsHD", adjust="fdr")
#head(topTab_fdr_GSE153380)
dim(topTab_fdr_GSE153380)
# 22586     6

## add Entrez IDs and save
topTab_fdr_GSE153380_annotated <- merge(as.data.frame(topTab_fdr_GSE153380),database_gene_symbols, by.x = 0, by.y = "id")
write.csv(topTab_fdr_GSE153380_annotated, file="./GSE153380_DEG_limma.csv")

## Create table of significant ones
topTab_fdr_sig_GSE153380 <- subset(topTab_fdr_GSE153380, topTab_fdr_GSE153380$adj.P.Val < 0.05 & (topTab_fdr_GSE153380$logFC < -1 | topTab_fdr_GSE153380$logFC > 1))
dim(topTab_fdr_sig_GSE153380)
# 206    6

## Add Entrez IDs and save
topTab_fdr_sig_GSE153380_annotated <- merge(as.data.frame(topTab_fdr_sig_GSE153380),database_gene_symbols, by.x = 0, by.y = "id")
write.csv(topTab_fdr_sig_GSE153380_annotated, file="./GSE153380_DEG_limma_sig.csv")
