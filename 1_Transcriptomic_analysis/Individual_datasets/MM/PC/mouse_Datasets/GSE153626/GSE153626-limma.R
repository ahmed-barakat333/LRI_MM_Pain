## Workflow protocols
## https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html
## https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html

# Needed packages
library(biomaRt)
library(edgeR)
library(limma)

# setting working directory 
setwd("W:/PhD/Projects/MM Pain/Computation MM TME/Cells/MM/mDatasets MM/GSE153626")

# reading count files
counts_GSE153626 <- read.csv("GSE153626_raw_counts.csv", header=T,row.names = 1)

# number of rows
nrow(counts_GSE153626)

# removing non-relevant samples (AL, MGUS,Plasmablasts)
counts_GSE153626[,c("B3463.BM.V","B3585.BM.V","B3698.BM.V","X4562.BM.VK12653","X4680.BM.VK12653")] <- NULL


# print number of reads which could be aligned to the different genes
print(colSums(counts_GSE153626))

# saving files
id_GSE153626 <- row.names(counts_GSE153626)
counts_GSE153626_save <- cbind("id" = id_GSE153626, counts_GSE153626)
write.csv(counts_GSE153626_save, file="./Unfiltered.Unnormalized.Annotated.Data.csv",row.names = FALSE)

# saving NetworkAnalyst files
class_GSE153626 <- data.frame("#CLASS","HD","HD","HD","MM","MM","MM","MM","MM",stringsAsFactors = FALSE) 
names(class_GSE153626) <- colnames(counts_GSE153626_save)
counts_GSE153626_save_NA <- rbind(class_GSE153626,counts_GSE153626_save)
colnames(counts_GSE153626_save_NA)[1] <- "#NAME"
write.table(counts_GSE153626_save_NA, file="./NA.UnFiltered.Unnormalized.Annotated.Data.txt", row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)

# getting ensembl_gene_id annotation
mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
database_gene_symbols <- getBM(
  attributes = c("mgi_symbol", "entrezgene_id"),
  filters = "mgi_symbol", values = rownames(counts_GSE153626),
  mart = mart, useCache = F)
database_gene_symbols<- na.omit(database_gene_symbols)
colnames(database_gene_symbols)[1] <- "id"

# filtering by mapping ids
ids_to_keep_GSE153626 <- database_gene_symbols$id
counts_GSE153626_filtered <- counts_GSE153626[ids_to_keep_GSE153626, ]
counts_GSE153626_filtered[is.na(counts_GSE153626_filtered)] <- 0
nrow(counts_GSE153626_filtered)
# 21487

# saving files with annotation and filtering genes (ID)
counts_GSE153626_save <- merge(counts_GSE153626_save,database_gene_symbols,by = "id")
counts_GSE153626_save <- counts_GSE153626_save[colnames(counts_GSE153626_save)[c(10,1:9)]]
counts_GSE153626_save$id <- NULL
write.csv(counts_GSE153626_save, file="./Filtered(ID).Unnormalized.Annotated.Data.csv",row.names = FALSE)

# saving NetworkAnalyst files
names(class_GSE153626) <- colnames(counts_GSE153626_save)
counts_GSE153626_save_NA <- rbind(class_GSE153626,counts_GSE153626_save)
colnames(counts_GSE153626_save_NA)[1] <- "#NAME"
write.table(counts_GSE153626_save_NA, file="./NA.Filtered(ID).Unnormalized.Annotated.Data.txt", row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)

# defining groups
gsms_GSE153626 <- "00011111"
sml_GSE153626 <- c()
conditions_GSE153626 <- c("HD","MM")
for (i in 1:nchar(gsms_GSE153626)) { sml_GSE153626[i] <- conditions_GSE153626[as.integer(substr(gsms_GSE153626,i,i))+1] }
fl_GSE153626 <- as.factor(sml_GSE153626)

# create DGEList object from the count matrix ??
merged_dge_GSE153626 <- DGEList(counts_GSE153626_filtered)
merged_dge_GSE153626$samples$group <- fl_GSE153626
print(merged_dge_GSE153626$samples)
dim(merged_dge_GSE153626)
# 21487   8

#  filter genes using edgeR, while keeping as many genes as possible with worthwhile counts 
table(rowSums(merged_dge_GSE153626$counts==0))

keep.exprs_GSE153626 <- filterByExpr(merged_dge_GSE153626, group=fl_GSE153626, min.count = 1)
merged_filtered_GSE153626 <- merged_dge_GSE153626[keep.exprs_GSE153626, keep.lib.sizes=FALSE]
dim(merged_filtered_GSE153626)

table(rowSums(merged_filtered_GSE153626$counts==0))


# saving files with annotation and filtering genes ID and counts
id_GSE153626 <- row.names(merged_filtered_GSE153626)
counts_GSE153626_save <- cbind("id" = id_GSE153626, merged_filtered_GSE153626$counts)
counts_GSE153626_save <- merge(counts_GSE153626_save,database_gene_symbols,by = "id")
counts_GSE153626_save <- counts_GSE153626_save[colnames(counts_GSE153626_save)[c(10,1:9)]]
counts_GSE153626_save$id <- NULL
write.csv(counts_GSE153626_save, file="./Filtered(ID,Counts).Unnormalized.Annotated.Data.csv",row.names = FALSE)

# saving NetworkAnalyst files
counts_GSE153626_save_NA <- rbind(class_GSE153626,counts_GSE153626_save)
colnames(counts_GSE153626_save_NA)[1] <- "#NAME"
write.table(counts_GSE153626_save_NA, file="./NA.Filtered(ID,Counts).Unnormalized.Annotated.Data.txt", row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)

# normalising gene expression distribution using edgeR
merged_norm_GSE153626 <- calcNormFactors(merged_filtered_GSE153626, method = "TMM")

# saving files with annotation and filtering genes (ID), counts and normalization
## https://support.bioconductor.org/p/103747/#103769
merged_norm_GSE153626_counts <- cpm(merged_norm_GSE153626)
counts_GSE153626_save <- cbind("id" = id_GSE153626, merged_norm_GSE153626_counts)
counts_GSE153626_save <- merge(counts_GSE153626_save,database_gene_symbols,by = "id")
counts_GSE153626_save <- counts_GSE153626_save[colnames(counts_GSE153626_save)[c(10,1:9)]]
counts_GSE153626_save$id <- NULL
write.csv(counts_GSE153626_save, file="./Filtered(ID,Counts).Normalized.Annotated.Data_cpm.csv",row.names = FALSE)

# saving NetworkAnalyst files
counts_GSE153626_save_NA <- rbind(class_GSE153626,counts_GSE153626_save)
colnames(counts_GSE153626_save_NA)[1] <- "#NAME"
write.table(counts_GSE153626_save_NA, file="./NA.Filtered(ID,Counts).Normalized.Annotated.Data_cpm.txt", row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)

# [optional] multidimensional scaling (MDS) plot
par(mfrow=c(1,1))
lcpm_norm_GSE153626 <- cpm(merged_norm_GSE153626, log=TRUE)
plotMDS(lcpm_norm_GSE153626, labels=fl_GSE153626, col = as.numeric(fl_GSE153626))
title(main="Sample groups")


# [optional] filter based on variance 
#filtered_var_GSE175384 <- varFilter(lcpm_norm_GSE175384, var.func=IQR, var.cutoff=0.5, filterByQuantile=TRUE)
#rows_to_keep_var_GSE175384 <- rownames(filtered_var_GSE175384)
#merged_norm_filtered_GSE175384 <- merged_norm_GSE175384[rows_to_keep_var_GSE175384,]
#dim(merged_norm_filtered_GSE175384)
# 6012   64

# creating a design matrix and contrasts
designMat.id_GSE153626 <- model.matrix(~0+fl_GSE153626)
colnames(designMat.id_GSE153626) <- levels(fl_GSE153626)
rownames(designMat.id_GSE153626) <- colnames(merged_norm_GSE153626)
print(designMat.id_GSE153626)
cont.matrix.id_GSE153626 <- makeContrasts(MMvsHD = MM-HD,levels=colnames(designMat.id_GSE153626))
print(cont.matrix.id_GSE153626)

# remove heteroscedascity from count data
par(mfrow=c(1,2))
v_GSE153626 <- voom(merged_norm_GSE153626, designMat.id_GSE153626, plot=TRUE)
# counts in v_GSE153626$E
# When operating on a DGEList-object, voom converts raw counts to log-CPM values by automatically extracting library sizes and normalisation factors from x itself. 
# Additional normalisation to log-CPM values can be specified within voom using the normalize.method argument.

# saving files with annotation and filtering genes (ID), counts  and normalization using voom
## see https://support.bioconductor.org/p/62541/
merged_norm_GSE153626_voom_counts <- data.frame(2^(v_GSE153626$E)*v_GSE153626$targets$norm.factors)
counts_GSE153626_save <- cbind("id" = id_GSE153626, merged_norm_GSE153626_voom_counts)
counts_GSE153626_save <- merge(counts_GSE153626_save,database_gene_symbols,by = "id")
counts_GSE153626_save <- counts_GSE153626_save[colnames(counts_GSE153626_save)[c(10,1:9)]]
counts_GSE153626_save$id <- NULL
write.csv(counts_GSE153626_save, file="./Filtered(ID,Counts).Normalized.Annotated.Data_voom.csv",row.names = FALSE)

## saving NetworkAnalyst files
counts_GSE153626_save_NA <- rbind(class_GSE153626,counts_GSE153626_save)
colnames(counts_GSE153626_save_NA)[1] <- "#NAME"
write.table(counts_GSE153626_save_NA, file="./NA.Filtered(ID,Counts).Normalized.Annotated.Data_voom.txt", row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)

## fitting linear models for comparisons of interest
vfit_GSE153626 <- lmFit(v_GSE153626, designMat.id_GSE153626)
vfit_GSE153626 <- contrasts.fit(vfit_GSE153626, contrasts=cont.matrix.id_GSE153626)
efit_GSE153626 <- eBayes(vfit_GSE153626)
plotSA(efit_GSE153626, main="Final model: Mean-variance trend")

# Linear modelling in limma is carried out using the lmFit and contrasts.fit functions originally written for application to microarrays. 
# The functions can be used for both microarray and RNA-seq data and fit a separate model to the expression values for each gene. Next, 
# empirical Bayes moderation is carried out by borrowing information across all the genes to obtain more precise estimates of gene-wise 
# variability (Smyth 2004). The model’s residual variances are plotted against average expression values in the next figure. It can be 
# seen from this plot that the variance is no longer dependent on the mean expression level.

## examining the number of DE genes
dt_fdr_GSE153626 = decideTests(efit_GSE153626, adjust.method	= "fdr", p.value=0.05)
summary(dt_fdr_GSE153626)
#MMvsHD
#Down     1909
#NotSig   9958
#Up       1388
dt_fdr_lfc_GSE153626 = decideTests(efit_GSE153626, adjust.method	= "fdr", p.value=0.05, lfc = 1)
summary(dt_fdr_lfc_GSE153626)
#MMvsHD
#Down     1672
#NotSig   10613
#Up       970
topTab_fdr_GSE153626 <- topTable (efit_GSE153626, number=nrow(efit_GSE153626), coef="MMvsHD", adjust="fdr")
#head(topTab_fdr_GSE153626)
dim(topTab_fdr_GSE153626)
# 13255     6

## add Entrez IDs and save
topTab_fdr_GSE153626_annotated <- merge(as.data.frame(topTab_fdr_GSE153626),database_gene_symbols, by.x = 0, by.y = "id")
write.csv(topTab_fdr_GSE153626_annotated, file="./GSE153626_DEG_limma.csv")

## Create table of significant ones
topTab_fdr_sig_GSE153626 <- subset(topTab_fdr_GSE153626, topTab_fdr_GSE153626$adj.P.Val < 0.05 & (topTab_fdr_GSE153626$logFC < -1 | topTab_fdr_GSE153626$logFC > 1))
dim(topTab_fdr_sig_GSE153626)
# 2642    6

## Add Entrez IDs and save
topTab_fdr_sig_GSE153626_annotated <- merge(as.data.frame(topTab_fdr_sig_GSE153626),database_gene_symbols, by.x = 0, by.y = "id")
write.csv(topTab_fdr_sig_GSE153626_annotated, file="./GSE153626_DEG_limma_sig.csv")
