## Workflow protocols
## https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html
## https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html

# Needed packages
library(biomaRt)
library(edgeR)
library(limma)

# setting working directory 
setwd("W:/PhD/Projects/MM Pain-Bone/Computation/Cells/MP/mDatasets MP/GSE176385")

# reading count files
counts_GSE176385 <- read.table("GSE176385_TAM_count_matrix.txt", header=T,row.names = 1)

# number of rows
nrow(counts_GSE176385)

# removing non-relevant samples (samples not in contact with MM)
counts_GSE176385[,c("Mf1" , "Mf3")] <- NULL


# print number of reads which could be aligned to the different genes
print(colSums(counts_GSE176385))

# stripping version number from Ensembl gene ID
library(stringr)
rownames(counts_GSE176385) <- str_replace(rownames(counts_GSE176385),
                        pattern = ".[0-9]+$",
                        replacement = "")

# saving files
id_GSE176385 <- row.names(counts_GSE176385)
counts_GSE176385_save <- cbind("id" = id_GSE176385, counts_GSE176385)
write.csv(counts_GSE176385_save, file="./Unfiltered.Unnormalized.Annotated.Data.csv",row.names = FALSE)

# saving NetworkAnalyst files
class_GSE176385 <- data.frame("#CLASS","MM","MM","HD","HD",stringsAsFactors = FALSE) 
names(class_GSE176385) <- colnames(counts_GSE176385_save)
counts_GSE176385_save_NA <- rbind(class_GSE176385,counts_GSE176385_save)
colnames(counts_GSE176385_save_NA)[1] <- "#NAME"
write.table(counts_GSE176385_save_NA, file="./NA.UnFiltered.Unnormalized.Annotated.Data.txt", row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)

# getting ensembl_gene_id annotation
mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
database_gene_symbols <- getBM(
  attributes = c("ensembl_gene_id", "entrezgene_id"),
  filters = "ensembl_gene_id", values = rownames(counts_GSE176385),
  mart = mart, useCache = F)
database_gene_symbols<- na.omit(database_gene_symbols)
colnames(database_gene_symbols)[1] <- "id"

# filtering by mapping ids
ids_to_keep_GSE176385 <- database_gene_symbols$id
counts_GSE176385_filtered <- counts_GSE176385[ids_to_keep_GSE176385, ]
counts_GSE176385_filtered[is.na(counts_GSE176385_filtered)] <- 0
nrow(counts_GSE176385_filtered)
# 25783

# saving files with annotation and filtering genes (ID)
counts_GSE176385_save <- merge(counts_GSE176385_save,database_gene_symbols,by = "id")
counts_GSE176385_save <- counts_GSE176385_save[colnames(counts_GSE176385_save)[c(6,1:5)]]
counts_GSE176385_save$id <- NULL
write.csv(counts_GSE176385_save, file="./Filtered(ID).Unnormalized.Annotated.Data.csv",row.names = FALSE)

# saving NetworkAnalyst files
names(class_GSE176385) <- colnames(counts_GSE176385_save)
counts_GSE176385_save_NA <- rbind(class_GSE176385,counts_GSE176385_save)
colnames(counts_GSE176385_save_NA)[1] <- "#NAME"
write.table(counts_GSE176385_save_NA, file="./NA.Filtered(ID).Unnormalized.Annotated.Data.txt", row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)

# defining groups
gsms_GSE176385 <- "1100"
sml_GSE176385 <- c()
conditions_GSE176385 <- c("HD","MM")
for (i in 1:nchar(gsms_GSE176385)) { sml_GSE176385[i] <- conditions_GSE176385[as.integer(substr(gsms_GSE176385,i,i))+1] }
fl_GSE176385 <- as.factor(sml_GSE176385)

# create DGEList object from the count matrix ??
merged_dge_GSE176385 <- DGEList(counts_GSE176385_filtered)
merged_dge_GSE176385$samples$group <- fl_GSE176385
print(merged_dge_GSE176385$samples)
dim(merged_dge_GSE176385)
# 25783     4

#  filter genes using edgeR, while keeping as many genes as possible with worthwhile counts 
table(rowSums(merged_dge_GSE176385$counts==0))

keep.exprs_GSE176385 <- filterByExpr(merged_dge_GSE176385, group=fl_GSE176385, min.count = 1)
merged_filtered_GSE176385 <- merged_dge_GSE176385[keep.exprs_GSE176385, keep.lib.sizes=FALSE]
dim(merged_filtered_GSE176385)

table(rowSums(merged_filtered_GSE176385$counts==0))


# saving files with annotation and filtering genes ID and counts
id_GSE176385 <- row.names(merged_filtered_GSE176385)
counts_GSE176385_save <- cbind("id" = id_GSE176385, merged_filtered_GSE176385$counts)
counts_GSE176385_save <- merge(counts_GSE176385_save,database_gene_symbols,by = "id")
counts_GSE176385_save <- counts_GSE176385_save[colnames(counts_GSE176385_save)[c(6,1:5)]]
counts_GSE176385_save$id <- NULL
write.csv(counts_GSE176385_save, file="./Filtered(ID,Counts).Unnormalized.Annotated.Data.csv",row.names = FALSE)

# saving NetworkAnalyst files
counts_GSE176385_save_NA <- rbind(class_GSE176385,counts_GSE176385_save)
colnames(counts_GSE176385_save_NA)[1] <- "#NAME"
write.table(counts_GSE176385_save_NA, file="./NA.Filtered(ID,Counts).Unnormalized.Annotated.Data.txt", row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)

# normalising gene expression distribution using edgeR
merged_norm_GSE176385 <- calcNormFactors(merged_filtered_GSE176385, method = "TMM")

# saving files with annotation and filtering genes (ID), counts and normalization
## https://support.bioconductor.org/p/103747/#103769
merged_norm_GSE176385_counts <- cpm(merged_norm_GSE176385)
counts_GSE176385_save <- cbind("id" = id_GSE176385, merged_norm_GSE176385_counts)
counts_GSE176385_save <- merge(counts_GSE176385_save,database_gene_symbols,by = "id")
counts_GSE176385_save <- counts_GSE176385_save[colnames(counts_GSE176385_save)[c(6,1:5)]]
counts_GSE176385_save$id <- NULL
write.csv(counts_GSE176385_save, file="./Filtered(ID,Counts).Normalized.Annotated.Data_cpm.csv",row.names = FALSE)

# saving NetworkAnalyst files
counts_GSE176385_save_NA <- rbind(class_GSE176385,counts_GSE176385_save)
colnames(counts_GSE176385_save_NA)[1] <- "#NAME"
write.table(counts_GSE176385_save_NA, file="./NA.Filtered(ID,Counts).Normalized.Annotated.Data_cpm.txt", row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)

# [optional] multidimensional scaling (MDS) plot
par(mfrow=c(1,1))
lcpm_norm_GSE176385 <- cpm(merged_norm_GSE176385, log=TRUE)
plotMDS(lcpm_norm_GSE176385, labels=fl_GSE176385, col = as.numeric(fl_GSE176385))
title(main="Sample groups")


# [optional] filter based on variance 
#filtered_var_GSE175384 <- varFilter(lcpm_norm_GSE175384, var.func=IQR, var.cutoff=0.5, filterByQuantile=TRUE)
#rows_to_keep_var_GSE175384 <- rownames(filtered_var_GSE175384)
#merged_norm_filtered_GSE175384 <- merged_norm_GSE175384[rows_to_keep_var_GSE175384,]
#dim(merged_norm_filtered_GSE175384)
# 6012   64

# creating a design matrix and contrasts
designMat.id_GSE176385 <- model.matrix(~0+fl_GSE176385)
colnames(designMat.id_GSE176385) <- levels(fl_GSE176385)
rownames(designMat.id_GSE176385) <- colnames(merged_norm_GSE176385)
print(designMat.id_GSE176385)
cont.matrix.id_GSE176385 <- makeContrasts(MMvsHD = MM-HD,levels=colnames(designMat.id_GSE176385))
print(cont.matrix.id_GSE176385)

# remove heteroscedascity from count data
par(mfrow=c(1,2))
v_GSE176385 <- voom(merged_norm_GSE176385, designMat.id_GSE176385, plot=TRUE)
# counts in v_GSE176385$E
# When operating on a DGEList-object, voom converts raw counts to log-CPM values by automatically extracting library sizes and normalisation factors from x itself. 
# Additional normalisation to log-CPM values can be specified within voom using the normalize.method argument.

# saving files with annotation and filtering genes (ID), counts  and normalization using voom
## see https://support.bioconductor.org/p/62541/
merged_norm_GSE176385_voom_counts <- data.frame(2^(v_GSE176385$E)*v_GSE176385$targets$norm.factors)
counts_GSE176385_save <- cbind("id" = id_GSE176385, merged_norm_GSE176385_voom_counts)
counts_GSE176385_save <- merge(counts_GSE176385_save,database_gene_symbols,by = "id")
counts_GSE176385_save <- counts_GSE176385_save[colnames(counts_GSE176385_save)[c(6,1:5)]]
counts_GSE176385_save$id <- NULL
write.csv(counts_GSE176385_save, file="./Filtered(ID,Counts).Normalized.Annotated.Data_voom.csv",row.names = FALSE)

## saving NetworkAnalyst files
counts_GSE176385_save_NA <- rbind(class_GSE176385,counts_GSE176385_save)
colnames(counts_GSE176385_save_NA)[1] <- "#NAME"
write.table(counts_GSE176385_save_NA, file="./NA.Filtered(ID,Counts).Normalized.Annotated.Data_voom.txt", row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)

## fitting linear models for comparisons of interest
vfit_GSE176385 <- lmFit(v_GSE176385, designMat.id_GSE176385)
vfit_GSE176385 <- contrasts.fit(vfit_GSE176385, contrasts=cont.matrix.id_GSE176385)
efit_GSE176385 <- eBayes(vfit_GSE176385)
plotSA(efit_GSE176385, main="Final model: Mean-variance trend")

# Linear modelling in limma is carried out using the lmFit and contrasts.fit functions originally written for application to microarrays. 
# The functions can be used for both microarray and RNA-seq data and fit a separate model to the expression values for each gene. Next, 
# empirical Bayes moderation is carried out by borrowing information across all the genes to obtain more precise estimates of gene-wise 
# variability (Smyth 2004). The model’s residual variances are plotted against average expression values in the next figure. It can be 
# seen from this plot that the variance is no longer dependent on the mean expression level.

## examining the number of DE genes
dt_fdr_GSE176385 = decideTests(efit_GSE176385, adjust.method	= "fdr", p.value=0.05)
summary(dt_fdr_GSE176385)
#MMvsHD
#Down     160
#NotSig   13762
#Up       117
dt_fdr_lfc_GSE176385 = decideTests(efit_GSE176385, adjust.method	= "fdr", p.value=0.05, lfc = 1)
summary(dt_fdr_lfc_GSE176385)
#MMvsHD
#Down     160
#NotSig   13762
#Up       117
topTab_fdr_GSE176385 <- topTable (efit_GSE176385, number=nrow(efit_GSE176385), coef="MMvsHD", adjust="fdr")
#head(topTab_fdr_GSE176385)
dim(topTab_fdr_GSE176385)
# 14039     6

## add Entrez IDs and save
topTab_fdr_GSE176385_annotated <- merge(as.data.frame(topTab_fdr_GSE176385),database_gene_symbols, by.x = 0, by.y = "id")
write.csv(topTab_fdr_GSE176385_annotated, file="./GSE176385_DEG_limma.csv")

## Create table of significant ones
topTab_fdr_sig_GSE176385 <- subset(topTab_fdr_GSE176385, topTab_fdr_GSE176385$adj.P.Val < 0.05 & (topTab_fdr_GSE176385$logFC < -1 | topTab_fdr_GSE176385$logFC > 1))
dim(topTab_fdr_sig_GSE176385)
# 277    6

## Add Entrez IDs and save
topTab_fdr_sig_GSE176385_annotated <- merge(as.data.frame(topTab_fdr_sig_GSE176385),database_gene_symbols, by.x = 0, by.y = "id")
write.csv(topTab_fdr_sig_GSE176385_annotated, file="./GSE176385_DEG_limma_sig.csv")
