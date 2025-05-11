## Workflow protocols
## https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html
## https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html

# Needed packages
library(biomaRt)
library(edgeR)
library(limma)

# setting working directory 
setwd("~/Work/PhD/Projects/MoD TrD LRI MM-Pain/Computational Analysis/Ligand_Receptor/Data/MM/PC/hDatasets/GSE175384")

# reading count files
counts_GSE175384 <- read.table("GSE175384_Counts.txt", header=T)

# number of rows
nrow(counts_GSE175384)

# removing non-relevant samples (AL, MGUS,Plasmablasts)
counts_GSE175384[,c("RNAseq1","RNAseq2","RNAseq3","RNAseq4","RNAseq5","RNAseq6","RNAseq7","RNAseq8","RNAseq9","RNAseq10","RNAseq11","RNAseq12","RNAseq13","RNAseq14","RNAseq15","RNAseq16","RNAseq17","RNAseq18","RNAseq19","RNAseq20","RNAseq21","RNAseq22","RNAseq23","RNAseq24","RNAseq25","RNAseq26","RNAseq27","RNAseq28","RNAseq29","RNAseq30",
          "RNAseq31","RNAseq32","RNAseq65","RNAseq66","RNAseq67","RNAseq68","RNAseq69","RNAseq70","RNAseq71","RNAseq72","RNAseq73","RNAseq74","RNAseq75","RNAseq76","RNAseq77","RNAseq78","RNAseq79")] <- NULL

# read design table
design <- read.csv("design.csv", stringsAsFactors = FALSE)

# reorder count table based on design table 
counts_GSE175384 <- counts_GSE175384[,match(design$sample, colnames(counts_GSE175384))]

# stop if not(all.equal(colnames(counts), design$sample))
all.equal(colnames(counts_GSE175384), design$sample)

# print number of reads which could be aligned to the different genes
print(colSums(counts_GSE175384))

# saving files
id_GSE175384 <- row.names(counts_GSE175384)
counts_GSE175384_save <- cbind("id" = id_GSE175384, counts_GSE175384)
write.csv(counts_GSE175384_save, file="./Unfiltered.Unnormalized.Annotated.Data.csv",row.names = FALSE)

# getting ensembl_gene_id annotation
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
database_gene_symbols <- getBM(
  attributes = c("hgnc_symbol", "entrezgene_id"),
  filters = "hgnc_symbol", values = rownames(counts_GSE175384),
  mart = mart, useCache = F)
database_gene_symbols<- na.omit(database_gene_symbols)
colnames(database_gene_symbols)[1] <- "id"

# filtering by mapping ids
ids_to_keep_GSE175384 <- database_gene_symbols$id
counts_GSE175384_filtered <- counts_GSE175384[ids_to_keep_GSE175384, ]
counts_GSE175384_filtered[is.na(counts_GSE175384_filtered)] <- 0
nrow(counts_GSE175384_filtered)
# 19535

# saving files with annotation and filtering genes (ID)
counts_GSE175384_save <- merge(counts_GSE175384_save,database_gene_symbols,by = "id")
counts_GSE175384_save <- counts_GSE175384_save[colnames(counts_GSE175384_save)[c(66,1:65)]]
# counts_GSE175384_save$id <- NULL
write.csv(counts_GSE175384_save, file="./Filtered(ID).Unnormalized.Annotated.Data.csv",row.names = FALSE)

# defining groups
gsms_GSE175384 <- "1111111111110000000000000000000000000000000011111111111111111111"
sml_GSE175384 <- c()
conditions_GSE175384 <- c("HD","MM")
for (i in 1:nchar(gsms_GSE175384)) { sml_GSE175384[i] <- conditions_GSE175384[as.integer(substr(gsms_GSE175384,i,i))+1] }
fl_GSE175384 <- as.factor(sml_GSE175384)

# create DGEList object from the count matrix ??
merged_dge_GSE175384 <- DGEList(counts_GSE175384_filtered)
merged_dge_GSE175384$samples$group <- fl_GSE175384
print(merged_dge_GSE175384$samples)
dim(merged_dge_GSE175384)
# 19535    64

#  filter genes using edgeR, while keeping as many genes as possible with worthwhile counts 
table(rowSums(merged_dge_GSE175384$counts==0))

keep.exprs_GSE175384 <- filterByExpr(merged_dge_GSE175384, group=fl_GSE175384, min.count = 1)
merged_filtered_GSE175384 <- merged_dge_GSE175384[keep.exprs_GSE175384, keep.lib.sizes=FALSE]
dim(merged_filtered_GSE175384)
# 11904    64

table(rowSums(merged_filtered_GSE175384$counts==0))


# saving files with annotation and filtering genes ID and counts
id_GSE175384 <- row.names(merged_filtered_GSE175384)
counts_GSE175384_save <- cbind("id" = id_GSE175384, merged_filtered_GSE175384$counts)
counts_GSE175384_save <- merge(counts_GSE175384_save,database_gene_symbols,by = "id")
counts_GSE175384_save <- counts_GSE175384_save[colnames(counts_GSE175384_save)[c(66,1:65)]]
#counts_GSE175384_save$id <- NULL
write.csv(counts_GSE175384_save, file="./Filtered(ID,Counts).Unnormalized.Annotated.Data.csv",row.names = FALSE)

# saving NetworkAnalyst files (cancelled)
counts_GSE175384_save_NA <- rbind(class_GSE175384,counts_GSE175384_save)
colnames(counts_GSE175384_save_NA)[1] <- "#NAME"
write.table(counts_GSE175384_save_NA, file="./NA.Filtered(ID,Counts).Unnormalized.Annotated.Data.txt", row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)

# normalising gene expression distribution using edgeR
merged_norm_GSE175384 <- calcNormFactors(merged_filtered_GSE175384, method = "TMM")

# saving files with annotation and filtering genes (ID), counts and normalization
## https://support.bioconductor.org/p/103747/#103769
merged_norm_GSE175384_counts <- cpm(merged_norm_GSE175384)
counts_GSE175384_save <- cbind("id" = id_GSE175384, merged_norm_GSE175384_counts)
counts_GSE175384_save <- merge(counts_GSE175384_save,database_gene_symbols,by = "id")
counts_GSE175384_save <- counts_GSE175384_save[colnames(counts_GSE175384_save)[c(66,1:65)]]
#counts_GSE175384_save$id <- NULL
write.csv(counts_GSE175384_save, file="./Filtered(ID,Counts).Normalized.Annotated.Data_cpm.csv",row.names = FALSE)

# saving NetworkAnalyst files (cancelled)
counts_GSE175384_save_NA <- rbind(class_GSE175384,counts_GSE175384_save)
colnames(counts_GSE175384_save_NA)[1] <- "#NAME"
write.table(counts_GSE175384_save_NA, file="./NA.Filtered(ID,Counts).Normalized.Annotated.Data_cpm.txt", row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)

# [optional] multidimensional scaling (MDS) plot
par(mfrow=c(1,1))
lcpm_norm_GSE175384 <- cpm(merged_norm_GSE175384, log=TRUE)
plotMDS(lcpm_norm_GSE175384, labels=fl_GSE175384, col = as.numeric(fl_GSE175384))
title(main="Sample groups")


# [optional] filter based on variance 
filtered_var_GSE175384 <- varFilter(lcpm_norm_GSE175384, var.func=IQR, var.cutoff=0.5, filterByQuantile=TRUE)
rows_to_keep_var_GSE175384 <- rownames(filtered_var_GSE175384)
merged_norm_filtered_GSE175384 <- merged_norm_GSE175384[rows_to_keep_var_GSE175384,]
dim(merged_norm_filtered_GSE175384)
# 6012   64

# creating a design matrix and contrasts
designMat.id_GSE175384 <- model.matrix(~0+fl_GSE175384)
colnames(designMat.id_GSE175384) <- levels(fl_GSE175384)
rownames(designMat.id_GSE175384) <- colnames(merged_norm_GSE175384)
print(designMat.id_GSE175384)
cont.matrix.id_GSE175384 <- makeContrasts(MMvsHD = MM-HD,levels=colnames(designMat.id_GSE175384))
print(cont.matrix.id_GSE175384)

# remove heteroscedascity from count data
par(mfrow=c(1,2))
v_GSE175384 <- voom(merged_norm_GSE175384, designMat.id_GSE175384, plot=TRUE)
# counts in v_GSE175384$E
# When operating on a DGEList-object, voom converts raw counts to log-CPM values by automatically extracting library sizes and normalisation factors from x itself. 
# Additional normalisation to log-CPM values can be specified within voom using the normalize.method argument.

# saving files with annotation and filtering genes (ID), counts  and normalization using voom
## see https://support.bioconductor.org/p/62541/
merged_norm_GSE175384_voom_counts <- data.frame(2^(v_GSE175384$E)*v_GSE175384$targets$norm.factors)
counts_GSE175384_save <- cbind("id" = id_GSE175384, merged_norm_GSE175384_voom_counts)
counts_GSE175384_save <- merge(counts_GSE175384_save,database_gene_symbols,by = "id")
counts_GSE175384_save <- counts_GSE175384_save[colnames(counts_GSE175384_save)[c(66,1:65)]]
# counts_GSE175384_save$id <- NULL
write.csv(counts_GSE175384_save, file="./Filtered(ID,Counts).Normalized.Annotated.Data_voom.csv",row.names = FALSE)

## saving NetworkAnalyst files (cancelled)
counts_GSE175384_save_NA <- rbind(class_GSE175384,counts_GSE175384_save)
colnames(counts_GSE175384_save_NA)[1] <- "#NAME"
write.table(counts_GSE175384_save_NA, file="./NA.Filtered(ID,Counts).Normalized.Annotated.Data_voom.txt", row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)

## fitting linear models for comparisons of interest
vfit_GSE175384 <- lmFit(v_GSE175384, designMat.id_GSE175384)
vfit_GSE175384 <- contrasts.fit(vfit_GSE175384, contrasts=cont.matrix.id_GSE175384)
efit_GSE175384 <- eBayes(vfit_GSE175384)
plotSA(efit_GSE175384, main="Final model: Mean-variance trend")

# Linear modelling in limma is carried out using the lmFit and contrasts.fit functions originally written for application to microarrays. 
# The functions can be used for both microarray and RNA-seq data and fit a separate model to the expression values for each gene. Next, 
# empirical Bayes moderation is carried out by borrowing information across all the genes to obtain more precise estimates of gene-wise 
# variability (Smyth 2004). The model’s residual variances are plotted against average expression values in the next figure. It can be 
# seen from this plot that the variance is no longer dependent on the mean expression level.

## examining the number of DE genes
dt_fdr_GSE175384 = decideTests(efit_GSE175384, adjust.method	= "fdr", p.value=0.05)
summary(dt_fdr_GSE175384)
#MMvsHD
#Down     1782
#NotSig   7393
#Up       2795
dt_fdr_lfc_GSE175384 = decideTests(efit_GSE175384, adjust.method	= "fdr", p.value=0.05, lfc = 1)
summary(dt_fdr_lfc_GSE175384)
#MMvsHD
#Down     1066
#NotSig   9446
#Up       1458
topTab_fdr_GSE175384 <- topTable (efit_GSE175384, number=nrow(efit_GSE175384), coef="MMvsHD", adjust="fdr")
#head(topTab_fdr_GSE175384)
dim(topTab_fdr_GSE175384)
# 11970     6

## add Entrez IDs and save
topTab_fdr_GSE175384_annotated <- merge(as.data.frame(topTab_fdr_GSE175384),database_gene_symbols, by.x = 0, by.y = "id")
write.csv(topTab_fdr_GSE175384_annotated, file="./GSE175384_DEG_limma.csv")

## Create table of significant ones
topTab_fdr_sig_GSE175384 <- subset(topTab_fdr_GSE175384, topTab_fdr_GSE175384$adj.P.Val < 0.05 & (topTab_fdr_GSE175384$logFC < -1 | topTab_fdr_GSE175384$logFC > 1))
dim(topTab_fdr_sig_GSE175384)
# 2524    6

## Add Entrez IDs and save
topTab_fdr_sig_GSE175384_annotated <- merge(as.data.frame(topTab_fdr_sig_GSE175384),database_gene_symbols, by.x = 0, by.y = "id")
write.csv(topTab_fdr_sig_GSE175384_annotated, file="./GSE175384_DEG_limma_sig.csv")
