## Workflow protocols
## https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html
## https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html

# Needed packages
library(biomaRt)
library(edgeR)
library(limma)

# setting working directory 
setwd("W:/PhD/Projects/MM Pain/Computation MM TME/Cells/MDSC/Datasets MDSC/GSE150021")

# reading count files
counts_GSE150021 <- read.table("GSE150021_MDSC_RNAseq.txt", header=T)

# number of rows
nrow(counts_GSE150021)

# print number of reads which could be aligned to the different genes
print(colSums(counts_GSE150021))

# saving files
id_GSE150021 <- row.names(counts_GSE150021)
counts_GSE150021_save <- cbind("id" = id_GSE150021, counts_GSE150021)
write.csv(counts_GSE150021_save, file="./Unfiltered.Unnormalized.Annotated.Data.csv",row.names = FALSE)

# saving NetworkAnalyst files
class_GSE150021 <- data.frame("#CLASS","HD","HD","HD","HD","HD","MM","MM","MM","MM","HD","HD","HD","HD","HD"
                              ,"MM","MM","MM","MM","MM","HD","HD","HD","HD","HD","MM","MM","MM","MM","MM","MM","MM",
                              "HD","HD","HD","MM","MM","MM","MM","MM","MM","HD","HD","HD","HD","HD","HD",stringsAsFactors = FALSE) 
names(class_GSE150021) <- colnames(counts_GSE150021_save)
counts_GSE150021_save_NA <- rbind(class_GSE150021,counts_GSE150021_save)
colnames(counts_GSE150021_save_NA)[1] <- "#NAME"
write.table(counts_GSE150021_save_NA, file="./NA.UnFiltered.Unnormalized.Annotated.Data.txt", row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)

# getting ensembl_gene_id annotation
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
database_gene_symbols <- getBM(
  attributes = c("hgnc_symbol", "entrezgene_id"),
  filters = "hgnc_symbol", values = rownames(counts_GSE150021),
  mart = mart, useCache = F)
database_gene_symbols<- na.omit(database_gene_symbols)
colnames(database_gene_symbols)[1] <- "id"

# filtering by mapping ids
ids_to_keep_GSE150021 <- database_gene_symbols$id
counts_GSE150021_filtered <- counts_GSE150021[ids_to_keep_GSE150021, ]
counts_GSE150021_filtered[is.na(counts_GSE150021_filtered)] <- 0
nrow(counts_GSE150021_filtered)
# 21057

# saving files with annotation and filtering genes (ID)
counts_GSE150021_save <- merge(counts_GSE150021_save,database_gene_symbols,by = "id")
counts_GSE150021_save <- counts_GSE150021_save[colnames(counts_GSE150021_save)[c(48,1:47)]]
counts_GSE150021_save$id <- NULL
write.csv(counts_GSE150021_save, file="./Filtered(ID).Unnormalized.Annotated.Data.csv",row.names = FALSE)

# saving NetworkAnalyst files
names(class_GSE150021) <- colnames(counts_GSE150021_save)
counts_GSE150021_save_NA <- rbind(class_GSE150021,counts_GSE150021_save)
colnames(counts_GSE150021_save_NA)[1] <- "#NAME"
write.table(counts_GSE150021_save_NA, file="./NA.Filtered(ID).Unnormalized.Annotated.Data.txt", row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)

# defining groups
gsms_GSE150021 <- "0000011110000011111000001111111000111111000000"
sml_GSE150021 <- c()
conditions_GSE150021 <- c("HD","MM")
for (i in 1:nchar(gsms_GSE150021)) { sml_GSE150021[i] <- conditions_GSE150021[as.integer(substr(gsms_GSE150021,i,i))+1] }
fl_GSE150021 <- as.factor(sml_GSE150021)

# create DGEList object from the count matrix ??
merged_dge_GSE150021 <- DGEList(counts_GSE150021_filtered)
merged_dge_GSE150021$samples$group <- fl_GSE150021
print(merged_dge_GSE150021$samples)
dim(merged_dge_GSE150021)
# 21057    46

#  filter genes using edgeR, while keeping as many genes as possible with worthwhile counts 
table(rowSums(merged_dge_GSE150021$counts==0))

keep.exprs_GSE150021 <- filterByExpr(merged_dge_GSE150021, group=fl_GSE150021, min.count = 1)
merged_filtered_GSE150021 <- merged_dge_GSE150021[keep.exprs_GSE150021, keep.lib.sizes=FALSE]
dim(merged_filtered_GSE150021)

table(rowSums(merged_filtered_GSE150021$counts==0))


# saving files with annotation and filtering genes ID and counts
id_GSE150021 <- row.names(merged_filtered_GSE150021)
counts_GSE150021_save <- cbind("id" = id_GSE150021, merged_filtered_GSE150021$counts)
counts_GSE150021_save <- merge(counts_GSE150021_save,database_gene_symbols,by = "id")
counts_GSE150021_save <- counts_GSE150021_save[colnames(counts_GSE150021_save)[c(48,1:47)]]
counts_GSE150021_save$id <- NULL
write.csv(counts_GSE150021_save, file="./Filtered(ID,Counts).Unnormalized.Annotated.Data.csv",row.names = FALSE)

# saving NetworkAnalyst files
names(class_GSE150021) <- names(counts_GSE150021_save)
counts_GSE150021_save_NA <- rbind(class_GSE150021,counts_GSE150021_save)
colnames(counts_GSE150021_save_NA)[1] <- "#NAME"
write.table(counts_GSE150021_save_NA, file="./NA.Filtered(ID,Counts).Unnormalized.Annotated.Data.txt", row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)

# normalising gene expression distribution using edgeR
merged_norm_GSE150021 <- calcNormFactors(merged_filtered_GSE150021, method = "TMM")

# saving files with annotation and filtering genes (ID), counts and normalization
## https://support.bioconductor.org/p/103747/#103769
merged_norm_GSE150021_counts <- cpm(merged_norm_GSE150021)
counts_GSE150021_save <- cbind("id" = id_GSE150021, merged_norm_GSE150021_counts)
counts_GSE150021_save <- merge(counts_GSE150021_save,database_gene_symbols,by = "id")
counts_GSE150021_save <- counts_GSE150021_save[colnames(counts_GSE150021_save)[c(48,1:47)]]
counts_GSE150021_save$id <- NULL
write.csv(counts_GSE150021_save, file="./Filtered(ID,Counts).Normalized.Annotated.Data_cpm.csv",row.names = FALSE)

# saving NetworkAnalyst files
counts_GSE150021_save_NA <- rbind(class_GSE150021,counts_GSE150021_save)
colnames(counts_GSE150021_save_NA)[1] <- "#NAME"
write.table(counts_GSE150021_save_NA, file="./NA.Filtered(ID,Counts).Normalized.Annotated.Data_cpm.txt", row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)

# [optional] multidimensional scaling (MDS) plot
par(mfrow=c(1,1))
lcpm_norm_GSE150021 <- cpm(merged_norm_GSE150021, log=TRUE)
plotMDS(lcpm_norm_GSE150021, labels=fl_GSE150021, col = as.numeric(fl_GSE150021))
title(main="Sample groups")


# [optional] filter based on variance 
require(genefilter)
filtered_var_GSE150021 <- varFilter(lcpm_norm_GSE150021, var.func=IQR, var.cutoff=0.5, filterByQuantile=TRUE)
rows_to_keep_var_GSE150021 <- rownames(filtered_var_GSE150021)
merged_norm_filtered_GSE150021 <- merged_norm_GSE150021[rows_to_keep_var_GSE150021,]
dim(merged_norm_filtered_GSE150021)
# 5277   46

# creating a design matrix and contrasts
designMat.id_GSE150021 <- model.matrix(~0+fl_GSE150021)
colnames(designMat.id_GSE150021) <- levels(fl_GSE150021)
rownames(designMat.id_GSE150021) <- colnames(merged_norm_GSE150021)
print(designMat.id_GSE150021)
cont.matrix.id_GSE150021 <- makeContrasts(MMvsHD = MM-HD,levels=colnames(designMat.id_GSE150021))
print(cont.matrix.id_GSE150021)

# remove heteroscedascity from count data
par(mfrow=c(1,2))
v_GSE150021 <- voom(merged_norm_GSE150021, designMat.id_GSE150021, plot=TRUE)
# counts in v_GSE150021$E
# When operating on a DGEList-object, voom converts raw counts to log-CPM values by automatically extracting library sizes and normalisation factors from x itself. 
# Additional normalisation to log-CPM values can be specified within voom using the normalize.method argument.

# saving files with annotation and filtering genes (ID), counts  and normalization using voom
## see https://support.bioconductor.org/p/62541/
merged_norm_GSE150021_voom_counts <- data.frame(2^(v_GSE150021$E)*v_GSE150021$targets$norm.factors)
counts_GSE150021_save <- cbind("id" = id_GSE150021, merged_norm_GSE150021_voom_counts)
counts_GSE150021_save <- merge(counts_GSE150021_save,database_gene_symbols,by = "id")
counts_GSE150021_save <- counts_GSE150021_save[colnames(counts_GSE150021_save)[c(48,1:47)]]
counts_GSE150021_save$id <- NULL
write.csv(counts_GSE150021_save, file="./Filtered(ID,Counts).Normalized.Annotated.Data_voom.csv",row.names = FALSE)

## saving NetworkAnalyst files
counts_GSE150021_save_NA <- rbind(class_GSE150021,counts_GSE150021_save)
colnames(counts_GSE150021_save_NA)[1] <- "#NAME"
write.table(counts_GSE150021_save_NA, file="./NA.Filtered(ID,Counts).Normalized.Annotated.Data_voom.txt", row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)

## fitting linear models for comparisons of interest
## blocking factor method https://www.biostars.org/p/54565/
require(Biobase)
my.targets <-read.AnnotatedDataFrame(file.path("targets.csv"), 
                                     header = TRUE, row.names = 1, 
                                     sep=";") 
corfit <- duplicateCorrelation(v_GSE150021, designMat.id_GSE150021, block=my.targets$Stage)
vfit_GSE150021 <- lmFit(v_GSE150021, designMat.id_GSE150021, block=my.targets$Stage, correlation=corfit$consensus.correlation)
vfit_GSE150021 <- contrasts.fit(vfit_GSE150021, contrasts=cont.matrix.id_GSE150021)
efit_GSE150021 <- eBayes(vfit_GSE150021)
plotSA(efit_GSE150021, main="Final model: Mean-variance trend")

# Linear modelling in limma is carried out using the lmFit and contrasts.fit functions originally written for application to microarrays. 
# The functions can be used for both microarray and RNA-seq data and fit a separate model to the expression values for each gene. Next, 
# empirical Bayes moderation is carried out by borrowing information across all the genes to obtain more precise estimates of gene-wise 
# variability (Smyth 2004). The model’s residual variances are plotted against average expression values in the next figure. It can be 
# seen from this plot that the variance is no longer dependent on the mean expression level.

## examining the number of DE genes
dt_fdr_GSE150021 = decideTests(efit_GSE150021, adjust.method	= "fdr", p.value=0.05)
summary(dt_fdr_GSE150021)
#MMvsHD
#Down     162
#NotSig   10256
#Up       137
dt_fdr_lfc_GSE150021= decideTests(efit_GSE150021, adjust.method	= "fdr", p.value=0.05, lfc = 1)
summary(dt_fdr_lfc_GSE150021)
#MMvsHD
#Down     130
#NotSig   10349
#Up       76
topTab_fdr_GSE150021 <- topTable (efit_GSE150021, number=nrow(efit_GSE150021), coef="MMvsHD", adjust="fdr")
#head(topTab_fdr_GSE150021)
dim(topTab_fdr_GSE150021)
# 10555        6

## add Entrez IDs and save
topTab_fdr_GSE150021_annotated <- merge(as.data.frame(topTab_fdr_GSE150021),database_gene_symbols, by.x = 0, by.y = "id")
write.csv(topTab_fdr_GSE150021_annotated, file="./GSE150021_DEG_limma.csv")

## Create table of significant ones
topTab_fdr_sig_GSE150021 <- subset(topTab_fdr_GSE150021, topTab_fdr_GSE150021$adj.P.Val < 0.05 & (topTab_fdr_GSE150021$logFC < -1 | topTab_fdr_GSE150021$logFC > 1))
dim(topTab_fdr_sig_GSE150021)
# 206    6

## Add Entrez IDs and save
topTab_fdr_sig_GSE150021_annotated <- merge(as.data.frame(topTab_fdr_sig_GSE150021),database_gene_symbols, by.x = 0, by.y = "id")
write.csv(topTab_fdr_sig_GSE150021_annotated, file="./GSE175384_DEG_limma_sig.csv")
