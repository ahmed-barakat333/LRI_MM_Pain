# needed packages
library(edgeR)
library(RColorBrewer)
library(limma)

## Workflow based on 
### see https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html
### see also https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html

# setting working directoroy 
setwd("~/Work/PhD/Projects/MoD TrD LRI MM-Pain/Computational Analysis/Ligand_Receptor/Data/MM/PC/hDatasets/GSE76425")

## read data in
counts_1 <- read.csv("data/GSM2011243_BMPC49.txt", sep="\t", header = TRUE)
counts_1 <- counts_1[ , c("gid","id", "count")]
counts_1 <- counts_1[complete.cases(counts_1), ] 
colnames(counts_1)[3] <- "GSM2011243_BMPC49"

counts_2 <- read.csv("data/GSM2011244_BMPC56.txt", sep="\t", header = TRUE)
counts_2 <- counts_2[ , c("gid","id", "count")]
counts_2 <- counts_2[complete.cases(counts_2), ] 
colnames(counts_2)[3] <- "GSM2011244_BMPC56"

counts_3 <- read.csv("data/GSM2011245_BMPC57.txt", sep="\t", header = TRUE)
counts_3 <- counts_3[ , c("gid","id", "count")]
counts_3 <- counts_3[complete.cases(counts_3), ] 
colnames(counts_3)[3] <- "GSM2011245_BMPC57"

counts_4 <- read.csv("data/GSM2011246_BMPC65.txt", sep="\t", header = TRUE)
counts_4 <- counts_4[ , c("gid","id", "count")]
counts_4 <- counts_4[complete.cases(counts_4), ] 
colnames(counts_4)[3] <- "GSM2011246_BMPC65"

counts_5 <- read.csv("data/GSM2011247_PDLen1_A.txt", sep="\t", header = TRUE)
counts_5 <- counts_5[ , c("gid","id", "count")]
counts_5 <- counts_5[complete.cases(counts_5), ] 
colnames(counts_5)[3] <- "GSM2011247_PDLen1_A"

counts_6 <- read.csv("data/GSM2011252_PDLen2_A.txt", sep="\t", header = TRUE)
counts_6 <- counts_6[ , c("gid","id", "count")]
counts_6 <- counts_6[complete.cases(counts_6), ] 
colnames(counts_6)[3] <- "GSM2011252_PDLen2_A"

counts_7 <- read.csv("data/GSM2011257_PDLen3_A.txt", sep="\t", header = TRUE)
counts_7 <- counts_7[ , c("gid","id", "count")]
counts_7 <- counts_7[complete.cases(counts_7), ] 
colnames(counts_7)[3] <- "GSM2011257_PDLen3_A"

counts_8 <- read.csv("data/GSM2011262_PDLen4_A.txt", sep="\t", header = TRUE)
counts_8 <- counts_8[ , c("gid","id", "count")]
counts_8 <- counts_8[complete.cases(counts_8), ] 
colnames(counts_8)[3] <- "GSM2011262_PDLen4_A"

counts_9 <- read.csv("data/GSM2011267_PDLen5_A.txt", sep="\t", header = TRUE)
counts_9 <- counts_9[ , c("gid","id", "count")]
counts_9 <- counts_9[complete.cases(counts_9), ] 
colnames(counts_9)[3] <- "GSM2011267_PDLen5_A"

counts_10 <- read.csv("data/GSM2011271_PDLen6_A.txt", sep="\t", header = TRUE)
counts_10 <- counts_10[ , c("gid","id", "count")]
counts_10 <- counts_10[complete.cases(counts_10), ] 
colnames(counts_10)[3] <- "GSM2011271_PDLen6_A"

counts_11 <- read.csv("data/GSM2011276_PDLen7_A.txt", sep="\t", header = TRUE)
counts_11 <- counts_11[ , c("gid","id", "count")]
counts_11 <- counts_11[complete.cases(counts_11), ] 
colnames(counts_11)[3] <- "GSM2011276_PDLen7_A"

counts_GSE76425 <- Reduce(function(...) merge(..., by=c ("id","gid"), all=TRUE), list(counts_1, counts_2, counts_3,counts_4,counts_5,counts_6,counts_7,counts_8,counts_9,counts_10,counts_11))
remove(counts_1,counts_2,counts_3,counts_4,counts_5,counts_6,counts_7,counts_8,counts_9,counts_10,counts_11)

# remove genes without symbol

# take gene entrez id with max count
require(dplyr)
counts_GSE76425 <- data.frame(counts_GSE76425 %>% 
                              group_by(gid) %>% 
                              summarise_all(funs(max)))


## saving files
write.csv(counts_GSE76425, file="results/Unfiltered.Unnormalized.Unannotated.Data.csv",row.names = FALSE)

## create dataframe with entrez id and gene symbol
symbol_GSE76425 <- counts_GSE76425[,c(1,2)]
colnames(symbol_GSE76425) <- c("id","symbol")

## move gene names/ids from the first column to the row names
genes_GSE76425 <- counts_GSE76425[,1]
counts_GSE76425 <- counts_GSE76425[,-c(1,2)]
row.names(counts_GSE76425) <- genes_GSE76425

counts_GSE76425[is.na(counts_GSE76425)] <- 0
colSums(counts_GSE76425)
# read depths range from ~25M - ~54M

## define groups
gsms_GSE76425 <- "00001111111"
sml_GSE76425 <- c()
conditions_GSE76425 <- c("HD","MM")
for (i in 1:nchar(gsms_GSE76425)) { sml_GSE76425[i] <- conditions_GSE76425[as.integer(substr(gsms_GSE76425,i,i))+1] }
fl_GSE76425 <- as.factor(sml_GSE76425)

## create DGEList object from the count matrix
merged_dge_GSE76425 <- DGEList(counts_GSE76425)
merged_dge_GSE76425$samples$group <- fl_GSE76425
print(merged_dge_GSE76425$samples)
dim(merged_dge_GSE76425)
#[1] 23495    11

## is there a need to filter lowly expressed genes? 
## output number of genes with 0 counts for a number samples
table(rowSums(merged_dge_GSE76425$counts==0))
## yes, there is a need to filter lowly expressed genes! 1100 genes have 0 counts in 10 samples

## automatic way to filter genes using edgeR, while keeping as many genes as possible with worthwhile counts 
keep.exprs_GSE76425 <- filterByExpr(merged_dge_GSE76425, group=fl_GSE76425, min.count = 1)
merged_filtered_GSE76425 <- merged_dge_GSE76425[keep.exprs_GSE76425, keep.lib.sizes=FALSE]
dim(merged_filtered_GSE76425)
#[1] 20437    11
table(rowSums(merged_filtered_GSE76425$counts==0))
# -> less samples with in which a gene is 0

## saving files with annotation and filtering genes expression
id_GSE76425 <- row.names(merged_filtered_GSE76425)
counts_GSE76425_save <- cbind("id" = id_GSE76425, merged_filtered_GSE76425$counts)
counts_GSE76425_save <- merge(symbol_GSE76425,counts_GSE76425_save, by="id")
write.csv(counts_GSE76425_save, file="results/Filtered(Counts).Unnormalized.Annotated.Data.csv",row.names = FALSE)

## saving NetworkAnalyst files (cancelled)
names(class_GSE76425) <- colnames(counts_GSE76425_save)
counts_GSE76425_save_NA <- rbind(class_GSE76425,counts_GSE76425_save)
colnames(counts_GSE76425_save_NA)[1] <- "#NAME"
write.table(counts_GSE76425_save_NA, file="results/NA.Filtered(Counts).Unnormalized.Annotated.Data.txt", row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)

## normalising gene expression distribution
merged_norm_GSE76425 <- calcNormFactors(merged_filtered_GSE76425, method = "TMM")
#print(merged_norm_GSE76425$samples)

## saving files with annotation and filtering genes expression and normalization
## see https://support.bioconductor.org/p/103747/#103769
merged_norm_GSE76425_counts <- cpm(merged_norm_GSE76425)
counts_GSE76425_save <- cbind("id" = id_GSE76425, merged_norm_GSE76425_counts)
counts_GSE76425_save <- merge(symbol_GSE76425,counts_GSE76425_save, by="id")
write.csv(counts_GSE76425_save, file="results/Filtered(Counts).Normalized.Annotated.Data_cpm.csv",row.names = FALSE)

## saving NetworkAnalyst files (cancelled)
counts_GSE76425_save_NA <- rbind(class_GSE76425,counts_GSE76425_save)
colnames(counts_GSE76425_save_NA)[1] <- "#NAME"
write.table(counts_GSE76425_save_NA, file="results/NA.Filtered(Counts).Normalized.Annotated.Data_cpm.txt", row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)

# [optional] multidimensional scaling (MDS) plot
par(mfrow=c(1,1))
lcpm_norm_GSE76425 <- cpm(merged_norm_GSE76425, log=TRUE)
plotMDS(lcpm_norm_GSE76425, labels=fl_GSE76425, col = as.numeric(fl_GSE76425))
title(main="Sample groups")
# [optional] 

# [optional] filter based on variance
filtered_var_GSE76425 <- varFilter(lcpm_norm_GSE76425, var.func=IQR, var.cutoff=0.5, filterByQuantile=TRUE)
rows_to_keep_var_GSE76425 <- rownames(filtered_var_GSE76425)
merged_norm_filtered_GSE76425 <- merged_norm_GSE76425[rows_to_keep_var_GSE76425,]
dim(merged_norm_filtered_GSE76425)
# 10220    11

## creating a design matrix and contrasts
designMat.id_GSE76425 <- model.matrix(~0+fl_GSE76425)
colnames(designMat.id_GSE76425) <- levels(fl_GSE76425)
rownames(designMat.id_GSE76425) <- colnames(merged_norm_GSE76425)
print(designMat.id_GSE76425)
cont.matrix.id_GSE76425 <- makeContrasts(MMvsHD = MM-HD,levels=colnames(designMat.id_GSE76425))
print(cont.matrix.id_GSE76425)

## remove heteroscedascity from count data
par(mfrow=c(1,2))
v_GSE76425 <- voom(merged_norm_GSE76425, designMat.id_GSE76425, plot=TRUE)
# When operating on a DGEList-object, voom converts raw counts to log-CPM values by automatically extracting library sizes and normalisation factors from x itself. 
# Additional normalisation to log-CPM values can be specified within voom using the normalize.method argument.

## saving files with annotation and filtering genes (ID) and expression and normalization
## see https://support.bioconductor.org/p/62541/
merged_norm_GSE76425_voom_counts <- data.frame(2^(v_GSE76425$E)*v_GSE76425$targets$norm.factors)
counts_GSE76425_save <- cbind("id" = id_GSE76425, merged_norm_GSE76425_voom_counts)
counts_GSE76425_save <- merge(symbol_GSE76425,counts_GSE76425_save, by="id")
write.csv(counts_GSE76425_save, file="results/Filtered(Counts).Normalized.Annotated.Data_voom.csv",row.names = FALSE)

## saving NetworkAnalyst files (cancelled)
counts_GSE76425_save_NA <- rbind(class_GSE76425,counts_GSE76425_save)
colnames(counts_GSE76425_save_NA)[1] <- "#NAME"
write.table(counts_GSE76425_save_NA, file="results/NA.Filtered(Counts).Normalized.Annotated.Data_voom.txt", row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)

## fitting linear models for comparisons of interest
vfit_GSE76425 <- lmFit(v_GSE76425, designMat.id_GSE76425)
vfit_GSE76425 <- contrasts.fit(vfit_GSE76425, contrasts=cont.matrix.id_GSE76425)
efit_GSE76425 <- eBayes(vfit_GSE76425)
plotSA(efit_GSE76425, main="Final model: Mean-variance trend")
# Linear modelling in limma is carried out using the lmFit and contrasts.fit functions originally written for application to microarrays. 
# The functions can be used for both microarray and RNA-seq data and fit a separate model to the expression values for each gene. Next, 
# empirical Bayes moderation is carried out by borrowing information across all the genes to obtain more precise estimates of gene-wise 
# variability (Smyth 2004). The model’s residual variances are plotted against average expression values in the next figure. It can be 
# seen from this plot that the variance is no longer dependent on the mean expression level.

## examining the number of DE genes
dt_fdr_GSE76425 = decideTests(efit_GSE76425, adjust.method	= "fdr", p.value=0.05)
summary(dt_fdr_GSE76425)
#MMvsHD
#Down     5224
#NotSig  10988
#Up       4225
dt_fdr_lfc_GSE76425 = decideTests(efit_GSE76425, adjust.method	= "fdr", p.value=0.05, lfc = 1)
summary(dt_fdr_lfc_GSE76425)
#Down     4495
#NotSig  12662
#Up       3280
topTab_fdr_GSE76425 <- topTable (efit_GSE76425, number=nrow(efit_GSE76425), coef="MMvsHD", adjust="fdr")
dim(topTab_fdr_GSE76425)
# 20437     6

## save
write.csv(topTab_fdr_GSE76425, file="results/GSE76425_DEG_limma.csv")

## Create table of significant ones
topTab_fdr_sig_GSE76425 <- subset(topTab_fdr_GSE76425, topTab_fdr_GSE76425$adj.P.Val < 0.05 & (topTab_fdr_GSE76425$logFC < -1 | topTab_fdr_GSE76425$logFC > 1))
dim(topTab_fdr_sig_GSE76425)
# 7775    6

## Save
write.csv(topTab_fdr_sig_GSE76425, file="results/GSE76425_DEG_limma_sig.csv")
