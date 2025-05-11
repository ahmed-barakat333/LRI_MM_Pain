# setting working directoroy 
setwd("~/Work/PhD/Projects/MM Pain/Computational Analysis/Ligand_Receptor/Data/MM/WBM/hDatasets/GSE118985")

# load series and platform data from GEO
require(GEOquery)
gset <- getGEO("GSE118985", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

samples_to_keep <- c("TT","Normal")
gset <- gset[,gset@phenoData@data[["patient group:ch1"]] %in% samples_to_keep]

require(Biobase)
my.targets <-read.AnnotatedDataFrame(file.path("data/targets_2.csv"), 
                                     header = TRUE, row.names = 1, 
                                     sep=";") 

# check that sample annotations in gset are the same as my.targets file
hd_gset <- gset[,gset@phenoData@data[["patient group:ch1"]] %in% "Normal"]
hd_gset <- sampleNames(hd_gset)

hd_targets <- read.csv("data/targets_2.csv", sep = ";")
hd_targets <- hd_targets[hd_targets$Group == 'HD',]
hd_targets <- hd_targets$FileName

sum(hd_gset %in% hd_targets)

mm_gset <- gset[,gset@phenoData@data[["patient group:ch1"]] %in% "TT"]
mm_gset <- sampleNames(mm_gset)

mm_targets <- read.csv("data/targets_2.csv", sep = ";")
mm_targets <- mm_targets[mm_targets$Group == 'MM',]
mm_targets <- mm_targets$FileName

sum(mm_gset %in% mm_targets)

# plots of normalized data 
require(ggplot2)
require(ggrepel)
plotPCA3 <- function (datos, labels, factor, title, scale,colores, size = 1.5, glineas = 0.25) {
  data <- prcomp(t(datos),scale=scale)
  # plot adjustments
  dataDf <- data.frame(data$x)
  Group <- factor
  loads <- round(data$sdev^2/sum(data$sdev^2)*100,1)
  # main plot
  p1 <- ggplot(dataDf,aes(x=PC1, y=PC2)) +
    theme_classic() +
    geom_hline(yintercept = 0, color = "gray70") +
    geom_vline(xintercept = 0, color = "gray70") +
    geom_point(aes(color = Group), alpha = 0.55, size = 3) +
    coord_cartesian(xlim = c(min(data$x[,1])-5,max(data$x[,1])+5)) +
    scale_fill_discrete(name = "Group")
  
  # avoiding labels superposition
  p1 + geom_text_repel(aes(y = PC2 + 0.25, label = labels),segment.size = 0.25, size = size) + 
    labs(x = c(paste("PC1",loads[1],"%")),y=c(paste("PC2",loads[2],"%"))) +  
    ggtitle(paste("Principal Component Analysis for: ",title,sep=" "))+ 
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values=colores)
}  
plotPCA3(exprs(gset), labels = my.targets$ShortName, factor = my.targets$Group, 
         title="Normalized data", scale = FALSE, size = 3, 
         colores = c("blue", "red"))

boxplot(exprs(gset), cex.axis=0.5, las=2,  which="all",
        main="Distribution of normalized intensity values")

# saving files
write.csv(exprs(gset), file="./results/Normalized.Unfiltered.Unannotated.Data.csv")

# probe id annotation
annot <- fData(gset)
annotated<-function(topTab, anotTable)
{
  topTab <- cbind(PROBEID=rownames(topTab), topTab)
  myProbes <- rownames(topTab)
  geneAnots <- anotTable[, c("ID", "ENTREZ_GENE_ID","Gene Symbol")]
  annotatedTopTab<- merge(x=geneAnots, y=topTab, by.x="ID", by.y="PROBEID")
  return(annotatedTopTab)
}
eset_rma_unfiltered_annotated <- annotated(exprs(gset),anotTable=annot)


# saving files
write.csv(eset_rma_unfiltered_annotated, file="./results/Normalized.Unfiltered.Annotated.Data.csv",row.names = FALSE)

# NetworkAnalyst formatting  (cancelled)
eset_rma_unfiltered_annotated[,c("ID","GENE_SYMBOL")] <- NULL
class <- data.frame("#CLASS","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","HD","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","HD","MM","MM","MM","MM","MM","MM","MM","HD") 
names(class) <- colnames(eset_rma_unfiltered_annotated)
eset_rma_unfiltered_annotated <- rbind(class,eset_rma_unfiltered_annotated)
colnames(eset_rma_unfiltered_annotated) <- c("#NAME","GSM3569443","GSM3569444","GSM3569445","GSM3569446","GSM3569447","GSM3569448","GSM3569449","GSM3569450","GSM3569451","GSM3569452","GSM3569453","GSM3569454","GSM3569455","GSM3569456","GSM3569457","GSM3569458","GSM3569459","GSM3569460","GSM3569461","GSM3569462","GSM3569463","GSM3569464","GSM3569465","GSM3569466","GSM3569467","GSM3569468","GSM3569469","GSM3569470","GSM3569471","GSM3569472","GSM3569473","GSM3569474","GSM3569475","GSM3569476","GSM3569477","GSM3569478","GSM3569479","GSM3569480","GSM3569481","GSM3569482","GSM3569483","GSM3569484","GSM3569485","GSM3569486","GSM3569487","GSM3569488","GSM3569489","GSM3569490")

# saving files
write.table(eset_rma_unfiltered_annotated, file="./results/NA.Normalized.Unfiltered.Annotated.Data.txt", row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)

# detecting variable genes
sds <- apply (exprs(gset), 1, sd)
sdsO<- sort(sds)
plot(1:length(sdsO), sdsO, main="Distribution of variability for all genes",
     sub="Vertical lines represent 90% and 95% percentiles",
     xlab="Gene index (from least to most variable)", ylab="Standard deviation")
abline(v=length(sds)*c(0.9,0.95))

# filtering genes (ID)
require(genefilter)

# filter based on ID available -> remove all annotations without an ENTREZ ID 
annot_filtered <- annot[!is.na(annot$ENTREZ_GENE_ID), ]
rows_to_keep <- rownames(annot_filtered)
gset_filtered <- gset[rows_to_keep,]

# filter duplicated entries based on variance
eset <- exprs(gset_filtered)
eset_variance <- rowQ(eset, ceiling(0.75 * ncol(eset))) - rowQ(eset, floor(0.25 * ncol(eset)))
entries_order_var <- order(eset_variance, decreasing=TRUE)
dup_gb_list_var <- duplicated(annot_filtered$ENTREZ_GENE_ID[entries_order_var])
gset_filtered_unique <- gset_filtered[entries_order_var,][!dup_gb_list_var, ]

# saving files
write.csv(exprs(gset_filtered_unique), file="./results/Normalized.Filtered(ID).Unannotated.Data.csv")

# probe id annotation
eset_rma_filtered.id_annotated<-annotated(exprs(gset_filtered_unique),anotTable=annot)

# saving files
eset_rma_filtered.id_annotated[5721,3]
eset_rma_filtered.id_annotated[5721,3] <- gsub("AGO2 /// CASC7 /// CASC7","AGO2 /// CASC7",eset_rma_unfiltered_annotated[5721,3] )
library(tidyverse)
eset_rma_filtered.id_annotated <- eset_rma_filtered.id_annotated %>% 
  separate_rows(ENTREZ_GENE_ID, `Gene Symbol`, sep = "///") 

write.csv(eset_rma_filtered.id_annotated, file="./results/Normalized.Filtered(ID).Annotated.Data.csv",row.names = FALSE)

# NetworkAnalyst formatting  (cancelled)
eset_rma_filtered.id_annotated[,c("ID","GENE_SYMBOL")] <- NULL
class <- data.frame("#CLASS","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","HD","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","HD","MM","MM","MM","MM","MM","MM","MM","HD") 
names(class) <- colnames(eset_rma_filtered.id_annotated)
eset_rma_filtered.id_annotated <- rbind(class,eset_rma_filtered.id_annotated)
colnames(eset_rma_filtered.id_annotated) <- c("#NAME","GSM3569443","GSM3569444","GSM3569445","GSM3569446","GSM3569447","GSM3569448","GSM3569449","GSM3569450","GSM3569451","GSM3569452","GSM3569453","GSM3569454","GSM3569455","GSM3569456","GSM3569457","GSM3569458","GSM3569459","GSM3569460","GSM3569461","GSM3569462","GSM3569463","GSM3569464","GSM3569465","GSM3569466","GSM3569467","GSM3569468","GSM3569469","GSM3569470","GSM3569471","GSM3569472","GSM3569473","GSM3569474","GSM3569475","GSM3569476","GSM3569477","GSM3569478","GSM3569479","GSM3569480","GSM3569481","GSM3569482","GSM3569483","GSM3569484","GSM3569485","GSM3569486","GSM3569487","GSM3569488","GSM3569489","GSM3569490")

# saving files
write.table(eset_rma_filtered.id_annotated, file="./results/NA.Normalized.Filtered(ID).Annotated.Data.txt", row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)

## differential gene experssion
# design matrix
require(limma)
designMat.id<- model.matrix(~0+my.targets$Group, pData(gset_filtered_unique))
colnames(designMat.id) <- c("HD","MM")
print(designMat.id)

# contrast matrix
cont.matrix.id <- makeContrasts (MMvsHD = MM-HD,levels=designMat.id)
print(cont.matrix.id)

# linear model fit
require(limma)
fit.id<-lmFit(gset_filtered_unique, designMat.id)
fit.main.id<-contrasts.fit(fit.id, cont.matrix.id)
fit.main.id<-eBayes(fit.main.id)
class(fit.main.id)

# top Tabs
topTab_MMvsHD.id <- topTable (fit.main.id, number=nrow(fit.main.id), coef="MMvsHD", adjust="fdr" ) 
head(topTab_MMvsHD.id)

# probe id annotation
annotatedTopTable <- function(topTab, anotTable)
{
  topTab <- cbind(PROBEID=rownames(topTab), topTab)
  myProbes <- rownames(topTab)
  geneAnots <- anotTable[, c("ID", "ENTREZ_GENE_ID","Gene Symbol")]
  annotatedTopTab<- merge(x=geneAnots, y=topTab, by.x="ID", by.y="PROBEID")
  return(annotatedTopTab)
}

topAnnotated_MMvsHD.id <- annotatedTopTable(topTab_MMvsHD.id, anotTable=annot)

write.csv(topAnnotated_MMvsHD.id, file="./results/DEG.ID.csv")

# filtering genes (variance/ ID)
filtered_var <- varFilter(exprs(gset_filtered_unique), var.func=IQR, var.cutoff=0.75, filterByQuantile=TRUE)
rows_to_keep_var <- rownames(filtered_var)
gset_filtered_unique_var <- gset[rows_to_keep_var,]

# plots of normalized filtered data
plotPCA3(exprs(gset_filtered_unique_var), labels = my.targets$ShortName, factor = my.targets$Group, 
         title="Normalized filtered data", scale = FALSE, size = 3, 
         colores = c("blue", "red"))

# saving files
write.csv(exprs(gset_filtered_unique_var), file="./results/Normalized.Filtered(ID,Var).Unannotated.Data.csv")

# probe id annotation
eset_rma_filtered.id.var_annotated<-annotated(exprs(gset_filtered_unique_var),anotTable=annot)

# saving files
write.csv(eset_rma_filtered.id.var_annotated, file="./results/Normalized.Filtered(ID,Var).Annotated.Data.csv",row.names = FALSE)

# NetworkAnalyst formatting  (cancelled)
eset_rma_filtered.id.var_annotated[,c("ID","GENE_SYMBOL")] <- NULL
class <- data.frame("#CLASS","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","HD","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","HD","MM","MM","MM","MM","MM","MM","MM","HD") 
names(class) <- colnames(eset_rma_filtered.id.var_annotated)
eset_rma_filtered.id.var_annotated <- rbind(class,eset_rma_filtered.id.var_annotated)
colnames(eset_rma_filtered.id.var_annotated) <- c("#NAME","GSM3569443","GSM3569444","GSM3569445","GSM3569446","GSM3569447","GSM3569448","GSM3569449","GSM3569450","GSM3569451","GSM3569452","GSM3569453","GSM3569454","GSM3569455","GSM3569456","GSM3569457","GSM3569458","GSM3569459","GSM3569460","GSM3569461","GSM3569462","GSM3569463","GSM3569464","GSM3569465","GSM3569466","GSM3569467","GSM3569468","GSM3569469","GSM3569470","GSM3569471","GSM3569472","GSM3569473","GSM3569474","GSM3569475","GSM3569476","GSM3569477","GSM3569478","GSM3569479","GSM3569480","GSM3569481","GSM3569482","GSM3569483","GSM3569484","GSM3569485","GSM3569486","GSM3569487","GSM3569488","GSM3569489","GSM3569490")

# saving files
write.table(eset_rma_filtered.id.var_annotated, file="./results/NA.Normalized.Filtered(ID,Var).Annotated.Data.txt", row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)

## differential gene experssion
# design matrix
require(limma)
designMat.id.var<- model.matrix(~0+my.targets$Group, pData(gset_filtered_unique_var))
colnames(designMat.id.var) <- c("HD","MM")
print(designMat.id.var)

# contrast matrix
cont.matrix.id.var <- makeContrasts (MMvsHD = MM-HD,levels=designMat.id.var)
print(cont.matrix.id.var)

# linear model fit
require(limma)
fit.id.var<-lmFit(gset_filtered_unique_var, designMat.id.var)
fit.main.id.var<-contrasts.fit(fit.id.var, cont.matrix.id.var)
fit.main.id.var<-eBayes(fit.main.id.var)
class(fit.main.id.var)

# top Tabs
topTab_MMvsHD.id.var <- topTable (fit.main.id.var, number=nrow(fit.main.id.var), coef="MMvsHD", adjust="fdr" ) 
head(topTab_MMvsHD.id.var)

# probe id annotation
topAnnotated_MMvsHD.id.var <- annotatedTopTable(topTab_MMvsHD.id.var,anotTable=annot)

write.csv(topAnnotated_MMvsHD.id.var, file="./results/DEG.ID.Var.csv")

#____________________________________________________________________________________________
#____________________________________________________________________________________________
## solve multiple ids 

library(tidyverse)
library(org.Hs.eg.db)
library(stringr)

DEG.ID <- read.csv("results/GSE118985 DEG.ID.csv")

sum(str_count(DEG.ID$ENTREZ_GENE_ID.x, '///') != str_count(DEG.ID$Gene.Symbol.1, '///'))
d <- which(str_count(DEG.ID$ENTREZ_GENE_ID.x, '///') != str_count(DEG.ID$Gene.Symbol.1, '///'))

#____________________________________________________________________________________________

DEG.ID <- DEG.ID %>% 
  separate_rows(ENTREZ_GENE_ID.x, sep = "///") %>% 
  mutate(across(where(is.character), str_trim))

mappings <- select(org.Hs.eg.db, DEG.ID$ENTREZ_GENE_ID.x,  'SYMBOL', 'ENTREZID')
colnames(mappings)[1] <- "ENTREZ_GENE_ID.x"
d <- merge(DEG.ID, mappings, by = "ENTREZ_GENE_ID.x" )

#____________________________________________________________________________________________

DEG.ID.SYMBOL <- DEG.ID %>% 
  separate_rows(Gene.Symbol.1, sep = "///") %>% 
  mutate(across(where(is.character), str_trim))

write.csv(DEG.ID.SYMBOL, file="./results/DEG.ID.SYMBOL.csv",row.names = FALSE)

DEG.ID.ENTREZID <- DEG.ID %>% 
  separate_rows(ENTREZ_GENE_ID.x, sep = "///") %>% 
  mutate(across(where(is.character), str_trim))

write.csv(DEG.ID.ENTREZID, file="./results/DEG.ID.ENTREZID.csv",row.names = FALSE)
