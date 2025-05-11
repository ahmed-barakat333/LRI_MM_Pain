# setting working directoroy 
setwd("~/Work/PhD/Projects/MM Pain/Computational Analysis/Ligand_Receptor/Data/MM/ADCY/mouse_Datasets/GSE143269")

# reading .CEL files
require(oligo)
celFiles <- list.celfiles("./data", full.names = TRUE)

require(Biobase)
my.targets <-read.AnnotatedDataFrame(file.path("./data","targets.csv"), 
                                     header = TRUE, row.names = 1, 
                                     sep=";") 

rawData<-read.celfiles(celFiles, phenoData = my.targets)

# quality control of raw data 
require(arrayQualityMetrics)
arrayQualityMetrics(rawData, outdir = file.path("./results",
                                                "QCDir.Raw"), force=TRUE)
# plots of raw data 
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
plotPCA3(exprs(rawData), labels = my.targets$ShortName, factor = my.targets$Group, 
         title="Raw data", scale = FALSE, size = 3, 
         colores = c("blue","red"))

boxplot(rawData, cex.axis=0.5, las=2,  which="all",
        col = c(rep("blue", 4), rep("red", 3)),
        main="Distribution of raw intensity values")

# data normalization
eset_rma<-rma(rawData)        

# quality control of normalized data
arrayQualityMetrics(eset_rma, outdir = file.path("./results",
                                                 "QCDir.Norm"), force=TRUE)
# plots of normalized data
plotPCA3(exprs(eset_rma), labels = my.targets$ShortName, factor = my.targets$Group, 
         title="Normalized data", scale = FALSE, size = 3, 
         colores = c("blue","red"))

boxplot(eset_rma, cex.axis=0.5, las=2,  which="all" ,
        col = c(rep("blue", 4), rep("red", 3)),
        main="Boxplot for arrays intensity: Normalized Data")

# saving files
write.csv(exprs(eset_rma), file="./results/Normalized.Unfiltered.Unannotated.Data.csv")

# probe id annotation
require(clariomsmousehttranscriptcluster.db)
annotation(eset_rma) <- "clariomsmousehttranscriptcluster.db"
annotated<-function(topTab, anotPackage)
{
  topTab <- cbind(PROBEID=rownames(topTab), topTab)
  myProbes <- rownames(topTab)
  thePackage <- eval(parse(text = anotPackage))
  geneAnots <- select(thePackage, myProbes, c("ENTREZID","SYMBOL"))
  annotatedTopTab<- merge(x=geneAnots, y=topTab, by.x="PROBEID", by.y="PROBEID")
  return(annotatedTopTab)
}
eset_rma_unfiltered_annotated <- annotated(exprs(eset_rma),anotPackage="clariomsmousehttranscriptcluster.db")

# saving files
write.csv(eset_rma_unfiltered_annotated, file="./results/Normalized.Unfiltered.Annotated.Data.csv",row.names = FALSE)

# NetworkAnalyst formatting (cancelled)
eset_rma_unfiltered_annotated[,c("PROBEID","SYMBOL")] <- NULL
class <- data.frame("#CLASS", "HD", "MM", "HD", "MM", "MM", "HD", "HD") 
colnames(class)<-colnames(eset_rma_unfiltered_annotated)
eset_rma_unfiltered_annotated <- rbind(class,eset_rma_unfiltered_annotated)
colnames(eset_rma_unfiltered_annotated)[1] <- "#NAME"

# saving files
write.table(eset_rma_unfiltered_annotated, file="./results/NA.Normalized.Unfiltered.Annotated.Data.txt", row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)

# detecting variable genes
sds <- apply (exprs(eset_rma), 1, sd)
sdsO<- sort(sds)
plot(1:length(sdsO), sdsO, main="Distribution of variability for all genes",
     sub="Vertical lines represent 90% and 95% percentiles",
     xlab="Gene index (from least to most variable)", ylab="Standard deviation")
abline(v=length(sds)*c(0.9,0.95))

# filtering genes (ID)
require(genefilter)
filtered.id <- nsFilter(eset_rma, 
                     require.entrez = TRUE, remove.dupEntrez = TRUE,
                     var.filter=FALSE, var.func=IQR, var.cutoff=0.75, 
                     filterByQuantile=TRUE, feature.exclude = "^AFFX")
names(filtered.id)
class(filtered.id$eset)
print(filtered.id$filter.log)
eset_filtered.id <-filtered.id$eset

# saving files
write.csv(exprs(eset_filtered.id), file="./results/Normalized.Filtered(ID).Unannotated.Data.csv")

# probe id annotation
eset_rma_filtered.id_annotated<-annotated(exprs(eset_filtered.id),anotPackage="clariomsmousehttranscriptcluster.db")

# saving files
write.csv(eset_rma_filtered.id_annotated, file="./results/Normalized.Filtered(ID).Annotated.Data.csv",row.names = FALSE)

# NetworkAnalyst formatting (cancelled)
eset_rma_filtered.id_annotated[,c("PROBEID","SYMBOL")] <- NULL
class <- data.frame("#CLASS", "MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","HD",	"HD",	"HD") 
colnames(class)<-colnames(eset_rma_filtered.id_annotated)
eset_rma_filtered.id_annotated <- rbind(class,eset_rma_filtered.id_annotated)
colnames(eset_rma_filtered.id_annotated)[1] <- "#NAME"

# saving files
write.table(eset_rma_filtered.id_annotated, file="./results/NA.Normalized.Filtered(ID).Annotated.Data.txt", row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)

## differential gene experssion
# design matrix
require(limma)
designMat.id<- model.matrix(~0+Group, pData(eset_filtered.id))
colnames(designMat.id) <- c("HD","MM")
print(designMat.id)

# contrast matrix
cont.matrix.id <- makeContrasts (MMvsHD = MM-HD,levels=designMat.id)
print(cont.matrix.id)

# linear model fit
require(limma)
fit.id<-lmFit(eset_filtered.id, designMat.id)
fit.main.id<-contrasts.fit(fit.id, cont.matrix.id)
fit.main.id<-eBayes(fit.main.id)
class(fit.main.id)

# top Tabs
topTab_MMvsHD.id <- topTable (fit.main.id, number=nrow(fit.main.id), coef="MMvsHD", adjust="fdr" ) 
head(topTab_MMvsHD.id)

# probe id annotation
annotatedTopTable <- function(topTab, anotPackage)
{
  topTab <- cbind(PROBEID=rownames(topTab), topTab)
  myProbes <- rownames(topTab)
  thePackage <- eval(parse(text = anotPackage))
  geneAnots <- select(thePackage, myProbes, c("SYMBOL", "ENTREZID", "GENENAME"))
  annotatedTopTab<- merge(x=geneAnots, y=topTab, by.x="PROBEID", by.y="PROBEID")
  return(annotatedTopTab)
}

topAnnotated_MMvsHD.id <- annotatedTopTable(topTab_MMvsHD.id,anotPackage="clariomsmousehttranscriptcluster.db")

write.csv(topAnnotated_MMvsHD.id, file="./results/DEG.ID.csv")

# filtering genes (variance/ ID)
filtered.id.var <- nsFilter(eset_rma, 
                     require.entrez = TRUE, remove.dupEntrez = TRUE,
                     var.filter=TRUE, var.func=IQR, var.cutoff=0.75, 
                     filterByQuantile=TRUE, feature.exclude = "^AFFX")
names(filtered.id.var)
class(filtered.id.var$eset)
print(filtered.id.var$filter.log)
eset_filtered.id.var <-filtered.id.var$eset

# plots of normalized filtered data
plotPCA3(exprs(eset_filtered.id.var), labels = my.targets$ShortName, factor = my.targets$Group, 
         title="Normalized filtered data", scale = FALSE, size = 3, 
         colores = c("blue","red"))

# saving files
write.csv(exprs(eset_filtered.id.var), file="./results/Normalized.Filtered(ID,Var).Unannotated.Data.csv")

# probe id annotation
eset_rma_filtered.id.var_annotated<-annotated(exprs(eset_filtered.id.var),anotPackage="clariomsmousehttranscriptcluster.db")

# saving files
write.csv(eset_rma_filtered.id.var_annotated, file="./results/Normalized.Filtered(ID,Var).Annotated.Data.csv",row.names = FALSE)

# NetworkAnalyst formatting
eset_rma_filtered.id.var_annotated[,c("PROBEID","SYMBOL")] <- NULL
class <- data.frame("#CLASS", "MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","HD",	"HD",	"HD") 
colnames(class)<-colnames(eset_rma_filtered.id.var_annotated)
eset_rma_filtered.id.var_annotated <- rbind(class,eset_rma_filtered.id.var_annotated)
colnames(eset_rma_filtered.id.var_annotated)[1] <- "#NAME"

# saving files
write.table(eset_rma_filtered.id.var_annotated, file="./results/NA.Normalized.Filtered(ID,Var).Annotated.Data.txt", row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)

## differential gene experssion
# design matrix
require(limma)
designMat.id.var<- model.matrix(~0+Group, pData(eset_filtered.id.var))
colnames(designMat.id.var) <- c("HD","MM")
print(designMat.id.var)

# contrast matrix
cont.matrix.id.var <- makeContrasts (MMvsHD = MM-HD,levels=designMat.id.var)
print(cont.matrix.id.var)

# linear model fit
require(limma)
fit.id.var<-lmFit(eset_filtered.id.var, designMat.id.var)
fit.main.id.var<-contrasts.fit(fit.id.var, cont.matrix.id.var)
fit.main.id.var<-eBayes(fit.main.id.var)
class(fit.main.id.var)

# top Tabs
topTab_MMvsHD.id.var <- topTable (fit.main.id.var, number=nrow(fit.main.id.var), coef="MMvsHD", adjust="fdr" ) 
head(topTab_MMvsHD.id.var)

# probe id annotation
topAnnotated_MMvsHD.id.var <- annotatedTopTable(topTab_MMvsHD.id.var,anotPackage="clariomsmousehttranscriptcluster.db")

write.csv(topAnnotated_MMvsHD.id.var, file="./results/DEG.ID.Var.csv")



# ortholog mapping from mouse to human of topAnnotated_MMvsHD.id
musGenes <- topAnnotated_MMvsHD.id$SYMBOL

convertMouseGeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl",host="http://dec2021.archive.ensembl.org")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl",host="http://dec2021.archive.ensembl.org")
  
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  humanx <- unique(genesV2[,c(1,2)])
  
  return(humanx)
}

gene_symbol_mapped.hs <- convertMouseGeneList(musGenes)


colnames(topAnnotated_MMvsHD.id)[2] <- "MGI.symbol"
topAnnotated_MMvsHD.id_hs.mapped <- merge(topAnnotated_MMvsHD.id,gene_symbol_mapped.hs, by = "MGI.symbol")


# check for duplicate symbols
sum(duplicated(topAnnotated_MMvsHD.id_hs.mapped$HGNC.symbol))
sum(duplicated(topAnnotated_MMvsHD.id_hs.mapped$MGI.symbol))
## There are many to many mapping

# summarize multiple MGI.symbol mapping to same HGNC.symbol 
library(dplyr)
topAnnotated_MMvsHD.id_hs.mapped <- topAnnotated_MMvsHD.id_hs.mapped[, -c(2,4)]
topAnnotated_MMvsHD.id_hs.mapped <- topAnnotated_MMvsHD.id_hs.mapped %>% 
  group_by(HGNC.symbol) %>% 
  mutate(MGI.symbol = paste0(MGI.symbol, collapse = ","), ENTREZID = paste0(ENTREZID, collapse = ",")) %>% 
  slice_min(order_by = adj.P.Val) %>% as.data.frame()

write.csv(topAnnotated_MMvsHD.id_hs.mapped,"results/DEG.ID.HS.Mapped.csv")


# ortholog mapping from mouse to human of topAnnotated_MMvsHD.id.var
musGenes <- topAnnotated_MMvsHD.id.var$SYMBOL


gene_symbol_mapped.hs <- convertMouseGeneList(musGenes)


#gene_symbol_mapped.hs <- convert_mouse_to_human(musGenes)


colnames(topAnnotated_MMvsHD.id.var)[2] <- "MGI.symbol"
topAnnotated_MMvsHD.id.var_hs.mapped <- merge(topAnnotated_MMvsHD.id.var,gene_symbol_mapped.hs, by = "MGI.symbol")


# check for duplicate symbols
sum(duplicated(topAnnotated_MMvsHD.id.var_hs.mapped$HGNC.symbol))
sum(duplicated(topAnnotated_MMvsHD.id.var_hs.mapped$MGI.symbol))
## There are many to many mapping

# summarize multiple MGI.symbol mapping to same HGNC.symbol 
topAnnotated_MMvsHD.id.var_hs.mapped <- topAnnotated_MMvsHD.id.var_hs.mapped[, -c(2,4)]
topAnnotated_MMvsHD.id.var_hs.mapped <- topAnnotated_MMvsHD.id.var_hs.mapped %>% 
  group_by(HGNC.symbol) %>% 
  mutate(MGI.symbol = paste0(MGI.symbol, collapse = ","), ENTREZID = paste0(ENTREZID, collapse = ",")) %>% 
  slice_min(order_by = adj.P.Val) %>% as.data.frame()

write.csv(topAnnotated_MMvsHD.id.var_hs.mapped,"results/DEG.ID.Var.HS.Mapped.csv")




























