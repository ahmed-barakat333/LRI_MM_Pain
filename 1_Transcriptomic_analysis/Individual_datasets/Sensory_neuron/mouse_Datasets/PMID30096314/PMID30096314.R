# set working directory
setwd("~/Work/PhD/Projects/Comput MM Pain/Computational Analysis/Ligand_Receptor/Data/Sensory_Neuron/Neuron/mDatasets/PMID30096314")

# load package to read loom file
library(loomR)

# read the file
lfile <- connect(filename = "data/l1_drg.loom", mode = "r", skip.validate = T)
lfile

# content of the file
lfile[["matrix"]]

lfile[["col_attrs"]]

lfile[["row_attrs"]] 

# prepare dataframe of gene cell matrix
all_cells.genes <- data.frame(lfile[["matrix"]][,])
dim(x = all_cells.genes)

gene_symbol <- lfile[["row_attrs/Gene"]][]
head(x = gene_symbol)
gene_symbol <- data.frame(as.list(gene_symbol))


cell_class <- lfile[["col_attrs/Class"]][]
head(x = cell_class)
cell_class <- data.frame(as.list(cell_class))


row.names(all_cells.genes) <- make.names(cell_class, unique=TRUE)
colnames(all_cells.genes ) <- gene_symbol

# remove not relevant cells (vascular, glia, astrocyte, blood)
a <- all_cells.genes[-grep("PeripheralGlia",rownames(all_cells.genes)),]
b <- a[-grep("Vascular",rownames(a)),]
c <- b[-grep("Astrocytes",rownames(b)),]
d <- c[-grep("Blood",rownames(c)),]
e <- d[-grep("Oligos",rownames(d)),]
PMID30096314_raw <- e[-grep("Excluded",rownames(e)),]
remove(a,b,c,d,e,all_cells.genes,cell_class,gene_symbol,lfile)

# sum counts per gene
PMID30096314_raw <- data.frame(t(PMID30096314_raw))
PMID30096314_raw$raw_counts <- rowSums(PMID30096314_raw)
PMID30096314_raw <- PMID30096314_raw[1738]

# check for dupicate symbols
sum(duplicated(rownames(PMID30096314_raw)))

# get ensembl id from symbol 
library(EnsDb.Mmusculus.v79)
ensembl_id <- ensembldb::select(EnsDb.Mmusculus.v79, keys= rownames(PMID30096314_raw), keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))
PMID30096314_raw$SYMBOL <- rownames(PMID30096314_raw)
PMID30096314_raw <- merge(PMID30096314_raw, ensembl_id, by="SYMBOL")

# get gene length
library (EDASeq)
ensembl_list <- geneIDs$GENEID
genelength <- getGeneLengthAndGCContent(PMID30096314_raw$GENEID, "mmu")

## write gene length file
write.csv(genelength,file="results/PMID30096314_genelength.csv",row.names = TRUE)

## read gene length file 
genelength <- read.csv("data/PMID30096314_genelength.csv", header = TRUE)

# combine files
colnames(genelength)[1] <- "GENEID"
PMID30096314_raw<- merge(PMID30096314_raw, genelength, by="GENEID")
PMID30096314_raw$gc <- NULL
PMID30096314_raw <- na.omit(PMID30096314_raw)

# summarize multiple ensembl id mapping to same symbol 
require(dplyr)
PMID30096314_raw <- PMID30096314_raw %>% 
                    group_by(SYMBOL) %>% 
                    mutate(GENEID = paste0(GENEID, collapse = ", ")) 


# taking average gene length
PMID30096314_raw <- data.frame(PMID30096314_raw %>% 
                               group_by(GENEID,SYMBOL) %>% 
                               summarise_all(funs(mean)))


# calculate tpm counts
PMID30096314_tpm <- PMID30096314_raw
PMID30096314_tpm$tpm_counts <- PMID30096314_raw$raw_counts/PMID30096314_raw$length
PMID30096314_tpm$tpm_counts <-  t(t(PMID30096314_tpm$tpm_counts) * 1e6 / sum(PMID30096314_tpm$tpm_counts))

# filer based on hard threshold (tpm >= 0.5)
PMID30096314_tpm_filtered <- PMID30096314_tpm[PMID30096314_tpm$tpm_counts >= 0.5, ] 

# sort genes based on tpm counts
PMID30096314_tpm_filtered_sorted <- PMID30096314_tpm_filtered[order(-PMID30096314_tpm_filtered$tpm_counts), ]

# check for duplicate symbols
sum(duplicated(PMID30096314_tpm_filtered_sorted$SYMBOL))

# write file
PMID30096314_tpm_filtered_sorted[3:4] <- list(NULL)
write.csv(PMID30096314_tpm_filtered_sorted,file="PMID30096314_tpm_filtered_sorted.csv",row.names = FALSE)


