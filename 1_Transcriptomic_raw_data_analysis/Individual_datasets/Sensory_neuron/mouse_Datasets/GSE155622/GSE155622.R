# set working directory
setwd("~/Work/PhD/Projects/Comput MM Pain/Computational Analysis/Ligand_Receptor/Data/Sensory_Neuron/Neuron/mDatasets/GSE155622")

# read data
GSE155622_fpkm <- read.table("data/GSE155622_Smart-seq2_FPKM.txt", header = TRUE)

# drop non-relevant columns 
GSE155622_fpkm <- GSE155622_fpkm[,c("X","S721_FPKM","S722_FPKM","S723_FPKM","S724_FPKM",
                                  "S726_FPKM","S731_FPKM","S732_FPKM", "S734_FPKM",
                                  "S736_FPKM","S738_FPKM","S739_FPKM","S740_FPKM",
                                  "S751_FPKM","S752_FPKM","S753_FPKM","S755_FPKM",
                                  "S756_FPKM","S757_FPKM","S758_FPKM")] 

# sum counts from different samples 
GSE155622_fpkm$fpkm_counts <- rowSums(GSE155622_fpkm[,-1])

# drop samples
GSE155622_fpkm <-  GSE155622_fpkm[,c("X","fpkm_counts")]

# check for duplicate symbols
sum(duplicated(GSE155622_fpkm$X))

# move gene id to rownames 
rownames(GSE155622_fpkm)<- GSE155622_fpkm$X
GSE155622_fpkm$X <- NULL

# normalize from fpkm to tpm
library(GeoTcgaData)
GSE155622_tpm <- data.frame(fpkmToTpm_matrix(GSE155622_fpkm))
colnames(GSE155622_tpm) <- "tpm_counts"

# get ensembl id from symbol 
library(EnsDb.Mmusculus.v79)
ensembl_id <- ensembldb::select(EnsDb.Mmusculus.v79, keys= rownames(GSE155622_tpm), keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))
GSE155622_tpm$SYMBOL <- rownames(GSE155622_tpm)
GSE155622_tpm <- merge(GSE155622_tpm, ensembl_id, by="SYMBOL")

# summarize multiple ensembl id mapping to same symbol 
require(dplyr)
GSE155622_tpm <- GSE155622_tpm %>% 
                 group_by(SYMBOL) %>% 
                 mutate(GENEID = paste0(GENEID, collapse = ", ")) 

# remove duplicates
GSE155622_tpm <- GSE155622_tpm[!duplicated(GSE155622_tpm), ]


# filer based on hard threshold (tpm >= 0.5)
GSE155622_tpm_filtered <- GSE155622_tpm[GSE155622_tpm$tpm_counts >= 0.5, ] 

# sort genes based on tpm counts
GSE155622_tpm_filtered_sorted <- GSE155622_tpm_filtered[order(-GSE155622_tpm_filtered$tpm_counts), ]

# check for dupicate symbols
sum(duplicated(GSE155622_tpm_filtered_sorted$SYMBOL))

# write file
GSE155622_tpm_filtered_sorted <- GSE155622_tpm_filtered_sorted[,c(3,1,2)]
write.csv(GSE155622_tpm_filtered_sorted,file="results/GSE155622_tpm_filtered_sorted.csv",row.names = FALSE)
