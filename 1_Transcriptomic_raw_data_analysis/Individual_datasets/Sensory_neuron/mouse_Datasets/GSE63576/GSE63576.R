# set working directory
setwd("~/Work/PhD/Projects/Comput MM Pain/Computational Analysis/Ligand_Receptor/Data/Sensory_Neuron/Neuron/mDatasets/GSE63576")

# read data
GSE63576_fpkm <- read.table("data/GSE63576_DRG_neurons_RNA_Seq.txt", header = TRUE, row.names = 1)

# sum fpkm counts from different samples
GSE63576_fpkm$fpkm_counts <- rowSums(GSE63576_fpkm)

# drop columns 
GSE63576_fpkm[1:204] <- list(NULL) 

# check for duplicate symbols and ensembl id
sum(duplicated(rownames(GSE63576_fpkm)))

# calculate tpm counts from fpkm 
library(GeoTcgaData)
GSE63576_tpm <- data.frame(fpkmToTpm_matrix(GSE63576_fpkm))
colnames(GSE63576_tpm) <- "tpm_counts"

# get ensembl id from symbol 
library(EnsDb.Mmusculus.v79)
ensembl_id <- ensembldb::select(EnsDb.Mmusculus.v79, keys= rownames(GSE63576_tpm), keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))
GSE63576_tpm$SYMBOL <- rownames(GSE63576_tpm)
GSE63576_tpm <- merge(GSE63576_tpm, ensembl_id, by="SYMBOL")

# summarize multiple ensembl id mapping to same symbol 
require(dplyr)
GSE63576_tpm <- GSE63576_tpm %>% 
                group_by(SYMBOL) %>% 
                mutate(GENEID = paste0(GENEID, collapse = ", ")) 

# remove duplicates
GSE63576_tpm <- GSE63576_tpm[!duplicated(GSE63576_tpm), ]

# filer based on hard threshold (tpm >= 0.5)
GSE63576_tpm_filtered <- GSE63576_tpm[GSE63576_tpm$tpm_counts >= 0.5, ] 

# sort genes based on tpm counts
GSE63576_tpm_filtered_sorted <- GSE63576_tpm_filtered[order(-GSE63576_tpm_filtered$tpm_counts), ]

# check for dupicate symbols
sum(duplicated(GSE63576_tpm_filtered_sorted$SYMBOL))

# write file
GSE63576_tpm_filtered_sorted <- GSE63576_tpm_filtered_sorted[,c(3,1,2)]
write.csv(GSE63576_tpm_filtered_sorted,file="results/GSE63576_tpm_filtered_sorted.csv",row.names = FALSE)
