# set working directory
setwd("~/Work/PhD/Projects/Comput MM Pain/Computational Analysis/Ligand_Receptor/Data/Sensory_Neuron/Neuron/mDatasets/GSE62424")

# read data
GSE62424_all_fpkm <- read.csv("data/GSE62424_FPKM_all.csv", header = TRUE)

# drop non-relevant columns 
GSE62424_all_fpkm[2:11] <- list(NULL) 

# check for dupicate symbols
sum(duplicated(GSE62424_all_fpkm$gene_short_name))

# normalize from fpkm to tpm
GSE62424_all_tpm <- GSE62424_all_fpkm
GSE62424_all_tpm$MACS.1_tpm_counts <- GSE62424_all_tpm$MACS.1 / sum(GSE62424_all_tpm$MACS.1) * 1e6
GSE62424_all_tpm$MACS.2_tpm_counts <- GSE62424_all_tpm$MACS.2 / sum(GSE62424_all_tpm$MACS.2) * 1e6
GSE62424_all_tpm$MACS.3_tpm_counts <- GSE62424_all_tpm$MACS.3 / sum(GSE62424_all_tpm$MACS.3) * 1e6
GSE62424_all_tpm$MACS.4_tpm_counts <- GSE62424_all_tpm$MACS.4 / sum(GSE62424_all_tpm$MACS.4) * 1e6
GSE62424_all_tpm[2:5] <- list(NULL)

# sum count of same gene symbol
require(dplyr)
GSE62424_all_tpm <- data.frame(GSE62424_all_tpm %>% 
                                 group_by(gene_short_name) %>% 
                                 summarise_all(funs(sum)))


# get ensembl id from symbol 
library(EnsDb.Mmusculus.v79)
ensembl_id <- ensembldb::select(EnsDb.Mmusculus.v79, keys= GSE62424_all_tpm$gene_short_name, keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))
colnames(GSE62424_all_tpm)[1] <- "SYMBOL"
GSE62424_all_tpm <- merge(GSE62424_all_tpm, ensembl_id, by="SYMBOL")

# summarize multiple ensembl id mapping to same symbol 
require(dplyr)
GSE62424_all_tpm <- GSE62424_all_tpm %>% 
                    group_by(SYMBOL) %>% 
                    mutate(GENEID = paste0(GENEID, collapse = ", ")) 

# remove duplicates
GSE62424_all_tpm <- GSE62424_all_tpm[!duplicated(GSE62424_all_tpm), ]

# summarize (average) 4 samples into one
GSE62424_all_tpm$tpm_counts <- rowMeans(GSE62424_all_tpm[,-c(1,6)])
GSE62424_all_tpm[2:5] <- list(NULL) 

#  filer based on hard threshold (tpm >= 0.5)
GSE62424_all_tpm_filtered <- GSE62424_all_tpm[GSE62424_all_tpm$tpm_counts >= 0.5,] 

# sort based on tpm counts
GSE62424_all_tpm_filtered_sorted <- GSE62424_all_tpm_filtered[order(-GSE62424_all_tpm_filtered$tpm_counts), ]

# write file
GSE62424_all_tpm_filtered_sorted <- GSE62424_all_tpm_filtered_sorted[,c(2,1,3)]
write.csv(GSE62424_all_tpm_filtered_sorted,file="results/GSE62424_all_tpm_filtered_sorted.csv",row.names = FALSE)
