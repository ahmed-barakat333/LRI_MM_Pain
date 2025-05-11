# set working directory
setwd("~/Work/PhD/Projects/MM Pain-Bone/Computational Analysis/Sensory Data/Neuron/mDatasets/GSE100175")

# read data
GSE100175_all_fpkm <- read.csv("data/cufflinks_w_expression_cutoff.csv", header = TRUE)

# drop non-relevant columns 
GSE100175_all_fpkm[,c(2,4,5,6,7,8,9,10,11,12,13,14,15,16,23,24,25,26)] <- list(NULL) 

# check for duplicated Ensembl IDs
sum(duplicated(GSE100175_all_fpkm$EnsemblID))

# extract duplicated entries to review them
duplicated_genes_id <- GSE100175_all_fpkm$EnsemblID[duplicated(GSE100175_all_fpkm$EnsemblID)]
GSE100175_all_fpkm_duplicated_id <- GSE100175_all_fpkm[GSE100175_all_fpkm$EnsemblID %in% duplicated_genes_id, ]
# since duplicated entries have the same gene symbol and Ensembl ID, we can keep the one with the highest counts

# check for duplicated symbols
sum(duplicated(GSE100175_all_fpkm$gene_short_name.1))

# extract duplicated entries to review them
duplicated_genes_symbol <- GSE100175_all_fpkm$gene_short_name.1[duplicated(GSE100175_all_fpkm$gene_short_name.1)]
GSE100175_all_fpkm_duplicated_symbol <- GSE100175_all_fpkm[GSE100175_all_fpkm$gene_short_name.1 %in% duplicated_genes_symbol, ]
# since many duplicated entries are just 0-s, we can keep the one with the highest counts

# keep ensembl id with max count 
require(dplyr)
GSE100175_all_fpkm <- data.frame(GSE100175_all_fpkm %>% 
                                          group_by(EnsemblID) %>% 
                                          summarise_all(funs(max)))

# move gene id to rownames 
GSE100175_all_fpkm$cobmined_id <- paste(GSE100175_all_fpkm$EnsemblID, "/", GSE100175_all_fpkm$gene_short_name.1)
rownames(GSE100175_all_fpkm) <- GSE100175_all_fpkm$cobmined_id
GSE100175_all_fpkm[,c(1,2,9)] <- list(NULL)

# normalize from fpkm to tpm
library(GeoTcgaData)
GSE100175_all_tpm <- data.frame(fpkmToTpm_matrix(GSE100175_all_fpkm))

# put gene id into columns
GSE100175_all_tpm$combined_id <- rownames(GSE100175_all_tpm)
library(dplyr)
library(tidyr)
GSE100175_all_tpm <- GSE100175_all_tpm %>%
                     separate(combined_id, c("GENEID", "SYMBOL"), "/")

# summarize multiple ensembl id mapping to same symbol 
GSE100175_all_tpm <- GSE100175_all_tpm %>% 
                     group_by(SYMBOL) %>% 
                     mutate(GENEID = paste0(GENEID, collapse = ", ")) 


# sum counts of same gene symbol 
GSE100175_all_tpm <- data.frame(GSE100175_all_tpm %>% 
                               group_by(GENEID,SYMBOL) %>% 
                               summarise_all(funs(sum)))

# check for dupicate symbols
sum(duplicated(GSE100175_all_tpm$SYMBOL))

# subset the experiment into 6 libraries and filer based on hard threshold (tpm >= 0.5)
GSE100175_1_tpm <- GSE100175_all_tpm[,c(1,2,3)]
colnames(GSE100175_1_tpm)[3] <- "tpm_counts"
GSE100175_1_tpm_filtered <- GSE100175_1_tpm[GSE100175_1_tpm$tpm_counts >= 0.5,] 
GSE100175_1_tpm_filtered_sorted <- GSE100175_1_tpm_filtered[order(-GSE100175_1_tpm_filtered$tpm_counts), ]

GSE100175_2_tpm <- GSE100175_all_tpm[,c(1,2,4)]
colnames(GSE100175_2_tpm)[3] <- "tpm_counts"
GSE100175_2_tpm_filtered <- GSE100175_2_tpm[GSE100175_2_tpm$tpm_counts >= 0.5,] 
GSE100175_2_tpm_filtered_sorted <- GSE100175_2_tpm_filtered[order(-GSE100175_2_tpm_filtered$tpm_counts), ]

GSE100175_3_tpm <- GSE100175_all_tpm[,c(1,2,5)]
colnames(GSE100175_3_tpm)[3] <- "tpm_counts"
GSE100175_3_tpm_filtered <- GSE100175_3_tpm[GSE100175_3_tpm$tpm_counts >= 0.5,] 
GSE100175_3_tpm_filtered_sorted <- GSE100175_3_tpm_filtered[order(-GSE100175_3_tpm_filtered$tpm_counts), ]

GSE100175_4_tpm <- GSE100175_all_tpm[,c(1,2,6)]
colnames(GSE100175_4_tpm)[3] <- "tpm_counts"
GSE100175_4_tpm_filtered <- GSE100175_4_tpm[GSE100175_4_tpm$tpm_counts >= 0.5,] 
GSE100175_4_tpm_filtered_sorted <- GSE100175_4_tpm_filtered[order(-GSE100175_4_tpm_filtered$tpm_counts), ]

GSE100175_5_tpm <- GSE100175_all_tpm[,c(1,2,7)]
colnames(GSE100175_5_tpm)[3] <- "tpm_counts"
GSE100175_5_tpm_filtered <- GSE100175_5_tpm[GSE100175_5_tpm$tpm_counts >= 0.5,] 
GSE100175_5_tpm_filtered_sorted <- GSE100175_5_tpm_filtered[order(-GSE100175_5_tpm_filtered$tpm_counts), ]

GSE100175_6_tpm <- GSE100175_all_tpm[,c(1,2,8)]
colnames(GSE100175_6_tpm)[3] <- "tpm_counts"
GSE100175_6_tpm_filtered <- GSE100175_6_tpm[GSE100175_6_tpm$tpm_counts >= 0.5,] 
GSE100175_6_tpm_filtered_sorted <- GSE100175_6_tpm_filtered[order(-GSE100175_6_tpm_filtered$tpm_counts), ]


# write files
write.csv(GSE100175_1_tpm_filtered_sorted,file="results/GSE100175_1_tpm_filtered_sorted.csv",row.names = FALSE)
write.csv(GSE100175_2_tpm_filtered_sorted,file="results/GSE100175_2_tpm_filtered_sorted.csv",row.names = FALSE)
write.csv(GSE100175_3_tpm_filtered_sorted,file="results/GSE100175_3_tpm_filtered_sorted.csv",row.names = FALSE)
write.csv(GSE100175_4_tpm_filtered_sorted,file="results/GSE100175_4_tpm_filtered_sorted.csv",row.names = FALSE)
write.csv(GSE100175_5_tpm_filtered_sorted,file="results/GSE100175_5_tpm_filtered_sorted.csv",row.names = FALSE)
write.csv(GSE100175_6_tpm_filtered_sorted,file="results/GSE100175_6_tpm_filtered_sorted.csv",row.names = FALSE)
