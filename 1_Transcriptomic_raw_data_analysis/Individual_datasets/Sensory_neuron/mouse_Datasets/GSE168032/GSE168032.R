# set working directory
setwd("~/Work/PhD/Projects/Comput MM Pain/Computational Analysis/Ligand_Receptor/Data/Sensory_Neuron/Neuron/mDatasets/GSE168032")

# read data
GSE168032_tpm <- read.table("data/GSE168032_tpm_counts.txt", header = TRUE)

# drop non-relevant columns 
GSE168032_tpm[,c(3,4)] <- list(NULL) 

# sum counts of different samples
GSE168032_tpm$tpm_counts <- rowMeans(GSE168032_tpm[3:55])

# drop samples 
GSE168032_tpm <- GSE168032_tpm[,c(1,2,56)]

# check for duplicate symbols and ensembl id
sum(duplicated(GSE168032_tpm$gene_symbol))
sum(duplicated(GSE168032_tpm$ensembl))

# extract duplicated entries to review them
duplicated_genes_id <- GSE168032_tpm$gene_symbol[duplicated(GSE168032_tpm$gene_symbol)]
GSE168032_tpm_duplicated_id <- GSE168032_tpm[GSE168032_tpm$gene_symbol %in% duplicated_genes_id, ]
### since duplicated entries do NOT have the same gene symbol and Ensembl ID, we wil sum up counts

# summarize multiple ensembl id mapping to same symbol 
require(dplyr)
GSE168032_tpm <- GSE168032_tpm %>% 
  group_by(gene_symbol) %>% 
  mutate(ensembl = paste0(ensembl, collapse = ", ")) 

# sum counts for same gene
require(dplyr)
GSE168032_tpm <- data.frame(GSE168032_tpm %>% 
                               group_by(gene_symbol, ensembl) %>% 
                               summarise_all(funs(sum)))

# filer based on hard threshold (tpm > 0.5)
GSE168032_tpm_filtered <- GSE168032_tpm[GSE168032_tpm$tpm_counts >= 0.5, ] 

# order genes based on tpm counts
GSE168032_tpm_filtered_sorted <- GSE168032_tpm_filtered[order(-GSE168032_tpm_filtered$tpm_counts), ]

# check for dupicate symbols
sum(duplicated(GSE168032_tpm_filtered_sorted$gene_symbol))

# write file
GSE168032_tpm_filtered_sorted <- GSE168032_tpm_filtered_sorted[,c(2,1,3)]
colnames(GSE168032_tpm_filtered_sorted) <- c("GENEID","SYMBOL","tpm_counts")
write.csv(GSE168032_tpm_filtered_sorted,file="results/GSE168032_tpm_filtered_sorted.csv",row.names = FALSE)

