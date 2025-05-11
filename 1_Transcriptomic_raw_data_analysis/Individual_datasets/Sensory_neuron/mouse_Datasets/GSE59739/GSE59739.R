# set working directory
setwd("~/Work/PhD/Projects/Comput MM Pain/Computational Analysis/Ligand_Receptor/Data/Sensory_Neuron/Neuron/mDatasets/GSE59739")

# read data
library("readxl")
GSE59739_tpm.c <- read_excel("data/GSE59739.xlsx")

# drop non-neuronal samples  
GSE59739_tpm.c[633:708] <- list(NULL) 

# drop non-relevant columns  
GSE59739_tpm.c[2:10] <- list(NULL) 

# drop non-relevant rows  
GSE59739_tpm.c <- GSE59739_tpm.c[-c(1:11), ] 

# write file 
write.csv(GSE59739_tpm.c,file="data/GSE59739_tpm.c.csv",row.names = FALSE)

# read file 
GSE59739_tpm.c <- read.csv("data/GSE59739_tpm.c.csv",row.names = 1)

# convert dataframe from character  to numerical, this code remove rownames 
GSE59739_tpm <- as.data.frame(sapply(GSE59739_tpm.c, as.numeric))
rownames(GSE59739_tpm) <- rownames(GSE59739_tpm.c)

# sum counts per gene
GSE59739_tpm$tpm_counts <- rowMeans(GSE59739_tpm)

# drop samples 
GSE59739_tpm <- GSE59739_tpm[623]

# check for duplicate symbols
sum(duplicated(rownames(GSE59739_tpm)))

# make a column for gene symbol
GSE59739_tpm$SYMBOL <- rownames(GSE59739_tpm)

# get ensembl id from symbol 
library(EnsDb.Mmusculus.v79)
ensembl_id <- ensembldb::select(EnsDb.Mmusculus.v79, keys= rownames(GSE59739_tpm), keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))
GSE59739_tpm <- merge(GSE59739_tpm, ensembl_id, by="SYMBOL")

# summarize multiple ensembl id mapping to same symbol 
require(dplyr)
GSE59739_tpm <- GSE59739_tpm %>% 
                group_by(SYMBOL) %>% 
                mutate(GENEID = paste0(GENEID, collapse = ", ")) 

# remove duplicates
GSE59739_tpm <- GSE59739_tpm[!duplicated(GSE59739_tpm), ]

# filer based on hard threshold (tpm >= 0.5)
GSE59739_tpm_filtered <- GSE59739_tpm[GSE59739_tpm$tpm_counts >= 0.5, ] 

# sort genes based on tpm counts
GSE59739_tpm_filtered_sorted <- GSE59739_tpm_filtered[order(-GSE59739_tpm_filtered$tpm_counts), ]

# check for dupicate symbols
sum(duplicated(GSE59739_tpm_filtered_sorted$SYMBOL))

# write file
GSE59739_tpm_filtered_sorted <- GSE59739_tpm_filtered_sorted[,c(3,1,2)]
write.csv(GSE59739_tpm_filtered_sorted,file="results/GSE59739_tpm_filtered_sorted.csv",row.names = FALSE)

