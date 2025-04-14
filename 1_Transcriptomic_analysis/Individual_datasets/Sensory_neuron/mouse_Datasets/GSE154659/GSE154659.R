# set working directory
setwd("~/Work/PhD/Projects/Comput MM Pain/Computational Analysis/Ligand_Receptor/Data/Sensory_Neuron/Neuron/mDatasets/GSE154659")

# read data
filename <- file.choose()
rds_file <- readRDS(filename)
summ <- summary(rds_file)
GSE154659_raw <- data.frame(SYMBOL      = rownames(rds_file)[summ$i],
                           sample  = colnames(rds_file)[summ$j],
                           value = summ$x)

# drop non-relevant samples 
GSE154659_raw <- dplyr::filter(GSE154659_raw, grepl('Naive', sample))
GSE154659_raw$sample <- NULL

# check for duplicate symbols and ensembl id
sum(duplicated(GSE154659_raw$SYMBOL))

# write file 
colnames(GSE154659_raw)[2] <- "raw_counts"
write.csv(GSE154659_raw,file="data/GSE154659_raw.csv",row.names = FALSE)

# read file 
GSE154659_raw <- read.csv("data/GSE154659_raw.csv",row.names = 1)

# get ensembl id from symbol 
library(EnsDb.Mmusculus.v79)
ensembl_id <- ensembldb::select(EnsDb.Mmusculus.v79, keys= rownames(GSE154659_raw), keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))
GSE154659_raw$SYMBOL <- rownames(GSE154659_raw)
GSE154659_raw <- merge(GSE154659_raw, ensembl_id, by="SYMBOL")

# get gene length
library (EDASeq)
ensembl_list <- GSE154659_raw$GENEID
genelength <- getGeneLengthAndGCContent(GSE154659_raw$GENEID, "mmu")

## write gene length file
write.csv(genelength,file="data/GSE154659_genelength.csv",row.names = TRUE)

## read gene length file 
genelength <- read.csv("data/GSE154659_genelength.csv", header = TRUE)

# combine files
colnames(genelength)[1] <- "GENEID"
GSE154659_raw<- merge(GSE154659_raw, genelength, by="GENEID")
GSE154659_raw$gc <- NULL
GSE154659_raw <- na.omit(GSE154659_raw)

# summarize multiple ensembl id mapping to same symbol 
require(dplyr)
GSE154659_raw <- GSE154659_raw %>% 
  group_by(SYMBOL) %>% 
  mutate(GENEID = paste0(GENEID, collapse = ", ")) 


# taking average gene length
GSE154659_raw <- data.frame(GSE154659_raw %>% 
                                 group_by(GENEID,SYMBOL) %>% 
                                 summarise_all(funs(mean)))


# calculate tpm counts
GSE154659_tpm <- GSE154659_raw
GSE154659_tpm$tpm_counts <- GSE154659_raw$raw_counts/GSE154659_raw$length
GSE154659_tpm$tpm_counts <-  t(t(GSE154659_tpm$tpm_counts) * 1e6 / sum(GSE154659_tpm$tpm_counts))

# filer based on hard threshold (tpm >= 0.5)
GSE154659_tpm_filtered <- GSE154659_tpm[GSE154659_tpm$tpm_counts >= 0.5, ] 

# sort genes based on tpm counts
GSE154659_tpm_filtered_sorted <- GSE154659_tpm_filtered[order(-GSE154659_tpm_filtered$tpm_counts), ]

# check for dupicate symbols
sum(duplicated(GSE154659_tpm_filtered_sorted$SYMBOL))

# write file
GSE154659_tpm_filtered_sorted[3:4] <- list(NULL)
write.csv(GSE154659_tpm_filtered_sorted,file="results/GSE154659_tpm_filtered_sorted.csv",row.names = FALSE)


