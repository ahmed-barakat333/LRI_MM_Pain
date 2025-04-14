# set working directory
setwd("~/Work/PhD/Projects/Comput MM Pain/Computational Analysis/Ligand_Receptor/Data/Sensory_Neuron/Neuron/mDatasets/GSE139088")

# read data
all_content <- readLines("data/GSM4130750_WT_1.csv")
skip_third <- all_content[-3]
GSE139088_raw.c <- read.csv(textConnection(skip_third), header = TRUE, stringsAsFactors = FALSE)

# drop non-relevant columns (unknown cells) and barcode row
GSE139088_raw.c["Unknown"] <- list(NULL) 
GSE139088_raw.c <- GSE139088_raw.c[-1, ]

# make rownames from first column, there are two rows with gene 1-Mar and two rows with gene 2-Mar so I used make unique (i dont know another way)
rownames(GSE139088_raw.c) <- make.names(GSE139088_raw.c[,1], unique = TRUE)
GSE139088_raw.c$Cell.identity <- NULL

# convert dataframe from character  to numerical, this code remove rownames 
GSE139088_raw <- as.data.frame(sapply(GSE139088_raw.c, as.numeric))
rownames(GSE139088_raw) <- rownames(GSE139088_raw.c)

# sum counts per gene
GSE139088_raw$raw_counts <- rowSums(GSE139088_raw)

# make a column for gene symbol
GSE139088_raw$SYMBOL <- rownames(GSE139088_raw)

# drop samples 
GSE139088_raw <- GSE139088_raw[,c(11139,11140)]

# write file
GSE139088_raw <- GSE139088_raw[,c(2,1)]
write.csv(GSE139088_raw,file="data/GSE139088_raw.csv",row.names = FALSE)

# read file
GSE139088_raw <- read.csv("data/GSE139088_raw.csv")

# get ensembl id from symbol 
library(EnsDb.Mmusculus.v79)
ensembl_id <- ensembldb::select(EnsDb.Mmusculus.v79, keys= GSE139088_raw$SYMBOL, keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))
GSE139088_raw <- merge(GSE139088_raw, ensembl_id, by="SYMBOL")

# get gene length
library (EDASeq)
ensembl_list <- GSE139088_raw$GENEID
genelength <- getGeneLengthAndGCContent(GSE139088_raw$GENEID, "mmu")

## write gene length file
write.csv(genelength,file="data/GSE139088_genelength.csv",row.names = TRUE)

## read gene length file 
genelength <- read.csv("data/GSE139088_genelength.csv", header = TRUE)

# combine files
colnames(genelength)[1] <- "GENEID"
GSE139088_raw<- merge(GSE139088_raw, genelength, by="GENEID")
GSE139088_raw$gc <- NULL
GSE139088_raw <- na.omit(GSE139088_raw)

# summarize multiple ensembl id mapping to same symbol 
require(dplyr)
GSE139088_raw <- GSE139088_raw %>% 
  group_by(SYMBOL) %>% 
  mutate(GENEID = paste0(GENEID, collapse = ", ")) 


# taking average gene length
GSE139088_raw <- data.frame(GSE139088_raw %>% 
                                 group_by(GENEID,SYMBOL) %>% 
                                 summarise_all(funs(mean)))


# calculate tpm counts
GSE139088_tpm <- GSE139088_raw
GSE139088_tpm$tpm_counts <- GSE139088_tpm$raw_counts/GSE139088_tpm$length
GSE139088_tpm$tpm_counts <-  t(t(GSE139088_tpm$tpm_counts) * 1e6 / sum(GSE139088_tpm$tpm_counts))

# filer based on hard threshold (tpm >= 0.5)
GSE139088_tpm_filtered <- GSE139088_tpm[GSE139088_tpm$tpm_counts >= 0.5, ] 

# sort genes based on tpm counts
GSE139088_tpm_filtered_sorted <- GSE139088_tpm_filtered[order(-GSE139088_tpm_filtered$tpm_counts), ]

# check for dupicate symbols
sum(duplicated(GSE139088_tpm_filtered_sorted$SYMBOL))

# write file
GSE139088_tpm_filtered_sorted[3:4] <- list(NULL)
write.csv(GSE139088_tpm_filtered_sorted,file="GSE139088_tpm_filtered_sorted.csv",row.names = FALSE)


