# set working directory
setwd("~/Work/PhD/Projects/Comput MM Pain/Computational Analysis/Ligand_Receptor/Data/Sensory_Neuron/Neuron/mDatasets/GSE131230")

# read the file
GSE131230_raw <- read.csv("data/GSE131230_counts_official.csv", row.names = 1)

# sum counts from different samples 
GSE131230_raw$raw_counts <- rowSums(GSE131230_raw)
GSE131230_raw <- GSE131230_raw[27]

# get symbol from ensembl id
library(EnsDb.Mmusculus.v79)
ensembl_id <- ensembldb::select(EnsDb.Mmusculus.v79, keys= rownames(GSE131230_raw), keytype = "GENEID", columns = c("GENEID","SYMBOL"))
GSE131230_raw$GENEID <- rownames(GSE131230_raw)
GSE131230_raw <- merge(GSE131230_raw, ensembl_id, by="GENEID")

# get gene length
library (EDASeq)
ensembl_list <- GSE131230_raw$GENEID
genelength <- getGeneLengthAndGCContent(GSE131230_raw$GENEID, "mmu")

## write gene length file
write.csv(genelength,file="data/GSE131230_genelength.csv",row.names = TRUE)

## read gene length file 
genelength <- read.csv("data/GSE131230_genelength.csv", header = TRUE)

# combine files
colnames(genelength)[1] <- "GENEID"
GSE131230_raw<- merge(GSE131230_raw, genelength, by="GENEID")
GSE131230_raw$gc <- NULL
GSE131230_raw <- na.omit(GSE131230_raw)

# calculate tpm counts
GSE131230_tpm <- GSE131230_raw
GSE131230_tpm$tpm_counts <- GSE131230_raw$raw_counts/GSE131230_raw$length
GSE131230_tpm$tpm_counts <-  t(t(GSE131230_tpm$tpm_counts) * 1e6 / sum(GSE131230_tpm$tpm_counts))
GSE131230_tpm[,c(2,4)] <- list(NULL)

# summarize multiple ensembl id mapping to same symbol 
require(dplyr)
GSE131230_tpm <- GSE131230_tpm %>% 
  group_by(SYMBOL) %>% 
  mutate(GENEID = paste0(GENEID, collapse = ", ")) 


# sum counts of same gene symbol 
GSE131230_tpm <- data.frame(GSE131230_tpm %>% 
                                 group_by(GENEID,SYMBOL) %>% 
                                 summarise_all(funs(sum)))

# check for dupicate symbols
sum(duplicated(GSE131230_tpm$SYMBOL))

# filer based on hard threshold (tpm >= 0.5)
GSE131230_tpm_filtered <- GSE131230_tpm[GSE131230_tpm$tpm_counts >= 0.5, ] 

# sort genes based on tpm counts
GSE131230_tpm_filtered_sorted <- GSE131230_tpm_filtered[order(-GSE131230_tpm_filtered$tpm_counts), ]

# delete rows with empty symbol
GSE131230_tpm_filtered_sorted <- GSE131230_tpm_filtered_sorted[-which(GSE131230_tpm_filtered_sorted$SYMBOL == ""), ]

# write file
write.csv(GSE131230_tpm_filtered_sorted,file="GSE131230_tpm_filtered_sorted.csv",row.names = FALSE)


