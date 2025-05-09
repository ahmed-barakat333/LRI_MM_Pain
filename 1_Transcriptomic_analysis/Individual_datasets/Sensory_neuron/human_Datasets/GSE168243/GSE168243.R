# set working directory
setwd("~/Work/PhD/Projects/Comput MM Pain/Computational Analysis/Ligand_Receptor/Data/Sensory_Neuron/Neuron/hDatasets /GSE168243")

# read  files
GSM5134537_prep1_raw <- read.csv("data/GSM5134537_prep1.csv", row.names = 1)

GSM5134538_prep2_raw <- read.csv("data/GSM5134538_prep2.csv", row.names = 1)

GSM5134539_prep3_raw <- read.csv("data/GSM5134539_prep3.csv", row.names = 1)

GSM5134540_prep4_raw <- read.csv("data/GSM5134540_prep4.csv", row.names = 1)

GSM5134541_prep5_raw <- read.csv("data/GSM5134541_prep5.csv", row.names = 1)

GSM5134542_prep6_raw <- read.csv("data/GSM5134542_prep6.csv", row.names = 1)

# sum counts from different cells 
GSM5134537_prep1_raw$raw_counts <- rowSums(GSM5134537_prep1_raw)
GSM5134537_prep1_raw <- GSM5134537_prep1_raw[213]

GSM5134538_prep2_raw$raw_counts <- rowSums(GSM5134538_prep2_raw)
GSM5134538_prep2_raw <- GSM5134538_prep2_raw[153]

GSM5134539_prep3_raw$raw_counts <- rowSums(GSM5134539_prep3_raw)
GSM5134539_prep3_raw <- GSM5134539_prep3_raw[771]

GSM5134540_prep4_raw$raw_counts <- rowSums(GSM5134540_prep4_raw)
GSM5134540_prep4_raw <- GSM5134540_prep4_raw[282]

GSM5134541_prep5_raw$raw_counts <- rowSums(GSM5134541_prep5_raw)
GSM5134541_prep5_raw <- GSM5134541_prep5_raw[81]

GSM5134542_prep6_raw$raw_counts <- rowSums(GSM5134542_prep6_raw)
GSM5134542_prep6_raw <- GSM5134542_prep6_raw[343]

# get ensembl id from symbol 
library(EnsDb.Hsapiens.v86)
ensembl_id <- ensembldb::select(EnsDb.Hsapiens.v86, keys= rownames(GSM5134537_prep1_raw), keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))

GSM5134537_prep1_raw$SYMBOL <- rownames(GSM5134537_prep1_raw)
GSM5134537_prep1_raw <- merge(GSM5134537_prep1_raw, ensembl_id, by="SYMBOL")

GSM5134538_prep2_raw$SYMBOL <- rownames(GSM5134538_prep2_raw)
GSM5134538_prep2_raw <- merge(GSM5134538_prep2_raw, ensembl_id, by="SYMBOL")

GSM5134539_prep3_raw$SYMBOL <- rownames(GSM5134539_prep3_raw)
GSM5134539_prep3_raw <- merge(GSM5134539_prep3_raw, ensembl_id, by="SYMBOL")

GSM5134540_prep4_raw$SYMBOL <- rownames(GSM5134540_prep4_raw)
GSM5134540_prep4_raw <- merge(GSM5134540_prep4_raw, ensembl_id, by="SYMBOL")

GSM5134541_prep5_raw$SYMBOL <- rownames(GSM5134541_prep5_raw)
GSM5134541_prep5_raw <- merge(GSM5134541_prep5_raw, ensembl_id, by="SYMBOL")

GSM5134542_prep6_raw$SYMBOL <- rownames(GSM5134542_prep6_raw)
GSM5134542_prep6_raw <- merge(GSM5134542_prep6_raw, ensembl_id, by="SYMBOL")

# get gene length
library (EDASeq)
ensembl_list <- ensembl_id$GENEID
genelength <- getGeneLengthAndGCContent(ensembl_id$GENEID, "hsa")

## write gene length file
write.csv(genelength,file="data/GSE168243_genelength.csv",row.names = TRUE)

## read gene length file 
genelength <- read.csv("data/GSE168243_genelength.csv", header = TRUE)

# combine files
colnames(genelength)[1] <- "GENEID"
genelength$gc <- NULL

GSM5134537_prep1_raw <- merge(GSM5134537_prep1_raw, genelength, by="GENEID")
GSM5134537_prep1_raw <- na.omit(GSM5134537_prep1_raw)

GSM5134538_prep2_raw <- merge(GSM5134538_prep2_raw, genelength, by="GENEID")
GSM5134538_prep2_raw <- na.omit(GSM5134538_prep2_raw)

GSM5134539_prep3_raw <- merge(GSM5134539_prep3_raw, genelength, by="GENEID")
GSM5134539_prep3_raw <- na.omit(GSM5134539_prep3_raw)

GSM5134540_prep4_raw <- merge(GSM5134540_prep4_raw, genelength, by="GENEID")
GSM5134540_prep4_raw <- na.omit(GSM5134540_prep4_raw)

GSM5134541_prep5_raw <- merge(GSM5134541_prep5_raw, genelength, by="GENEID")
GSM5134541_prep5_raw <- na.omit(GSM5134541_prep5_raw)

GSM5134542_prep6_raw <- merge(GSM5134542_prep6_raw, genelength, by="GENEID")
GSM5134542_prep6_raw <- na.omit(GSM5134542_prep6_raw)


# summarize multiple ensembl id mapping to same symbol 
require(dplyr)

GSM5134537_prep1_raw <- GSM5134537_prep1_raw %>% 
  group_by(SYMBOL) %>% 
  mutate(GENEID = paste0(GENEID, collapse = ", ")) 

GSM5134538_prep2_raw <- GSM5134538_prep2_raw %>% 
  group_by(SYMBOL) %>% 
  mutate(GENEID = paste0(GENEID, collapse = ", ")) 

GSM5134539_prep3_raw <- GSM5134539_prep3_raw %>% 
  group_by(SYMBOL) %>% 
  mutate(GENEID = paste0(GENEID, collapse = ", ")) 

GSM5134540_prep4_raw <- GSM5134540_prep4_raw %>% 
  group_by(SYMBOL) %>% 
  mutate(GENEID = paste0(GENEID, collapse = ", ")) 

GSM5134541_prep5_raw <- GSM5134541_prep5_raw %>% 
  group_by(SYMBOL) %>% 
  mutate(GENEID = paste0(GENEID, collapse = ", ")) 

GSM5134542_prep6_raw <- GSM5134542_prep6_raw %>% 
  group_by(SYMBOL) %>% 
  mutate(GENEID = paste0(GENEID, collapse = ", ")) 


# taking average gene length
GSM5134537_prep1_raw <- data.frame(GSM5134537_prep1_raw %>% 
                              group_by(GENEID,SYMBOL) %>% 
                              summarise_all(funs(mean)))

GSM5134538_prep2_raw <- data.frame(GSM5134538_prep2_raw %>% 
                                     group_by(GENEID,SYMBOL) %>% 
                                     summarise_all(funs(mean)))

GSM5134539_prep3_raw <- data.frame(GSM5134539_prep3_raw %>% 
                                     group_by(GENEID,SYMBOL) %>% 
                                     summarise_all(funs(mean)))

GSM5134540_prep4_raw <- data.frame(GSM5134540_prep4_raw %>% 
                                     group_by(GENEID,SYMBOL) %>% 
                                     summarise_all(funs(mean)))

GSM5134541_prep5_raw <- data.frame(GSM5134541_prep5_raw %>% 
                                     group_by(GENEID,SYMBOL) %>% 
                                     summarise_all(funs(mean)))

GSM5134542_prep6_raw <- data.frame(GSM5134542_prep6_raw %>% 
                                     group_by(GENEID,SYMBOL) %>% 
                                     summarise_all(funs(mean)))


# calculate tpm counts
GSM5134537_prep1_tpm <- GSM5134537_prep1_raw
GSM5134537_prep1_tpm$tpm_counts <- GSM5134537_prep1_raw$raw_counts/GSM5134537_prep1_raw$length
GSM5134537_prep1_tpm$tpm_counts <-  t(t(GSM5134537_prep1_tpm$tpm_counts) * 1e6 / sum(GSM5134537_prep1_tpm$tpm_counts))

GSM5134538_prep2_tpm <- GSM5134538_prep2_raw
GSM5134538_prep2_tpm$tpm_counts <- GSM5134538_prep2_raw$raw_counts/GSM5134538_prep2_raw$length
GSM5134538_prep2_tpm$tpm_counts <-  t(t(GSM5134538_prep2_tpm$tpm_counts) * 1e6 / sum(GSM5134538_prep2_tpm$tpm_counts))

GSM5134539_prep3_tpm <- GSM5134539_prep3_raw
GSM5134539_prep3_tpm$tpm_counts <- GSM5134539_prep3_raw$raw_counts/GSM5134539_prep3_raw$length
GSM5134539_prep3_tpm$tpm_counts <-  t(t(GSM5134539_prep3_tpm$tpm_counts) * 1e6 / sum(GSM5134539_prep3_tpm$tpm_counts))

GSM5134540_prep4_tpm <- GSM5134540_prep4_raw
GSM5134540_prep4_tpm$tpm_counts <- GSM5134540_prep4_raw$raw_counts/GSM5134540_prep4_raw$length
GSM5134540_prep4_tpm$tpm_counts <-  t(t(GSM5134540_prep4_tpm$tpm_counts) * 1e6 / sum(GSM5134540_prep4_tpm$tpm_counts))

GSM5134541_prep5_tpm <- GSM5134541_prep5_raw
GSM5134541_prep5_tpm$tpm_counts <- GSM5134541_prep5_raw$raw_counts/GSM5134541_prep5_raw$length
GSM5134541_prep5_tpm$tpm_counts <-  t(t(GSM5134541_prep5_tpm$tpm_counts) * 1e6 / sum(GSM5134541_prep5_tpm$tpm_counts))

GSM5134542_prep6_tpm <- GSM5134542_prep6_raw
GSM5134542_prep6_tpm$tpm_counts <- GSM5134542_prep6_raw$raw_counts/GSM5134542_prep6_raw$length
GSM5134542_prep6_tpm$tpm_counts <-  t(t(GSM5134542_prep6_tpm$tpm_counts) * 1e6 / sum(GSM5134542_prep6_tpm$tpm_counts))


# filer based on hard threshold (tpm >= 0.5)
GSM5134537_prep1_tpm_filtered <- GSM5134537_prep1_tpm[GSM5134537_prep1_tpm$tpm_counts >= 0.5, ] 

GSM5134538_prep2_tpm_filtered <- GSM5134538_prep2_tpm[GSM5134538_prep2_tpm$tpm_counts >= 0.5, ] 

GSM5134539_prep3_tpm_filtered <- GSM5134539_prep3_tpm[GSM5134539_prep3_tpm$tpm_counts >= 0.5, ] 

GSM5134540_prep4_tpm_filtered <- GSM5134540_prep4_tpm[GSM5134540_prep4_tpm$tpm_counts >= 0.5, ] 

GSM5134541_prep5_tpm_filtered <- GSM5134541_prep5_tpm[GSM5134541_prep5_tpm$tpm_counts >= 0.5, ] 

GSM5134542_prep6_tpm_filtered <- GSM5134542_prep6_tpm[GSM5134542_prep6_tpm$tpm_counts >= 0.5, ] 

# sort genes based on tpm counts
GSM5134537_prep1_tpm_filtered_sorted <- GSM5134537_prep1_tpm_filtered[order(-GSM5134537_prep1_tpm_filtered$tpm_counts), ]

GSM5134538_prep2_tpm_filtered_sorted <- GSM5134538_prep2_tpm_filtered[order(-GSM5134538_prep2_tpm_filtered$tpm_counts), ]

GSM5134539_prep3_tpm_filtered_sorted <- GSM5134539_prep3_tpm_filtered[order(-GSM5134539_prep3_tpm_filtered$tpm_counts), ]

GSM5134540_prep4_tpm_filtered_sorted <- GSM5134540_prep4_tpm_filtered[order(-GSM5134540_prep4_tpm_filtered$tpm_counts), ]

GSM5134541_prep5_tpm_filtered_sorted <- GSM5134541_prep5_tpm_filtered[order(-GSM5134541_prep5_tpm_filtered$tpm_counts), ]

GSM5134542_prep6_tpm_filtered_sorted <- GSM5134542_prep6_tpm_filtered[order(-GSM5134542_prep6_tpm_filtered$tpm_counts), ]

# check for dupicate symbols
sum(duplicated(GSM5134537_prep1_tpm_filtered_sorted$SYMBOL))

sum(duplicated(GSM5134538_prep2_tpm_filtered_sorted$SYMBOL))

sum(duplicated(GSM5134539_prep3_tpm_filtered_sorted$SYMBOL))

sum(duplicated(GSM5134540_prep4_tpm_filtered_sorted$SYMBOL))

sum(duplicated(GSM5134541_prep5_tpm_filtered_sorted$SYMBOL))

sum(duplicated(GSM5134542_prep6_tpm_filtered_sorted$SYMBOL))

# write file
GSM5134537_prep1_tpm_filtered_sorted[3:4] <- list(NULL)
write.csv(GSM5134537_prep1_tpm_filtered_sorted,file="GSM5134537_prep1_tpm_filtered_sorted.csv",row.names = FALSE)

GSM5134538_prep2_tpm_filtered_sorted[3:4] <- list(NULL)
write.csv(GSM5134538_prep2_tpm_filtered_sorted,file="GSM5134538_prep2_tpm_filtered_sorted.csv",row.names = FALSE)

GSM5134539_prep3_tpm_filtered_sorted[3:4] <- list(NULL)
write.csv(GSM5134539_prep3_tpm_filtered_sorted,file="GSM5134539_prep3_tpm_filtered_sorted.csv",row.names = FALSE)

GSM5134540_prep4_tpm_filtered_sorted[3:4] <- list(NULL)
write.csv(GSM5134540_prep4_tpm_filtered_sorted,file="GSM5134540_prep4_tpm_filtered_sorted.csv",row.names = FALSE)

GSM5134541_prep5_tpm_filtered_sorted[3:4] <- list(NULL)
write.csv(GSM5134541_prep5_tpm_filtered_sorted,file="GSM5134541_prep5_tpm_filtered_sorted.csv",row.names = FALSE)

GSM5134542_prep6_tpm_filtered_sorted[3:4] <- list(NULL)
write.csv(GSM5134542_prep6_tpm_filtered_sorted,file="GSM5134542_prep6_tpm_filtered_sorted.csv",row.names = FALSE)

# prepare datasets for meta-analysis
GSM5134537_prep1 <- GSM5134537_prep1_tpm_filtered_sorted$SYMBOL

GSM5134538_prep2 <- GSM5134538_prep2_tpm_filtered_sorted$SYMBOL

GSM5134539_prep3 <- GSM5134539_prep3_tpm_filtered_sorted$SYMBOL

GSM5134540_prep4 <- GSM5134540_prep4_tpm_filtered_sorted$SYMBOL

GSM5134541_prep5 <- GSM5134541_prep5_tpm_filtered_sorted$SYMBOL

GSM5134542_prep6 <- GSM5134542_prep6_tpm_filtered_sorted$SYMBOL

# rank meta-analysis
library(RobustRankAggreg)

glist <- list(GSM5134537_prep1 = GSM5134537_prep1, GSM5134538_prep2 = GSM5134538_prep2, GSM5134539_prep3 = GSM5134539_prep3, GSM5134540_prep4 = GSM5134540_prep4, GSM5134541_prep5 = GSM5134541_prep5, GSM5134542_prep6 = GSM5134542_prep6)

hs_neuron_ma <- aggregateRanks(glist = glist)

# write file
write.csv(hs_neuron_ma,file="HS_Sensory_Neuron_MA.csv",row.names = FALSE)


