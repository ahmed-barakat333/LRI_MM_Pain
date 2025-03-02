# set working directory 
setwd("~/Work/PhD/Projects/Comput MM Pain/Computational Analysis/Ligand_Receptor/Data/Sensory_Neuron/Neuron/mMA ")

# load data (18 datasets)
GSE59739 <- read.csv("data/GSE59739_tpm_filtered_sorted.csv")
GSE59739 <- trimws(GSE59739$SYMBOL)

GSE62424 <- read.csv("data/GSE62424_all_tpm_filtered_sorted.csv")
GSE62424 <- trimws(GSE62424$SYMBOL)

GSE63576 <- read.csv("data/GSE63576_tpm_filtered_sorted.csv")
GSE63576 <- trimws(GSE63576$SYMBOL)

GSE100175 <- read.csv("data/GSE100175_all_tpm_filtered_sorted.csv")
GSE100175 <- trimws(GSE100175$SYMBOL)

GSE131230 <- read.csv("data/GSE131230_tpm_filtered_sorted.csv")
GSE131230 <- trimws(GSE131230$SYMBOL)

GSE139088 <- read.csv("data/GSE139088_tpm_filtered_sorted.csv")
GSE139088 <- trimws(GSE139088$SYMBOL)

GSE154659 <- read.csv("data/GSE154659_tpm_filtered_sorted.csv")
GSE154659 <- trimws(GSE154659$SYMBOL)

GSE155622 <- read.csv("data/GSE155622_tpm_filtered_sorted.csv")
GSE155622 <- trimws(GSE155622$SYMBOL)

GSE168032 <- read.csv("data/GSE168032_tpm_filtered_sorted.csv")
GSE168032 <- trimws(GSE168032$SYMBOL)

PMID30096314 <- read.csv("data/PMID30096314_tpm_filtered_sorted.csv")
PMID30096314 <- trimws(PMID30096314$SYMBOL)

# rank meta-analysis
library(RobustRankAggreg)

glist <- list(GSE59739 = GSE59739, GSE62424 = GSE62424, GSE63576 = GSE63576, 
              GSE100175 = GSE100175,GSE131230 = GSE131230, GSE139088 = GSE139088,
              GSE154659 = GSE154659, GSE155622 = GSE155622, GSE168032 = GSE168032, PMID30096314 = PMID30096314)

mus_neuron_ma <- aggregateRanks(glist = glist)

# compare mean gene ranks to rra score
library(tidyverse)

GSE59739 <- read.csv("data/GSE59739_tpm_filtered_sorted.csv")
GSE59739$rank <- 1:nrow(GSE59739)
GSE59739 <- GSE59739[,c(2,4)]

GSE62424 <- read.csv("data/GSE62424_all_tpm_filtered_sorted.csv")
GSE62424$rank <- 1:nrow(GSE62424)
GSE62424 <- GSE62424[,c(2,4)]


GSE63576 <- read.csv("data/GSE63576_tpm_filtered_sorted.csv")
GSE63576$rank <- 1:nrow(GSE63576)
GSE63576 <- GSE63576[,c(2,4)]

GSE100175 <- read.csv("data/GSE100175_all_tpm_filtered_sorted.csv")
GSE100175$rank <- 1:nrow(GSE100175)
GSE100175 <- GSE100175[,c(2,4)]

GSE131230 <- read.csv("data/GSE131230_tpm_filtered_sorted.csv")
GSE131230$rank <- 1:nrow(GSE131230)
GSE131230 <- GSE131230[,c(2,4)]

GSE139088 <- read.csv("data/GSE139088_tpm_filtered_sorted.csv")
GSE139088$rank <- 1:nrow(GSE139088)
GSE139088 <- GSE139088[,c(2,4)]

GSE154659 <- read.csv("data/GSE154659_tpm_filtered_sorted.csv")
GSE154659$rank <- 1:nrow(GSE154659)
GSE154659 <- GSE154659[,c(2,4)]

GSE155622 <- read.csv("data/GSE155622_tpm_filtered_sorted.csv")
GSE155622$rank <- 1:nrow(GSE155622)
GSE155622 <- GSE155622[,c(2,4)]

GSE168032 <- read.csv("data/GSE168032_tpm_filtered_sorted.csv")
GSE168032$rank <- 1:nrow(GSE168032)
GSE168032 <- GSE168032[,c(2,4)]

PMID30096314 <- read.csv("data/PMID30096314_tpm_filtered_sorted.csv")
PMID30096314$rank <- 1:nrow(PMID30096314)
PMID30096314 <- PMID30096314[,c(2,4)]

df_list <- list(GSE59739, GSE62424, GSE63576, 
               GSE100175,GSE131230, GSE139088,
               GSE154659, GSE155622, GSE168032, PMID30096314) 
df2 <- df_list %>% reduce(full_join, by='SYMBOL')

df2 <- subset(df2, SYMBOL %in% mus_neuron_ma$Name)
df2$averageRank <- rowMeans(df2[,2:11], na.rm=TRUE)
df2$median = apply(df2, 1, median, na.rm=T)

df3 <- merge(df2[,c(1,12)], mus_neuron_ma, by.x='SYMBOL', by.y='Name')
df3$significance <- ifelse(df3$Score < 0.05, 'significant', 'non-significant')

ggplot(df3, aes(x=averageRank, color=significance)) +
  geom_histogram(fill="white")

ggplot(df3, aes(x=significance, y=averageRank, color=significance)) +
  geom_boxplot()

plot(df3$averageRank, df3$Score)


# ortholog mapping from mouse to human
musGenes <- mus_neuron_ma$Name

convertMouseGeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  humanx <- unique(genesV2[,c(1,2)])

  return(humanx)
}

genes <- convertMouseGeneList(musGenes)

colnames(mus_neuron_ma)[1] <- "MGI.symbol"
hs_mapped_neuron_ma <- merge(genes,mus_neuron_ma, by = "MGI.symbol")

# check for duplicate symbols
sum(duplicated(hs_mapped_neuron_ma$HGNC.symbol))
sum(duplicated(hs_mapped_neuron_ma$MGI.symbol))
## There are many to many mapping

# summarize multiple MGI.symbol mapping to same HGNC.symbol 
require(dplyr)
hs_mapped_neuron_ma <- hs_mapped_neuron_ma %>% 
  group_by(HGNC.symbol) %>% 
  mutate(MGI.symbol = paste0(MGI.symbol, collapse = ", ")) 


# taking average score
hs_mapped_neuron_ma <- data.frame(hs_mapped_neuron_ma %>% 
                                 group_by(MGI.symbol,HGNC.symbol) %>% 
                                 summarise_all(funs(mean)))


hs_mapped_neuron_ma$MGI.symbol <- NULL 

# write file
write.csv(mus_neuron_ma,file="results/MUS_Sensory_Neuron_Transcript_MA.csv",row.names = FALSE)
write.csv(hs_mapped_neuron_ma,file="results/MUS_Mapped_Sensory_Neuron_Transcript_MA.csv",row.names = FALSE)
