# set working directory 
setwd("~/Documents/GitHub/MM_Pain_Network/PC")

####################### Import packages and functions
library(readr)
library(vsn)
library(limma)
library(progeny)
library(dorothea)
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(reshape)
library(gridExtra)
library(grid)
library(cowplot)
library(ggrepel)
library(hexbin)
library(visNetwork)
library(CARNIVAL)
library(OmnipathR)
library(scales)
library(data.table)


source("scripts/support_functions.R")
source("scripts/assignPROGENyScores.r")
source("scripts/generateTFList.r")
source("scripts/carnival_visNetwork.r")
source("scripts/support_enrichment.r")
source("scripts/support_networks.r")

############################# Import pc normalized data files and prepare targets file

GSE6477 <- read.csv("data/pc/GSE6477 Normalized.Filtered(ID).Annotated.Data.csv")

GSE6477_dataset <- GSE6477[,3:119]
names(GSE6477_dataset)[1] <- "Symbol"
sum(duplicated(GSE6477_dataset$Symbol))

targets_6477 <- as.data.frame(matrix(NA,116,2))
names(targets_6477) <- c("sample","condition")
targets_6477$sample <- names(GSE6477[,-c(1:3)])
targets_6477$condition <- c("MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM", "MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM", "MM",	"MM",	"MM",	"MM",	"MM", "MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",
                            "MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM", "MM", "MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"HD" ,"MM", "MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",	
                            "MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"HD",	"MM",	"MM",	"MM",	"MM",	"MM","MM",	"MM",	"MM",	"MM",	"MM",	"HD",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",	"MM",
                            "MM",	"MM",	"HD",	"HD",	"HD",	"HD",	"HD",	"HD",	"HD",	"HD",	"HD",	"HD",	"HD"	,"HD")

###
GSE39754 <- read.csv("data/pc/GSE39754 Normalized.Filtered(ID).Annotated.Data.csv")

GSE39754_dataset <- GSE39754[,3:179]
names(GSE39754_dataset)[1] <- "Symbol"
sum(duplicated(GSE39754_dataset$Symbol))


targets_GSE39754 <- as.data.frame(matrix(NA,176,2))
names(targets_GSE39754) <- c("sample","condition")
targets_GSE39754$sample <- names(GSE39754[,-c(1:3)])
targets_GSE39754$condition <- c("HD","HD","HD","HD","HD","HD","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM")

###
GSE47552 <- read.csv("data/pc/GSE47552 Normalized.Filtered(ID).Annotated.Data.csv")

GSE47552_dataset <- GSE47552[,3:49]
names(GSE47552_dataset)[1] <- "Symbol"
sum(duplicated(GSE47552_dataset$Symbol))

targets_GSE47552 <- as.data.frame(matrix(NA,46,2))
names(targets_GSE47552) <- c("sample","condition")
targets_GSE47552$sample <- names(GSE47552[,-c(1:3)])
targets_GSE47552$condition <- c("MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","HD","HD","HD","HD","HD")

###
GSE76425 <- read.csv("data/pc/GSE76425_Filtered(Counts).Normalized.Annotated.Data_voom.csv")

GSE76425_dataset <- GSE76425[,2:13]
names(GSE76425_dataset)[1] <- "Symbol"
sum(duplicated(GSE76425_dataset$Symbol))
GSE76425_dataset <-  GSE76425_dataset[!(is.na(GSE76425_dataset$Symbol) | GSE76425_dataset$Symbol=="-"), ]
sum(duplicated(GSE76425_dataset$Symbol))

targets_GSE76425 <- as.data.frame(matrix(NA,11,2))
names(targets_GSE76425) <- c("sample","condition")
targets_GSE76425$sample <- names(GSE76425[,-c(1:2)])
targets_GSE76425$condition <- c("HD","HD","HD","HD","MM","MM","MM","MM","MM","MM","MM")

###
GSE116294  <- read.csv("data/pc/GSE116294 Normalized.Filtered(ID).Annotated.Data.csv")

GSE116294_dataset <- GSE116294[,3:57]
names(GSE116294_dataset)[1] <- "Symbol"
sum(duplicated(GSE116294_dataset$Symbol))

targets_GSE116294 <- as.data.frame(matrix(NA,54,2))
names(targets_GSE116294) <- c("sample","condition")
targets_GSE116294$sample <- names(GSE116294[,-c(1:3)])
targets_GSE116294$condition <- c("HD","HD","HD","HD","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM")
  
###
GSE125361 <- read.csv("data/pc/GSE125361 Normalized.Filtered(ID).Annotated.Data.csv")

GSE125361_dataset <- GSE125361[,3:51]
names(GSE125361_dataset)[1] <- "Symbol"
GSE125361_dataset <-  GSE125361_dataset[!(is.na(GSE125361_dataset$Symbol) | GSE125361_dataset$Symbol==""), ]
sum(duplicated(GSE125361_dataset$Symbol))

targets_GSE125361 <- as.data.frame(matrix(NA,48,2))
names(targets_GSE125361) <- c("sample","condition")
targets_GSE125361$sample <- names(GSE125361[,-c(1:3)])
targets_GSE125361$condition <- c("MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","HD","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","HD","MM","MM","MM","MM","MM","MM","MM","HD")

###
GSE153380 <- read.csv("data/pc/GSE153380_Filtered(ID,Counts).Normalized.Annotated.Data_voom.csv")

GSE153380_dataset <- GSE153380[,-c(1,3)]
names(GSE153380_dataset)[1] <- "Symbol"
sum(duplicated(GSE153380_dataset$Symbol))
GSE153380_dataset <-  GSE153380_dataset[!(is.na(GSE153380_dataset$Symbol) | GSE153380_dataset$Symbol==""), ]
sum(duplicated(GSE153380_dataset$Symbol))
GSE153380_dataset <- data.frame(GSE153380_dataset %>% 
                                  group_by(Symbol) %>% 
                                  summarise_all(funs(max)))
sum(duplicated(GSE153380_dataset$Symbol))


targets_GSE153380 <- as.data.frame(matrix(NA,31,2))
names(targets_GSE153380) <- c("sample","condition")
targets_GSE153380$sample <- names(GSE153380[,-c(1:3)])
targets_GSE153380$condition <- c("MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","HD","HD","HD")

###
GSE175384 <- read.csv("data/pc/GSE175384_Filtered(ID,Counts).Normalized.Annotated.Data_voom.csv")

GSE175384_dataset <- GSE175384 [,-1]
names(GSE175384_dataset)[1] <- "Symbol"
sum(duplicated(GSE175384_dataset$Symbol))
GSE175384_dataset <- data.frame(GSE175384_dataset %>% 
                                  group_by(Symbol) %>% 
                                  summarise_all(funs(max)))
sum(duplicated(GSE175384_dataset$Symbol))

targets_GSE175384  <- as.data.frame(matrix(NA,64,2))
names(targets_GSE175384) <- c("sample","condition")
targets_GSE175384$sample <- names(GSE175384 [,-c(1:2)])
targets_GSE175384$condition <- c("MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","HD","HD","HD","HD","HD","HD","HD","HD","HD","HD","HD","HD","HD","HD","HD","HD","HD","HD","HD","HD","HD","HD","HD","HD","HD","HD","HD","HD","HD","HD","HD","HD",
                                 "MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM","MM")



##
rm(GSE6477, GSE116294,GSE125361,GSE153380, GSE175384, GSE39754, GSE47552, GSE76425)

##
targets <- list(targets_6477,targets_GSE116294,targets_GSE125361,targets_GSE153380,targets_GSE175384,targets_GSE39754,targets_GSE47552,targets_GSE76425)
targets <- rbindlist(targets)

rm(targets_6477,targets_GSE116294,targets_GSE125361,targets_GSE153380,targets_GSE175384,targets_GSE39754,targets_GSE47552,targets_GSE76425)

normalized_counts <- Reduce(function(...) merge(..., all=TRUE, by="Symbol"), list(GSE116294_dataset,GSE125361_dataset,GSE153380_dataset,GSE175384_dataset,GSE39754_dataset,GSE47552_dataset,GSE6477_dataset,GSE76425_dataset))
rownames(normalized_counts) <- normalized_counts$Symbol
sum(duplicated(normalized_counts$Symbol))
normalized_counts$Symbol <- NULL

rm(GSE116294_dataset,GSE125361_dataset,GSE153380_dataset,GSE175384_dataset,GSE39754_dataset,GSE47552_dataset,GSE6477_dataset,GSE76425_dataset)

# load pc meta-analysis data
load("data/data_string_PC_BMSC_220613.RData")

pc_meta.analysis <- merged_analysis_PC_8DS_ID_all[,c(2,4)]
sum(duplicated(pc_meta.analysis$hgnc_symbol))
pc_meta.analysis <- data.frame(pc_meta.analysis %>% 
                              group_by(hgnc_symbol) %>% 
                              summarise_all(funs(max)))

pc_meta.analysis_symbol <- pc_meta.analysis$hgnc_symbol

# filter normalized counts file by genes in pc meta-analysis
normalized_counts <- normalized_counts[pc_meta.analysis_symbol, ] 
normalized_counts$Symbol <- rownames(normalized_counts)
sum(duplicated(normalized_counts$Symbol))

# write files
write.csv(targets,"data/pc/targets.csv", row.names = FALSE)
write.csv(normalized_counts,"data/pc/normalized_counts.csv", row.names = FALSE)
write.csv(pc_meta.analysis,"data/pc/pc_meta.analysis.csv", row.names = FALSE)

################### 02_Pathway_activity_with_Progeny

## We read the normalised counts and the experimental design 
Normalised_counts <- read.csv("data/pc/normalized_counts.csv")
Experimental_design <- read.csv("data/pc/targets.csv")

## We read the results from the differential analysis. 
ttop_KOvsWT <- read.csv("data/pc/pc_meta.analysis.csv")

Normalised_counts_matrix <- Normalised_counts %>% 
  dplyr::mutate_if(~ any(is.na(.x)),~ ifelse(is.na(.x),0,.x)) %>% 
  tibble::column_to_rownames(var = "Symbol") %>% 
  as.matrix()

ttop_KOvsWT_matrix <- ttop_KOvsWT %>% 
  dplyr::select(hgnc_symbol, effectSize) %>% 
  dplyr::filter(!is.na(effectSize)) %>% 
  column_to_rownames(var = "hgnc_symbol") %>%
  as.matrix()

## Pathway activity with Progeny

PathwayActivity_counts <- progeny(Normalised_counts_matrix, scale=TRUE, 
                                  organism="Human", top = 100)
Activity_counts <- as.vector(PathwayActivity_counts)

paletteLength <- 100
myColor <- 
  colorRampPalette(c("darkblue", "whitesmoke","indianred"))(paletteLength)

progenyBreaks <- c(seq(min(Activity_counts), 0, 
                       length.out=ceiling(paletteLength/2) + 1),
                   seq(max(Activity_counts)/paletteLength, 
                       max(Activity_counts), 
                       length.out=floor(paletteLength/2)))

progeny_hmap <- pheatmap(t(PathwayActivity_counts),fontsize=14, 
                         fontsize_row = 10, fontsize_col = 10, 
                         color=myColor, breaks = progenyBreaks, 
                         main = "PROGENy (100)", angle_col = 45,
                         treeheight_col = 0,  border_color = NA)

PathwayActivity_zscore <- progeny(ttop_KOvsWT_matrix, 
                                  scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = TRUE) %>%
  t()
colnames(PathwayActivity_zscore) <- "NES"

PathwayActivity_zscore_df <- as.data.frame(PathwayActivity_zscore) %>% 
  rownames_to_column(var = "Pathway") %>%
  dplyr::arrange(NES) %>%
  dplyr::mutate(Pathway = factor(Pathway))

ggplot(PathwayActivity_zscore_df,aes(x = reorder(Pathway, NES), y = NES)) + 
  geom_bar(aes(fill = NES), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
          element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Pathways")

prog_matrix <- getModel("Human", top=100) %>% 
  as.data.frame()  %>%
  tibble::rownames_to_column("GeneID")

ttop_KOvsWT_df <- ttop_KOvsWT_matrix %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("GeneID")

scat_plots <- progeny::progenyScatter(df = ttop_KOvsWT_df, 
                                      weight_matrix = prog_matrix, 
                                      statName = "t_values", verbose = FALSE)

plot(scat_plots[[1]]$`MAPK`) 

PathwayActivity_CARNIVALinput <- progeny(ttop_KOvsWT_matrix, 
                                         scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = FALSE) %>%
  t () %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "Pathway") 
colnames(PathwayActivity_CARNIVALinput)[2] <- "score"
write_csv(PathwayActivity_CARNIVALinput, 
          "results/PathwayActivity_CARNIVALinput.csv")

################# 04_TranscriptionFactor_activity_with_Dorothea

## Getting Started
## We read the normalised counts and the experimental design 
Normalised_counts <- read.csv("data/pc/normalized_counts.csv")
Experimental_design <- read.csv("data/pc/targets.csv")

## We read the results from the differential analysis. 
ttop_KOvsWT <- read.csv("data/pc/pc_meta.analysis.csv")

Normalised_counts_matrix <- Normalised_counts %>% 
  dplyr::mutate_if(~ any(is.na(.x)),~ ifelse(is.na(.x),0,.x)) %>% 
  tibble::column_to_rownames(var = "Symbol") %>% 
  as.matrix()

ttop_KOvsWT_matrix <- ttop_KOvsWT %>% 
  dplyr::select(hgnc_symbol, effectSize) %>% 
  dplyr::filter(!is.na(effectSize)) %>% 
  column_to_rownames(var = "hgnc_symbol") %>%
  as.matrix()

## Transcription Factor activity with DoRothEA

## We load Dorothea Regulons
data(dorothea_hs, package = "dorothea")
regulons <- dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B","C"))

tf_activities_stat <- dorothea::run_viper(ttop_KOvsWT_matrix, regulons,
                                          options =  list(minsize = 5, eset.filter = FALSE, 
                                                          cores = 1, verbose = FALSE, nes = TRUE))

tf_activities_stat_top25 <- tf_activities_stat %>%
  as.data.frame() %>% 
  rownames_to_column(var = "GeneID") %>%
  dplyr::rename(NES = "effectSize") %>%
  dplyr::top_n(25, wt = abs(NES)) %>%
  dplyr::arrange(NES) %>% 
  dplyr::mutate(GeneID = factor(GeneID))

ggplot(tf_activities_stat_top25,aes(x = reorder(GeneID, NES), y = NES)) + 
  geom_bar(aes(fill = NES), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
          element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Transcription Factors")

targets_IRF2 <- regulons$target[regulons$tf == "IRF2"]
volcano_nice(as.data.frame(ttop_KOvsWT[ttop_KOvsWT$ID %in% targets_IRF2,]), 
             FCIndex = 2, pValIndex = 5, IDIndex = 1,nlabels = 20, label = TRUE, 
             straight = FALSE) 

tf_activities_CARNIVALinput<- tf_activities_stat %>%
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "TF") 
write_csv(tf_activities_CARNIVALinput, "results/TFActivity_CARNIVALinput.csv")

tf_activities_counts <- 
  dorothea::run_viper(Normalised_counts_matrix, regulons,
                      options =  list(minsize = 5, eset.filter = FALSE, 
                                      cores = 1, verbose = FALSE, method = c("scale")))

tf_activities_counts_filter <- tf_activities_counts %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "GeneID") %>%
  dplyr::filter(GeneID %in% tf_activities_stat_top25$GeneID) %>%
  column_to_rownames(var = "GeneID") %>%
  as.matrix()
tf_activities_vector <- as.vector(tf_activities_counts_filter)

paletteLength <- 100
myColor <- 
  colorRampPalette(c("darkblue", "whitesmoke","indianred"))(paletteLength)

dorotheaBreaks <- c(seq(min(tf_activities_vector), 0, 
                        length.out=ceiling(paletteLength/2) + 1),
                    seq(max(tf_activities_vector)/paletteLength, 
                        max(tf_activities_vector), 
                        length.out=floor(paletteLength/2)))
dorothea_hmap <- pheatmap(tf_activities_counts_filter,
                          fontsize=14, fontsize_row = 8, fontsize_col = 8, 
                          color=myColor, breaks = dorotheaBreaks,
                          main = "Dorothea ABC", angle_col = 45,
                          treeheight_col = 0,  border_color = NA)


################ 05_network_reconstruction_with_CARNIVAL

## Getting Started

## We read the normalised counts and the experimental design 
tf_activities <- read_csv("results/TFActivity_CARNIVALinput.csv")
PathwayActivity <- read_csv("results/PathwayActivity_CARNIVALinput.csv")


## Getting the scaffold network from Omnipath

omniR <- import_omnipath_interactions()

# signed and directed
omnipath_sd <- omniR %>% dplyr::filter(consensus_direction == 1 &
                                         (consensus_stimulation == 1 | 
                                            consensus_inhibition == 1
                                         ))

# changing 0/1 criteria in consensus_stimulation/inhibition to -1/1
omnipath_sd$consensus_stimulation[which( omnipath_sd$consensus_stimulation == 0)] = -1
omnipath_sd$consensus_inhibition[which( omnipath_sd$consensus_inhibition == 1)] = -1
omnipath_sd$consensus_inhibition[which( omnipath_sd$consensus_inhibition == 0)] = 1

# check consistency on consensus sign and select only those in a SIF format
sif <- omnipath_sd[,c('source_genesymbol', 'consensus_stimulation', 'consensus_inhibition', 'target_genesymbol')] %>%
  dplyr::filter(consensus_stimulation==consensus_inhibition) %>%
  unique.data.frame()

sif$consensus_stimulation <- NULL
colnames(sif) <- c('source', 'interaction', 'target')

# remove complexes
sif$source <- gsub(":", "_", sif$source)
sif$target <- gsub(":", "_", sif$target)

sif <- sif[!grepl("_",sif$source),]
sif <- sif[!grepl("_",sif$target),]

#save SIF
write_tsv(sif, "results/omnipath_carnival.tsv")

## Transcription Factor and pathway activities for CARNIVAL

# dorothea for CARNIVAL
tf_activities_carnival <- data.frame(tf_activities, stringsAsFactors = F)
rownames(tf_activities_carnival) <- tf_activities$TF
tf_activities_carnival$TF <- NULL
tfList = generateTFList(tf_activities_carnival, top=50, access_idx = 1)

# progeny for CARNIVAL
load(file = system.file("progenyMembers.RData",package="CARNIVAL"))

PathwayActivity_carnival <- data.frame(PathwayActivity, stringsAsFactors = F)
rownames(PathwayActivity_carnival) <- PathwayActivity_carnival$Pathway
PathwayActivity_carnival$Pathway <- NULL
progenylist = assignPROGENyScores(progeny = t(PathwayActivity_carnival), 
                                  progenyMembers = progenyMembers, 
                                  id = "gene", 
                                  access_idx = 1)

## Running CARNIVAL
load("data/data_string_PC_BMSC_220613.RData")

# create a subset dataframe that contains DEG of MM-pc extracellular ligands 
mm_pc_ligands_ex <- merged_analysis_PC_8DS_ID_annot[merged_analysis_PC_8DS_ID_annot$ligand == TRUE & merged_analysis_PC_8DS_ID_annot$num_sources_extracell > 0, ]
mm_pc_ligands_ex <- mm_pc_ligands_ex[,c(2,4)]

mm_pc_ligands_ex <- data.frame(mm_pc_ligands_ex %>% 
                                 group_by(hgnc_symbol) %>% 
                                 summarise_all(funs(max)))

rownames(mm_pc_ligands_ex) <- mm_pc_ligands_ex$hgnc_symbol
mm_pc_ligands_ex$hgnc_symbol <- NULL
mm_pc_ligands_ex$effectSize[mm_pc_ligands_ex$effectSize < 0] <- -1
mm_pc_ligands_ex$effectSize[mm_pc_ligands_ex$effectSize > 0] <- 1

mm_pc_ligands_ex <-as.data.frame(t(mm_pc_ligands_ex))

# get initial nodes
iniMTX = base::setdiff(sif$source, sif$target)
iniciators = base::data.frame(base::matrix(data = NaN, nrow = 1, ncol = length(iniMTX)), stringsAsFactors = F)
colnames(iniciators) = iniMTX

iniciators = mm_pc_ligands_ex

# run carnival
carnival_result = runCARNIVAL( inputObj= iniciators,
                               measObj = tfList$effectSize, 
                               netObj = sif, 
                               weightObj = progenylist$score, 
                               # solverPath = "/Applications/CPLEX_Studio_Community221/cplex/bin/x86-64_osx/cplex", 
                               # solverPath = "cbc-osx/cbc", 
                               solverPath = "/Library/gurobi952/macos_universal2/bin/gurobi_cl", 
                               # solver = "cplex",
                               # solver = "cbc",
                               solver = 'gurobi',
                               timelimit=3600,
                               mipGAP=0,
                               poolrelGAP=0 )

#transform to data.frame
carnival_result$weightedSIF <- data.frame(carnival_result$weightedSIF, stringsAsFactors = F)
carnival_result$weightedSIF$Sign <- as.numeric(carnival_result$weightedSIF$Sign)
carnival_result$weightedSIF$Weight <- as.numeric(carnival_result$weightedSIF$Weight)
carnival_result$weightedSIF$Weight <- as.numeric(carnival_result$weightedSIF$Weight)

carnival_result$nodesAttributes <- data.frame(carnival_result$nodesAttributes, stringsAsFactors = F)
carnival_result$nodesAttributes$ZeroAct <- as.numeric(carnival_result$nodesAttributes$ZeroAct)
carnival_result$nodesAttributes$UpAct <- as.numeric(carnival_result$nodesAttributes$UpAct)
carnival_result$nodesAttributes$DownAct <- as.numeric(carnival_result$nodesAttributes$DownAct)
carnival_result$nodesAttributes$AvgAct <- as.numeric(carnival_result$nodesAttributes$AvgAct)

# filer weightedSIF based on weight
df_weight <- carnival_result[["weightedSIF"]]
df_weight <- df_weight[df_weight$Weight > 0,]
carnival_result[["weightedSIF"]] <- df_weight

df_attributes <- carnival_result[["nodesAttributes"]]
df_attributes <- df_attributes[df_attributes$ZeroAct == 0,]
carnival_result[["nodesAttributes"]] <- df_attributes



saveRDS(carnival_result,"results/pc_carnival_result.rds")

# visualization (not really work)
visNet = carnival_visNet(evis = carnival_result$weightedSIF,
                         nvis = carnival_result$nodesAttributes)
visSave(visNet, file = paste0('carnival_visualization_visNetwork.html'), selfcontained = TRUE)


# export to cytoscape 
SIF <- as.data.frame(carnival_result$weightedSIF)
ATT <- as.data.frame(carnival_result$nodesAttributes)
ATT <- left_join(ATT, merged_analysis_PC_8DS_ID_all[ , c("hgnc_symbol","effectSize", "fisherFDRMin")], by = c("Node"="hgnc_symbol"), all.x = TRUE)

ATT[is.na(ATT)] <- " "

write_csv(SIF, file = "results/PC_SIF.csv")
write_csv(ATT, file = "results/PC_ATT.csv")


################ 06_analysis_CARNIVAL_results

#read CARNIVAL results
carnival_result = readRDS("results/pc_carnival_result.rds")
carnival_sample_resolution = readRDS("results/carnival_sample_resolution.rds")
pkn = read_tsv("results/omnipath_carnival.tsv")

## Enrichment Analysis

# Load pathways
pathways = gmt_to_csv("data/c2.cp.v7.1.symbols.gmt")

# Extract nodes and background
nodes_carnival = extractCARNIVALnodes(carnival_result)

# Run GSA hyper Geometric test
sig_pathways <- runGSAhyper(genes = nodes_carnival$sucesses, 
                            universe = nodes_carnival$bg, gsc = loadGSC(pathways))
sig_pathways_df <- as.data.frame(sig_pathways$resTab)  %>% 
  tibble::rownames_to_column(var = "pathway") 

#data for plotting
PathwaysSelect <- sig_pathways_df %>%
  dplyr::select(pathway, `p-value`, `Adjusted p-value`) %>%
  dplyr::filter(`Adjusted p-value` <= 0.01) %>%
  dplyr::rename(pvalue = `p-value`, AdjPvalu = `Adjusted p-value`) %>% 
  dplyr::mutate(pathway = as.factor(pathway))

PathwaysSelect <- data.frame(t(apply(PathwaysSelect, 1, function(r){
  aux = unlist(strsplit( sub("_",";", r["pathway"]), ";" ))
  r["pathway"] = gsub("_", " ", aux[2])
  return(c(r, "source" = aux[1]))
})))

colnames(PathwaysSelect) = c("pathway", "pvalue", "AdjPvalu", "source")
PathwaysSelect$AdjPvalu = as.numeric(PathwaysSelect$AdjPvalu)

ggdata = PathwaysSelect %>% 
  dplyr::filter(AdjPvalu <= 0.05) %>% 
  dplyr::group_by(source) %>% 
  dplyr::arrange(AdjPvalu) %>%
  dplyr::slice(1:5)


# Visualize top results
ggplot(ggdata, aes(y = reorder(pathway, AdjPvalu), x = -log10(AdjPvalu)), color = source) + 
  geom_bar(stat = "identity") +
  facet_grid(source ~ ., scales="free_y") +
  scale_x_continuous(
    expand = c(0.01, 0.01),
    limits = c(0, ceiling(max(-log10(PathwaysSelect$AdjPvalu)))),
    breaks = seq(floor(min(-log10(PathwaysSelect$AdjPvalu))), ceiling(max(-log10(PathwaysSelect$AdjPvalu))), 1),
    labels = math_format(10^-.x)
  ) +
  annotation_logticks(sides = "bt") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.y = element_text(size = 6)) +
  ylab("")
