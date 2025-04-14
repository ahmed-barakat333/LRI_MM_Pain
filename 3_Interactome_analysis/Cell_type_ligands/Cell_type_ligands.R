# set working directory
setwd("~/Work/PhD/Projects/MoD TrD LRI MM-Pain/Computational Analysis/Ligand_Receptor/Cells_Ligands")

# load packages
library(readxl)
library(tidyverse)
library(UniProtKeywords)
library(org.Hs.eg.db)
library(qpcR)
library(reshape2)
library(data.table)
library(stringr)
library(cowplot)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(gdata)
library(scales)

#___________________________________________________________________________________________
#___________________________________________________________________________________________

# get gene GO keywords annotations

biological_process_kw = load_keyword_genesets(9606, as_table = TRUE,  category = "Biological process")
biological_process_kw <- biological_process_kw %>% 
                          group_by(gene) %>% 
                          mutate(keyword = paste0(keyword, collapse = " - ")) %>% 
                          distinct()
colnames(biological_process_kw)[1] <- "biological_process_keyword"

molecular_function_kw = load_keyword_genesets(9606, as_table = TRUE,  category = "Molecular function")
molecular_function_kw <- molecular_function_kw %>% 
                          group_by(gene) %>% 
                          mutate(keyword = paste0(keyword, collapse = " - ")) %>% 
                          distinct()
colnames(molecular_function_kw)[1] <- "molecular_function_keyword"

#___________________________________________________________________________________________
#___________________________________________________________________________________________

# load ligands data
pc_ligands <- read_excel("data/PC_meta_analysis_annotated_220408.xlsx")
pc_ligands <- subset(pc_ligands, pc_sn_nb_r == 'TRUE' & ligand == 'TRUE')
pc_ligands <- dplyr::select(pc_ligands, hgnc_symbol, entrez_id, effectSize)
#write.csv(pc_ligands,"results/PC_Ligands.csv", row.names = FALSE)
colnames(pc_ligands)[2] <- "gene"
pc_ligands$gene <- as.character(pc_ligands$gene)
pc_ligands <- merge(pc_ligands,molecular_function_kw, by="gene", all.x = TRUE)
pc_ligands <- merge(pc_ligands,biological_process_kw, by="gene", all.x = TRUE)
pc_ligands$gene <- NULL
colnames(pc_ligands)[1] <- "Plasma Cells"


bmsc_ligands <- read_excel("data/BMSC_meta_analysis_annotated_220317.xlsx")
bmsc_ligands <- subset(bmsc_ligands, bmsc_sn_nb_r == 'TRUE' & ligand == 'TRUE')
bmsc_ligands <- dplyr::select(bmsc_ligands, hgnc_symbol, entrez_id, effectSize)
#write.csv(bmsc_ligands,"results/BMSC_Ligands.csv", row.names = FALSE)
colnames(bmsc_ligands)[2] <- "gene"
bmsc_ligands$gene <- as.character(bmsc_ligands$gene)
bmsc_ligands <- merge(bmsc_ligands,molecular_function_kw, by="gene", all.x = TRUE)
bmsc_ligands <- merge(bmsc_ligands,biological_process_kw, by="gene", all.x = TRUE)
bmsc_ligands$gene <- NULL
colnames(bmsc_ligands)[1] <- "Bone Marrow Stromal Cells"

ec_ligands <- read.table("data/EC_DEG_analysis_annotated_220622.tsv", sep = '\t', header = TRUE)
ec_ligands <- subset(ec_ligands, ec_sn_nb_r == 'true' & ligand == 'true')
ec_ligands <- dplyr::select(ec_ligands, hgnc_symbol, ENTREZID, logFC)
#write.csv(ec_ligands,"results/EC_Ligands.csv", row.names = FALSE)
colnames(ec_ligands)[2] <- "gene"
ec_ligands$gene <- as.character(ec_ligands$gene)
ec_ligands <- merge(ec_ligands,molecular_function_kw, by="gene", all.x = TRUE)
ec_ligands <- merge(ec_ligands,biological_process_kw, by="gene", all.x = TRUE)
ec_ligands$gene <- NULL
colnames(ec_ligands)[1] <- "Endothelial Cells"

ocy_ligands <- read.table("data/OC_DEG_analysis_annotated_220518.tsv", sep = '\t', header = TRUE)
ocy_ligands <- subset(ocy_ligands, oc_sn_nb_r == 'true' & ligand == 'true')
ocy_ligands <- dplyr::select(ocy_ligands, hgnc_symbol, ENTREZID, logFC)
#write.csv(ocy_ligands,"results/OCY_Ligands.csv", row.names = FALSE)
colnames(ocy_ligands)[2] <- "gene"
ocy_ligands$gene <- as.character(ocy_ligands$gene)
ocy_ligands <- merge(ocy_ligands,molecular_function_kw, by="gene", all.x = TRUE)
ocy_ligands <- merge(ocy_ligands,biological_process_kw, by="gene", all.x = TRUE)
ocy_ligands$gene <- NULL
colnames(ocy_ligands)[1] <- "Osteocytes"

opc_ligands <- read.table("data/opc_DEG_analysis_annotated_221215.tsv", sep = '\t', header = TRUE)
opc_ligands <- subset(opc_ligands, sn_nb_r == 'true' & ligand == 'true')
opc_ligands <- dplyr::select(opc_ligands, hgnc_symbol, ENTREZID, logFC)
#write.csv(opc_ligands,"results/OPC_Ligands.csv", row.names = FALSE)
colnames(opc_ligands)[2] <- "gene"
opc_ligands$gene <- as.character(opc_ligands$gene)
opc_ligands <- merge(opc_ligands,molecular_function_kw, by="gene", all.x = TRUE)
opc_ligands <- merge(opc_ligands,biological_process_kw, by="gene", all.x = TRUE)
opc_ligands$gene <- NULL
colnames(opc_ligands)[1] <- "Osteogenic Precursor Cells"

adcy_ligands <- read.table("data/ac_DEG_analysis_annotated_230102.tsv", sep = '\t', header = TRUE)
adcy_ligands <- subset(adcy_ligands, sn_nb_r == 'true' & ligand == 'true')
adcy_ligands$ENTREZID <- mapIds(org.Hs.eg.db, adcy_ligands$hgnc_symbol,  'ENTREZID','SYMBOL')
#adcy_ligands$ENTREZID <- adcy_ligands$ENTREZID$ENTREZID
adcy_ligands <- dplyr::select(adcy_ligands,hgnc_symbol, ENTREZID, logFC)
#write.csv(adcy_ligands,"results/ADCY_Ligands.csv", row.names = FALSE)
colnames(adcy_ligands)[2] <- "gene"
adcy_ligands$gene <- as.character(adcy_ligands$gene)
adcy_ligands <- merge(adcy_ligands,molecular_function_kw, by="gene", all.x = TRUE)
adcy_ligands <- merge(adcy_ligands,biological_process_kw, by="gene", all.x = TRUE)
adcy_ligands$gene <- NULL
colnames(adcy_ligands)[1] <- "Adipocytes"

nkc_ligands <- read.table("data/NKC_DEG_analysis_annotated_220626.tsv", sep = '\t', header = TRUE)
nkc_ligands <- subset(nkc_ligands, nkc_sn_nb_r == 'true' & ligand == 'true')
nkc_ligands <- dplyr::select(nkc_ligands, hgnc_symbol, ENTREZID, logFC)
#write.csv(nkc_ligands,"results/NKC_Ligands.csv", row.names = FALSE)
colnames(nkc_ligands)[2] <- "gene"
nkc_ligands$gene <- as.character(nkc_ligands$gene)
nkc_ligands <- merge(nkc_ligands,molecular_function_kw, by="gene", all.x = TRUE)
nkc_ligands <- merge(nkc_ligands,biological_process_kw, by="gene", all.x = TRUE)
nkc_ligands$gene <- NULL
colnames(nkc_ligands)[1] <- "Natural Killer Cells"

neut_ligands <- read.table("data/Neutrophil_DEG_analysis_annotated_230102.tsv", sep = '\t', header = TRUE)
neut_ligands <- subset(neut_ligands, sn_nb_r == 'true' & ligand == 'true')
neut_ligands <- dplyr::select(neut_ligands,hgnc_symbol,entrezgene_id, logFC)
#write.csv(neut_ligands,"results/NEUT_Ligands.csv", row.names = FALSE)
colnames(neut_ligands)[2] <- "gene"
neut_ligands$gene <- as.character(neut_ligands$gene)
neut_ligands <- merge(neut_ligands,molecular_function_kw, by="gene", all.x = TRUE)
neut_ligands <- merge(neut_ligands,biological_process_kw, by="gene", all.x = TRUE)
neut_ligands$gene <- NULL
colnames(neut_ligands)[1] <- "Neutrophils"

mp_ligands <- read.table("data/MPH_DEG_analysis_annotated_230104.tsv", sep = '\t', header = TRUE)
mp_ligands <- subset(mp_ligands, sn_nb_r == 'true' & ligand == 'true')
mp_ligands <- dplyr::select(mp_ligands,hgnc_symbol,ENTREZID.Hs, effectSize)
#write.csv(mp_ligands,"results/MP_Ligands.csv", row.names = FALSE)
colnames(mp_ligands)[2] <- "gene"
mp_ligands$gene <- as.character(mp_ligands$gene)
mp_ligands <- merge(mp_ligands,molecular_function_kw, by="gene", all.x = TRUE)
mp_ligands <- merge(mp_ligands,biological_process_kw, by="gene", all.x = TRUE)
mp_ligands$gene <- NULL
colnames(mp_ligands)[1] <- "Macrophages"

treg_ligands <- read.table("data/rtc_GSE109533_DEG_limma_uniprot_annot_230102.tsv", sep = '\t', header = TRUE)
treg_ligands <- subset(treg_ligands, lrpair == 'TRUE' & ligand == 'TRUE')
treg_ligands <- dplyr::select(treg_ligands,hgnc_symbol,ENTREZID.Hs, logFC)
#write.csv(treg_ligands,"results/Treg_Ligands.csv", row.names = FALSE)
colnames(treg_ligands)[2] <- "gene"
treg_ligands$gene <- as.character(treg_ligands$gene)
treg_ligands <- merge(treg_ligands,molecular_function_kw, by="gene", all.x = TRUE)
treg_ligands <- merge(treg_ligands,biological_process_kw, by="gene", all.x = TRUE)
treg_ligands$gene <- NULL
colnames(treg_ligands)[1] <- "Regulatory T Cells"

hsc_ligands <- read.table("data/hst_DEG_analysis_annotated_221215.tsv", sep = '\t', header = TRUE)
hsc_ligands <- subset(hsc_ligands, sn_nb_r == 'true' & ligand == 'true')
hsc_ligands <- dplyr::select(hsc_ligands,hgnc_symbol, ENTREZID, logFC)
hsc_ligands <- hsc_ligands[c(1,2),]
#write.csv(hsc_ligands,"results/HSC_Ligands.csv", row.names = FALSE)
colnames(hsc_ligands)[2] <- "gene"
hsc_ligands$gene <- as.character(hsc_ligands$gene)
hsc_ligands <- merge(hsc_ligands,molecular_function_kw, by="gene", all.x = TRUE)
hsc_ligands <- merge(hsc_ligands,biological_process_kw, by="gene", all.x = TRUE)
hsc_ligands$gene <- NULL
colnames(hsc_ligands)[1] <- "Hematopoietic Stem Cells"

wbm_ligands <- read.table("data/wbm_DEG_analysis_annotated_230109.tsv", sep = '\t', header = TRUE)
wbm_ligands <- subset(wbm_ligands, sn_nb_r == 'true' & ligand == 'true')
wbm_ligands <- dplyr::select(wbm_ligands, hgnc_symbol, ENTREZ_GENE_ID.x, logFC)
#write.csv(wbm_ligands,"results/WBM_Ligands.csv", row.names = FALSE)
colnames(wbm_ligands)[2] <- "gene"
wbm_ligands$gene <- as.character(wbm_ligands$gene)
wbm_ligands <- merge(wbm_ligands,molecular_function_kw, by="gene", all.x = TRUE)
wbm_ligands <- merge(wbm_ligands,biological_process_kw, by="gene", all.x = TRUE)
wbm_ligands$gene <- NULL
colnames(wbm_ligands)[1] <- "Whole Bone Marrow"

#___________________________________________________________________________________________
#___________________________________________________________________________________________

# prepare dataset for ggplot

ligands <- data.frame(qpcR:::cbind.na(pc_ligands$`Plasma Cells`,bmsc_ligands$`Bone Marrow Stromal Cells`,ec_ligands$`Endothelial Cells`,ocy_ligands$Osteocytes,
                           opc_ligands$`Osteogenic Precursor Cells`, adcy_ligands$Adipocytes,nkc_ligands$`Natural Killer Cells`, neut_ligands$Neutrophils,mp_ligands$Macrophages,
                           treg_ligands$`Regulatory T Cells`,hsc_ligands$`Hematopoietic Stem Cells`,wbm_ligands$`Whole Bone Marrow`))
colnames(ligands) <- c("Plasma Cells","Bone Marrow Stromal Cells", "Endothelial Cells", "Osteocytes","Osteogenic Precursor Cells",
                       "Adipocytes", "Natural Killer Cells", "Neutrophils","Macrophages","Regulatory T Cells", "Hematopoietic Stem Cells","Whole Bone Marrow")

write.csv(ligands,"results/MM_Cells_Ligands_v1.csv", row.names = FALSE)


ligands_no <- length(Reduce(union, list(pc_ligands$`Plasma Cells`,bmsc_ligands$`Bone Marrow Stromal Cells`,ec_ligands$`Endothelial Cells`,ocy_ligands$Osteocytes,
                                        opc_ligands$`Osteogenic Precursor Cells`, adcy_ligands$Adipocytes,nkc_ligands$`Natural Killer Cells`, neut_ligands$Neutrophils,mp_ligands$Macrophages,
                                        treg_ligands$`Regulatory T Cells`,hsc_ligands$`Hematopoietic Stem Cells`,wbm_ligands$`Whole Bone Marrow`)))

ligands <- Reduce(union, list(pc_ligands$`Plasma Cells`,bmsc_ligands$`Bone Marrow Stromal Cells`,ec_ligands$`Endothelial Cells`,ocy_ligands$Osteocytes,
                              opc_ligands$`Osteogenic Precursor Cells`, adcy_ligands$Adipocytes,nkc_ligands$`Natural Killer Cells`, neut_ligands$Neutrophils,mp_ligands$Macrophages,
                              treg_ligands$`Regulatory T Cells`,hsc_ligands$`Hematopoietic Stem Cells`,wbm_ligands$`Whole Bone Marrow`))

write.csv(ligands,"results/MM_Ligands.csv", row.names = FALSE)

cells_ligands <- data.frame(matrix(NA, nrow = 348, ncol = 12))

colnames(cells_ligands) <- c("Plasma Cells","Bone Marrow Stromal Cells", "Endothelial Cells", "Osteocytes","Osteogenic Precursor Cells",
                             "Adipocytes", "Natural Killer Cells", "Neutrophils","Macrophages","Regulatory T Cells", "Hematopoietic Stem Cells","Whole Bone Marrow")

rownames(cells_ligands) <- Reduce(union, list(pc_ligands$`Plasma Cells`,bmsc_ligands$`Bone Marrow Stromal Cells`,ec_ligands$`Endothelial Cells`,ocy_ligands$Osteocytes,
                                              opc_ligands$`Osteogenic Precursor Cells`, adcy_ligands$Adipocytes,nkc_ligands$`Natural Killer Cells`, neut_ligands$Neutrophils,mp_ligands$Macrophages,
                                              treg_ligands$`Regulatory T Cells`,hsc_ligands$`Hematopoietic Stem Cells`,wbm_ligands$`Whole Bone Marrow`))

cells_ligands$ligand <- rownames(cells_ligands)

cells_ligands$`Plasma Cells` <- ifelse(cells_ligands$ligand %in% pc_ligands$`Plasma Cells`, 1, NA)
cells_ligands$`Bone Marrow Stromal Cells` <- ifelse(cells_ligands$ligand %in% bmsc_ligands$`Bone Marrow Stromal Cells`, 1, NA)
cells_ligands$`Endothelial Cells` <- ifelse( cells_ligands$ligand %in% ec_ligands$`Endothelial Cells`, 1, NA)
cells_ligands$Osteocytes <- ifelse( cells_ligands$ligand %in% ocy_ligands$Osteocytes, 1, NA)
cells_ligands$`Osteogenic Precursor Cells` <- ifelse( cells_ligands$ligand %in% opc_ligands$`Osteogenic Precursor Cells`, 1, NA)
cells_ligands$Adipocytes <- ifelse( cells_ligands$ligand %in% adcy_ligands$Adipocytes, 1, NA)
cells_ligands$`Natural Killer Cells` <- ifelse( cells_ligands$ligand %in% nkc_ligands$`Natural Killer Cells`, 1, NA)
cells_ligands$Neutrophils <- ifelse( cells_ligands$ligand %in% neut_ligands$Neutrophils, 1, NA)
cells_ligands$Macrophages <- ifelse( cells_ligands$ligand %in% mp_ligands$Macrophages, 1, NA)
cells_ligands$`Regulatory T Cells` <- ifelse( cells_ligands$ligand %in% treg_ligands$`Regulatory T Cells`, 1, NA)
cells_ligands$`Hematopoietic Stem Cells` <- ifelse( cells_ligands$ligand %in% hsc_ligands$`Hematopoietic Stem Cells`, 1, NA)
cells_ligands$`Whole Bone Marrow` <- ifelse( cells_ligands$ligand %in% wbm_ligands$`Whole Bone Marrow`, 1, NA)
write.csv(cells_ligands,"results/MM_Cells_Ligands_v2.csv", row.names = TRUE)


cells_ligands$`Plasma Cells` <- ifelse(cells_ligands$ligand %in% pc_ligands$`Plasma Cells`, pc_ligands$effectSize[match(cells_ligands$ligand, pc_ligands$`Plasma Cells`)], NA)
cells_ligands$`Bone Marrow Stromal Cells` <- ifelse(cells_ligands$ligand %in% bmsc_ligands$`Bone Marrow Stromal Cells`, bmsc_ligands$effectSize[match(cells_ligands$ligand, bmsc_ligands$`Bone Marrow Stromal Cells`)], NA)
cells_ligands$`Endothelial Cells` <- ifelse( cells_ligands$ligand %in% ec_ligands$`Endothelial Cells`, ec_ligands$logFC[match(cells_ligands$ligand,ec_ligands$`Endothelial Cells`)], NA)
cells_ligands$Osteocytes <- ifelse( cells_ligands$ligand %in% ocy_ligands$Osteocytes, ocy_ligands$logFC[match(cells_ligands$ligand,ocy_ligands$Osteocytes)], NA)
cells_ligands$`Osteogenic Precursor Cells` <- ifelse( cells_ligands$ligand %in% opc_ligands$`Osteogenic Precursor Cells`, opc_ligands$logFC[match(cells_ligands$ligand,opc_ligands$`Osteogenic Precursor Cells`)], NA)
cells_ligands$Adipocytes <- ifelse( cells_ligands$ligand %in% adcy_ligands$Adipocytes, adcy_ligands$logFC[match(cells_ligands$ligand,adcy_ligands$Adipocytes)], NA)
cells_ligands$`Natural Killer Cells` <- ifelse( cells_ligands$ligand %in% nkc_ligands$`Natural Killer Cells`, nkc_ligands$logFC[match(cells_ligands$ligand,nkc_ligands$`Natural Killer Cells`)], NA)
cells_ligands$Neutrophils <- ifelse( cells_ligands$ligand %in% neut_ligands$Neutrophils, neut_ligands$logFC[match(cells_ligands$ligand,neut_ligands$Neutrophils)], NA)
cells_ligands$Macrophages <- ifelse( cells_ligands$ligand %in% mp_ligands$Macrophages, mp_ligands$effectSize[match(cells_ligands$ligand,mp_ligands$Macrophages)], NA)
cells_ligands$`Regulatory T Cells` <- ifelse( cells_ligands$ligand %in% treg_ligands$`Regulatory T Cells`, treg_ligands$logFC[match(cells_ligands$ligand,treg_ligands$`Regulatory T Cells`)], NA)
cells_ligands$`Hematopoietic Stem Cells` <- ifelse( cells_ligands$ligand %in% hsc_ligands$`Hematopoietic Stem Cells`, hsc_ligands$logFC[match(cells_ligands$ligand,hsc_ligands$`Hematopoietic Stem Cells`)], NA)
cells_ligands$`Whole Bone Marrow` <- ifelse( cells_ligands$ligand %in% wbm_ligands$`Whole Bone Marrow`, wbm_ligands$logFC[match(cells_ligands$ligand,wbm_ligands$`Whole Bone Marrow`)], NA)


cells_ligands$sum <- rowSums(!is.na(cells_ligands[,c(1:12)]))

cells_ligands <- cells_ligands[order(-cells_ligands$sum),]

cells_ligands$ligand <- reorder.factor(cells_ligands$ligand, new.order=cells_ligands$ligand)
cells_ligands <- cells_ligands %>%
  arrange(ligand)

cells_ligands$sum <- NULL
cells_ligands$ligand <- as.character(cells_ligands$ligand)

cells_ligands_long <- melt(setDT(cells_ligands), id.vars = c("ligand"), variable.name = "Cells")
cells_ligands_long <- cells_ligands_long[complete.cases(cells_ligands_long), ]
cells_ligands_long$gene <- mapIds(org.Hs.eg.db, cells_ligands_long$ligand,  'ENTREZID','SYMBOL')
cells_ligands_long <- merge(cells_ligands_long,molecular_function_kw, by="gene", all.x = TRUE)
cells_ligands_long <- merge(cells_ligands_long,biological_process_kw, by="gene", all.x = TRUE)

cells_ligands_long <- cells_ligands_long[order(match(cells_ligands_long$ligand, cells_ligands$ligand)),]


tiff("results/cells_ligands.tiff", units="in", width=8, height=50, res=100)

cells_ligands_plot <- ggplot(cells_ligands_long, aes(x= Cells, y = ligand, color=value)) +
                        geom_point(size=3,shape =15) +
                        scale_color_gradientn(colors = c("blue", "white", "red"),limits = c(-5, 5),breaks = c(-5, 0,5), labels = c("-5 or less","0","5 or more"),oob=squish) +
                        theme_classic() +
                        xlab("") +
                        ylab("") +
                        #theme(plot.title = element_text(face="bold")) +
                        #ggtitle("") +
                        theme(axis.title = element_text(size=12)) + 
                        theme(axis.text.x = element_text( color="#000000", size=10, angle=90, hjust = 0.95, vjust = 0.4),
                              axis.text.y = element_text(color="#000000", size=10)) +
                        ylim(rev(unique(unlist(as.vector(cells_ligands_long[,2])))[1:50]))+
                        theme(legend.position="bottom") +
                        guides(colour = guide_colorbar(title.position = "top"))+
                        labs(color="Effect Size/ Fold Change") + 
                        theme(legend.title = element_text(size = 15),
                              legend.key.size = unit(0.5, "cm"),
                              legend.key.width = unit(1,"cm"),
                              legend.text=element_text(size=15))
  

cells_ligands_plot                      
dev.off()

cells_ligands_plot_1_hor <- ggplot(cells_ligands_long, aes(x= ligand , y = Cells, color=value)) +
  geom_point(size=3,shape =15) +
  scale_color_gradientn(colors = c("blue", "white", "red"),limits = c(-5, 5),breaks = c(-5, 0,5), labels = c("-5 or less","0","5 or more"),oob=squish) +
  theme_classic() +
  xlab("Ligand name") +
  ylab("Cell type") +
  #theme(plot.title = element_text(face="bold")) +
  #ggtitle("") +
  theme(axis.title = element_text(size=12)) + 
  theme(axis.text.x = element_text( color="#000000", size=10, angle=90, hjust = 0.95, vjust = 0.4),
        axis.text.y = element_text(color="#000000", size=10)) +
  xlim(rev(unique(unlist(as.vector(cells_ligands_long[,2])))[50:1])) +
  theme(legend.position="none") +
  scale_y_discrete(limits=rev) 

cells_ligands_plot_1_hor   

cells_ligands_plot_1 <- ggplot(cells_ligands_long, aes(x= Cells, y = ligand, color=value)) +
  geom_point(size=3,shape =15) +
  scale_color_gradientn(colors = c("blue", "white", "red"),limits = c(-5, 5),breaks = c(-5, 0,5), labels = c("-5 or less","0","5 or more"),oob=squish) +
  theme_classic() +
  xlab("") +
  ylab("") +
  #theme(plot.title = element_text(face="bold")) +
  #ggtitle("") +
  theme(axis.title = element_text(size=12)) + 
  theme(axis.text.x = element_text( color="#000000", size=10, angle=90, hjust = 0.95, vjust = 0.4),
        axis.text.y = element_text(color="#000000", size=10)) +
  ylim(rev(unique(unlist(as.vector(cells_ligands_long[,2])))[1:50])) +
  theme(legend.position="none") 

cells_ligands_plot_1   

cells_ligands_plot_2 <- ggplot(cells_ligands_long, aes(x= Cells, y = ligand, color=value)) +
                        geom_point(size=3,shape =15) +
                        scale_color_gradientn(colors = c("blue", "white", "red"),limits = c(-5, 5),breaks = c(-5, 0,5), labels = c("-5 or less","0","5 or more"),oob=squish) +
                        theme_classic() +
                        xlab("") +
                        ylab("") +
                        #theme(plot.title = element_text(face="bold")) +
                        #ggtitle("") +
                        theme(axis.title = element_text(size=12)) + 
                        theme(axis.text.x = element_text( color="#000000", size=10, angle=90, hjust = 0.95, vjust = 0.4),
                              axis.text.y = element_text(color="#000000", size=10)) +
                        ylim(rev(unique(unlist(as.vector(cells_ligands_long[,2])))[51:100])) +
                        theme(legend.position="none") 
  
cells_ligands_plot_2   

cells_ligands_plot_3 <- ggplot(cells_ligands_long, aes(x= Cells, y = ligand, color=value)) +
                        geom_point(size=3,shape =15) +
                        scale_color_gradientn(colors = c("blue", "white", "red"),limits = c(-5, 5),breaks = c(-5, 0,5), labels = c("-5 or less","0","5 or more"),oob=squish) +
                        theme_classic() +
                        xlab("") +
                        ylab("") +
                        #theme(plot.title = element_text(face="bold")) +
                        #ggtitle("") +
                        theme(axis.title = element_text(size=12)) + 
                        theme(axis.text.x = element_text( color="#000000", size=10, angle=90, hjust = 0.95, vjust = 0.4),
                              axis.text.y = element_text(color="#000000", size=10)) +
                        ylim(rev(unique(unlist(as.vector(cells_ligands_long[,2])))[101:150])) +
                        theme(legend.position="none") 

cells_ligands_plot_3       

cells_ligands_plot_4 <- ggplot(cells_ligands_long, aes(x= Cells, y = ligand, color=value)) +
                        geom_point(size=3,shape =15) +
                        scale_color_gradientn(colors = c("blue", "white", "red"),limits = c(-5, 5),breaks = c(-5, 0,5), labels = c("-5 or less","0","5 or more"),oob=squish) +
                        theme_classic() +
                        xlab("") +
                        ylab("") +
                        #theme(plot.title = element_text(face="bold")) +
                        #ggtitle("") +
                        theme(axis.title = element_text(size=12)) + 
                        theme(axis.text.x = element_text( color="#000000", size=10, angle=90, hjust = 0.95, vjust = 0.4),
                              axis.text.y = element_text(color="#000000", size=10)) +
                        ylim(rev(unique(unlist(as.vector(cells_ligands_long[,2])))[151:200])) +
                        theme(legend.position="none") 

cells_ligands_plot_4       

cells_ligands_plot_5 <- ggplot(cells_ligands_long, aes(x= Cells, y = ligand, color=value)) +
                        geom_point(size=3,shape =15) +
                        scale_color_gradientn(colors = c("blue", "white", "red"),limits = c(-5, 5),breaks = c(-5, 0,5), labels = c("-5 or less","0","5 or more"),oob=squish) +
                        theme_classic() +
                        xlab("") +
                        ylab("") +
                        #theme(plot.title = element_text(face="bold")) +
                        #ggtitle("") +
                        theme(axis.title = element_text(size=12)) + 
                        theme(axis.text.x = element_text( color="#000000", size=10, angle=90, hjust = 0.95, vjust = 0.4),
                              axis.text.y = element_text(color="#000000", size=10)) +
                        ylim(rev(unique(unlist(as.vector(cells_ligands_long[,2])))[201:250])) +
                        theme(legend.position="none") 

cells_ligands_plot_5       

cells_ligands_plot_6 <- ggplot(cells_ligands_long, aes(x= Cells, y = ligand, color=value)) +
                        geom_point(size=3,shape =15) +
                        scale_color_gradientn(colors = c("blue", "white", "red"),limits = c(-5, 5),breaks = c(-5, 0,5), labels = c("-5 or less","0","5 or more"),oob=squish) +
                        theme_classic() +
                        xlab("") +
                        ylab("") +
                        #theme(plot.title = element_text(face="bold")) +
                        #ggtitle("") +
                        theme(axis.title = element_text(size=12)) + 
                        theme(axis.text.x = element_text( color="#000000", size=10, angle=90, hjust = 0.95, vjust = 0.4),
                              axis.text.y = element_text(color="#000000", size=10)) +
                        ylim(rev(unique(unlist(as.vector(cells_ligands_long[,2])))[251:300])) +
                        theme(legend.position="none") 

cells_ligands_plot_6       

cells_ligands_plot_7 <- ggplot(cells_ligands_long, aes(x= Cells, y = ligand, color=value)) +
  geom_point(size=3,shape =15) +
  scale_color_gradientn(colors = c("blue", "white", "red"),limits = c(-5, 5),breaks = c(-5, 0,5), labels = c("-5 or less","0","5 or more"),oob=squish) +
  theme_classic() +
  xlab("") +
  ylab("") +
  #theme(plot.title = element_text(face="bold")) +
  #ggtitle("") +
  theme(axis.title = element_text(size=12)) + 
  theme(axis.text.x = element_text( color="#000000", size=10, angle=90, hjust = 0.95, vjust = 0.4),
        axis.text.y = element_text(color="#000000", size=10)) +
  ylim(rev(unique(unlist(as.vector(cells_ligands_long[,2])))[301:348])) +
  theme(legend.position="none") 

cells_ligands_plot_7       
#________________________________________________________________________________
#________________________________________________________________________________

# Draw Only FC Legend without plot

tiff("results/fc_legends.tiff", units="in", width=6, height=3, res=300)
legend <- get_legend(cells_ligands_plot_1)                    

grid.newpage()                              

grid.draw(legend) 
dev.off()
#___________________________________________________________________________________________
#___________________________________________________________________________________________
# prepare molecular function annotation 

cells_ligands_long_mf <- cells_ligands_long[,c(2,5)]
cells_ligands_long_mf$x <- "Molecular function keyword"
cells_ligands_long_mf$molecular_function_keyword <- gsub("\\Cytokine.*", "Cytokine", cells_ligands_long_mf$molecular_function_keyword)
cells_ligands_long_mf$molecular_function_keyword <- gsub("\\Growth factor.*", "Growth factor", cells_ligands_long_mf$molecular_function_keyword)
cells_ligands_long_mf$molecular_function_keyword <- gsub("\\Hydrolase.*", "Enzyme", cells_ligands_long_mf$molecular_function_keyword)
cells_ligands_long_mf$molecular_function_keyword <- gsub("\\Acyltransferase.*", "Enzyme", cells_ligands_long_mf$molecular_function_keyword)
cells_ligands_long_mf$molecular_function_keyword <- gsub("\\Kinase.*", "Enzyme", cells_ligands_long_mf$molecular_function_keyword)
cells_ligands_long_mf$molecular_function_keyword <- gsub("\\Dioxygenase.*", "Enzyme", cells_ligands_long_mf$molecular_function_keyword)
cells_ligands_long_mf$molecular_function_keyword <- gsub("\\Allosteric.*", "Enzyme", cells_ligands_long_mf$molecular_function_keyword)
cells_ligands_long_mf$molecular_function_keyword <- gsub("\\Glycosidase.*", "Enzyme", cells_ligands_long_mf$molecular_function_keyword)
cells_ligands_long_mf$molecular_function_keyword <- gsub("\\Protease inhibitor.*", "Enzyme inhibitor", cells_ligands_long_mf$molecular_function_keyword)
cells_ligands_long_mf$molecular_function_keyword <- gsub("\\Metalloenzyme inhibitor.*", "Enzyme inhibitor", cells_ligands_long_mf$molecular_function_keyword)
cells_ligands_long_mf$molecular_function_keyword <- gsub("\\Heparin-binding - Enzyme inhibitor.*", "Enzyme inhibitor", cells_ligands_long_mf$molecular_function_keyword)
cells_ligands_long_mf$molecular_function_keyword <- gsub("\\Endorphin.*", "Neuropeptide", cells_ligands_long_mf$molecular_function_keyword)
cells_ligands_long_mf$molecular_function_keyword <- gsub("\\Activator.*", "Miscellaneous", cells_ligands_long_mf$molecular_function_keyword)
cells_ligands_long_mf$molecular_function_keyword <- gsub("\\Blood.*", "Miscellaneous", cells_ligands_long_mf$molecular_function_keyword)
cells_ligands_long_mf$molecular_function_keyword <- gsub("\\Vasoactive.*", "Miscellaneous", cells_ligands_long_mf$molecular_function_keyword)
cells_ligands_long_mf$molecular_function_keyword <- gsub("\\Host cell.*", "Miscellaneous", cells_ligands_long_mf$molecular_function_keyword)
cells_ligands_long_mf$molecular_function_keyword <- gsub("\\Heparin-binding.*", "Miscellaneous", cells_ligands_long_mf$molecular_function_keyword)
cells_ligands_long_mf$molecular_function_keyword <- gsub("\\G-protein coupled receptor.*", "Miscellaneous", cells_ligands_long_mf$molecular_function_keyword)
cells_ligands_long_mf$molecular_function_keyword <- gsub("\\Chaperone.*", "Miscellaneous", cells_ligands_long_mf$molecular_function_keyword)
cells_ligands_long_mf$molecular_function_keyword <- gsub("\\Guanine-nucleotide.*", "Miscellaneous", cells_ligands_long_mf$molecular_function_keyword)
cells_ligands_long_mf$molecular_function_keyword <- gsub("\\Antioxidant.*", "Miscellaneous", cells_ligands_long_mf$molecular_function_keyword)
cells_ligands_long_mf$molecular_function_keyword <- gsub("\\Ribonucleoprotein.*", "Miscellaneous", cells_ligands_long_mf$molecular_function_keyword)
cells_ligands_long_mf$molecular_function_keyword <- gsub("\\Calcium channel.*", "Miscellaneous", cells_ligands_long_mf$molecular_function_keyword)
cells_ligands_long_mf$molecular_function_keyword <- str_replace(cells_ligands_long_mf$molecular_function_keyword, "Phospholipase A2 inhibitor", "Enzyme")
cells_ligands_long_mf$molecular_function_keyword <- str_replace(cells_ligands_long_mf$molecular_function_keyword, "Decarboxylase - Lyase", "Enzyme")
cells_ligands_long_mf$molecular_function_keyword <- str_replace(cells_ligands_long_mf$molecular_function_keyword, "Neuropeptide - Neurotransmitter - Opioid peptide", "Neuropeptide")
cells_ligands_long_mf$molecular_function_keyword <- str_replace(cells_ligands_long_mf$molecular_function_keyword, "BMiscellaneous", "Miscellaneous")


tiff("results/mf.tiff", units="in", width=7, height=50, res=100)

mf_plot <- ggplot(cells_ligands_long_mf, aes(y = ligand, x=x,color=molecular_function_keyword)) +
            geom_point(shape =15, size=3) +
            scale_color_brewer(palette="PuOr", na.value="grey") + 
            theme_classic() +
            theme(legend.position="bottom") +
            theme(legend.direction="vertical") +
            ylim(rev(unique(unlist(as.vector(cells_ligands_long[,2])))[1:50]))+
            theme(legend.text = element_text(size=15))+
            theme(legend.title = element_text(size=15)) +
            theme(legend.position="bottom")+ 
            theme(axis.text.x = element_text(color="#000000", size=12, angle = 90,vjust=.5, hjust=1),
                  axis.text.y=element_blank(),
                  axis.title.x=element_blank(),
                  axis.title.y=element_blank(),legend.position="bottom", axis.ticks.y = element_blank(),
                  panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank(),plot.background=element_blank(), axis.line=element_blank()) +                                        
            guides(color = guide_legend(override.aes = list(size = 5),title="Molecular function keyword"))



mf_plot
dev.off()

#transpose data frame
cells_ligands_long_mf_t <- transpose(cells_ligands_long_mf)

#redefine row and column names
rownames(df_t) <- colnames(df)

mf_plot_1_hor <- ggplot(cells_ligands_long_mf, aes(y = ligand, x=x,color=molecular_function_keyword)) +
            geom_point(shape =15, size=3) +
            scale_color_brewer(palette="PuOr", na.value="grey") + 
            theme_classic() +
            theme(legend.position="bottom") +
            guides(color=guide_legend(title="Molecular function keyword")) +
            theme(legend.direction="vertical") +
            ylim(rev(unique(unlist(as.vector(cells_ligands_long[,2])))[50:1]))+
            theme(legend.text = element_text(size=15))+
            theme(legend.title = element_text(size=15)) +
            theme(legend.position="none")+ 
            theme(axis.text.y = element_text(color="#000000", size=12, angle = 0,vjust=.5, hjust=1),
                  axis.text.x=element_blank(),
                  axis.title.y=element_blank(),
                  axis.title.x=element_blank(),legend.position="none", axis.ticks.x = element_blank(),
                  panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank(),plot.background=element_blank(), axis.line=element_blank()) + coord_flip()



mf_plot_1_hor

mf_plot_1 <- ggplot(cells_ligands_long_mf, aes(y = ligand, x=x,color=molecular_function_keyword)) +
  geom_point(shape =15, size=3) +
  scale_color_brewer(palette="PuOr", na.value="grey") + 
  theme_classic() +
  theme(legend.position="bottom") +
  guides(color=guide_legend(title="Molecular function keyword")) +
  theme(legend.direction="vertical") +
  ylim(rev(unique(unlist(as.vector(cells_ligands_long[,2])))[1:50]))+
  theme(legend.text = element_text(size=15))+
  theme(legend.title = element_text(size=15)) +
  theme(legend.position="none")+ 
  theme(axis.text.x = element_text(color="#000000", size=12, angle = 90,vjust=.5, hjust=1),
        axis.text.y=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none", axis.ticks.y = element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank(), axis.line=element_blank())



mf_plot_1

mf_plot_2 <- ggplot(cells_ligands_long_mf, aes(y = ligand, x=x,color=molecular_function_keyword)) +
            geom_point(shape =15, size=3) +
            scale_color_brewer(palette="PuOr", na.value="grey") + 
            theme_classic() +
            theme(legend.position="bottom") +
            guides(color=guide_legend(title="Molecular function keyword")) +
            theme(legend.direction="vertical") +
            ylim(rev(unique(unlist(as.vector(cells_ligands_long[,2])))[51:100]))+
            theme(legend.text = element_text(size=15))+
            theme(legend.title = element_text(size=15)) +
            theme(legend.position="none")+ 
            theme(axis.text.x = element_text(color="#000000", size=12, angle = 90,vjust=.5, hjust=1),
                  axis.text.y=element_blank(),
                  axis.title.x=element_blank(),
                  axis.title.y=element_blank(),legend.position="none", axis.ticks.y = element_blank(),
                  panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank(),plot.background=element_blank(), axis.line=element_blank())



mf_plot_2

mf_plot_3 <- ggplot(cells_ligands_long_mf, aes(y = ligand, x=x,color=molecular_function_keyword)) +
            geom_point(shape =15, size=3) +
            scale_color_brewer(palette="PuOr", na.value="grey") + 
            theme_classic() +
            theme(legend.position="bottom") +
            guides(color=guide_legend(title="Molecular function keyword")) +
            theme(legend.direction="vertical") +
            ylim(rev(unique(unlist(as.vector(cells_ligands_long[,2])))[101:150]))+
            theme(legend.text = element_text(size=15))+
            theme(legend.title = element_text(size=15)) +
            theme(legend.position="none")+ 
            theme(axis.text.x = element_text(color="#000000", size=12, angle = 90,vjust=.5, hjust=1),
                  axis.text.y=element_blank(),
                  axis.title.x=element_blank(),
                  axis.title.y=element_blank(),legend.position="none", axis.ticks.y = element_blank(),
                  panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank(),plot.background=element_blank(), axis.line=element_blank())



mf_plot_3

mf_plot_4 <- ggplot(cells_ligands_long_mf, aes(y = ligand, x=x,color=molecular_function_keyword)) +
            geom_point(shape =15, size=3) +
            scale_color_brewer(palette="PuOr", na.value="grey") + 
            theme_classic() +
            theme(legend.position="bottom") +
            guides(color=guide_legend(title="Molecular function keyword")) +
            theme(legend.direction="vertical") +
            ylim(rev(unique(unlist(as.vector(cells_ligands_long[,2])))[151:200]))+
            theme(legend.text = element_text(size=15))+
            theme(legend.title = element_text(size=15)) +
            theme(legend.position="none")+ 
            theme(axis.text.x = element_text(color="#000000", size=12, angle = 90,vjust=.5, hjust=1),
                  axis.text.y=element_blank(),
                  axis.title.x=element_blank(),
                  axis.title.y=element_blank(),legend.position="none", axis.ticks.y = element_blank(),
                  panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank(),plot.background=element_blank(), axis.line=element_blank())



mf_plot_4

mf_plot_5 <- ggplot(cells_ligands_long_mf, aes(y = ligand, x=x,color=molecular_function_keyword)) +
            geom_point(shape =15, size=3) +
            scale_color_brewer(palette="PuOr", na.value="grey") + 
            theme_classic() +
            theme(legend.position="bottom") +
            guides(color=guide_legend(title="Molecular function keyword")) +
            theme(legend.direction="vertical") +
            ylim(rev(unique(unlist(as.vector(cells_ligands_long[,2])))[201:250]))+
            theme(legend.text = element_text(size=15))+
            theme(legend.title = element_text(size=15)) +
            theme(legend.position="none")+ 
            theme(axis.text.x = element_text(color="#000000", size=12, angle = 90,vjust=.5, hjust=1),
                  axis.text.y=element_blank(),
                  axis.title.x=element_blank(),
                  axis.title.y=element_blank(),legend.position="none", axis.ticks.y = element_blank(),
                  panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank(),plot.background=element_blank(), axis.line=element_blank())



mf_plot_5

mf_plot_6 <- ggplot(cells_ligands_long_mf, aes(y = ligand, x=x,color=molecular_function_keyword)) +
            geom_point(shape =15, size=3) +
            scale_color_brewer(palette="PuOr", na.value="grey") + 
            theme_classic() +
            theme(legend.position="bottom") +
            guides(color=guide_legend(title="Molecular function keyword")) +
            theme(legend.direction="vertical") +
            ylim(rev(unique(unlist(as.vector(cells_ligands_long[,2])))[251:300]))+
            theme(legend.text = element_text(size=15))+
            theme(legend.title = element_text(size=15)) +
            theme(legend.position="none")+ 
            theme(axis.text.x = element_text(color="#000000", size=12, angle = 90,vjust=.5, hjust=1),
                  axis.text.y=element_blank(),
                  axis.title.x=element_blank(),
                  axis.title.y=element_blank(),legend.position="none", axis.ticks.y = element_blank(),
                  panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank(),plot.background=element_blank(), axis.line=element_blank())



mf_plot_6

mf_plot_7 <- ggplot(cells_ligands_long_mf, aes(y = ligand, x=x,color=molecular_function_keyword)) +
            geom_point(shape =15, size=3) +
            scale_color_brewer(palette="PuOr", na.value="grey") + 
            theme_classic() +
            theme(legend.position="bottom") +
            guides(color=guide_legend(title="Molecular function keyword")) +
            theme(legend.direction="vertical") +
            ylim(rev(unique(unlist(as.vector(cells_ligands_long[,2])))[301:348]))+
            theme(legend.text = element_text(size=15))+
            theme(legend.title = element_text(size=15)) +
            theme(legend.position="none")+ 
            theme(axis.text.x = element_text(color="#000000", size=12, angle = 90,vjust=.5, hjust=1),
                  axis.text.y=element_blank(),
                  axis.title.x=element_blank(),
                  axis.title.y=element_blank(),legend.position="none", axis.ticks.y = element_blank(),
                  panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank(),plot.background=element_blank(), axis.line=element_blank())



mf_plot_7
#________________________________________________________________________________
#________________________________________________________________________________

# Draw Only MF Legend without plot

tiff("results/mf_legends.tiff", units="in", width=7, height=5, res=300)
#legend <- get_legend(mf_plot)                    
legend <- cowplot::get_plot_component(mf_plot, 'guide-box-bottom', return_all = TRUE)

grid.newpage()                              

grid.draw(legend) 
dev.off()
#___________________________________________________________________________________________
#___________________________________________________________________________________________
# prepare molecular function annotation 

cells_ligands_long_bp <- cells_ligands_long[,c(2,6)]
cells_ligands_long_bp$x <- "Biological process keyword"
cells_ligands_long_bp$biological_process_keyword <- gsub("\\Cell adhesion - Neurogenesis.*", "Neurogenesis", cells_ligands_long_bp$biological_process_keyword)
cells_ligands_long_bp$biological_process_keyword <- gsub("\\Cell adhesion.*", "Cell adhesion", cells_ligands_long_bp$biological_process_keyword)
cells_ligands_long_bp$biological_process_keyword <- gsub("\\Chemotaxis.*", "Chemotaxis", cells_ligands_long_bp$biological_process_keyword)
cells_ligands_long_bp$biological_process_keyword <- gsub("\\Immunity.*", "Immunity", cells_ligands_long_bp$biological_process_keyword)
cells_ligands_long_bp$biological_process_keyword <- gsub("\\Complement.*", "Immunity", cells_ligands_long_bp$biological_process_keyword)
cells_ligands_long_bp$biological_process_keyword <- gsub("\\Inflammatory.*", "Immunity", cells_ligands_long_bp$biological_process_keyword)
cells_ligands_long_bp$biological_process_keyword <- gsub("\\Antiviral.*", "Immunity", cells_ligands_long_bp$biological_process_keyword)
cells_ligands_long_bp$biological_process_keyword <- str_replace(cells_ligands_long_bp$biological_process_keyword,"B-cell activation", "Immunity")
cells_ligands_long_bp$biological_process_keyword <- gsub("\\Chondrogenesis.*", "Chondrogenesis - Osteogenesis", cells_ligands_long_bp$biological_process_keyword)
cells_ligands_long_bp$biological_process_keyword <- str_replace(cells_ligands_long_bp$biological_process_keyword,"Wnt signaling pathway", "Chondrogenesis - Osteogenesis")
cells_ligands_long_bp$biological_process_keyword <- str_replace(cells_ligands_long_bp$biological_process_keyword,"Notch signaling pathway", "Chondrogenesis - Osteogenesis")
cells_ligands_long_bp$biological_process_keyword <- str_replace(cells_ligands_long_bp$biological_process_keyword,"Endocytosis - Osteogenesis", "Chondrogenesis - Osteogenesis")
cells_ligands_long_bp$biological_process_keyword <- gsub("\\Collagen degradation.*", "Chondrogenesis - Osteogenesis", cells_ligands_long_bp$biological_process_keyword)
cells_ligands_long_bp$biological_process_keyword <- gsub("\\Biomineralization.*", "Chondrogenesis - Osteogenesis", cells_ligands_long_bp$biological_process_keyword)
cells_ligands_long_bp$biological_process_keyword <- gsub("\\Angiogenesis.*", "Angiogenesis", cells_ligands_long_bp$biological_process_keyword)

cells_ligands_long_bp$biological_process_keyword <- gsub("\\Carbohydrate.*", "Metabolism", cells_ligands_long_bp$biological_process_keyword)
cells_ligands_long_bp$biological_process_keyword <- gsub("\\Fatty.*", "Metabolism", cells_ligands_long_bp$biological_process_keyword)
cells_ligands_long_bp$biological_process_keyword <- gsub("\\Lipid.*", "Metabolism", cells_ligands_long_bp$biological_process_keyword)
cells_ligands_long_bp$biological_process_keyword <- gsub("\\Gluconeogenesis.*", "Metabolism", cells_ligands_long_bp$biological_process_keyword)
cells_ligands_long_bp$biological_process_keyword <- gsub("\\Glycolysis.*", "Metabolism", cells_ligands_long_bp$biological_process_keyword)

cells_ligands_long_bp$biological_process_keyword <- gsub("\\Acute.*", "Miscellaneous", cells_ligands_long_bp$biological_process_keyword)
cells_ligands_long_bp$biological_process_keyword <- gsub("\\Apoptosis.*", "Miscellaneous", cells_ligands_long_bp$biological_process_keyword)
cells_ligands_long_bp$biological_process_keyword <- gsub("\\Biological rhythms.*", "Miscellaneous", cells_ligands_long_bp$biological_process_keyword)
cells_ligands_long_bp$biological_process_keyword <- gsub("\\Transcription.*", "Miscellaneous", cells_ligands_long_bp$biological_process_keyword)
cells_ligands_long_bp$biological_process_keyword <- gsub("\\Stress.*", "Miscellaneous", cells_ligands_long_bp$biological_process_keyword)
cells_ligands_long_bp$biological_process_keyword <- gsub("\\Ion.*", "Miscellaneous", cells_ligands_long_bp$biological_process_keyword)
cells_ligands_long_bp$biological_process_keyword <- gsub("\\Growth.*", "Miscellaneous", cells_ligands_long_bp$biological_process_keyword)
cells_ligands_long_bp$biological_process_keyword <- gsub("\\Fertilization.*", "Miscellaneous", cells_ligands_long_bp$biological_process_keyword)
cells_ligands_long_bp$biological_process_keyword <- gsub("\\Exocytosis.*", "Miscellaneous", cells_ligands_long_bp$biological_process_keyword)
cells_ligands_long_bp$biological_process_keyword <- gsub("\\ER-.*", "Miscellaneous", cells_ligands_long_bp$biological_process_keyword)
cells_ligands_long_bp$biological_process_keyword <- gsub("\\Copper.*", "Miscellaneous", cells_ligands_long_bp$biological_process_keyword)
cells_ligands_long_bp$biological_process_keyword <- gsub("\\Blood coagulation.*", "Miscellaneous", cells_ligands_long_bp$biological_process_keyword)
cells_ligands_long_bp$biological_process_keyword <- gsub("\\Calcium channel.*", "Miscellaneous", cells_ligands_long_bp$biological_process_keyword)
cells_ligands_long_bp$biological_process_keyword <- gsub("\\Cell cycle.*", "Miscellaneous", cells_ligands_long_bp$biological_process_keyword)
cells_ligands_long_bp$biological_process_keyword <- gsub("\\Cilium.*", "Miscellaneous", cells_ligands_long_bp$biological_process_keyword)
cells_ligands_long_bp$biological_process_keyword <- str_replace(cells_ligands_long_bp$biological_process_keyword, "Plasminogen activation", "Miscellaneous")
cells_ligands_long_bp$biological_process_keyword <- str_replace(cells_ligands_long_bp$biological_process_keyword, "Neurotransmitter degradation", "Miscellaneous")
cells_ligands_long_bp$biological_process_keyword <- str_replace(cells_ligands_long_bp$biological_process_keyword, "Sensory transduction", "Miscellaneous")
cells_ligands_long_bp$biological_process_keyword <- str_replace(cells_ligands_long_bp$biological_process_keyword, "Neurotransmitter biosynthesis", "Miscellaneous")
cells_ligands_long_bp$biological_process_keyword <- str_replace(cells_ligands_long_bp$biological_process_keyword, "Leukotriene biosynthesis", "Miscellaneous")
cells_ligands_long_bp$biological_process_keyword <- str_replace(cells_ligands_long_bp$biological_process_keyword, "Autophagy - Chemotaxis", "Chemotaxis")


cells_ligands_long_bp$biological_process_keyword <- str_replace(cells_ligands_long_bp$biological_process_keyword, "BMiscellaneous", "Miscellaneous")
cells_ligands_long_bp$biological_process_keyword <- str_replace(cells_ligands_long_bp$biological_process_keyword,"BChondrogenesis - Osteogenesis", "Chondrogenesis - Osteogenesis")


tiff("results/bp.tiff", units="in", width=7, height=50, res=100)

bp_plot <- ggplot(cells_ligands_long_bp, aes(y = ligand, x=x,color=biological_process_keyword)) +
            geom_point(shape =15, size=3) +
            theme_classic() +
            scale_color_brewer(palette="PiYG", na.value="grey", direction=-1) + 
            theme(legend.position="bottom") +
            theme(legend.direction="vertical") +
            ylim(rev(unique(unlist(as.vector(cells_ligands_long[,2])))[1:50]))+
            theme(legend.text = element_text(size=15))+
            theme(legend.title = element_text(size=15)) +
            theme(legend.position="none") + 
            theme(axis.text.x = element_text(color="#000000", size=12, angle = 90,vjust=.5, hjust=1),
            axis.text.y=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),legend.position="bottom", axis.ticks.y = element_blank(),
            panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),plot.background=element_blank(), axis.line=element_blank()) +
            guides(color = guide_legend(override.aes = list(size = 5),title="Biological process keyword"))



bp_plot
dev.off()

bp_plot_1 <- ggplot(cells_ligands_long_bp, aes(y = ligand, x=x,color=biological_process_keyword)) +
            geom_point(shape =15, size=3) +
            theme_classic() +
            scale_color_brewer(palette="PiYG", na.value="grey", direction=-1) + 
            theme(legend.position="bottom") +
            guides(color=guide_legend(title="Biological process keyword")) +
            theme(legend.direction="vertical") +
            ylim(rev(unique(unlist(as.vector(cells_ligands_long[,2])))[1:50]))+
            theme(legend.text = element_text(size=15))+
            theme(legend.title = element_text(size=15)) +
            theme(legend.position="none")+ 
            theme(axis.text.x = element_text(color="#000000", size=12, angle = 90,vjust=.5, hjust=1),
                  axis.text.y=element_blank(),
                  axis.title.x=element_blank(),
                  axis.title.y=element_blank(),legend.position="none", axis.ticks.y = element_blank(),
                  panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank(),plot.background=element_blank(), axis.line=element_blank())
bp_plot_1

bp_plot_1_hor <- ggplot(cells_ligands_long_bp, aes(y = ligand, x=x,color=biological_process_keyword)) +
  geom_point(shape =15, size=3) +
  theme_classic() +
  scale_color_brewer(palette="PiYG", na.value="grey", direction=-1) + 
  theme(legend.position="bottom") +
  guides(color=guide_legend(title="Biological process keyword")) +
  theme(legend.direction="vertical") +
  ylim(rev(unique(unlist(as.vector(cells_ligands_long[,2])))[1:50]))+
  theme(legend.text = element_text(size=15))+
  theme(legend.title = element_text(size=15)) +
  theme(legend.position="none")+ 
  theme(axis.text.y = element_text(color="#000000", size=12, angle = 0,vjust=.5, hjust=1),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),legend.position="none", axis.ticks.x = element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank(), axis.line=element_blank()) + 
  coord_flip()

bp_plot_1_hor

bp_plot_2 <- ggplot(cells_ligands_long_bp, aes(y = ligand, x=x,color=biological_process_keyword)) +
            geom_point(shape =15, size=3) +
            theme_classic() +
            scale_color_brewer(palette="PiYG", na.value="grey", direction=-1) + 
            theme(legend.position="bottom") +
            guides(color=guide_legend(title="Biological process keyword")) +
            theme(legend.direction="vertical") +
            ylim(rev(unique(unlist(as.vector(cells_ligands_long[,2])))[51:100]))+
            theme(legend.text = element_text(size=15))+
            theme(legend.title = element_text(size=15)) +
            theme(legend.position="none")+ 
            theme(axis.text.x = element_text(color="#000000", size=12, angle = 90,vjust=.5, hjust=1),
                  axis.text.y=element_blank(),
                  axis.title.x=element_blank(),
                  axis.title.y=element_blank(),legend.position="none", axis.ticks.y = element_blank(),
                  panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank(),plot.background=element_blank(), axis.line=element_blank())



bp_plot_2

bp_plot_3 <- ggplot(cells_ligands_long_bp, aes(y = ligand, x=x,color=biological_process_keyword)) +
              geom_point(shape =15, size=3) +
              theme_classic() +
              scale_color_brewer(palette="PiYG", na.value="grey", direction=-1) + 
              theme(legend.position="bottom") +
              guides(color=guide_legend(title="Biological process keyword")) +
              theme(legend.direction="vertical") +
              ylim(rev(unique(unlist(as.vector(cells_ligands_long[,2])))[101:150]))+
              theme(legend.text = element_text(size=15))+
              theme(legend.title = element_text(size=15)) +
              theme(legend.position="none")+ 
              theme(axis.text.x = element_text(color="#000000", size=12, angle = 90,vjust=.5, hjust=1),
                    axis.text.y=element_blank(),
                    axis.title.x=element_blank(),
                    axis.title.y=element_blank(),legend.position="none", axis.ticks.y = element_blank(),
                    panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                    panel.grid.minor=element_blank(),plot.background=element_blank(), axis.line=element_blank())



bp_plot_3

bp_plot_4 <- ggplot(cells_ligands_long_bp, aes(y = ligand, x=x,color=biological_process_keyword)) +
            geom_point(shape =15, size=3) +
            theme_classic() +
            scale_color_brewer(palette="PiYG", na.value="grey", direction=-1) + 
            theme(legend.position="bottom") +
            guides(color=guide_legend(title="Biological process keyword")) +
            theme(legend.direction="vertical") +
            ylim(rev(unique(unlist(as.vector(cells_ligands_long[,2])))[151:200]))+
            theme(legend.text = element_text(size=15))+
            theme(legend.title = element_text(size=15)) +
            theme(legend.position="none")+ 
            theme(axis.text.x = element_text(color="#000000", size=12, angle = 90,vjust=.5, hjust=1),
                  axis.text.y=element_blank(),
                  axis.title.x=element_blank(),
                  axis.title.y=element_blank(),legend.position="none", axis.ticks.y = element_blank(),
                  panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank(),plot.background=element_blank(), axis.line=element_blank())



bp_plot_4

bp_plot_5 <- ggplot(cells_ligands_long_bp, aes(y = ligand, x=x,color=biological_process_keyword)) +
  geom_point(shape =15, size=3) +
  theme_classic() +
  scale_color_brewer(palette="PiYG", na.value="grey", direction=-1) + 
  theme(legend.position="bottom") +
  guides(color=guide_legend(title="Biological process keyword")) +
  theme(legend.direction="vertical") +
  ylim(rev(unique(unlist(as.vector(cells_ligands_long[,2])))[201:250]))+
  theme(legend.text = element_text(size=15))+
  theme(legend.title = element_text(size=15)) +
  theme(legend.position="none")+ 
  theme(axis.text.x = element_text(color="#000000", size=12, angle = 90,vjust=.5, hjust=1),
        axis.text.y=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none", axis.ticks.y = element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank(), axis.line=element_blank())

bp_plot_5

bp_plot_6 <- ggplot(cells_ligands_long_bp, aes(y = ligand, x=x,color=biological_process_keyword)) +
            geom_point(shape =15, size=3) +
            theme_classic() +
            scale_color_brewer(palette="PiYG", na.value="grey", direction=-1) + 
            theme(legend.position="bottom") +
            guides(color=guide_legend(title="Biological process keyword")) +
            theme(legend.direction="vertical") +
            ylim(rev(unique(unlist(as.vector(cells_ligands_long[,2])))[251:300]))+
            theme(legend.text = element_text(size=15))+
            theme(legend.title = element_text(size=15)) +
            theme(legend.position="none")+ 
            theme(axis.text.x = element_text(color="#000000", size=12, angle = 90,vjust=.5, hjust=1),
                  axis.text.y=element_blank(),
                  axis.title.x=element_blank(),
                  axis.title.y=element_blank(),legend.position="none", axis.ticks.y = element_blank(),
                  panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank(),plot.background=element_blank(), axis.line=element_blank())



bp_plot_6

bp_plot_7 <- ggplot(cells_ligands_long_bp, aes(y = ligand, x=x,color=biological_process_keyword)) +
  geom_point(shape =15, size=3) +
  theme_classic() +
  scale_color_brewer(palette="PiYG", na.value="grey", direction=-1) + 
  theme(legend.position="bottom") +
  guides(color=guide_legend(title="Biological process keyword")) +
  theme(legend.direction="vertical") +
  ylim(rev(unique(unlist(as.vector(cells_ligands_long[,2])))[301:348]))+
  theme(legend.text = element_text(size=15))+
  theme(legend.title = element_text(size=15)) +
  theme(legend.position="none")+ 
  theme(axis.text.x = element_text(color="#000000", size=12, angle = 90,vjust=.5, hjust=1),
        axis.text.y=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none", axis.ticks.y = element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank(), axis.line=element_blank())



bp_plot_7



#________________________________________________________________________________
#________________________________________________________________________________

# Draw Only MF Legend without plot

tiff("results/bp_legends.tiff", units="in", width=7, height=5, res=300)
#legend <- get_legend(bp_plot)                    
legend <- cowplot::get_plot_component(bp_plot, 'guide-box-bottom', return_all = TRUE)
grid.newpage()                              

grid.draw(legend) 
dev.off()

#___________________________________________________________________________________________
#___________________________________________________________________________________________
tiff("results/cells_ligands_go.annotations_1_100_hor.tiff", units="in", width=10, height=5, res=300)

final_plot_1_100_hor <- plot_grid(bp_plot_1_hor,NULL,mf_plot_1_hor,NULL,cells_ligands_plot_1_hor, 
                        
                              align="v",ncol = 1, rel_heights  = c(0.3,-0.2,0.3,-0.1,1), axis = "l")
final_plot_1_100_hor
dev.off()

tiff("results/cells_ligands_go.annotations_1_100.tiff", units="in", width=10, height=10, res=300)

final_plot_1_100 <- plot_grid(bp_plot_1,NULL,mf_plot_1,NULL,cells_ligands_plot_1, 
                        NULL,
                        bp_plot_2,NULL,mf_plot_2,NULL,cells_ligands_plot_2,
                        align="h",nrow = 1, rel_widths = c(0.3,-0.2,0.3,-0.1,1,
                                                           0,
                                                           0.3,-0.2,0.3,-0.1,1), axis = "b")
final_plot_1_100
dev.off()

tiff("results/cells_ligands_go.annotations_101_200.tiff", units="in", width=10, height=10, res=300)

final_plot_101_200 <- plot_grid(bp_plot_3,NULL,mf_plot_3,NULL,cells_ligands_plot_3, 
                              NULL,
                              bp_plot_4,NULL,mf_plot_4,NULL,cells_ligands_plot_4,
                              align="h",nrow = 1, rel_widths = c(0.3,-0.2,0.3,-0.1,1,
                                                                 0,
                                                                 0.3,-0.2,0.3,-0.1,1), axis = "b")
final_plot_101_200
dev.off()

tiff("results/cells_ligands_go.annotations_201_300.tiff", units="in", width=10, height=10, res=300)

final_plot_201_300 <- plot_grid(bp_plot_5,NULL,mf_plot_5,NULL,cells_ligands_plot_5, 
                                NULL,
                                bp_plot_6,NULL,mf_plot_6,NULL,cells_ligands_plot_6,
                                align="h",nrow = 1, rel_widths = c(0.3,-0.2,0.3,-0.1,1,
                                                                   0,
                                                                   0.3,-0.2,0.3,-0.1,1), axis = "b")
final_plot_201_300
dev.off()

tiff("results/cells_ligands_go.annotations_301_348.tiff", units="in", width=5, height=10, res=300)

final_plot_301_348 <- plot_grid(bp_plot_7,NULL,mf_plot_7,NULL,cells_ligands_plot_7, 
                                align="h",nrow = 1, rel_widths = c(0.3,-0.2,0.3,-0.1,1), axis = "b")
final_plot_301_348
dev.off()

#___________________________________________________________________________________________
#___________________________________________________________________________________________
# frequencty of GO annotations 

ligands_go <- cells_ligands_long[!duplicated(cells_ligands_long[,'ligand']),]

ligands_mo_freq <- ligands_go %>% 
  group_by(molecular_function_keyword) %>% 
  count()


ligands_bp_freq <- ligands_go %>% 
  group_by(biological_process_keyword) %>% 
  count()


#___________________________________________________________________________________________
#___________________________________________________________________________________________
# load receptors data

pc_receptors <- read_excel("data/PC_meta_analysis_SensoryNeurons_neighbors_annotated_220408.xlsx")
pc_receptors <- subset(pc_receptors, pc_sn_nb_r == 'TRUE' & ligand == 'TRUE' & receptor =="TRUE")
pc_receptors <- dplyr::select(pc_receptors, c(hgnc_symbol.y, entrez_id.y, score_RRA_human,  target..family.y))
pc_receptors <- pc_receptors[!duplicated(pc_receptors[,'hgnc_symbol.y']),]
pc_receptors$entrez_id.y <- mapIds(org.Hs.eg.db, pc_receptors$hgnc_symbol.y,  'ENTREZID','SYMBOL')
write.csv(pc_receptors,"results/PC_Receptors.csv", row.names = FALSE)
colnames(pc_receptors)[2] <- "gene"
pc_receptors$gene <- as.character(pc_receptors$gene)
pc_receptors <- merge(pc_receptors,molecular_function_kw, by="gene", all.x = TRUE)
pc_receptors <- merge(pc_receptors,biological_process_kw, by="gene", all.x = TRUE)
pc_receptors$gene <- NULL
colnames(pc_receptors)[1] <- "Plasma Cells"

pc_receptors_freq <- pc_receptors %>% 
  group_by(target..family.y) %>% 
  count()


bmsc_receptors <- read_excel("data/BMSC_meta_analysis_SensoryNeurons_neighbors_annotated_220317.xlsx")
bmsc_receptors <- subset(bmsc_receptors, bmsc_sn_nb_r == 'TRUE' & ligand == 'TRUE' & receptor =="TRUE")
bmsc_receptors <- dplyr::select(bmsc_receptors, c(hgnc_symbol.y, entrez_id.y, score_RRA_human,  target..family.y))
bmsc_receptors <- bmsc_receptors[!duplicated(bmsc_receptors[,'hgnc_symbol.y']),]
bmsc_receptors$entrez_id.y <- mapIds(org.Hs.eg.db, bmsc_receptors$hgnc_symbol.y,  'ENTREZID','SYMBOL')
write.csv(bmsc_receptors,"results/BMSC_Receptors.csv", row.names = FALSE)
colnames(bmsc_receptors)[2] <- "gene"
bmsc_receptors$gene <- as.character(bmsc_receptors$gene)
bmsc_receptors <- merge(bmsc_receptors,molecular_function_kw, by="gene", all.x = TRUE)
bmsc_receptors <- merge(bmsc_receptors,biological_process_kw, by="gene", all.x = TRUE)
bmsc_receptors$gene <- NULL
colnames(bmsc_receptors)[1] <- "Bone Marrow Stromal Cells"


bmsc_receptors_freq <- bmsc_receptors %>% 
  group_by(target..family.y) %>% 
  count()


#___________________________________________________________________________________________
#___________________________________________________________________________________________
# load experimental ligands data
exp_ligands_receptors <- read_excel("data/PC_meta_analysis_SensoryNeurons_neighbors_annotated_220408.xlsx")
exp_ligands_filter <- c("IGF1","MIF", "NRG2", "WNT5A", "THBS1", "FGF7", "SEMA6A")

exp_ligands_receptors <- exp_ligands_receptors[exp_ligands_receptors$hgnc_symbol.x %in% exp_ligands_filter ,]

#___________________________________________________________________________________________
#___________________________________________________________________________________________
# save files
cells_ligands_long_mf <- cells_ligands_long[,c(2,5)]
cells_ligands_long_mf <- cells_ligands_long_mf[!duplicated(cells_ligands_long_mf), ]
cells_ligands_long_bp <- cells_ligands_long[,c(2,6)]
cells_ligands_long_bp <- cells_ligands_long_bp[!duplicated(cells_ligands_long_bp), ]
cells_ligands_long_mf_bp <- merge(cells_ligands_long_mf, cells_ligands_long_bp, by="ligand")

write.csv(cells_ligands_long_mf_bp, "results/Ligand_MF.BP.Keywords.csv", row.names = FALSE)


#___________________________________________________________________________________________
#___________________________________________________________________________________________


# load ligands data
pc_ligands <- read_excel("data/PC_meta_analysis_annotated_220408.xlsx")
pc_ligands <- subset(pc_ligands, pc_sn_nb_r == 'TRUE' & ligand == 'TRUE')
pc_ligands <- dplyr::select(pc_ligands, hgnc_symbol, effectSize, fisherFDRMin)
colnames(pc_ligands) <- c("hgnc_symbol", "effectSize/logFC", "adj.p.value")
pc_ligands$cell_type <- "Plasma cells"
pc_ligands <- pc_ligands[c(4,1,2,3)]


bmsc_ligands <- read_excel("data/BMSC_meta_analysis_annotated_220317.xlsx")
bmsc_ligands <- subset(bmsc_ligands, bmsc_sn_nb_r == 'TRUE' & ligand == 'TRUE')
bmsc_ligands <- dplyr::select(bmsc_ligands, hgnc_symbol, effectSize, fisherFDRMin)
colnames(bmsc_ligands) <- c("hgnc_symbol", "effectSize/logFC", "adj.p.value")
bmsc_ligands$cell_type <- "Bone marrow stromal cells"
bmsc_ligands <- bmsc_ligands[c(4,1,2,3)]

ec_ligands <- read.table("data/EC_DEG_analysis_annotated_220622.tsv", sep = '\t', header = TRUE)
ec_ligands <- subset(ec_ligands, ec_sn_nb_r == 'true' & ligand == 'true')
ec_ligands <- dplyr::select(ec_ligands, hgnc_symbol, logFC, adj.P.Val)
colnames(ec_ligands) <- c("hgnc_symbol", "effectSize/logFC", "adj.p.value")
ec_ligands$cell_type <- "Endothelial cells"
ec_ligands <- ec_ligands[c(4,1,2,3)]

ocy_ligands <- read.table("data/OC_DEG_analysis_annotated_220518.tsv", sep = '\t', header = TRUE)
ocy_ligands <- subset(ocy_ligands, oc_sn_nb_r == 'true' & ligand == 'true')
ocy_ligands <- dplyr::select(ocy_ligands, hgnc_symbol, logFC, adj.P.Val)
colnames(ocy_ligands) <- c("hgnc_symbol", "effectSize/logFC", "adj.p.value")
ocy_ligands$cell_type <- "Osteocytes"
ocy_ligands <- ocy_ligands[c(4,1,2,3)]

opc_ligands <- read.table("data/opc_DEG_analysis_annotated_221215.tsv", sep = '\t', header = TRUE)
opc_ligands <- subset(opc_ligands, sn_nb_r == 'true' & ligand == 'true')
opc_ligands <- dplyr::select(opc_ligands, hgnc_symbol, logFC, adj.P.Val)
colnames(opc_ligands) <- c("hgnc_symbol", "effectSize/logFC", "adj.p.value")
opc_ligands$cell_type <- "Osteogenic precursor cells"
opc_ligands <- opc_ligands[c(4,1,2,3)]

adcy_ligands <- read.table("data/ac_DEG_analysis_annotated_230102.tsv", sep = '\t', header = TRUE)
adcy_ligands <- subset(adcy_ligands, sn_nb_r == 'true' & ligand == 'true')
adcy_ligands$ENTREZID <- mapIds(org.Hs.eg.db, adcy_ligands$hgnc_symbol,  'ENTREZID','SYMBOL')
#adcy_ligands$ENTREZID <- adcy_ligands$ENTREZID$ENTREZID
adcy_ligands <- dplyr::select(adcy_ligands, hgnc_symbol, logFC, adj.P.Val)
colnames(adcy_ligands) <- c("hgnc_symbol", "effectSize/logFC", "adj.p.value")
adcy_ligands$cell_type <- "Adipocytes"
adcy_ligands <- adcy_ligands[c(4,1,2,3)]


nkc_ligands <- read.table("data/NKC_DEG_analysis_annotated_220626.tsv", sep = '\t', header = TRUE)
nkc_ligands <- subset(nkc_ligands, nkc_sn_nb_r == 'true' & ligand == 'true')
nkc_ligands <- dplyr::select(nkc_ligands, hgnc_symbol, logFC, adj.P.Val)
colnames(nkc_ligands) <- c("hgnc_symbol", "effectSize/logFC", "adj.p.value")
nkc_ligands$cell_type <- "Natural killer cells"
nkc_ligands <- nkc_ligands[c(4,1,2,3)]


neut_ligands <- read.table("data/Neutrophil_DEG_analysis_annotated_230102.tsv", sep = '\t', header = TRUE)
neut_ligands <- subset(neut_ligands, sn_nb_r == 'true' & ligand == 'true')
neut_ligands <- dplyr::select(neut_ligands, hgnc_symbol, logFC, adj.P.Val)
colnames(neut_ligands) <- c("hgnc_symbol", "effectSize/logFC", "adj.p.value")
neut_ligands$cell_type <- "Neutrophils"
neut_ligands <- neut_ligands[c(4,1,2,3)]

mp_ligands <- read.table("data/MPH_DEG_analysis_annotated_230104.tsv", sep = '\t', header = TRUE)
mp_ligands <- subset(mp_ligands, sn_nb_r == 'true' & ligand == 'true')
mp_ligands <- dplyr::select(mp_ligands, hgnc_symbol, effectSize, fisherFDRMin)
colnames(mp_ligands) <- c("hgnc_symbol", "effectSize/logFC", "adj.p.value")
mp_ligands$cell_type <- "Macrophages"
mp_ligands <- mp_ligands[c(4,1,2,3)]

treg_ligands <- read.table("data/rtc_GSE109533_DEG_limma_uniprot_annot_230102.tsv", sep = '\t', header = TRUE)
treg_ligands <- subset(treg_ligands, lrpair == 'TRUE' & ligand == 'TRUE')
treg_ligands <- dplyr::select(treg_ligands, hgnc_symbol, logFC, adj.P.Val)
colnames(treg_ligands) <- c("hgnc_symbol", "effectSize/logFC", "adj.p.value")
treg_ligands$cell_type <- "Regulatory T cells"
treg_ligands <- treg_ligands[c(4,1,2,3)]

hsc_ligands <- read.table("data/hst_DEG_analysis_annotated_221215.tsv", sep = '\t', header = TRUE)
hsc_ligands <- subset(hsc_ligands, sn_nb_r == 'true' & ligand == 'true')
hsc_ligands <- dplyr::select(hsc_ligands, hgnc_symbol, logFC, adj.P.Val)
colnames(hsc_ligands) <- c("hgnc_symbol", "effectSize/logFC", "adj.p.value")
hsc_ligands$cell_type <- "Hematopoietic stem cells"
hsc_ligands <- hsc_ligands[c(4,1,2,3)]
hsc_ligands <- hsc_ligands[complete.cases(hsc_ligands), ]

wbm_ligands <- read.table("data/wbm_DEG_analysis_annotated_230109.tsv", sep = '\t', header = TRUE)
wbm_ligands <- subset(wbm_ligands, sn_nb_r == 'true' & ligand == 'true')
wbm_ligands <- dplyr::select(wbm_ligands, hgnc_symbol, logFC, adj.P.Val)
colnames(wbm_ligands) <- c("hgnc_symbol", "effectSize/logFC", "adj.p.value")
wbm_ligands$cell_type <- "Whole bone marrow"
wbm_ligands <- wbm_ligands[c(4,1,2,3)]


cells_ligands <- rbind(pc_ligands, bmsc_ligands, ec_ligands, ocy_ligands, opc_ligands, adcy_ligands,
                       nkc_ligands, neut_ligands, mp_ligands, treg_ligands, hsc_ligands, wbm_ligands)


write.csv(cells_ligands,"results/MM_Cells_Ligands_v3.csv", row.names = FALSE)


cells_ligands_mf_bp <- merge(cells_ligands, cells_ligands_long_mf_bp, by.x="hgnc_symbol", by.y="ligand")[]


write.csv(cells_ligands_mf_bp,"results/MM_Cells_Ligands_MF.BP.Keywords.csv", row.names = FALSE)











