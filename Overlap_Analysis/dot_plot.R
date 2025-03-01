# set working directory
setwd("~/Work/PhD/Projects/Comput MM Pain/Computational Analysis/Dot_Plot")

# load packages
library(data.table)
library(stats)
library(tidyverse)
library(reshape2)
library(readxl)
library(GeneOverlap)
library(plotly)
library(ggdendro)
library(gdata)
library(cowplot)

#_________________________________________________________________________________________________
#_________________________________________________________________________________________________


mm_cells_order <- c("Plasma Cells","Bone Marrow Stromal Cells", "Endothelial Cellss", "Osteocytess","Osteogenic Precursor Cells",
                    "Adipocytess", "Natural Killer Cells", "Neutrophils", "Macrophagess",
                    "Regulatory T Cellss", "Hematopoietic Stem Cells","Whole Bone Marrow")



#_________________________________________________________________________________________________
#_________________________________________________________________________________________________


# assess sensory neuron human and mouse transcriptome and proteome 

## load data
neuron_human_transcript <- read.csv("data/sensory_neuron/HS_Sensory_Neuron_MA.csv")
neuron_human_transcript <- neuron_human_transcript[-2]
neuron_human_transcript <- as.vector(neuron_human_transcript)

neuron_mouse_transcript <- read.csv("data/sensory_neuron/HS_Mapped_Sensory_Neuron_MA.csv")
neuron_mouse_transcript <- neuron_mouse_transcript[-2]
neuron_mouse_transcript <- as.vector(neuron_mouse_transcript)

### test overlap

neuron_human_mouse_transcript <- newGeneOverlap(neuron_human_transcript$Name,
                                                neuron_mouse_transcript$HGNC.symbol,
                                                genome.size=length(neuron_human_transcript$Name)+length(neuron_mouse_transcript$HGNC.symbol))

neuron_human_mouse_transcript <- testGeneOverlap(neuron_human_mouse_transcript)

neuron_human_mouse_transcript

print(neuron_human_mouse_transcript)

# assess sensory neuron human and mouse protome 

## load data
neuron_human_prot <- read_excel("data/sensory_neuron/41598_2018_31189_MOESM2_ESM.xlsx", skip = 1, col_names = TRUE)
neuron_human_prot <- neuron_human_prot[ -c(1, 3) ]
neuron_human_prot <- separate_rows(neuron_human_prot, "Gene names", sep = ";")
neuron_human_prot <- as.vector(neuron_human_prot)

neuron_mouse_prot <- read_excel("data/sensory_neuron/pain_2021_01_19_zahn_pain-d-20-01021_sdc2.xlsx", skip = 0, col_names = TRUE)
neuron_mouse_prot <- neuron_mouse_prot[2]
neuron_mouse_prot <- as.vector(neuron_mouse_prot)

### test overlap

neuron_human_mouse_prot <- newGeneOverlap(neuron_human_prot$Name,
                                          neuron_mouse_prot$HGNC.symbol,
                                          genome.size=length(neuron_human_prot$Name)+length(neuron_mouse_prot$HGNC.symbol))

neuron_human_mouse_prot <- testGeneOverlap(neuron_human_mouse_prot)

neuron_human_mouse_prot

print(neuron_human_mouse_prot)


# assess sensory neuron overlapped human+mouse transcriptome and human proteome overlap

## load  data
neuron_human_prot <- read_excel("data/sensory_neuron/41598_2018_31189_MOESM2_ESM.xlsx", skip = 1, col_names = TRUE)
neuron_human_prot <- neuron_human_prot[ -c(1, 3) ]
neuron_human_prot <- separate_rows(neuron_human_prot, "Gene names", sep = ";")
neuron_human_prot <- as.vector(neuron_human_prot)

neuron_human.mouse_transcript <- read.csv("data/sensory_neuron/HM_Sensory_Neuron_overlapped.csv")
neuron_human.mouse_transcript <- neuron_human.mouse_transcript[ -c(2, 3) ]
neuron_human.mouse_transcript <- as.vector(neuron_human.mouse_transcript)

### test overlap

neuron_transcript_human.prot <- newGeneOverlap(neuron_human_prot$`Gene names`,
                                               neuron_human.mouse_transcript$Name)
                         
neuron_transcript_human.prot <- testGeneOverlap(neuron_transcript_human.prot)

neuron_transcript_human.prot

print(neuron_transcript_human.prot)

#_________________________________________________________________________________________________
#_________________________________________________________________________________________________

# assess mm cells pathways overlap 

## load data
pc_pathways <- read.csv("data/pathways/PC_Pathways.csv")
pc_pathways$pathway <- with(pc_pathways, ifelse(NES > 0 , paste0(pathway, "_UP"), paste0(pathway, "_DOWN")))
colnames(pc_pathways)[1] <- "Plasma Cells"
pc_pathways <- pc_pathways[1]
pc_pathways <- as.vector(pc_pathways)

bmsc_pathways <- read.csv("data/pathways/BMSC_Pathways.csv")
bmsc_pathways$pathway <- with(bmsc_pathways, ifelse(NES > 0 , paste0(pathway, "_UP"), paste0(pathway, "_DOWN")))
colnames(bmsc_pathways)[1] <- "Bone Marrow Stromal Cells"
bmsc_pathways <- bmsc_pathways[1]
bmsc_pathways <- as.vector(bmsc_pathways)

ec_pathways <- read.csv("data/pathways/EC_Pathways.csv")
ec_pathways$pathway <- with(ec_pathways, ifelse(NES > 0 , paste0(pathway, "_UP"), paste0(pathway, "_DOWN")))
colnames(ec_pathways)[1] <- "Endothelial Cellss"
ec_pathways <- ec_pathways[1]
ec_pathways <- as.vector(ec_pathways)

ocy_pathways <- read.csv("data/pathways/OCY_Pathways.csv")
ocy_pathways$pathway <- with(ocy_pathways, ifelse(NES > 0 , paste0(pathway, "_UP"), paste0(pathway, "_DOWN")))
colnames(ocy_pathways)[1] <- "Osteocytess"
ocy_pathways <- ocy_pathways[1]
ocy_pathways <- as.vector(ocy_pathways)

opc_pathways <- read.csv("data/pathways/opc_Pathways.csv")
opc_pathways$pathway <- with(opc_pathways, ifelse(NES > 0 , paste0(pathway, "_UP"), paste0(pathway, "_DOWN")))
colnames(opc_pathways)[1] <- "Osteogenic Precursor Cells"
opc_pathways <- opc_pathways[1]
opc_pathways <- as.vector(opc_pathways)

adcy_pathways <- read.csv("data/pathways/ADCY_Pathways.csv")
adcy_pathways$pathway <- with(adcy_pathways, ifelse(NES > 0 , paste0(pathway, "_UP"), paste0(pathway, "_DOWN")))
colnames(adcy_pathways)[1] <- "Adipocytess"
adcy_pathways <- adcy_pathways[1]
adcy_pathways <- as.vector(adcy_pathways)

nkc_pathways <- read.csv("data/pathways/NKC_Pathways.csv")
nkc_pathways$pathway <- with(nkc_pathways, ifelse(NES > 0 , paste0(pathway, "_UP"), paste0(pathway, "_DOWN")))
colnames(nkc_pathways)[1] <- "Natural Killer Cells"
nkc_pathways <- nkc_pathways[1]
nkc_pathways <- as.vector(nkc_pathways)

neut_pathways <- read.csv("data/pathways/NEUT_Pathways.csv")
neut_pathways$pathway <- with(neut_pathways, ifelse(NES > 0 , paste0(pathway, "_UP"), paste0(pathway, "_DOWN")))
colnames(neut_pathways)[1] <- "Neutrophils"
neut_pathways <- neut_pathways[1]
neut_pathways <- as.vector(neut_pathways)

mp_pathways <- read.csv("data/pathways/MP_Pathways.csv")
mp_pathways$pathway <- with(mp_pathways, ifelse(NES > 0 , paste0(pathway, "_UP"), paste0(pathway, "_DOWN")))
colnames(mp_pathways)[1] <- "Macrophagess"
mp_pathways <- mp_pathways[1]
mp_pathways <- as.vector(mp_pathways)

treg_pathways <- read.csv("data/pathways/Treg_Pathways.csv")
treg_pathways$pathway <- with(treg_pathways, ifelse(NES > 0 , paste0(pathway, "_UP"), paste0(pathway, "_DOWN")))
colnames(treg_pathways)[1] <- "Regulatory T Cellss"
treg_pathways <- treg_pathways[1]
treg_pathways <- as.vector(treg_pathways)

hsc_pathways <- read.csv("data/pathways/HSC_Pathways.csv")
hsc_pathways$pathway <- with(hsc_pathways, ifelse(NES > 0 , paste0(pathway, "_UP"), paste0(pathway, "_DOWN")))
colnames(hsc_pathways)[1] <- "Hematopoietic Stem Cells"
hsc_pathways <- hsc_pathways[1]
hsc_pathways <- as.vector(hsc_pathways)

wbm_pathways <- read.csv("data/pathways/WBM_Pathways.csv")
wbm_pathways$pathway <- with(wbm_pathways, ifelse(NES > 0 , paste0(pathway, "_UP"), paste0(pathway, "_DOWN")))
colnames(wbm_pathways)[1] <- "Whole Bone Marrow"
wbm_pathways <- wbm_pathways[1]
wbm_pathways <- as.vector(wbm_pathways)


### test overlap
pathways <- list("Plasma Cells" = pc_pathways$`Plasma Cells`,"Bone Marrow Stromal Cells" =bmsc_pathways$`Bone Marrow Stromal Cells`,"Endothelial Cellss"  =ec_pathways$`Endothelial Cellss`, 
                 "Osteocytess"=ocy_pathways$Osteocytess, "Osteogenic Precursor Cells"=opc_pathways$`Osteogenic Precursor Cells`,"Adipocytess"=adcy_pathways$Adipocytess, "Natural Killer Cells"=nkc_pathways$`Natural Killer Cells`, "Neutrophils"=neut_pathways$Neutrophils, 
                 "Macrophagess"=mp_pathways$Macrophagess,"Regulatory T Cellss"=treg_pathways$`Regulatory T Cellss`,"Hematopoietic Stem Cells"=hsc_pathways$`Hematopoietic Stem Cells`,
                  "Whole Bone Marrow"=wbm_pathways$`Whole Bone Marrow`)




pathways.gom.obj <- newGOM(pathways, pathways)
save(pathways.gom.obj, file = "data/pathways/pathways.gom.obj.RData")
load("data/pathways/pathways.gom.obj.RData")

pathways.pval <- data.frame(getMatrix(pathways.gom.obj, name="pval"))
pathways.pval$Cell_1 <- rownames(pathways.pval)
pathways.pval <- melt(setDT(pathways.pval), id.vars = c("Cell_1"), variable.name = "Cell_2")
pathways.pval$Cell_2 <- gsub("\\.", " ", pathways.pval$Cell_2)
colnames(pathways.pval)[3] <- "p_value"
pathways.pval <- filter(pathways.pval, Cell_1 != Cell_2)
#pathways.pval$Cells<- with(pathways.pval, paste(pmin(Cell_1, Cell_2),"_" ,pmax(Cell_1, Cell_2)))
#pathways.pval <- pathways.pval[!duplicated(pathways.pval$Cells),]
#pathways.pval$Cells <- NULL
pathways.pval$p_adj <- p.adjust(pathways.pval$p_value, method = "fdr")
pathways.pval$log10_p_adj <- -log10(pathways.pval$p_adj)
#pathways.pval$log10_p_adj[!is.finite(pathways.pval$log10_p_adj)] <- 0
pathways.pval$log10_p_adj_max <- with(pathways.pval, ifelse(log10_p_adj < 100, log10_p_adj, 100))


pathways.odds <- data.frame(getMatrix(pathways.gom.obj, name="odds.ratio"))
pathways.odds$Cell_1 <- rownames(pathways.odds)
pathways.odds <- melt(setDT(pathways.odds), id.vars = c("Cell_1"), variable.name = "Cell_2")
pathways.odds$Cell_2 <- gsub("\\.", " ", pathways.odds$Cell_2)
colnames(pathways.odds)[3] <- "odds_ratio"
pathways.odds <- filter(pathways.odds, Cell_1 != Cell_2)
#pathways.odds$Cells<- with(pathways.odds, paste(pmin(Cell_1, Cell_2),"_" ,pmax(Cell_1, Cell_2)))
#pathways.odds <- pathways.odds[!duplicated(pathways.odds$Cells),]
#pathways.odds$Cells <- NULL
pathways.odds$odds_ratio_max <- with(pathways.odds, ifelse(odds_ratio < 150, odds_ratio, 150))

pathways.odds.pval <- merge(pathways.pval, pathways.odds, by=c("Cell_1","Cell_2"))



pathways.odds.pval$Cell_1 <- reorder.factor(pathways.odds.pval$Cell_1, new.order=mm_cells_order)
pathways.odds.pval$Cell_2 <- reorder.factor(pathways.odds.pval$Cell_2, new.order=mm_cells_order)
pathways.odds.pval <- pathways.odds.pval %>%
                                        arrange(Cell_1,Cell_2)


#tiff("results/pathways_dot_plot.tiff", units="in", width=14, height=13, res=100)

#dot_plot(data.to.plot = pathways.odds.pval, col_var = pathways.odds.pval$log10_p_adj_max,
#         size_var = pathways.odds.pval$odds_ratio_max, shape.scale = 16, cols.use= c("#d3d3d3","#FFFF00","#FF0000"),
#         size.breaks.values= c(1, 25, 50, 75, 100), display_max_sizes = FALSE,color.breaks.number = 5,
#         dend_y_var= c("log10_p_adj","odds_ratio"),
#         x.lab.pos="bottom", y.lab.pos="left", color.breaks.values = c(0,5,50,100),
#         size_legend="Odds Ratio", col_legend = expression("-log"['10']*"P"['adj']),
#         x.lab.size.factor=1, y.lab.size.factor=1, 
#         dist_method =  "euclidean")

#dev.off()

#_________________________________________________________________________________________________


tiff("results/pathways_dot_plot.tiff", units="in", width=7, height=6, res=300)


ggplot(pathways.odds.pval, aes(x= Cell_1, y = Cell_2, color = odds_ratio, size =  log10_p_adj)) +
  scale_color_gradientn(colors = c( "white", "#ffe0b3", "#ffad33","#ff9900", "#cc7a00"),breaks = c(1, 50, 100,200), limits=c(1,200)) +
  scale_size_continuous(limits = c(1.5, NA))  +
  geom_point() +
  theme_classic() +
  xlab("") +
  ylab("") +
  #theme(plot.title = element_text(face="bold")) +
  #ggtitle("") +
  theme(axis.title = element_text(size=12)) + 
  theme(axis.text.x = element_text( color="#000000", size=10, angle=90, hjust = 0.95, vjust = 0.4),
        axis.text.y = element_text(color="#000000", size=10)) +
  #theme( axis.line = element_line(colour = "#000000", linewidth = 1, linetype = "solid")) +
  labs(color="Odds Ratio") +
  labs(size= expression("-log"['10']*"P"['adj'])) +
  guides(colour = guide_colourbar(order = 1))
  #theme(legend.title = element_text(colour="#000000", size=10, face="bold"))


dev.off()
#_________________________________________________________________________________________________

## plot dendogram
pathways.odds.pval_clust <- pathways.odds.pval %>% 
  select(-p_value, -p_adj, -log10_p_adj_max, -log10_p_adj,-odds_ratio_max) %>%  # drop unused columns to faciliate widening
  pivot_wider(names_from = Cell_2, values_from = odds_ratio) %>% 
  data.frame() # make df as tibbles -> matrix annoying
row.names(pathways.odds.pval_clust) <- pathways.odds.pval_clust$Cell_1  # put gene in `row`
pathways.odds.pval_clust <- pathways.odds.pval_clust[,-1] #drop gene column as now in rows

pathways.odds.pval_clust[is.na(pathways.odds.pval_clust)] <- 0

#is.na(pathways.odds.pval_clust)<-sapply(pathways.odds.pval_clust, is.infinite)
#pathways.odds.pval_clust[is.na(pathways.odds.pval_clust)]<-max(pathways.odds.pval_clust[pathways.odds.pval_clust!=max(pathways.odds.pval_clust)])

pathways.odds.pval_clust <- hclust(dist(pathways.odds.pval_clust %>% as.matrix())) # hclust with distance matrix
dhc <- as.dendrogram(pathways.odds.pval_clust)
data <- dendro_data(dhc, type = "rectangle")

#ggplotly(ggdendrogram(pathways.odds.pval_clust, rotate = TRUE, size = 10))

tiff("results/pathways_dendo_plot.tiff", units="in", width=8, height=7, res=300)
ggplot(data$segments) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend),size=0.5)+
  geom_text(data = data$labels, aes(x, y, label = label), hjust = -0.03, angle = 0, size = 4, fontface=1) +
  coord_flip() + 
  scale_y_reverse() +
  theme_classic () +
  theme(axis.title = element_text()) +
  ylim(300, -300) +
  xlab("Multiple Myeloma Cells") +
  theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.y=element_blank(), axis.text.x=element_blank(), axis.title.y=element_text(size=12), axis.title.x=element_blank()) 
  #theme(axis.text.x = element_text(face="bold", color="#000000", size=10, angle=180, hjust = 0.95, vjust = 0.4))
dev.off()

#_________________________________________________________________________________________________
#_________________________________________________________________________________________________
# assess mm cells ligands overlap 

## load data
pc_ligands <- read_excel("data/ligands/PC_meta_analysis_annotated_220408.xlsx")
pc_ligands <- subset(pc_ligands, pc_sn_nb_r == 'TRUE' & ligand == 'TRUE')
pc_ligands <- select(pc_ligands, hgnc_symbol, fisherFDRMin, effectSize)
colnames(pc_ligands)[1] <- "Plasma Cells"
pc_ligands <- pc_ligands[1]
write.csv(pc_ligands,"data/ligands/PC_Ligands.csv", row.names = FALSE)
pc_ligands <- as.vector(pc_ligands)


bmsc_ligands <- read_excel("data/ligands/BMSC_meta_analysis_annotated_220317.xlsx")
bmsc_ligands <- subset(bmsc_ligands, bmsc_sn_nb_r == 'TRUE' & ligand == 'TRUE')
bmsc_ligands <- select(bmsc_ligands, hgnc_symbol, fisherFDRMin, effectSize)
colnames(bmsc_ligands)[1] <- "Bone Marrow Stromal Cells"
bmsc_ligands <- bmsc_ligands[1]
write.csv(bmsc_ligands,"data/ligands/BMSC_Ligands.csv", row.names = FALSE)
bmsc_ligands <- as.vector(bmsc_ligands)

ec_ligands <- read.table("data/ligands/EC_DEG_analysis_annotated_220622.tsv", sep = '\t', header = TRUE)
ec_ligands <- subset(ec_ligands, ec_sn_nb_r == 'true' & ligand == 'true')
ec_ligands <- select(ec_ligands, hgnc_symbol, adj.P.Val, logFC)
colnames(ec_ligands)[1] <- "Endothelial Cells"
ec_ligands <- ec_ligands[1]
write.csv(ec_ligands,"data/ligands/EC_Ligands.csv", row.names = FALSE)
ec_ligands <- as.vector(ec_ligands)

ocy_ligands <- read.table("data/ligands/OC_DEG_analysis_annotated_220518.tsv", sep = '\t', header = TRUE)
ocy_ligands <- subset(ocy_ligands, oc_sn_nb_r == 'true' & ligand == 'true')
ocy_ligands <- select(ocy_ligands, hgnc_symbol, adj.P.Val, logFC)
colnames(ocy_ligands)[1] <- "Osteocytes"
ocy_ligands <- ocy_ligands[1]
write.csv(ocy_ligands,"data/ligands/OCY_Ligands.csv", row.names = FALSE)
ocy_ligands <- as.vector(ocy_ligands)

nkc_ligands <- read.table("data/ligands/NKC_DEG_analysis_annotated_220626.tsv", sep = '\t', header = TRUE)
nkc_ligands <- subset(nkc_ligands, nkc_sn_nb_r == 'true' & ligand == 'true')
nkc_ligands <- select(nkc_ligands, hgnc_symbol, adj.P.Val, logFC)
colnames(nkc_ligands)[1] <- "Natural Killer Cell"
nkc_ligands <- nkc_ligands[1]
write.csv(nkc_ligands,"data/ligands/NKC_Ligands.csv", row.names = FALSE)
nkc_ligands <- as.vector(nkc_ligands)

### test overlap
ligands <- list("Plasma Cells" = pc_ligands$`Plasma Cells`,"Bone Marrow Stromal Cells" =bmsc_ligands$`Bone Marrow Stromal Cells`,"Endothelial Cells"  =ec_ligands$`Endothelial Cells`, 
                "Osteocytes"=ocy_ligands$Osteocytes,"Natural Killer Cell"=nkc_ligands$`Natural Killer Cell`)




ligands.gom.obj <- newGOM(ligands, ligands)
save(ligands.gom.obj, file = "data/ligands/ligands.gom.obj.RData")
load("data/ligands/ligands.gom.obj.RData")

ligands.pval <- data.frame(getMatrix(ligands.gom.obj, name="pval"))
ligands.pval$Cell_1 <- rownames(ligands.pval)
ligands.pval <- melt(setDT(ligands.pval), id.vars = c("Cell_1"), variable.name = "Cell_2")
ligands.pval$Cell_2 <- gsub("\\.", " ", ligands.pval$Cell_2)
colnames(ligands.pval)[3] <- "p_value"
ligands.pval <- filter(ligands.pval, Cell_1 != Cell_2)
#ligands.pval$Cells<- with(ligands.pval, paste(pmin(Cell_1, Cell_2),"_" ,pmax(Cell_1, Cell_2)))
#ligands.pval <- ligands.pval[!duplicated(ligands.pval$Cells),]
#ligands.pval$Cells <- NULL
ligands.pval$p_adj <- p.adjust(ligands.pval$p_value, method = "fdr")
ligands.pval$log10_p_adj <- -log10(ligands.pval$p_adj)
#ligands.pval$log10_p_adj[!is.finite(ligands.pval$log10_p_adj)] <- 0
ligands.pval$log10_p_adj_max <- with(ligands.pval, ifelse(log10_p_adj < 100, log10_p_adj, 100))


ligands.odds <- data.frame(getMatrix(ligands.gom.obj, name="odds.ratio"))
ligands.odds$Cell_1 <- rownames(ligands.odds)
ligands.odds <- melt(setDT(ligands.odds), id.vars = c("Cell_1"), variable.name = "Cell_2")
ligands.odds$Cell_2 <- gsub("\\.", " ", ligands.odds$Cell_2)
colnames(ligands.odds)[3] <- "odds_ratio"
ligands.odds <- filter(ligands.odds, Cell_1 != Cell_2)
#ligands.odds$Cells<- with(ligands.odds, paste(pmin(Cell_1, Cell_2),"_" ,pmax(Cell_1, Cell_2)))
#ligands.odds <- ligands.odds[!duplicated(ligands.odds$Cells),]
#ligands.odds$Cells <- NULL
ligands.odds$odds_ratio_max <- with(ligands.odds, ifelse(odds_ratio < 150, odds_ratio, 150))

ligands.odds.pval <- merge(ligands.pval, ligands.odds, by=c("Cell_1","Cell_2"))

ligands.odds.pval$Cell_1 <- reorder.factor(ligands.odds.pval$Cell_1, new.order=mm_cells_order)
ligands.odds.pval$Cell_2 <- reorder.factor(ligands.odds.pval$Cell_2, new.order=mm_cells_order)
ligands.odds.pval <- ligands.odds.pval %>%
  arrange(Cell_1,Cell_2)


#tiff("results/ligands_dot_plot.tiff", units="in", width=12, height=10, res=100)

#dot_plot(data.to.plot = ligands.odds.pval, col_var = ligands.odds.pval$log10_p_adj,
#         size_var = ligands.odds.pval$odds_ratio, shape.scale = 20, cols.use= c("#d3d3d3","#FFFF00","#FF0000"),
#         size.breaks.values= c(40, 100, 140), display_max_sizes = FALSE,color.breaks.number = 5,
#         dend_y_var= c("log10_p_adj","odds_ratio"),
#        x.lab.pos="bottom", y.lab.pos="left", color.breaks.values = c(0,5,10,30),
#         size_legend="Odds Ratio", col_legend = expression("-log"['10']*"P"['adj']),
#         x.lab.size.factor=1, y.lab.size.factor=1, 
#         dist_method =  "euclidean")
#dev.off()


ligands.inter <- getNestedList(ligands.gom.obj, name="intersection")

tiff("results/ligands_dot_plot.tiff", units="in", width=8, height=7, res=100)


ggplot(ligands.odds.pval, aes(x= Cell_1, y = Cell_2, color = log10_p_adj, size = odds_ratio)) +
  scale_color_gradientn(colors = c("#ffe6e6", "#ff6666", "#ff0000", "#b30000"),limits = c(1.5, NA),breaks = c(1, 5, 20, 30),oob = scales::squish) +
  scale_size_continuous(limits = c(1, NA))  +
  geom_point() +
  theme_classic() +
  xlab("Multiple Myeloma Cells") +
  ylab("Multiple Myeloma Cells") +
  #theme(plot.title = element_text(face="bold")) +
  #ggtitle("") +
  theme(axis.title = element_text(size=12)) + 
  theme(axis.text.x = element_text( color="#000000", size=10, angle=90, hjust = 0.95, vjust = 0.4),
        axis.text.y = element_text(color="#000000", size=10)) +
  #theme( axis.line = element_line(colour = "#000000", linewidth = 1, linetype = "solid")) +
  labs(color=expression("-log"['10']*"P"['adj'])) +
  labs(size="Odds Ratio") 

#theme(legend.title = element_text(colour="#000000", size=10, face="bold"))

dev.off()
#_________________________________________________________________________________________________

## plot dendogram
ligands.odds.pval_clust <- ligands.odds.pval %>% 
  select(-p_value, -p_adj, -log10_p_adj_max, -log10_p_adj,-odds_ratio_max) %>%  # drop unused columns to faciliate widening
  pivot_wider(names_from = Cell_2, values_from = odds_ratio) %>% 
  data.frame() # make df as tibbles -> matrix annoying
row.names(ligands.odds.pval_clust) <- ligands.odds.pval_clust$Cell_1  # put gene in `row`
ligands.odds.pval_clust <- ligands.odds.pval_clust[,-1] #drop gene column as now in rows

ligands.odds.pval_clust[is.na(ligands.odds.pval_clust)] <- 0

#is.na(ligands.odds.pval_clust)<-sapply(ligands.odds.pval_clust, is.infinite)
#ligands.odds.pval_clust[is.na(ligands.odds.pval_clust)]<-max(ligands.odds.pval_clust[ligands.odds.pval_clust!=max(ligands.odds.pval_clust)])

ligands.odds.pval_clust <- hclust(dist(ligands.odds.pval_clust %>% as.matrix())) # hclust with distance matrix
dhc <- as.dendrogram(ligands.odds.pval_clust)
data <- dendro_data(dhc, type = "rectangle")

#ggplotly(ggdendrogram(ligands.odds.pval_clust, rotate = TRUE, size = 10))

tiff("results/ligands_dendo_plot.tiff", units="in", width=8, height=7, res=100)
ggplot(data$segments) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend),size=0.5)+
  geom_text(data = data$labels, aes(x, y, label = label), hjust = -0.03, angle = 0, size = 4, fontface=1) +
  coord_flip() + 
  scale_y_reverse() +
  theme_classic () +
  theme(axis.title = element_text()) +
  ylim(300, -300) +
  xlab("Multiple Myeloma Cells") +
  theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.y=element_blank(), axis.text.x=element_blank(), axis.title.y=element_text(size=12), axis.title.x=element_blank()) 
#theme(axis.text.x = element_text(face="bold", color="#000000", size=10, angle=180, hjust = 0.95, vjust = 0.4))
dev.off()

#_________________________________________________________________________________________________
#_________________________________________________________________________________________________
## assess mm cells kinases overlap 

## load data
pc_kinases <- read.csv("data/kinases/PC_Kinases.csv")
pc_kinases <- subset(pc_kinases, Hypergeometric.p.value < 0.05)
pc_kinases <- select(pc_kinases, Protein.Kinase, Hypergeometric.p.value)
colnames(pc_kinases)[1] <- "Plasma Cells"
pc_kinases <- pc_kinases[1]
#write.csv(pc_kinases,"data/kinases/PC_kinases.csv", row.names = FALSE)
pc_kinases <- as.vector(pc_kinases)


bmsc_kinases <- read.csv("data/kinases/BMSC_Kinases.csv")
bmsc_kinases <- subset(bmsc_kinases, Hypergeometric.p.value < 0.05)
bmsc_kinases <- select(bmsc_kinases, Protein.Kinase, Hypergeometric.p.value)
colnames(bmsc_kinases)[1] <- "Bone Marrow Stromal Cells"
bmsc_kinases <- bmsc_kinases[1]
#write.csv(bmsc_kinases,"data/kinases/bmsc_kinases.csv", row.names = FALSE)
bmsc_kinases <- as.vector(bmsc_kinases)


ec_kinases <- read.csv("data/kinases/EC_Kinases.csv")
ec_kinases <- subset(ec_kinases, Hypergeometric.p.value < 0.05)
ec_kinases <- select(ec_kinases, Protein.Kinase, Hypergeometric.p.value)
colnames(ec_kinases)[1] <- "Endothelial Cells"
ec_kinases <- ec_kinases[1]
#write.csv(ec_kinases,"data/kinases/ec_kinases.csv", row.names = FALSE)
ec_kinases <- as.vector(ec_kinases)

ocy_kinases <- read.csv("data/kinases/OCY_Kinases.csv")
ocy_kinases <- subset(ocy_kinases, Hypergeometric.p.value < 0.05)
ocy_kinases <- select(ocy_kinases, Protein.Kinase, Hypergeometric.p.value)
colnames(ocy_kinases)[1] <- "Osteocytes"
ocy_kinases <- ocy_kinases[1]
#write.csv(ocy_kinases,"data/kinases/ocy_kinases.csv", row.names = FALSE)
ocy_kinases <- as.vector(ocy_kinases)

nkc_kinases <- read.csv("data/kinases/NKC_Kinases.csv")
nkc_kinases <- subset(nkc_kinases, Hypergeometric.p.value < 0.05)
nkc_kinases <- select(nkc_kinases, Protein.Kinase, Hypergeometric.p.value)
colnames(nkc_kinases)[1] <- "Natural Killer Cell"
nkc_kinases <- nkc_kinases[1]
#write.csv(nkc_kinases,"data/kinases/nkc_kinases.csv", row.names = FALSE)
nkc_kinases <- as.vector(nkc_kinases)

#_________________________________________________________________________________________________

### test overlap
kinases <- list("Plasma Cells" = pc_kinases$`Plasma Cells`,"Bone Marrow Stromal Cells" =bmsc_kinases$`Bone Marrow Stromal Cells`,"Endothelial Cells"  =ec_kinases$`Endothelial Cells`, 
                "Osteocytes"=ocy_kinases$Osteocytes,"Natural Killer Cell"=nkc_kinases$`Natural Killer Cell`)




kinases.gom.obj <- newGOM(kinases, kinases)
save(kinases.gom.obj, file = "data/kinases/kinases.gom.obj.RData")
load("data/kinases/kinases.gom.obj.RData")

kinases.pval <- data.frame(getMatrix(kinases.gom.obj, name="pval"))
kinases.pval$Cell_1 <- rownames(kinases.pval)
kinases.pval <- melt(setDT(kinases.pval), id.vars = c("Cell_1"), variable.name = "Cell_2")
kinases.pval$Cell_2 <- gsub("\\.", " ", kinases.pval$Cell_2)
colnames(kinases.pval)[3] <- "p_value"
kinases.pval <- filter(kinases.pval, Cell_1 != Cell_2)
#kinases.pval$Cells<- with(kinases.pval, paste(pmin(Cell_1, Cell_2),"_" ,pmax(Cell_1, Cell_2)))
#kinases.pval <- kinases.pval[!duplicated(kinases.pval$Cells),]
#kinases.pval$Cells <- NULL
kinases.pval$p_adj <- p.adjust(kinases.pval$p_value, method = "fdr")
kinases.pval$log10_p_adj <- -log10(kinases.pval$p_adj)
#kinases.pval$log10_p_adj[!is.finite(kinases.pval$log10_p_adj)] <- 0
#kinases.pval$log10_p_adj_max <- with(kinases.pval, ifelse(log10_p_adj < 100, log10_p_adj, 100))


kinases.odds <- data.frame(getMatrix(kinases.gom.obj, name="odds.ratio"))
kinases.odds$Cell_1 <- rownames(kinases.odds)
kinases.odds <- melt(setDT(kinases.odds), id.vars = c("Cell_1"), variable.name = "Cell_2")
kinases.odds$Cell_2 <- gsub("\\.", " ", kinases.odds$Cell_2)
colnames(kinases.odds)[3] <- "odds_ratio"
kinases.odds <- filter(kinases.odds, Cell_1 != Cell_2)
#kinases.odds$Cells<- with(kinases.odds, paste(pmin(Cell_1, Cell_2),"_" ,pmax(Cell_1, Cell_2)))
#kinases.odds <- kinases.odds[!duplicated(kinases.odds$Cells),]
#kinases.odds$Cells <- NULL
#kinases.odds$odds_ratio_max <- with(kinases.odds, ifelse(odds_ratio < 150, odds_ratio, 150))

kinases.odds.pval <- merge(kinases.pval, kinases.odds, by=c("Cell_1","Cell_2"))

kinases.odds.pval$Cell_1 <- reorder.factor(kinases.odds.pval$Cell_1, new.order=mm_cells_order)
kinases.odds.pval$Cell_2 <- reorder.factor(kinases.odds.pval$Cell_2, new.order=mm_cells_order)
kinases.odds.pval <- kinases.odds.pval %>%
  arrange(Cell_1,Cell_2)


#tiff("results/kinases_dot_plot.tiff", units="in", width=12, height=10, res=100)
#dot_plot(data.to.plot = kinases.odds.pval, col_var = kinases.odds.pval$log10_p_adj,
#         size_var = kinases.odds.pval$odds_ratio, shape.scale = 20, cols.use= c("#d3d3d3","#FFFF00","#FF0000"),
#         size.breaks.values= c(40, 100, 140), display_max_sizes = FALSE,color.breaks.number = 5,
#         dend_y_var= c("log10_p_adj","odds_ratio"),
#        x.lab.pos="bottom", y.lab.pos="left", color.breaks.values = c(0,5,10,30),
#         size_legend="Odds Ratio", col_legend = expression("-log"['10']*"P"['adj']),
#         x.lab.size.factor=1, y.lab.size.factor=1, 
#         dist_method =  "euclidean")
#dev.off()


kinases.inter <- getNestedList(kinases.gom.obj, name="intersection")
#_________________________________________________________________________________________________

## plot overlap
tiff("results/kinases_dot_plot.tiff", units="in", width=8, height=7, res=100)

ggplot(kinases.odds.pval, aes(x= Cell_1, y = Cell_2, color = log10_p_adj, size = odds_ratio)) +
  scale_color_gradientn(colors = c("#ffe6e6", "#ff6666", "#ff0000", "#b30000"),limits = c(150, 300),breaks = c(150, 200,300)) +
  #scale_size(limits = size.limits) +
  geom_point() +
  theme_classic() +
  xlab("Multiple Myeloma Cells") +
  ylab("Multiple Myeloma Cells") +
  theme(plot.title = element_text(face="bold")) +
  #ggtitle("") +
  theme(axis.title = element_text(size=12)) + 
  theme(axis.text.x = element_text( color="#000000", size=10, angle=90, hjust = 0.95, vjust = 0.4),
        axis.text.y = element_text( color="#000000", size=10)) +
  theme( axis.line = element_line(colour = "#000000", linewidth = 0.5, linetype = "solid")) +
  labs(color=expression("-log"['10']*"P"['adj'])) +
  labs(size="Odds Ratio") +
  guides(colour = guide_colourbar(order = 1))
#theme(legend.title = element_text(colour="#000000", size=10, face="bold"))


dev.off()
#_________________________________________________________________________________________________

## plot dendogram
kinases.odds.pval_clust <- kinases.odds.pval %>% 
  select(-p_value, -p_adj, -log10_p_adj) %>%  # drop unused columns to faciliate widening
  pivot_wider(names_from = Cell_2, values_from = odds_ratio) %>% 
  data.frame() # make df as tibbles -> matrix annoying
row.names(kinases.odds.pval_clust) <- kinases.odds.pval_clust$Cell_1  # put gene in `row`
kinases.odds.pval_clust <- kinases.odds.pval_clust[,-1] #drop gene column as now in rows

kinases.odds.pval_clust[is.na(kinases.odds.pval_clust)] <- 0

#is.na(kinases.odds.pval_clust)<-sapply(kinases.odds.pval_clust, is.infinite)
#kinases.odds.pval_clust[is.na(kinases.odds.pval_clust)]<-max(kinases.odds.pval_clust[kinases.odds.pval_clust!=max(kinases.odds.pval_clust)])

kinases.odds.pval_clust <- hclust(dist(kinases.odds.pval_clust %>% as.matrix())) # hclust with distance matrix
dhc <- as.dendrogram(kinases.odds.pval_clust)
data <- dendro_data(dhc, type = "rectangle")

#ggplotly(ggdendrogram(kinases.odds.pval_clust, rotate = TRUE, size = 10))

tiff("results/kinases_dendo_plot.tiff", units="in", width=8, height=7, res=100)
ggplot(data$segments) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend),size=0.5)+
  geom_text(data = data$labels, aes(x, y, label = label), hjust = -0.03, angle = 0, size = 4, fontface=1) +
  coord_flip() + 
  scale_y_reverse() +
  theme_classic () +
  theme(axis.title = element_text()) +
  ylim(9000, -9000) +
  xlab("Multiple Myeloma Cells") +
  theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.y=element_blank(), axis.text.x=element_blank(), axis.title.y=element_text(size=12), axis.title.x=element_blank()) 
#theme(axis.text.x = element_text(face="bold", color="#000000", size=10, angle=180, hjust = 0.95, vjust = 0.4))
dev.off()

#_________________________________________________________________________________________________

## check kinases across mm cells 

## load data
pc_kinases <- read.csv("data/kinases/PC_Kinases.csv")
pc_kinases <- subset(pc_kinases, Hypergeometric.p.value < 0.05)
pc_kinases <- select(pc_kinases, Protein.Kinase, Hypergeometric.p.value)
pc_kinases$log10_p_value <- -log10(pc_kinases$Hypergeometric.p.value)
pc_kinases$Hypergeometric.p.value <- NULL
colnames(pc_kinases)[2] <- "Plasma Cells"


bmsc_kinases <- read.csv("data/kinases/BMSC_Kinases.csv")
bmsc_kinases <- subset(bmsc_kinases, Hypergeometric.p.value < 0.05)
bmsc_kinases <- select(bmsc_kinases, Protein.Kinase, Hypergeometric.p.value)
bmsc_kinases$log10_p_value <- -log10(bmsc_kinases$Hypergeometric.p.value)
bmsc_kinases$Hypergeometric.p.value <- NULL
colnames(bmsc_kinases)[2] <- "Bone Marrow Stromal Cells"


ec_kinases <- read.csv("data/kinases/EC_Kinases.csv")
ec_kinases <- subset(ec_kinases, Hypergeometric.p.value < 0.05)
ec_kinases <- select(ec_kinases, Protein.Kinase, Hypergeometric.p.value)
ec_kinases$log10_p_value <- -log10(ec_kinases$Hypergeometric.p.value)
ec_kinases$Hypergeometric.p.value <- NULL
colnames(ec_kinases)[2] <- "Endothelial Cells"

ocy_kinases <- read.csv("data/kinases/OCY_Kinases.csv")
ocy_kinases <- subset(ocy_kinases, Hypergeometric.p.value < 0.05)
ocy_kinases <- select(ocy_kinases, Protein.Kinase, Hypergeometric.p.value)
ocy_kinases$log10_p_value <- -log10(ocy_kinases$Hypergeometric.p.value)
ocy_kinases$Hypergeometric.p.value <- NULL
colnames(ocy_kinases)[2] <- "Osteocytes"

# adcy_kinases <- read.csv("data/kinases/adcy_Kinases.csv")
# adcy_kinases <- subset(adcy_kinases, Hypergeometric.p.value < 0.05)
# adcy_kinases <- select(adcy_kinases, Protein.Kinase, Hypergeometric.p.value)
# adcy_kinases$log10_p_value <- -log10(adcy_kinases$Hypergeometric.p.value)
# adcy_kinases$Hypergeometric.p.value <- NULL
# colnames(adcy_kinases)[2] <- "Adipocytes"

mp_kinases <- read.csv("data/kinases/MP_kinases.csv")
mp_kinases <- subset(mp_kinases, Hypergeometric.p.value < 0.05)
mp_kinases <- select(mp_kinases, Protein.Kinase, Hypergeometric.p.value)
mp_kinases$log10_p_value <- -log10(mp_kinases$Hypergeometric.p.value)
mp_kinases$Hypergeometric.p.value <- NULL
colnames(mp_kinases)[2] <- "Macrophages"

treg_kinases <- read.csv("data/kinases/treg_kinases.csv")
treg_kinases <- subset(treg_kinases, Hypergeometric.p.value < 0.05)
treg_kinases <- select(treg_kinases, Protein.Kinase, Hypergeometric.p.value)
treg_kinases$log10_p_value <- -log10(treg_kinases$Hypergeometric.p.value)
treg_kinases$Hypergeometric.p.value <- NULL
colnames(treg_kinases)[2] <- "Regulatory T Cells"

wbm_kinases <- read.csv("data/kinases/wbm_kinases.csv")
wbm_kinases <- subset(wbm_kinases, Hypergeometric.p.value < 0.05)
wbm_kinases <- select(wbm_kinases, Protein.Kinase, Hypergeometric.p.value)
wbm_kinases$log10_p_value <- -log10(wbm_kinases$Hypergeometric.p.value)
wbm_kinases$Hypergeometric.p.value <- NULL
colnames(wbm_kinases)[2] <- "Whole Bone Marrow"

## merge datasets
kinases <- list(pc_kinases, bmsc_kinases, ec_kinases, ocy_kinases,mp_kinases,treg_kinases,wbm_kinases)
kinases <- Reduce(function(x, y) merge(x, y, all=TRUE), kinases)

kinases$sum.na <- rowSums(is.na(kinases))
hist(kinases$sum.na)
kinases <- kinases[kinases$sum.na < 1,]
kinases$sum.na <- NULL
#kinases[is.na(kinases)] <- 1

kinases$mean <- rowMeans(kinases[,c(2:8)])
#kinases$max <- do.call(pmax, kinases[2:9])
#kinases$min <- do.call(pmin, kinases[2:9])
#kinases$range <- kinases$max - kinases$min

kinases <- kinases[order(kinases$mean),]

kinases$Protein.Kinase <- reorder.factor(kinases$Protein.Kinase, new.order=kinases$Protein.Kinase)
kinases <- kinases %>%
  arrange(Protein.Kinase)
kinases$mean <- NULL

kinases_long <- melt(setDT(kinases), id.vars = c("Protein.Kinase"), variable.name = "Cell")
colnames(kinases_long)[3] <- "log10_p_value"

tiff("results/cells_kinases_dot_plot.tiff", units="in", width=8, height=15, res=100)

# remove repeated protein kinases names 

#kinases_long <-subset(kinases_long, Protein.Kinase!="ERK1" & Protein.Kinase!="ERK2" & Protein.Kinase!="JNK1" & Protein.Kinase!="JNK2")

cells_kinases_dot_plot <- 
ggplot(kinases_long, aes(x= Cell, y = Protein.Kinase, size = log10_p_value)) +
  geom_point(colour="grey") +
  scale_size_continuous(range = c(1,6), limits = c(1,30), breaks = c(1,10,30)) +
  theme_classic() +
  xlab("") +
  ylab("") +
  #theme(plot.title = element_text(size = 15)) +
  ggtitle("") +
  theme(axis.title = element_text(size=14)) + 
  theme(axis.text.x = element_text( color="#000000", size=14, angle=90, hjust = 0.95, vjust = 0.4),
        axis.text.y = element_text( color="#000000", size=14)) +
  theme( axis.line = element_line(colour = "#000000", linewidth = 0.5, linetype = "solid")) +
  labs(size=expression("-log"['10']*" p-value")) +
  theme(legend.title = element_text(colour="#000000", size=14)) +
  theme(legend.position="right") 
cells_kinases_dot_plot

dev.off()

#_________________________________________________________________________________________________
#_________________________________________________________________________________________________

## assess mm cells tf overlap 

## load data
pc_tf <- read.csv("data/tf/PC_tf.csv")
pc_tf <- subset(pc_tf, Hypergeometric.p.value < 0.05)
pc_tf <- select(pc_tf, Transcription.Factor, Hypergeometric.p.value)
colnames(pc_tf)[1] <- "Plasma Cells"
pc_tf <- pc_tf[1]
#write.csv(pc_tf,"data/tf/PC_tf.csv", row.names = FALSE)
pc_tf <- as.vector(pc_tf)


bmsc_tf <- read.csv("data/tf/BMSC_tf.csv")
bmsc_tf <- subset(bmsc_tf, Hypergeometric.p.value < 0.05)
bmsc_tf <- select(bmsc_tf, Transcription.Factor, Hypergeometric.p.value)
colnames(bmsc_tf)[1] <- "Bone Marrow Stromal Cells"
bmsc_tf <- bmsc_tf[1]
#write.csv(bmsc_tf,"data/tf/bmsc_tf.csv", row.names = FALSE)
bmsc_tf <- as.vector(bmsc_tf)


ec_tf <- read.csv("data/tf/EC_tf.csv")
ec_tf <- subset(ec_tf, Hypergeometric.p.value < 0.05)
ec_tf <- select(ec_tf, Transcription.Factor, Hypergeometric.p.value)
colnames(ec_tf)[1] <- "Endothelial Cells"
ec_tf <- ec_tf[1]
#write.csv(ec_tf,"data/tf/ec_tf.csv", row.names = FALSE)
ec_tf <- as.vector(ec_tf)

ocy_tf <- read.csv("data/tf/OCY_tf.csv")
ocy_tf <- subset(ocy_tf, Hypergeometric.p.value < 0.05)
ocy_tf <- select(ocy_tf, Transcription.Factor, Hypergeometric.p.value)
colnames(ocy_tf)[1] <- "Osteocytes"
ocy_tf <- ocy_tf[1]
#write.csv(ocy_tf,"data/tf/ocy_tf.csv", row.names = FALSE)
ocy_tf <- as.vector(ocy_tf)

nkc_tf <- read.csv("data/tf/NKC_tf.csv")
nkc_tf <- subset(nkc_tf, Hypergeometric.p.value < 0.05)
nkc_tf <- select(nkc_tf, Transcription.Factor, Hypergeometric.p.value)
colnames(nkc_tf)[1] <- "Natural Killer Cell"
nkc_tf <- nkc_tf[1]
#write.csv(nkc_tf,"data/tf/nkc_tf.csv", row.names = FALSE)
nkc_tf <- as.vector(nkc_tf)

#_________________________________________________________________________________________________

### test overlap
tf <- list("Plasma Cells" = pc_tf$`Plasma Cells`,"Bone Marrow Stromal Cells" =bmsc_tf$`Bone Marrow Stromal Cells`,"Endothelial Cells"  =ec_tf$`Endothelial Cells`, 
           "Osteocytes"=ocy_tf$Osteocytes,"Natural Killer Cell"=nkc_tf$`Natural Killer Cell`)




tf.gom.obj <- newGOM(tf, tf)
save(tf.gom.obj, file = "data/tf/tf.gom.obj.RData")
load("data/tf/tf.gom.obj.RData")

tf.pval <- data.frame(getMatrix(tf.gom.obj, name="pval"))
tf.pval$Cell_1 <- rownames(tf.pval)
tf.pval <- melt(setDT(tf.pval), id.vars = c("Cell_1"), variable.name = "Cell_2")
tf.pval$Cell_2 <- gsub("\\.", " ", tf.pval$Cell_2)
colnames(tf.pval)[3] <- "p_value"
tf.pval <- filter(tf.pval, Cell_1 != Cell_2)
#tf.pval$Cells<- with(tf.pval, paste(pmin(Cell_1, Cell_2),"_" ,pmax(Cell_1, Cell_2)))
#tf.pval <- tf.pval[!duplicated(tf.pval$Cells),]
#tf.pval$Cells <- NULL
tf.pval$p_adj <- p.adjust(tf.pval$p_value, method = "fdr")
tf.pval$log10_p_adj <- -log10(tf.pval$p_adj)
#tf.pval$log10_p_adj[!is.finite(tf.pval$log10_p_adj)] <- 0
#tf.pval$log10_p_adj_max <- with(tf.pval, ifelse(log10_p_adj < 100, log10_p_adj, 100))


tf.odds <- data.frame(getMatrix(tf.gom.obj, name="odds.ratio"))
tf.odds$Cell_1 <- rownames(tf.odds)
tf.odds <- melt(setDT(tf.odds), id.vars = c("Cell_1"), variable.name = "Cell_2")
tf.odds$Cell_2 <- gsub("\\.", " ", tf.odds$Cell_2)
colnames(tf.odds)[3] <- "odds_ratio"
tf.odds <- filter(tf.odds, Cell_1 != Cell_2)
#tf.odds$Cells<- with(tf.odds, paste(pmin(Cell_1, Cell_2),"_" ,pmax(Cell_1, Cell_2)))
#tf.odds <- tf.odds[!duplicated(tf.odds$Cells),]
#tf.odds$Cells <- NULL
#tf.odds$odds_ratio_max <- with(tf.odds, ifelse(odds_ratio < 150, odds_ratio, 150))

tf.odds.pval <- merge(tf.pval, tf.odds, by=c("Cell_1","Cell_2"))

tf.odds.pval$Cell_1 <- reorder.factor(tf.odds.pval$Cell_1, new.order=mm_cells_order)
tf.odds.pval$Cell_2 <- reorder.factor(tf.odds.pval$Cell_2, new.order=mm_cells_order)
tf.odds.pval <- tf.odds.pval %>%
  arrange(Cell_1,Cell_2)


is.na(tf.odds.pval)<-sapply(tf.odds.pval, is.infinite)
tf.odds.pval[is.na(tf.odds.pval)] <- 11000


#tiff("results/tf_dot_plot.tiff", units="in", width=12, height=10, res=100)
#dot_plot(data.to.plot = tf.odds.pval, col_var = tf.odds.pval$log10_p_adj,
#         size_var = tf.odds.pval$odds_ratio, shape.scale = 20, cols.use= c("#d3d3d3","#FFFF00","#FF0000"),
#         size.breaks.values= c(40, 100, 140), display_max_sizes = FALSE,color.breaks.number = 5,
#         dend_y_var= c("log10_p_adj","odds_ratio"),
#        x.lab.pos="bottom", y.lab.pos="left", color.breaks.values = c(0,5,10,30),
#         size_legend="Odds Ratio", col_legend = expression("-log"['10']*"P"['adj']),
#         x.lab.size.factor=1, y.lab.size.factor=1, 
#         dist_method =  "euclidean")
#dev.off()


tf.inter <- getNestedList(tf.gom.obj, name="intersection")
#_________________________________________________________________________________________________

## plot overlap
tiff("results/tf_dot_plot.tiff", units="in", width=8, height=7, res=100)


ggplot(tf.odds.pval, aes(x= Cell_1, y = Cell_2, color = log10_p_adj, size = odds_ratio)) +
  scale_color_gradientn(colors = c("#ffe6e6", "#ff6666", "#ff0000", "#b30000"),limits = c(1, 30),breaks = c(1, 15,30)) +
  #scale_size(limits = size.limits) +
  geom_point() +
  theme_classic() +
  xlab("Multiple Myeloma Cells") +
  ylab("Multiple Myeloma Cells") +
  theme(plot.title = element_text(face="bold")) +
  #ggtitle("") +
  theme(axis.title = element_text(size=12)) + 
  theme(axis.text.x = element_text( color="#000000", size=10, angle=90, hjust = 0.95, vjust = 0.4),
        axis.text.y = element_text( color="#000000", size=10)) +
  theme( axis.line = element_line(colour = "#000000", linewidth = 0.5, linetype = "solid")) +
  labs(color=expression("-log"['10']*"P"['adj'])) +
  labs(size="Odds Ratio") +
  guides(colour = guide_colourbar(order = 1))
#theme(legend.title = element_text(colour="#000000", size=10, face="bold"))



dev.off()
#_________________________________________________________________________________________________

## plot dendogram
tf.odds.pval_clust <- tf.odds.pval %>% 
  select(-p_value, -p_adj, -log10_p_adj) %>%  # drop unused columns to faciliate widening
  pivot_wider(names_from = Cell_2, values_from = odds_ratio) %>% 
  data.frame() # make df as tibbles -> matrix annoying
row.names(tf.odds.pval_clust) <- tf.odds.pval_clust$Cell_1  # put gene in `row`
tf.odds.pval_clust <- tf.odds.pval_clust[,-1] #drop gene column as now in rows

tf.odds.pval_clust[is.na(tf.odds.pval_clust)] <- 0

#is.na(tf.odds.pval_clust)<-sapply(tf.odds.pval_clust, is.infinite)
#tf.odds.pval_clust[is.na(tf.odds.pval_clust)]<-max(tf.odds.pval_clust[tf.odds.pval_clust!=max(tf.odds.pval_clust)])

tf.odds.pval_clust <- hclust(dist(tf.odds.pval_clust %>% as.matrix())) # hclust with distance matrix
dhc <- as.dendrogram(tf.odds.pval_clust)
data <- dendro_data(dhc, type = "rectangle")

#ggplotly(ggdendrogram(tf.odds.pval_clust, rotate = TRUE, size = 10))

tiff("results/tf_dendo_plot.tiff", units="in", width=8, height=7, res=100)
ggplot(data$segments) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend),size=0.5)+
  geom_text(data = data$labels, aes(x, y, label = label), hjust = -0.03, angle = 0, size = 4, fontface=1) +
  coord_flip() + 
  scale_y_reverse() +
  theme_classic () +
  theme(axis.title = element_text()) +
  ylim(19000, -19000) +
  xlab("Multiple Myeloma Cells") +
  theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.y=element_blank(), axis.text.x=element_blank(), axis.title.y=element_text(size=12), axis.title.x=element_blank()) 
#theme(axis.text.x = element_text(face="bold", color="#000000", size=10, angle=180, hjust = 0.95, vjust = 0.4))
dev.off()

#_________________________________________________________________________________________________

## check tf across mm cells 

## load data
pc_tf <- read.csv("data/tf/PC_tf.csv")
pc_tf <- subset(pc_tf, Hypergeometric.p.value < 0.05)
pc_tf <- select(pc_tf, Transcription.Factor, Hypergeometric.p.value)
pc_tf$log10_p_value <- -log10(pc_tf$Hypergeometric.p.value)
pc_tf$Hypergeometric.p.value <- NULL
colnames(pc_tf)[2] <- "Plasma Cells"


bmsc_tf <- read.csv("data/tf/BMSC_tf.csv")
bmsc_tf <- subset(bmsc_tf, Hypergeometric.p.value < 0.05)
bmsc_tf <- select(bmsc_tf, Transcription.Factor, Hypergeometric.p.value)
bmsc_tf$log10_p_value <- -log10(bmsc_tf$Hypergeometric.p.value)
bmsc_tf$Hypergeometric.p.value <- NULL
colnames(bmsc_tf)[2] <- "Bone Marrow Stromal Cells"


ec_tf <- read.csv("data/tf/EC_tf.csv")
ec_tf <- subset(ec_tf, Hypergeometric.p.value < 0.05)
ec_tf <- select(ec_tf, Transcription.Factor, Hypergeometric.p.value)
ec_tf$log10_p_value <- -log10(ec_tf$Hypergeometric.p.value)
ec_tf$Hypergeometric.p.value <- NULL
colnames(ec_tf)[2] <- "Endothelial Cells"

ocy_tf <- read.csv("data/tf/OCY_tf.csv")
ocy_tf <- subset(ocy_tf, Hypergeometric.p.value < 0.05)
ocy_tf <- select(ocy_tf, Transcription.Factor, Hypergeometric.p.value)
ocy_tf$log10_p_value <- -log10(ocy_tf$Hypergeometric.p.value)
ocy_tf$Hypergeometric.p.value <- NULL
colnames(ocy_tf)[2] <- "Osteocytes"

# adcy_tf <- read.csv("data/tf/adcy_tf.csv")
# adcy_tf <- subset(adcy_tf, Hypergeometric.p.value < 0.05)
# adcy_tf <- select(adcy_tf, Transcription.Factor, Hypergeometric.p.value)
# adcy_tf$log10_p_value <- -log10(adcy_tf$Hypergeometric.p.value)
# adcy_tf$Hypergeometric.p.value <- NULL
# colnames(adcy_tf)[2] <- "Adipocytes"

mp_tf <- read.csv("data/tf/MP_tf.csv")
mp_tf <- subset(mp_tf, Hypergeometric.p.value < 0.05)
mp_tf <- select(mp_tf, Transcription.Factor, Hypergeometric.p.value)
mp_tf$log10_p_value <- -log10(mp_tf$Hypergeometric.p.value)
mp_tf$Hypergeometric.p.value <- NULL
colnames(mp_tf)[2] <- "Macrophages"

treg_tf <- read.csv("data/tf/treg_tf.csv")
treg_tf <- subset(treg_tf, Hypergeometric.p.value < 0.05)
treg_tf <- select(treg_tf, Transcription.Factor, Hypergeometric.p.value)
treg_tf$log10_p_value <- -log10(treg_tf$Hypergeometric.p.value)
treg_tf$Hypergeometric.p.value <- NULL
treg_tf <- treg_tf %>% 
  group_by(Transcription.Factor) %>% 
  slice_max(order_by = log10_p_value) %>% as.data.frame()
colnames(treg_tf)[2] <- "Regulatory T Cells"


wbm_tf <- read.csv("data/tf/wbm_tf.csv")
wbm_tf <- subset(wbm_tf, Hypergeometric.p.value < 0.05)
wbm_tf <- select(wbm_tf, Transcription.Factor, Hypergeometric.p.value)
wbm_tf$log10_p_value <- -log10(wbm_tf$Hypergeometric.p.value)
wbm_tf$Hypergeometric.p.value <- NULL
colnames(wbm_tf)[2] <- "Whole Bone Marrow"

## merge datasets
tf <- list(pc_tf, bmsc_tf, ec_tf, ocy_tf,mp_tf,treg_tf,wbm_tf)
tf <- Reduce(function(x, y) merge(x, y, all=TRUE), tf)

tf$sum.na <- rowSums(is.na(tf))
hist(tf$sum.na)
#tf <- tf[tf$sum.na < 7,]
tf$sum.na <- NULL
tf[is.na(tf)] <- 0

tf$mean <- rowMeans(tf[,c(2:8)])
#tf$max <- do.call(pmax, tf[2:9])
#tf$min <- do.call(pmin, tf[2:9])
#tf$range <- tf$max - tf$min

tf <- tf[order(tf$mean),]

tf$Transcription.Factor <- reorder.factor(tf$Transcription.Factor, new.order=tf$Transcription.Factor)
tf <- tf %>%
  arrange(Transcription.Factor)
tf$mean <- NULL

tf_long <- melt(setDT(tf), id.vars = c("Transcription.Factor"), variable.name = "Cell")
colnames(tf_long)[3] <- "log10_p_value"
tf_long[tf_long==0] <- NA
tf_long<-tf_long[complete.cases(tf_long),]

tiff("results/cells_tf_dot_plot.tiff", units="in", width=8, height=7.25, res=100)

cells_tf_dot_plot <-
ggplot(tf_long, aes(x= Cell, y = Transcription.Factor, size = log10_p_value)) +
  geom_point(colour="grey") +
  scale_size_continuous(range = c(1,6), limits = c(1,30), breaks = c(1,10,30)) +
  theme_classic() +
  xlab("") +
  ylab("") +
  #theme(plot.title = element_text(size = 15)) +
  ggtitle("") +
  theme(axis.title = element_text(size=14)) + 
  theme(axis.text.x = element_text( color="#000000", size=14, angle=90, hjust = 0.95, vjust = 0.4),
        axis.text.y = element_text( color="#000000", size=14)) +
  theme( axis.line = element_line(colour = "#000000", linewidth = 0.5, linetype = "solid")) +
  labs(size=expression("-log"['10']*" p-value")) +
  theme(legend.title = element_text(colour="#000000", size=14)) +
  theme(legend.position="none") 

cells_tf_dot_plot
dev.off()


#________________________________________________________________________________
#________________________________________________________________________________

## plot final plot 

tiff("results/cells_kinases_tf_plot.tiff", units="in", width=10, height=15, res=300)

cells_kinases_tf_plot <- plot_grid(cells_tf_dot_plot, NULL, cells_kinases_dot_plot, 
                        rel_widths=c(1.1,0.2,1.7),nrow=1, align = "h",axis = "b")

cells_kinases_tf_plot
dev.off()

#________________________________________________________________________________
#________________________________________________________________________________

## save plots

saveRDS(cells_kinases_tf_plot, file="results/cells_kinases_tf_plot.RData")

cells_kinases_tf_plot <- readRDS("results/cells_kinases_tf_plot.RData")
cells_kinases_tf_plot







