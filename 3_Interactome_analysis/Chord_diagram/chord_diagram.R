# set working directory
setwd("~/Work/PhD/Projects/MoD TrD LRI MM-Pain/Computational Analysis/Ligand_Receptor/Chord_Diagram")

# load packages
library(readxl)
library(tidyverse)
library(circlize)
library(UniProtKeywords)
library(ComplexHeatmap)
library(org.Hs.eg.db)



## annotation of ligands and receptors with keywords from uniprot database
biological_process_kw = load_keyword_genesets(9606, as_table = TRUE,  category = "Biological process")
molecular_function_kw = load_keyword_genesets(9606, as_table = TRUE,  category = "Molecular function")

## load and process data

bmsc <- read_excel("data/BMSC_Chord_1.xlsx")
bmsc <- bmsc[bmsc$num_sources_extracell.y > 1,]
bmsc$entrez_id.y <- select(org.Hs.eg.db, bmsc$hgnc_symbol.y,  'ENTREZID','SYMBOL')
bmsc$entrez_id.y <- bmsc$entrez_id.y$ENTREZID
bmsc$hgnc_symbol.y <- gsub("\\ITG.*", "ITGs", bmsc$hgnc_symbol.y)

bmsc_ligand_receptor <- bmsc[,c(4,5,11,14,20,21,22,27,37,41)]
bmsc_ligand_receptor$fisherFDRMin <- -log10(bmsc_ligand_receptor$fisherFDRMin)
bmsc_ligand_receptor$score_RRA_human <- -log10(bmsc_ligand_receptor$score_RRA_human)

bmsc_ligand_receptor_unique <- subset(bmsc_ligand_receptor, !hgnc_symbol.y %in% hgnc_symbol.x)
bmsc_ligand_receptor_unique$target..family.x[bmsc_ligand_receptor_unique$hgnc_symbol.x=='FGF7'|bmsc_ligand_receptor_unique$hgnc_symbol.x=='VEGFB'|bmsc_ligand_receptor_unique$hgnc_symbol.x=='PGF'] <- 'Growth factor'
bmsc_ligand_receptor_unique$target..family.x[bmsc_ligand_receptor_unique$hgnc_symbol.x=='CXCL12'|bmsc_ligand_receptor_unique$hgnc_symbol.x=='MIF'|bmsc_ligand_receptor_unique$hgnc_symbol.x=='CCL28'] <- 'Cytokine'


## set chord plot options

circos.clear()

tiff("results/bmsc_chord.tiff", height = 30, width = 30, units='cm', res = 300)

circos.par(start.degree = 0, clock.wise = TRUE, canvas.xlim=c(-1.05, 1.05),  
           canvas.ylim=c(-1.05, 1.05), track.margin = c(0.01, 0.05), track.height = 0.05)

par(cex = 1)

## create an empty dataframe to store target family colors 
family_colors <- data.frame(matrix(NA, nrow = 7, ncol = 2))
colnames(family_colors)[1] <- "target..family" 
colnames(family_colors)[2] <- "color"
family_colors[1,] <- c(NA,"#cccccc")
family_colors[2,] <- c("Enzyme","#bf8040")
family_colors[3,] <- c("IC","#e0ca3a")
family_colors[4,] <- c("Kinase","#0099ff")
family_colors[5,] <- c("GPCR","#00ffff")
family_colors[6,] <- c("Growth factor","#d966ff")
family_colors[7,] <- c("Cytokine","#ff0080")

## add ligand/ receptor family color

bmsc_ligand_receptor_unique <- left_join(bmsc_ligand_receptor_unique,family_colors, by= c("target..family.y" = "target..family"))
colnames(bmsc_ligand_receptor_unique)[11] <- "family_color_code_r"

bmsc_ligand_receptor_unique <- left_join(bmsc_ligand_receptor_unique,family_colors, by= c("target..family.x" = "target..family"))
colnames(bmsc_ligand_receptor_unique)[12] <- "family_color_code_l"

## color mapping functions for different variables

col_fun <- colorRamp2(c(max(bmsc_ligand_receptor_unique$effectSize), min(bmsc_ligand_receptor_unique$effectSize)), c("red", "white"))

for(i in 1:nrow(bmsc_ligand_receptor_unique)) {      
  bmsc_ligand_receptor_unique$effectsize_color_code_l[[i]] <- col_fun(bmsc_ligand_receptor_unique$effectSize[[i]])
}  

col_fun2 <- colorRamp2(c(max(bmsc_ligand_receptor_unique$fisherFDRMin), min(bmsc_ligand_receptor_unique$effectSize)), c("yellow", "white"))

for(i in 1:nrow(bmsc_ligand_receptor_unique)) {      
  bmsc_ligand_receptor_unique$fisherFDRMin_color_code_l[[i]] <- col_fun2(bmsc_ligand_receptor_unique$fisherFDRMin[[i]])
}  


col_fun3 <- colorRamp2(c(max(bmsc_ligand_receptor_unique$score_RRA_human), min(bmsc_ligand_receptor_unique$score_RRA_human)), c("orange", "white"))

for(i in 1:nrow(bmsc_ligand_receptor_unique)) {      
  bmsc_ligand_receptor_unique$rra_color_code_r[[i]] <- col_fun3(bmsc_ligand_receptor_unique$score_RRA_human[[i]])
}  

col_fun4 <- colorRamp2(c(max(bmsc_ligand_receptor_unique$num_sources_pain.y), min(bmsc_ligand_receptor_unique$num_sources_pain.y)), c("green", "white"))

for(i in 1:nrow(bmsc_ligand_receptor_unique)) {      
  bmsc_ligand_receptor_unique$pain_color_code_r[[i]] <- col_fun4(bmsc_ligand_receptor_unique$num_sources_pain.y[[i]])
}  


ligand_family_color <- bmsc_ligand_receptor_unique %>% distinct(hgnc_symbol.x, .keep_all=TRUE)
ligand_family_color <- unlist(ligand_family_color$family_color_code_l)

receptor_family_color <- bmsc_ligand_receptor_unique %>% distinct(hgnc_symbol.y, .keep_all=TRUE)
receptor_family_color <- unlist(receptor_family_color$family_color_code_r)

family_color <- c(ligand_family_color,receptor_family_color)

## plot chord diagram

chordDiagram(bmsc_ligand_receptor_unique[,c(1,6)], annotationTrack =  "grid", 
             big.gap = 10, annotationTrackHeight = c(0.03, 0.01), grid.col = family_color, 
               preAllocateTracks = 3)


## add tracks 

ligand_effectsize_color <- bmsc_ligand_receptor_unique %>% distinct(hgnc_symbol.x, .keep_all=TRUE)
ligand_effectsize_color <- unlist(ligand_effectsize_color$effectsize_color_code_l)

receptor_expression_color <- bmsc_ligand_receptor_unique %>% distinct(hgnc_symbol.y, .keep_all=TRUE)
receptor_expression_color <- unlist(receptor_expression_color$rra_color_code_r)

effectsize_expression_color <- c(ligand_effectsize_color,receptor_expression_color)

circos.track(track.index = 3, panel.fun = function(x, y) {
}, bg.col = effectsize_expression_color, bg.border = "black", bg.lwd = 0.2)



ligand_fisherfdrmin_color <- bmsc_ligand_receptor_unique %>% distinct(hgnc_symbol.x, .keep_all=TRUE)
ligand_fisherfdrmin_color <- unlist(ligand_fisherfdrmin_color$fisherFDRMin_color_code_l)

receptor_pain_color <- bmsc_ligand_receptor_unique %>% distinct(hgnc_symbol.y, .keep_all=TRUE)
receptor_pain_color <- unlist(receptor_pain_color$pain_color_code_r)

fisher_pain_color <- c(ligand_fisherfdrmin_color,receptor_pain_color)

circos.track(track.index = 2, panel.fun = function(x, y) {
}, bg.col = fisher_pain_color, bg.border = "black", bg.lwd = 0.2)


circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA, track.height = max(strwidth(unlist(dimnames(bmsc_ligand_receptor_unique)))))


text(0, -1.1, "Ligands", cex = 1.5)
text(0, 1.1, "Receptors", cex = 1.5)

## add legend

# continuous
lgd_links_effectsize = Legend(at = c( 0, 2), col_fun = col_fun,
                              title_position = "topleft", title =  expression("Effect size"), grid_height = unit(100, "cm"),
                              title_gp = gpar(fontsize = 14), labels_gp = gpar(fontsize = 12),
                              legend_height = unit(2, "cm"), grid_width = unit(0.5, "cm"))
lgd_links_fisherfdr = Legend(at = c( 1, 4), col_fun = col_fun2, 
                            title_position = "topleft", title = expression("-log"['10']*"Fisher FDR"), grid_height = unit(100, "cm"),
                            title_gp = gpar(fontsize = 14), labels_gp = gpar(fontsize = 12),
                            legend_height = unit(2, "cm"), grid_width = unit(0.5, "cm"))
lgd_links_rrascore = Legend(at = c(0, 10), col_fun = col_fun3, 
                             title_position = "topleft", title = expression("-log"['10']*"RRA score"), grid_height = unit(100, "cm"),
                             title_gp = gpar(fontsize = 14), labels_gp = gpar(fontsize = 12),
                             legend_height = unit(2, "cm"), grid_width = unit(0.5, "cm"))
lgd_links_pain = Legend(at = c(0, 10), col_fun = col_fun4, 
                            title_position = "topleft", title = expression("Pain-related annotation score"), grid_height = unit(100, "cm"),
                            title_gp = gpar(fontsize = 14), labels_gp = gpar(fontsize = 12),
                            legend_height = unit(2, "cm"), grid_width = unit(0.5, "cm"))


lgd_list_ligand = packLegend(lgd_links_effectsize, lgd_links_fisherfdr)
lgd_list_ligand
draw(lgd_list_ligand, x = unit(4, "mm"), y = unit(10, "mm"), just = c("left", "bottom"))

lgd_list_receptor = packLegend(lgd_links_pain, lgd_links_rrascore)
lgd_list_receptor
draw(lgd_list_receptor, x = unit(4, "mm"), y = unit(240, "mm"), just = c("left", "bottom"))

# discrete
lgd_lines = Legend(labels =  c("NA","Enzyme","Ion channel", "Kinase", "GPCR", "Growth factor", "Cytokine"), 
                   legend_gp = gpar(fill = c("#cccccc","#bf8040","#e0ca3a","#0099ff","#00ffff","#d966ff","#ff0080")), title_position = "topleft", 
                   title = expression("Ligand/ receptor family"),
                   title_gp = gpar(fontsize = 14), labels_gp = gpar(fontsize = 12),
                   legend_height = unit(3, "cm"), grid_width = unit(0.5, "cm"))

lgd_list_class = packLegend(lgd_lines)
lgd_list_class

draw(lgd_list_class, x = unit(300, "mm"), y = unit(45, "mm"), just = c("right", "top"))

dev.off() 

circos.clear()

#_____________________________________________________________________________________
#_____________________________________________________________________________________
## load and process data

pc <- read_excel("data/pc_Chord_1.xlsx")
pc <- pc[pc$num_sources_extracell.y > 1,]
pc$entrez_id.y <- select(org.Hs.eg.db, pc$hgnc_symbol.y,  'ENTREZID','SYMBOL')
pc$entrez_id.y <- pc$entrez_id.y$ENTREZID
pc$hgnc_symbol.y <- gsub("\\ITG.*", "ITGs", pc$hgnc_symbol.y)

pc_ligand_receptor <- pc[,c(4,5,15,18,31,32,33,42,52,58)]
pc_ligand_receptor$fisherFDRMin <- -log10(pc_ligand_receptor$fisherFDRMin)
pc_ligand_receptor$score_RRA_human <- -log10(pc_ligand_receptor$score_RRA_human)

pc_ligand_receptor_unique <- subset(pc_ligand_receptor, !hgnc_symbol.y %in% hgnc_symbol.x)
pc_ligand_receptor_unique$target..family.x[pc_ligand_receptor_unique$hgnc_symbol.x=='IGF1'|pc_ligand_receptor_unique$hgnc_symbol.x=='VEGFB'|pc_ligand_receptor_unique$hgnc_symbol.x=='NRG3'|pc_ligand_receptor_unique$hgnc_symbol.x=='NRG2'|pc_ligand_receptor_unique$hgnc_symbol.x=='HGF'] <- 'Growth factor'
pc_ligand_receptor_unique$target..family.x[pc_ligand_receptor_unique$hgnc_symbol.x=='MIF'] <- 'Cytokine'


## set chord plot options

circos.clear()

tiff("results/pc_chord.tiff", height = 30, width = 30, units='cm', res = 300)

circos.par(start.degree = 0, clock.wise = TRUE, canvas.xlim=c(-1.05, 1.05),  
           canvas.ylim=c(-1.05, 1.05), track.margin = c(0.01, 0.05), track.height = 0.05)

par(cex = 1)

## create an empty dataframe to store target family colors 
family_colors <- data.frame(matrix(NA, nrow = 7, ncol = 2))
colnames(family_colors)[1] <- "target..family" 
colnames(family_colors)[2] <- "color"
family_colors[1,] <- c(NA,"#cccccc")
family_colors[2,] <- c("Enzyme","#bf8040")
family_colors[3,] <- c("IC","#e0ca3a")
family_colors[4,] <- c("Kinase","#0099ff")
family_colors[5,] <- c("GPCR","#00ffff")
family_colors[6,] <- c("Growth factor","#d966ff")
family_colors[7,] <- c("Cytokine","#ff0080")

## add ligand/ receptor family color

pc_ligand_receptor_unique <- left_join(pc_ligand_receptor_unique,family_colors, by= c("target..family.y" = "target..family"))
colnames(pc_ligand_receptor_unique)[11] <- "family_color_code_r"

pc_ligand_receptor_unique <- left_join(pc_ligand_receptor_unique,family_colors, by= c("target..family.x" = "target..family"))
colnames(pc_ligand_receptor_unique)[12] <- "family_color_code_l"

## color mapping functions for different variables

col_fun <- colorRamp2(c(max(pc_ligand_receptor_unique$effectSize), min(pc_ligand_receptor_unique$effectSize)), c("red", "white"))

for(i in 1:nrow(pc_ligand_receptor_unique)) {      
  pc_ligand_receptor_unique$effectsize_color_code_l[[i]] <- col_fun(pc_ligand_receptor_unique$effectSize[[i]])
}  

col_fun2 <- colorRamp2(c(max(pc_ligand_receptor_unique$fisherFDRMin), min(pc_ligand_receptor_unique$effectSize)), c("yellow", "white"))

for(i in 1:nrow(pc_ligand_receptor_unique)) {      
  pc_ligand_receptor_unique$fisherFDRMin_color_code_l[[i]] <- col_fun2(pc_ligand_receptor_unique$fisherFDRMin[[i]])
}  


col_fun3 <- colorRamp2(c(max(pc_ligand_receptor_unique$score_RRA_human), min(pc_ligand_receptor_unique$score_RRA_human)), c("orange", "white"))

for(i in 1:nrow(pc_ligand_receptor_unique)) {      
  pc_ligand_receptor_unique$rra_color_code_r[[i]] <- col_fun3(pc_ligand_receptor_unique$score_RRA_human[[i]])
}  

col_fun4 <- colorRamp2(c(max(pc_ligand_receptor_unique$num_sources_pain.y), min(pc_ligand_receptor_unique$num_sources_pain.y)), c("green", "white"))

for(i in 1:nrow(pc_ligand_receptor_unique)) {      
  pc_ligand_receptor_unique$pain_color_code_r[[i]] <- col_fun4(pc_ligand_receptor_unique$num_sources_pain.y[[i]])
}  


ligand_family_color <- pc_ligand_receptor_unique %>% distinct(hgnc_symbol.x, .keep_all=TRUE)
ligand_family_color <- unlist(ligand_family_color$family_color_code_l)

receptor_family_color <- pc_ligand_receptor_unique %>% distinct(hgnc_symbol.y, .keep_all=TRUE)
receptor_family_color <- unlist(receptor_family_color$family_color_code_r)

family_color <- c(ligand_family_color,receptor_family_color)

## plot chord diagram

chordDiagram(pc_ligand_receptor_unique[,c(1,6)], annotationTrack =  "grid", 
             big.gap = 10, annotationTrackHeight = c(0.03, 0.01), grid.col = family_color, 
             preAllocateTracks = 3)


## add tracks 

ligand_effectsize_color <- pc_ligand_receptor_unique %>% distinct(hgnc_symbol.x, .keep_all=TRUE)
ligand_effectsize_color <- unlist(ligand_effectsize_color$effectsize_color_code_l)

receptor_expression_color <- pc_ligand_receptor_unique %>% distinct(hgnc_symbol.y, .keep_all=TRUE)
receptor_expression_color <- unlist(receptor_expression_color$rra_color_code_r)

effectsize_expression_color <- c(ligand_effectsize_color,receptor_expression_color)

circos.track(track.index = 3, panel.fun = function(x, y) {
}, bg.col = effectsize_expression_color, bg.border = "black", bg.lwd = 0.2)



ligand_fisherfdrmin_color <- pc_ligand_receptor_unique %>% distinct(hgnc_symbol.x, .keep_all=TRUE)
ligand_fisherfdrmin_color <- unlist(ligand_fisherfdrmin_color$fisherFDRMin_color_code_l)

receptor_pain_color <- pc_ligand_receptor_unique %>% distinct(hgnc_symbol.y, .keep_all=TRUE)
receptor_pain_color <- unlist(receptor_pain_color$pain_color_code_r)

fisher_pain_color <- c(ligand_fisherfdrmin_color,receptor_pain_color)

circos.track(track.index = 2, panel.fun = function(x, y) {
}, bg.col = fisher_pain_color, bg.border = "black", bg.lwd = 0.2)


circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA, track.height = max(strwidth(unlist(dimnames(pc_ligand_receptor_unique)))))


text(0, -1.1, "Ligands", cex = 1.5)
text(0, 1.1, "Receptors", cex = 1.5)

## add legend

# continuous
lgd_links_effectsize = Legend(at = c( 0, 2), col_fun = col_fun,
                              title_position = "topleft", title =  expression("Effect size"), grid_height = unit(100, "cm"),
                              title_gp = gpar(fontsize = 14), labels_gp = gpar(fontsize = 12),
                              legend_height = unit(2, "cm"), grid_width = unit(0.5, "cm"))
lgd_links_fisherfdr = Legend(at = c( 5, 25), col_fun = col_fun2, 
                             title_position = "topleft", title = expression("-log"['10']*"fisher FDR"), grid_height = unit(100, "cm"),
                             title_gp = gpar(fontsize = 14), labels_gp = gpar(fontsize = 12),
                             legend_height = unit(2, "cm"), grid_width = unit(0.5, "cm"))
lgd_links_rrascore = Legend(at = c(0, 11), col_fun = col_fun3, 
                            title_position = "topleft", title = expression("-log"['10']*"RRA score"), grid_height = unit(100, "cm"),
                            title_gp = gpar(fontsize = 14), labels_gp = gpar(fontsize = 12),
                            legend_height = unit(2, "cm"), grid_width = unit(0.5, "cm"))
lgd_links_pain = Legend(at = c(0, 13), col_fun = col_fun4, 
                        title_position = "topleft", title = expression("Pain-related annotation score"), grid_height = unit(100, "cm"),
                        title_gp = gpar(fontsize = 14), labels_gp = gpar(fontsize = 12),
                        legend_height = unit(2, "cm"), grid_width = unit(0.5, "cm"))


lgd_list_ligand = packLegend(lgd_links_effectsize, lgd_links_fisherfdr)
lgd_list_ligand
draw(lgd_list_ligand, x = unit(4, "mm"), y = unit(10, "mm"), just = c("left", "bottom"))

lgd_list_receptor = packLegend(lgd_links_pain, lgd_links_rrascore)
lgd_list_receptor
draw(lgd_list_receptor, x = unit(4, "mm"), y = unit(240, "mm"), just = c("left", "bottom"))

# discrete
lgd_lines = Legend(labels =  c("NA","Enzyme","Ion channel", "Kinase", "GPCR", "Growth factor", "Cytokine"), 
                   legend_gp = gpar(fill = c("#cccccc","#bf8040","#e0ca3a","#0099ff","#00ffff","#d966ff","#ff0080")), title_position = "topleft", 
                   title = expression("Ligand/ receptor family"),
                   title_gp = gpar(fontsize = 14), labels_gp = gpar(fontsize = 12),
                   legend_height = unit(3, "cm"), grid_width = unit(0.5, "cm"))

lgd_list_class = packLegend(lgd_lines)
lgd_list_class

draw(lgd_list_class, x = unit(300, "mm"), y = unit(45, "mm"), just = c("right", "top"))

dev.off() 

circos.clear()

#_____________________________________________________________________________________
#_____________________________________________________________________________________

## get GO annotations 

library(GO.db)
library(AnnotationHub)
ah <- AnnotationHub()
orgs <- subset(ah, ah$rdataclass == "OrgDb")
orgdb <- query(orgs, "Homo sapiens")[[1]]

select(orgdb, keys="COL4A5", columns="ENSEMBL", keytype="SYMBOL")
go <- select(orgdb, "COL4A5", "GO", "SYMBOL")
go <- go[go$ONTOLOGY == "MF",]
go2 <- select(GO.db, go$GO, c("TERM","DEFINITION"), "GOID")

ligands <- unique(bmsc_ligand_receptor_unique$hgnc_symbol.x)

mylist <- list()

for (i in (ligands)) {
  select(orgdb, keys=i, columns="ENSEMBL", keytype="SYMBOL")
  go <- select(orgdb, i, "GO", "SYMBOL")
  go <- go[go$ONTOLOGY == "BP",]
  go2 <- select(GO.db, go$GO, c("TERM","DEFINITION"), "GOID")
  go2 <- go2[,2]
  mylist[[i]] <- go2
}
df <- data.frame(do.call("cbind",mylist))


dat.overlap = Reduce(intersect, mylist)


#retrieved <- AnnotationDbi::select(org.Hs.eg.db, keytype="SYMBOL", keys="THBS1", columns="GOALL")














