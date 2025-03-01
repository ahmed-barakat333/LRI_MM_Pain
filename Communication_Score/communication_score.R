# set working directory
setwd("~/Work/PhD/Projects/Comput MM Pain/Computational Analysis/Communication_Score")

# load data
load("data/data_string_PC_BMSC_220613.RData")

# create an empty dataframe to store p values
p_value <- data.frame(matrix(NA, nrow = 3, ncol = 4))
rownames(p_value) <- c("0.4 cutoff","0.7 cutoff","0.9 cutoff")
colnames(p_value) <- c("BMSC(All)_Receptor(All)","BMSC(All)_Receptor(Pain)","PC(All)_Receptor(All)","PC(All)_Receptor(Pain)")

####BMSC(All)_Receptor(All)#### 
# filter string interactions 0.4
string_human_physical_uniprot_0.4 <- string_human_physical_uniprot[string_human_physical_uniprot$score >= 400,]

# check how many sensory neuron are annotated as receptors and have at least 1 extracellular resource
nrow(proteins_sensory_neurons_annot_transc[proteins_sensory_neurons_annot_transc$receptor == TRUE & proteins_sensory_neurons_annot_transc$num_sources_extracell > 0, ])
#[1] 786

# check how many DEG of MM-BMSC genes are annotated as ligand and have at least 1 extracellular resource
nrow(merged_analysis_BMSC_5DS_ID_proteome_annot[merged_analysis_BMSC_5DS_ID_proteome_annot$ligand == TRUE & merged_analysis_BMSC_5DS_ID_proteome_annot$num_sources_extracell > 0, ])
#[1] 104

# create a subset dataframe that contains DEG of MM-BMSC extracellular ligands 
mm_bmsc_ligands_ex <- merged_analysis_BMSC_5DS_ID_proteome_annot[merged_analysis_BMSC_5DS_ID_proteome_annot$ligand == TRUE & merged_analysis_BMSC_5DS_ID_proteome_annot$num_sources_extracell > 0, ]

# get degree_0.4 (no. of interactions) distribution of DEG of MM-BMSC ligands - 0.4
mm_bmsc_ligands_ex$degree_0.4 <- NA

for(i in 1:nrow(mm_bmsc_ligands_ex)) {      
  protein <- mm_bmsc_ligands_ex[i, ]
  protein_degree_0.4 <-subset(string_human_physical_uniprot_0.4, 
                              string_human_physical_uniprot_0.4$protein1 %in% 
                                protein$uniprotswissprot)
  
  mm_bmsc_ligands_ex$degree_0.4[[i]] <- nrow(protein_degree_0.4)
  
}

# plot the distribution
mm_bmsc_degree_0.4 <-  as.numeric(unlist(mm_bmsc_ligands_ex$degree_0.4, use.names = FALSE))
hist(mm_bmsc_degree_0.4)

# sum of MM-BMSC degree_0.4
sum(mm_bmsc_degree_0.4)
#[1] 3013

# check how many all BMSC genes are annotated as ligand and have at least 1 extracellular resource
nrow(merged_analysis_BMSC_5DS_ID_all[merged_analysis_BMSC_5DS_ID_all$ligand == TRUE & merged_analysis_BMSC_5DS_ID_all$num_sources_extracell > 0, ])
#[1] 763

# create a subset dataframe that contains all BMSC extracellular ligands to sample from
all_bmsc_ligands_ex <- merged_analysis_BMSC_5DS_ID_all[merged_analysis_BMSC_5DS_ID_all$ligand == TRUE & merged_analysis_BMSC_5DS_ID_all$num_sources_extracell > 0, ]
 
# get degree_0.4 (no. of interactions) distribution of all BMSC ligands - 0.4
all_bmsc_ligands_ex$degree_0.4 <- NA
for(i in 1:nrow(all_bmsc_ligands_ex)) {      
  protein <- all_bmsc_ligands_ex[i, ]
  protein_degree_0.4 <-subset(string_human_physical_uniprot_0.4, 
                              string_human_physical_uniprot_0.4$protein1 %in% 
                                protein$uniprotswissprot)
  
  all_bmsc_ligands_ex$degree_0.4[[i]] <- nrow(protein_degree_0.4)
  
}

# plot the distribution
all_bmsc_degree_0.4 <- as.numeric(unlist(all_bmsc_ligands_ex$degree_0.4, use.names = FALSE))
hist(all_bmsc_degree_0.4)

# sum of all BMSC degree_0.4
sum(all_bmsc_degree_0.4)
#[1] 17551

# density function for sampling same degree distribution
set.seed(3333)
bmsc.dens.obs_0.4 <- density(mm_bmsc_ligands_ex$degree, n=length(all_bmsc_ligands_ex$degree))

# calculate p value for observed interactions for BMSC at 0.4 cutoff

## observed interactions for DEG from MM-BMSC
mm_bmsc_104_receptor_interactions_0.4 <- subset(string_human_physical_uniprot_0.4, 
                                                string_human_physical_uniprot_0.4$protein1 %in% 
                                                  mm_bmsc_ligands_ex$uniprotswissprot 
                                                & string_human_physical_uniprot_0.4$protein2 %in% 
                                                  proteins_sensory_neurons_annot_transc[proteins_sensory_neurons_annot_transc$receptor == TRUE & proteins_sensory_neurons_annot_transc$num_sources_extracell > 0 ,]$uniprot_id)
nrow(mm_bmsc_104_receptor_interactions_0.4)
#[1] 678

## create FOR loop for bootstrapping 
bmsc_104_1e4_receptor_interactions_0.4 <- list()
set.seed(3333)
for(i in 1:10000){
  sample <- all_bmsc_ligands_ex[sample(nrow(all_bmsc_ligands_ex), 104, replace = TRUE, prob=bmsc.dens.obs_0.4$y), ]
  interactions_sample <-subset(string_human_physical_uniprot_0.4, 
                               string_human_physical_uniprot_0.4$protein1 %in% 
                                 sample$uniprotswissprot 
                               & string_human_physical_uniprot_0.4$protein2 %in% 
                                 proteins_sensory_neurons_annot_transc[proteins_sensory_neurons_annot_transc$receptor == TRUE & proteins_sensory_neurons_annot_transc$num_sources_extracell > 0,]$uniprot_id)
  
  bmsc_104_1e4_receptor_interactions_0.4[[i]] <- nrow(interactions_sample)
}


## plot null (normal) sampling distribution
bmsc_104_1e4_receptor_interactions_0.4 <- unlist(bmsc_104_1e4_receptor_interactions_0.4, use.names = FALSE)
hist(bmsc_104_1e4_receptor_interactions_0.4)

write.csv(bmsc_104_1e4_receptor_interactions_0.4,file="data/bmsc_104_1e4_receptor_interactions_0.4.csv",row.names = FALSE)
bmsc_104_1e4_receptor_interactions_0.4 <- as.vector(t(read.csv("data/bmsc_104_1e4_receptor_interactions_0.4.csv")))

## add observed interaction 
abline(v=678, col="blue")

## calculate p value
pvalue_bmsc_104_1e4_receptor_interactions_0.4 <- mean(bmsc_104_1e4_receptor_interactions_0.4>=678)
#[1] 0.047

## store the p value
p_value[1,1] <- mean(bmsc_104_1e4_receptor_interactions_0.4>=678)

# calculate p value for observed interactions for BMSC at 0.7 cutoff

## filter string interactions
string_human_physical_uniprot_0.7 <- string_human_physical_uniprot[string_human_physical_uniprot$score >= 700,]

# get degree_0.7 (no. of interactions) distribution of DEG of MM-BMSC ligands - 0.7
mm_bmsc_ligands_ex$degree_0.7 <- NA

for(i in 1:nrow(mm_bmsc_ligands_ex)) {      
  protein <- mm_bmsc_ligands_ex[i, ]
  protein_degree_0.7 <-subset(string_human_physical_uniprot_0.7, 
                              string_human_physical_uniprot_0.7$protein1 %in% 
                                protein$uniprotswissprot)
  
  mm_bmsc_ligands_ex$degree_0.7[[i]] <- nrow(protein_degree_0.7)
  
}

# plot the distribution
mm_bmsc_degree_0.7 <-  as.numeric(unlist(mm_bmsc_ligands_ex$degree_0.7, use.names = FALSE))
hist(mm_bmsc_degree_0.7)

# sum of MM-BMSC degree_0.7
sum(mm_bmsc_degree_0.7)
#[1] 1145

# get degree_0.7 (no. of interactions) distribution of all BMSC ligands - 0.7
all_bmsc_ligands_ex$degree_0.7 <- NA
for(i in 1:nrow(all_bmsc_ligands_ex)) {      
  protein <- all_bmsc_ligands_ex[i, ]
  protein_degree_0.7 <-subset(string_human_physical_uniprot_0.7, 
                              string_human_physical_uniprot_0.7$protein1 %in% 
                                protein$uniprotswissprot)
  
  all_bmsc_ligands_ex$degree_0.7[[i]] <- nrow(protein_degree_0.7)
  
}

# plot the distribution
all_bmsc_degree_0.7 <- as.numeric(unlist(all_bmsc_ligands_ex$degree_0.7, use.names = FALSE))
hist(all_bmsc_degree_0.7)

# sum of all BMSC degree_0.7
sum(all_bmsc_degree_0.7)
#[1] 6683

# density function 
set.seed(3333)
bmsc.dens.obs_0.7 <- density(mm_bmsc_ligands_ex$degree_0.7, n=length(all_bmsc_ligands_ex$degree_0.7), adjust = 2)

## create FOR loop for bootstrapping 
bmsc_104_1e4_receptor_interactions_0.7 <- list()
set.seed(3333)
for(i in 1:10000){
  sample <- all_bmsc_ligands_ex[sample(nrow(all_bmsc_ligands_ex), 104, replace = TRUE, prob=bmsc.dens.obs_0.7$y), ]
  interactions_sample <-subset(string_human_physical_uniprot_0.7, 
                            string_human_physical_uniprot_0.7$protein1 %in% 
                           sample$uniprotswissprot 
                         & string_human_physical_uniprot_0.7$protein2 %in% 
                           proteins_sensory_neurons_annot_transc[proteins_sensory_neurons_annot_transc$receptor == TRUE & proteins_sensory_neurons_annot_transc$num_sources_extracell > 0,]$uniprot_id)
   
  bmsc_104_1e4_receptor_interactions_0.7[[i]] <- nrow(interactions_sample)
    }


## plot null (normal) sampling distribution
bmsc_104_1e4_receptor_interactions_0.7 <- unlist(bmsc_104_1e4_receptor_interactions_0.7, use.names = FALSE)
hist(bmsc_104_1e4_receptor_interactions_0.7)

write.csv(bmsc_104_1e4_receptor_interactions_0.7,file="data/bmsc_104_1e4_receptor_interactions_0.7.csv",row.names = FALSE)
bmsc_104_1e4_receptor_interactions_0.7 <- as.vector(t(read.csv("data/bmsc_104_1e4_receptor_interactions_0.7.csv")))

## observed interactions for DEG from MM-BMSC
mm_bmsc_104_receptor_interactions_0.7 <- subset(string_human_physical_uniprot_0.7, 
                                        string_human_physical_uniprot_0.7$protein1 %in% 
                                          mm_bmsc_ligands_ex$uniprotswissprot 
                                        & string_human_physical_uniprot_0.7$protein2 %in% 
                                          proteins_sensory_neurons_annot_transc[proteins_sensory_neurons_annot_transc$receptor == TRUE & proteins_sensory_neurons_annot_transc$num_sources_extracell > 0 ,]$uniprot_id)
nrow(mm_bmsc_104_receptor_interactions_0.7)
#[1] 390

## add observed interaction 
abline(v=390, col="blue")

## calculate p value
pvalue_bmsc_104_1e4_receptor_interactions_0.7 <- mean(bmsc_104_1e4_receptor_interactions_0.7>=390)
#[1] 0.0046

## store the p value
p_value[2,1] <- mean(bmsc_104_1e4_receptor_interactions_0.7>=390)

# calculate p value for observed interactions for BMSC at 0.9 cutoff ######

## filter string interactions
string_human_physical_uniprot_0.9 <- string_human_physical_uniprot[string_human_physical_uniprot$score >= 900,]

# get degree_0.9 (no. of interactions) distribution of DEG of MM-BMSC ligands - 0.9
mm_bmsc_ligands_ex$degree_0.9 <- NA

for(i in 1:nrow(mm_bmsc_ligands_ex)) {      
  protein <- mm_bmsc_ligands_ex[i, ]
  protein_degree_0.9 <-subset(string_human_physical_uniprot_0.9, 
                              string_human_physical_uniprot_0.9$protein1 %in% 
                                protein$uniprotswissprot)
  
  mm_bmsc_ligands_ex$degree_0.9[[i]] <- nrow(protein_degree_0.9)
  
}

# plot the distribution
mm_bmsc_degree_0.9 <-  as.numeric(unlist(mm_bmsc_ligands_ex$degree_0.9, use.names = FALSE))
hist(mm_bmsc_degree_0.9)

# sum of MM-BMSC degree_0.9
sum(mm_bmsc_degree_0.9)
#[1] 732

# get degree_0.9 (no. of interactions) distribution of all BMSC ligands - 0.9
all_bmsc_ligands_ex$degree_0.9 <- NA
for(i in 1:nrow(all_bmsc_ligands_ex)) {      
  protein <- all_bmsc_ligands_ex[i, ]
  protein_degree_0.9 <-subset(string_human_physical_uniprot_0.9, 
                              string_human_physical_uniprot_0.9$protein1 %in% 
                                protein$uniprotswissprot)
  
  all_bmsc_ligands_ex$degree_0.9[[i]] <- nrow(protein_degree_0.9)
  
}

# plot the distribution
all_bmsc_degree_0.9 <- as.numeric(unlist(all_bmsc_ligands_ex$degree_0.9, use.names = FALSE))
hist(all_bmsc_degree_0.9)

# sum of all BMSC degree_0.9
sum(all_bmsc_degree_0.9)
#[1] 4093

# density function 
set.seed(3333)
bmsc.dens.obs_0.9 <- density(mm_bmsc_ligands_ex$degree_0.9, n=length(all_bmsc_ligands_ex$degree_0.9), adjust = 2)

## create FOR loop for bootstrapping 
bmsc_104_1e4_receptor_interactions_0.9 <- list()
set.seed(3333)
for(i in 1:10000){
  sample <- all_bmsc_ligands_ex[sample(nrow(all_bmsc_ligands_ex), 104, replace = TRUE, prob=bmsc.dens.obs_0.9$y), ]
  interactions_sample <-subset(string_human_physical_uniprot_0.9, 
                               string_human_physical_uniprot_0.9$protein1 %in% 
                                 sample$uniprotswissprot 
                               & string_human_physical_uniprot_0.9$protein2 %in% 
                                 proteins_sensory_neurons_annot_transc[proteins_sensory_neurons_annot_transc$receptor == TRUE & proteins_sensory_neurons_annot_transc$num_sources_extracell > 0,]$uniprot_id)
  
  bmsc_104_1e4_receptor_interactions_0.9[[i]] <- nrow(interactions_sample)
}

## plot null (normal) sampling distribution
bmsc_104_1e4_receptor_interactions_0.9 <- unlist(bmsc_104_1e4_receptor_interactions_0.9, use.names = FALSE)
hist(bmsc_104_1e4_receptor_interactions_0.9)

write.csv(bmsc_104_1e4_receptor_interactions_0.9,file="data/bmsc_104_1e4_receptor_interactions_0.9.csv",row.names = FALSE)
bmsc_104_1e4_receptor_interactions_0.9 <- as.vector(t(read.csv("data/bmsc_104_1e4_receptor_interactions_0.9.csv")))

## observed interactions for DEG from MM-BMSC
mm_bmsc_104_receptor_interactions_0.9 <- subset(string_human_physical_uniprot_0.9, 
                                   string_human_physical_uniprot_0.9$protein1 %in% 
                                     mm_bmsc_ligands_ex$uniprotswissprot 
                                   & string_human_physical_uniprot_0.9$protein2 %in% 
                                     proteins_sensory_neurons_annot_transc[proteins_sensory_neurons_annot_transc$receptor == TRUE & proteins_sensory_neurons_annot_transc$num_sources_extracell > 0 ,]$uniprot_id)
nrow(mm_bmsc_104_receptor_interactions_0.9)
#[1] 288

## add observed interaction 
abline(v=288, col="blue")

## calculate p value
pvalue_bmsc_104_1e4_receptor_interactions_0.9 <- mean(bmsc_104_1e4_receptor_interactions_0.9 >= 288)
#[1] 0

## store the p value
p_value[3,1] <- mean(bmsc_104_1e4_receptor_interactions_0.9 >= 288)

# plot ggplot histogram
library(ggplot2)
library(gridExtra)
library(grid)

tiff("test.tiff", units="in", width=5, height=5, res=300)

bmsc_0.4 <- ggplot() + aes(as.numeric(bmsc_104_1e4_receptor_interactions_0.4))+ geom_histogram(aes(y=after_stat(density)), binwidth=3, colour="darkblue", fill="lightblue")+
  geom_density(alpha=0.2, fill="#FF6666") +
  ylim(0.0000, 0.015) +
  annotate("text", x = 760, y=0.010, label = "p-value = 0.0268") +
  labs(title="0.4 cutoff",x="Number of interactions", y = "Density")+
  theme_classic()+
  geom_vline(aes(xintercept=678),
             color="red", linetype="dashed", linewidth=1)

bmsc_0.4

bmsc_0.7 <- ggplot() + aes(as.numeric(bmsc_104_1e4_receptor_interactions_0.7))+ geom_histogram(aes(y=after_stat(density)), binwidth=3, colour="darkblue", fill="lightblue")+
  geom_density(alpha=0.2, fill="#FF6666") +
  ylim(0.0000, 0.015) +
  annotate("text", x = 440, y=0.010, label = "p-value = 0.0044") +
  labs(title="0.7 cutoff",x="Number of interactions", y = "Density")+
  theme_classic()+
  geom_vline(aes(xintercept=390),
             color="red", linetype="dashed", linewidth=1)

bmsc_0.7

bmsc_0.9 <- ggplot() + aes(as.numeric(bmsc_104_1e4_receptor_interactions_0.9))+ geom_histogram(aes(y=after_stat(density)), binwidth=3, colour="darkblue", fill="lightblue")+
  geom_density(alpha=0.2, fill="#FF6666") +
  ylim(0.0000, 0.015) +
  annotate("text", x = 320, y=0.010, label = "p-value = 0.001") +
  labs(title="0.9 cutoff",x="Number of interactions", y = "Density")+
  theme_classic() +
  geom_vline(aes(xintercept=288),
             color="red", linetype="dashed", linewidth=1)

bmsc_0.9

grid.arrange(bmsc_0.4, bmsc_0.7, bmsc_0.9, ncol=3)

bmsc <- arrangeGrob(bmsc_0.4, bmsc_0.7, bmsc_0.9, ncol=3)
ggsave(file="results/bmsc_all.png", bmsc) 


tiff("results/bmsc_all_0.7.tiff", units="in", width=6, height=5, res=300)

bmsc_0.7 <- ggplot() + aes(as.numeric(bmsc_104_1e4_receptor_interactions_0.7))+ geom_histogram(aes(y=after_stat(density)), binwidth=3, colour="darkblue", fill="lightblue")+
  geom_density(alpha=0.2, fill="#FF6666") +
  ylim(0.0000, 0.015) +
  annotate("text", x = 425, y=0.010, label = "p-value = 0.0044") +
  labs(title="",x="Number of interactions", y = "Density")+
  theme_classic()+
  geom_vline(aes(xintercept=390),
             color="red", linetype="dashed", linewidth=1)

bmsc_0.7

dev.off()

#_____________________________________________________________________________________
#_____________________________________________________________________________________


####PC(All)_Receptor(All)####

# check how many DEG of MM-PC genes are annotated as ligand and have at least 1 extracellular resource
nrow(merged_analysis_PC_8DS_ID_annot[merged_analysis_PC_8DS_ID_annot$ligand == TRUE & merged_analysis_PC_8DS_ID_annot$num_sources_extracell > 0, ])
#[1] 178

# create a subset dataframe that contains DEG of MM-pc extracellular ligands 
mm_pc_ligands_ex <- merged_analysis_PC_8DS_ID_annot[merged_analysis_PC_8DS_ID_annot$ligand == TRUE & merged_analysis_PC_8DS_ID_annot$num_sources_extracell > 0, ]

# get degree_0.4 (no. of interactions) distribution of DEG of MM-pc ligands - 0.4
mm_pc_ligands_ex$degree_0.4 <- NA

for(i in 1:nrow(mm_pc_ligands_ex)) {      
  protein <- mm_pc_ligands_ex[i, ]
  protein_degree_0.4 <-subset(string_human_physical_uniprot_0.4, 
                              string_human_physical_uniprot_0.4$protein1 %in% 
                                protein$uniprotswissprot)
  
  mm_pc_ligands_ex$degree_0.4[[i]] <- nrow(protein_degree_0.4)
  
}

# plot the distribution
mm_pc_degree_0.4 <-  as.numeric(unlist(mm_pc_ligands_ex$degree_0.4, use.names = FALSE))
hist(mm_pc_degree_0.4)

# sum of MM-pc degree_0.4
sum(mm_pc_degree_0.4)
#[1] 3927

# check how many all PC genes are annotated as ligand and have at least 1 extracellular resource
nrow(merged_analysis_PC_8DS_ID_all[merged_analysis_PC_8DS_ID_all$ligand == TRUE & merged_analysis_PC_8DS_ID_all$num_sources_extracell > 0, ])
#[1] 790

# create a subset dataframe that contains all of pc extracellular ligands 
all_pc_ligands_ex <- merged_analysis_PC_8DS_ID_all[merged_analysis_PC_8DS_ID_all$ligand == TRUE & merged_analysis_PC_8DS_ID_all$num_sources_extracell > 0, ]

# get degree_0.4 (no. of interactions) distribution of all pc ligands - 0.4
all_pc_ligands_ex$degree_0.4 <- NA
for(i in 1:nrow(all_pc_ligands_ex)) {      
  protein <- all_pc_ligands_ex[i, ]
  protein_degree_0.4 <-subset(string_human_physical_uniprot_0.4, 
                              string_human_physical_uniprot_0.4$protein1 %in% 
                                protein$uniprotswissprot)
  
  all_pc_ligands_ex$degree_0.4[[i]] <- nrow(protein_degree_0.4)
  
}

# plot the distribution
all_pc_degree_0.4 <- as.numeric(unlist(all_pc_ligands_ex$degree_0.4, use.names = FALSE))
hist(all_pc_degree_0.4)

# sum of all pc degree_0.4
sum(all_pc_degree_0.4)
#[1] 17989

# density function 
set.seed(3333)
pc.dens.obs_0.4 <- density(mm_pc_ligands_ex$degree_0.4, n=length(all_pc_ligands_ex$degree_0.4), adjust = 0.5)


# check how many sensory neuron genes are annotated as receptor and have at least 1 extracellular resource
nrow(proteins_sensory_neurons_annot_transc[proteins_sensory_neurons_annot_transc$receptor == TRUE & proteins_sensory_neurons_annot_transc$num_sources_extracell > 0, ])
#[1] 786

# calculate p value for observed interactions for pc at 0.4 cutoff

## filter string interactions
string_human_physical_uniprot_0.4 <- string_human_physical_uniprot[string_human_physical_uniprot$score >= 400,]

## create FOR loop for bootstrapping 
pc_178_1e4_receptor_interactions_0.4 <- list()
set.seed(3333)
for(i in 1:10000){
  sample <- all_pc_ligands_ex[sample(nrow(all_pc_ligands_ex), 178, replace = TRUE, prob=pc.dens.obs_0.4$y), ]
  interactions_sample <-subset(string_human_physical_uniprot_0.4, 
                               string_human_physical_uniprot_0.4$protein1 %in% 
                                 sample$uniprotswissprot 
                               & string_human_physical_uniprot_0.4$protein2 %in% 
                                 proteins_sensory_neurons_annot_transc[proteins_sensory_neurons_annot_transc$receptor == TRUE & proteins_sensory_neurons_annot_transc$num_sources_extracell > 0,]$uniprot_id)
  
  pc_178_1e4_receptor_interactions_0.4[[i]] <- nrow(interactions_sample)
}


## plot null (normal) sampling distribution
pc_178_1e4_receptor_interactions_0.4 <- unlist(pc_178_1e4_receptor_interactions_0.4, use.names = FALSE)
hist(pc_178_1e4_receptor_interactions_0.4)

write.csv(pc_178_1e4_receptor_interactions_0.4,file="data/pc_178_1e4_receptor_interactions_0.4.csv",row.names = FALSE)
pc_178_1e4_receptor_interactions_0.4 <- as.vector(t(read.csv("data/pc_178_1e4_receptor_interactions_0.4.csv")))

## observed interactions for DEG from MM-pc
mm_pc_178_receptor_interactions_0.4 <- subset(string_human_physical_uniprot_0.4, 
                                 string_human_physical_uniprot_0.4$protein1 %in% 
                                   mm_pc_ligands_ex$uniprotswissprot 
                                 & string_human_physical_uniprot_0.4$protein2 %in% 
                                   proteins_sensory_neurons_annot_transc[proteins_sensory_neurons_annot_transc$receptor == TRUE & proteins_sensory_neurons_annot_transc$num_sources_extracell > 0 ,]$uniprot_id)
nrow(mm_pc_178_receptor_interactions_0.4)
#[1] 1071

## add observed interaction 
abline(v=1071, col="blue")

## calculate p value
pvalue_pc_178_1e4_receptor_interactions_0.4 <- mean(pc_178_1e4_receptor_interactions_0.4>=1071)
#[1] 0

## store the p value
p_value[1,3] <- mean(pc_178_1e4_receptor_interactions_0.4>=1071)

# calculate p value for observed interactions for pc at 0.7 cutoff

## filter string interactions
string_human_physical_uniprot_0.7 <- string_human_physical_uniprot[string_human_physical_uniprot$score >= 700,]

# get degree_0.7 (no. of interactions) distribution of DEG of MM-pc ligands - 0.7
mm_pc_ligands_ex$degree_0.7 <- NA

for(i in 1:nrow(mm_pc_ligands_ex)) {      
  protein <- mm_pc_ligands_ex[i, ]
  protein_degree_0.7 <-subset(string_human_physical_uniprot_0.7, 
                              string_human_physical_uniprot_0.7$protein1 %in% 
                                protein$uniprotswissprot)
  
  mm_pc_ligands_ex$degree_0.7[[i]] <- nrow(protein_degree_0.7)
  
}

# plot the distribution
mm_pc_degree_0.7 <-  as.numeric(unlist(mm_pc_ligands_ex$degree_0.7, use.names = FALSE))
hist(mm_pc_degree_0.7)

# sum of MM-pc degree_0.7
sum(mm_pc_degree_0.7)
#[1] 1574

# get degree_0.7 (no. of interactions) distribution of all pc ligands - 0.7
all_pc_ligands_ex$degree_0.7 <- NA
for(i in 1:nrow(all_pc_ligands_ex)) {      
  protein <- all_pc_ligands_ex[i, ]
  protein_degree_0.7 <-subset(string_human_physical_uniprot_0.7, 
                              string_human_physical_uniprot_0.7$protein1 %in% 
                                protein$uniprotswissprot)
  
  all_pc_ligands_ex$degree_0.7[[i]] <- nrow(protein_degree_0.7)
  
}

# plot the distribution
all_pc_degree_0.7 <- as.numeric(unlist(all_pc_ligands_ex$degree_0.7, use.names = FALSE))
hist(all_pc_degree_0.7)

# sum of all pc degree_0.7
sum(all_pc_degree_0.7)
#[1] 6842

# density function 
set.seed(3333)
pc.dens.obs_0.7 <- density(mm_pc_ligands_ex$degree_0.7, n=length(all_pc_ligands_ex$degree_0.7))

## create FOR loop for bootstrapping 
pc_178_1e4_receptor_interactions_0.7 <- list()
set.seed(3333)
for(i in 1:10000){
  sample <- all_pc_ligands_ex[sample(nrow(all_pc_ligands_ex), 178, replace = TRUE,  prob=pc.dens.obs_0.7$y), ]
  interactions_sample <-subset(string_human_physical_uniprot_0.7, 
                               string_human_physical_uniprot_0.7$protein1 %in% 
                                 sample$uniprotswissprot 
                               & string_human_physical_uniprot_0.7$protein2 %in% 
                                 proteins_sensory_neurons_annot_transc[proteins_sensory_neurons_annot_transc$receptor == TRUE & proteins_sensory_neurons_annot_transc$num_sources_extracell > 0,]$uniprot_id)
  
  pc_178_1e4_receptor_interactions_0.7[[i]] <- nrow(interactions_sample)
}


## plot null (normal) sampling distribution
pc_178_1e4_receptor_interactions_0.7 <- unlist(pc_178_1e4_receptor_interactions_0.7, use.names = FALSE)
hist(pc_178_1e4_receptor_interactions_0.7)

write.csv(pc_178_1e4_receptor_interactions_0.7,file="data/pc_178_1e4_receptor_interactions_0.7.csv",row.names = FALSE)
pc_178_1e4_receptor_interactions_0.7 <- as.vector(t(read.csv("data/pc_178_1e4_receptor_interactions_0.7.csv")))

## observed interactions for DEG from MM-pc
mm_pc_178_receptor_interactions_0.7 <- subset(string_human_physical_uniprot_0.7, 
                                 string_human_physical_uniprot_0.7$protein1 %in% 
                                   mm_pc_ligands_ex$uniprotswissprot 
                                 & string_human_physical_uniprot_0.7$protein2 %in% 
                                   proteins_sensory_neurons_annot_transc[proteins_sensory_neurons_annot_transc$receptor == TRUE & proteins_sensory_neurons_annot_transc$num_sources_extracell > 0 ,]$uniprot_id)
nrow(mm_pc_178_receptor_interactions_0.7)
#[1] 527

## add observed interaction 
abline(v=527, col="blue")

## calculate p value
pvalue_pc_178_1e4_receptor_interactions_0.7 <- mean(pc_178_1e4_receptor_interactions_0.7>=527)
#[1] 0

## store the p value
p_value[2,3] <- mean(pc_178_1e4_receptor_interactions_0.7>=527)

# calculate p value for observed interactions for pc at 0.9 cutoff ######

## filter string interactions
string_human_physical_uniprot_0.9 <- string_human_physical_uniprot[string_human_physical_uniprot$score >= 900,]

# get degree_0.9 (no. of interactions) distribution of DEG of MM-pc ligands - 0.9
mm_pc_ligands_ex$degree_0.9 <- NA

for(i in 1:nrow(mm_pc_ligands_ex)) {      
  protein <- mm_pc_ligands_ex[i, ]
  protein_degree_0.9 <-subset(string_human_physical_uniprot_0.9, 
                              string_human_physical_uniprot_0.9$protein1 %in% 
                                protein$uniprotswissprot)
  
  mm_pc_ligands_ex$degree_0.9[[i]] <- nrow(protein_degree_0.9)
  
}

# plot the distribution
mm_pc_degree_0.9 <-  as.numeric(unlist(mm_pc_ligands_ex$degree_0.9, use.names = FALSE))
hist(mm_pc_degree_0.9)

# sum of MM-pc degree_0.9
sum(mm_pc_degree_0.9)
#[1] 953

# get degree_0.9 (no. of interactions) distribution of all pc ligands - 0.9
all_pc_ligands_ex$degree_0.9 <- NA
for(i in 1:nrow(all_pc_ligands_ex)) {      
  protein <- all_pc_ligands_ex[i, ]
  protein_degree_0.9 <-subset(string_human_physical_uniprot_0.9, 
                              string_human_physical_uniprot_0.9$protein1 %in% 
                                protein$uniprotswissprot)
  
  all_pc_ligands_ex$degree_0.9[[i]] <- nrow(protein_degree_0.9)
  
}

# plot the distribution
all_pc_degree_0.9 <- as.numeric(unlist(all_pc_ligands_ex$degree_0.9, use.names = FALSE))
hist(all_pc_degree_0.9)

# sum of all pc degree_0.9
sum(all_pc_degree_0.9)
#[1] 4168

# density function 
set.seed(3333)
pc.dens.obs_0.9 <- density(mm_pc_ligands_ex$degree_0.9, n=length(all_pc_ligands_ex$degree_0.9),  adjust = 1.5)


## create FOR loop for bootstrapping 
pc_178_1e4_receptor_interactions_0.9 <- list()
set.seed(3333)
for(i in 1:10000){
  sample <- all_pc_ligands_ex[sample(nrow(all_pc_ligands_ex), 178, replace = TRUE, pc.dens.obs_0.9$y), ]
  interactions_sample <-subset(string_human_physical_uniprot_0.9, 
                               string_human_physical_uniprot_0.9$protein1 %in% 
                                 sample$uniprotswissprot 
                               & string_human_physical_uniprot_0.9$protein2 %in% 
                                 proteins_sensory_neurons_annot_transc[proteins_sensory_neurons_annot_transc$receptor == TRUE & proteins_sensory_neurons_annot_transc$num_sources_extracell > 0,]$uniprot_id)
  
  pc_178_1e4_receptor_interactions_0.9[[i]] <- nrow(interactions_sample)
}


## plot null (normal) sampling distribution
pc_178_1e4_receptor_interactions_0.9 <- unlist(pc_178_1e4_receptor_interactions_0.9, use.names = FALSE)
hist(pc_178_1e4_receptor_interactions_0.9)

write.csv(pc_178_1e4_receptor_interactions_0.9,file="data/pc_178_1e4_receptor_interactions_0.9.csv",row.names = FALSE)
pc_178_1e4_receptor_interactions_0.9 <- as.vector(t(read.csv("data/pc_178_1e4_receptor_interactions_0.9.csv")))

## observed interactions for DEG from MM-pc
mm_pc_178_receptor_interactions_0.9 <- subset(string_human_physical_uniprot_0.9, 
                                 string_human_physical_uniprot_0.9$protein1 %in% 
                                   mm_pc_ligands_ex$uniprotswissprot 
                                 & string_human_physical_uniprot_0.9$protein2 %in% 
                                   proteins_sensory_neurons_annot_transc[proteins_sensory_neurons_annot_transc$receptor == TRUE & proteins_sensory_neurons_annot_transc$num_sources_extracell > 0 ,]$uniprot_id)
nrow(mm_pc_178_receptor_interactions_0.9)
#[1] 388

## add observed interaction 
abline(v=388, col="blue")

## calculate p value
pvalue_pc_178_1e4_receptor_interactions_0.9 <- mean(pc_178_1e4_receptor_interactions_0.9>=388)
#[1]0

## store the p value
p_value[3,3] <- mean(pc_178_1e4_receptor_interactions_0.9>=388)

# plot ggplot histogram

tiff("test.tiff", units="in", width=15, height=10, res=300)

pc_0.4 <- ggplot() + aes(as.numeric(pc_178_1e4_receptor_interactions_0.4))+ geom_histogram(aes(y=after_stat(density)), binwidth=3, colour="darkblue", fill="lightblue")+
  geom_density(alpha=0.2, fill="#FF6666") +
  ylim(0.0000, 0.017) +
  annotate("text", x = 850, y=0.010, label = "p-value=0") +
  labs(title="0.4 cutoff",x="Number of interactions", y = "Density")+
  theme_classic()+
  geom_vline(aes(xintercept=1071),
                   color="red", linetype="dashed", linewidth=1)

pc_0.4

pc_0.7 <- ggplot() + aes(as.numeric(pc_178_1e4_receptor_interactions_0.7))+ geom_histogram(aes(y=after_stat(density)), binwidth=3, colour="darkblue", fill="lightblue")+
  geom_density(alpha=0.2, fill="#FF6666") +
  ylim(0.0000, 0.017) +
  annotate("text", x = 470, y=0.010, label = "p-value=0") +
  labs(title="0.7 cutoff",x="Number of interactions", y = "Density")+
  theme_classic()+
  geom_vline(aes(xintercept=527),
                   color="red", linetype="dashed", linewidth=1)

pc_0.7

pc_0.9 <- ggplot() + aes(as.numeric(pc_178_1e4_receptor_interactions_0.9))+ geom_histogram(aes(y=after_stat(density)), binwidth=3, colour="darkblue", fill="lightblue")+
  geom_density(alpha=0.2, fill="#FF6666") +
  ylim(0.0000, 0.017) +
  annotate("text", x = 325, y=0.010, label = "p-value=0") +
  labs(title="0.9 cutoff",x="Number of interactions", y = "Density")+
  theme_classic() +
  geom_vline(aes(xintercept=388),
              color="red", linetype="dashed", linewidth=1)

pc_0.9

grid.arrange(pc_0.4, pc_0.7, pc_0.9, ncol=3)

g <- arrangeGrob(pc_0.4, pc_0.7, pc_0.9, ncol=3)
ggsave(file="results/pc_all.png", g) 

tiff("results/pc_all_0.7.tiff", units="in", width=6, height=5, res=300)

pc_0.7 <- ggplot() + aes(as.numeric(pc_178_1e4_receptor_interactions_0.7))+ geom_histogram(aes(y=after_stat(density)), binwidth=3, colour="darkblue", fill="lightblue")+
  geom_density(alpha=0.2, fill="#FF6666") +
  ylim(0.0000, 0.015) +
  annotate("text", x = 425, y=0.01, label = "p-value = 0") +
  labs(title="",x="Number of interactions", y = "Density")+
  theme_classic()+
  geom_vline(aes(xintercept=527),
             color="red", linetype="dashed", linewidth=1)

pc_0.7

dev.off()
#_____________________________________________________________________________________
#_____________________________________________________________________________________

####BMSC(All)_Receptor(Pain)#### 

# check how many DEG of MM-BMSC genes are annotated as ligand and have at least 1 extracellular resource
nrow(merged_analysis_BMSC_5DS_ID_proteome_annot[merged_analysis_BMSC_5DS_ID_proteome_annot$ligand == TRUE & merged_analysis_BMSC_5DS_ID_proteome_annot$num_sources_extracell > 0, ])
#[1] 104

# check how many all BMSC genes are annotated as ligand and have at least 1 extracellular resource
nrow(merged_analysis_BMSC_5DS_ID_all[merged_analysis_BMSC_5DS_ID_all$ligand == TRUE & merged_analysis_BMSC_5DS_ID_all$num_sources_extracell > 0, ])
#[1] 763

# check how many sensory neuron genes are annotated as receptor and have at least 1 extracellular resource and have at least 1 "pain" annotation
nrow(proteins_sensory_neurons_annot_transc[proteins_sensory_neurons_annot_transc$receptor == TRUE & proteins_sensory_neurons_annot_transc$num_sources_extracell > 0 & proteins_sensory_neurons_annot_transc$num_sources_pain > 0, ])
#[1] 384

# create a subset dataframe that contains all BMSC extracellular ligands to sample from
all_bmsc_ligands_ex <- merged_analysis_BMSC_5DS_ID_all[merged_analysis_BMSC_5DS_ID_all$ligand == TRUE & merged_analysis_BMSC_5DS_ID_all$num_sources_extracell > 0, ]

# calculate p value for observed interactions for BMSC at 0.4 cutoff

## filter string interactions
string_human_physical_uniprot_0.4 <- string_human_physical_uniprot[string_human_physical_uniprot$score >= 400,]

## create FOR loop for bootstrapping 
bmsc_104_1e4_pain_receptor_interactions_0.4 <- list()
set.seed(3333)
for(i in 1:10000){
  sample <- all_bmsc_ligands_ex[sample(nrow(all_bmsc_ligands_ex), 104, replace = TRUE, prob=bmsc.dens.obs_0.4$y), ]
  interactions_sample <-subset(string_human_physical_uniprot_0.4, 
                               string_human_physical_uniprot_0.4$protein1 %in% 
                                 sample$uniprotswissprot 
                               & string_human_physical_uniprot_0.4$protein2 %in% 
                                 proteins_sensory_neurons_annot_transc[proteins_sensory_neurons_annot_transc$receptor == TRUE & proteins_sensory_neurons_annot_transc$num_sources_extracell > 0 & proteins_sensory_neurons_annot_transc$num_sources_pain > 0,]$uniprot_id)
  
  bmsc_104_1e4_pain_receptor_interactions_0.4[[i]] <- nrow(interactions_sample)
}


## plot null (normal) sampling distribution
bmsc_104_1e4_pain_receptor_interactions_0.4 <- unlist(bmsc_104_1e4_pain_receptor_interactions_0.4, use.names = FALSE)
hist(bmsc_104_1e4_pain_receptor_interactions_0.4)

write.csv(bmsc_104_1e4_pain_receptor_interactions_0.4,file="data/bmsc_104_1e4_pain_receptor_interactions_0.4.csv",row.names = FALSE)
bmsc_104_1e4_pain_receptor_interactions_0.4 <- as.vector(t(read.csv("data/bmsc_104_1e4_pain_receptor_interactions_0.4.csv")))

## observed interactions for DEG from MM-BMSC
mm_bmsc_104_pain_receptor_interactions_0.4 <- subset(string_human_physical_uniprot_0.4, 
                                   string_human_physical_uniprot_0.4$protein1 %in% 
                                     merged_analysis_BMSC_5DS_ID_proteome_annot[merged_analysis_BMSC_5DS_ID_proteome_annot$ligand == TRUE & merged_analysis_BMSC_5DS_ID_proteome_annot$num_sources_extracell > 0, ]$uniprotswissprot 
                                   & string_human_physical_uniprot_0.4$protein2 %in% 
                                     proteins_sensory_neurons_annot_transc[proteins_sensory_neurons_annot_transc$receptor == TRUE & proteins_sensory_neurons_annot_transc$num_sources_extracell > 0 & proteins_sensory_neurons_annot_transc$num_sources_pain > 0,]$uniprot_id)
nrow(mm_bmsc_104_pain_receptor_interactions_0.4)
#[1] 411

## add observed interaction 
abline(v=411, col="blue")

## calculate p value
pvalue_bmsc_104_1e4_pain_receptor_interactions_0.4 <- mean(bmsc_104_1e4_pain_receptor_interactions_0.4>=411)
#[1] 0.021

## store the p value
p_value[1,2] <- mean(bmsc_104_1e4_pain_receptor_interactions_0.4>=411)

# calculate p value for observed interactions for BMSC at 0.7 cutoff

## filter string interactions
string_human_physical_uniprot_0.7 <- string_human_physical_uniprot[string_human_physical_uniprot$score >= 700,]

## create FOR loop for bootstrapping 
bmsc_104_1e4_pain_receptor_interactions_0.7 <- list()
set.seed(3333)
for(i in 1:10000){
  sample <- all_bmsc_ligands_ex[sample(nrow(all_bmsc_ligands_ex), 104, replace = TRUE, prob=bmsc.dens.obs_0.7$y), ]
  interactions_sample <-subset(string_human_physical_uniprot_0.7, 
                               string_human_physical_uniprot_0.7$protein1 %in% 
                                 sample$uniprotswissprot 
                               & string_human_physical_uniprot_0.7$protein2 %in% 
                                 proteins_sensory_neurons_annot_transc[proteins_sensory_neurons_annot_transc$receptor == TRUE & proteins_sensory_neurons_annot_transc$num_sources_extracell > 0 & proteins_sensory_neurons_annot_transc$num_sources_pain > 0,]$uniprot_id)
  
  bmsc_104_1e4_pain_receptor_interactions_0.7[[i]] <- nrow(interactions_sample)
}


## plot null (normal) sampling distribution
bmsc_104_1e4_pain_receptor_interactions_0.7 <- unlist(bmsc_104_1e4_pain_receptor_interactions_0.7, use.names = FALSE)
hist(bmsc_104_1e4_pain_receptor_interactions_0.7)

write.csv(bmsc_104_1e4_pain_receptor_interactions_0.7,file="data/bmsc_104_1e4_pain_receptor_interactions_0.7.csv",row.names = FALSE)
bmsc_104_1e4_pain_receptor_interactions_0.7 <- as.vector(t(read.csv("data/bmsc_104_1e4_pain_receptor_interactions_0.7.csv")))

## observed interactions for DEG from MM-BMSC
mm_bmsc_104_pain_receptor_interactions_0.7 <- subset(string_human_physical_uniprot_0.7, 
                                   string_human_physical_uniprot_0.7$protein1 %in% 
                                     merged_analysis_BMSC_5DS_ID_proteome_annot[merged_analysis_BMSC_5DS_ID_proteome_annot$ligand == TRUE & merged_analysis_BMSC_5DS_ID_proteome_annot$num_sources_extracell > 0, ]$uniprotswissprot 
                                   & string_human_physical_uniprot_0.7$protein2 %in% 
                                     proteins_sensory_neurons_annot_transc[proteins_sensory_neurons_annot_transc$receptor == TRUE & proteins_sensory_neurons_annot_transc$num_sources_extracell > 0 & proteins_sensory_neurons_annot_transc$num_sources_pain > 0,]$uniprot_id)
nrow(mm_bmsc_104_pain_receptor_interactions_0.7)
#[1] 218

## add observed interaction 
abline(v=218, col="blue")

## calculate p value
pvalue_bmsc_104_1e4_pain_receptor_interactions_0.7 <- mean(bmsc_104_1e4_pain_receptor_interactions_0.7>=218)
#[1] 0.0349

## store the p value
p_value[2,2] <- mean(bmsc_104_1e4_pain_receptor_interactions_0.7>=218)

# calculate p value for observed interactions for BMSC at 0.9 cutoff ######

## filter string interactions
string_human_physical_uniprot_0.9 <- string_human_physical_uniprot[string_human_physical_uniprot$score >= 900,]

## create FOR loop for bootstrapping 
bmsc_104_1e4_pain_receptor_interactions_0.9 <- list()
set.seed(3333)
for(i in 1:10000){
  sample <- all_bmsc_ligands_ex[sample(nrow(all_bmsc_ligands_ex), 104, replace = TRUE, prob=bmsc.dens.obs_0.9$y), ]
  interactions_sample <-subset(string_human_physical_uniprot_0.9, 
                               string_human_physical_uniprot_0.9$protein1 %in% 
                                 sample$uniprotswissprot 
                               & string_human_physical_uniprot_0.9$protein2 %in% 
                                 proteins_sensory_neurons_annot_transc[proteins_sensory_neurons_annot_transc$receptor == TRUE & proteins_sensory_neurons_annot_transc$num_sources_extracell > 0 & proteins_sensory_neurons_annot_transc$num_sources_pain > 0,]$uniprot_id)
  
  bmsc_104_1e4_pain_receptor_interactions_0.9[[i]] <- nrow(interactions_sample)
}


## plot null (normal) sampling distribution
bmsc_104_1e4_pain_receptor_interactions_0.9 <- unlist(bmsc_104_1e4_pain_receptor_interactions_0.9, use.names = FALSE)
hist(bmsc_104_1e4_pain_receptor_interactions_0.9)

write.csv(bmsc_104_1e4_pain_receptor_interactions_0.9,file="data/bmsc_104_1e4_pain_receptor_interactions_0.9.csv",row.names = FALSE)
bmsc_104_1e4_pain_receptor_interactions_0.9 <- as.vector(t(read.csv("data/bmsc_104_1e4_pain_receptor_interactions_0.9.csv")))

## observed interactions for DEG from MM-BMSC
mm_bmsc_104_pain_receptor_interactions_0.9 <- subset(string_human_physical_uniprot_0.9, 
                                   string_human_physical_uniprot_0.9$protein1 %in% 
                                     merged_analysis_BMSC_5DS_ID_proteome_annot[merged_analysis_BMSC_5DS_ID_proteome_annot$ligand == TRUE & merged_analysis_BMSC_5DS_ID_proteome_annot$num_sources_extracell > 0, ]$uniprotswissprot 
                                   & string_human_physical_uniprot_0.9$protein2 %in% 
                                     proteins_sensory_neurons_annot_transc[proteins_sensory_neurons_annot_transc$receptor == TRUE & proteins_sensory_neurons_annot_transc$num_sources_extracell > 0 & proteins_sensory_neurons_annot_transc$num_sources_pain > 0,]$uniprot_id)
nrow(mm_bmsc_104_pain_receptor_interactions_0.9)
#[1] 159

## add observed interaction 
abline(v=159, col="blue")

## calculate p value
pvalue_bmsc_104_1e4_pain_receptor_interactions_0.9 <- mean(bmsc_104_1e4_pain_receptor_interactions_0.9 >= 159)
#[1] 0.0015

## store the p value
p_value[3,2] <- mean(bmsc_104_1e4_pain_receptor_interactions_0.9 >= 159)

# plot ggplot histogram

tiff("test.tiff", units="in", width=15, height=8, res=300)

bmsc_0.4_pain <- ggplot() + aes(as.numeric(bmsc_104_1e4_pain_receptor_interactions_0.4))+ geom_histogram(aes(y=after_stat(density)), binwidth=3, colour="darkblue", fill="lightblue")+
                  geom_density(alpha=0.2, fill="#FF6666") +
                  ylim(0.0000, 0.028) +
                  annotate("text", x = 450, y=0.01, label = "p-value = 0.021") +
                  labs(title="0.4 cutoff",x="Number of interactions", y = "Density")+
                  theme_classic()+
                  geom_vline(aes(xintercept=411),
                             color="red", linetype="dashed", linewidth=1)

bmsc_0.4_pain

bmsc_0.7_pain <- ggplot() + aes(as.numeric(bmsc_104_1e4_pain_receptor_interactions_0.7))+ geom_histogram(aes(y=after_stat(density)), binwidth=3, colour="darkblue", fill="lightblue")+
                  geom_density(alpha=0.2, fill="#FF6666") +
                  ylim(0.0000, 0.028) +
                  annotate("text", x = 250, y=0.01, label = "p-value = 0.0349") +
                  labs(title="0.7 cutoff",x="Number of interactions", y = "Density")+
                  theme_classic()+
                  geom_vline(aes(xintercept=218),
                             color="red", linetype="dashed", linewidth=1)

bmsc_0.7_pain

bmsc_0.9_pain <- ggplot() + aes(as.numeric(bmsc_104_1e4_pain_receptor_interactions_0.9))+ geom_histogram(aes(y=after_stat(density)), binwidth=3, colour="darkblue", fill="lightblue")+
                geom_density(alpha=0.22, fill="#FF6666") +
                ylim(0.0000, 0.028) +
                annotate("text", x = 150, y=0.01, label = "p-value = 0.0015") +
                labs(title="0.9 cutoff",x="Number of interactions", y = "Density")+
                theme_classic() +
                geom_vline(aes(xintercept=159),
                           color="red", linetype="dashed", linewidth=1)

bmsc_0.9_pain


grid.arrange(bmsc_0.4_pain, bmsc_0.7_pain, bmsc_0.9_pain, ncol=3)

bmsc <- arrangeGrob(bmsc_0.4_pain, bmsc_0.7_pain, bmsc_0.9_pain, ncol=3)
ggsave(file="results/bmsc_all_pain.png", bmsc) 


tiff("results/bmsc_all_pain_0.7.tiff", units="in", width=11.5, height=3, res=300)

bmsc_0.7_pain <- ggplot() + aes(as.numeric(bmsc_104_1e4_pain_receptor_interactions_0.7))+ geom_histogram(aes(y=after_stat(density)), binwidth=3, colour="darkblue", fill="lightblue")+
  geom_density(alpha=0.2, fill="#FF6666") +
  ylim(0.0000, 0.028) +
  annotate("text", x = 250, y=0.01, label = "p-value = 0.0349") +
  labs(title="",x="Number of interactions", y = "Density")+
  theme_classic()+
  geom_vline(aes(xintercept=218),
             color="red", linetype="dashed", linewidth=1)

bmsc_0.7_pain

dev.off()
#_____________________________________________________________________________________
#_____________________________________________________________________________________

####PC(All)_Receptor(Pain)####

# check how many DEG of MM-PC genes are annotated as ligand and have at least 1 extracellular resource
nrow(merged_analysis_PC_8DS_ID_annot[merged_analysis_PC_8DS_ID_annot$ligand == TRUE & merged_analysis_PC_8DS_ID_annot$num_sources_extracell > 0, ])
#[1] 178

# check how many all PC genes are annotated as ligand and have at least 1 extracellular resource
nrow(merged_analysis_PC_8DS_ID_all[merged_analysis_PC_8DS_ID_all$ligand == TRUE & merged_analysis_PC_8DS_ID_all$num_sources_extracell > 0, ])
#[1] 790

# check how many sensory neuron genes are annotated as receptor and have at least 1 extracellular resource
nrow(proteins_sensory_neurons_annot_transc[proteins_sensory_neurons_annot_transc$receptor == TRUE & proteins_sensory_neurons_annot_transc$num_sources_extracell > 0 & proteins_sensory_neurons_annot_transc$num_sources_pain > 0, ])
#[1] 384

# create a subset dataframe that contains all PC extracellular ligands to sample from
all_pc_ligands_ex <- merged_analysis_PC_8DS_ID_all[merged_analysis_PC_8DS_ID_all$ligand == TRUE & merged_analysis_PC_8DS_ID_all$num_sources_extracell > 0, ]

# calculate p value for observed interactions for pc at 0.4 cutoff

## filter string interactions
string_human_physical_uniprot_0.4 <- string_human_physical_uniprot[string_human_physical_uniprot$score >= 400,]

## create FOR loop for bootstrapping 
pc_178_1e4_pain_receptor_interactions_0.4 <- list()
set.seed(3333)
for(i in 1:10000){
  sample <- all_pc_ligands_ex[sample(nrow(all_pc_ligands_ex), 178, replace = TRUE, prob=pc.dens.obs_0.4$y), ]
  interactions_sample <-subset(string_human_physical_uniprot_0.4, 
                               string_human_physical_uniprot_0.4$protein1 %in% 
                                 sample$uniprotswissprot 
                               & string_human_physical_uniprot_0.4$protein2 %in% 
                                 proteins_sensory_neurons_annot_transc[proteins_sensory_neurons_annot_transc$receptor == TRUE & proteins_sensory_neurons_annot_transc$num_sources_extracell > 0 & proteins_sensory_neurons_annot_transc$num_sources_pain > 0 ,]$uniprot_id)
  
  pc_178_1e4_pain_receptor_interactions_0.4[[i]] <- nrow(interactions_sample)
}


## plot null (normal) sampling distribution
pc_178_1e4_pain_receptor_interactions_0.4 <- unlist(pc_178_1e4_pain_receptor_interactions_0.4, use.names = FALSE)
hist(pc_178_1e4_pain_receptor_interactions_0.4)

write.csv(pc_178_1e4_pain_receptor_interactions_0.4,file="data/pc_178_1e4_pain_receptor_interactions_0.4.csv",row.names = FALSE)
pc_178_1e4_pain_receptor_interactions_0.4 <- as.vector(t(read.csv("data/pc_178_1e4_pain_receptor_interactions_0.4.csv")))

## observed interactions for DEG from MM-pc
mm_pc_178_pain_receptor_interactions_0.4 <- subset(string_human_physical_uniprot_0.4, 
                                 string_human_physical_uniprot_0.4$protein1 %in% 
                                   merged_analysis_PC_8DS_ID_annot[merged_analysis_PC_8DS_ID_annot$ligand == TRUE & merged_analysis_PC_8DS_ID_annot$num_sources_extracell > 0, ]$uniprotswissprot 
                                 & string_human_physical_uniprot_0.4$protein2 %in% 
                                   proteins_sensory_neurons_annot_transc[proteins_sensory_neurons_annot_transc$receptor == TRUE & proteins_sensory_neurons_annot_transc$num_sources_extracell > 0 & proteins_sensory_neurons_annot_transc$num_sources_pain > 0,]$uniprot_id)
nrow(mm_pc_178_pain_receptor_interactions_0.4)
#[1] 582

## add observed interaction 
abline(v=582, col="blue")

## calculate p value
pvalue_pc_178_1e4_pain_receptor_interactions_0.4 <- mean(pc_178_1e4_pain_receptor_interactions_0.4>=582)
#[1] 3e-04

## store the p value
p_value[1,4] <- mean(pc_178_1e4_pain_receptor_interactions_0.4>=582)

# calculate p value for observed interactions for pc at 0.7 cutoff

## filter string interactions
string_human_physical_uniprot_0.7 <- string_human_physical_uniprot[string_human_physical_uniprot$score >= 700,]

## create FOR loop for bootstrapping 
pc_178_1e4_pain_receptor_interactions_0.7 <- list()
set.seed(3333)
for(i in 1:10000){
  sample <- all_pc_ligands_ex[sample(nrow(all_pc_ligands_ex), 178, replace = TRUE,  prob=pc.dens.obs_0.7$y), ]
  interactions_sample <-subset(string_human_physical_uniprot_0.7, 
                               string_human_physical_uniprot_0.7$protein1 %in% 
                                 sample$uniprotswissprot 
                               & string_human_physical_uniprot_0.7$protein2 %in% 
                                 proteins_sensory_neurons_annot_transc[proteins_sensory_neurons_annot_transc$receptor == TRUE & proteins_sensory_neurons_annot_transc$num_sources_extracell > 0 & proteins_sensory_neurons_annot_transc$num_sources_pain > 0,]$uniprot_id)
  
  pc_178_1e4_pain_receptor_interactions_0.7[[i]] <- nrow(interactions_sample)
}


## plot null (normal) sampling distribution
pc_178_1e4_pain_receptor_interactions_0.7 <- unlist(pc_178_1e4_pain_receptor_interactions_0.7, use.names = FALSE)
hist(pc_178_1e4_pain_receptor_interactions_0.7)

write.csv(pc_178_1e4_pain_receptor_interactions_0.7,file="data/pc_178_1e4_pain_receptor_interactions_0.7.csv",row.names = FALSE)
pc_178_1e4_pain_receptor_interactions_0.7 <- as.vector(t(read.csv("data/pc_178_1e4_pain_receptor_interactions_0.7.csv")))


## observed interactions for DEG from MM-pc
mm_pc_178_pain_receptor_interactions_0.7 <- subset(string_human_physical_uniprot_0.7, 
                                 string_human_physical_uniprot_0.7$protein1 %in% 
                                   merged_analysis_PC_8DS_ID_annot[merged_analysis_PC_8DS_ID_annot$ligand == TRUE & merged_analysis_PC_8DS_ID_annot$num_sources_extracell > 0, ]$uniprotswissprot 
                                 & string_human_physical_uniprot_0.7$protein2 %in% 
                                   proteins_sensory_neurons_annot_transc[proteins_sensory_neurons_annot_transc$receptor == TRUE & proteins_sensory_neurons_annot_transc$num_sources_extracell > 0 & proteins_sensory_neurons_annot_transc$num_sources_pain > 0,]$uniprot_id)
nrow(mm_pc_178_pain_receptor_interactions_0.7)
#[1] 303

## add observed interaction 
abline(v=303, col="blue")

## calculate p value
pvalue_pc_178_1e4_pain_receptor_interactions_0.7 <- mean(pc_178_1e4_pain_receptor_interactions_0.7>=303)
#[1] 5e-04

## store the p value
p_value[2,4] <-  mean(pc_178_1e4_pain_receptor_interactions_0.7>=303)

# calculate p value for observed interactions for pc at 0.9 cutoff ######

## filter string interactions
string_human_physical_uniprot_0.9 <- string_human_physical_uniprot[string_human_physical_uniprot$score >= 900,]

## create FOR loop for bootstrapping 
pc_178_1e4_pain_receptor_interactions_0.9 <- list()
set.seed(3333)
for(i in 1:10000){
  sample <- all_pc_ligands_ex[sample(nrow(all_pc_ligands_ex), 178, replace = TRUE, pc.dens.obs_0.9$y), ]
  interactions_sample <-subset(string_human_physical_uniprot_0.9, 
                               string_human_physical_uniprot_0.9$protein1 %in% 
                                 sample$uniprotswissprot 
                               & string_human_physical_uniprot_0.9$protein2 %in% 
                                 proteins_sensory_neurons_annot_transc[proteins_sensory_neurons_annot_transc$receptor == TRUE & proteins_sensory_neurons_annot_transc$num_sources_extracell > 0 & proteins_sensory_neurons_annot_transc$num_sources_pain > 0,]$uniprot_id)
  
  pc_178_1e4_pain_receptor_interactions_0.9[[i]] <- nrow(interactions_sample)
}


## plot null (normal) sampling distribution
pc_178_1e4_pain_receptor_interactions_0.9 <- unlist(pc_178_1e4_pain_receptor_interactions_0.9, use.names = FALSE)
hist(pc_178_1e4_pain_receptor_interactions_0.9)

write.csv(pc_178_1e4_pain_receptor_interactions_0.9,file="data/pc_178_1e4_pain_receptor_interactions_0.9.csv",row.names = FALSE)
pc_178_1e4_pain_receptor_interactions_0.9 <- as.vector(t(read.csv("data/pc_178_1e4_pain_receptor_interactions_0.9.csv")))

## observed interactions for DEG from MM-pc
mm_pc_178_pain_receptor_interactions_0.9 <- subset(string_human_physical_uniprot_0.9, 
                                 string_human_physical_uniprot_0.9$protein1 %in% 
                                   merged_analysis_PC_8DS_ID_annot[merged_analysis_PC_8DS_ID_annot$ligand == TRUE & merged_analysis_PC_8DS_ID_annot$num_sources_extracell > 0, ]$uniprotswissprot 
                                 & string_human_physical_uniprot_0.9$protein2 %in% 
                                   proteins_sensory_neurons_annot_transc[proteins_sensory_neurons_annot_transc$receptor == TRUE & proteins_sensory_neurons_annot_transc$num_sources_extracell > 0 & proteins_sensory_neurons_annot_transc$num_sources_pain > 0,]$uniprot_id)
nrow(mm_pc_178_pain_receptor_interactions_0.9)
#[1] 196

## add observed interaction 
abline(v=196, col="blue")

## calculate p value
pvalue_pc_178_1e4_pain_receptor_interactions_0.9 <- mean(pc_178_1e4_pain_receptor_interactions_0.9 >= 196)
#[1] 0.0049

## store the p value
p_value[3,4] <- mean(pc_178_1e4_pain_receptor_interactions_0.9 >= 196)

# plot ggplot histogram

tiff("test.tiff", units="in", width=15, height=8, res=300)

pc_0.4_pain <- ggplot() + aes(as.numeric(pc_178_1e4_pain_receptor_interactions_0.4))+ geom_histogram(aes(y=after_stat(density)), binwidth=3, colour="darkblue", fill="lightblue")+
          geom_density(alpha=0.2, fill="#FF6666") +
          ylim(0.0000, 0.025) +
          annotate("text", x = 500, y=0.010, label = "p-value=0.0003") +
          labs(title="0.4 cutoff",x="Number of interactions", y = "Density")+
          theme_classic()+
          geom_vline(aes(xintercept=582),
                     color="red", linetype="dashed", linewidth=1)

pc_0.4_pain

pc_0.7_pain <- ggplot() + aes(as.numeric(pc_178_1e4_pain_receptor_interactions_0.7))+ geom_histogram(aes(y=after_stat(density)), binwidth=3, colour="darkblue", fill="lightblue")+
          geom_density(alpha=0.2, fill="#FF6666") +
          ylim(0.0000, 0.025) +
          annotate("text", x = 300, y=0.010, label = "p-value=0.0005") +
          labs(title="0.7 cutoff",x="Number of interactions", y = "Density")+
          theme_classic()+
          geom_vline(aes(xintercept=303),
                     color="red", linetype="dashed", linewidth=1)

pc_0.7_pain

pc_0.9_pain <- ggplot() + aes(as.numeric(pc_178_1e4_pain_receptor_interactions_0.9))+ geom_histogram(aes(y=after_stat(density)), binwidth=3, colour="darkblue", fill="lightblue")+
          geom_density(alpha=0.2, fill="#FF6666") +
          ylim(0.0000, 0.025) +
          annotate("text", x = 200, y=0.010, label = "p-value=0.0049") +
          labs(title="0.9 cutoff",x="Number of interactions", y = "Density")+
          theme_classic() +
          geom_vline(aes(xintercept=196),
                     color="red", linetype="dashed", linewidth=1)

pc_0.9_pain

grid.arrange(pc_0.4_pain, pc_0.7_pain, pc_0.9_pain, ncol=3)

g <- arrangeGrob(pc_0.4_pain, pc_0.7_pain, pc_0.9_pain, ncol=3)
ggsave(file="results/pc_all_pain.png", g) 

tiff("results/pc_all_pain_0.7.tiff", units="in", width=11.5, height=3, res=300)

pc_0.7_pain <- ggplot() + aes(as.numeric(pc_178_1e4_pain_receptor_interactions_0.7))+ geom_histogram(aes(y=after_stat(density)), binwidth=3, colour="darkblue", fill="lightblue")+
              geom_density(alpha=0.2, fill="#FF6666") +
              ylim(0.0000, 0.025) +
              annotate("text", x = 280, y=0.010, label = "p-value=0.0005") +
              labs(title="",x="Number of interactions", y = "Density")+
              theme_classic()+
              geom_vline(aes(xintercept=303),
              color="red", linetype="dashed", linewidth=1)

pc_0.7_pain

dev.off()
#_____________________________________________________________________________________
#_____________________________________________________________________________________


write.csv(p_value,file="results/p-value.csv")
