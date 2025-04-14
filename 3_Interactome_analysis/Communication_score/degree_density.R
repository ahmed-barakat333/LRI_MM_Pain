# the aim of this code is to try the density function to control the random sampling and show if it is working 

# set working directory
setwd("~/Work/PhD/Projects/MM Pain/Computational Analysis/Communication_Score")

# load data
load("data/data_string_PC_BMSC_220613.RData")

## bmsc ligands

# load string 0.4 data
string_human_physical_uniprot_0.4 <- string_human_physical_uniprot[string_human_physical_uniprot$score >= 400,]

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

# create a subset dataframe that contains all of BMSC extracellular ligands 
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

# density function without adjustment
set.seed(3333)
bmsc.dens.obs_0.4 <- density(mm_bmsc_ligands_ex$degree_0.4, n=length(all_bmsc_ligands_ex$degree_0.4))
plot(bmsc.dens.obs_0.4$y)

# if the density functions works, the distribution of random sample degree_0.4 sum should be centered around 3000 (MM-BMSC degree_0.4)

# random sampling with density function - checking the sum of degree for all samples
degree_0.4_weighted <- list()
for(i in 1:10000) {      
  sample_weighted <- sample(all_bmsc_ligands_ex$degree_0.4, length(mm_bmsc_ligands_ex$degree_0.4), replace=TRUE, prob=bmsc.dens.obs_0.4$y)
  
  degree_0.4_weighted[[i]] <- sum(sample_weighted)
}
degree_0.4_weighted <-  as.numeric(unlist(degree_0.4_weighted, use.names = FALSE))
hist(degree_0.4_weighted)

# random sampling without density function - checking the sum of degree for all samples
degree_0.4 <- list()
for(i in 1:10000) {      
  sample <- sample(all_bmsc_ligands_ex$degree_0.4, length(mm_bmsc_ligands_ex$degree_0.4), replace=TRUE)
  
  degree_0.4[[i]] <- sum(sample)
}
degree_0.4 <- as.numeric(unlist(degree_0.4, use.names = FALSE))
hist(degree_0.4)

# random sampling without the density function produce samples that have lower degree_0.4 sum ~ 2200 compared to MM-BMSC degree_0.4 3000

# repeat the same procedure but applying it on dataframe
degree_0.4_weighted <- list()
for(i in 1:10000) {      
  sample_weighted <- all_bmsc_ligands_ex[sample(nrow(all_bmsc_ligands_ex), 104, replace = TRUE, prob=bmsc.dens.obs_0.4$y), ]
  
  degree_0.4_weighted[[i]] <- sum(sample_weighted$degree_0.4)
}
degree_0.4_weighted <-  as.numeric(unlist(degree_0.4_weighted, use.names = FALSE))
hist(degree_0.4_weighted)


degree_0.4<- list()
for(i in 1:10000) {      
  sample <- all_bmsc_ligands_ex[sample(nrow(all_bmsc_ligands_ex), 104, replace = TRUE), ]
  
  degree_0.4[[i]] <- sum(sample$degree_0.4)
}
degree_0.4 <-  as.numeric(unlist(degree_0.4, use.names = FALSE))
hist(degree_0.4)

######################

# load string 0.7 data
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

# density function with adjustment
set.seed(3333)
bmsc.dens.obs_0.7 <- density(mm_bmsc_ligands_ex$degree_0.7, n=length(all_bmsc_ligands_ex$degree_0.7), adjust = 2)
plot(bmsc.dens.obs_0.7$y)

# if the density functions works, the distribution of random sample degree_0.7 sum should be centered around 1200 (MM-BMSC degree_0.7)

# random sampling with density function - checking the sum of degree for all samples
degree_0.7_weighted <- list()
for(i in 1:10000) {      
  sample_weighted <- sample(all_bmsc_ligands_ex$degree_0.7, length(mm_bmsc_ligands_ex$degree_0.7), replace=TRUE, prob=bmsc.dens.obs_0.7$y)
  
  degree_0.7_weighted[[i]] <- sum(sample_weighted)
}
degree_0.7_weighted <-  as.numeric(unlist(degree_0.7_weighted, use.names = FALSE))
hist(degree_0.7_weighted)

# random sampling without density function - checking the sum of degree for all samples
degree_0.7 <- list()
for(i in 1:10000) {      
  sample <- sample(all_bmsc_ligands_ex$degree_0.7, length(mm_bmsc_ligands_ex$degree_0.7), replace=TRUE)
  
  degree_0.7[[i]] <- sum(sample)
}
degree_0.7 <- as.numeric(unlist(degree_0.7, use.names = FALSE))
hist(degree_0.7)

# repeat the same procedure but applying it on dataframe
degree_0.7_weighted <- list()
for(i in 1:10000) {      
  sample_weighted <- all_bmsc_ligands_ex[sample(nrow(all_bmsc_ligands_ex), 104, replace = TRUE, prob=bmsc.dens.obs_0.7$y), ]
  
  degree_0.7_weighted[[i]] <- sum(sample_weighted$degree_0.7)
}
degree_0.7_weighted <-  as.numeric(unlist(degree_0.7_weighted, use.names = FALSE))
hist(degree_0.7_weighted)


degree_0.7<- list()
for(i in 1:10000) {      
  sample <- all_bmsc_ligands_ex[sample(nrow(all_bmsc_ligands_ex), 104, replace = TRUE), ]
  
  degree_0.7[[i]] <- sum(sample$degree_0.7)
}
degree_0.7 <-  as.numeric(unlist(degree_0.7, use.names = FALSE))
hist(degree_0.7)

# random sampling without the density function produce samples that have lower degree_0.7 ~ 800 compared to MM-BMSC degree_0.7 1200

#############################

# load string 0.9 data
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

# density function with adjustment
set.seed(3333)
bmsc.dens.obs_0.9 <- density(mm_bmsc_ligands_ex$degree_0.9, n=length(all_bmsc_ligands_ex$degree_0.9), adjust = 3.5)
plot(bmsc.dens.obs_0.9$y)

# if the density functions works, the distribution of degree_0.9 should be centered around 700 (MM-BMSC degree_0.9)

# random sampling with density function - checking the sum of degree for all samples
degree_0.9_weighted <- list()
for(i in 1:10000) {      
  sample_weighted <- sample(all_bmsc_ligands_ex$degree_0.9, length(mm_bmsc_ligands_ex$degree_0.9), replace=TRUE, prob=bmsc.dens.obs_0.9$y)
  
  degree_0.9_weighted[[i]] <- sum(sample_weighted)
}
degree_0.9_weighted <-  as.numeric(unlist(degree_0.9_weighted, use.names = FALSE))
hist(degree_0.9_weighted)

# random sampling without density function - checking the sum of degree for all samples
degree_0.9 <- list()
for(i in 1:10000) {      
  sample <- sample(all_bmsc_ligands_ex$degree_0.9, length(mm_bmsc_ligands_ex$degree_0.9), replace=TRUE)
  
  degree_0.9[[i]] <- sum(sample)
}
degree_0.9 <- as.numeric(unlist(degree_0.9, use.names = FALSE))
hist(degree_0.9)

# repeat the same procedure but applying it on dataframe
degree_0.9_weighted <- list()
for(i in 1:10000) {      
  sample_weighted <- all_bmsc_ligands_ex[sample(nrow(all_bmsc_ligands_ex), 104, replace = TRUE, prob=bmsc.dens.obs_0.9$y), ]
  
  degree_0.9_weighted[[i]] <- sum(sample_weighted$degree_0.9)
}
degree_0.9_weighted <-  as.numeric(unlist(degree_0.9_weighted, use.names = FALSE))
hist(degree_0.9_weighted)


degree_0.9<- list()
for(i in 1:10000) {      
  sample <- all_bmsc_ligands_ex[sample(nrow(all_bmsc_ligands_ex), 104, replace = TRUE), ]
  
  degree_0.9[[i]] <- sum(sample$degree_0.9)
}
degree_0.9 <-  as.numeric(unlist(degree_0.9, use.names = FALSE))
hist(degree_0.9)

# random sampling without the density functins produce samples that have lower degree_0.9 ~ 500 compared to MM-BMSC degree_0.9 700

### pc ligands

# load string 0.4 data
string_human_physical_uniprot_0.4 <- string_human_physical_uniprot[string_human_physical_uniprot$score >= 400,]

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

# density function with adjustment
set.seed(3333)
pc.dens.obs_0.4 <- density(mm_pc_ligands_ex$degree_0.4, n=length(all_pc_ligands_ex$degree_0.4), adjust = 0.5)
plot(pc.dens.obs_0.4$y)

# if the density functions works, the distribution of degree_0.4 should be centered around 3927 (MM-pc degree_0.4)

# random sampling with density function - checking the sum of degree for all samples
degree_0.4_weighted <- list()
for(i in 1:10000) {      
  sample_weighted <- sample(all_pc_ligands_ex$degree_0.4, length(mm_pc_ligands_ex$degree_0.4), replace=TRUE, prob=pc.dens.obs_0.4$y)
  
  degree_0.4_weighted[[i]] <- sum(sample_weighted)
}
degree_0.4_weighted <-  as.numeric(unlist(degree_0.4_weighted, use.names = FALSE))
hist(degree_0.4_weighted)

# random sampling without density function - checking the sum of degree for all samples
degree_0.4 <- list()
for(i in 1:10000) {      
  sample <- sample(all_pc_ligands_ex$degree_0.4, length(mm_pc_ligands_ex$degree_0.4), replace=TRUE)
  
  degree_0.4[[i]] <- sum(sample)
}
degree_0.4 <- as.numeric(unlist(degree_0.4, use.names = FALSE))
hist(degree_0.4)

# random sampling with and without the density function produce samples that have similar degree_0.4 ~ 4500 compared to MM-pc degree_0.4 3900

# repeat the same procedure but applying it on dataframe
degree_0.4_weighted <- list()
for(i in 1:10000) {      
  sample_weighted <- all_pc_ligands_ex[sample(nrow(all_pc_ligands_ex), 178, replace = TRUE, prob=pc.dens.obs_0.4$y), ]
  
  degree_0.4_weighted[[i]] <- sum(sample_weighted$degree_0.4)
}
degree_0.4_weighted <-  as.numeric(unlist(degree_0.4_weighted, use.names = FALSE))
hist(degree_0.4_weighted)


degree_0.4<- list()
for(i in 1:10000) {      
  sample <- all_pc_ligands_ex[sample(nrow(all_pc_ligands_ex), 178, replace = TRUE), ]
  
  degree_0.4[[i]] <- sum(sample$degree_0.4)
}
degree_0.4 <-  as.numeric(unlist(degree_0.4, use.names = FALSE))
hist(degree_0.4)

######################

# load string 0.7 data
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

# density function without adjustment
set.seed(3333)
pc.dens.obs_0.7 <- density(mm_pc_ligands_ex$degree_0.7, n=length(all_pc_ligands_ex$degree_0.7))
plot(pc.dens.obs_0.7$y)

# if the density functions works, the distribution of degree_0.7 should be centered around 1574 (MM-pc degree_0.7)

# random sampling without density function - checking the sum of degree for all samples
degree_0.7_weighted <- list()
for(i in 1:10000) {      
  sample_weighted <- sample(all_pc_ligands_ex$degree_0.7, length(mm_pc_ligands_ex$degree_0.7), replace=TRUE, prob=pc.dens.obs_0.7$y)
  
  degree_0.7_weighted[[i]] <- sum(sample_weighted)
}
degree_0.7_weighted <-  as.numeric(unlist(degree_0.7_weighted, use.names = FALSE))
hist(degree_0.7_weighted)

# random sampling without density function - checking the sum of degree for all samples
degree_0.7 <- list()
for(i in 1:10000) {      
  sample <- sample(all_pc_ligands_ex$degree_0.7, length(mm_pc_ligands_ex$degree_0.7), replace=TRUE)
  
  degree_0.7[[i]] <- sum(sample)
}
degree_0.7 <- as.numeric(unlist(degree_0.7, use.names = FALSE))
hist(degree_0.7)

# random sampling with and without the density function produce samples that have similar degree_0.7 ~ 1500 compared to MM-pc degree_0.7 1574

# repeat the same procedure but applying it on dataframe
degree_0.7_weighted <- list()
for(i in 1:10000) {      
  sample_weighted <- all_pc_ligands_ex[sample(nrow(all_pc_ligands_ex), 178, replace = TRUE, prob=pc.dens.obs_0.7$y), ]
  
  degree_0.7_weighted[[i]] <- sum(sample_weighted$degree_0.7)
}
degree_0.7_weighted <-  as.numeric(unlist(degree_0.7_weighted, use.names = FALSE))
hist(degree_0.7_weighted)


degree_0.7<- list()
for(i in 1:10000) {      
  sample <- all_pc_ligands_ex[sample(nrow(all_pc_ligands_ex), 178, replace = TRUE), ]
  
  degree_0.7[[i]] <- sum(sample$degree_0.7)
}
degree_0.7 <-  as.numeric(unlist(degree_0.7, use.names = FALSE))
hist(degree_0.7)


#############################

# load string 0.9 data
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

# density function with adjustment
set.seed(3333)
pc.dens.obs_0.9 <- density(mm_pc_ligands_ex$degree_0.9, n=length(all_pc_ligands_ex$degree_0.9), adjust = 1.5)
plot(pc.dens.obs_0.9$y)

# if the density functions works, the distribution of degree_0.9 should be centered around 953 (MM-pc degree_0.9)

# random sampling with density function - checking the sum of degree for all samples
degree_0.9_weighted <- list()
for(i in 1:10000) {      
  sample_weighted <- sample(all_pc_ligands_ex$degree_0.9, length(mm_pc_ligands_ex$degree_0.9), replace=TRUE, prob=pc.dens.obs_0.9$y)
  
  degree_0.9_weighted[[i]] <- sum(sample_weighted)
}
degree_0.9_weighted <-  as.numeric(unlist(degree_0.9_weighted, use.names = FALSE))
hist(degree_0.9_weighted)

# random sampling without density function - checking the sum of degree for all samples
degree_0.9 <- list()
for(i in 1:10000) {      
  sample <- sample(all_pc_ligands_ex$degree_0.9, length(mm_pc_ligands_ex$degree_0.9), replace=TRUE)
  
  degree_0.9[[i]] <- sum(sample)
}
degree_0.9 <- as.numeric(unlist(degree_0.9, use.names = FALSE))
hist(degree_0.9)

# random sampling with and without the density function produce samples that have similar degree_0.9 ~ 900 compared to MM-pc degree_0.9 953

# repeat the same procedure but applying it on dataframe
degree_0.9_weighted <- list()
for(i in 1:10000) {      
  sample_weighted <- all_pc_ligands_ex[sample(nrow(all_pc_ligands_ex), 178, replace = TRUE, prob=pc.dens.obs_0.9$y), ]
  
  degree_0.9_weighted[[i]] <- sum(sample_weighted$degree_0.9)
}
degree_0.9_weighted <-  as.numeric(unlist(degree_0.9_weighted, use.names = FALSE))
hist(degree_0.9_weighted)


degree_0.9<- list()
for(i in 1:10000) {      
  sample <- all_pc_ligands_ex[sample(nrow(all_pc_ligands_ex), 178, replace = TRUE), ]
  
  degree_0.9[[i]] <- sum(sample$degree_0.9)
}
degree_0.9 <-  as.numeric(unlist(degree_0.9, use.names = FALSE))
hist(degree_0.9)


## bmsc receptors

# load string 0.4 data
string_human_physical_uniprot_0.4 <- string_human_physical_uniprot[string_human_physical_uniprot$score >= 400,]

# create a subset dataframe that contains sensory neuron annotated as receptors and have at least 1 extracellular resource
all_sensory_neuron_receptors_ex <- proteins_sensory_neurons_annot_transc[proteins_sensory_neurons_annot_transc$receptor == TRUE & proteins_sensory_neurons_annot_transc$num_sources_extracell > 0,]

# get degree_0.4 (no. of interactions) distribution of all receptors 
all_sensory_neuron_receptors_ex$degree_0.4 <- NA
for(i in 1:nrow(all_sensory_neuron_receptors_ex)) {      
  protein <- all_sensory_neuron_receptors_ex[i, ]
  protein_degree_0.4 <-subset(string_human_physical_uniprot_0.4, 
                              string_human_physical_uniprot_0.4$protein1 %in% 
                                protein$uniprot_id)
  
  all_sensory_neuron_receptors_ex$degree_0.4[[i]] <- nrow(protein_degree_0.4)
  
}

## observed interactions for DEG from MM-BMSC 0.4
mm_bmsc_104_receptor_interactions_0.4 <- subset(string_human_physical_uniprot_0.4, 
                                                string_human_physical_uniprot_0.4$protein1 %in% 
                                                  merged_analysis_BMSC_5DS_ID_proteome_annot[merged_analysis_BMSC_5DS_ID_proteome_annot$ligand == TRUE & merged_analysis_BMSC_5DS_ID_proteome_annot$num_sources_extracell > 0, ]$uniprotswissprot 
                                                & string_human_physical_uniprot_0.4$protein2 %in% 
                                                  proteins_sensory_neurons_annot_transc[proteins_sensory_neurons_annot_transc$receptor == TRUE & proteins_sensory_neurons_annot_transc$num_sources_extracell > 0 ,]$uniprot_id)
nrow(mm_bmsc_104_receptor_interactions_0.4)
#[1] 678

# create a subset dataframe of receptors that DEG from MM-BMSC interact with at 0.4
mm_bmsc_104_receptor_0.4 <- subset(all_sensory_neuron_receptors_ex, 
                                   all_sensory_neuron_receptors_ex$uniprot_id%in% 
                                     mm_bmsc_104_receptor_interactions_0.4$protein2) 


nrow(mm_bmsc_104_receptor_0.4)
#[1] 331

# plot the distribution
mm_bmsc_receptor_degree_0.4 <-  as.numeric(unlist(mm_bmsc_104_receptor_0.4$degree_0.4, use.names = FALSE))
hist(mm_bmsc_receptor_degree_0.4)

# sum of MM-BMSC degree_0.4
sum(mm_bmsc_receptor_degree_0.4)
#[1] 17899


# plot the distribution
all_receptor_degree_0.4 <- as.numeric(unlist(all_sensory_neuron_receptors_ex$degree_0.4, use.names = FALSE))
hist(all_receptor_degree_0.4)

# sum of all BMSC degree_0.4
sum(all_receptor_degree_0.4)
#[1] 25715

# density function 
set.seed(3333)
bmsc.receptor.dens.obs_0.4 <- density(mm_bmsc_104_receptor_0.4$degree_0.4, n=length(all_sensory_neuron_receptors_ex$degree_0.4), adjust = 0.5)
plot(bmsc.receptor.dens.obs_0.4$y)

# if the density functions works, the distribution of random samples degree_0.4 sum should be centered around 17899 (MM-BMSC degree_0.4)

# random sampling with density function - checking the sum of degree for all samples
receptor_degree_0.4_weighted <- list()
for(i in 1:10000) {      
  sample_weighted <- sample(all_sensory_neuron_receptors_ex$degree_0.4, length(mm_bmsc_104_receptor_0.4$degree_0.4), replace=TRUE, prob=bmsc.receptor.dens.obs_0.4$y)
  
  receptor_degree_0.4_weighted[[i]] <- sum(sample_weighted)
}
receptor_degree_0.4_weighted <-  as.numeric(unlist(receptor_degree_0.4_weighted, use.names = FALSE))
hist(receptor_degree_0.4_weighted)

# random sampling without density function - checking the sum of degree for all samples
degree_0.4 <- list()
for(i in 1:10000) {      
  sample <- sample(all_sensory_neuron_receptors_ex$degree_0.4, length(mm_bmsc_104_receptor_0.4$degree_0.4), replace=TRUE)

  degree_0.4[[i]] <- sum(sample)
}
degree_0.4 <- as.numeric(unlist(degree_0.4, use.names = FALSE))
hist(degree_0.4)

# random sampling without the density function produce samples that have lower degree_0.4 sum ~ 2200 compared to MM-BMSC degree_0.4 3000



