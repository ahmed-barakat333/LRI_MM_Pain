# LRI_MM_Pain

This is a repo to reproduce results and figures from the manuscript "A ligand-receptor interactome of the bone tumor microenvironment in multiple myeloma bone pain".

The data and code are organized into folders:

1- CARNIVAL_Network_Analysis_BMSCs: Data and code used to construct signaling network of multiple myeloma - bone marrow stromal cells (BMSCs). Please note to run this code, a Gurobi solver is needed (https://www.gurobi.com/solutions/gurobi-optimizer/). 

2- Cells_Ligands: Data and code used to visualize ligands of each cell type of multiple myeloma bone marrow, in addition to annotate each ligand with "biological process" and "molecular function" keywords. 

3- Chord_Diagram: Data and code used to visualize some of the interactions between multiple myeloma - plasma cells (PCs) and bone marrow stromal cells (BMSCs) ligands with seonsory neuron receptors. 

4- Communication_Score: Data and code used to the significance of interactions between multiple myeloma - plasma cells (PCs) and bone marrow stromal cells (BMSCs) ligands with seonsory neuron receptors.

5- Functional_Analysis: Data and code used to perform enrichment analysis of differentially expressed genes of mutiple myeloma cell types. 

6- Ligand_Receptor_Annotations: A file containing final results of ligand-receptor pairs of each multple myeloma cell type with senosry neurons. 

7- Overlap_Analysis: Data and code used to assess and visualize the overlap between activated pathways in multple myeloma cell types. 

8- RRHO_Analysis: Data and code used to run Rank Rank Hypergeometric Overlap (RRHO) analysis of differentially expressed genes of multiple myeloma cell types. 

9- Transcriptomics_Analysis_Sensory_Neuron: Data and code used to analyse senosry neuron datasets. This folder is orgnaized into the following subfolders 
  a- human_Datasets: Data and code to analyse the human dataset (GSE168243).
  
  b- human_mouse_overlap:  Data and code to get genes that overlap in their expression between mouse and human datasets. 
  
  c- mouse_Datasets: Data and code to analyse different mouse datasets. Due to size constrains, this folder was uploaded to Google Drive, and the link has been provided. 
  
  d- mouse_Meta-analysis: Data and code to run rank-aggregaton meta-analysis of different mouse datasets.

10- Transcriptomics_Analysis_MM: Data and code used to analyse multiple myeloma cell-types datasets. This folder is organized into subfolders based on cell type:
