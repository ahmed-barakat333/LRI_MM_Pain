# LRI_MM_Pain

This is a repo to reproduce results and figures from the manuscript "A ligand-receptor interactome of the bone tumor microenvironment in multiple myeloma bone pain".

The data and code are organized into folders based on the analysis type. All folders have a common structure (.R file, "data" folder, and "results" folder) so the ".R" file should read data files from the "data" folder and output analysis results in the "results" folder.

> Please note that due to size constrains, different subfolders located inside the main folders are compressed to save space.

The contents of this repository include:
---------------------------------------

1- **CARNIVAL_Network_Analysis_BMSCs**: Data and code used to construct signaling network of multiple myeloma - bone marrow stromal cells (BMSCs). Please note, to run this code, a Gurobi solver is needed (https://www.gurobi.com/solutions/gurobi-optimizer/). 

2- **CARNIVAL_Network_Analysis_PCs**: Data and code used to construct signaling network of multiple myeloma - plasma cells (PCs). Please note, to run this code, a Gurobi solver is needed (https://www.gurobi.com/solutions/gurobi-optimizer/). 

3- **Cells_Ligands**: Data and code used to visualize ligands of each cell type of multiple myeloma bone marrow, in addition to annotate each ligand with "biological process" and "molecular function" keywords. 

4- **Chord_Diagram**: Data and code used to visualize some of the interactions between multiple myeloma - plasma cells (PCs) and bone marrow stromal cells (BMSCs) ligands with seonsory neuron receptors. 

5- **Communication_Score**: Data and code used to the significance of interactions between multiple myeloma - plasma cells (PCs) and bone marrow stromal cells (BMSCs) ligands with seonsory neuron receptors.

6- **Functional_Analysis**: Data and code used to perform enrichment analysis of differentially expressed genes of mutiple myeloma cell types. 

7- **Overlap_Analysis**: Data and code used to assess and visualize the overlap between activated pathways in multiple myeloma cell types. 

8- **RRHO_Analysis**: Data and code used to run Rank Rank Hypergeometric Overlap (RRHO) analysis of differentially expressed genes of multiple myeloma cell types. 

9- **Transcriptomics_Analysis_MM**: Data and code used to analyse multiple myeloma cell-types transcriptomic (microarray, bulk RNA-seq, meta-analysis) datasets. This folder is organized into subfolders based on cell type:

 - ADCY: Data and code used to analyse adipocytes data
 - BMSC: Data and code used to analyse bone marrow stromal cells data
 - EC: Data and code used to analyse endothelial cells data
 - HSC: Data and code used to analyse hematopoietic stem cell data
 - MP: Data and code used to analyse macrophages data
 - NEUT: Data and code used to analyse neutrophils data
 - NKC: Data and code used to analyse natural killer cells data
 - OCY: Data and code used to analyse osteocytes data
 - OPC: Data and code used to analyse osteogenic precursor cells data
 - PC: Data and code used to analyse plasma cells data
 - Treg: Data and code used to analyse regulatory T cells data
 - WBM: Data and code used to analyse whole bone marrow data

 > Please note that due to size and other constraints, the raw data files (.CEL files and raw count matrices) could not be uploaded.
 > A link for downloading the raw data files has been provided. This link is provided inside a "data_link.txt" file located inside the folder of each dataset. 

10- **Transcriptomics_Analysis_Sensory_Neuron**: Data and code used to analyse senosry neuron transcriptomic datasets. This folder is orgnaized into the following subfolders:

  - human_Datasets: Data and code to analyse the human dataset (GSE168243).
  - human_mouse_overlap:  Data and code to get genes that overlap in their expression between mouse and human datasets. 
  - mouse_Datasets: Data and code to analyse different mouse datasets.
      > Please note that due to size and other constraints, the raw data files (raw count matrices) could not be uploaded.
      > A link for downloading the raw data files has been provided. This link is provided inside a "data_link.txt" file located inside the folder of each dataset. 
  - mouse_Meta-analysis: Data and code to run rank-aggregaton meta-analysis of different mouse datasets.


