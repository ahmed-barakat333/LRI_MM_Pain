# LRI_MM_Pain

This is a repo to reproduce results and figures from the manuscript "A ligand-receptor interactome of the bone tumor microenvironment in multiple myeloma bone pain".

The data and code are organized into folders based on the analysis type. Folders have a common structure (.R filescript, a “data” folder, and a “results”  folder) so the ".R" file should read data files from the "data" folder and output analysis results in the "results" folder.

The contents of this repository include:
---------------------------------------

1- **Transcriptomic_analysis**: Data and code used to analyse transcriptomic data of MM cell types and senosry neurons. This folder has the following subfolders:

 - Individual_datasets: Data and code used to analyse individual datasets of MM cell types and senosry neurons. This folder is subdivided into the following subfolders:
     - MM. This folder contains folders of individual datasets of MM cell types (ADCY, BMSC, EC, HSC, MP, NEUT, NKC, OCY, OPC, PC, Treg, WBM).
     - Sensory_neuron. This folder contains individual datasets of senosry neurons.
 - Meta-analysis: Data and code used to perform meta-analysis of MM cell types and senosry neurons. This folder is subdivided into the following subfolders:
     - MM. This folder contains folders of meta-analysis of individual datasets of some MM cell types (human_BMSC, human_PC, mouse_MP, mouse_PC).
     - Sensory_neuron. This folder contains folders of meta-analysis of individual datasets of senosry neurons.
  
 2- **Signature_analysis**: Data and code used to analyse the gene signature of MM cell types. Signature analysis is based on Rank Rank Hypergeometric Overlap (RRHO) and pathway enrichment (fgsea). Overlap of enriched pathways, transcription factors, and kinases were assessed between MM cell types were assessed using Fisher's exact test. 
 


