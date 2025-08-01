# LRI_MM_Pain

This is a repo to reproduce results and figures from the manuscript ["A ligand-receptor interactome of the bone tumor microenvironment in multiple myeloma bone pain"](https://journals.lww.com/pain/pages/articleviewer.aspx?year=9900&issue=00000&article=00979&type=Abstract&fbclid=IwZXh0bgNhZW0CMTEAAR66BrmDwcRfD8HMMF4EO69k6ezkz9FHtRBkykMh45LiwGN0_fd_hnGXqzK07A_aem_Rz6LXf7HYHsPPYlGnbhuVA).

The data and code are organized into folders based on the analysis type. Folders have a common structure (.R script, a `data` folder, and a `results`  folder) so the `.R` file should read data files from the `data` folder and output analysis results in the `results` folder.

The contents of this repository include:
---------------------------------------

1- **Transcriptomic_raw_data_analysis**: Data and code used to analyse raw transcriptomic data (microarray and RNA-seq) of multiple myeloma (MM) cell types and senosry neurons. This folder has the following subfolders:

 - `Individual_datasets`: Data and code used to analyse individual datasets of MM cell types and senosry neurons. This folder is subdivided into the following subfolders:
     - `MM`: This folder contains folders of individual datasets of MM cell types (`ADCY`, `BMSC`, `EC`, `HSC`, `MP`, `NEUT`, `NKC`, `OCY`, `OPC`, `PC`, `Treg`, `WBM`).
     - `Sensory_neuron`: This folder contains individual datasets of senosry neurons for mouse and human.
 - `Meta-analysis`: Data and code used to perform meta-analysis of MM cell types and senosry neurons. This folder is subdivided into the following subfolders:
     - `MM`: This folder contains folders of meta-analysis of individual datasets of some MM cell types (`human_BMSC`, `human_PC`, `mouse_MP`, `mouse_PC`).
     - `Sensory_neuron`: This folder contains folders of meta-analysis of mouse datasets of senosry neurons and a comparison with the human dataset.
  
2- **Signature_analysis**: Data and code used to analyse the gene signature of MM cell types. Signature analysis is based on Rank Rank Hypergeometric Overlap (RRHO) and pathway enrichment (fgsea). Overlap of enriched pathways, transcription factors, and kinases were assessed between MM cell types using Fisher's exact test. 
 
3- **Interactome_analysis**: Data and code used to analyse and visualize protein-protein interactions between MM cell types and sensory neurons. This folder has the following subfolders:
 - `Cell_type_ligands`: Data and code used to visualize ligands of each MM cell type, and annotate each ligand with "biological process" and "molecular function" Uniprot keywords.
 - `Chord_diagram`: Data and code used to visualize some of the interactions between ligands of MM - plasma cells (PCs) and bone marrow stromal cells (BMSCs) with sensory neuron receptors.
 - `Communication_score`: Data and code used to analyse the significance of interactions between ligands of MM - PCs and BMSCs with sensory neuron receptors.

4- **Network_analysis**:  Data and code used to construct a signaling network of MM PCs and BMSCs using [CARNIVAL](https://saezlab.github.io/CARNIVAL/) R package.

