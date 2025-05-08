#Transcriptomic analysis

This folder contains the code used to analyse transcriptomic data of multiple myeloma (MM) cell types and senosry neurons. 

## Individual datasets

### MM (Multiple Myeloma cell types)

This folder contains folders of individual datasets of MM cell types (ADCY, BMSC, EC, HSC, MP, NEUT, NKC, OCY, OPC, PC, Treg, WBM). Depending on the cell type, there are datasets from human, mouse, or both organisms.  
For each dataset, there is a folder with the R script that was used to process the data and a `data` directory. The `data` directory does not contain the actual data from GEO (due to size restrictions on github), but a file with a link from where the data can be downloaded. For microarray datasets, a `target.csv` file is also included in the `data` folder that contains meta information used during the analysis. For RNA-seq datasets, the `data` folder contains a `design.csv` file, which contains the meta information.  

The follwoing datasets were processed and used in the analysis: 
- ADCY: GSE143269 (mouse)
- BMSC: GSE36474, GSE46053_Coculture, GSE78235, GSE80608, GSE87073
- EC: GSE14230 (human)
- HSC: GSE24870  (human)
- MP: GSE143025, GSE176385 (mouse)
- NEUT: GSE150021 (human)
- NKC: GSE27838 (human)
- OCY: GSE27372 (human)
- OPC: GSE87073 (human)
- PC: 
	- GSE116294, GSE125361, GSE153380, GSE175384, GSE39754, GSE47552, GSE6477, GSE76425 (human)
	- GSE111921, GSE153626 (mouse)
- Treg: GSE109533 (mouse)
- WBM: GSE118985 (human)

The results files for each dataset can be found in the OSF repository: 

### Sensory neurons

This folder contains individual datasets of senosry neurons.

## Meta-analysis


