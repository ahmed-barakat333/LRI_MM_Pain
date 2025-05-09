# Transcriptomic analysis

This folder contains the code used to analyse transcriptomic data of multiple myeloma (MM) cell types and sensory neurons. 

## Individual datasets

### MM (Multiple Myeloma cell types)

This folder contains folders of individual datasets of MM cell types (ADCY, BMSC, EC, HSC, MP, NEUT, NKC, OCY, OPC, PC, Treg, WBM). Depending on the cell type, there are datasets from human, mouse, or both organisms.  
For each dataset, there is a folder with the R script that was used to process the data and a `data` directory. The `data` directory does not contain the actual data from GEO (due to size restrictions on github), but a file with a link from where the data can be downloaded. For microarray datasets, a `target.csv` file is also included in the `data` folder that contains meta information used during the analysis. For RNA-seq datasets, the `data` folder contains a `design.csv` file, which contains the meta information. For the cell types with more than 1 dataset available, the result files used for meta analysis can also be found in `/Meta-analysis/MM/`. All results files are available in the [OSF repository](https://osf.io/3jwys/files/osfstorage#) under `Processed data`.  

The following datasets were processed and used in the analysis: 
- ADCY: [GSE143269](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE143269) (mouse)
- BMSC: [GSE36474](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE36474), [GSE46053_Coculture](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE46053), [GSE78235](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE78235), [GSE80608](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80608), [GSE87073](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87073) (human)
- EC: [GSE14230](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE14230) (human)
- HSC: [GSE24870](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE24870
)  (human)
- MP: [GSE143025](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE143025), [GSE176385](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176385) (mouse)
- NEUT: [GSE150021](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150021) (human)
- NKC: [GSE27838](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE27838) (human)
- OCY: [GSE27372](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE27372) (human)
- OPC: [GSE87073](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87073) (human)
- PC: 
	- [GSE116294](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE116294), [GSE125361](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125361), [GSE153380](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153380), [GSE175384](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE175384), [GSE39754](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE39754), [GSE47552](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47552), [GSE6477](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE6477), [GSE76425](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76425) (human)
	- [GSE111921](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111921), [GSE153626](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153626) (mouse)
- Treg: [GSE109533](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109533) (mouse)
- WBM: [GSE118985](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118985) (human)

### Sensory_neuron

This folder contains individual datasets of sensory neurons divided by organism (mouse or human). For each dataset, there is a folder with the R script that was used to process the data and a `data` directory. The `data` directory does not contain the actual data from GEO (due to size restrictions on github), but a file with a link from where the data can be downloaded. The `results` folder contains the processed file. 

The following datasets were processed and used in the analysis: 
- Human: [GSE168243](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE168243)
- Mouse: [GSE100175](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100175), [GSE131230](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131230), [GSE139088](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139088), [GSE154659](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154659), [GSE155622](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE155622), [GSE168032](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE168032), [GSE59739](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE59739), [GSE62424](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62424), [GSE63576](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63576), [PMID30096314](http://mousebrain.org/)

## Meta-analysis

### MM 

This folder contains the meta-analysis for BMSC, human PC, mouse PC, and MP cells. Each of these has a dedicated folder with an R script, `data` folder with all needed input files (coming from the transcriptomics analysis of the individual datasets), and a `results` folder with the file generated by MetaIntegrator.  

### Sensory_neuron

This folder contains the meta-analysis of the sensory neuron data for mouse (`mouse_Meta-analysis`) as well as a comparison of the mouse meta-analysis with the human dataset (`human_mouse_overlap`). Each folder contains the respective R script, `data` folder, and a `results` folder.   

