Data and code used to analyse the gene signature of MM cell types. 
Signature analysis is based on:

 -   [Rank Rank Hypergeometric Overlap](https://academic.oup.com/nar/article/38/17/e169/1033168) (RRHO): a threshold-free, statistical method to compare two ranked gene lists to identify regions of significant overlap.

 -  [Pathway enrichment analysis](https://www.biorxiv.org/content/10.1101/060012v3.abstract) (fgsea): gene set enrichment analysis (GSEA) to determine whether predefined gene sets show statistically significant, concordant differences between MM vs control.

-  Assessment of the overlap of enriched pathways, transcription factors, and kinases between MM cell types using Fisher's exact test. This statistical test evaluates whether the observed overlap between two sets of enriched features (e.g., pathways shared between two cell types) is greater than expected by chance, providing a measure of similarity or shared biology between cell types.

   - Enriched pathways were retrieved from fgsea analysis. Enriched transcription factors and kinases were retreved from [X2K](https://maayanlab.cloud/X2K/) online tool. 

