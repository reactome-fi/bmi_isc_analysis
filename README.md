This repo hosts the code for a collobration project between Dr. Missy Wong's group and Dr. Guanming Wu's group at OHSU for analyzing BMI1+ intestinal stem cells during development, homeostasis and damage repair.

The current version of the manuscript is in the [manuscript folder](https://github.com/reactome-fi/bmi_isc_analysis/tree/main/manuscript).

Raw 10x counts and FASTA files can be found on GEO ([GSE212798](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212798)) and BioProject ([PRJNA877285](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA877285)).

The resources folder contains static files to help with reproducing the results where some variability is always observed like CytoTRACE results and UMAP coordinates. DEG lists and Reactome pathway-gene sets are also present as well as the exported conda environment used for generating the results in this manuscript.

Python and R scripts are in their respective folders and should output results within their directories. 

The release was updated so that a DOI could be generated automatically via Zenodo. The generated DOI is: 10.5281/zenodo.14957091.


