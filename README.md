## REMP: Repetitive Element Methylation Prediction

*REMP* provides machine learning-based tools to predict DNA methylation of locus-specific repetitive elements (RE) by learning surrounding genetic and epigenetic information. These tools provide genomewide and single-base resolution of DNA methylation prediction on RE that are difficult to measure directly using array-based or sequencing-based platforms, which enables epigenome-wide association study (EWAS) and differentially methylated region (DMR) analysis on RE. 

Supporting methylation platforms:
+ ** Illumina 450k BeadChip Array
+ ** Illumina EPIC BeadChip Array
+ ** Illumina methylation sequencing platforms (e.g., TruSeq Methyl Capture EPIC)

## Installation 

REMP is available in Bioconductor repository
```r
## try http:// if https:// URLs are not supported
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")

## To get the latest version of REMP, please make sure you are using the latest Bioconductor

## Install REMP
BiocManager::install("REMP")
```

## Citation:

Zheng Y, Joyce BT, Liu L, Zhang Z, Kibbe WA, Zhang W, Hou L. Prediction of genome-wide DNA methylation in repetitive elements. Nucleic Acids Research. 2017;45(15):8697-711. [PubMed PMID: 28911103](https://www.ncbi.nlm.nih.gov/pubmed/28911103); [PMCID: PMC5587781](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5587781/). [http://dx.doi.org/10.1093/nar/gkx587](http://dx.doi.org/10.1093/nar/gkx587).

## Contact package maintainer:
Yinan Zheng 

Email: y-zheng@northwestern.edu
