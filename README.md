## REMP: Repetitive Element Methylation Prediction
[![GitHub release](https://img.shields.io/badge/release-v1.2.5-blue.svg)](https://www.bioconductor.org/packages/release/bioc/html/REMP.html)

*REMP* provides machine learning-based tools to predict DNA methylation of locus-specific repetitive elements (RE) by learning surrounding genetic and epigenetic information. These tools provide genomewide and single-base resolution of DNA methylation prediction on RE that are difficult to measure using array-based or sequencing-based platforms, which enables epigenome-wide association study (EWAS) and differentially methylated region (DMR) analysis on RE. 

## Installation 

REMP is available in Bioconductor repository
```r

## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")

## To get the latest REMP, please first make sure you are using the latest Bioconductor
biocLite("BiocUpgrade")

## Install REMP
biocLite("REMP")

```

## Citation:

Zheng Y, Joyce BT, Liu L, Zhang Z, Kibbe WA, Zhang W, Hou L. Prediction of genome-wide DNA methylation in repetitive elements. Nucleic Acids Research. 2017;45(15):8697-711. [PubMed PMID: 28911103](https://www.ncbi.nlm.nih.gov/pubmed/28911103); [PMCID: PMC5587781](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5587781/). [http://dx.doi.org/10.1093/nar/gkx587](http://dx.doi.org/10.1093/nar/gkx587).

## Contact package maintainer:
Yinan Zheng 

Email: y-zheng@northwestern.edu
