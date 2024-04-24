The repository contains the code for the computational studies in [The fibro-adipogenic progenitor APOD+DCN+LUM+ cell population in aggressive carcinomas](https://doi.org/10.1007/s10555-024-10181-y).

# Adaptive Attractor
An adaptive method for finding co-expression signatures in single-cell/bulk RNA-seq datasets. 

- [Overview](#Overview)
- [Tutorials](#Tutorials)
  - [Quick start](#Quickstart)
  - [Examples](#Examples)
  - [Parameters](#Parameters)
- [Seed selection](#Seed-selection)
  - [Demonstration](#Demonstration)
- [Description of original attractor algorithm (fixed exponent parameter)](#Description-of-original-attractor-algorithm)


## Overview
An adaptive version of attractor algorithm. 

`cafr::findAttractor` finds a converged attractor based on the seed gene provided. The `findAttractor.adaptive` gradually decreased the exponent parameter to maximize the strength of the Nth-ranked genes in the converged attractor.


The attractor algorithm was first proposed for identifying co-expression signatures from bulk expression values in samples [[Ref.1](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002920)]. A detailed description and application of the attractor algorithm on single-cell RNA-seq data can be found in [[Ref.2](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009228#pcbi.1009228.ref020)]. Users can refer to our recent manuscript [[Ref.3](https://doi.org/10.1093/bioinformatics/btae283)] for a detailed description and application of adaptive attractor algorithm.

## Tutorials
### Quick start
Please see [Seed selection](#Seed-selection) for more information about the seed. <br /> 
```R
source("./findAttractor.adaptive.r")
data <- as.matrix(data)  # a normalized expression matrix with genes in the rows, cells in the columns.
seed <- "LUM"
attr <- findAttractor.adaptive(data, seed)
```

### Examples
R script for generating Supplementary Table 1 in our recent manuscript [link]. 
```R
## install cafr package 
devtools::install_github("weiyi-bitw/cafr")

## load metadata
# Data website: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129363
metadata <- data.table::fread("./GSE129363/raw/GSE129363_Discovery_Cohort_CellAnnotation.txt")

## normalize the count matrix 
library("Seurat")
adata <- Read10X("./GSE129363/raw") #including barcodes.tsv, genes.tsv, matrix.mtx
adata <- CreateSeuratObject(counts = adata[, metadata$CellID])    
adata <- NormalizeData(adata)
data_all <- adata[["RNA"]]@data; rm(adata); gc()

## select one sample
sampleID <- "Ate-05-SAT"
data <- data_all[, metadata[which(metadata$SampleName == sampleID), ]$CellID]
  
## apply adaptive attractor algorithm
source("./findAttractor.adaptive.r")
data <- as.matrix(data)  #normalized matrix for each sample
seed <- "LUM"
attr <- findAttractor.adaptive(data, seed, exponent.max = 5, exponent.min = 2, step.large = 1, step.small = 0.1, dominantThreshold = 0.2, seedTopN = 50, compareNth = 10, verbose = FALSE)

## results
head(attr$attractor.final)
attr$final.exponent

```


### Parameters
`data` An expression matrix with genes in the rows, samples in the columns. <br />
`seed` Seed gene. Please see [Seed selection](#Seed-selection) for more information. <br />
`exponent.max` The maximum exponent of the mutual information, used to create weight vector for metagenes. <br />
`exponent.min` The minimum exponent of the mutual information, used to create weight vector for metagenes. <br />
`step.large` Decreasing step of the first round of scanning. <br />
`step.small` Decreasing step of the second round of scanning. <br />
`dominantThreshold` Threshold of dominance. <br />
`seedTopN` Threshold of divergence. <br />
`compareNth` It will maxmize the weight of the top Nth gene.<br />
`verbose` When TRUE, it will show the top 20 genes of the metagene in each iteration.<br />
`epsilon` Threshold of convergence.<br />
`minimize` Default FALSE,when minimize = TRUE, the algorithm will minize the difference of metagene[first.idx] and metagene[second.idx] instead of maxmizing  the weight of the top Nth gene. <br />
`first.idx` Default first.idx = 1, the first gene of the metagene in each iteration.  <br />
`second.idx` Default second.idx = 2, the second gene of the metagene in each iteration.  <br />

## Seed selection

Users can choose some general markers as seed genes, such as:

| Fibroblast  | Pericyte | Macrophage | Endothelial | Mitotic |
| ----------- | -------- |----------- | ----------- | ------- |
| LUM         | RGS5     | AIF1       | PECAM1      | TOP2A   |
| DCN         |          |            |             |         |   


If the dataset has many cells of one cell type (e.g. fibroblasts) and the attractor exponent (`a`) is fixed, then identical attractors will be found using different general markers for one cell type, like DCN, LUM, and COL1A1. 

### Demonstration
Here, we applied the attractor finding algorithm using the `findAttractor` function implemented in the cafr (v0.312) R package [[Ref.1](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002920)] with the general fibroblastic marker gene LUM and DCN as seed.<br />
![](https://github.com/LingyiC/adaptiveAttractor/blob/main/others/LUM_a3.gif)
![](https://github.com/LingyiC/adaptiveAttractor/blob/main/others/DCN_a3.gif)

## Description of original attractor algorithm
Detailed descriptions can be found in [Ref.1](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002920). 

The exponent (`a`) is fixed. For the analysis of UMI based (e.g. 10x) and full-length-based (e.g. Smart-seq2) datasets, we would suggest using `a = 3` and `a = 5`, respectively.

```R
## install cafr package 
devtools::install_github("weiyi-bitw/cafr")
cafr::findAttractor(data, vec, a = 3, epsilon = 1e-7)
```

