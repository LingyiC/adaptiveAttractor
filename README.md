# Adaptive Attractor
The adaptive version of attractor algorithm. 

The adaptive attractor algorithm gradually decreased the exponent parameter to maximize the strength of the Nth-ranked genes in the converged attractor.

Users can refer to our manuscript [link] for a detailed description and application.

## Example
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
data <- data  #normalized matrix for each sample
seed <- "LUM"
attr <- findAttractor.adaptive(data, seed, exponent.max = 5, exponent.min = 2, step.large = 1, step.small = 0.1, dominantThreshold = 0.2, seedTopN = 50, compareNth = 10, verbose = FALSE)

## results
head(attr$attractor.final)
attr$final.exponent

```

## Parameters
`data` An expression matrix with genes in the rows, samples in the columns. <br />
`seed` Seed gene. <br />
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



