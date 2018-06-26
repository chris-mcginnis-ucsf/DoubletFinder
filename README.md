# DoubletFinder
DoubletFinder is an R package that predicts doublets in single-cell RNA sequencing data. 

DoubletFinder is implemented to interface with Seurat (https://satijalab.org/seurat/).

For more information, check out our preprint (https://www.biorxiv.org/content/early/2018/06/20/352484)

## Installation
In R/RStudio:

devtools::install_github('chris-mcginnis-ucsf/doubletFinder')

## Dependencies
DoubletFinder requires the Seurat (>= 2.0) and Matrix R packages

## Usage

DoubletFinder can be broken up into 4 steps:

(1) Generate artificial doublets from existing scRNA-seq data 

(2) Perform PCA on the merged real-artificial data

(3) Calculate the euclidean distance matrix in PC space and find each cell's proportion of artificial nearest neighbors (pANN)

(4) Rank order and threshold pANN values according to the expected number of doublets

![Alt text](Users/cmcginnis/Desktop/DDFig1B?raw=true "DoubletFinder Worklfow")

DoubletFinder takes as an input a fully-processed Seurat object (i.e., after NormalizeData, FindVariableGenes, ScaleData, RunPCA, and RunTSNE ahve all been run).  

