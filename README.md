# DoubletFinder
DoubletFinder is an R package that predicts doublets in single-cell RNA sequencing data. 

DoubletFinder is implemented to interface with Seurat >= 2.0 (https://satijalab.org/seurat/) 

For more information, check out our preprint (https://www.biorxiv.org/content/early/2018/06/20/352484)

## Installation (in R/RStudio)

devtools::install_github('chris-mcginnis-ucsf/doubletFinder')

## Dependencies
DoubletFinder requires the Seurat (>= 2.0) and Matrix R packages

## Usage

DoubletFinder can be broken up into 4 steps:

(1) Generate artificial doublets from existing scRNA-seq data 

(2) Perform PCA on the merged real-artificial data

(3) Calculate the euclidean distance matrix in PC space and find each cell's proportion of artificial nearest neighbors (pANN)

(4) Rank order and threshold pANN values according to the expected number of doublets

![alternativetext](DF.screenshots/Workflow2.png)

DoubletFinder takes the following arguments:

seu ~ This is a fully-processed Seurat object (i.e., after NormalizeData, FindVariableGenes, ScaleData, RunPCA, and RunTSNE ahve all been run).

proportion.artificial ~ This defines the number of generated artificial doublets, expressed as a proportion of the merged real-artificial data. Default is set to 25%, based on previous optimization (see McGinnis, Murrow and Gartner 2018, BioRxiv).

proportion.NN ~ This defines the PC neighborhood size sued to calculate pANN, expressed as a proportion of the merged real-artificial data. Default is set to 1%, based on previous optimization (see McGinnis, Murrow and Gartner 2018, BioRxiv).

expected.doublets ~ This defines the pANN threshold used to make final doublet/singlet predictions. This value can best be estimated from cell loading densities into the 10X/Drop-Seq device.

## Application to Cell Hashing and Demuxlet data

![alternativetext](DF.screenshots/Results2.png)

## References

1.	Stoeckius M, Zheng S, Houck-Loomis B, Hao S, Yeung BZ, Smibert P, Satija R. Cell "hashing" with barcoded antibodies enables multiplexing and doublet detection for single cell genomics. 2017. Preprint. bioRxiv doi: 10.1101/237693.

2.  Kang HM, Subramaniam M, Targ S, Nguyen M, Maliskova L, McCarthy E. Multiplexed droplet single-cell RNA-sequencing using natural genetic variation. Nat Biotechnol. 2018; 36(1):89-94. 


## SessionInfo

R version 3.4.1 (2017-06-30)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS High Sierra 10.13.2

Matrix products: default

BLAS: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRblas.0.dylib

LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats     [2] graphics  [3] grDevices utils     
[4] datasets  [5] methods   [6] base     

other attached packages:
[1] doubletFinder_1.0.0  [2] devtools_1.13.5   [3] cowplot_0.9.2       
[4] ggplot2_2.2.1        [5] Matrix_1.2-13      




