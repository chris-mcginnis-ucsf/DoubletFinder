# DoubletFinder

DoubletFinder is an R package that predicts doublets in single-cell RNA sequencing data. 

DoubletFinder is implemented to interface with Seurat >= 2.0 (https://satijalab.org/seurat/) 

For more information, check out our preprint (https://www.biorxiv.org/content/early/2018/06/20/352484)

## DoubletFinder V2.0 (11/28/2018) 

New Features:
1. Increased computational efficiency during pANN computation
2. Implemented strategy for determining optimal pK values for any scRNA-seq data using pN-pK parameter sweeps and mean-variance-normalized bimodality coefficient (BCmvn)
3. Included vignette describing 'best-practices' for applying DoubletFinder to scRNA-seq data generated without sample multiplexing

## Installation (in R/RStudio)

devtools::install_github('chris-mcginnis-ucsf/DoubletFinder')

## Dependencies

DoubletFinder requires the following R packages: 
* Seurat (>= 2.0) 
* Matrix (1.2.14) 
* fields (9.6) 
* KernSmooth (2.23-15)
* modes (0.7.0)
* ROCR (1.0-7)
* NOTE: These package versions were used in the bioRxiv paper, but other versions may work, as well.

# DoubletFinder Overview

DoubletFinder can be broken up into 4 steps:

(1) Generate artificial doublets from existing scRNA-seq data 

(2) Pre-process merged real-artificial data

(3) Perform PCA and use the PC distance matrix to find each cell's proportion of artificial k nearest neighbors (pANN)

(4) Rank order and threshold pANN values according to the expected number of doublets

![alternativetext](DF.screenshots/DoubletFinderOverview.png)

DoubletFinder takes the following arguments:

seu ~ This is a fully-processed Seurat object (i.e., after NormalizeData, FindVariableGenes, ScaleData, RunPCA, and RunTSNE have all been run).

pN ~ This defines the number of generated artificial doublets, expressed as a proportion of the merged real-artificial data. Default is set to 25%, based on observation that DoubletFinder performance is largely pN-invariant (see McGinnis, Murrow and Gartner 2018, BioRxiv).

pK ~ This defines the PC neighborhood size used to compute pANN, expressed as a proportion of the merged real-artificial data. No default is set, as pK should be adjusted for each scRNA-seq dataset. Optimal pK values should be estimated using the strategy described below.

nExp ~ This defines the pANN threshold used to make final doublet/singlet predictions. This value can best be estimated from cell loading densities into the 10X/Drop-Seq device, and adjusted according to the estimated proportion of homotypic doublets.

## Application to Cell Hashing and Demuxlet data

DoubletFinder successfully recapitulates ground-truth doublet classifications determined using antibody-barcode sample multiplexing (Cell Hashing) and SNP deconvolution (Demuxlet). DoubletFinder identifies false-negative Demuxlet classifications caused by doublets formed from cells with identical SNP profiles. DoubletFinder is insensitive to homotypic doublets -- i.e., doublets dervied from transcriptionally-similar cell states. 

![alternativetext](DF.screenshots/Results_Demux.png)
![alternativetext](DF.screenshots/Results_Hashing.png)

# 'Best-Practices' for scRNA-seq data generated without sample multiplexing

## pK Selection

ROC analysis across pN-pK parameter sweeps for Cell Hashing and Demuxlet datasets demonstrate that DoubletFinder performance is largely invariant of pN value selection:

![alternativetext](DF.screenshots/ParamSweep_Schematic.png)
![alternativetext](DF.screenshots/ParamSweep_HeatMap.png)

ROC analysis across pN-pK parameter sweeps for simulated scRNA-seq data with (I) Variable numbers of cell states and (II) Variable magnitudes of transcriptional heterogeneity demonstrates that (I) Optimal pK value selection depends on the total number of cell states and (II) DoubletFinder performance suffers when applied to transcriptionally-homogenous data. Simulated data was generated using a strategy similar to as described in Wolock, Lopex & Klein, 2018.

![alternativetext](DF.screenshots/Simulation_Schematic.png)
![alternativetext](DF.screenshots/Results_Simulation.png)

Simulated and sample-multiplexed data are unique in that ground-truth doublet classifications can be leveraged to characterize how DoubletFinder parameters must be 'fit' to distinct scRNA-seq datasets. However, doublets remain unknown in real-world contexts -- which is likely why you are interested in DoubletFinder, at all!

To maximize the accuracy of DoubletFinder predictions, we sought a ground-truth-agnostic metric that coincides with pK values that maximize AUC in Cell Hashing and Demuxlet data. Mean-variance normalized bimodality coefficient (BCmvn) achieves this goal, featuring a single, easily-discernible maximum at pK values that optimize AUC. 

![alternativetext](DF.screenshots/BCmvn.png)

BCmvn distributions also feature a single maximum for scRNA-seq datasets generated without sample-multiplexing (e.g., Mouse pancreas, Byrnes et al., 2018, Nature Communcations; Mouse kidney, Park et al., 2018, Science), enabling pK selection.

## Doublet Number Estimation

DoubletFinder is sensitive to heterotypic doublets -- i.e., doublets formed from transcriptionally-distinct cell states -- but is insensitive to homotypic doublets -- i.e., doublets formed from transcriptionally-similar cell states. In our original manuscript, we suggested using DoubletFinder to predict the number of doublets expected from Poisson statistical estimates realting to the droplet microfluidics cell loading density. However, Poisson estimates are agnostic of homotypic doublets, and will thus invariably overestimate the number of *detectable* doublets.

To address this issue, we suggest users utilize literature-supported cell type annotations to model the proportion of homotypic doublets present in their data. As an example, we present an analysis of mouse kidney scRNA-seq data (Park et al., 2018, Science):

![alternativetext](DF.screenshots/HomotypicAdjustment.png)

Notably, it is conceivable that literature-suppoted cell type annotations may not accurately recapitulate the magnitude of transcriptional divergence necessary for DoubletFinder sensitivity. For example, nominally-homogenous cells (e.g., CD4+ T-cells) may exist along a spectrum of gene expression states (e.g., distinct anatomical locations, disease states, naive/Tregs/Th17 cells, etc.), and doublets formed by cell sub-types may be detectable by DoubletFinder. Thus, we consider doublet number estimates based on Poisson statistics with and without homotypic doublet proportion adjustment to 'bookend' the real detectable doublet rate. 

## Example code for 'real-world' applications

```R
## Pre-process Seurat object -------------------------------------------------------------------------------------------------
seu_kidney <- CreateSeuratObject(kidney.data)
seu_kidney <- NormalizeData(seu_kidney)
seu_kidney <- ScaleData(seu_kidney, vars.to.regress = "nUMI")
seu_kidney <- FindVariableGenes(seu_kidney, x.low.cutoff = 0.0125, y.cutoff = 0.25, do.plot=FALSE)
seu_kidney <- RunPCA(seu_kidney, pc.genes = seu_kidney@var.genes, pcs.print = 0)
seu_kidney <- RunTSNE(seu_kidney, dims.use = 1:10, verbose=TRUE)

## pK Identification ---------------------------------------------------------------------------------------------------------
sweep.res.list_kidney <- doubletFinder_ParamSweep(seu_kidney)
sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
bcmvn_kidney <- find.pK(sweep.stats_kidney)

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.075*length(seu_kidney@cell.names))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
seu_kidney <- doubletFinder(seu_kidney, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE)
seu_kidney <- doubletFinder(seu_kidney, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913")

## Plot results --------------------------------------------------------------------------------------------------------------
seu_kidney@meta.data[,"DF_hi.lo"] <- seu_kidney@meta.data$DF.classifications_0.25_0.09_913
seu_kidney@meta.data$DF_hi.lo[which(seu_kidney@meta.data$DF_hi.lo == "Doublet" & seu_kidney@meta.data$DF.classifications_0.25_0.09_473 == "Singlet")] <- "Doublet_lo"
seu_kidney@meta.data$DF_hi.lo[which(seu_kidney@meta.data$DF_hi.lo == "Doublet")] <- "Doublet_hi"
TSNEPlot(seu_kidney, group.by="DF_hi.lo", plot.order=c("Doublet_hi","Doublet_lo","Singlet"), colors.use=c("black","gold","red"))
```

![alternativetext](DF.screenshots/DFkidney_low.vs.high.png)

## References

1.	Stoeckius M, Zheng S, Houck-Loomis B, Hao S, Yeung BZ, Smibert P, Satija R. Cell "hashing" with barcoded antibodies enables multiplexing and doublet detection for single cell genomics. 2017. Preprint. bioRxiv doi: 10.1101/237693.

2.  Kang HM, Subramaniam M, Targ S, Nguyen M, Maliskova L, McCarthy E, Wan E, Wong S, Byrnes L, Lanata CM, Gate RE, Mostafavi S, Marson A, Zaitlen N, Criswell LA, Ye JC. Multiplexed droplet single-cell RNA-sequencing using natural genetic variation. Nat Biotechnol. 2018; 36(1):89-94. 

3.  Wolock SL, Lopez R, Klein AM. Scrublet: computational identification of cell doublets in single-cell transcriptomic data. bioRxiv doi: 10.1101/357368.

4.  Park J, Shrestha R, Qiu C, Kondo A, Huang S, Werth M, Li M, Barasch J, Suszták K. Single-cell transcriptomics of the mouse kidney reveals potential cellular targets of kidney disease. Science. 2018; 360(6390):758-63.

5.  Byrnes LE, Wong DM, Subramaniam M, Meyer NP, Gilchrist CL, Knox SM, Tward AD, Ye CJ, Sneddon JB. Lineage dynamics of murine pancreatic development at single-cell resolution. Nat Commun. 2018; 9:3922.

