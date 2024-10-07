# **Comparison of IDH Mutant and IDH Wild Type samples, including the comparison of pan-glioma subtypes**
 
## **Description**
This project is centered around the analysis of Low Grade Glioma (LGG) gene expression. The project is split into three parts:
1. Interpretation of the gene expression dataset using RStudio
2. Visualization of the functional enrichment analysis results using RStudio

You will need the following packages in RStudio:
1. BiocManager: Installs, updates, and manages Bioconductor packages.
2. TCGAbiolinks:nProvides tools for integrating TCGA (The Cancer Genome Atlas) data, facilitating data retrieval and analysis.
3. limma: Linear Models for Microarray and RNA-seq data analysis, focusing on differential expression and gene set testing.
4. edgeR:nEmpirical analysis of digital gene expression data, specifically for differential expression analysis of RNA-seq data.
5. EDASeq: Exploratory Data Analysis for sequencing data, providing methods for quality control, normalization, and visualization.
6. gplots: For R plotting functions, including heatmaps, boxplots, and volcano plots.
7. sesameData:nProvides example datasets for use with the sesame package, which performs gene set enrichment analysis.
8. SummarizedExperiment: Stores and manipulates summarized assay data together with row and column metadata.
9. biomaRt: Provides an interface to BioMart databases, enabling retrieval of genomic data and annotations.

**These packages can be installed using install.packages() and loaded using library()**

## Preprocessing of the Dataset
1. Dataset is uploaded directly onto RStudio using TCGABiolink, then preproccessed by: Subsetting into the designated groups 
2. Normalization
3. Filteration with a quantile cut off 0.25

## Differential Expressionn Analysis
1. Mutant and WT Samples are grouped together
2. Treatment level DEA is performed
3. Ouput is used as data to generate heatmaps comparing the clustering of males vs females in both tumor types

## Functional Enrichment Analysis
1. The up and downregulated genes for all 4 subsets are extracted individually
2. The data undergoes functional enrichment analysis (particularly gene ontology), selecting the top 5 pathways.
3. Results are visualized in a barplot form and the PDF is saved on the computed
