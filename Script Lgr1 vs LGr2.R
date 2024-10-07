# Load Important Packages
library("BiocManager")
library("TCGAbiolinks")
library("limma")
library("edgeR")
library("EDASeq")
library("gplots")
library("sesameData")
library("SummarizedExperiment")

# Get an overview on  LGG data
getProjectSummary("TCGA-LGG")

# Query the LGG Data
LGG <- GDCquery(project = "TCGA-LGG",
                data.category = "Transcriptome Profiling",
                data.type = "Gene Expression Quantification")

# Download the dataset
GDCdownload(LGG)

# Prepare the dataset
LGG.Data <- GDCprepare(LGG)

# Check the output of the data for LGG
head(LGG.Data)
View(LGG.Data)

# Explore the data
LGG.Data$gender
LGG.Data$barcode
LGG.Data$patient
LGG$
  
# Check the Subsets of the metadata and add them in one data frame
LGG.metaData <- data.frame("IDH_status" = LGG.Data$paper_IDH.status,
                             "Expression_cluster" = LGG.Data$paper_Pan.Glioma.RNA.Expression.Cluster,
                             "Barcode" = LGG.Data$barcode)

# Extract the raw data from the prepared dataset
LGG.Raw <- assays(LGG.Data)

#Select unstranded data
dim(LGG.Raw$unstranded)
View(LGG.Raw$unstranded)

# Subset Meta Data
SelectedBarcodes2 <- c(subset(LGG.metaData, IDH_status == "Mutant" & Expression_cluster == "LGr1")$Barcode,
                      subset(LGG.metaData, IDH_status == "WT" & Expression_cluster == "LGr1")$Barcode,
                      subset(LGG.metaData, IDH_status == "Mutant" & Expression_cluster == "LGr2")$Barcode,
                      subset(LGG.metaData, IDH_status == "WT" & Expression_cluster == "LGr2")$Barcode)

#Retrieve the unstranded data for the selected barcodes
Selected_LGG_Data2 <- LGG.Raw$unstranded[, c(SelectedBarcodes2)]
dim(Selected_LGG_Data2)
View(Selected_LGG_Data2)

#Normalize Data
LGG_normalized2 <- TCGAanalyze_Normalization(tabDF = Selected_LGG_Data2, geneInfo = geneInfoHT, method = "geneLength")

# Then Filter
LGG_filtered2 <- TCGAanalyze_Filtering(tabDF = LGG_normalized2,
                                      method = "quantile",
                                      qnt.cut = 0.25)

View(LGG_filtered2)
dim(LGG_filtered2)


# Create matrices for IDH Mutant samples
mat1_Mutant_LGr1a <- LGG_filtered2[, colnames(LGG_filtered2) %in% subset(LGG.metaData, IDH_status == "Mutant" & Expression_cluster == "LGr1")$Barcode]
mat2_Mutant_LGr2 <- LGG_filtered2[, colnames(LGG_filtered2) %in% subset(LGG.metaData, IDH_status == "Mutant" & Expression_cluster == "LGr2")$Barcode]

# Create matrices for IDH Wild type samples 
mat1_Wild_LGr1a <- LGG_filtered2[, colnames(LGG_filtered2) %in% subset(LGG.metaData, IDH_status == "WT" & Expression_cluster == "LGr1")$Barcode]
mat2_Wild_LGr2 <- LGG_filtered2[, colnames(LGG_filtered2) %in% subset(LGG.metaData, IDH_status == "WT" & Expression_cluster == "LGr2")$Barcode]


# Perform DEA for IDH Mutant Samples
results_IDH_mutant2 <- TCGAanalyze_DEA(mat1 = mat1_Mutant_LGr1a,
                                      mat2 = mat2_Mutant_LGr2,
                                      Cond1type = "Mutant LGr1",
                                      Cond2type = "Mutant LGr2",
                                      pipeline = "edgeR",
                                      fdr.cut = 0.01,
                                      logFC.cut = 2)

# Perform DEA for IDH Wild type Samples
results_IDH_WT2 <- TCGAanalyze_DEA(mat1 = mat1_Wild_LGr1a,
                                  mat2 = mat2_Wild_LGr2,
                                  Cond1type = "WT LGr1",
                                  Cond2type = "WT LGr42",
                                  pipeline = "edgeR",
                                  fdr.cut = 0.01,
                                  logFC.cut = 2)

#Volcano Plot
plot(results_IDH_mutant2$logFC, -log10(results_IDH_mutant2$FDR))
plot(results_IDH_WT2$logFC, -log10(results_IDH_WT2$FDR))

#DEA with expression levels
results_Mut.level2 <- TCGAanalyze_LevelTab(results_IDH_mutant2,"Mutant LGr1","Mutant LGr2", mat1_Mutant_LGr1a, mat2_Mutant_LGr2)
results_WT.level2 <- TCGAanalyze_LevelTab(results_IDH_WT2,"WT LGr1","WT LGr2", mat1_Wild_LGr1a, mat2_Wild_LGr2)

head(results_Mut.level)
head(results_WT.level)

# Create vectors for the barcodes
IDH_Mut1_barcodes2 <- subset(LGG.metaData, IDH_status == "Mutant" & Expression_cluster == "LGr1")$Barcode
IDH_Mut2_barcodes <- subset(LGG.metaData, IDH_status == "Mutant" & Expression_cluster == "LGr2")$Barcode
IDH_WT1_barcodes2 <- subset(LGG.metaData, IDH_status == "WT" & Expression_cluster == "LGr1")$Barcode
IDH_WT2_barcodes <- subset(LGG.metaData, IDH_status == "WT" & Expression_cluster == "LGr2")$Barcode

# Combine the barcodes for heatmap data
selected_barcodes_Mutant2 <- c(IDH_Mut1_barcodes2, IDH_Mut2_barcodes)
selected_barcodes_WT2 <- c(IDH_WT1_barcodes2, IDH_WT2_barcodes)

# Subset the heatmap data
heat.data.Mut2 <- LGG_filtered2[rownames(results_Mut.level2), selected_barcodes_Mutant2]
heat.data.WT2 <- LGG_filtered2[rownames(results_WT.level2), selected_barcodes_WT2]

#View number of columns of heatmap data for IDH Mutant samples
ncol(heat.data.Mut2) # To ensure that the total length of sample_colors is equal to ncol
colnames(heat.data.Mut)

# Create a color vector for the IDH mutant samples
sample_colorsMT2 <- c(rep("midnightblue", 115), rep("red4", 74)) #Check the number of the barcodes for LGr 1 and 2 then adjust

#Visualization of DEA for IDH Mutant Samples
heatmap.2(x = as.matrix(heat.data.Mut2),
          col = hcl.colors(10, palette = 'Blue-Red 3'),  # Corrected the closing parenthesis here
          Rowv = FALSE, Colv = TRUE,
          scale = 'row',
          sepcolor = 'black',
          trace = "none",
          key = TRUE,
          dendrogram = "col",
          cexRow = 0.5, cexCol = 1,
          main = "Heatmap of IDH Mutation Expression Level LGr1 vs LGr2",
          na.color = 'black',
          ColSideColors = sample_colorsMT2)

#Add a Legend
legend("topright", legend = c("LGr1", "LGr2"),
       fill = c("midnightblue", "red4"),
       border = "black",  # Black outline
       bty = "n")  # No box around the legend

#View number of columns of heatmap data for IDH Mutant samples
ncol(heat.data.WT2) # To ensure that the total length of sample_colors is equal to ncol

# Create a color vector for the IDH mutant samplesamples
sample_colorsWT2 <- c(rep("midnightblue", 7), rep("red4", 12)) # Check the number of the barcodes for LGr 1 and 2 then adjust.

#Visualization of DEA for IDH Wild type Samples
heatmap.2(x = as.matrix(heat.data.WT2),
          col = hcl.colors(10, palette = 'Blue-Red 3'),  # Custom color palette
          Rowv = FALSE, Colv = TRUE,
          scale = 'row',
          sepcolor = 'black',
          trace = "none",
          key = TRUE,
          dendrogram = "col",
          cexRow = 0.5, cexCol = 1,
          main = "Heatmap of IDH Wild Type Expression Level LGr1 vs LGr2",
          na.color = 'black',
          ColSideColors = sample_colorsWT2)

#Add a Legend
legend("topright", legend = c("LGr1", "LGr2"),
       fill = c("midnightblue", "red4"),
       border = "black",  # Black outline
       bty = "n")

# Upregulated and downregulated genes for LGr 1 & 2 IDH mutation
upreg_LGr1_Mut2 <- rownames(subset(results_Mut.level2, logFC > 2 & `Mutant LGr1` > 0))
upreg_LGr2_Mut <- rownames(subset(results_Mut.level2, logFC > 2 & `Mutant LGr2` > 0))
downreg_LGr1_Mut2 <- rownames(subset(results_Mut.level2, logFC < -2 & `Mutant LGr1` > 0))
downreg_LGr2_Mut <- rownames(subset(results_Mut.level2, logFC < -2 & `Mutant LGr2` > 0))

# Upregulated and downregulated genes for LGr 1 & 2 IDH WT
upreg_LGr1_WT2 <- rownames(subset(results_WT.level2, logFC > 2 & `WT LGr1` > 0))
upreg_LGr2_WT <- rownames(subset(results_WT.level2, logFC > 2 & `WT LGr2` > 0))
downreg_LGr1_WT2 <- rownames(subset(results_WT.level2, logFC < -2 & `WT LGr1` > 0))
downreg_LGr2_WT <- rownames(subset(results_WT.level2, logFC < -2 & `WT LGr2` > 0))


library("biomaRt")
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# Convert Ensemble IDs to gene symbols
up_LGr1_Mut_symbols2 <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = upreg_LGr1_Mut2,
  mart = mart
)

upreg_LGr2_Mut_symbols <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = upreg_LGr2_Mut,
  mart = mart
)

downreg_LGr1_Mut_symbols2 <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = downreg_LGr1_Mut2,
  mart = mart
)

downreg_LGr2_Mut_symbols <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = downreg_LGr2_Mut,
  mart = mart
)

upreg_LGr1_WT_symbols2 <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = upreg_LGr1_WT2,
  mart = mart
)

upreg_LGr2_WT_symbols <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = upreg_LGr2_WT,
  mart = mart
)

downreg_LGr1_WT_symbols2 <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = downreg_LGr1_WT2,
  mart = mart
)

downreg_LGr2_WT_symbols <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = downreg_LGr2_WT,
  mart = mart
)

dev.off()

# Functional Enrichment Analysis

up.EA.Muta.LGr1a <- TCGAanalyze_EAcomplete(TFname = "Mutant LGr1", RegulonList = up_LGr1_Mut_symbols2$hgnc_symbol)
up.EA.Muta.LGr2 <- TCGAanalyze_EAcomplete(TFname = "Mutant LGr2", RegulonList = upreg_LGr2_Mut_symbols$hgnc_symbol)
up.EA.WT.LGr1a <- TCGAanalyze_EAcomplete(TFname = "WT LGr1", RegulonList = upreg_LGr1_WT_symbols2$hgnc_symbol)
up.EA.WT.LGr2 <- TCGAanalyze_EAcomplete(TFname = "WT LGr2", RegulonList = upreg_LGr2_WT_symbols$hgnc_symbol)
down.EA.Muta.LGr1a <- TCGAanalyze_EAcomplete(TFname = "Mutant LGr1 (Downregulated)", RegulonList = downreg_LGr1_Mut_symbols2$hgnc_symbol)
down.EA.Muta.LGr2 <- TCGAanalyze_EAcomplete(TFname = "Mutant LGr2 (Downregulated)", RegulonList = downreg_LGr2_Mut_symbols$hgnc_symbol)
down.EA.WT.LGr1a <- TCGAanalyze_EAcomplete(TFname = "WT LGr1 (Downregulated)", RegulonList = downreg_LGr1_WT_symbols2$hgnc_symbol)
down.EA.WT.LGr2 <- TCGAanalyze_EAcomplete(TFname = "WT LGr2 (Downregulated)", RegulonList = downreg_LGr2_WT_symbols$hgnc_symbol)

#Visualization of the enrichment analysis
UP.Muta.LGr1a <- TCGAvisualize_EAbarplot(tf = rownames(up.EA.Muta.LGr1a$ResBP),
                                        GOBPTab = up.EA.Muta.LGr1a$ResBP,
                                        GOCCTab = up.EA.Muta.LGr1a$ResCC,
                                        GOMFTab = up.EA.Muta.LGr1a$ResMF,
                                        PathTab = up.EA.Muta.LGr1a$ResPat,
                                        nRGTab = up.EA.Muta.LGr1a,
                                        nBar = 10,
                                        text.size = 2,
                                        fig.width = 30,
                                        fig.height = 15)


UP.Muta.LGr2 <- TCGAvisualize_EAbarplot(tf = rownames(up.EA.Muta.LGr2$ResBP),
                                        GOBPTab = up.EA.Muta.LGr2$ResBP,
                                        GOCCTab = up.EA.Muta.LGr2$ResCC,
                                        GOMFTab = up.EA.Muta.LGr2$ResMF,
                                        PathTab = up.EA.Muta.LGr2$ResPat,
                                        nRGTab = up.EA.Muta.LGr2,
                                        nBar = 10,
                                        text.size = 2,
                                        fig.width = 30,
                                        fig.height = 15)


UP.WT.LGr1a <- TCGAvisualize_EAbarplot(tf = rownames(up.EA.WT.LGr1a$ResBP),
                                      GOBPTab = up.EA.WT.LGr1a$ResBP,
                                      GOCCTab = up.EA.WT.LGr1a$ResCC,
                                      GOMFTab = up.EA.WT.LGr1a$ResMF,
                                      PathTab = up.EA.WT.LGr1a$ResPat,
                                      nRGTab = up.EA.WT.LGr1a,
                                      nBar = 10,
                                      text.size = 2,
                                      fig.width = 30,
                                      fig.height = 15)


UP.WT.LGr2 <- TCGAvisualize_EAbarplot(tf = rownames(up.EA.WT.LGr2$ResBP),
                                      GOBPTab = up.EA.WT.LGr2$ResBP,
                                      GOCCTab = up.EA.WT.LGr2$ResCC,
                                      GOMFTab = up.EA.WT.LGr2$ResMF,
                                      PathTab = up.EA.WT.LGr2$ResPat,
                                      nRGTab = up.EA.WT.LGr2,
                                      nBar = 10,
                                      text.size = 2,
                                      fig.width = 30,
                                      fig.height = 15)


DOWN.Muta.LGr1a <- TCGAvisualize_EAbarplot(tf = rownames(down.EA.Muta.LGr1a$ResBP),
                                          GOBPTab = down.EA.Muta.LGr1a$ResBP,
                                          GOCCTab = down.EA.Muta.LGr1a$ResCC,
                                          GOMFTab = down.EA.Muta.LGr1a$ResMF,
                                          PathTab = down.EA.Muta.LGr1a$ResPat,
                                          nRGTab = down.EA.Muta.LGr1a,
                                          nBar = 10,
                                          text.size = 2,
                                          fig.width = 30,
                                          fig.height = 15)


DOWN.Muta.LGr2 <- TCGAvisualize_EAbarplot(tf = rownames(down.EA.Muta.LGr2$ResBP),
                                          GOBPTab = down.EA.Muta.LGr2$ResBP,
                                          GOCCTab = down.EA.Muta.LGr2$ResCC,
                                          GOMFTab = down.EA.Muta.LGr2$ResMF,
                                          PathTab = down.EA.Muta.LGr2$ResPat,
                                          nRGTab = down.EA.Muta.LGr2,
                                          nBar = 10,
                                          text.size = 2,
                                          fig.width = 30,
                                          fig.height = 15)


DOWN.WT.LGr1a <- TCGAvisualize_EAbarplot(tf = rownames(down.EA.WT.LGr1a$ResBP),
                                        GOBPTab = down.EA.WT.LGr1a$ResBP,
                                        GOCCTab = down.EA.WT.LGr1a$ResCC,
                                        GOMFTab = down.EA.WT.LGr1a$ResMF,
                                        PathTab = down.EA.WT.LGr1a$ResPat,
                                        nRGTab = down.EA.WT.LGr1a,
                                        nBar = 10,
                                        text.size = 2,
                                        fig.width = 30,
                                        fig.height = 15)


DOWN.WT.LGr2 <- TCGAvisualize_EAbarplot(tf = rownames(down.EA.WT.LGr2$ResBP),
                                        GOBPTab = down.EA.WT.LGr2$ResBP,
                                        GOCCTab = down.EA.WT.LGr2$ResCC,
                                        GOMFTab = down.EA.WT.LGr2$ResMF,
                                        PathTab = down.EA.WT.LGr2$ResPat,
                                        nRGTab = down.EA.WT.LGr2,
                                        nBar = 10,
                                        text.size = 2,
                                        fig.width = 30,
                                        fig.height = 15)
