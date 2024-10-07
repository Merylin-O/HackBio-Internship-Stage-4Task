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
SelectedBarcodes <- c(subset(LGG.metaData, IDH_status == "Mutant" & Expression_cluster == "LGr1")$Barcode,
                      subset(LGG.metaData, IDH_status == "WT" & Expression_cluster == "LGr1")$Barcode,
                      subset(LGG.metaData, IDH_status == "Mutant" & Expression_cluster == "LGr4")$Barcode,
                      subset(LGG.metaData, IDH_status == "WT" & Expression_cluster == "LGr4")$Barcode)

#Retrieve the unstranded data for the selected barcodes
Selected_LGG_Data <- LGG.Raw$unstranded[, c(SelectedBarcodes)]
dim(Selected_LGG_Data)
View(Selected_LGG_Data)

#Normalize Data
LGG_normalized <- TCGAanalyze_Normalization(tabDF = Selected_LGG_Data, geneInfo = geneInfoHT, method = "geneLength")

# Then Filter
LGG_filtered <- TCGAanalyze_Filtering(tabDF = LGG_normalized,
                                      method = "quantile",
                                      qnt.cut = 0.25)

View(LGG_filtered)
dim(LGG_filtered)


# Create matrices for IDH Mutant samples
mat1_Mutant_LGr1 <- LGG_filtered[, colnames(LGG_filtered) %in% subset(LGG.metaData, IDH_status == "Mutant" & Expression_cluster == "LGr1")$Barcode]
mat2_Mutant_LGr4 <- LGG_filtered[, colnames(LGG_filtered) %in% subset(LGG.metaData, IDH_status == "Mutant" & Expression_cluster == "LGr4")$Barcode]

# Create matrices for IDH Wild type samples 
mat1_Wild_LGr1 <- LGG_filtered[, colnames(LGG_filtered) %in% subset(LGG.metaData, IDH_status == "WT" & Expression_cluster == "LGr1")$Barcode]
mat2_Wild_LGr4 <- LGG_filtered[, colnames(LGG_filtered) %in% subset(LGG.metaData, IDH_status == "WT" & Expression_cluster == "LGr4")$Barcode]


# Perform DEA for IDH Mutant Samples
results_IDH_mutant <- TCGAanalyze_DEA(mat1 = mat1_Mutant_LGr1,
                                      mat2 = mat2_Mutant_LGr4,
                                      Cond1type = "Mutant LGr1",
                                      Cond2type = "Mutant LGr4",
                                      pipeline = "edgeR",
                                      fdr.cut = 0.01,
                                      logFC.cut = 2)

# Perform DEA for IDH Wild type Samples
results_IDH_WT <- TCGAanalyze_DEA(mat1 = mat1_Wild_LGr1,
                                  mat2 = mat2_Wild_LGr4,
                                  Cond1type = "WT LGr1",
                                  Cond2type = "WT LGr4",
                                  pipeline = "edgeR",
                                  fdr.cut = 0.01,
                                  logFC.cut = 2)

#Volcano Plot
plot(results_IDH_mutant$logFC, -log10(results_IDH_mutant$FDR))
plot(results_IDH_WT$logFC, -log10(results_IDH_WT$FDR))

#DEA with expression levels
results_Mut.level <- TCGAanalyze_LevelTab(results_IDH_mutant,"Mutant LGr1","Mutant LGr4", mat1_Mutant_LGr1, mat2_Mutant_LGr4)
results_WT.level <- TCGAanalyze_LevelTab(results_IDH_WT,"WT LGr1","WT LGr4", mat1_Wild_LGr1, mat2_Wild_LGr4)

head(results_Mut.level)
head(results_WT.level)

# Create vectors for the barcodes
IDH_Mut1_barcodes <- subset(LGG.metaData, IDH_status == "Mutant" & Expression_cluster == "LGr1")$Barcode
IDH_Mut4_barcodes <- subset(LGG.metaData, IDH_status == "Mutant" & Expression_cluster == "LGr4")$Barcode
IDH_WT1_barcodes <- subset(LGG.metaData, IDH_status == "WT" & Expression_cluster == "LGr1")$Barcode
IDH_WT4_barcodes <- subset(LGG.metaData, IDH_status == "WT" & Expression_cluster == "LGr4")$Barcode

# Combine the barcodes for heatmap data
selected_barcodes_Mutant <- c(IDH_Mut1_barcodes, IDH_Mut4_barcodes)
selected_barcodes_WT <- c(IDH_WT1_barcodes, IDH_WT4_barcodes)

# Subset the heatmap data
heat.data.Mut <- LGG_filtered[rownames(results_Mut.level), selected_barcodes_Mutant]
heat.data.WT <- LGG_filtered[rownames(results_WT.level), selected_barcodes_WT]

#View number of columns of heatmap data for IDH Mutant samples
ncol(heat.data.Mut) # To ensure that the total length of sample_colors is equal to ncol
colnames(heat.data.Mut)

# Create a color vector for the IDH mutant samples
sample_colorsM <- c(rep("midnightblue", 115), rep("red4", 5)) #Check the number of the barcodes for LGr 1 and 4 then adjust. midgnight blue is for LGr1 and red4 is for LGr

#Visualization of DEA for IDH Mutant Samples
heatmap.2(x = as.matrix(heat.data.Mut),
          col = hcl.colors(10, palette = 'Blue-Red 3'),  # Corrected the closing parenthesis here
          Rowv = FALSE, Colv = TRUE,
          scale = 'row',
          sepcolor = 'black',
          trace = "none",
          key = TRUE,
          dendrogram = "col",
          cexRow = 0.5, cexCol = 1,
          main = "Heatmap of IDH Mutation Expression Level",
          na.color = 'black',
          ColSideColors = sample_colorsM)

#Add a Legend
legend("topright", legend = c("LGr1", "LGr4"),
       fill = c("midnightblue", "red4"),
       border = "black",  # Black outline
       bty = "n")  # No box around the legend

#View number of columns of heatmap data for IDH Mutant samples
ncol(heat.data.WT) # To ensure that the total length of sample_colors is equal to ncol

# Create a color vector for the IDH mutant samplesamples
sample_colorsWT <- c(rep("midnightblue", 7), rep("red4", 71)) # Check the number of the barcodes for LGr 1 and 4 then adjust. midgnight blue is for LGr1 and red4 is for LGr

#Visualization of DEA for IDH Wild type Samples
heatmap.2(x = as.matrix(heat.data.WT),
          col = hcl.colors(10, palette = 'Blue-Red 3'),  # Custom color palette
          Rowv = FALSE, Colv = TRUE,
          scale = 'row',
          sepcolor = 'black',
          trace = "none",
          key = TRUE,
          dendrogram = "col",
          cexRow = 0.5, cexCol = 1,
          main = "Heatmap of IDH Wild Type Expression Level",
          na.color = 'black',
          ColSideColors = sample_colorsWT)

#Add a Legend
legend("topright", legend = c("LGr1", "LGr4"),
       fill = c("midnightblue", "red4"),
       border = "black",  # Black outline
       bty = "n")

# Upregulated and downregulated genes for LGr 1 & 4 IDH mutation
upreg_LGr1_Mut <- rownames(subset(results_Mut.level, logFC > 2 & `Mutant LGr1` > 0))
upreg_LGr4_Mut <- rownames(subset(results_Mut.level, logFC > 2 & `Mutant LGr4` > 0))
downreg_LGr1_Mut <- rownames(subset(results_Mut.level, logFC < -2 & `Mutant LGr1` > 0))
downreg_LGr4_Mut <- rownames(subset(results_Mut.level, logFC < -2 & `Mutant LGr4` > 0))

# Upregulated and downregulated genes for LGr 1 & 4 IDH WT
upreg_LGr1_WT <- rownames(subset(results_WT.level, logFC > 2 & `WT LGr1` > 0))
upreg_LGr4_WT <- rownames(subset(results_WT.level, logFC > 2 & `WT LGr4` > 0))
downreg_LGr1_WT <- rownames(subset(results_WT.level, logFC < -2 & `WT LGr1` > 0))
downreg_LGr4_WT <- rownames(subset(results_WT.level, logFC < -2 & `WT LGr4` > 0))


library("biomaRt")
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# Convert Ensemble IDs to gene symbols
up_LGr1_Mut_symbols <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = upreg_LGr1_Mut,
  mart = mart
)

upreg_LGr4_Mut_symbols <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = upreg_LGr4_Mut,
  mart = mart
)

downreg_LGr1_Mut_symbols <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = downreg_LGr1_Mut,
  mart = mart
)

downreg_LGr4_Mut_symbols <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = downreg_LGr4_Mut,
  mart = mart
)

upreg_LGr1_WT_symbols <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = upreg_LGr1_WT,
  mart = mart
)

upreg_LGr4_WT_symbols <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = upreg_LGr4_WT,
  mart = mart
)

downreg_LGr1_WT_symbols <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = downreg_LGr1_WT,
  mart = mart
)

downreg_LGr4_WT_symbols <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = downreg_LGr4_WT,
  mart = mart
)

dev.off()

# Functional Enrichment Analysis

up.EA.Muta.LGr1 <- TCGAanalyze_EAcomplete(TFname = "Mutant LGr1", RegulonList = up_LGr1_Mut_symbols$hgnc_symbol)
up.EA.Muta.LGr4 <- TCGAanalyze_EAcomplete(TFname = "Mutant LGr4", RegulonList = upreg_LGr4_Mut_symbols$hgnc_symbol)
up.EA.WT.LGr1 <- TCGAanalyze_EAcomplete(TFname = "WT LGr1", RegulonList = upreg_LGr1_WT_symbols$hgnc_symbol)
up.EA.WT.LGr4 <- TCGAanalyze_EAcomplete(TFname = "WT LGr4", RegulonList = upreg_LGr4_WT_symbols$hgnc_symbol)
down.EA.Muta.LGr1 <- TCGAanalyze_EAcomplete(TFname = "Mutant LGr1 (Downregulated)", RegulonList = downreg_LGr1_Mut_symbols$hgnc_symbol)
down.EA.Muta.LGr4 <- TCGAanalyze_EAcomplete(TFname = "Mutant LGr4 (Downregulated)", RegulonList = downreg_LGr4_Mut_symbols$hgnc_symbol)
down.EA.WT.LGr1 <- TCGAanalyze_EAcomplete(TFname = "WT LGr1 (Downregulated)", RegulonList = downreg_LGr1_WT_symbols$hgnc_symbol)
down.EA.WT.LGr4 <- TCGAanalyze_EAcomplete(TFname = "WT LGr4 (Downregulated)", RegulonList = downreg_LGr4_WT_symbols$hgnc_symbol)

#Visualization of the enrichment analysis
UP.Muta.LGr1 <- TCGAvisualize_EAbarplot(tf = rownames(up.EA.Muta.LGr1$ResBP),
                                        GOBPTab = up.EA.Muta.LGr1$ResBP,
                                        GOCCTab = up.EA.Muta.LGr1$ResCC,
                                        GOMFTab = up.EA.Muta.LGr1$ResMF,
                                        PathTab = up.EA.Muta.LGr1$ResPat,
                                        nRGTab = up.EA.Muta.LGr1,
                                        nBar = 10,
                                        text.size = 2,
                                        fig.width = 30,
                                        fig.height = 15)

UP.Muta.LGr4 <- TCGAvisualize_EAbarplot(tf = rownames(up.EA.Muta.LGr4$ResBP),
                                        GOBPTab = up.EA.Muta.LGr4$ResBP,
                                        GOCCTab = up.EA.Muta.LGr4$ResCC,
                                        GOMFTab = up.EA.Muta.LGr4$ResMF,
                                        PathTab = up.EA.Muta.LGr4$ResPat,
                                        nRGTab = up.EA.Muta.LGr4,
                                        nBar = 10,
                                        text.size = 2,
                                        fig.width = 30,
                                        fig.height = 15)

UP.WT.LGr1 <- TCGAvisualize_EAbarplot(tf = rownames(up.EA.WT.LGr1$ResBP),
                                      GOBPTab = up.EA.WT.LGr1$ResBP,
                                      GOCCTab = up.EA.WT.LGr1$ResCC,
                                      GOMFTab = up.EA.WT.LGr1$ResMF,
                                      PathTab = up.EA.WT.LGr1$ResPat,
                                      nRGTab = up.EA.WT.LGr1,
                                      nBar = 10,
                                      text.size = 2,
                                      fig.width = 30,
                                      fig.height = 15)

UP.WT.LGr4 <- TCGAvisualize_EAbarplot(tf = rownames(up.EA.WT.LGr4$ResBP),
                                      GOBPTab = up.EA.WT.LGr4$ResBP,
                                      GOCCTab = up.EA.WT.LGr4$ResCC,
                                      GOMFTab = up.EA.WT.LGr4$ResMF,
                                      PathTab = up.EA.WT.LGr4$ResPat,
                                      nRGTab = up.EA.WT.LGr4,
                                      nBar = 10,
                                      text.size = 2,
                                      fig.width = 30,
                                      fig.height = 15)

DOWN.Muta.LGr1 <- TCGAvisualize_EAbarplot(tf = rownames(down.EA.Muta.LGr1$ResBP),
                                          GOBPTab = down.EA.Muta.LGr1$ResBP,
                                          GOCCTab = down.EA.Muta.LGr1$ResCC,
                                          GOMFTab = down.EA.Muta.LGr1$ResMF,
                                          PathTab = down.EA.Muta.LGr1$ResPat,
                                          nRGTab = down.EA.Muta.LGr1,
                                          nBar = 10,
                                          text.size = 2,
                                          fig.width = 30,
                                          fig.height = 15)

DOWN.Muta.LGr4 <- TCGAvisualize_EAbarplot(tf = rownames(down.EA.Muta.LGr4$ResBP),
                                          GOBPTab = down.EA.Muta.LGr4$ResBP,
                                          GOCCTab = down.EA.Muta.LGr4$ResCC,
                                          GOMFTab = down.EA.Muta.LGr4$ResMF,
                                          PathTab = down.EA.Muta.LGr4$ResPat,
                                          nRGTab = down.EA.Muta.LGr4,
                                          nBar = 10,
                                          text.size = 2,
                                          fig.width = 30,
                                          fig.height = 15)

DOWN.WT.LGr1 <- TCGAvisualize_EAbarplot(tf = rownames(down.EA.WT.LGr1$ResBP),
                                        GOBPTab = down.EA.WT.LGr1$ResBP,
                                        GOCCTab = down.EA.WT.LGr1$ResCC,
                                        GOMFTab = down.EA.WT.LGr1$ResMF,
                                        PathTab = down.EA.WT.LGr1$ResPat,
                                        nRGTab = down.EA.WT.LGr1,
                                        nBar = 10,
                                        text.size = 2,
                                        fig.width = 30,
                                        fig.height = 15)

DOWN.WT.LGr4 <- TCGAvisualize_EAbarplot(tf = rownames(down.EA.WT.LGr4$ResBP),
                                        GOBPTab = down.EA.WT.LGr4$ResBP,
                                        GOCCTab = down.EA.WT.LGr4$ResCC,
                                        GOMFTab = down.EA.WT.LGr4$ResMF,
                                        PathTab = down.EA.WT.LGr4$ResPat,
                                        nRGTab = down.EA.WT.LGr4,
                                        nBar = 10,
                                        text.size = 2,
                                        fig.width = 30,
                                        fig.height = 15)
