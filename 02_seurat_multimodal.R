suppressMessages(library(tidyverse))
suppressMessages(library(pacman))
suppressMessages(library(data.table))
suppressMessages(library(Seurat))
suppressMessages(library(patchwork))
options(stringsAsFactors = F)
rm(list = ls())
# Load in the RNA UMI matrix

# Note that this dataset also contains ~5% of mouse cells, which we can use as negative controls
# for the protein measurements. For this reason, the gene expression matrix has HUMAN_ or MOUSE_
# appended to the beginning of each gene.
inFile = "~/BioFiles/GSE100866_CBMC/GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz"
cbmc.rna <- as.sparse(read.csv(file = inFile, sep = ",", 
                               header = TRUE, row.names = 1))
# To make life a bit easier going forward, we're going to discard all but the top 100 most
# highly expressed mouse genes, and remove the 'HUMAN_' from the CITE-seq prefix
cbmc.rna <- CollapseSpeciesExpressionMatrix(cbmc.rna)
# Load in the ADT UMI matrix
inFile = "~/BioFiles/GSE100866_CBMC/GSE100866_CBMC_8K_13AB_10X-ADT_umi.csv.gz"
cbmc.adt <- as.sparse(read.csv(file = inFile, sep = ",", 
                               header = TRUE, row.names = 1))
# Note that since measurements were made in the same cells, the two matrices have identical
# column names
all.equal(colnames(cbmc.rna), colnames(cbmc.adt))

#---------Seurat: 添加 RNA 和蛋白质数据--------
# creates a Seurat object based on the scRNA-seq data
cbmc <- CreateSeuratObject(counts = cbmc.rna)
# We can see that by default, the cbmc object contains an assay storing RNA measurement
Assays(cbmc)
# create a new assay to store ADT information
adt_assay <- CreateAssayObject(counts = cbmc.adt)
# add this assay to the previously created Seurat object
cbmc[["ADT"]] <- adt_assay
Assays(cbmc)
# Extract a list of features measured in the ADT assay
rownames(cbmc[["ADT"]])
# List the current default assay
DefaultAssay(cbmc)
# Switch the default to ADT
DefaultAssay(cbmc) <- "ADT"
DefaultAssay(cbmc)
#-------------Cluster cell-----------------
# Note that all operations below are performed on the RNA assay Set and verify that the default
# assay is RNA
DefaultAssay(cbmc) <- "RNA"
DefaultAssay(cbmc)
# perform visualization and clustering steps
cbmc <- NormalizeData(cbmc)
cbmc <- FindVariableFeatures(cbmc)
cbmc <- ScaleData(cbmc)
cbmc <- RunPCA(cbmc, verbose = FALSE)
cbmc <- FindNeighbors(cbmc, dims = 1:30)
cbmc <- FindClusters(cbmc, resolution = 0.8, verbose = FALSE)
cbmc <- RunUMAP(cbmc, dims = 1:30)
DimPlot(cbmc, label = TRUE)

#---------------并排查看的多种模式---------
# Normalize ADT data,
DefaultAssay(cbmc) <- "ADT"
cbmc <- NormalizeData(cbmc, normalization.method = "CLR", margin = 2)
DefaultAssay(cbmc) <- "RNA"

# Note that the following command is an alternative but returns the same result
cbmc <- NormalizeData(cbmc, normalization.method = "CLR", margin = 2, assay = "ADT")

# Now, we will visualize CD14 levels for RNA and protein By setting the default assay, we can
# visualize one or the other
DefaultAssay(cbmc) <- "ADT"
p1 <- FeaturePlot(cbmc, "CD19", cols = c("lightgrey", "darkgreen")) + ggtitle("CD19 protein")
DefaultAssay(cbmc) <- "RNA"
p2 <- FeaturePlot(cbmc, "CD19") + ggtitle("CD19 RNA")
# place plots side-by-side
p1 | p2
# for the RNA and protein assays
Key(cbmc[["RNA"]])
Key(cbmc[["ADT"]])
# Now, we can include the key in the feature name, which overrides the default assay
p1 <- FeaturePlot(cbmc, "adt_CD19", cols = c("lightgrey", "darkgreen")) + ggtitle("CD19 protein")
p2 <- FeaturePlot(cbmc, "rna_CD19") + ggtitle("CD19 RNA")
p1 | p2

#-----------识别细胞表面标记-----------
# surface
VlnPlot(cbmc, "adt_CD19")
# we can also identify alternative protein and RNA markers for this cluster through differential
# expression
adt_markers <- FindMarkers(cbmc, ident.1 = 5, assay = "ADT")
rna_markers <- FindMarkers(cbmc, ident.1 = 5, assay = "RNA")
head(adt_markers)
head(rna_markers)
#-----------更多可视化---------
# Draw ADT scatter plots (like biaxial plots for FACS). Note that you can even 'gate' cells if
# desired by using HoverLocator and FeatureLocator
FeatureScatter(cbmc, feature1 = "adt_CD19", feature2 = "adt_CD3")
# view relationship between protein and RNA
FeatureScatter(cbmc, feature1 = "adt_CD3", feature2 = "rna_CD3E")
FeatureScatter(cbmc, feature1 = "adt_CD4", feature2 = "adt_CD8")
# number in cells, which significantly reduces 'drop-out' in ADT data
FeatureScatter(cbmc, feature1 = "adt_CD4", feature2 = "adt_CD8", slot = "counts")
#----------------10X 多模态数据---------------
pbmc10k.data <- Read10X(data.dir = "../data/pbmc10k/filtered_feature_bc_matrix/")
rownames(x = pbmc10k.data[["Antibody Capture"]]) <- gsub(pattern = "_[control_]*TotalSeqB", replacement = "", 
                                                         x = rownames(x = pbmc10k.data[["Antibody Capture"]]))

pbmc10k <- CreateSeuratObject(counts = pbmc10k.data[["Gene Expression"]], min.cells = 3, min.features = 200)
pbmc10k <- NormalizeData(pbmc10k)
pbmc10k[["ADT"]] <- CreateAssayObject(pbmc10k.data[["Antibody Capture"]][, colnames(x = pbmc10k)])
pbmc10k <- NormalizeData(pbmc10k, assay = "ADT", normalization.method = "CLR")

plot1 <- FeatureScatter(pbmc10k, feature1 = "adt_CD19", feature2 = "adt_CD3", pt.size = 1)
plot2 <- FeatureScatter(pbmc10k, feature1 = "adt_CD4", feature2 = "adt_CD8a", pt.size = 1)
plot3 <- FeatureScatter(pbmc10k, feature1 = "adt_CD3", feature2 = "CD3E", pt.size = 1)
(plot1 + plot2 + plot3) & NoLegend()