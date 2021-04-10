library(Seurat)
library(SeuratData)
library(tidyverse)
library(cowplot)
library(patchwork)
rm(list = ls())
# InstallData("panc8")
# install.packages("../tmp/panc8.SeuratData_3.0.2.tar.gz", repos = NULL
library(panc8.SeuratData)
data("panc8")
pancreas.list <- SplitObject(panc8, split.by = "tech")
# pancreas.list <- pancreas.list[c("celseq", "celseq2", "fluidigmc1", "smartseq2")]
for (i in 1:length(pancreas.list)) {
  pancreas.list[[i]] <- NormalizeData(pancreas.list[[i]], verbose = FALSE)
  pancreas.list[[i]] <- FindVariableFeatures(pancreas.list[[i]], selection.method = "vst", nfeatures = 2000, 
                                             verbose = FALSE)
}

reference.list <- pancreas.list[c("celseq", "celseq2", "smartseq2")]
pancreas.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)

# switch to integrated assay. The variable features of this assay are automatically set during
# IntegrateData
DefaultAssay(pancreas.integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)
pancreas.integrated <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:30, verbose = FALSE)
p1 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "tech")
p2 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE) + 
  NoLegend()
p1 + p2

# 细胞分类------------
pancreas.query <- pancreas.list[["indrop"]]
pancreas.anchors <- FindTransferAnchors(reference = pancreas.integrated, query = pancreas.query, 
                                        dims = 1:30)
predictions <- TransferData(anchorset = pancreas.anchors, refdata = pancreas.integrated$celltype, 
                            dims = 1:30)
pancreas.query <- AddMetaData(pancreas.query, metadata = predictions)
pancreas.query$prediction.match <- pancreas.query$predicted.id == pancreas.query$celltype
table(pancreas.query$prediction.match)
table(pancreas.query$predicted.id)
VlnPlot(pancreas.query, c("REG1A", "PPY", "SST", "GHRL", "VWF", "SOX10"), group.by = "predicted.id")

pancreas.integrated <- RunUMAP(pancreas.integrated, dims = 1:30, reduction = "pca", return.model = TRUE)
pancreas.query <- MapQuery(anchorset = pancreas.anchors, reference = pancreas.integrated, query = pancreas.query, 
                           refdata = list(celltype = "celltype"), reference.reduction = "pca", reduction.model = "umap")

p1 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "celltype", label = TRUE, label.size = 3, 
              repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(pancreas.query, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE, 
              label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
p1 + p2
