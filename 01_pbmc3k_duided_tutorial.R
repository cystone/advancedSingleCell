suppressMessages(library(tidyverse))
suppressMessages(library(pacman))
suppressMessages(library(data.table))
suppressMessages(library(Seurat))
suppressMessages(library(patchwork))
options(stringsAsFactors = F)
rm(list = ls())

inDir = '~/BioFiles/pbmc3k/'
# Load the PBMC dataset
inFile = paste0(inDir, 'filtered_gene_bc_matrices/hg19/')
pbmc.data <- Read10X(data.dir = inFile)
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc
# Lets examine a few genes in the first thirty cells
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]

#-----标准预处理工作流程-----------
#下面的步骤包含了 Seurat 的 scRNA-seq 数据的标准预处理流程。
#这些代表了基于 QC 指标的单元的选择和筛选、数据规范化和缩放以及高度可变特征的检测。

#-----------QC和选择细胞-----------
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc[["percent.ribo"]] <- PercentageFeatureSet(pbmc, pattern = "^RP[SL]")
# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nFeature_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(pbmc, feature1 = "nFeature_RNA", feature2 = "percent.ribo")
plot1 + plot2
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#----------规范化数据----------
#从数据集中删除不需要的单元格后，下一步是规范化数据。默认情况下，
#我们使用一个全局缩放标准化方法“ LogNormalize”，
#该方法通过总表达式对每个单元格的特征表达式度量值进行标准化，
#将其乘以一个比例因子(默认为10,000) ，并对结果进行 log-transforms。
#规范化值存储在 pbmc[["RNA"]]@data。
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
#---------特征选择----------
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#-----------归一化数据--------
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
# pbmc <- SCTransform(pbmc,variable.features.n = 3000, vars.to.regress = "percent.mt")
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)

#--确定数据集的维数------
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)
ElbowPlot(pbmc)

#---------细胞聚集---------
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

#---------非线性降维----------
pbmc <- RunUMAP(pbmc, dims = 1:10)
set.seed(123)
DimPlot(pbmc, reduction = "umap")
saveRDS(pbmc, file = "./processData/01_pbmc_tutorial.rds")

#-------聚类生物标志物--------
# find all markers of cluster 1
cluster1.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
head(cluster1.markers, n = 5)
# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)
# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

# 修拉有几个差异表达式的测试，可以通过 test.use 参数设置(详见我们的 DE vignette)。
#例如，ROC 测试返回任何单个标记的“分类能力”(范围从0-random 到1-perfect)。
cluster1.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
# you can plot raw counts as well
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", 
                               "CD8A"))
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

FeaturePlot(pbmc, features = c("FCGR3A", "MS4A7", "FCER1A", "CST3", "PPBP"))
#------------cell type identity-------------
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono", 
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

saveRDS(pbmc, file = "./processData/01_pbmc3k_final.rds")