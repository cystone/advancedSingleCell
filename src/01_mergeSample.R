library(Seurat)
library(tidyverse)
library(patchwork)
library(data.table)
dir.create('./output/01_mergeSample/cluster1', recursive=T)
dir.create('./output/01_mergeSample/cluster1', recursive=T)
dir.create('./output/01_mergeSample/cluster1', recursive=T)
set.seed(123)  #设置随机数种子，使结果可重复

#---------样本合并-----------
##使用目录向量合并
dir = c('./data/raw/GSE139324_RAW/GSM4138110', 
        './data/raw/GSE139324_RAW/GSM4138111',
        './data/raw/GSE139324_RAW/GSM4138128',
        './data/raw/GSE139324_RAW/GSM4138129',
        './data/raw/GSE139324_RAW/GSM4138148',
        './data/raw/GSE139324_RAW/GSM4138149',
        './data/raw/GSE139324_RAW/GSM4138162',
        './data/raw/GSE139324_RAW/GSM4138163',
        './data/raw/GSE139324_RAW/GSM4138168',
        './data/raw/GSE139324_RAW/GSM4138169')

for(i in dir){
  ind = str_split(i, '/',simplify=T)[5]
  indP = paste0(ind, '.')
  if(!dir.exists(i)){dir.create(i, recursive=T)}
  filern = list.files('./data/raw/GSE139324_RAW/', pattern=indP,include.dirs = F)
  if(length(filern == 3)){
    filname = str_split(filern, '_',simplify=T)[,5]
    file.rename(paste0('./data/raw/GSE139324_RAW/', filern),
                paste0('./data/raw/GSE139324_RAW/', ind, '/',filname))
  }else(print(paste0('please check ', ind)))
  file.rename(paste0('./data/raw/GSE139324_RAW/',ind, '/genes.tsv.gz'),
              paste0('./data/raw/GSE139324_RAW/', ind, '/features.tsv.gz'))
  
  }

names(dir) = c('HNC01PBMC', 'HNC01TIL', 'HNC10PBMC', 'HNC10TIL', 'HNC20PBMC', 
               'HNC20TIL', 'PBMC1', 'PBMC2', 'Tonsil1', 'Tonsil2')
counts <- Read10X(data.dir = dir)


scRNA1 = CreateSeuratObject(counts, min.cells=1)
dim(scRNA1)   #查看基因数和细胞总数
#[1] 23603 19750 
table(scRNA1@meta.data$orig.ident)  #查看每个样本的细胞数
#HNC01PBMC  HNC01TIL HNC10PBMC  HNC10TIL HNC20PBMC  HNC20TIL     PBMC1     PBMC2   Tonsil1   Tonsil2 
#     1725      1298      1750      1384      1530      1148      2445      2436      3325      2709

#使用merge函数合并seurat对象
scRNAlist <- list()
#以下代码会把每个样本的数据创建一个seurat对象，并存放到列表scRNAlist里
for(i in 1:length(dir)){
  counts <- Read10X(data.dir = dir[i])
  scRNAlist[[i]] <- CreateSeuratObject(counts, min.cells=1)
}
#使用merge函数讲10个seurat对象合并成一个seurat对象
scRNA2 <- merge(scRNAlist[[1]], y=c(scRNAlist[[2]], scRNAlist[[3]], 
                                    scRNAlist[[4]], scRNAlist[[5]], scRNAlist[[6]], scRNAlist[[7]], 
                                    scRNAlist[[8]], scRNAlist[[9]], scRNAlist[[10]]))
#dim(scRNA2)
# [1] 23603 19750
table(scRNA2@meta.data$orig.ident)
#HNC01PBMC  HNC01TIL HNC10PBMC  HNC10TIL HNC20PBMC  HNC20TIL     PBMC1     PBMC2   Tonsil1   Tonsil2 
#     1725      1298      1750      1384      1530      1148      2445      2436      3325      2709


scRNA1 <- NormalizeData(scRNA1)
scRNA1 <- FindVariableFeatures(scRNA1, selection.method = "vst")
scRNA1 <- ScaleData(scRNA1, features = VariableFeatures(scRNA1))
scRNA1 <- RunPCA(scRNA1, features = VariableFeatures(scRNA1))
plot1 <- DimPlot(scRNA1, reduction = "pca", group.by="orig.ident")
plot2 <- ElbowPlot(scRNA1, ndims=30, reduction="pca") 
plotc <- plot1+plot2
ggsave("./output/01_mergeSample/cluster1/pca.png", plot = plotc, width = 8, height = 4)
print(c("请选择哪些pc轴用于后续分析？示例如下：","pc.num=1:15"))
#选取主成分
pc.num=1:30

##细胞聚类
scRNA1 <- FindNeighbors(scRNA1, dims = pc.num) 
scRNA1 <- FindClusters(scRNA1, resolution = 0.5)
table(scRNA1@meta.data$seurat_clusters)
metadata <- scRNA1@meta.data
cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)
write.csv(cell_cluster,'./output/01_mergeSample/cluster1/cell_cluster.csv',row.names = F)

##非线性降维
#tSNE
scRNA1 = RunTSNE(scRNA1, dims = pc.num)
embed_tsne <- Embeddings(scRNA1, 'tsne')   #提取tsne图坐标
write.csv(embed_tsne,'./output/01_mergeSample/cluster1/embed_tsne.csv')
#group_by_cluster
plot1 = DimPlot(scRNA1, reduction = "tsne", label=T) 
ggsave("./output/01_mergeSample/cluster1/tSNE.png", plot = plot1, width = 8, height = 7)
#group_by_sample
plot2 = DimPlot(scRNA1, reduction = "tsne", group.by='orig.ident') 
ggsave("./output/01_mergeSample/cluster1/tSNE_sample.png", plot = plot2, width = 8, height = 7)
#combinate
plotc <- plot1+plot2
ggsave("./output/01_mergeSample/cluster1/tSNE_cluster_sample.png", plot = plotc, width = 10, height = 5)

#UMAP
scRNA1 <- RunUMAP(scRNA1, dims = pc.num)
embed_umap <- Embeddings(scRNA1, 'umap')   #提取umap图坐标
write.csv(embed_umap,'./output/01_mergeSample/cluster1/embed_umap.csv') 
#group_by_cluster
plot3 = DimPlot(scRNA1, reduction = "umap", label=T) 
ggsave("./output/01_mergeSample/cluster1/UMAP.png", plot = plot3, width = 8, height = 7)
#group_by_sample
plot4 = DimPlot(scRNA1, reduction = "umap", group.by='orig.ident')
ggsave("./output/01_mergeSample/cluster1/UMAP.png", plot = plot4, width = 8, height = 7)
#combinate
plotc <- plot3+plot4
ggsave("./output/01_mergeSample/cluster1/UMAP_cluster_sample.png", plot = plotc, width = 10, height = 5)

#合并tSNE与UMAP
plotc <- plot2+plot4+ plot_layout(guides = 'collect')
ggsave("./output/01_mergeSample/cluster1/tSNE_UMAP.png", plot = plotc, width = 10, height = 5)

##scRNA2对象的降维聚类参考scRNA1的代码



#------------数据集合并---------------

#scRNAlist是之前代码运行保存好的seurat对象列表，保存了10个样本的独立数据
#数据整合之前要对每个样本的seurat对象进行数据标准化和选择高变基因
for (i in 1:length(scRNAlist)) {
  scRNAlist[[i]] <- NormalizeData(scRNAlist[[i]])
  scRNAlist[[i]] <- FindVariableFeatures(scRNAlist[[i]], selection.method = "vst")
}
##以VariableFeatures为基础寻找锚点，运行时间较长
scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist)
##利用锚点整合数据，运行时间较长
scRNA3 <- IntegrateData(anchorset = scRNA.anchors)
dim(scRNA3)
#[1]  2000 19750    
#有没有发现基因数据只有2000个了？这是因为seurat整合数据时只用2000个高变基因。
#降维聚类的代码省略

#-----------数据质控------------
##==数据质控==#
scRNA <- scRNA3  #以后的分析使用整合的数据进行
##meta.data添加信息
dir.create("./output/01_mergeSample/QC", recursive=T)
proj_name <- data.frame(proj_name=rep("demo2",ncol(scRNA)))
rownames(proj_name) <- row.names(scRNA@meta.data)
scRNA <- AddMetaData(scRNA, proj_name)

##切换数据集
DefaultAssay(scRNA) <- "RNA"

##计算线粒体和红细胞基因比例
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
#计算红细胞比例
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB.genes, rownames(scRNA@assays$RNA)) 
HB.genes <- rownames(scRNA@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
scRNA[["percent.HB"]]<-PercentageFeatureSet(scRNA, features=HB.genes) 
#head(scRNA@meta.data)
col.num <- length(levels(as.factor(scRNA@meta.data$orig.ident)))

##绘制小提琴图
#所有样本一个小提琴图用group.by="proj_name"，每个样本一个小提琴图用group.by="orig.ident"
violin <-VlnPlot(scRNA, group.by = "proj_name",  
                 features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
                 cols =rainbow(col.num), 
                 pt.size = 0.01, #不需要显示点，可以设置pt.size = 0
                 ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
ggsave("./output/01_mergeSample/QC/vlnplot_before_qc.pdf", plot = violin, width = 12, height = 6) 
ggsave("./output/01_mergeSample/QC/vlnplot_before_qc.png", plot = violin, width = 12, height = 6)  
plot1 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.HB")
pearplot <- CombinePlots(plots = list(plot1, plot2, plot3), nrow=1, legend="none") 
ggsave("./output/01_mergeSample/QC/pearplot_before_qc.pdf", plot = pearplot, width = 12, height = 5) 
ggsave("./output/01_mergeSample/QC/pearplot_before_qc.png", plot = pearplot, width = 12, height = 5)

##设置质控标准
print(c("请输入允许基因数和核糖体比例，示例如下：", "minGene=500", "maxGene=4000", "pctMT=20"))
minGene=500
maxGene=3000
pctMT=10

##数据质控
scRNA <- subset(scRNA, subset = nFeature_RNA > minGene & nFeature_RNA < maxGene & percent.mt < pctMT)
col.num <- length(levels(as.factor(scRNA@meta.data$orig.ident)))
violin <-VlnPlot(scRNA, group.by = "proj_name",
                 features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
                 cols =rainbow(col.num), 
                 pt.size = 0.1, 
                 ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
ggsave("./output/01_mergeSample/QC/vlnplot_after_qc.pdf", plot = violin, width = 12, height = 6) 
ggsave("./output/01_mergeSample/QC/vlnplot_after_qc.png", plot = violin, width = 12, height = 6)

#-----------------细胞类型鉴定-------------------
##==鉴定细胞类型==##
library(SingleR)
dir.create("./output/01_mergeSample/CellType")
refdata <- MonacoImmuneData()
testdata <- GetAssayData(scRNA, slot="data")
clusters <- scRNA@meta.data$seurat_clusters
#使用Monaco参考数据库鉴定
cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.main, 
                    method = "cluster", clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")
celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
write.csv(celltype,"./output/01_mergeSample/CellType/celltype_Monaco.csv",row.names = F)
scRNA@meta.data$celltype_Monaco = "NA"
for(i in 1:nrow(celltype)){
  scRNA@meta.data[which(scRNA@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype_Monaco'] <- celltype$celltype[i]}
p1 = DimPlot(scRNA, group.by="celltype_Monaco", repel=T, label=T, label.size=5, reduction='tsne')
p2 = DimPlot(scRNA, group.by="celltype_Monaco", repel=T, label=T, label.size=5, reduction='umap')
p3 = p1+p2+ plot_layout(guides = 'collect')
ggsave("./output/01_mergeSample/CellType/tSNE_celltype_Monaco.png", p1, width=7 ,height=6)
ggsave("./output/01_mergeSample/CellType/UMAP_celltype_Monaco.png", p2, width=7 ,height=6)
ggsave("./output/01_mergeSample/CellType/celltype_Monaco.png", p3, width=10 ,height=5)
#使用DICE参考数据库鉴定
refdata <- DatabaseImmuneCellExpressionData()
# load('~/database/SingleR_ref/ref_DICE_1561s.RData')
# refdata <- ref_DICE
testdata <- GetAssayData(scRNA, slot="data")
clusters <- scRNA@meta.data$seurat_clusters
cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.main, 
                    method = "cluster", clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")
celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
write.csv(celltype,"./output/01_mergeSample/CellType/celltype_DICE.csv",row.names = F)
scRNA@meta.data$celltype_DICE = "NA"
for(i in 1:nrow(celltype)){
  scRNA@meta.data[which(scRNA@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype_DICE'] <- celltype$celltype[i]}
p4 = DimPlot(scRNA, group.by="celltype_DICE", repel=T, label=T, label.size=5, reduction='tsne')
p5 = DimPlot(scRNA, group.by="celltype_DICE", repel=T, label=T, label.size=5, reduction='umap')
p6 = p3+p4+ plot_layout(guides = 'collect')
ggsave("./output/01_mergeSample/CellType/tSNE_celltype_DICE.png", p4, width=7 ,height=6)
ggsave("./output/01_mergeSample/CellType/UMAP_celltype_DICE.png", p5, width=7 ,height=6)
ggsave("./output/01_mergeSample/CellType/celltype_DICE.png", p6, width=10 ,height=5)
#对比两种数据库鉴定的结果
p8 = p1+p4
ggsave("./output/01_mergeSample/CellType/Monaco_DICE.png", p8, width=12 ,height=5)

##保存数据
saveRDS(scRNA,'./output/01_mergeSample/scRNA.rds')