library(Seurat)
library(tidyverse)
library(patchwork)
library(SCENIC)
rm(list=ls())

BiocManager::install('arrow')

setwd('~/Desktop/advancedSingleCell')
##==分析准备==##
dir.create("./output/02_SCENIC")
dir.create("./output/02_SCENIC/int")
scRNA <- readRDS("./output/01_mergeSample/scRNA.rds")
##准备细胞meta信息
cellInfo <- data.frame(scRNA@meta.data)
colnames(cellInfo)[which(colnames(cellInfo)=="orig.ident")] <- "sample"
colnames(cellInfo)[which(colnames(cellInfo)=="seurat_clusters")] <- "cluster"
colnames(cellInfo)[which(colnames(cellInfo)=="celltype_Monaco")] <- "celltype"
cellInfo <- cellInfo[,c("sample","cluster","celltype")]
saveRDS(cellInfo, file="./output/02_SCENIC/cellInfo.Rds")


##准备表达矩阵
#为了节省计算资源，随机抽取1000个细胞的数据子集
subcell <- sample(colnames(scRNA),1000)
scRNAsub <- scRNA[,subcell]
saveRDS(scRNAsub, "./output/02_SCENIC/scRNAsub.rds")
exprMat <- as.matrix(scRNAsub@assays$RNA@counts)
##设置分析环境
mydbDIR <- "~/database/cisTarget_databases/"
mydbs <- c("hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather",
           "hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather")
names(mydbs) <- c("500bp", "10kb")
scenicOptions <- initializeScenic(org="hgnc", 
                                  nCores=8,
                                  dbDir=mydbDIR, 
                                  dbs = mydbs,
                                  datasetTitle = "HNSCC")
saveRDS(scenicOptions, "./output/02_SCENIC/scenicOptions.rds")

#----------共表达网络计算-------

##==转录调控网络推断==##
##基因过滤
#过滤标准是基因表达量之和>细胞数*3%，且在1%的细胞中表达
genesKept <- geneFiltering(exprMat, scenicOptions, 
                           minCountsPerGene = 3 * 0.01 * ncol(exprMat), 
                           minSamples = ncol(exprMat) * 0.01)
exprMat_filtered <- exprMat[genesKept, ]
##计算相关性矩阵
runCorrelation(exprMat_filtered, scenicOptions)
##TF-Targets相关性回归分析
exprMat_filtered_log <- log2(exprMat_filtered+1)
runGenie3(exprMat_filtered_log, scenicOptions, nParts = 20)
#这一步消耗的计算资源非常大，个人电脑需要几个小时的运行时间
