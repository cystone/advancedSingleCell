library(Seurat)

library(loomR)
library(SeuratDisk)
inFile = paste0('~/BioFiles/immuneCellAtlas/',
                'cc95ff89-2e68-4a08-a234-480eca21ce79.homo_sapiens.loom')


bm280k.loom <- connect(filename = inFile, mode = "r+")
bm280k.loom
bm280k =  as.Seurat(bm280k.loom)
bm280k.loom$close_all()
# bm280k.data <- Read10X_h5("../data/ica_bone_marrow_h5.h5")
# bm280k <- CreateSeuratObject(counts = bm280k.data, min.cells = 100, min.features = 500)
bm280k.list <- SplitObject(bm280k, split.by = "orig.ident")
bm280k.list <- lapply(X = bm280k.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = bm280k.list)
bm280k.list <- lapply(X = bm280k.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

anchors <- FindIntegrationAnchors(object.list = bm280k.list, reference = c(1, 2), reduction = "rpca", 
                                  dims = 1:50)
bm280k.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)

bm280k.integrated <- ScaleData(bm280k.integrated, verbose = FALSE)
bm280k.integrated <- RunPCA(bm280k.integrated, verbose = FALSE)
bm280k.integrated <- RunUMAP(bm280k.integrated, dims = 1:50)

DimPlot(bm280k.integrated, group.by = "orig.ident")