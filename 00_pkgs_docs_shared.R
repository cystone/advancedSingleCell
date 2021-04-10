suppressMessages(library(tidyverse))
suppressMessages(library(pacman))
suppressMessages(library(data.table))

wkPath <- c('./result', './processData')
for(i in wkPath){
  wkPathi = i
  # wkPathi = paste0(sectionName, '/', i)
  #每一个子项目都含plot、result、input
  if (!dir.exists(wkPathi)) dir.create(wkPathi, recursive=T)
}
rm(list=c('i', 'wkPathi', 'wkPath'))

# install.packages('umap')
# BiocManager::install("glmGamPoi")
# remotes::install_github('satijalab/seurat-data')
