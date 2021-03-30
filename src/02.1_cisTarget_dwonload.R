##1, For human:
dbFiles1 <- c("https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather",
              "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather")
##2, For mouse:
dbFiles2 <- c("https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc9nr/gene_based/mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather",
              "https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc9nr/gene_based/mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather")
# mc9nr: Motif collection version 9: 24k motifs

##4, download
dir.create("~/database/cisTarget_databases");   #创建一个文件夹保存数据库
setwd("~/database/cisTarget_databases")
#如果3个参考数据库都想下载，每次设置变量dbFiles后，都要运行以下代码
dbFiles = c(dbFiles1, dbFiles2)
for(featherURL in dbFiles){
  download.file(featherURL, destfile=basename(featherURL)) # saved in current dir
}
# mc9nr: Motif collection version 9: 24k motifs