library(Seurat)
#library(SeuratDisk)

#### Define args
args = commandArgs(trailingOnly=TRUE)

### Load 10x data
illumina <- Read10X(args[1])
illumina <- CreateSeuratObject(counts=illumina, project='illumina', min.cells=3)
raw_counts = as.data.frame(illumina@assays$RNA@counts)
write.csv(raw_counts, args[2])
