cp /data/cephfs-2/unmirrored/projects/lipomap/LPS_Atlas/GSE221493_sc_LPS_atlas.rds /data/cephfs-1/home/users/rauertc_c/scratch/scRNAseq_LPS/temp
cp /data/cephfs-2/unmirrored/projects/lipomap/flex/flex-watss_n8.im.rds /data/cephfs-1/home/users/rauertc_c/scratch/scRNAseq_LPS/temp


library(Seurat)
library(SeuratData)
library(ggplot2)

lipo.ref <- readRDS("/data/cephfs-2/unmirrored/projects/lipomap/flex/flex-watss_n8.im.rds")

lipo.query <- readRDS("/data/cephfs-2/unmirrored/projects/lipomap/LPS_Atlas/GSE221493_sc_LPS_atlas.rds")

lipo.ref <- RunUMAP(lipo.ref, reduction = "integrated.cca", dims = 1:30)
DimPlot(lipo.ref, group.by = c("tech", "celltype"))

