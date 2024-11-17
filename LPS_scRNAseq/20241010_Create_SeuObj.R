library(Seurat)
expression_matrix <- ReadMtx(
  mtx = "/data/cephfs-2/unmirrored/projects/lipomap/LPS_Atlas/GSE221493_sc_LPS_atlas_raw_counts_matrix.mtx", features = "/data/cephfs-2/unmirrored/projects/lipomap/LPS_Atlas/GSE221493_sc_LPS_atlas_raw_counts_features.tsv",
  cells = "/data/cephfs-2/unmirrored/projects/lipomap/LPS_Atlas/GSE221493_sc_LPS_atlas_raw_counts_barcodes.tsv"
)
seurat_object <- CreateSeuratObject(counts = expression_matrix)
saveRDS(seurat_object,"/data/cephfs-2/unmirrored/projects/lipomap/LPS_Atlas/GSE221493_sc_LPS_atlas.rds")

seuratd_object<-readRDS("/data/cephfs-2/unmirrored/projects/lipomap/LPS_Atlas/GSE221493_sc_LPS_atlas.rds")

#ColSums to determine whether only integers
Cols <- colSums(as.matrix(seurat_object@assays$RNA$counts))
floors <- floor(Cols)-Cols
unique(floors)
# [1] 0

# Visualize QC metrics as a violin plot
pdf(file="/data/cephfs-1/home/users/rauertc_c/work/LPS_Atlas/20241010_Plots/20241010_QC_Plots.pdf")
VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()