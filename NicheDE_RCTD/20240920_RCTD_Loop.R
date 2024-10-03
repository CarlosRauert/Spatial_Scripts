library(Seurat)
library(future)
plan("multisession", workers = 12)
library(ggplot2)
library(data.table)
library(devtools)
library(spacexr)
#library(pbmcapply)

Args <- commandArgs(trailingOnly=TRUE)
Region <- Args[1]

reference <- readRDS("/data/cephfs-1/home/users/rauertc_c/work/NicheDE/RCTD/RCTD_Reference_Flex_Wat_Im.rds")

RCTD_Region <- function(xR){
  print(xR)
  xenium.obj_region_xR <- readRDS(paste0("/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/Niche_DE_Rang20/Xenium_Resegments_SeuObj/20240923_R",xR,"_Resegment.rds"))
  query.counts <- GetAssayData(xenium.obj_region_xR, assay = "Xenium", slot = "counts")[, Cells(xenium.obj_region_xR[["fov"]])]
  coords <- GetTissueCoordinates(xenium.obj_region_xR[["fov"]], which = "centroids")
  rownames(coords) <- coords$cell
  coords$cell <- NULL
  query <- SpatialRNA(coords, query.counts, colSums(query.counts))
  # run RCTD with many cores
  # hist(colSums(query.counts), breaks = 'FD', xlim = c(0,100))
  RCTD <- create.RCTD(query, reference, max_cores = 12, MAX_MULTI_TYPES = 2, UMI_min = 10, CELL_MIN_INSTANCE = 3)
  RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
  saveRDS(RCTD, file = paste0('/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/Niche_DE_Rang20/RCTD/20240923_Region',xR,'_RCTD.rds'))  # 
  # annotations.df <- RCTD@results$results_df
  # annotations <- annotations.df$first_type
  # names(annotations) <- rownames(annotations.df)
  # xenium.obj_region_xR$predicted.celltype <- annotations
  # saveRDS(xenium.obj_region_xR, file = paste0('/Users/duboisf/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/Liposarcoma/data/Xenium/data/20240603rctd/20240603region',xR,'_preNiche.rds'))
  # # ImageDimPlot(xenium.obj_region_xR, group.by = "predicted.celltype", size = 1.5, cols = "polychrome",
  # #              dark.background = F) + ggtitle("Cell type")
  # keep.cells <- Cells(xenium.obj_region_xR)[!is.na(xenium.obj_region_xR$predicted.celltype)] # tumor cells?
  # xenium.obj_region_xR <- subset(xenium.obj_region_xR, cells = keep.cells)
  # 
  # xenium.obj_region_xR <- BuildNicheAssay(object = xenium.obj_region_xR, fov = "fov", group.by = "predicted.celltype",
  #                               niches.k = 5, neighbors.k = 200)
  # 
  # pdf(file = paste0('/Users/duboisf/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/Liposarcoma/data/Xenium/data/20240603rctd/20240603region',xR,'_preNiche.pdf'), 
  #      paper = "a4r" )
  #   print(ImageDimPlot(xenium.obj_region_xR, group.by = "predicted.celltype", size = 0.2, cols = "polychrome",
  #                                 dark.background = F) + ggtitle("Cell type") )
  #   print(ImageDimPlot(xenium.obj_region_xR, group.by = "niches", size = 0.2, dark.background = F) + ggtitle("Niches") )
  # dev.off()
  # saveRDS(xenium.obj_region_xR, file = paste0('/Users/duboisf/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/Liposarcoma/data/Xenium/data/20240603rctd/20240603region',xR,'_postNiche.rds'))
}

print(paste0("starting RCTD for Region ",Region))
RCTD_Region(Region)
print("RCTD Done")