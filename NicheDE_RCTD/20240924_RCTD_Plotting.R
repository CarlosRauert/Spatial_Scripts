library(spacexr)
library(Matrix)

for (i in 1:9){
  #Read in RCTD Object
  RCTD_xR <- readRDS(paste0("/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/Niche_DE_Rang20/RCTD/20240923_Region",i,"_RCTD.rds"))
  cell_type_names <- RCTD_xR@cell_type_info$info[[2]]
  dir.create(paste0("/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/Niche_DE_Rang20/RCTD/20240924_RCTD_Plots_Region",i))
  plot_all_cell_types(RCTD_xR@results$results_df, RCTD_xR@spatialRNA@coords, cell_type_names, resultsdir = paste0("/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/Niche_DE_Rang20/RCTD/20240924_RCTD_Plots_Region",i)) 
}