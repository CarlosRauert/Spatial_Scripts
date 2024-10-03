library(Seurat)

# Create Seurat Objects for Resegmented Outs

R1 <- LoadXenium("/data/cephfs-2/unmirrored/projects/lipomap/Xenium_data2/reg1_rang2_ne20/outs/")
R2 <- LoadXenium("/data/cephfs-2/unmirrored/projects/lipomap/Xenium_data2/reg2_rang2_ne20/outs/")
R3 <- LoadXenium("/data/cephfs-2/unmirrored/projects/lipomap/Xenium_data2/reg3_rang2_ne20/outs/")
R4 <- LoadXenium("/data/cephfs-2/unmirrored/projects/lipomap/Xenium_data2/reg4_rang2_ne20/outs/")
R5 <- LoadXenium("/data/cephfs-2/unmirrored/projects/lipomap/Xenium_data2/reg5_rang2_ne20/outs/")
R6 <- LoadXenium("/data/cephfs-2/unmirrored/projects/lipomap/Xenium_data2/reg6_rang2_ne20/outs/")
R7 <- LoadXenium("/data/cephfs-2/unmirrored/projects/lipomap/Xenium_data2/reg7_rang2_ne20/outs/")

# Save Objects as RDS 

saveRDS(R1, "Xenium_Resegments_SeuObj/20240923_R1_Resegment.rds")
saveRDS(R2, "Xenium_Resegments_SeuObj/20240923_R2_Resegment.rds")
saveRDS(R3, "Xenium_Resegments_SeuObj/20240923_R3_Resegment.rds")
saveRDS(R4, "Xenium_Resegments_SeuObj/20240923_R4_Resegment.rds")
saveRDS(R5, "Xenium_Resegments_SeuObj/20240923_R5_Resegment.rds")
saveRDS(R6, "Xenium_Resegments_SeuObj/20240923_R6_Resegment.rds")
saveRDS(R7, "Xenium_Resegments_SeuObj/20240923_R7_Resegment.rds")

# Split Region 7

R7_Unsplit <- R7

R7_new_cells <- read.csv("20242309_Split_R7_Rang20/Region_7_cells_stats.csv", skip=2)
R7_new <- subset(R7_Unsplit, cells=R7_new_cells$Cell.ID)
R7 <- R7_new

R8_new_cells <- read.csv("20242309_Split_R7_Rang20/Region_8_cells_stats.csv", skip=2)
R8 <- subset(R7_Unsplit, cells=R8_new_cells$Cell.ID)

R9_new_cells <- read.csv("20242309_Split_R7_Rang20/Region_9_cells_stats.csv", skip=2)
R9 <- subset(R7_Unsplit, cells=R9_new_cells$Cell.ID)

# Save Split Seurat Objects

saveRDS(R7, "Xenium_Resegments_SeuObj/20240923_R7_Resegment.rds")
saveRDS(R8, "Xenium_Resegments_SeuObj/20240923_R8_Resegment.rds")
saveRDS(R9, "Xenium_Resegments_SeuObj/20240923_R9_Resegment.rds")

