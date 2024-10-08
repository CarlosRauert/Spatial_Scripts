library(nicheDE)
library(data.table)
library(spacexr)
library(Matrix)
library(Seurat) 

# 6. Perform Niche-DE ####

NDE_obj_m <- readRDS('/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/Niche_DE_Rang20/20241008_NicheDE/20241008NDEobj_m_postEffniche.rds')

NDE_obj_m = niche_DE(NDE_obj_m,num_cores = 32, outfile = "",C = 30, M = 10, gamma = 0.8,print = T, Int = T, batch = T,self_EN = F,G = 1)

saveRDS(NDE_obj, '/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/Niche_DE_Rang20/20241008_NicheDE/20241008NDEobj_m_postNDE.rds')

print("NicheDE Done")