library(Seurat)
library(spacexr)
library(nicheDE)
library(data.table)

Args <- commandArgs(trailingOnly=TRUE)
xR <- Args[1]
print(paste0("NicheDE Downstream Analysis for Region", xR))

# Read in NicheDE Object

NDE <- readRDS(paste0("/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/Niche_DE_Rang20/20240923_NicheDE/20240924_NicheDE_Post_pVal_Region",xR,".rds"))

# Get Niche De Interaction Genes for ASPCs

