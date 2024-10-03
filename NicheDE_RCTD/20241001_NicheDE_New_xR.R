library(Seurat)
library(spacexr)
library(nicheDE)
library(data.table)

xR <- 1

Args <- commandArgs(trailingOnly=TRUE)
xR <- Args[1]
print(paste0("NicheDE for Region", xR))

# read in reference dataset

reference=readRDS("/data/cephfs-1/home/users/rauertc_c/work/NicheDE/RCTD/RCTD_Reference_Flex_Wat_Im.rds")
Ref_Raw <- readRDS("/data/cephfs-2/unmirrored/projects/lipomap/flex/flex-wat.im.rds")
#get cell types of reference dataset
cell_types = unique(Ref_Raw$wat.im_pred)
#cell_types = as.character(unique(reference@cell_types))

# read in Seurat Object of Xenium Sample

SeuratObj <- readRDS(paste0("/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/Niche_DE_Rang20/Xenium_Resegments_SeuObj/20240923_R",xR,"_Resegment.rds"))


# read in RCTD Object

RCTD_xR <- readRDS(paste0("/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/Niche_DE_Rang20/RCTD/20240923_Region",xR,"_RCTD.rds"))

###Step 4: Reformat results into a matrix 

# Create Dummy variable Matrix

RCTD_OnlyName_Celltype <- as.data.frame(rownames(RCTD_xR@results$results_df))
RCTD_OnlyName_Celltype$Ident <- RCTD_xR@results$results_df[,2]

for (i in 1:length(cell_types)){
  TypeI <- cell_types[i]
  RCTD_OnlyName_Celltype[,2+i] <- ifelse(RCTD_OnlyName_Celltype$Ident == cell_types[i], 1, 0)
  colnames(RCTD_OnlyName_Celltype)[2+i] <- cell_types[i]
}

DummyMatrix <- RCTD_OnlyName_Celltype[,-c(2)]
colnames(DummyMatrix)[1] <- "Barcode"
rownames(DummyMatrix)<-DummyMatrix$Barcode
DummyMatrix<-DummyMatrix[,-c(1)]
DummyMatrix <- as.matrix(DummyMatrix)

# Create Average Expression Profile Matrix

query.counts<- t(GetAssayData(Ref_Raw, assay = "RNA", slot = "counts"))
query.counts_onlyXen <- query.counts[,which(colnames(query.counts) %in% rownames(SeuratObj))]
predicted.celltype <- as.data.table(Ref_Raw$wat.im_pred, names(Ref_Raw$wat.im_pred))
aveExpL= CreateLibraryMatrix(as.data.frame(query.counts_onlyXen),as.data.frame(predicted.celltype))

# Create NicheDE Object

Counts_Mat <- t(GetAssayData(SeuratObj, assay = "Xenium", slot = "counts"))
Counts_Mat <- as.matrix(Counts_Mat)
Counts_Mat_Cropped <- subset(Counts_Mat, rownames(Counts_Mat) %in% rownames(DummyMatrix))

Coordinates_Mat <- GetTissueCoordinates(SeuratObj)
rownames(Coordinates_Mat)<-Coordinates_Mat[,3]
Coordinates_Mat <- Coordinates_Mat[,-c(3)]
Coordinates_Mat_Cropped <- subset(Coordinates_Mat, rownames(Coordinates_Mat) %in% rownames(DummyMatrix))

AvgDistRectangle <- function(CoordsMat){
  NCells <- nrow(CoordsMat)
  X_Span <- max(CoordsMat$x)-min(CoordsMat$x)
  Y_Span <- max(CoordsMat$y)-min(CoordsMat$y)
  Area <- X_Span*Y_Span
  AreaPerCell <- Area / NCells
  CellsDist <- sqrt(AreaPerCell)
  return(CellsDist)
}

#rm(Ref_Raw)

Kernel <- floor(AvgDistRectangle(Coordinates_Mat_Cropped))

Kernel.lower <- floor(Kernel/10)
Kernel.upper <- floor(Kernel*3)

NicheDE_Obj <- CreateNicheDEObject(Counts_Mat_Cropped, Coordinates_Mat_Cropped, aveExpL, DummyMatrix, c(Kernel.lower,Kernel,Kernel.upper))

NicheDE_Obj = CalculateEffectiveNicheLargeScale(NicheDE_Obj,batch_size = 1000, cutoff = 0.05)

# Perform Niche-DE

NicheDE_Obj = niche_DE(NicheDE_Obj, outfile=paste0("/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/Niche_DE_Rang20/20241001_NicheDE/20240924_NicheDE_Out_Region",xR,".txt"), num_cores = 16, batch = T)

# Get P Values

NicheDE_Obj = get_niche_DE_pval(NicheDE_Obj,pos = T)
NicheDE_Obj = get_niche_DE_pval(NicheDE_Obj,pos = F)

saveRDS(NicheDE_Obj,paste0("/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/Niche_DE_Rang20/20240923_NicheDE/20240924_NicheDE_Post_pVal_Region",xR,".rds"))

print("done")

