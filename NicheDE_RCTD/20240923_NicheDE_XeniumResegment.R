library(Seurat)
library(spacexr)
library(nicheDE)
library(data.table)

# read in reference dataset

xR <-1

reference=readRDS("/data/cephfs-1/home/users/rauertc_c/work/NicheDE/RCTD/RCTD_Reference_Flex_Wat_Im.rds")
Ref_Raw <- readRDS("/data/cephfs-2/unmirrored/projects/lipomap/flex/flex-wat.im.rds")
#get cell types of reference dataset
cell_types = unique(Ref_Raw$wat.im_pred)

# read in Seurat Object of Xenium Sample

SeuratObj <- readRDS("/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/Niche_DE_Rang20/Xenium_Resegments_SeuObj/20240923_R1_Resegment.rds")


# read in RCTD Object

RCTD_R1 <- readRDS("RCTD/20240923_Region1_RCTD.rds")

###Step 4: Reformat results into a matrix 

# Create Dummy variable Matrix

RCTD_R1_OnlyName_Celltype <- as.data.frame(rownames(RCTD_R1@results$results_df))
RCTD_R1_OnlyName_Celltype$Ident <- RCTD_R1@results$results_df[,2]

for (i in 1:length(cell_types)){
  TypeI <- cell_types[i]
  RCTD_R1_OnlyName_Celltype[,2+i] <- ifelse(RCTD_R1_OnlyName_Celltype$Ident == cell_types[i], 1, 0)
  colnames(RCTD_R1_OnlyName_Celltype)[2+i] <- cell_types[i]
}

DummyMatrix <- RCTD_R1_OnlyName_Celltype[,-c(2)]
colnames(DummyMatrix)[1] <- "Barcode"
rownames(DummyMatrix)<-DummyMatrix$Barcode
DummyMatrix<-DummyMatrix[,-c(1)]
DummyMatrix <- as.matrix(DummyMatrix)

#saveRDS(DummyMatrix, paste0("/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/Niche_DE_Rang20/20240923_NicheDE/Dummy_Matrixes/Reg",xR,"DummyMat.rds"))

# Create Average Expression Profile Matrix

query.counts<- t(GetAssayData(Ref_Raw, assay = "RNA", slot = "counts"))
predicted.celltype <- as.data.table(Ref_Raw$wat.im_pred, names(Ref_Raw$wat.im_pred))
aveExpL= CreateLibraryMatrix(as.data.frame(query.counts),as.data.frame(predicted.celltype))

# Create NicheDE Object

Counts_Mat <- t(GetAssayData(SeuratObj, assay = "Xenium", slot = "counts"))
Counts_Mat <- as.matrix(Counts_Mat)
Counts_Mat_Cropped <- subset(Counts_Mat, rownames(Counts_Mat) %in% rownames(DummyMatrix))

Coordinates_Mat <- GetTissueCoordinates(SeuratObj)
rownames(Coordinates_Mat)<-Coordinates_Mat[,3]
Coordinates_Mat <- Coordinates_Mat[,-c(3)]
Coordinates_Mat_Cropped <- subset(Coordinates_Mat, rownames(Coordinates_Mat) %in% rownames(DummyMatrix))

# Estimate Kernel Size

AvgDistRectangle <- function(CoordsMat){
  NCells <- nrow(CoordsMat)
  X_Span <- max(CoordsMat$x)-min(CoordsMat$x)
  Y_Span <- max(CoordsMat$y)-min(CoordsMat$y)
  Area <- X_Span*Y_Span
  AreaPerCell <- Area / NCells
  CellsDist <- sqrt(AreaPerCell)
  return(CellsDist)
}

rm(Ref_Raw)

Kernel <- floor(AvgDistRectangle(Coordinates_Mat_Cropped))

Kernel.lower <- floor(Kernel/10)
Kernel.upper <- floor(Kernel*3)

NicheDE_Obj <- CreateNicheDEObject(Counts_Mat_Cropped, Coordinates_Mat_Cropped, aveExpL, DummyMatrix, c(Kernel.lower,Kernel,Kernel.upper))

# Calculate Effective Niche

NicheDE_Obj = CalculateEffectiveNicheLargeScale(NicheDE_Obj,batch_size = 1000, cutoff = 0.05)

# Perform Niche-DE

NicheDE_Obj = niche_DE(NicheDE_Obj, outfile=paste0("/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/Niche_DE_Rang20/20240923_NicheDE/20240924_NicheDE_PostNicheDE_Region",xR,".rds")

