library(nicheDE)
library(data.table)
library(spacexr)
library(Matrix)
library(Seurat) 

# 1. Get Avg Expression from Reference ####

# read Flex_Wat_IM

Flex_Wat_IM <- readRDS("/data/cephfs-2/unmirrored/projects/lipomap/flex/flex-wat.im.rds")

# make sure Idents are Cell_types
# determine Assay

Idents(Flex_Wat_IM) <- Flex_Wat_IM@meta.data$wat.im_pred

# Create Avg Expression Matrix:

query.counts<- GetAssayData(Flex_Wat_IM, assay = "RNA", slot = "counts")
  
Flex_Wat_IM_predicted.celltype <- as.data.table(names(Flex_Wat_IM$wat.im_pred ))
Flex_Wat_IM_predicted.celltype[, celltype:= Flex_Wat_IM$wat.im_pred]
aveExpL= CreateLibraryMatrix(as.data.frame(t(query.counts)),as.data.frame(Flex_Wat_IM_predicted.celltype))

# 2. Make Deconvolution Matrix from RCTD ####

xR <- 1
dumLst <- lapply(1:9, function(xR){
    rctd_obj <- readRDS(paste0("/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/Niche_DE_Rang20/RCTD/20240923_Region",xR,"_RCTD.rds"))
    # Extract the cell type weights (proportions for each spot)
    weights_matrix <- rctd_obj@results$weights
    dummy_matrix <- apply(weights_matrix, 1, function(x) {
        dummy <- rep(0, length(x))  # Create a vector of zeros
        dummy[which.max(x)] <- 1    # Assign 1 to the index with the maximum proportion
        return(dummy)
    })
    # Transpose the dummy matrix to match the original dimensions (spots x cell types) 
    dummy_matrix <- t(dummy_matrix)
    colnames(dummy_matrix) <- colnames(weights_matrix)
    CT_in <- rownames(aveExpL)[rownames(aveExpL)%in%colnames(dummy_matrix)]
    dummy_matrix <- dummy_matrix[, CT_in]
    return(dummy_matrix)
})

# 3. Create Niche-DE Objects ####

xR <-1
NDE_obj_list <- lapply(c(1,2,3,4,5,6,7,9), function(xR){
    # Read in Seurat Object for Counts & Coords
    Seu_xR <- readRDS(paste0("/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/Niche_DE_Rang20/Xenium_Resegments_SeuObj/20240923_R",xR,"_Resegment.rds"))
    
    # Get Coords
    coords<- GetTissueCoordinates(Seu_xR[["fov"]], which = "centroids")
    coords_NDE <- coords[, 1:2]
    row.names(coords_NDE) <- coords$cell

    # Get Dummy Mat
    Dummy_NDE <- dumLst[[xR]]
    Dummy_NDE <- Dummy_NDE[rownames(Dummy_NDE) %in% rownames(coords_NDE),]

    coords_NDE <- coords_NDE[rownames(coords_NDE) %in% rownames(Dummy_NDE),]

    # Get Counts
    Seu.counts<- GetAssayData(Seu_xR, assay = "Xenium", slot = "counts")[, Cells(Seu_xR[["fov"]])]
    counts_NDE <- t(Seu.counts)
    counts_NDE <- counts_NDE[rownames(counts_NDE) %in% rownames(Dummy_NDE),]
    
    # Crop aveExpL
    aveExpL_crop <- aveExpL[rownames(aveExpL) %in% colnames(Dummy_NDE),]

    #Create NicheDE object
    NDE_obj <- CreateNicheDEObject(counts_NDE,coords_NDE,
                            aveExpL_crop,Dummy_NDE,
                            sigma = c(1,100,250))
})

# Only Region 3

NDE_obj_R3 <- NDE_obj_list[[3]]
NDE_obj_R3 <- CalculateEffectiveNicheLargeScale(NDE_obj_R3,batch_size = 1000, cutoff = 0.05)
saveRDS(NDE_obj_R3, '/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/Niche_DE_Rang20/20241007_NicheDE/20241007NDEobj_R3_postEffniche.rds')
NDE_obj_R3 = niche_DE(NDE_obj_R3,num_cores = 16, outfile = "",C = 30, M = 10, gamma = 0.8,print = T, Int = T, batch = T,self_EN = F,G = 1)
saveRDS(NDE_obj_R3, '/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/Niche_DE_Rang20/20241007_NicheDE/20241007NDEobj_R3_postNDE.rds')