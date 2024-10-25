library(nicheDE)
library(data.table)
library(spacexr)
library(Matrix)
library(Seurat) 

#dir.create('/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/Niche_DE_Rang20/20241017_NicheDE')

# 1. Make Deconvolution Matrices & Average Expression Matrices from RCTD and Seurat Object ####

xR <- 1
MatLst <- lapply(1:9, function(xR){
    print(xR)
    
    rctd_obj <- readRDS(paste0("/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/Niche_DE_Rang20/RCTD/20240923_Region",xR,"_RCTD.rds"))
    Seu_xR <- readRDS(paste0("/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/Niche_DE_Rang20/Xenium_Resegments_SeuObj/20240923_R",xR,"_Resegment.rds"))
    
    # Create Library expression Matrix
    Seu.counts<- GetAssayData(Seu_xR, assay = "Xenium", slot = "counts")[, Cells(Seu_xR[["fov"]])]
    RCTD_Predicts <- as.data.table(rctd_obj@results$results_df, keep.rownames=TRUE)
    RCTD_Predicts <- RCTD_Predicts[, .(rn, first_type)]
    # Exclude Cells that didn't receive RCTD Celltype.
    Seu.counts_Cropped <- Seu.counts[,which(colnames(Seu.counts)%in%RCTD_Predicts$rn)]
    aveExpL= CreateLibraryMatrix(as.data.frame(t(Seu.counts_Cropped)),as.data.frame(RCTD_Predicts))
    aveExpL_sort<-aveExpL[order(rownames(aveExpL)),]

    # Extract the cell type weights (proportions for each spot) and create Dummy Matrix
    weights_matrix <- rctd_obj@results$weights
    dummy_matrix <- apply(weights_matrix, 1, function(x) {
        dummy <- rep(0, length(x))  # Create a vector of zeros
        dummy[which.max(x)] <- 1    # Assign 1 to the index with the maximum proportion
        return(dummy)
    })
    
    # Transpose the dummy matrix to match the original dimensions (spots x cell types) 
    dummy_matrix <- t(dummy_matrix)
    colnames(dummy_matrix) <- colnames(weights_matrix)
    dummy_matrix <- dummy_matrix[, order(colnames(dummy_matrix))]
    return(list(dummy_matrix, aveExpL_sort, RCTD_Predicts))
})

# 2. Create Niche-DE Objects ####

xR <-9
NDE_obj_list <- lapply(c(1,2,3,4,5,6,7,9), function(xR){
    # Read in Seurat Object for Counts / Coords
    Seu_xR <- readRDS(paste0("/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/Niche_DE_Rang20/Xenium_Resegments_SeuObj/20240923_R",xR,"_Resegment.rds"))

    # Get Dummy Mat
    Dummy_NDE <- MatLst[[xR]][[1]]

    # Get Average Expression Matrix
    aveExpL <- MatLst[[xR]][[2]]
    
    # Get Coords
    coords<- GetTissueCoordinates(Seu_xR[["fov"]], which = "centroids")
    coords_NDE <- coords[, 1:2]
    row.names(coords_NDE) <- coords$cell
    coords_NDE <- coords_NDE[rownames(coords_NDE) %in% rownames(Dummy_NDE),]

    # Get Counts
    Seu.counts<- GetAssayData(Seu_xR, assay = "Xenium", slot = "counts")[, Cells(Seu_xR[["fov"]])]
    counts_NDE <- t(Seu.counts)
    counts_NDE <- counts_NDE[rownames(counts_NDE) %in% rownames(Dummy_NDE),]

    #Create NicheDE object
    NDE_obj <- CreateNicheDEObject(counts_NDE,coords_NDE,
                            aveExpL,Dummy_NDE,
                            sigma = c(1,6,9))
})

# 3. merge with other objects ####
NDE_obj_m <- MergeObjects(NDE_obj_list)

# Remove MatLst and NDE_obj_list
#rm(list(MatLst, NDE_obj_list))

# 4. Calculate Effective Niche ####

NDE_obj_m = CalculateEffectiveNicheLargeScale(NDE_obj_m,batch_size = 1000, cutoff = 0.05)
saveRDS(NDE_obj_m, '/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/Niche_DE_Rang20/20241017_NicheDE/20241017NDEobj_m_postEffniche.rds')

# 5. Perform Niche-DE ####

NDE_obj_m <- readRDS('/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/Niche_DE_Rang20/20241017_NicheDE/20241017NDEobj_m_postEffniche.rds')

#NDE_obj_m = niche_DE_no_parallel(NDE_obj_m,C = 30, M = 10, gamma = 0.8,print = T, Int = T, batch = T,self_EN = F)
NDE_obj_m = niche_DE(NDE_obj_m,num_cores = 32, outfile = "",C = 30, M = 10, gamma = 0.8,print = T, Int = T, batch = T,self_EN = F,G = 1)

saveRDS(NDE_obj_m, '/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/Niche_DE_Rang20/20241008_NicheDE/20241017NDEobj_m_postNDE_1_4_7.rds')

DE_genes_T_preA = get_niche_DE_genes(NDE_obj_m,'I',index='ASPC',niche = 't_cell',positive = T,alpha = 0.05)
DE_genes_preA_T = get_niche_DE_genes(NDE_obj_m,'I',index='t_cell',niche = 'ASPC',positive = T,alpha = 0.05)
DE_genes_preA_preA = get_niche_DE_genes(NDE_obj_m,'CT',index='ASPC',niche = 'ASPC',positive = T,alpha = 0.05)
DE_genes_T_B = get_niche_DE_genes(NDE_obj_m,'I',index='b_cell',niche = 't_cell',positive = T,alpha = 0.05)
DE_genes_B_T = get_niche_DE_genes(NDE_obj_m,'I',index='b_cell',niche = 't_cell',positive = T,alpha = 0.05)