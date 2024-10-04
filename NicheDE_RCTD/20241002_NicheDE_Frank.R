library(nicheDE)
library(data.table)
library(spacexr)
library(Matrix)
library(Seurat) 

xR = 2
NDE_obj_L <- lapply(1:4, function(xR){
  print(xR)
  postnicheSeuratXenium<- readRDS(paste0('/Users/duboisf/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/Liposarcoma/data/Xenium/data/20240603rctd/20240603region',xR,'_postNiche.rds'))
  Idents(postnicheSeuratXenium) <- postnicheSeuratXenium$predicted.celltype 
  #get average expression profile matrix 
  #create library matrix
  query.counts<- GetAssayData(postnicheSeuratXenium, assay = "Xenium", slot = "counts")[, Cells(postnicheSeuratXenium[["fov"]])]
  
  postnicheSeuratXenium_predicted.celltype <- as.data.table(names(postnicheSeuratXenium$predicted.celltype ))
  postnicheSeuratXenium_predicted.celltype[, celltype:= postnicheSeuratXenium$predicted.celltype]
  aveExpL= CreateLibraryMatrix(t(query.counts),postnicheSeuratXenium_predicted.celltype)
  View(L)
  
  # deconv_mat
  RCTD <- readRDS(paste0('/Users/duboisf/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/Liposarcoma/data/Xenium/data/20240603rctd/20240603region',xR,'_RCTD.rds'))
  RCTD_results <- RCTD@results[["results_df"]]
  RCTD_weights <- RCTD@results[["weights"]] # = deconvMat
  
  #dummy deconv_mat
  postnicheSeuratXenium_predicted.celltype[, val := 1]
  colnames(postnicheSeuratXenium_predicted.celltype)[1] = 'ID'
  dummyDeconv<- dcast(postnicheSeuratXenium_predicted.celltype, formula = ID ~ celltype, fill = 0, value.var = 'val')
  dummyDeconv_mat <- as.matrix(dummyDeconv[, -1])
  rownames(dummyDeconv_mat) <- dummyDeconv$ID
  #make niche-DE object
  coords<- GetTissueCoordinates(postnicheSeuratXenium[["fov"]], which = "centroids")
  coords_NDE <- coords[, 1:2]
  row.names(coords_NDE) <- coords$cell
  
  idx<- sapply(colnames(RCTD_weights), function(xr){
    which(rownames(aveExpL)  == xr)
  })
  aveExpL_sort <- aveExpL[idx, ]
  
  counts_NDE <- t(query.counts)
  RCTD_weights_mat <- as.matrix(RCTD_weights)
  # NDE_obj= CreateNicheDEObject(counts_NDE,coords_NDE,
  #                              aveExpL_sort,dummyDeconv_mat,
  #                              sigma = c(1,100,250))
  
  NDE_obj= CreateNicheDEObject(counts_NDE,coords_NDE,
                                   aveExpL_sort,RCTD_weights_mat,
                                sigma = c(1,100,250))
  return(NDE_obj)
})


# merge with other objects
NDE_obj_m <- MergeObjects(NDE_obj_L)

# EffectiveNiche
NDE_obj = CalculateEffectiveNicheLargeScale(NDE_obj_m,batch_size = 1000, cutoff = 0.05)
saveRDS(NDE_obj, '/Users/duboisf/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/Liposarcoma/data/Xenium/data/20240617nicheDE/20240617NDEobj_postEffniche.rds')

#Perform Niche-DE
NDE_obj = niche_DE(NDE_obj,num_cores = 10, outfile = "",C = 30, M = 10, gamma = 0.8,print = T, Int = T, batch = T,self_EN = F,G = 1)
saveRDS(NDE_obj, '/Users/duboisf/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/Liposarcoma/data/Xenium/data/20240617nicheDE/20240617NDEobj_postNDE.rds')

res_cell_type_level <- NDE_obj@niche_DE_pval_pos[["cell_type_level"]]
# Remove rows where all elements are NA using rowSums
res_cell_type_level_clean <- res_cell_type_level[rowSums(is.na(res_cell_type_level)) != ncol(res_cell_type_level), ]
res_cell_type_level_clean_idx <- res_cell_type_level_clean < 0.05
res_cell_type_level_clean_sig <- res_cell_type_level_clean[rowSums(res_cell_type_level_clean_idx, na.rm = T) > 0, ]
idxC <- colnames(res_cell_type_level_clean_sig) %in% c("B cell", "dendritic cell" , "fat cell", "immature NK T cell", "macrophage" , "monocyte", "preadipocyte",  "T cell")
res_cell_type_level_clean_idx <- res_cell_tydpe_level_clean < 1e-12
res_cell_type_level_clean_sig2 <- res_cell_type_level_clean[rowSums(res_cell_type_level_clean_idx[, idxC], na.rm = T) > 0, idxC]
write.csv2(res_cell_type_level_clean_sig2, '/Users/duboisf/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/Liposarcoma/data/Xenium/data/20240617nicheDE/20240617spatialDEgenes_by_cell.csv', row.names = F)

# res_cell_type_Interaxtion_level <- NDE_obj@niche_DE_pval_pos[["interaction_level"]]
# res_cell_type_Interaxtion_level_predad <- res_cell_type_Interaxtion_level[13, , ]
# # clean_matrix <- res_cell_type_Interaxtion_level[!apply(res_cell_type_Interaxtion_level, 1, function(row) all(is.na(row))), ,]
DE_genes_T_preA = get_niche_DE_genes(NDE_obj,'I',index='preadipocyte',niche = 'T cell',positive = T,alpha = 0.05)
DE_genes_preA_T = get_niche_DE_genes(NDE_obj,'I',index='T cell',niche = 'preadipocyte',positive = T,alpha = 0.05)
DE_genes_preA_preA = get_niche_DE_genes(NDE_obj,'CT',index='preadipocyte',niche = 'preadipocyte',positive = T,alpha = 0.05)
DE_genes_T_B = get_niche_DE_genes(NDE_obj,'I',index='B cell',niche = 'T cell',positive = T,alpha = 0.05)
DE_genes_B_T = get_niche_DE_genes(NDE_obj,'I',index='T cell',niche = 'B cell',positive = T,alpha = 0.05)
# #Load enichr package
library(enrichR)
# #run pathway enrichment analysis
tum_processes = enrichr(DE_genes_preA_preA[,1],databases = 'Reactome_2022')
tum_processes_dt <- as.data.table(tum_processes$Reactome_2022)
tum_processes_dt_sig <- tum_processes_dt[Adjusted.P.value < 0.05]
# #run pathway enrichment analysis
T_preA_processes = enrichr(DE_genes_T_preA[,1],databases = 'Reactome_2022')
T_preA_processes_dt <- as.data.table(T_preA_processes$Reactome_2022)
T_preA_processes_dt_sig <- T_preA_processes_dt[Adjusted.P.value < 0.05]
write.csv2(T_preA_processes_dt_sig, '/Users/duboisf/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/Liposarcoma/data/Xenium/data/20240617nicheDE/20240617T_preA_interactions.csv', row.names = F)


preA_T_processes = enrichr(DE_genes_preA_T[,1],databases = 'Reactome_2022')
preA_T_processes_dt <- as.data.table(preA_T_processes$Reactome_2022)
preA_T_processes_dt_sig <- preA_T_processes_dt[Adjusted.P.value < 0.05]
write.csv2(preA_T_processes_dt_sig, '/Users/duboisf/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/Liposarcoma/data/Xenium/data/20240617nicheDE/20240617pread_t_interactions.csv', row.names = F)

#Ligand Receptor Analysis
data("niche_net_ligand_target_matrix")
data("ramilowski_ligand_receptor_list")
Tcell_tumor_LR = niche_LR_cell(NDE_obj,ligand_cell = 'T cell',receptor_cell = 'preadipocyte',
                               ligand_target_matrix = niche_net_ligand_target_matrix,
                               lr_mat = ramilowski_ligand_receptor_list,K = 25,M = 50,alpha = 0.05,truncation_value = 3)
# no ligand-receptor pairs to report

