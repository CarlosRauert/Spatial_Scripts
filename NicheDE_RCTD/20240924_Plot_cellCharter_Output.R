library(data.table)
library(Seurat)
library(qvalue)
library(ggplot2)
library(ggrepel)
library(parallel)
library(poolr)
library(dplyr)

adata.obs <- fread("/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20240923cellcharterXenium_auto_k2_k20/20240924_cc_adata_LS_Xenium_k11clstrd_adata.obs_.csv")
countMtx <- fread('/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20240923cellcharterXenium_auto_k2_k20/20240924_cc_adata_LS_Xenium_k11clstrd_CountMat_.csv')
genes_dt <- fread('/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20240923cellcharterXenium_auto_k2_k20/20240924_cc_LS_Xenium_k11clstrd_adata.var_.csv')

#countMtx[, cell_ID := adata.obs$cell_id]
colnames(countMtx) = genes_dt$V1
transcripts_per_cell <- rowSums(countMtx)
countMtx_fdt <- countMtx[-1, ]
countMtx_fdt[, cell_ID := adata.obs$cell_id]

#Region1_WDLS <- adata.obs[(adata.obs$sample == 1) & (adata.obs$cluster_cellcharter != 3), ]
#Region2_WDLS <- adata.obs[(adata.obs$sample == 2) & (adata.obs$cluster_cellcharter == 7), ]
#Region3_WDLS <- adata.obs[(adata.obs$sample == 3) & (adata.obs$cluster_cellcharter == 7), ]
#Region4_WDLS <- adata.obs[(adata.obs$sample == 4) & (adata.obs$cluster_cellcharter == 7), ]
#Region5_WDLS <- adata.obs[(adata.obs$sample == 5) & (adata.obs$cluster_cellcharter == 5), ]
#Region6_WDLS <- adata.obs[(adata.obs$sample == 6) & (adata.obs$cluster_cellcharter == 7), ]
#Region7_WDLS <- adata.obs[(adata.obs$sample == 7) & (adata.obs$cluster_cellcharter == 5), ]
#Region8_WDLS <- adata.obs[(adata.obs$sample == 8) & (adata.obs$cluster_cellcharter == 5), ]
#Region9_WDLS <- adata.obs[(adata.obs$sample == 9) & (adata.obs$cluster_cellcharter == 5), ]

#WDLS <- rbind(Region1_WDLS, Region2_WDLS, Region2_WDLS, Region4_WDLS, Region5_WDLS, Region6_WDLS, Region7_WDLS, Region8_WDLS, Region9_WDLS)

#Region2_DDLS <- adata.obs[(adata.obs$sample == 2) & (adata.obs$cluster_cellcharter ==1), ]
#Region3_DDLS <- adata.obs[(adata.obs$sample == 3) & (adata.obs$cluster_cellcharter ==2), ]
#Region4_DDLS <- adata.obs[(adata.obs$sample == 4) & (adata.obs$cluster_cellcharter ==4), ]
#Region5_DDLS <- adata.obs[(adata.obs$sample == 5) & (adata.obs$cluster_cellcharter >= 8), ]
#Region6_DDLS <- adata.obs[(adata.obs$sample == 6) & (adata.obs$cluster_cellcharter ==8), ]
#Region7_DDLS <- adata.obs[(adata.obs$sample == 7) & (adata.obs$cluster_cellcharter ==6), ]
#Region8_DDLS <- adata.obs[(adata.obs$sample == 8) & (adata.obs$cluster_cellcharter ==7), ]
#Region9_DDLS <- adata.obs[(adata.obs$sample == 9) & (adata.obs$cluster_cellcharter <= 1), ]

#DDLS <- rbind(Region2_DDLS, Region3_DDLS, Region4_DDLS, Region5_DDLS, Region6_DDLS, Region7_DDLS, Region8_DDLS, Region9_DDLS)

filter(adata.obs, )

#region loop ####

#1?


#adata.obs$subtype <- mclapply(1:length(adata.obs$cell_id),function(xCell){ 
#  if (adata.obs$sample[xCell]==2){
#    if (adata.obs$cluster_cellcharter[xCell] %in% c(1)){
#      return("DDLS")
#    }
#    else if ((adata.obs$cluster_cellcharter[xCell] %in% c(7,5))){
#      return("WDLS")
#    }
#    else {
#      return("both / healthy")
#    }
#  }
#  if (adata.obs$sample[xCell]==3){
#    if (adata.obs$cluster_cellcharter[xCell] %in% c(0,1,2)){
#      return("DDLS")
#    }
#    else if ((adata.obs$cluster_cellcharter[xCell] %in% c(7,5))){
#      return("WDLS")
#    }
#    else {
#      return("both / healthy")
#    }
#  }
#  if (adata.obs$sample[xCell]==4){
#    if (adata.obs$cluster_cellcharter[xCell] %in% c(4)){
#      return("DDLS")
#    }
#    else if ((adata.obs$cluster_cellcharter[xCell] %in% c(7,5))){
#      return("WDLS")
#    }
#    else {
#      return("both / healthy")
#    }
#  }    
#  if (adata.obs$sample[xCell]==5){
#    if (adata.obs$cluster_cellcharter[xCell] %in% c(0,1,3,8,10)){
#      return("DDLS")
#    }
#    else if ((adata.obs$cluster_cellcharter[xCell] %in% c(5))){
#      return("WDLS")
#    }
#    else {
#      return("both / healthy")
#    }
#  }
#  if (adata.obs$sample[xCell]==6){
#    if (adata.obs$cluster_cellcharter[xCell] %in% c(4,9)){
#      return("DDLS")
#    }
#    else if ((adata.obs$cluster_cellcharter[xCell] %in% c(5,7))){
#      return("WDLS")
#    }
#    else {
#      return("both / healthy")
#    }
#  }
#  if (adata.obs$sample[xCell]==7){
#    if (adata.obs$cluster_cellcharter[xCell] %in% c(6)){
#      return("DDLS")
#    }
#    else if ((adata.obs$cluster_cellcharter[xCell] %in% c(5,7))){
#      return("WDLS")
#    }
#    else {
#      return("both / healthy")
#    }
#  }
#  if (adata.obs$sample[xCell]==8){
#    if (adata.obs$cluster_cellcharter[xCell] %in% c(7,0)){
#      return("DDLS")
#    }
#    else if ((adata.obs$cluster_cellcharter[xCell] %in% c(5))){
#      return("WDLS")
#    }
#    else {
#      return("both / healthy")
#    }
#  }
#  if (adata.obs$sample[xCell]==9){
#    if (adata.obs$cluster_cellcharter[xCell] %in% c(2,0,1)){
#      return("DDLS")
#    }
#    else if ((adata.obs$cluster_cellcharter[xCell] %in% c(5,7))){
#      return("WDLS")
#    }
#    else {
#      return("both / healthy")
#    }
#  }
 # }, mc.cores=10)

#unique(adata.obs$subtype)

resall_dt <- rbindlist(lapply(2:9, function(xR){
  print(xR)
  adata.obs_xR <- adata.obs[sample == xR]
  #adata.obs_xR[grepl(pattern = 'DDLS', x = spatial_cluster), subtype := 'DDLS']
  #adata.obs_xR[grepl(pattern = 'WDLS', x = spatial_cluster), subtype := 'WDLS']
  adata.obs_xRf <- adata.obs_xR[ !is.na(subtype)]
  #loop over celltypes
  xCT = unique(adata.obs_xRf$cell_types)[1]
  fresCTs <- rbindlist(mclapply(unique(adata.obs_xRf$cell_types), function(xCT){
    print(xCT)
    confmatCTx = matrix(c(nrow(adata.obs_xRf[subtype == 'DDLS' & cell_types == xCT]), nrow(adata.obs_xRf[subtype == 'DDLS' & cell_types != xCT]), 
                          nrow(adata.obs_xRf[subtype == 'WDLS' & cell_types == xCT]), nrow(adata.obs_xRf[subtype == 'WDLS' & cell_types != xCT])), 
                        nrow = 2, ncol = 2, byrow = T)
    colnames(confmatCTx) = c('xCT', 'n_xCT')
    rownames(confmatCTx) = c('DDLS', 'WDLS')
    resxG <- fisher.test(confmatCTx)
    resdt <- as.data.table(xCT)
    resdt[, DDLS_share := nrow(adata.obs_xRf[subtype == 'DDLS' & cell_types == xCT]) /nrow(adata.obs_xRf[subtype == 'DDLS'])]
    resdt[, WDLS_share := nrow(adata.obs_xRf[subtype == 'WDLS' & cell_types == xCT]) /nrow(adata.obs_xRf[subtype == 'WDLS'])]
    resdt[, p_val:= resxG$p.value]
    resdt[, OR := resxG$estimate]
    resdt[, CI_l := resxG$conf.int[1]]
    resdt[, CI_H := resxG$conf.int[2]]
    return(resdt)
  }, mc.cores = 10))
  fresCTs[p_val < 0, p_val :=0]
  fresCTs[p_val > 1, p_val :=1]
  fresCTs[, q_val := qvalue(p_val, pi0 = 1)$qvalue]
  # write.csv2(WXres, paste0('/Users/duboisf/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/Liposarcoma/data/Xenium/data/20240810cc_nicheDiffexp/20240811_',xR,'DDLSvsWDLS', '.csv'),quote=FALSE, row.names=FALSE, col.names = TRUE)
  fresCTs[, region := xR]
  fresCTs[,log2FC := log2(DDLS_share/ WDLS_share), by = DDLS_share]
  fresCTs[,neglog10q := log10(q_val)*-1]
  fresCTs[q_val < 0.05 & (log2FC < 0) ,colr := 'blue']
  fresCTs[q_val < 0.05 & (log2FC > 0) ,colr := 'red']
  fresCTs[,labelid := xCT]
  #lower_x <- min(!is.na(fresCTs$log2FC))-1
  #upper_x <- max(!is.na(fresCTs$log2FC))+1
  #lower_y <- min(!is.na(fresCTs$neglog10q))-1 + xlim(lower_x, upper_x) + ylim(lower_y,upper_y))
  #upper_y <- max(!is.na(fresCTs$neglog10q))+1
  pdf(file = paste0('/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20240923cellcharterXenium_auto_k2_k20/Diff_WDvsDD/Region',xR,'_DDLSvsWDLS.pdf'), width = 18, height = 11)
  print(ggplot(fresCTs, aes(x=log2FC, y=neglog10q)) + geom_point(aes(col=colr, size = DDLS_share)) + geom_text_repel(aes(label=labelid)) + theme_classic() + theme(text = element_text(size=24)) + expand_limits(x = c(min(fresCTs$log2FC)-1,max(fresCTs$log2FC)+1), y=(c(max(fresCTs$neglog10q)+10))))    # force=3 
  dev.off()
  return(fresCTs)
}))

MeansDT <- function(DT){
  DF <- as.data.frame(DT)[2:11,]
  DF[!is.finite(m)] <- NA
  colMeans(m, na.rm=TRUE)
  #Means_List <- append(Means_List, DT[1,12])
  #if (Means_List[8] < 0.05 & Means_List[10]<0){
  #  Means_List <- append(Means_List,"blue")
  #}
  #else if (Means_List[8] < 0.05 & Means_List[10]>0){
  #  Means_List <- append(Means_List,"blue")
  #}
  #else{
  #  Means_List <- append(Means_List,NA)
  #}
}

MeansDF <- data.frame(xCT=rep("", 1), DDLS_share=rep(0, 1), WDLS_share=rep(0,1), p_val=rep(0,1), OR=rep(0,1), CI_l=rep(0,1), CI_H=rep(0,1), q_val=rep(0,1), region=rep(0,1), log2FC=rep(0,1), neglog10q=rep(0,1), colr=rep("", 1), labelid=rep("",1), stringsAsFactors=FALSE)

for (i in 1:length(unique(resall_dt$xCT))){
  MeansDF <- rbind(MeansDF, MeansList(resall_dt[,unique(resall_dt$xCT)[i]]))
}


print(MeansList(resall_dt[resall_dt$xCT=="macrophage",]))

write.csv2(resall_dt, file = '/Users/duboisf/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/Liposarcoma/data/Xenium/data/20240811nicheDiffCelltypes_DWDLS/20240811regionsAllresall_dt_int_DDLSvsWDLS.csv'
           , row.names = F)
# Fisher's methods for consitent changes #####
xCTm= unique(resall_dt$xCT)[1]
res_merged <- rbindlist(lapply(unique(resall_dt$xCT), function(xCTm) {
  print(xCTm)
  xmGresTbl <- resall_dt[xCT == xCTm]
  if ((sum(xmGresTbl$OR > 1) == 7) | (sum(xmGresTbl$OR < 1) == 7) ) {
    merged_res_xG <- xmGresTbl[1, 1]
    merged_res_xG[, meanShare_DDLS := mean(xmGresTbl$DDLS_share)]
    merged_res_xG[, meanShare_WDLS := mean(xmGresTbl$WDLS_share)]
    merged_res_xG[, merged_pval := fisher(p = xmGresTbl$p_val)$p]
    merged_res_xG[, median_OR := median(xmGresTbl$OR)]
    merged_res_xG[, mean_log2FC := log2(meanShare_DDLS/meanShare_WDLS)]
    return(merged_res_xG)
  } else  {
    return(data.table())
  }
}))
res_merged[merged_pval < 0, merged_pval :=0]
res_merged[merged_pval > 1, merged_pval :=1]
res_merged[, q_val := qvalue(merged_pval, pi0 = 1)$qvalue]
write.csv2(res_merged, file = '/Users/duboisf/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/Liposarcoma/data/Xenium/data/20240811nicheDiffCelltypes_DWDLS/20240811regions_mergedFisher_DDLSvsWDLS.csv'
           , row.names = F)
# pseudobulk####
bulktests <- rbindlist(lapply(unique(resall_dt$xCT), function(xmCT) {
  print(xmCT)
  xmGresTbl <- resall_dt[xCT == xmCT]
  resxG <- wilcox.test(xmGresTbl$DDLS_share, xmGresTbl$WDLS_share, exact = T, conf.int = T, paired = T)
  resdt <- as.data.table(xmCT)
  resdt[, meanShare_Niche_DDLS := mean(xmGresTbl$DDLS_share)]
  resdt[, meanShare_Niche_WDLS := mean(xmGresTbl$WDLS_share)]
  resdt[, sdShare_Niche_DDLS := sd(xmGresTbl$DDLS_share)]
  resdt[, sdShare_Niche_WDLS := sd(xmGresTbl$DDLS_share)]
  resdt[, p_val:= resxG$p.value]
  resdt[, estimate:= resxG$estimate]
  resdt[, conf_est_l:= resxG$conf.int[1]]
  resdt[, conf_est_h:= resxG$conf.int[2]]
  return(resdt)
}))
bulktests[p_val < 0, p_val :=0]
bulktests[p_val > 1, p_val :=1]
bulktests[, q_val := qvalue(p_val, pi0 = 1)$qvalue]
write.csv2(bulktests, file = '/Users/duboisf/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/Liposarcoma/data/Xenium/data/20240811nicheDiffCelltypes_DWDLS/20240811regionsAll_bulktests_DDLSvsWDLS.csv'
           , row.names = F)
bulktests[, log2FC := log2(meanShare_Niche_DDLS/meanShare_Niche_WDLS)]
bulktests[,neglog10q := log10(q_val)*-1]
bulktests[q_val < 0.1, labelid := xmCT]
ggplot(bulktests[q_val < 0.2], aes(x=log2FC, y=neglog10q)) + geom_point(aes(size = meanShare_Niche_DDLS)) + geom_text_repel(aes(label=labelid), force = 5) + theme_classic() + theme(text = element_text(size=24))