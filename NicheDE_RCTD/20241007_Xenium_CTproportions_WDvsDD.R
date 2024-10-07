library(data.table)
library(Seurat)
library(qvalue)
library(ggplot2)
library(ggrepel)
library(parallel)
library(poolr)
library(dplyr)

adata.obs <- fread("/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20240923cellcharterXenium_auto_k2_k20/20240924_cc_adata_LS_Xenium_k11k5clstrd_adata.obs_.csv")
countMtx <- fread('/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20240923cellcharterXenium_auto_k2_k20/20240924_cc_adata_LS_Xenium_k11k5clstrd_CountMat_.csv')
genes_dt <- fread('/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20240923cellcharterXenium_auto_k2_k20/20240924_cc_LS_Xenium_k11k5clstrd_adata.var_.csv')

#countMtx[, cell_ID := adata.obs$cell_id]
colnames(countMtx) = genes_dt$V1
transcripts_per_cell <- rowSums(countMtx)
countMtx_fdt <- countMtx[-1, ]
countMtx_fdt[, cell_ID := adata.obs$cell_id]

resall_dt <- rbindlist(lapply(1:9, function(xR){
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
  pdf(file = paste0('/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20241004_CellCharter_Xenium/Diff_WDvsDD/Region',xR,'_DDLSvsWDLS.pdf'), width = 18, height = 11)
  print(ggplot(fresCTs, aes(x=log2FC, y=neglog10q)) + geom_point(aes(col=colr, size = DDLS_share)) + geom_text_repel(aes(label=labelid)) + theme_classic() + theme(text = element_text(size=24)) + expand_limits(x = c(min(fresCTs$log2FC)-1,max(fresCTs$log2FC)+1), y=(c(max(fresCTs$neglog10q)+10))))    # force=3 
  dev.off()
  return(fresCTs)
}))

write.csv2(resall_dt, file = '/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20241004_CellCharter_Xenium/Diff_WDvsDD/20241004regionsAllresall_dt_DDLSvsWDLS.csv'
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
write.csv2(res_merged, file = '/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20241004_CellCharter_Xenium/Diff_WDvsDD/20241004regions_mergedFisher_DDLSvsWDLS.csv'
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
write.csv2(bulktests, file = '/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20241004_CellCharter_Xenium/Diff_WDvsDD/20241004regionsAll_bulktests_DDLSvsWDLS.csv'
           , row.names = F)
bulktests[, log2FC := log2(meanShare_Niche_DDLS/meanShare_Niche_WDLS)]
bulktests[,neglog10q := log10(q_val)*-1]
bulktests[q_val < 0.1, labelid := xmCT]

pdf(file = '/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20241004_CellCharter_Xenium/Diff_WDvsDD/Bulktests.pdf')
print(ggplot(bulktests[q_val < 0.2], aes(x=log2FC, y=neglog10q)) + geom_point(aes(size = meanShare_Niche_DDLS)) + geom_text_repel(aes(label=labelid), force = 5) + theme_classic() + theme(text = element_text(size=24)))
dev.off()

# Exclude Region 1 and 8 because they do not include both WDLS and DDLS parts.

resall_dt_read <- fread("/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20241004_CellCharter_Xenium/Diff_WDvsDD/20241004regionsAllresall_dt_DDLSvsWDLS.csv", sep=";", colClasses=c("xCT"="character", "DDLS_share"="numeric", "WDLS_share"="numeric", "p_val"="numeric", "OR"="numeric", "CI_l"="numeric", "CI_H"="numeric", "q_val"="numeric", "region"="integer", "log2FC"="numeric", "neglog10q"="numeric", "colr"="character"))

head(resall_dt_read)

resall_dt <- resall_dt[!(region %in% c(1,8))]

unique(resall_dt$region)

unique(resall_dt$xCT)

# Fisher's methods for consitent changes #####
xCTm= unique(resall_dt$xCT)[1]
xCTm = "adipocyte"
res_merged <- rbindlist(lapply(unique(resall_dt$xCT), function(xCTm) {
  print(xCTm)
  xmGresTbl <- resall_dt[xCT == xCTm]
  if ((sum(xmGresTbl$OR > 1) == 7) | (sum(xmGresTbl$OR < 1) == 7) ) {
    merged_res_xG <- xmGresTbl[1, 1]
    merged_res_xG[, labelid := xmGresTbl[1, 1]]
    merged_res_xG[, meanShare_DDLS := mean(as.numeric(xmGresTbl$DDLS_share))]
    merged_res_xG[, meanShare_WDLS := mean(as.numeric(xmGresTbl$WDLS_share))]
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
write.csv2(res_merged, file = '/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20241004_CellCharter_Xenium/Diff_WDvsDD/20241004regions_mergedFisher_DDLSvsWDLS_no1_8.csv'
           , row.names = F)

lapply(unique(resall_dt$xCT), function(xCTm) {
  xmGresTbl <- resall_dt[xCT == xCTm]
  sumORabove1 <- sum(xmGresTbl$OR > 1)
  sumORbelow1 <- sum(xmGresTbl$OR < 1)
  write(paste0("For Cell Type ",xCTm,", OR was above 1 in ",sumORabove1," Regions, and below 1 in ",sumORbelow1," Regions."), file='/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20241004_CellCharter_Xenium/Diff_WDvsDD/20241004_Regions_OR_Sum.txt', sep="\n", append=TRUE)
})

# Plot Fisher`s Method Results.

res_dt_ints <- resall_dt[xCT %in% c(res_merged$xG, 'ASPC')]

pdf(file = paste0('/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20241004_CellCharter_Xenium/Diff_WDvsDD/20241004region_merged_DDLSvsWDLS_no1_8.pdf'))
print(ggplot(res_merged, aes(x=mean_log2FC, y=q_val)) + geom_point(aes(size = meanShare_DDLS)) + geom_text_repel(aes(label=labelid)) + theme_classic() + theme(text = element_text(size=24)))    # force=3 
dev.off()

# pseudobulk####
xmCT <- "adipocyte"
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
write.csv2(bulktests, file = '/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20241004_CellCharter_Xenium/Diff_WDvsDD/20241004regionsAll_bulktests_DDLSvsWDLS_no1_8.csv'
           , row.names = F)
bulktests[, log2FC := log2(meanShare_Niche_DDLS/meanShare_Niche_WDLS)]
bulktests[,neglog10q := log10(q_val)*-1]
bulktests[q_val < 0.1, labelid := xmCT]

pdf(file = '/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20241004_CellCharter_Xenium/Diff_WDvsDD/Bulktests_no1_8.pdf', height=3, width=10)
print(ggplot(bulktests[q_val < 0.2], aes(x=log2FC, y=q_val)) + geom_point(aes(size = meanShare_Niche_DDLS)) + geom_text_repel(aes(label=labelid), force = 5) + theme_classic() + theme(text = element_text(size=24)))
dev.off()