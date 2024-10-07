library(data.table)
library(Seurat)
library(qvalue)
library(ggplot2)
library(ggrepel)
library(parallel)
library(poolr)

adata.obs <- fread('/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20240923cellcharterXenium_auto_k2_k20/20240924_cc_adata_LS_Xenium_k11k5clstrd_adata.obs_.csv')
countMtx <- fread('/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20240923cellcharterXenium_auto_k2_k20/20240924_cc_adata_LS_Xenium_k11k5clstrd_CountMat_.csv')
countMtx[1:10, 1:10]
table(adata.obs$sample)
colnames(adata.obs)
table(adata.obs$spatial_cluster_k11)

genes_dt <- fread('/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20240923cellcharterXenium_auto_k2_k20/20240924_cc_LS_Xenium_k11k5clstrd_adata.var_.csv')
#countMtx[, cell_ID := adata.obs$cell_id]
colnames(countMtx) = genes_dt$V1
transcripts_per_cell <- rowSums(countMtx)
countMtx_fdt <- countMtx[-1, ]
countMtx_fdt[, cell_ID := adata.obs$cell_id]
colSums(countMtx[-1, ])

# Perform Wilcoxon Rank Sum Test for each Region individually.
# In this Version Region 1 (No DDLS Part) and Region 8 (Either Only DDLS or WDLS to be determined) were excluded.

#region loop
xR = 2
resall_dt <- rbindlist(lapply(c(2,3,4,5,6,7,9), function(xR){
  print(xR)
  adata.obs_xR <- adata.obs[sample == xR]
  # get preadipocytes
  Preadipo.counts_xR_WDLS <- countMtx_fdt[cell_ID %in% adata.obs_xR[cell_types == 'ASPC' & subtype == 'WDLS', cell_id]]
  Preadipo.counts_xR_DDLS <- countMtx_fdt[cell_ID %in% adata.obs_xR[cell_types == 'ASPC' & grepl(pattern = 'DDLS', x = subtype), cell_id]]
  bulkExp_WDLS <- colMeans(Preadipo.counts_xR_WDLS[, -378])
  WDLSintgenes <- names(bulkExp_WDLS[bulkExp_WDLS > 0.01])
  bulkExp_DDLS <- colMeans(Preadipo.counts_xR_DDLS[, -378])
  DDLSintgenes <- names(bulkExp_DDLS[bulkExp_DDLS > 0.01])
  intgenes <- unique(c(DDLSintgenes, WDLSintgenes))
  #loop over expressed genes
  xG = intgenes[1]
  WXres <- rbindlist(mclapply(intgenes, function(xG){
    print(xG)
    idxG <- colnames(countMtx_fdt) == xG
    Niche_DDLS_counts <- unlist(Preadipo.counts_xR_DDLS[, ..idxG])
    Niche_WDLS_counts <- unlist(Preadipo.counts_xR_WDLS[, ..idxG])
    resxG <- wilcox.test(Niche_DDLS_counts, Niche_WDLS_counts, exact = T, conf.int = T)
    resdt <- as.data.table(xG)
    resdt[, medianCount_Niche_DDLS := median(Niche_DDLS_counts)]
    resdt[, medianCount_Niche_WDLS := median(Niche_WDLS_counts)]
    resdt[, meanCount_Niche_DDLS := mean(Niche_DDLS_counts)]
    resdt[, meanCount_Niche_WDLS := mean(Niche_WDLS_counts)]
    resdt[, p_val:= resxG$p.value]
    resdt[, estimate:= resxG$estimate]
    return(resdt)
  }, mc.cores = 4))
  WXres_old <- WXres
  WXres[p_val < 0, p_val :=0]
  WXres[p_val > 1, p_val :=1]
  WXres[, q_val := qvalue(p_val, pi0 = 1)$qvalue]
  write.csv2(WXres, paste0('/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20241004_CellCharter_Xenium/ASPC_WDvsDD_no1_8/20241004_',xR,'DDLSvsWDLS', '.csv'),quote=FALSE, row.names=FALSE, col.names = TRUE)
  WXres[, region := xR]
  WXres[,log2FC := log2(medianCount_Niche_DDLS/ medianCount_Niche_WDLS), by = medianCount_Niche_WDLS]
  #Compute Log2FC of Mean Expression
  WXres[,log2FC_mean := log2(meanCount_Niche_DDLS/ meanCount_Niche_WDLS)]
  WXres[,neglog10q := log10(q_val)*-1]
  # Step 1: Calculate the range of non-infinite values
  neglog10q_range <- range(WXres[is.finite(neglog10q), neglog10q])  # Get range of non-Inf y values
  neglog10q_padding <- 0.1 * diff(neglog10q_range)  # Add some padding to the range
  # Step 2: Replace Inf with values slightly beyond the non-Inf range for plotting
  WXres[, neglog10q_plot := ifelse(neglog10q == Inf, neglog10q_range[2] + neglog10q_padding,  # Plot Inf slightly above max
                    ifelse(neglog10q == -Inf, neglog10q_range[1] - neglog10q_padding, neglog10q))]  # Plot -Inf slightly below min
  WXres[q_val < 0.05 & (log2FC < 0) ,colr := 'blue']
  WXres[q_val < 0.05 & (log2FC > 0) ,colr := 'red']
  WXres[q_val < 0.05 & (log2FC_mean < 0) ,colr_mean := 'blue']
  WXres[q_val < 0.05 & (log2FC_mean > 0) ,colr_mean := 'red']
  WXres[q_val < 0.1 ,labelid := xG]
  pdf(file = paste0('/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20241004_CellCharter_Xenium/ASPC_WDvsDD_no1_8/20241004_',xR,'DDLSvsWDLS', '.pdf'), width = 16, height = 9)
  print(ggplot(WXres, aes(x=log2FC, y=neglog10q_plot)) + geom_point(aes(col=colr, size = medianCount_Niche_DDLS)) + geom_text_repel(aes(label=labelid), force = 3) + coord_cartesian(ylim = neglog10q_range) + theme_classic() + theme(text = element_text(size=24)) )
  dev.off()
  pdf(file = paste0('/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20241004_CellCharter_Xenium/ASPC_WDvsDD_no1_8/20241004_',xR,'DDLSvsWDLS_meanFC', '.pdf'), width = 16, height = 9)
  print(ggplot(WXres, aes(x=log2FC_mean, y=neglog10q_plot)) + geom_point(aes(col=colr_mean, size = meanCount_Niche_DDLS)) + geom_text_repel(aes(label=labelid), force = 3) + coord_cartesian(ylim = neglog10q_range) + theme_classic() + theme(text = element_text(size=24)) )
  dev.off()
  return(WXres)
}))

resall_dt[, inc := .N, by = xG]
resall_dt_int <- resall_dt[inc == 7]

write.csv2(resall_dt, "/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20241004_CellCharter_Xenium/ASPC_WDvsDD_no1_8/20241004_Resall.csv")

# Fisher's methods for consitent changesd

# For mean Count FC

xmG = 'MDM2'
res_merged <- rbindlist(lapply(unique(resall_dt_int$xG), function(xmG) {
  print(xmG)
  xmGresTbl <- resall_dt_int[xG == xmG]
  if ((sum(xmGresTbl$estimate > 0) == 7) | (sum(xmGresTbl$estimate < 0) == 7) ) {
    merged_res_xG <- xmGresTbl[1, 1:3]
    merged_res_xG[, meanCount_Niche_DDLS := mean(xmGresTbl$meanCount_Niche_DDLS)]
    merged_res_xG[, meanCount_Niche_WDLS := mean(xmGresTbl$meanCount_Niche_WDLS)]
    merged_res_xG[, merged_pval := fisher(p = xmGresTbl$p_val)$p]
    merged_res_xG[, median_estimate := median(xmGresTbl$estimate)]
    merged_res_xG[, mean_log2FC := log2(meanCount_Niche_DDLS/meanCount_Niche_WDLS)]
    return(merged_res_xG)
  } else  {
    return(data.table())
  }
}))
res_merged[merged_pval < 0, merged_pval :=0]
res_merged[merged_pval > 1, merged_pval :=1]
res_merged[, q_val := qvalue(merged_pval, pi0 = 1)$qvalue]
write.csv2(res_merged, file = '/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20241004_CellCharter_Xenium/ASPC_WDvsDD_no1_8/20241004_mergedFisher_DDLSvsWDLS_mean_log2FC.csv'
           , row.names = F)


# For Median Count FC

xmG = 'MDM2'
res_merged_median <- rbindlist(lapply(unique(resall_dt_int$xG), function(xmG) {
  print(xmG)
  xmGresTbl <- resall_dt_int[xG == xmG]
  if ((sum(xmGresTbl$estimate > 0) == 7) | (sum(xmGresTbl$estimate < 0) == 7) ) {
    merged_res_xG <- xmGresTbl[1, 1:3]
    merged_res_xG[, medianCount_Niche_DDLS := median(xmGresTbl$medianCount_Niche_DDLS)]
    merged_res_xG[, medianCount_Niche_WDLS := median(xmGresTbl$medianCount_Niche_WDLS)]
    merged_res_xG[, merged_pval := fisher(p = xmGresTbl$p_val)$p]
    merged_res_xG[, median_estimate := median(xmGresTbl$estimate)]
    merged_res_xG[, median_log2FC := median(xmGresTbl$log2FC)]
    return(merged_res_xG)
  } else  {
    return(data.table())
  }
}))
res_merged_median[merged_pval < 0, merged_pval :=0]
res_merged_median[merged_pval > 1, merged_pval :=1]
res_merged_median[, q_val := qvalue(merged_pval, pi0 = 1)$qvalue]
write.csv2(res_merged_median, file = '/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20241004_CellCharter_Xenium/ASPC_WDvsDD_no1_8/20241004_mergedFisher_DDLSvsWDLS_median_log2FC.csv'
           , row.names = F)

# Plot Fisher's Method Results.

res_dt_mean_ints <- resall_dt_int[xG %in% c(res_merged$xG, 'MDM2', 'CDK4')]

pdf(file = paste0('/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20241004_CellCharter_Xenium/ASPC_WDvsDD_no1_8/20241004region_mergedDDLSvsWDLS_meanFC.pdf'), height=5, width=30)
print(ggplot(res_merged, aes(x=mean_log2FC, y=q_val)) + geom_point(aes(size = meanCount_Niche_DDLS)) + geom_text_repel(aes(label=xG), force = 3) + theme_classic() + theme(text = element_text(size=24)))
dev.off()                                                                                                                                                                                

res_dt_median_ints <- resall_dt_int[xG %in% c(res_merged_median$xG, 'MDM2', 'CDK4')]

pdf(file = paste0('/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20241004_CellCharter_Xenium/ASPC_WDvsDD_no1_8/20241004region_mergedDDLSvsWDLS_medianFC.pdf'), height=5, width=30)
print(ggplot(res_merged_median, aes(x=median_log2FC, y=q_val)) + geom_point(aes(size = medianCount_Niche_DDLS)) + geom_text_repel(aes(label=xG), force = 3) + theme_classic() + theme(text = element_text(size=24)))
dev.off() 

# pseudobulk####
bulktests <- rbindlist(lapply(unique(resall_dt_int$xG), function(xmG) {
  print(xmG)
  xmGresTbl <- resall_dt_int[xG == xmG]
  resxG <- wilcox.test(xmGresTbl$meanCount_Niche_WDLS, xmGresTbl$meanCount_Niche_DDLS, exact = T, conf.int = T, paired = T)
  resdt <- as.data.table(xmG)
  resdt[, meanCount_Niche_DDLS := mean(xmGresTbl$meanCount_Niche_DDLS)]
  resdt[, meanCount_Niche_WDLS := mean(xmGresTbl$meanCount_Niche_WDLS)]
  resdt[, sdCount_Niche_DDLS := sd(xmGresTbl$meanCount_Niche_DDLS)]
  resdt[, sdCount_Niche_WDLS := sd(xmGresTbl$meanCount_Niche_WDLS)]
  resdt[, p_val:= resxG$p.value]
  resdt[, estimate:= resxG$estimate]
  resdt[, conf_est_l:= resxG$conf.int[1]]
  resdt[, conf_est_h:= resxG$conf.int[2]]
  return(resdt)
}))
bulktests[p_val < 0, p_val :=0]
bulktests[p_val > 1, p_val :=1]
bulktests[, q_val := qvalue(p_val)$qvalue]
write.csv2(bulktests, file = '/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20241004_CellCharter_Xenium/ASPC_WDvsDD_no1_8/20241004regionsAll_bulktests_DDLSvsWDLS.csv'
           , row.names = F)
bulktests[, log2FC := log2(meanCount_Niche_DDLS/meanCount_Niche_WDLS)]
bulktests[,neglog10q := log10(q_val)*-1]
bulktests[q_val < 0.1, labelid := xmG]

pdf(file='/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20241004_CellCharter_Xenium/ASPC_WDvsDD_no1_8/20241004regionsAll_bulktests_DDLSvsWDLS.pdf', width=30)
print(ggplot(bulktests[q_val < 0.2], aes(x=log2FC, y=neglog10q)) + geom_point(aes(size = meanCount_Niche_DDLS)) + geom_text_repel(aes(label=xmG), force = 5) + theme_classic() + theme(text = element_text(size=24)))
dev.off()      

# indivCounts ####
idxA <- colnames(adata.obs) %in% c("cell_id" ,"subtype", "cell_types" ,"sample")
xG = 'PDGFRA'
Counts_int_dt <- rbindlist(mclapply(unique(c(bulktests[(meanCount_Niche_WDLS >1 | meanCount_Niche_DDLS > 1) & q_val < 0.1 & (log2FC >1 | log2FC < -1) , xmG], 
                                           res_merged[(mean_log2FC >1 | mean_log2FC < -1) & (meanCount_Niche_WDLS >1 | meanCount_Niche_DDLS > 1), xG]
                                           , 'MDM2')), function(xG){
  print(xG)
  idxG <- colnames(countMtx_fdt) == xG
  idxG[378] =T
  xG_counts <- countMtx_fdt[, ..idxG]
  xR_counts_xG <- cbind(adata.obs[, ..idxA], xG_counts)
  xR_counts_xG <- xR_counts_xG[cell_types == 'ASPC']
  #xR_counts_xG[grepl(pattern = 'DDLS', x = spatial_cluster), subtype := 'DDLS']
  #xR_counts_xG[grepl(pattern = 'WDLS', x = spatial_cluster), subtype := 'WDLS']
  xR_counts_xGf <- xR_counts_xG[sample %in% c(2, 3, 4, 5, 6, 7, 9) & (subtype %in% c("WDLS","DDLS"))]
  xR_counts_xGf[,region:=paste0("Region",sample)]
  colnames(xR_counts_xGf)[5] = 'Counts'
  pdf(file = paste0('/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20241004_CellCharter_Xenium/ASPC_WDvsDD_no1_8/DiffGenes/20241007_DDLSvsWDLS_', xG, '.pdf'), width = 9, height = 16)
  print(ggplot(xR_counts_xGf, aes(x=subtype, y=Counts)) + geom_jitter(aes(col=region), size = 0.05) + geom_boxplot(outliers = F, aes(fill=region)) + theme_classic() + theme(text = element_text(size=24)) + 
          labs(subtitle= xG) + ylim(0,25))
  dev.off()
  xR_counts_xGf[, xiG := xG]
  return(xR_counts_xGf)
}, mc.cores = 4))
                    
# run pathway enrichment analysis

