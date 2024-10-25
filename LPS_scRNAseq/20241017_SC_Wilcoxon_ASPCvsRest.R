library(data.table)
library(Seurat)
library(qvalue)
library(ggplot2)
library(ggrepel)
library(parallel)
library(poolr)

dir.create('/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20241017_SingleCell')



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

#region loop
xR = 2
resall_dt <- rbindlist(mclapply(1:9, function(xR){
  print(xR)
  adata.obs_xR <- adata.obs[sample == xR]
  # get preadipocytes
  Preadipo.counts_xR <- countMtx_fdt[cell_ID %in% adata.obs_xR[cell_types == 'ASPC', cell_id]]
  Rest.counts_xR <- countMtx_fdt[cell_ID %in% adata.obs_xR[cell_types != 'ASPC', cell_id]]
  bulkExp_ASPC <- colMeans(Preadipo.counts_xR[, -378])
  ASPCintgenes <- names(bulkExp_ASPC[bulkExp_ASPC > 0.01])
  bulkExp_Rest <- colMeans(Rest.counts_xR[, -378])
  Restintgenes <- names(bulkExp_Rest[bulkExp_Rest > 0.01])
  intgenes <- unique(c(ASPCintgenes, Restintgenes))
  #loop over expressed genes
  xG = intgenes[1]
  WXres <- rbindlist(mclapply(intgenes, function(xG){
    print(xG)
    idxG <- colnames(countMtx_fdt) == xG
    Niche_ASPC_counts <- unlist(Preadipo.counts_xR[, ..idxG])
    Niche_Rest_counts <- unlist(Rest.counts_xR[, ..idxG])
    resxG <- wilcox.test(Niche_ASPC_counts, Niche_Rest_counts, exact = T, conf.int = T)
    resdt <- as.data.table(xG)
    resdt[, medianCount_Niche_ASPC := median(Niche_ASPC_counts)]
    resdt[, medianCount_Niche_Rest := median(Niche_Rest_counts)]
    resdt[, meanCount_Niche_ASPC := mean(Niche_ASPC_counts)]
    resdt[, meanCount_Niche_Rest := mean(Niche_Rest_counts)]
    resdt[, p_val:= resxG$p.value]
    resdt[, estimate:= resxG$estimate]
    return(resdt)
  }, mc.cores = 4))
  WXres_old <- WXres
  WXres[p_val < 0, p_val :=0]
  WXres[p_val > 1, p_val :=1]
  WXres[, q_val := qvalue(p_val, pi0 = 1)$qvalue]
  write.csv2(WXres, paste0('/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20241004_CellCharter_Xenium/ASPCvsRest/20241017_',xR,'ASPCvsRest', '.csv'),quote=FALSE, row.names=FALSE, col.names = TRUE)
  WXres[, region := xR]
  WXres[,log2FC := log2(medianCount_Niche_ASPC/ medianCount_Niche_Rest), by = medianCount_Niche_ASPC]
  #Compute Log2FC of Mean Expression
  WXres[,log2FC_mean := log2(meanCount_Niche_ASPC/ meanCount_Niche_Rest)]
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
  pdf(file = paste0('/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20241004_CellCharter_Xenium/ASPCvsRest/20241017_',xR,'ASPCvsRest', '.pdf'), width = 16, height = 9)
  print(ggplot(WXres, aes(x=log2FC, y=neglog10q_plot)) + geom_point(aes(col=colr, size = medianCount_Niche_ASPC)) + geom_text_repel(aes(label=labelid), force = 3) + coord_cartesian(ylim = neglog10q_range) + theme_classic() + theme(text = element_text(size=24)) )
  dev.off()
  pdf(file = paste0('/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20241004_CellCharter_Xenium/ASPCvsRest/20241017_',xR,'ASPCvsRest_meanFC', '.pdf'), width = 16, height = 9)
  print(ggplot(WXres, aes(x=log2FC_mean, y=neglog10q_plot)) + geom_point(aes(col=colr_mean, size = meanCount_Niche_ASPC)) + geom_text_repel(aes(label=labelid), force = 3) + coord_cartesian(ylim = neglog10q_range) + theme_classic() + theme(text = element_text(size=24)) )
  dev.off()
  return(WXres)
}, mc.cores=10))

resall_dt[, inc := .N, by = xG]
resall_dt_int <- resall_dt[inc == 9]

write.csv2(resall_dt, "/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20241004_CellCharter_Xenium/ASPCvsRest/20241017_Resall.csv")

# Fisher's methods for consitent changes

# For mean Count FC

xmG = 'MDM2'
res_merged <- rbindlist(lapply(unique(resall_dt_int$xG), function(xmG) {
  print(xmG)
  xmGresTbl <- resall_dt_int[xG == xmG]
  if ((sum(xmGresTbl$estimate > 0) == 9) | (sum(xmGresTbl$estimate < 0) == 9) ) {
    merged_res_xG <- xmGresTbl[1, 1:3]
    merged_res_xG[, meanCount_Niche_ASPC := mean(xmGresTbl$meanCount_Niche_ASPC)]
    merged_res_xG[, meanCount_Niche_Rest := mean(xmGresTbl$meanCount_Niche_Rest)]
    merged_res_xG[, merged_pval := fisher(p = xmGresTbl$p_val)$p]
    merged_res_xG[, median_estimate := median(xmGresTbl$estimate)]
    merged_res_xG[, mean_log2FC := log2(meanCount_Niche_ASPC/meanCount_Niche_Rest)]
    return(merged_res_xG)
  } else  {
    return(data.table())
  }
}))
res_merged[merged_pval < 0, merged_pval :=0]
res_merged[merged_pval > 1, merged_pval :=1]
res_merged[, q_val := qvalue(merged_pval, pi0 = 1)$qvalue]
write.csv2(res_merged, file = '/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20241004_CellCharter_Xenium/ASPCvsRest/20241017_mergedFisher_ASPCvsRest_mean_log2FC.csv'
           , row.names = F)

 + theme(text = element_text(size=24))
# Plot Fisher's Method Results.

res_dt_mean_ints <- resall_dt_int[xG %in% c(res_merged$xG, 'MDM2', 'CDK4')]

pdf(file ='/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20241004_CellCharter_Xenium/ASPCvsRest/20241017region_mergedASPCvsRest_meanFC.pdf')
print(ggplot(res_merged, aes(x=mean_log2FC, y=q_val)) + geom_point(aes(size = meanCount_Niche_ASPC)) + geom_text_repel(aes(label=xG)) + theme_classic())
dev.off()                                                                                                                                                                                


# pseudobulk####
bulktests <- rbindlist(lapply(unique(resall_dt_int$xG), function(xmG) {
  print(xmG)
  xmGresTbl <- resall_dt_int[xG == xmG]
  resxG <- wilcox.test(xmGresTbl$meanCount_Niche_ASPC, xmGresTbl$meanCount_Niche_Rest, exact = T, conf.int = T, paired = T)
  resdt <- as.data.table(xmG)
  resdt[, meanCount_Niche_ASPC := mean(xmGresTbl$meanCount_Niche_ASPC)]
  resdt[, meanCount_Niche_Rest := mean(xmGresTbl$meanCount_Niche_Rest)]
  resdt[, sdCount_Niche_ASPC := sd(xmGresTbl$meanCount_Niche_ASPC)]
  resdt[, sdCount_Niche_Rest := sd(xmGresTbl$meanCount_Niche_Rest)]
  resdt[, p_val:= resxG$p.value]
  resdt[, estimate:= resxG$estimate]
  resdt[, conf_est_l:= resxG$conf.int[1]]
  resdt[, conf_est_h:= resxG$conf.int[2]]
  return(resdt)
}))
bulktests[p_val < 0, p_val :=0]
bulktests[p_val > 1, p_val :=1]
bulktests[, q_val := qvalue(p_val)$qvalue]
write.csv2(bulktests, file = '/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20241004_CellCharter_Xenium/ASPCvsRest/20241017_Bulk_ASPCvsRest_mean_log2FC.csv'
           , row.names = F)
bulktests[, log2FC := log2(meanCount_Niche_ASPC/meanCount_Niche_Rest)]
bulktests[,neglog10q := log10(q_val)*-1]
bulktests[q_val < 0.1, labelid := xmG]

pdf(file = paste0('/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20241004_CellCharter_Xenium/ASPCvsRest/20241017_Bulk_ASPCvsRest_meanFC.pdf'), width=15, height=15)
ggplot(bulktests[q_val < 0.1], aes(x=log2FC, y=neglog10q)) + geom_point(aes(size = meanCount_Niche_ASPC)) + geom_text_repel(aes(label=xG), force = 5) + theme_classic() + theme(text = element_text(size=24))
dev.off()


# indivCounts ####
idxA <- colnames(adata.obs) %in% c("cell_id" , "spatial_cluster", "cell_types" ,"sample")
xG = 'PDGFRA'
Counts_int_dt <- rbindlist(mclapply(unique(c(bulktests[(meanCount_Niche_ASPC >1 | meanCount_Niche_Rest > 1) & q_val < 0.1 & (log2FC >1 | log2FC < -1) , xmG], 
                                           res_merged[(mean_log2FC >1 | mean_log2FC < -1) & (meanCount_Niche_ASPC >1 | meanCount_Niche_Rest > 1), xG]
                                           , 'MDM2')), function(xG){
  print(xG)
  idxG <- colnames(countMtx_fdt) == xG
  idxG[378] =T
  xG_counts <- countMtx_fdt[, ..idxG]
  xR_counts_xG <- cbind(adata.obs[, ..idxA], xG_counts)
  #xR_counts_xG <- xR_counts_xG[cell_types == 'preadipocyte']
  xR_counts_xG[grepl(pattern = 'ASPC', x = cell_types), subtype := 'ASPC']
  xR_counts_xG[!grepl(pattern = 'ASPC', x = cell_types), subtype := 'Rest']
  xR_counts_xGf <- xR_counts_xG[sample %in% c('1',"2","3","4","5","6","7","8","9") & !is.na(subtype)]
  colnames(xR_counts_xGf)[5] = 'Gene'
  pdf(file = paste0('/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20241004_CellCharter_Xenium/ASPCvsRest/20241017_ASPCvsRest_', xG, '.pdf'), width = 9, height = 16)
  print(ggplot(xR_counts_xGf, aes(x=subtype, y=Gene)) + geom_jitter(aes(col=as.character(sample)), size = 0.05) + geom_boxplot(outliers = F, aes(fill=as.character(sample))) + theme_classic() + theme(text = element_text(size=24)) + 
          labs(subtitle= xG) + ylim(0,25))
  dev.off()
  xR_counts_xGf[, xiG := xG]
  return(xR_counts_xGf)
}, mc.cores = 10))