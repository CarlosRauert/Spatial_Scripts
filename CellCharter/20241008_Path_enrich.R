library(data.table)
library(Seurat)
library(qvalue)
library(ggplot2)
library(ggrepel)
library(parallel)
library(poolr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(data.table)
library(fgsea)
library(msigdbr)
library(ReactomePA)

# Get Bulk DE Genes
BulkTests <- fread('/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20241004_CellCharter_Xenium/ASPC_WDvsDD_no1_8/20241004regionsAll_bulktests_DDLSvsWDLS.csv')
Bulk_Sign <- subset(BulkTests, q_val<=0.2)
Bulk_Sign[, log2FC := log2(meanCount_Niche_DDLS/meanCount_Niche_WDLS)]
upregd_Bulk <- subset(Bulk_Sign, log2FC>0)
downregd_Bulk <- subset(Bulk_Sign, log2FC<0)
# Get Fisher DE Genes
Merged_Fisher_MeanFC <- fread('/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20241004_CellCharter_Xenium/ASPC_WDvsDD_no1_8/20241004_mergedFisher_DDLSvsWDLS_mean_log2FC.csv', sep=";")
Merged_Fisher_MeanFC[, log2FC := mean_log2FC]
Merged_Fisher_MeanFC[, xmG := xG]
upregd_Fisher <- subset(Merged_Fisher_MeanFC, mean_log2FC>0)
downregd_Fisher <- subset(Merged_Fisher_MeanFC, mean_log2FC<0)

SetsLst <- list(upregd_Bulk=upregd_Bulk,downregd_Bulk=downregd_Bulk, upregd_Fisher=upregd_Fisher, downregd_Fisher=downregd_Fisher)
# loop over Gene Sets
Sign = 1
Enri <- function(Sign){
  # Create Plot Title
  varname=names(SetsLst)[Sign]
  title=paste0(ifelse(grepl("upregd",varname),"Upregulated Pathways in DDLS ASPCs,\n","Downregulated Pathways in DDLS ASPCs,\n"),ifelse(grepl("Fisher",varname), "Test: Wilcoxon Rank Sum Test / Fisher's Method", "Test: Wilcoxon Signed Rank Test"))
  # Make Gene sets, sort by FC, covert names to IDs
  gene_names <- SetsLst[[Sign]]$xmG
  fold_changes <- SetsLst[[Sign]]$log2FC
  gene_list <- setNames(fold_changes, gene_names)
  gene_list <- sort(abs(gene_list), decreasing = TRUE)
  genes_entrez <- bitr(names(gene_list), fromType = "SYMBOL", 
                       toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  # Merge the Entrez IDs with fold changes
  gene_list_entrez <- gene_list[genes_entrez$SYMBOL]
  names(gene_list_entrez) <- genes_entrez$ENTREZID
  
  # Run Kegg Pathway analysis
  kegg_enrich <- enrichKEGG(gene = names(gene_list_entrez),
                          organism = "hsa",     # hsa = human
                          pvalueCutoff = 0.05)
  if (nrow(kegg_enrich)>0){
    # Save KEGG
    head(kegg_enrich)
    # Get gene names for each KEGG pathway
    kegg_enrich@result$geneID <- lapply(kegg_enrich@result$geneID, function(gene_ids) {
      # Get gene symbols from the Entrez IDs
      symbols <- genes_entrez$SYMBOL[match(unlist(strsplit(gene_ids, "/")), genes_entrez$ENTREZID)]
      paste(symbols, collapse = ", ")  # Join symbols with a comma
    })
    kegg_df <- as.data.frame(kegg_enrich)
    #write.csv2(kegg_df, paste0("/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20241004_CellCharter_Xenium/ASPC_WDvsDD_no1_8/20241010_DDLSvsWDLS_",varname,"_keggEnrich.csv"))
    kegg_df$GeneRatio <- sapply(kegg_df$GeneRatio, function(x) {
      ratio <- as.numeric(unlist(strsplit(x, "/")))
      ratio[1] / ratio[2]  # Get the ratio value
    })
    kegg_df <- kegg_df[order(-kegg_df$GeneRatio, kegg_df$p.adjust), ]
    kegg_top_results <- head(kegg_df, 20)
    kegg_top_results <- kegg_top_results[rev(rownames(kegg_top_results)), ]
    kegg_top_results$Description <- factor(kegg_top_results$Description, levels=kegg_top_results$Description)
    # Plot KEGG
    pdf(file=paste0("/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20241004_CellCharter_Xenium/ASPC_WDvsDD_no1_8/20241010_DDLSvsWDLS_",varname,"_keggEnrich.pdf"))
    print(
    ggplot(kegg_top_results, aes(x = Description, y = GeneRatio, color = p.adjust)) +
      geom_point() +                                # Change to geom_point for a dot plot
      coord_flip() +                                # Flip coordinates to make it horizontal
      labs(x = "KEGG Pathway",                           # Update axis labels as needed
          y = "Gene Ratio") + 
      geom_text_repel(aes(label = geneID),    # Repel text labels for clarity
                      size = 3, 
                      max.overlaps = Inf, 
                      nudge_x = 0.1, 
                      nudge_y = 0.1, 
                      direction = "both", 
                      box.padding = 0.5) +
      ggtitle(title) +                             # Add the title as before
      scale_size(range = c(3, 6)) +                # Adjust size range for the dots
      scale_color_gradient(low = "blue", high = "red") # Add gradient for p.adjust color mapping
      )
    dev.off()
  }
  
  # Run GO enrichment analysis
  go_enrich <- enrichGO(gene = names(gene_list_entrez), 
                      OrgDb = org.Hs.eg.db, 
                      keyType = "ENTREZID", 
                      ont = "BP",           # Ontology: BP = Biological Process
                      pvalueCutoff = 0.05)
  go_enrich_df <- as.data.frame(go_enrich)
  go_enrich_df$GeneRatio <- sapply(go_enrich_df$GeneRatio, function(x) {
    ratio <- as.numeric(unlist(strsplit(x, "/")))
    ratio[1] / ratio[2]  # Get the ratio value
  })
  write.csv(go_enrich_df, paste0("/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20241004_CellCharter_Xenium/ASPC_WDvsDD_no1_8/20241010_DDLSvsWDLS_",varname,"_enrichGO.csv"))
  go_enrich_df <- go_enrich_df[order(-go_enrich_df$GeneRatio, go_enrich_df$p.adjust), ]
  # Get gene names for each GO pathway
  gene_conversion <- bitr(names(gene_list), fromType = "SYMBOL", 
                        toType = "ENTREZID", 
                        OrgDb = org.Hs.eg.db)
  # Add gene symbols to the enrichment result
  # Split the Entrez IDs in the 'geneID' column
  go_enrich_df$geneID_split <- strsplit(go_enrich_df$geneID, "/")
  # Convert Entrez IDs to Gene Symbols
  go_enrich_df$gene_symbols <- sapply(go_enrich_df$geneID_split, function(ids) {
    symbols <- gene_conversion$SYMBOL[match(ids, gene_conversion$ENTREZID)]
    paste(symbols, collapse = ", ")  # Combine gene symbols as a string
  })
  # Plot results
  top_results <- head(go_enrich_df, 20)
  top_results <- top_results[rev(rownames(top_results)), ]
  top_results$Description <- factor(top_results$Description, levels=top_results$Description)

  pdf(file=paste0("/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20241004_CellCharter_Xenium/ASPC_WDvsDD_no1_8/20241010_DDLSvsWDLS_",varname,"_EnrichGO.pdf"), width=9)
  # Generate the dotplot with ggplot2
  print(
  ggplot(top_results, aes(x = Description, y = GeneRatio, color = p.adjust, size = GeneRatio)) +
    geom_point() +                                # Change to geom_point for a dot plot
    coord_flip() +                                # Flip coordinates to make it horizontal
    labs(x = "GO Term",                           # Update axis labels as needed
         y = "Gene Ratio") + 
    geom_text_repel(aes(label = gene_symbols),    # Repel text labels for clarity
                    size = 3, 
                    max.overlaps = Inf, 
                    nudge_x = 0.1, 
                    nudge_y = 0.1, 
                    direction = "both", 
                    box.padding = 0.5) +
    ggtitle(title) +                             # Add the title as before
    scale_size(range = c(3, 6)) +                # Adjust size range for the dots
    scale_color_gradient(low = "blue", high = "red") # Add gradient for p.adjust color mapping
  )
  dev.off()

  # Run GSEA Analysis
  pathways <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "KEGG")
  pathways_list <- split(pathways$entrez_gene, pathways$gs_name)
  fgsea_res <- fgsea(pathways = pathways_list, 
                   stats = gene_list_entrez,   # Named vector of log2 fold change
                   maxSize = 500,
                   scoreType="pos")
  # Filter By P Value
  fgsea_res_filtered <- fgsea_res[fgsea_res$padj < 0.05 ]
  # Sort By P Value
  sorted_results <- fgsea_res_filtered[order(-fgsea_res_filtered$NES, fgsea_res_filtered$pval), ]
  if (nrow(fgsea_res_filtered)>0){
    # Extract top 5 leading-edge genes for each pathway and concatenate them
    fgsea_res_filtered$leadingEdgeNames <- rep(NULL, times=nrow(fgsea_res_filtered))
    for (i in 1:length(fgsea_res_filtered$leadingEdge)){
     fgsea_res_filtered$leadingEdgeNames[[i]] <- list()
      for (j in 1:length(fgsea_res_filtered$leadingEdge[[i]])){
        fgsea_res_filtered$leadingEdgeNames[[i]][j]<- names(gene_list_orig)[which(names(gene_list_entrez)==fgsea_res_filtered$leadingEdge[[i]][j])]
      }
    }
    fgsea_res_filtered$leading_edge_geneNames <- sapply(fgsea_res_filtered$leadingEdgeNames, function(genes) {
      paste(head(genes, 5), collapse = ", ")  # Get top 5 genes and join them with commas
    })

    which(names(gene_list_entrez)==fgsea_res_filtered$leadingEdge[[1]][1])

    # Modify the pathway names to include leading edge genes (first 5)
    fgsea_res_filtered$pathway_with_genes <- paste0(fgsea_res_filtered$pathway, "\n(", fgsea_res_filtered$leading_edge_geneNames, ")")
    # Save Result
    write.csv2(fgsea_res_filtered,paste0("/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20241004_CellCharter_Xenium/ASPC_WDvsDD_no1_8/20241010_",varname,"_DDLSvsWDLS_GSEA.csv")) 
    # Create bar plot of NES (Normalized Enrichment Score)
    pdf(paste0("/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20241004_CellCharter_Xenium/ASPC_WDvsDD_no1_8/20241010_",varname,"_DDLSvsWDLS_GSEA.pdf"), width=12)
    print(ggplot(fgsea_res_filtered, aes(x = reorder(pathway_with_genes, NES), y = NES)) +
    geom_bar(stat = "identity", aes(fill = padj)) +  # Fill based on -log10 of padj
    coord_flip() +  # Flip for a horizontal bar plot
    labs(x = "Pathway (with top contributing genes)", 
         y = "Normalized Enrichment Score (NES)", 
         title = "top enriched Pathways in DDLS ASPC") +
    scale_fill_gradient(low = "blue", high = "red", name = "Adjusted P-value") +  # Color gradient
    theme_minimal() +
    theme(axis.text.y = element_text(size = 10)) +  # Adjust text size for clarity
    ggtitle(title)
    )
    dev.off()
  }
  
  # Run Reactome Pathway analysis
  reactome_results <- enrichPathway(gene = names(gene_list_entrez), 
                                  organism = "human",
                                  pvalueCutoff = 0.05,
                                  readable = TRUE)  # readable = TRUE returns gene symbols instead of Entrez IDs
  if (!is.null(reactome_results)){
    sorted_results <- as.data.frame(reactome_results)
    sorted_results <- reactome_results[order(-reactome_results$Count,reactome_results$p.adjust),]
    sorted_results <- sorted_results[rev(rownames(sorted_results)), ]
    sorted_results$Description <- factor(sorted_results$Description, levels=sorted_results$Description)
    write.csv(sorted_results, paste0("/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20241004_CellCharter_Xenium/ASPC_WDvsDD_no1_8/20241010_",varname,"_DDLSvsWDLS_Reactome.csv"))
    pdf(file=paste0("/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20241004_CellCharter_Xenium/ASPC_WDvsDD_no1_8/20241010_",varname,"_DDLSvsWDLS_Reactome.pdf"), width=11)
    print(ggplot(sorted_results, aes(x = Description, y = Count, color = p.adjust)) +
    geom_point() +                                # Change to geom_point for a dot plot
    coord_flip() +                                # Flip coordinates to make it horizontal
    labs(x = "Reactome Pathway",                           # Update axis labels as needed
         y = "Count") + 
    geom_text_repel(aes(label = geneID),    # Repel text labels for clarity
                    size = 3, 
                    max.overlaps = Inf, 
                    nudge_x = 0.1, 
                    nudge_y = 0.1, 
                    direction = "both", 
                    box.padding = 0.5) +
    ggtitle(title) +                             # Add the title as before
    scale_size(range = c(3, 6)) +                # Adjust size range for the dots
    scale_color_gradient(low = "blue", high = "red")) # Add gradient for p.adjust color mapping
    dev.off()
  }
}

lapply(1:4, Enri)