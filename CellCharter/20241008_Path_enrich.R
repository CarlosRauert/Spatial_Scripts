library(data.table)
library(Seurat)
library(qvalue)
library(ggplot2)
library(ggrepel)
library(parallel)
library(poolr)

# run pathway enrichment analysis

#install.packages("BiocManager")
#BiocManager::install("clusterProfiler")
#BiocManager::install("org.Hs.eg.db")  # for human gene annotation
#BiocManager::install("DOSE")          # for disease ontology enrichment (optional)
#BiocManager::install("enrichplot")

# Example: A named vector where names are genes and values are log2 fold changes

library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(data.table)

BulkTests <- fread('/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20241004_CellCharter_Xenium/ASPC_WDvsDD_no1_8/20241004regionsAll_bulktests_DDLSvsWDLS.csv')

# Bulk

Bulk_Sign <- subset(BulkTests, q_val<=0.2)
Bulk_Sign[, log2FC := log2(meanCount_Niche_DDLS/meanCount_Niche_WDLS)]

gene_names <- Bulk_Sign$xmG
fold_changes <- Bulk_Sign$log2FC
gene_list <- setNames(fold_changes, gene_names)
gene_list <- sort(gene_list, decreasing = TRUE)
gene_list_orig <- setNames(fold_changes, gene_names)
gene_list_orig <- sort(gene_list, decreasing = TRUE)

genes_entrez <- bitr(names(gene_list), fromType = "SYMBOL", 
                     toType = "ENTREZID", OrgDb = org.Hs.eg.db)
# Merge the Entrez IDs with fold changes
gene_list_entrez <- gene_list[genes_entrez$SYMBOL]
names(gene_list_entrez) <- genes_entrez$ENTREZID
kegg_enrich <- enrichKEGG(gene = names(gene_list_entrez),
                          organism = "hsa",     # hsa = human
                          pvalueCutoff = 0.05)
# View results
head(kegg_enrich)
write.csv2(kegg_enrich, "/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20241004_CellCharter_Xenium/ASPC_WDvsDD_no1_8/20241004_mergedFisher_DDLSvsWDLS_Bulk_keggEnrich.csv")

pdf(file="/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20241004_CellCharter_Xenium/ASPC_WDvsDD_no1_8/20241004_mergedFisher_DDLSvsWDLS_Bulk_keggEnrich.pdf")

# Get gene names for each KEGG pathway
kegg_enrich@result$geneID <- lapply(kegg_enrich@result$geneID, function(gene_ids) {
  # Get gene symbols from the Entrez IDs
  symbols <- genes_entrez$SYMBOL[match(unlist(strsplit(gene_ids, "/")), genes_entrez$ENTREZID)]
  paste(symbols, collapse = ", ")  # Join symbols with a comma
})

# Create the dotplot
dot_plot <- dotplot(kegg_enrich, showCategory = 10)

# Add gene names with ggrepel to avoid overlap
dot_plot + 
  geom_text_repel(aes(label = geneID), size = 3, max.overlaps = Inf, 
                  nudge_x = 0.1, nudge_y = 0.1, direction = "both", 
                  box.padding = 0.5) +
  theme(plot.margin = margin(1, 1, 1, 1, "cm"))  # Adjust margins

dev.off()

# Run GO enrichment analysis
go_enrich <- enrichGO(gene = names(gene_list_entrez), 
                      OrgDb = org.Hs.eg.db, 
                      keyType = "ENTREZID", 
                      ont = "BP",           # Ontology: BP = Biological Process
                      pvalueCutoff = 0.05)

# Plot results
write.csv2(go_enrich, "/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20241004_CellCharter_Xenium/ASPC_WDvsDD_no1_8/20241004_mergedFisher_DDLSvsWDLS_Bulk_enrichGO.csv")
pdf(file="/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20241004_CellCharter_Xenium/ASPC_WDvsDD_no1_8/20241004_mergedFisher_DDLSvsWDLS_Bulk_EnrichGO.pdf")
dotplot(go_enrich, showCategory = 10)  # Show top 10 categories
dev.off()

# Run Gene Set enrichment analysis (GSEA)

#BiocManager::install("fgsea")
#BiocManager::install("msigdbr")  # To download gene sets (e.g., KEGG, Reactome)
library(fgsea)
library(msigdbr)
pathways <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "KEGG")
pathways_list <- split(pathways$entrez_gene, pathways$gs_name)

gene_entrez <- bitr(names(gene_list), fromType = "SYMBOL", 
                    toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Rebuild the gene list using Entrez IDs
gene_list_entrez <- gene_list[gene_entrez$SYMBOL]
names(gene_list_entrez) <- gene_entrez$ENTREZID

fgsea_res <- fgsea(pathways = pathways_list, 
                   stats = gene_list_entrez,   # Named vector of log2 fold change
                   maxSize = 500)               # Number of permutations

# View results, sorted by adjusted p-value
head(fgsea_res[order(fgsea_res$padj), ])

# Prepare data for plotting
top_pathways <- fgsea_res[order(fgsea_res$padj), ][1:20, ]  # Select top 20 pathways by adjusted p-value

# Extract top 5 leading-edge genes for each pathway and concatenate them
top_pathways$leading_edge_geneNames <- sapply(top_pathways$leadingEdgeNames, function(genes) {
  paste(head(genes, 5), collapse = ", ")  # Get top 5 genes and join them with commas
})
top_pathways$leadingEdgeNames <- rep(NULL, times=nrow(top_pathways))
for (i in 1:length(top_pathways$leadingEdge)){
  top_pathways$leadingEdgeNames[[i]] <- list()
  for (j in 1:length(top_pathways$leadingEdge[[i]])){
    top_pathways$leadingEdgeNames[[i]][j]<- names(gene_list_orig)[which(names(gene_list_entrez)==top_pathways$leadingEdge[[i]][j])]
  }
}

which(names(gene_list_entrez)==top_pathways$leadingEdge[[1]][1])

# Modify the pathway names to include leading edge genes (first 5)
top_pathways$pathway_with_genes <- paste0(top_pathways$pathway, "\n(", top_pathways$leading_edge_geneNames, ")")

# Create bar plot of NES (Normalized Enrichment Score)
pdf("/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20241004_CellCharter_Xenium/ASPC_WDvsDD_no1_8/20241004_mergedFisher_DDLSvsWDLS_Bulk_GSEA.pdf", width=12)
print(ggplot(top_pathways, aes(x = reorder(pathway_with_genes, NES), y = NES)) +
  geom_bar(stat = "identity", aes(fill = padj)) +  # Fill based on -log10 of padj
  coord_flip() +  # Flip for a horizontal bar plot
  labs(x = "Pathway (with top contributing genes)", 
       y = "Normalized Enrichment Score (NES)", 
       title = "top enriched Pathways in DDLS ASPC") +
  scale_fill_gradient(low = "blue", high = "red", name = "Adjusted P-value") +  # Color gradient
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10))  # Adjust text size for clarity
)
dev.off()


#Fisher

Merged_Fisher_MeanFC <- fread('/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20241004_CellCharter_Xenium/ASPC_WDvsDD_no1_8/20241004_mergedFisher_DDLSvsWDLS_mean_log2FC.csv', sep=";")

gene_names <- Merged_Fisher_MeanFC$xG
fold_changes <- Merged_Fisher_MeanFC$mean_log2FC
gene_list <- setNames(fold_changes, gene_names)


genes_entrez <- bitr(names(gene_list), fromType = "SYMBOL", 
                     toType = "ENTREZID", OrgDb = org.Hs.eg.db)
# Merge the Entrez IDs with fold changes
gene_list_entrez <- gene_list[genes_entrez$SYMBOL]
names(gene_list_entrez) <- genes_entrez$ENTREZID
# get upregd and downreg genes
gene_list_entrez_up <- gene_list_entrez[which(gene_list_entrez>0)]
gene_list_entrez_down <- gene_list_entrez[which(gene_list_entrez<0)]

#upreg'd

kegg_enrich_up <- enrichKEGG(gene = names(gene_list_entrez),
                          organism = "hsa",     # hsa = human
                          pvalueCutoff = 0.05)
# View results
head(kegg_enrich_up)
write.csv2(kegg_enrich, "/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20241004_CellCharter_Xenium/ASPC_WDvsDD_no1_8/20241004_mergedFisher_DDLSvsWDLS_mean_log2FC_keggEnrich.csv")
pdf(file="/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20241004_CellCharter_Xenium/ASPC_WDvsDD_no1_8/20241004_mergedFisher_DDLSvsWDLS_mean_log2FC_keggEnrich.pdf")
dotplot(kegg_enrich, showCategory = 10)  # Show top 10 categories
dev.off()

# Run GO enrichment analysis
go_enrich <- enrichGO(gene = names(gene_list_entrez), 
                      OrgDb = org.Hs.eg.db, 
                      keyType = "ENTREZID",
                      ont="BP", 
                      pvalueCutoff = 0.05)

# Plot results
write.csv2(go_enrich, "/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20241004_CellCharter_Xenium/ASPC_WDvsDD_no1_8/20241004_mergedFisher_DDLSvsWDLS_mean_log2FC_enrichGO.csv")
pdf(file="/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20241004_CellCharter_Xenium/ASPC_WDvsDD_no1_8/20241004_mergedFisher_DDLSvsWDLS_mean_log2FC_EnrichGO.pdf")
dotplot(go_enrich, showCategory = 10)  # Show top 10 categories
dev.off()

# Run Gene Set enrichment analysis (GSEA)

BiocManager::install("fgsea")
BiocManager::install("msigdbr")  # To download gene sets (e.g., KEGG, Reactome)
library(fgsea)
library(msigdbr)
gene_list <- sort(gene_list, decreasing = TRUE)
pathways <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "KEGG")
pathways_list <- split(pathways$entrez_gene, pathways$gs_name)

gene_entrez <- bitr(names(gene_list), fromType = "SYMBOL", 
                    toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Rebuild the gene list using Entrez IDs
gene_list_entrez <- gene_list[gene_entrez$SYMBOL]
names(gene_list_entrez) <- gene_entrez$ENTREZID

fgsea_res <- fgsea(pathways = pathways_list, 
                   stats = gene_list_entrez,   # Named vector of log2 fold changes
                   maxSize = 500)               # Number of permutations

# View results, sorted by adjusted p-value
head(fgsea_res[order(fgsea_res$padj), ])

# Prepare data for plotting
top_pathways <- fgsea_res[order(fgsea_res$padj), ]  # Select top 20 pathways by adjusted p-value

# Extract top 5 leading-edge genes for each pathway and concatenate them
top_pathways$leading_edge_geneNames <- sapply(top_pathways$leadingEdgeNames, function(genes) {
  paste(head(genes, 5), collapse = ", ")  # Get top 5 genes and join them with commas
})
top_pathways$leadingEdgeNames <- rep(list(), times=nrow(top_pathways))
for (i in 1:length(top_pathways$leadingEdge)){
  top_pathways$leadingEdgeNames[[i]] <- list()
  for (j in 1:length(top_pathways$leadingEdge[[i]])){
    top_pathways$leadingEdgeNames[[i]][j]<- names(gene_list)[which(names(gene_list_entrez)==top_pathways$leadingEdge[[i]][j])]
  }
}

which(names(gene_list_entrez)==top_pathways$leadingEdge[[1]][1])

# Modify the pathway names to include leading edge genes (first 5)
top_pathways$pathway_with_genes <- paste0(top_pathways$pathway, "\n(", top_pathways$leading_edge_geneNames, ")")

# Create bar plot of NES (Normalized Enrichment Score)
pdf("/data/cephfs-2/unmirrored/projects/liposarcoma-wgs/20241004_CellCharter_Xenium/ASPC_WDvsDD_no1_8/20241004_mergedFisher_DDLSvsWDLS_Fischer_GSEA.pdf", width=12, height=15)
print(ggplot(top_pathways, aes(x = reorder(pathway_with_genes, NES), y = NES)) +
  geom_bar(stat = "identity", aes(fill = padj)) +  # Fill based on -log10 of padj
  coord_flip() +  # Flip for a horizontal bar plot
  labs(x = "Pathway (with top contributing genes)", 
       y = "Normalized Enrichment Score (NES)", 
       title = "top enriched Pathways in DDLS ASPC") +
  scale_fill_gradient(low = "blue", high = "red", name = "Adjusted P-value") +  # Color gradient
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10))  # Adjust text size for clarity
)
dev.off()
