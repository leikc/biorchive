run.gsea <- function(markers.df, org=org.Hs.eg.db, keytype="SYMBOL", ont="BP", padj="none", seed=100){
  orig.gene.list <- markers.df$avg_log2FC
  names(orig.gene.list) <- rownames(markers.df)
  gene.list<-na.omit(orig.gene.list)
  gene.list = sort(gene.list, decreasing = TRUE)
  
  gse <- gseGO(geneList=gene.list, 
               ont = ont, 
               keyType = keytype, 
               minGSSize = 3, 
               maxGSSize = 800, 
               pvalueCutoff = 0.05, 
               verbose = TRUE, 
               OrgDb = org, 
               pAdjustMethod = padj,
               seed = seed)
  
  return(list(gse=gse, gene.list=gene.list))
}

run.kegg <- function(markers.df, org=org.Hs.eg.db, keytype="SYMBOL", kegg_organism = "hsa",
                     padj="none", seed=100){
  orig.gene.list <- markers.df$avg_log2FC
  names(orig.gene.list) <- rownames(markers.df)
  gene.list<-na.omit(orig.gene.list)
  gene.list = sort(gene.list, decreasing = TRUE)
  
  ids<-bitr(names(gene.list), fromType = keytype, toType = "ENTREZID", OrgDb=org)
  dedup_ids = ids[!duplicated(ids[c(keytype)]),]
  df2 = markers.df[rownames(markers.df) %in% dedup_ids$SYMBOL,]
  df2 <- df2[dedup_ids$SYMBOL,]
  df2$Y = dedup_ids$ENTREZID
  
  kegg_gene_list <- df2$avg_log2FC
  names(kegg_gene_list) <- df2$Y
  kegg_gene_list<-na.omit(kegg_gene_list)
  kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
  
  R.utils::setOption("clusterProfiler.download.method","auto")
  kk2 <- gseKEGG(geneList=kegg_gene_list, use_internal_data = TRUE,
                 organism=kegg_organism,
                 minGSSize=3,
                 maxGSSize=800,
                 pvalueCutoff=0.05,
                 pAdjustMethod=padj,
                 keyType="ncbi-geneid", 
                 seed = seed) # No results with adjusted p, use no adjustment
  
  return(list(kegg=kk2, gene.list=kegg_gene_list, dedup.ids=dedup_ids))
}

run.monocle3 <- function(so, assay.name="SCT", ndim=50){
  cds <- new_cell_data_set(expression_data = so[[assay.name]]@data, 
                           cell_metadata = so@meta.data, 
                           gene_metadata = data.frame(gene_short_name=rownames(so[[assay.name]]@counts), row.names=rownames(so[[assay.name]]@counts)))
  cds <- preprocess_cds(cds, num_dim = ndim, norm_method = "none", method = "PCA")
  cds <- reduce_dimension(cds)
  cds <- cluster_cells(cds)
  cds <- learn_graph(cds)
  cds
}