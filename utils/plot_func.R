gsea.barplot <- function(gse, showCategory=10){
  res <- gse@result
  res$generatio <- sapply(res$core_enrichment, function(x) {length(strsplit(x, "/")[[1]])}, USE.NAMES = FALSE) / res$setSize
  df <- res %>% 
    group_by(NES < 0) %>% 
    slice_min(n=showCategory, order_by=qvalue, with_ties = FALSE) %>% 
    arrange(order_by=NES)
  df$Description <- factor(df$Description, levels=df$Description)
  
  ggplot(df, aes(y=factor(Description), x=NES, fill=generatio)) + 
    geom_bar(stat="identity") +
    theme_classic() +
    theme(axis.title.y = element_blank(), axis.text.y = element_text(size=12)) + 
    xlab("NES") +
    # scale_fill_manual(values=c("red","blue"))
    scale_fill_gradient(low = "blue", high="red")
}

plot.bayesspace <- function(sce.enhanced, sce, features, ncol=4, compare=FALSE) {
  enhanced.plots <- purrr::map(features, function(x) featurePlot(sce.enhanced, x))
  spot.plots <- purrr::map(features, function(x) featurePlot(sce, x))
  # options(repr.plot.height=15, repr.plot.width=20)
  
  if(!compare){
    patchwork::wrap_plots(c(enhanced.plots), ncol=ncol)
  } else {
    patchwork::wrap_plots(c(enhanced.plots, spot.plots), ncol=ncol)
  }
}