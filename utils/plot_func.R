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

fgsea.barplot <- function(fgseares, showCategory=10, label.size=4, color="pval"){
  res <- fgseares
  if(class(res$leadingEdge)=="character"){
    res$generatio <- sapply(strsplit(res$leadingEdge,","), length) / res$size
  } else { # list
    res$generatio <- sapply(res$leadingEdge, length) / res$size
  }
  df <- res %>% 
    filter(!is.na(NES)) %>%
    group_by(NES < 0) %>% 
    slice_max(n=showCategory, order_by=abs(NES), with_ties = FALSE) %>% 
    arrange(order_by=NES)
  df$pathway <- factor(df$pathway, levels=df$pathway)
  x.max <- max(abs(df$NES))
  ggplot(df, aes(y=factor(pathway), x=NES, fill=get(color))) + 
    geom_text(aes(label = pathway, x = ifelse(NES > 0, -0.05, 0.05), group=`NES < 0`, hjust=ifelse(NES > 0, 1,0)), position = position_dodge(width = 1), size = label.size) +
    geom_col() +
    geom_segment(aes(x=0, y=pathway, xend=ifelse(NES > 0, -0.025, 0.025), yend=pathway), size=0.5, color="black") + 
    theme_classic() +
    theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
          axis.line.y.left = element_blank(), axis.ticks.y = element_blank()) + 
    xlab("NES") +
    geom_vline(xintercept = 0) + 
    scale_fill_gradient(color, low = "green", high="red") +
    scale_x_continuous(limits = c(-x.max, x.max))
}