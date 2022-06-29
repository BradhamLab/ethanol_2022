# .libPaths("/projectnb2/bradham/RPackages/3.6.2")
suppressMessages(library(GSVA))
suppressMessages(library(rjson))
suppressMessages(library(Biobase))
suppressMessages(library(RColorBrewer))
suppressMessages(library(pheatmap))
library(grid)


#' Transform gene expression data to pathway enrichment. Perform
#' One vs all wilcoxon-sign rank test to determine enrichment
#' between clusters.
#'
#' @param X data matrix of gene expression values.
#' @param obs observation level data with cell clusters.
#' @param genelist list of genesets containing gene identifiers for pathways
#'
#' @return X: pathway matrix, results: one-vs-all results with p-values + fold changes
#' @export
#'
#' @examples
pathway_enrichment <- function(X, obs, genelists, cluster = "Cluster") {
  if (!all(rownames(obs) == colnames(X))) {
    stop("Row and column names do not match between expression and observation data")
  }
  pdata <- Biobase::AnnotatedDataFrame(obs)
  eset <- Biobase::ExpressionSet(as.matrix(X), phenoData = pdata)

  gsva_out <- GSVA::gsva(eset, genelists)
  clusters <- pData(gsva_out)[, cluster]
  one_vs_all <- do.call(rbind, lapply(unique(clusters), function(x) {
    cluster_cells <- exprs(gsva_out[, clusters == x])
    not_cluster <- exprs(gsva_out[, clusters != x])
    go.out <- sapply(row.names(cluster_cells), function(go) {
      wilcox.out <- wilcox.test(
        cluster_cells[go, ],
        not_cluster[go, ]
      )
      mean_enrich <- as.numeric(mean(cluster_cells[go, ]))
      outlist <- list(wilcox.out$p.value, mean_enrich)
      names(outlist) <- c(
        paste0(x, ".vs.all.pvalue"),
        paste0(x, ".mean.enrichment")
      )
      return(outlist)
    })
    return(go.out)
  }))

  one_vs_all <- as.data.frame(t(one_vs_all))
  for (i in unique(clusters)) {
    pval <- paste0(i, ".vs.all.pvalue")
    one_vs_all[, paste0(i, ".vs.all.p.adjust")] <- p.adjust(one_vs_all[, pval])
  }
  one_vs_all <- one_vs_all[sort(colnames(one_vs_all))]
  one_vs_all <- as.data.frame(apply(one_vs_all, 2, as.numeric))
  row.names(one_vs_all) <- row.names(gsva_out)
  return(list("X" = exprs(gsva_out), "results" = one_vs_all))
}


# main <- function(X_csv, obs_csv, genelist_json, out_X, out_results,
#                  cluster='Cluster') {
#   X <- read.table(X_csv, sep=',',
#                   header=TRUE, row.names=1, check.names=FALSE)
#   obs <- read.table(obs_csv, row.names=1, sep=',', header=TRUE, check.names=FALSE)
#   genelists <- rjson::fromJSON(file=genelist_json)
#   results <- pathway_enrichment(X, obs, genelists, cluster)
#   return(results)
# #   write.csv(results$X, out_X)
# #   write.csv(results$results, out_results)
#   }

# setwd('~/Data/scratch/single-cell/w_chlorate/annotated/')
# main("collapsed_X.csv", "obs.csv", "../../../../go_gene_list.json",
#      "../gsva_X.csv", "../gsva_wilcox.csv")

# setwd('~/Data/scratch/single-cell/w_chlorate/treated')

if (exists("snakemake")) {
  counts <- read.csv(snakemake@input$corrected, row.names = 1, check.names = FALSE)
  annos <- read.csv(snakemake@input$annos, row.names = 1, check.names = FALSE)
  keep <- row.names(counts) %in% annos[, snakemake@params$uniprot]
  fcounts <- counts[keep, ]
  coldata <- read.csv(snakemake@input$samples, row.names = 1, check.names = FALSE)
  coldata$HPF <- as.factor(coldata$HPF)
  go_lookup <- read.csv(
    snakemake@input$golookup,
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  genelists <- rjson::fromJSON(file = snakemake@input$genesets)
  comparisons <- snakemake@params$comparisons
  top_pathways <- list(length(seq(1, length(comparisons), 2)))
  i <- 1
  for (idx in seq(1, length(comparisons), 2)) {
    each <- comparisons[idx:(idx + 1)]
    print(paste0("Comparison: ", paste(each, collapse = " vs ")))
    samples <- row.names(coldata)[apply(coldata, 1, function(x) {
      x[snakemake@params$group] %in% each
    })]
    print(paste0("Samples ", paste(samples, collapse = ", ")))
    sub_counts <- fcounts[, samples]
    sub_coldata <- coldata[samples, ]
    results <- pathway_enrichment(
      sub_counts, sub_coldata,
      genelists, snakemake@params$compare
    )
    results$results$desc <- go_lookup[row.names(results$results), "desc"]
    prefix <- paste(each, collapse = "_vs_")
    write.csv(
      results$X,
      file.path(snakemake@params$outdir, paste0(prefix, "_gsva.csv"))
    )
    write.csv(
      results$results,
      file.path(snakemake@params$outdir, paste0(prefix, "_results.csv"))
    )

    # let's extract + plot top pathways
    ordered_path <- row.names(results$results)[order(results$results[, "Control.mean.enrichment"])]
    most_down <- ordered_path[1:20]
    most_up <- ordered_path[(nrow(results$results) - 19):nrow(results$results)]
    most_diff <- c(most_down, most_up)
    plot_mat <- results$X[most_diff, ]
    print(paste0("length of up-down: ", length(most_diff)))
    print(length(most_down))
    print(length(most_up))
    top <- data.frame(
      GO = most_diff,
      Comparison = rep(prefix, length(most_diff)),
      In.Control = c(rep("-", length(most_down)), rep("+", length(most_up))),
      Desc = go_lookup[most_diff, 1]
    )
    top_pathways[[i]] <- top
    i <- i + 1
    # # TODO: change hard coding
    colnames(plot_mat) <- paste(sub_coldata$Group, sub_coldata$Sample,
      sep = "-"
    )

    row.names(sub_coldata) <- colnames(plot_mat)
    # colors <- colorRampPalette(rev(brewer.pal(9, "PiYG")) )(255)
    new_row <- lapply(rownames(plot_mat), function(x) bquote(bold(.(x))))
    new_col <- lapply(colnames(plot_mat), function(x) bquote(bold(.(x))))

    p <- pheatmap(plot_mat,
      # col=colors,
      annotation_col = sub_coldata[c("HPF", "Treatment", "Sample")],
      labels_col = as.expression(new_col),
      fontsize_row = 8,
      show_rownames = TRUE,
      fontsize_col = 14,
      fontsize = 16
    )

    save_pheatmap_png <- function(x, filename, width = 1200, height = 1000, res = 150) {
      png(filename, width = width, height = height, res = res)
      grid::grid.newpage()
      # grid.draw(rectGrob(gp=gpar(fill="black", lwd=0)))
      grid::grid.draw(x$gtable)
      # grid.gedit("layout", gp = gpar(col = "white", text = ""))
      dev.off()
    }

    save_pheatmap_png(p, file.path(
      snakemake@params$plotdir,
      paste0(prefix, "_top_enrichment.png")
    ))
  }
  top_results <- do.call(rbind, top_pathways)
  write.csv(top_results, snakemake@output$top)
}