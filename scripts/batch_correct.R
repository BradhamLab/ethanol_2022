.libPaths("/projectnb2/bradham/RPackages/4.0.2")
# .libPaths("/projectnb2/bradham/RPackages/3.6.2")
suppressMessages(library(sva))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(pheatmap))
library(grid)


theme_black <- function(base_size = 12, base_family = "") {
  theme_grey(base_size = base_size, base_family = base_family) %+replace%

    theme(
      # Specify axis options
      axis.line = element_blank(),
      axis.text.x = element_text(size = base_size * 0.8, color = "white", lineheight = 0.9),
      axis.text.y = element_text(size = base_size * 0.8, color = "white", lineheight = 0.9),
      axis.ticks = element_line(color = "white", size = 0.2),
      axis.title.x = element_text(size = base_size, color = "white", margin = margin(0, 10, 0, 0)),
      axis.title.y = element_text(size = base_size, color = "white", angle = 90, margin = margin(0, 10, 0, 0)),
      axis.ticks.length = unit(0.3, "lines"),
      # Specify legend options
      legend.background = element_rect(color = "NA", fill = "black"),
      legend.key = element_rect(color = "white", fill = "black"),
      legend.key.size = unit(1.2, "lines"),
      legend.key.height = NULL,
      legend.key.width = NULL,
      legend.text = element_text(size = base_size * 0.8, color = "white"),
      legend.title = element_text(size = base_size * 0.8, face = "bold", hjust = 0, color = "white"),
      legend.position = "right",
      legend.text.align = NULL,
      legend.title.align = NULL,
      legend.direction = "vertical",
      legend.box = NULL,

      # Specify panel options
      panel.background = element_rect(fill = "black", color = NA),
      panel.border = element_rect(fill = NA, color = "white"),
      panel.grid.major = element_line(color = "grey35"),
      panel.grid.minor = element_line(color = "grey20"),
      #   panel.spacing = unit(0.5, "lines"),
      # Specify facetting options
      strip.background = element_rect(fill = "grey30", color = "grey10"),
      strip.text.x = element_text(size = base_size * 0.8, color = "white"),
      strip.text.y = element_text(size = base_size * 0.8, color = "white", angle = -90),
      # Specify plot options
      plot.background = element_rect(color = "black", fill = "black"),
      plot.title = element_text(size = base_size * 1.2, color = "white"),
      plot.margin = unit(rep(1, 4), "lines")
    )
}

if (exists("snakemake")) {
  samples <- read.csv(snakemake@input$samples, row.names = 1, check.names = FALSE)
  counts <- read.csv(snakemake@input$counts, row.names = 1, check.names = FALSE)
  covar_mat <- do.call(cbind, lapply(c("HPF", "Treatment"), function(x) {
    return(as.integer(as.factor(samples[, x])))
  }))
  batch <- as.integer(as.factor(samples[, snakemake@params$batch]))
  print("Batch correcting...")
  # corrected <- ComBat_seq(as.matrix(counts),
  #                         batch=batch,
  #                         group=NULL,
  #                         covar_mod=covar_mat)
  corrected <- ComBat(as.matrix(counts),
    batch = batch
  )

  print("Performing PCA...")
  pca <- prcomp(t(corrected))

  samples$PC1 <- pca$x[, 1]
  samples$PC2 <- pca$x[, 2]
  variance <- pca$sdev^2 / sum(pca$sdev^2)
  colors <- RColorBrewer::brewer.pal(6, "Paired")
  names(colors) <- c(
    "Ethanol-15", "Control-15", "Ethanol-21", "Control-21", "Ethanol-18", "Control-18"
  )
  ggplot(data = samples, aes(x = PC1, y = PC2, color = Group)) +
    geom_point(size = 4) +
    scale_color_manual(values = colors) +
    xlab(paste("PC1: ", round(variance[1], 2) * 100, "% variance")) +
    ylab(paste("PC2: ", round(variance[2], 2) * 100, "% variance")) + 
    theme_black()
  ggsave(snakemake@output$pca, width = 4.25, height = 3, units = "in")

  # calculate distance off of PCA
  print("Calculating distances...")
  dmat <- dist(pca$x)
  colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
  for (each in snakemake@params$cofactors) {
    samples[, each] <- as.factor(samples[, each])
  }
  p <- pheatmap(dmat,
    col = colors,
    annotation_col = samples[snakemake@params$cofactors]
  )
  save_pheatmap_png <- function(x, filename, width = 1200, height = 1000, res = 150) {
    png(filename, width = width, height = height, res = res)
    grid::grid.newpage()
    grid.draw(rectGrob(gp = gpar(fill = "black", lwd = 0)))
    grid::grid.draw(x$gtable)
    grid.gedit("layout", gp = gpar(col = "white", text = ""))
    dev.off()
  }

  save_pheatmap_png(p, snakemake@output$dist)
  write.csv(corrected, snakemake@output$corrected)
}