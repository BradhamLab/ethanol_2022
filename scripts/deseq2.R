# DeSeq2 is broken under bradham project libraries -- unsure what fixing will do
if (length(.libPaths()) > 1) {
    .libPaths(.libPaths()[2])
}
suppressMessages(library(DESeq2))
suppressMessages(library(ggplot2))
suppressMessages(library(pheatmap))
suppressMessages(library(RColorBrewer))
suppressMessages(library(gridExtra))
library(grid)
library(dplyr)
library(RColorBrewer)

theme_black = function(base_size = 12, base_family = "") {
  
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    
    theme(
      # Specify axis options
      axis.line = element_blank(),  
      axis.text.x = element_text(size = base_size*0.8, color = "white", lineheight = 0.9),  
      axis.text.y = element_text(size = base_size*0.8, color = "white", lineheight = 0.9),  
      axis.ticks = element_line(color = "white", size  =  0.2),  
      axis.title.x = element_text(size = base_size, color = "white", margin = margin(0, 10, 0, 0)),  
      axis.title.y = element_text(size = base_size, color = "white", angle = 90, margin = margin(0, 10, 0, 0)),  
      axis.ticks.length = unit(0.3, "lines"),   
      # Specify legend options
      legend.background = element_rect(color = NA, fill = "black"),  
      legend.key = element_rect(color = "white",  fill = "black"),  
      legend.key.size = unit(1.2, "lines"),  
      legend.key.height = NULL,  
      legend.key.width = NULL,      
      legend.text = element_text(size = base_size*0.8, color = "white"),  
      legend.title = element_text(size = base_size*0.8, face = "bold", hjust = 0, color = "white"),  
      legend.position = "right",  
      legend.text.align = NULL,  
      legend.title.align = NULL,  
      legend.direction = "vertical",  
      legend.box = NULL, 
      # Specify panel options
      panel.background = element_rect(fill = "black", color  =  NA),  
      panel.border = element_rect(fill = NA, color = "white"),  
      panel.grid.major = element_line(color = "grey35"),  
      panel.grid.minor = element_line(color = "grey20"),  
    #   panel.spacing = unit(0.5, "lines"),   
      # Specify facetting options
      strip.background = element_rect(fill = "grey30", color = "grey10"),  
      strip.text.x = element_text(size = base_size*0.8, color = "white"),  
      strip.text.y = element_text(size = base_size*0.8, color = "white",angle = -90),  
      # Specify plot options
      plot.background = element_rect(color = "black", fill = "black"),  
      plot.title = element_text(size = base_size*1.2, color = "white"),  
      plot.margin = unit(rep(1, 4), "lines")
      
    )
  
}


volcano_plot <- function(res, fc_threshold=1, p_threshold=0.05,
                         title="", dark=TRUE) {
    
    # res <- as.data.frame(res) %>%
    #     dplyr::mutate(
    #         Significant = dplyr::case_when(
    #             !is.na(padj) & padj < p_threshold & log2FoldChange > fc_threshold ~ "+",
    #             !is.na(padj) & padj < p_threshold & log2FoldChange < -fc_threshold ~ '-',
    #             TRUE ~ "n.s."
    #         ),
    #         log10p = -log10(pvalue),
    #         Significant = factor(Significant, levels=c("+", "-", "n.s."))
    #     )

    y_height <- min(res$log10p[res$Significant], na.rm=TRUE)
    xticks1 <- seq(0, max(ceiling(res$log2FoldChange), na.rm=TRUE), 2)
    xticks2 <- seq(0, min(floor(res$log2FoldChange), na.rm=TRUE), -2)
    xticks <- sort(c(-fc_threshold, fc_threshold, xticks1, xticks2))
    # yticks <- sort(c(seq(0, 1, 0.25), p_threshold))
    p <- ggplot(data=res, aes(x=log2FoldChange, y=log10p, color=Significant)) +
            geom_point(alpha=0.75) +
            geom_vline(xintercept=fc_threshold, linetype='dotted') +
            geom_vline(xintercept=-fc_threshold, linetype='dotted') +
            geom_hline(yintercept=y_height, linetype='dotted') +
            xlab("Log2FC") +
            ylab('-Log10(p)') +
            scale_x_continuous(breaks=xticks) +
            ggtitle(title) + 
            scale_color_manual(
                values = c("+" = "#EF8A62", "-" = "#67A9CF", "n.s." = "#D3D3D3")
            ) +
            theme(legend.position = 'none',
                  axis.text=element_text(size=14, face='bold'),
                  axis.title=element_text(size=18, face="bold"))
    if (dark) {
        p <- p +
            theme_black(22) +
            scale_color_manual(values=c('#ff6ec7', "green", 'gray')) 
            theme(legend.position = 'none',
                  axis.text=element_text(size=14, face='bold'),
                  axis.title=element_text(size=18, face="bold"))
    }
    return(p)
}


get_vsd_and_plot <- function(dds, plotdir, group) {
    vsd <- vst(dds, blind=FALSE)
    print("Plotting Principle Coordinates")
    plotPCA(vsd, intgroup=group)
    ggsave(
        file.path(plotdir, 'treatment-sample-pca.png'),
        width=9, height=7, dpi='retina'
    )

    # # plot dissimilarity matrix between samples
    # print("Calculating correlations")
    # sample_corr <- 1 - cor(assay(vsd))
    # # sample_corr_mat <- as.matrix(sample_corr)
    # rownames(sample_corr) <- colData(vsd)[, group]
    # colnames(sample_corr) <- colData(vsd)[, group]
    # colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

    # p <- pheatmap(
    #     sample_corr,
    #     col=colors,
    #     annotation_col=as.data.frame(colData(vsd))[c('HPF', 'Treatment', group)]
    # )
    # save_pheatmap_png <- function(x, filename, width=1200, height=1000, res = 150) {
    #     png(filename, width = width, height = height, res = res)
    #     grid::grid.newpage()
    #     # grid.draw(rectGrob(gp=gpar(fill="black", lwd=0)))
    #     grid::grid.draw(x$gtable)
    #     # grid.gedit("layout", gp = gpar(col = "white", text = ""))
    #     dev.off()
    # }
    # png(
    #     file.path(plotdir, "sample-diss-matrix.png"),
    #     width=1200,
    #     height=1000,
    #     res=150
    # )
    # dev.off()
    # save_pheatmap_png(p, file.path(plotdir, 'sample-diss-matrix.png'))
    warnings()
    return(assay(vsd))
}


treatment_de <- function(counts, coldata, comparison_col, reference, experimental, name, plotdir, batch='Batch', p_threshold=0.05, fc_threshold=1) {
    if (!dir.exists(plotdir)) {
        dir.create(plotdir)
    }

    # compare between treatments
    formula_str <- paste0('~ ', comparison_col)
    if (!is.null(batch)) {
        formula_str <- paste0(formula_str, ' + ', batch)
    }
    print(paste0('Formula: ', formula_str))
    design = as.formula(formula_str)
    dds <- DESeqDataSetFromMatrix(countData=counts,
                                  colData=coldata,
                                  design=design)
    dds <- DESeq(dds)
    print("PLATDOAIDFAD")
    print(plotdir)
    saveRDS(dds, 'deseq_dds.RDS')
    deseq_results <- lapply(experimental, function(condition) {
        condition_results <- results(dds, contrast=c(comparison_col, condition, reference))
        condition_results <- data.frame(condition_results)
        condition_results$Gene <- row.names(condition_results)
        condition_results$Comparision <- paste0(reference, "_vs_", condition)

        condition_results <- condition_results %>%
            dplyr::mutate(
                Significant = dplyr::case_when(
                    !is.na(padj) & padj < p_threshold & log2FoldChange > fc_threshold ~ "+",
                    !is.na(padj) & padj < p_threshold & log2FoldChange < -fc_threshold ~ '-',
                    TRUE ~ "n.s."
                ),
                log10p = -log10(pvalue),
                Significant = factor(Significant, levels=c("+", "-", "n.s."))
            )

        p <- volcano_plot(condition_results, dark=FALSE) +
            ggtitle(paste0(reference, "_vs_", condition))
        ggsave(
            file.path(plotdir, paste0(reference, '_vs_', condition, '.png')),
            plot = p,
            width=9, height=7, dpi='retina'
        )
        
        return(condition_results)
    })
    vsd <- get_vsd_and_plot(
        dds,
        plotdir,
        comparison_col
    )
    return(list('normalized' = vsd, 'results'=dplyr::bind_rows(deseq_results)))
}

main <- function(counts, coldata, comparisons, outdir,
                 plotdir,
                 batch=NULL) {
    out_list <- list()
    for (each in names(comparisons)) {
        print(each)
        # get comparison info
        column <- comparisons[[each]][["column"]]
        reference <- comparisons[[each]][["reference"]]
        experimental <- comparisons[[each]][["experimental"]]
        print(column)
        for (every in experimental) {
            print(every)
            if (!every %in% coldata[, column]) {
                stop(sprintf("%s not found in column %s", every, column))
            }
        }
        # subset data to samples of interest
        to_test <- c(reference, experimental)
        print(paste0("Comparison: ", each))
        samples <- row.names(coldata)[apply(coldata, 1, function(x) x[column] %in% to_test)]
        print(paste0("Samples ", paste(samples, collapse=', ')))
        sub_counts <- counts[ , samples]
        sub_coldata <- coldata[samples, ]
        # set reference values
        print("ORIGINAL")
        print(plotdir)
        out <- treatment_de(
            sub_counts,
            sub_coldata,
            column,
            reference,
            experimental,
            each,
            plotdir=file.path(plotdir, each, "plots", "deseq2"),
            batch
        )
        print("WRITING HERE:            !!!!!")
        print(outdir)
        write.csv(
            out$results,
            file.path(outdir, each, "deseq2", "results.csv")
        )
        write.csv(
            out$normalized,
            file.path(outdir, each, "deseq2", "normalized_counts.csv")
        )
        print("wrote all output")       
    }
}


if (exists("snakemake")) {
    coldata <- read.csv(snakemake@input$samples, row.names=1, check.names=FALSE)
    if (!length(intersect(c("Sample", "HPF", "Treatment"), colnames(coldata)))) {
        stop("Expected Sample and HPF columns. This is case sensitive.")
    }
    coldata$HPF <- as.factor(coldata$HPF)
    # load counts and match sample order
    counts <- read.csv(snakemake@input$counts,
                       row.names=1,
                       check.names=FALSE)
    # not ideal, fix later or not
    colnames(counts) <- gsub("^sample-", "", colnames(counts))
    counts <- counts[ , row.names(coldata)]
    # run deseq2 on desired comparisons
    print(snakemake@params$comparisons)
    main(
        counts,
        coldata,
        snakemake@params$comparisons,
        snakemake@params$outdir,
        snakemake@params$plotdir,
        snakemake@params$batch
    )    
}