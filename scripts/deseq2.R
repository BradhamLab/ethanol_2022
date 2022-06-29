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

treatment_de <- function(counts, coldata, treatment='Treatment',
                         control='Control', batch='Batch') {
    # set control value as 
    levels = c(control, setdiff(unique(coldata[, treatment]),
                                c(control)))
    coldata[, treatment] <- factor(coldata[, treatment],
                                   levels=levels)
    # compare between treatments
    formula_str <- paste0('~ ', treatment)
    if (!is.null(batch) || batch != '') {
        formula_str <- paste0(formula_str, ' + ', batch)
    }
    print(paste0('Formula: ', formula_str))
    design = as.formula(formula_str)
    dds <- DESeqDataSetFromMatrix(countData=counts,
                                  colData=coldata,
                                  design=design)
    dds <- DESeq(dds)
    res <- results(dds)
    # plot DE genes
    p <- volcano_plot(data.frame(res), dark=FALSE)
    return(list('res'=data.frame(res), 'plot'=p))
}

volcano_plot <- function(res, fc_threshold=1, p_threshold=0.05,
                         title="", dark=TRUE) {
    ggplot(data=res, aes(x='log2FoldChange', y='pvalue'))
    res$log10p <- -log10(res$pvalue)
    res$significant <- apply(res, 1, function(x) {
        if (is.na(x['padj']) || is.na(x['log2FoldChange'])) {
            return(FALSE)
        }
        if (x['padj'] <= p_threshold && abs(x['log2FoldChange']) >= fc_threshold) {
            return(TRUE)
        }
        return(FALSE)
    })
    y_height <- min(res$log10p[res$significant], na.rm=TRUE)
    xticks1 <- seq(0, max(ceiling(res$log2FoldChange), na.rm=TRUE), 2)
    xticks2 <- seq(0, min(floor(res$log2FoldChange), na.rm=TRUE), -2)
    xticks <- sort(c(-fc_threshold, fc_threshold, xticks1, xticks2))
    # yticks <- sort(c(seq(0, 1, 0.25), p_threshold))
    p <- ggplot(data=res, aes(x=log2FoldChange, y=log10p, col=significant)) +
            geom_point(alpha=0.75) +
            geom_vline(xintercept=fc_threshold, linetype='dotted') +
            geom_vline(xintercept=-fc_threshold, linetype='dotted') +
            geom_hline(yintercept=y_height, linetype='dotted') +
            xlab('Log2(FoldChange)') +
            ylab('-Log10(p)') +
            # scale_color_manual(values=c('black', 'red')) +
            scale_x_continuous(breaks=xticks) +
            ggtitle(title) + 
            theme(legend.position = 'none',
                  axis.text=element_text(size=14, face='bold'),
                  axis.title=element_text(size=18, face="bold"))
    if (dark) {
        p <- p +
            theme_black(22) +
            scale_color_manual(values=c('gray', '#ff6ec7')) 
            theme(legend.position = 'none',
                  axis.text=element_text(size=14, face='bold'),
                  axis.title=element_text(size=18, face="bold"))
    } else {
        p <- p + scale_color_manual(values=c('black', 'red')) 
    }
    return(p)
}

main <- function(counts, coldata, comparisons, col, outdir,
                 plotdir, treatment='Treatment', control='Control',
                 batch='Batch') {
    for (idx in seq(1, length(comparisons), 2)) {
        each <- comparisons[idx:(idx + 1)]
        print(paste0("Comparison: ", paste(each, collapse=' vs ')))
        samples <- row.names(coldata)[apply(coldata, 1, function(x) x[col] %in% each)]
        print(paste0("Samples ", paste(samples, collapse=', ')))
        sub_counts <- counts[ , samples]
        sub_coldata <- coldata[samples, ]
        out <- treatment_de(sub_counts, sub_coldata, treatment, control, batch)
        write.csv(out$res, file.path(outdir,
                                     paste0(paste(each, collapse='_vs_'), '.csv')))
        out$plot + ggtitle(paste(each, collapse=' vs '))
        ggsave(file.path(plotdir, paste0(paste(each, collapse='_vs_'), '.png')),
               width=9, height=7, dpi='retina')
    }
}

if (exists("snakemake")) {
    coldata <- read.csv(snakemake@input$samples, row.names=1, check.names=FALSE)
    coldata$HPF <- as.factor(coldata$HPF)
    # load counts and match sample order
    counts <- read.csv(snakemake@input$counts,
                       row.names=1,
                       check.names=FALSE)[ , row.names(coldata)]
    # run deseq2 on desired comparisons
    print(snakemake@params$comparisons)
    main(counts, coldata, snakemake@params$comparisons, snakemake@params$col,
         snakemake@params$outdir, snakemake@params$plotdir,
         snakemake@params$treatment, snakemake@params$control,
         snakemake@params$batch)

    # plot pca of all samples
    design = as.formula(paste0('~ ', snakemake@params$treatment))
    dds <- DESeqDataSetFromMatrix(countData=counts,
                                  colData=coldata,
                                  design=design)
    vsd <- vst(dds, blind=FALSE)
    print("Plotting Principle Coordinates")
    plotPCA(vsd, intgroup=c("Treatment", "HPF")) + theme_black(14)
    ggsave(file.path(snakemake@params$plotdir, 'treatment-sample-pca.png'),
           width=9, height=7, dpi='retina')


    # plot dissimilarity matrix between samples
    print("Calculating correlations")
    sample_corr <- 1 - cor(assay(vsd))
    # sample_corr_mat <- as.matrix(sample_corr)
    rownames(sample_corr) <- paste(vsd$Group, vsd$Sample, sep="-")
    colnames(sample_corr) <- paste(vsd$Group, vsd$Sample, sep="-")
    row.names(coldata) <- row.names(sample_corr)
    colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

    p <- pheatmap(sample_corr,
                  col=colors,
                  annotation_col=coldata[c('HPF', 'Treatment', 'Sample')])
    save_pheatmap_png <- function(x, filename, width=1200, height=1000, res = 150) {
        png(filename, width = width, height = height, res = res)
        grid::grid.newpage()
        grid.draw(rectGrob(gp=gpar(fill="black", lwd=0)))
        grid::grid.draw(x$gtable)
        grid.gedit("layout", gp = gpar(col = "white", text = ""))
        dev.off()
    }

    save_pheatmap_png(p, file.path(snakemake@params$plotdir, 'sample-diss-matrix.png'))
    warnings()
    write.csv(assay(vsd), snakemake@output$counts)
}