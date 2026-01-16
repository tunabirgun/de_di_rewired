# 01_differential_expression.R
# Differential Gene Expression Analysis (DESeq2) Pipeline
# Includes: Normalization, DE Testing, PCA, MA Plot, Volcano Plot, Heatmap

suppressPackageStartupMessages({
    library(DESeq2)
    library(EnhancedVolcano)
    library(ggplot2)
    library(ggrepel)
    library(dplyr)
    library(openxlsx)
    library(pheatmap)
})

# --- Configuration Dictionary ---
CONFIG <- list(
    PADJ_CUTOFF = 0.05,
    LOGFC_CUTOFF = 1.0,
    MIN_COUNTS = 10,
    COLORS = list(
        UP = "#FF0000",
        DOWN = "#000080",
        NS = "grey70",
        GROUP_1 = "#0072B2", # Control
        GROUP_2 = "#D55E00"  # Treatment
    )
)

# --- Visualization Functions ---

plot_custom_pca <- function(vsd, intgroup, output_base) {
    pcaData <- plotPCA(vsd, intgroup=intgroup, returnData=TRUE)
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    
    p <- ggplot(pcaData, aes(x=PC1, y=PC2, color=group)) +
        stat_ellipse(aes(fill=group), geom="polygon", alpha=0.15, level=0.95, show.legend=FALSE) +
        geom_point(size=4, shape=21, stroke=1) +
        labs(x=paste0("PC1 (", percentVar[1], "%)"), y=paste0("PC2 (", percentVar[2], "%)")) +
        scale_color_manual(values=c(CONFIG$COLORS$GROUP_1, CONFIG$COLORS$GROUP_2)) +
        scale_fill_manual(values=c(CONFIG$COLORS$GROUP_1, CONFIG$COLORS$GROUP_2)) +
        theme_classic(base_size=14)
    
    ggsave(paste0(output_base, "_PCA.png"), p, width=6, height=6, dpi=300)
    ggsave(paste0(output_base, "_PCA.svg"), p, width=6, height=6)
}

plot_custom_volcano <- function(res, output_base) {
    df <- as.data.frame(res)
    df$sig <- "NS"
    df$sig[df$padj < CONFIG$PADJ_CUTOFF & df$log2FoldChange > CONFIG$LOGFC_CUTOFF] <- "Up"
    df$sig[df$padj < CONFIG$PADJ_CUTOFF & df$log2FoldChange < -CONFIG$LOGFC_CUTOFF] <- "Down"
    
    p <- ggplot(df, aes(x=log2FoldChange, y=-log10(padj), color=sig)) +
        geom_point(alpha=0.7, size=1.5) +
        scale_color_manual(values=c("Up"=CONFIG$COLORS$UP, "Down"=CONFIG$COLORS$DOWN, "NS"=CONFIG$COLORS$NS)) +
        geom_vline(xintercept=c(-CONFIG$LOGFC_CUTOFF, CONFIG$LOGFC_CUTOFF), linetype="dashed") +
        geom_hline(yintercept=-log10(CONFIG$PADJ_CUTOFF), linetype="dashed") +
        theme_bw() + theme(panel.grid=element_blank())
    
    ggsave(paste0(output_base, "_Volcano.png"), p, width=6, height=6, dpi=300)
    ggsave(paste0(output_base, "_Volcano.svg"), p, width=6, height=6)
}

plot_deg_heatmap <- function(vsd, res, output_base) {
    # Select Top DEGs
    sig_genes <- rownames(res)[res$padj < CONFIG$PADJ_CUTOFF & abs(res$log2FoldChange) > CONFIG$LOGFC_CUTOFF]
    top_genes <- head(sig_genes[order(res[sig_genes, "padj"])], 50)
    
    if(length(top_genes) < 2) return()
    
    mat <- assay(vsd)[top_genes, ]
    mat_z <- t(scale(t(mat)))
    
    p_heat <- pheatmap(mat_z, 
                       cluster_rows=TRUE, cluster_cols=TRUE,
                       color=colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(255),
                       border_color=NA, silent=TRUE)
    
    png(paste0(output_base, "_Heatmap.png"), width=7*300, height=6*300, res=300)
    grid::grid.draw(p_heat$gtable)
    dev.off()
    
    svg(paste0(output_base, "_Heatmap.svg"), width=7, height=6)
    grid::grid.draw(p_heat$gtable)
    dev.off()
}

# --- Main Pipeline ---

run_deseq_analysis <- function(counts_file, metadata_file, comparison_name, out_dir) {
    if(!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE)
    
    # 1. Load Data
    counts <- read.table(counts_file, header=TRUE, row.names=1, sep="\t", check.names=FALSE)
    colData <- read.table(metadata_file, header=TRUE, row.names=1)
    
    # 2. Filter & Match
    common <- intersect(rownames(colData), colnames(counts))
    counts <- counts[, common]
    colData <- colData[common, , drop=FALSE]
    
    # 3. DESeq2 Object
    dds <- DESeqDataSetFromMatrix(countData=counts, colData=colData, design=~condition)
    dds <- dds[rowSums(counts(dds)) >= CONFIG$MIN_COUNTS, ]
    
    # 4. Run Analysis
    dds <- DESeq(dds)
    res <- results(dds)
    resLFC <- lfcShrink(dds, coef=resultsNames(dds)[2], type="apeglm") # Auto-select coef 2
    
    # 5. Save Tables
    res_df <- as.data.frame(resLFC) %>% arrange(padj)
    write.xlsx(res_df, file.path(out_dir, paste0(comparison_name, "_All_Genes.xlsx")), rowNames=TRUE)
    
    up <- res_df %>% filter(padj < CONFIG$PADJ_CUTOFF & log2FoldChange > CONFIG$LOGFC_CUTOFF)
    down <- res_df %>% filter(padj < CONFIG$PADJ_CUTOFF & log2FoldChange < -CONFIG$LOGFC_CUTOFF)
    
    write.xlsx(up, file.path(out_dir, paste0(comparison_name, "_Upregulated.xlsx")), rowNames=TRUE)
    write.xlsx(down, file.path(out_dir, paste0(comparison_name, "_Downregulated.xlsx")), rowNames=TRUE)
    
    # 6. Dimensions Reduction & Plots
    vsd <- vst(dds, blind=FALSE)
    
    plot_custom_pca(vsd, "condition", file.path(out_dir, comparison_name))
    plot_custom_volcano(res_df, file.path(out_dir, comparison_name))
    plot_deg_heatmap(vsd, res_df, file.path(out_dir, comparison_name))
}
