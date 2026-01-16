# 04_enrichment.R
# Functional Enrichment Analysis of Rewired Proteins

library(clusterProfiler)
library(enrichplot)
library(ggplot2)

# --- Functions ---

perform_enrichment <- function(gene_list, organism_db, title) {
    if(length(gene_list) < 5) return(NULL)
    
    # Run GO Enrichment
    ego <- enrichGO(
        gene = gene_list,
        OrgDb = organism_db,
        keyType = "GID",  # Adjust based on your ID format
        ont = "ALL",      # BP, MF, CC
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2
    )
    
    if(is.null(ego) || nrow(ego) == 0) return(NULL)
    
    # Generate Plot
    p <- dotplot(ego, showCategory=15) + ggtitle(paste("GO Enrichment:", title))
    
    return(list(
        result = as.data.frame(ego),
        plot = p
    ))
}

# --- Main Execution Template ---

run_enrichment_step <- function(rewired_summary_file, output_dir) {
    rewired_data <- read.csv(rewired_summary_file)
    # db <- AnnotationHub()[["AH12345"]] # Example loading
    
    # Analyze
    # res <- perform_enrichment(rewired_data$GeneID, db, "My Study")
    
    # Save
    # write.csv(res$result, file.path(output_dir, "GO_Results.csv"))
    # ggsave(file.path(output_dir, "GO_Dotplot.png"), res$plot)
}
