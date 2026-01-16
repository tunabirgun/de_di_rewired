# 03_rewired_analysis.R
# Extraction of "Rewired" Proteins (High Connectivity Change, Stable Expression)

library(dplyr)
library(AnnotationHub)

# --- Parameters ---
DEG_PERCENTILE <- 0.90   # Threshold for High Differential Degree
LOGFC_CUTOFF <- 1.0      # Threshold for Stable Expression

# --- Functions ---

# 1. Functional Annotation Lookups
fetch_annotations <- function(gene_ids, taxid) {
    # Returns vector of descriptions for given gene IDs
    ah <- AnnotationHub()
    # Logic to find OrgDb would depend on 'taxid'
    # Simplified here for clarity
    return(rep("Protein Description", length(gene_ids))) 
}

# 2. Filtering Logic
extract_rewired_proteins <- function(diff_interactions, de_genes) {
    # diff_interactions: Output from Step 02 (Filtered for Significant Q-values)
    # de_genes: List of differentially expressed genes with LogFC
    
    # A. Calculate Differential Degree
    # Count appearances of each protein in the significant interactions list
    nodes <- c(diff_interactions$ProteinA, diff_interactions$ProteinB)
    diff_degree <- as.data.frame(table(nodes))
    colnames(diff_degree) <- c("GeneID", "DiffDegree")
    
    # B. Merge with Expression Data
    merged <- merge(diff_degree, de_genes, by="GeneID", all.x=TRUE)
    merged$LogFC[is.na(merged$LogFC)] <- 0 # Assume 0 if not in DE list
    
    # C. Apply Filters
    # Threshold: Top 10% of Differential Degree
    degree_thresh <- quantile(merged$DiffDegree[merged$DiffDegree > 0], DEG_PERCENTILE, na.rm=TRUE)
    
    # "Rewired" = High Degree AND Low LogFC
    rewired_candidates <- merged %>% filter(
        DiffDegree >= degree_thresh & 
        abs(LogFC) < LOGFC_CUTOFF
    )
    
    return(rewired_candidates)
}

# --- Main Execution Template ---

run_rewiring_analysis <- function(interaction_file, de_file, output_file) {
    # Load Data
    ints <- read.xlsx(interaction_file) %>% filter(Classification %in% c("Activated", "Repressed"))
    de <- read.xlsx(de_file)
    
    # Analyze
    rewired <- extract_rewired_proteins(ints, de)
    
    # Annotate
    # rewired$Function <- fetch_annotations(rewired$GeneID)
    
    # Save
    write.csv(rewired, output_file, row.names=FALSE)
}
