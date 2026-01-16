# 01_preprocess_ids.R
# Standardization of Gene IDs for downstream analysis.

# --- Configuration ---
# Path to directory containing raw count matrices
COUNTS_DIR <- "raw_counts"
# Path to directory containing STRING alias files
MAPPING_DIR <- "string_data"
# Output directory
OUT_DIR <- "processed_data"

if(!dir.exists(OUT_DIR)) dir.create(OUT_DIR)

# --- Functions ---

standardize_ids <- function(count_file, alias_map) {
    # Load raw counts
    counts <- read.csv(count_file, row.names=1)
    
    # Map IDs (Row names) to standardized STRING/Gene IDs
    # Assumes alias_map matches user's gene format to target format
    mapped_ids <- alias_map[rownames(counts)]
    
    # Filter unmapped
    valid_rows <- !is.na(mapped_ids)
    clean_counts <- counts[valid_rows, ]
    rownames(clean_counts) <- mapped_ids[valid_rows]
    
    return(clean_counts)
}

# --- Main Execution ---

# Example usage for a list of datasets
# In a real run, define specific files here
input_files <- list.files(COUNTS_DIR, pattern="*.csv", full.names=TRUE)

for(f in input_files) {
    # Logic to select appropriate alias file would go here
    # processed_counts <- standardize_ids(f, loaded_alias_map)
    # write.csv(processed_counts, file.path(OUT_DIR, basename(f)))
}
