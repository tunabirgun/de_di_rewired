# 05_visualization.R
# Visualization of Differential and Rewired Results

library(UpSetR)
library(VennDiagram)
library(ggplot2)
library(grid)

# --- Plotting Functions ---

# 1. UpSet Plot
# Visualizes intersections between multiple sets (e.g., rewired proteins across conditions)
generate_upset_plot <- function(data_frame, group_col, id_col, output_file) {
    # Convert long-format data to list
    list_input <- split(data_frame[[id_col]], data_frame[[group_col]])
    
    # Generate Matrix
    bin_mat <- fromList(list_input)
    
    # Plot
    png(output_file, width=1200, height=800, res=100)
    upset(bin_mat, order.by = "freq", nsets = length(list_input))
    dev.off()
}

# 2. Venn Diagram
# Compares defined sets (e.g., Condition A vs Condition B)
generate_venn_diagram <- function(set_list, title, output_file) {
    # set_list: Named list of vectors (e.g., list(A=c(...), B=c(...)))
    
    vp <- venn.diagram(
        x = set_list,
        filename = NULL,
        category.names = names(set_list),
        main = title,
        fill = rainbow(length(set_list)),
        alpha = 0.5
    )
    
    png(output_file, width=800, height=800)
    grid.draw(vp)
    dev.off()
}

# --- Main Execution Template ---

run_visualization_suite <- function(rewired_file, out_dir) {
    data <- read.csv(rewired_file)
    
    # UpSet
    generate_upset_plot(data, "Dataset", "GeneID", file.path(out_dir, "UpSet_Rewired.png"))
    
    # Venn (Example: Comparing two groups defined by user)
    # group1 <- data$GeneID[data$Condition == "A"]
    # group2 <- data$GeneID[data$Condition == "B"]
    # generate_venn_diagram(list(A=group1, B=group2), "Comparison", file.path(out_dir, "Venn.png"))
}
