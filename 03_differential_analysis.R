# 02_differential_analysis.R
# Core Algorithm: Differential Interactome Analysis based on Expression States

library(dplyr)

# --- Parameters ---
THETA_P <- 0.10      # Lower bound probability (Repressed)
THETA_DELTA <- 0.20  # Minimum frequency difference
Q_CRIT_1 <- THETA_P
Q_CRIT_2 <- 1 - THETA_P

# --- Core Algorithm ---

# 1. Binarization
# Converts continuous expression counts to discrete states {-1, 0, 1}
# State 1: Expression > Mean + SD
# State -1: Expression < Mean - SD
binarize_expression <- function(exp_matrix) {
    avg <- rowMeans(exp_matrix)
    sd_val <- apply(exp_matrix, 1, sd)
    
    upper <- avg + sd_val
    lower <- avg - sd_val
    
    bin_matrix <- matrix(0, nrow=nrow(exp_matrix), ncol=ncol(exp_matrix))
    bin_matrix[exp_matrix >= upper] <- 1
    bin_matrix[exp_matrix <= lower] <- -1
    
    return(bin_matrix)
}

# 2. Q-Value Calculation
# Calculates the conditional probability of interaction occurrence in the treatment group
process_interactions <- function(bin_matrix, n_ctl, ppi_reference) {
    # Split Control and Case samples
    mat_ctl <- bin_matrix[, 1:n_ctl]
    mat_case <- bin_matrix[, (n_ctl + 1):ncol(bin_matrix)]
    n_case <- ncol(mat_case)
    
    # Calculate Co-occurrence Frequencies for State [1 1] (Both proteins high)
    # NN: Frequency in Control
    # NT: Frequency in Case (Treatment)
    
    calc_state_freq <- function(mat, idxA, idxB) {
        # Check if both Protein A and Protein B are in State 1 (High)
        rowSums((mat[idxA, ] == 1) & (mat[idxB, ] == 1))
    }
    
    NN11 <- calc_state_freq(mat_ctl, ppi_reference$idxA, ppi_reference$idxB) / n_ctl
    ND11 <- calc_state_freq(mat_case, ppi_reference$idxA, ppi_reference$idxB) / n_case
    
    # Calculate Q-Value (Skew Score)
    # Q = Freq_Case / (Freq_Control + Freq_Case)
    denom <- NN11 + ND11
    Q11 <- ifelse(denom == 0, 0.5, ND11 / denom)
    
    results <- data.frame(
        ProteinA = ppi_reference$ProteinA,
        ProteinB = ppi_reference$ProteinB,
        Freq_Control = NN11,
        Freq_Case = ND11,
        Q_Value = Q11
    )
    
    return(results)
}

# --- Wrapper for Processing a Dataset ---

run_differential_pipeline <- function(expression_data, ppi_db) {
    # expression_data: Matrix of counts (Controls first, then Cases)
    # ppi_db: Reference list of interactors
    
    # Process
    bin_mat <- binarize_expression(expression_data)
    results <- process_interactions(bin_mat, n_ctl=ncol(expression_data)/2, ppi_db) # Assuming balanced design for example
    
    # Classify
    results <- results %>% mutate(
        Classification = case_when(
            Q_Value >= Q_CRIT_2 & pmax(Freq_Control, Freq_Case) >= THETA_DELTA ~ "Activated",
            Q_Value <= Q_CRIT_1 & pmax(Freq_Control, Freq_Case) >= THETA_DELTA ~ "Repressed",
            TRUE ~ "Unchanged"
        )
    )
    
    return(results)
}
