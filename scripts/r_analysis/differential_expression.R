#!/usr/bin/env Rscript

# ==============================================================================
# Differential Expression Analysis using Limma and Wilcoxon Tests
# ==============================================================================
# Description: Performs differential expression analysis comparing early and late
#              stage tumors using both Limma (linear models) and Wilcoxon tests.
#              Generates combined results and visualizations.
# Author: Automated Refactoring
# Date: 2025-12-17
# ==============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library("limma")
  library("edgeR")
  library("dplyr")
  library("effsize")
  library("optparse")
})

# ==============================================================================
# Logger Functions
# ==============================================================================

log_message <- function(message, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] %s: %s\n", timestamp, level, message))
}

log_info <- function(message) {
  log_message(message, "INFO")
}

log_warning <- function(message) {
  log_message(message, "WARNING")
}

log_error <- function(message) {
  log_message(message, "ERROR")
}

log_success <- function(message) {
  log_message(message, "SUCCESS")
}

# ==============================================================================
# Help Documentation Function
# ==============================================================================

#' Display Help for Differential Expression Analysis Pipeline
#'
#' @description
#' Prints comprehensive documentation for the differential expression analysis
#' pipeline including parameters, methods, and usage examples.
#'
#' @return NULL (prints to console)
#' @export
#'
#' @examples
#' show_de_help()
#'
show_de_help <- function() {
  cat("
================================================================================
          Differential Expression Analysis - Limma and Wilcoxon
================================================================================

DESCRIPTION:
  Performs differential expression analysis comparing two experimental groups
  (early and late stage tumors) using two complementary statistical methods:

  1. Limma (Linear Models for Microarray Data):
     - Uses empirical Bayes moderated t-statistics
     - Excellent for small sample sizes
     - Incorporates prior distribution knowledge
     - Provides logFC, t-statistics, and adjusted p-values

  2. Wilcoxon Rank-Sum Test (Mann-Whitney U):
     - Non-parametric alternative to t-tests
     - No assumptions about distribution shape
     - Includes effect size (Cohen's d) estimation with effsize package

  The pipeline:
    1. Loads RNA-Seq expression and clinical data
    2. Maps clinical stages to early/late groups
    3. Filters samples by group membership
    4. Runs both statistical methods on all genes
    5. Adjusts p-values for multiple testing (BH correction)
    6. Identifies significant genes per method and in both methods
    7. Combines results for comprehensive analysis
    8. Saves results to multiple formats

--------------------------------------------------------------------------------
USAGE:
--------------------------------------------------------------------------------

  Rscript scripts/r_analysis/differential_expression.R \\
    --data-path <path_to_rnaseq_data> \\
    --clinical-path <path_to_clinical_data> \\
    --output-dir <output_directory> \\
    --alpha <significance_threshold> \\
    --group1 <first_group_name> \\
    --group2 <second_group_name> \\
    --description <analysis_description>

  # Example with default parameters:
  Rscript scripts/r_analysis/differential_expression.R

  # Example with custom parameters:
  Rscript scripts/r_analysis/differential_expression.R \\
    --data-path data/processed/my_data.csv \\
    --clinical-path data/processed/my_clinical.csv \\
    --alpha 0.05 \\
    --description my_analysis

--------------------------------------------------------------------------------
PARAMETERS:
--------------------------------------------------------------------------------

  --data-path              Path to RNA-Seq expression data (CSV format)
                           Expected format: genes as rows, samples as columns
                           Values should be normalized log2(counts)
                           Default: data/interim/20240930_preprocessed_KIRC/rnaseq_maha.csv

  --clinical-path          Path to clinical metadata (CSV format)
                           Must have column 'ajcc_pathologic_tumor_stage'
                           with values: 'Stage I', 'Stage II', 'Stage III', 'Stage IV'
                           or pre-mapped 'early' and 'late' values
                           Default: data/interim/20240930_preprocessed_KIRC/CuratedClinicalData.csv

  --output-dir             Base directory for output results
                           Creates subdirectories: <output-dir>_limma/ and <output-dir>_wilcox/
                           Default: data/processed/<date>_de_analysis

  --alpha                  Significance threshold for adjusted p-values
                           Range: 0.0 to 1.0
                           Common values: 0.01 (strict), 0.05 (standard)
                           Default: 0.01

  --group1                 Name of first group for contrast
                           Default: early

  --group2                 Name of second group for contrast
                           Default: late

  --description            Description suffix for output directory
                           Used in directory naming for organization
                           Default: preprocessed

  --equivalences           Path to gene equivalence mapping file (optional)
                           Format: tab-separated with Ensembl ID and gene name
                           Default: NULL (no remapping)

  --transpose              Logical flag - transpose data if needed
                           Set to TRUE if samples are rows, genes are columns
                           Default: auto-detect based on sample matching

--------------------------------------------------------------------------------
METHODS DESCRIPTION:
--------------------------------------------------------------------------------

  Limma Analysis:
    - Model: Linear model fit with group contrast
    - Shrinkage: Empirical Bayes moderation of t-statistics
    - P-value correction: Benjamini-Hochberg (FDR)
    - Outputs: logFC, AveExpr, t, P.Value, adj.P.Val, B

  Wilcoxon Test:
    - Test type: Two-sided rank-sum test
    - Effect size: Cohen's d with 95% confidence interval
    - P-value correction: Benjamini-Hochberg (FDR)
    - Outputs: W statistic, p-value, effect size, Cohen's d CI

  Consensus Results:
    - Genes significant in BOTH methods are considered most robust
    - All genes from either method are included in combined results
    - Enables cross-validation of findings

--------------------------------------------------------------------------------
OUTPUT FILES:
--------------------------------------------------------------------------------

  In <output-dir>_limma/:
    - sorted_table.csv                      All genes with statistics (Limma)
    - significant_genes.csv                 Limma-significant genes only
    - combined_results_all_genes.csv        All genes from both methods
    - combined_results_both_significant.csv Genes significant in both methods
    - summary_statistics.csv                Summary counts

  In <output-dir>_wilcox/:
    - wilcox_results.csv                    All genes with statistics (Wilcox)
    - significant_genes_wilcox.csv          Wilcoxon-significant genes only
    - combined_results_all_genes.csv        All genes from both methods
    - combined_results_both_significant.csv Genes significant in both methods
    - summary_statistics.csv                Summary counts

  Column Descriptions:
    logFC              Log2 fold change (Limma only)
    AveExpr            Average expression across samples (Limma only)
    t                  t-statistic (Limma only)
    P.Value            Unadjusted p-value
    adj.P.Val          Adjusted p-value (Benjamini-Hochberg FDR)
    B                  Log-odds of differential expression (Limma only)
    W                  Wilcoxon test statistic
    effect_size        Cohen's d estimate
    cohens_d_lower     Lower bound of 95% CI
    cohens_d_upper     Upper bound of 95% CI

--------------------------------------------------------------------------------
EXAMPLES:
--------------------------------------------------------------------------------

  # Example 1: Basic usage with default parameters
  Rscript scripts/r_analysis/differential_expression.R

  # Example 2: Custom data paths
  Rscript scripts/r_analysis/differential_expression.R \\
    --data-path data/my_rnaseq.csv \\
    --clinical-path data/my_clinical.csv

  # Example 3: Relaxed significance threshold for exploratory analysis
  Rscript scripts/r_analysis/differential_expression.R \\
    --alpha 0.05 \\
    --description exploratory

  # Example 4: Custom group names with stricter threshold
  Rscript scripts/r_analysis/differential_expression.R \\
    --group1 early_stage \\
    --group2 late_stage \\
    --alpha 0.001

  # Example 5: Save to custom directory
  Rscript scripts/r_analysis/differential_expression.R \\
    --output-dir results/my_analysis \\
    --description batch_corrected

--------------------------------------------------------------------------------
ERROR HANDLING:
--------------------------------------------------------------------------------

  - Missing input files: Error is logged and execution stops
  - Data dimension mismatch: Auto-transposes if needed; warns if forced
  - Invalid clinical stages: Error if expected stage values not found
  - Missing genes: Continuing with available genes; warnings logged
  - Failed test computation: Error logged per gene; analysis continues
  - Missing output directories: Automatically created with recursive=TRUE

--------------------------------------------------------------------------------
LOGGING:
--------------------------------------------------------------------------------

  All major steps are logged with timestamps:
    - Input file loading
    - Sample/gene counts
    - Clinical stage mapping
    - Test execution progress
    - Significant gene counts (per method)
    - Output file saves
    - Summary statistics

  Log levels: INFO, WARNING, ERROR, SUCCESS

--------------------------------------------------------------------------------
CITATION:
--------------------------------------------------------------------------------

  If you use Limma in published research, please cite:

  Ritchie, M. E., Phipson, B., Wu, D., et al. (2015).
  limma powers differential expression analyses for RNA-sequencing and
  microarray studies. Nucleic Acids Research, 43(7), e47.

  For Wilcoxon test methodology:

  Wilcoxon, F. (1945). Individual comparisons by ranking methods.
  Biometrics, 1(6), 80-83.

  For effect size (Cohen's d):

  Cohen, J. (1988). Statistical Power Analysis for the Behavioral Sciences.
  Lawrence Erlbaum Associates.

================================================================================
\n")
  invisible(NULL)
}

# ==============================================================================
# Command Line Arguments
# ==============================================================================

option_list <- list(
  make_option(
    c("--data-path"),
    type = "character",
    default = "data/interim/20240930_preprocessed_KIRC/rnaseq_maha.csv",
    help = "Path to RNA-Seq expression data (CSV format) [default: %default]",
    metavar = "character"
  ),
  make_option(
    c("--clinical-path"),
    type = "character",
    default = "data/interim/20240930_preprocessed_KIRC/CuratedClinicalData.csv",
    help = "Path to clinical metadata (CSV format) [default: %default]",
    metavar = "character"
  ),
  make_option(
    c("--output-dir"),
    type = "character",
    default = NULL,
    help = "Base output directory [default: data/processed/<date>_de_analysis]",
    metavar = "character"
  ),
  make_option(
    c("--alpha"),
    type = "numeric",
    default = 0.01,
    help = "Significance threshold for adjusted p-values [default: %default]",
    metavar = "numeric"
  ),
  make_option(
    c("--group1"),
    type = "character",
    default = "early",
    help = "Name of first group for contrast [default: %default]",
    metavar = "character"
  ),
  make_option(
    c("--group2"),
    type = "character",
    default = "late",
    help = "Name of second group for contrast [default: %default]",
    metavar = "character"
  ),
  make_option(
    c("--description"),
    type = "character",
    default = "preprocessed",
    help = "Description suffix for output directory [default: %default]",
    metavar = "character"
  ),
  make_option(
    c("--equivalences"),
    type = "character",
    default = NULL,
    help = "Path to gene equivalence mapping file (optional) [default: %default]",
    metavar = "character"
  ),
  make_option(
    c("--transpose"),
    type = "logical",
    default = FALSE,
    help = "Transpose data if needed (samples as rows) [default: %default]",
    metavar = "logical"
  )
)

opt_parser <- OptionParser(
  option_list = option_list,
  description = "\nDifferential Expression Analysis - Limma and Wilcoxon"
)
opt <- parse_args(opt_parser)

# ==============================================================================
# Initialize
# ==============================================================================

start_time <- Sys.time()

# Get today's date without dashes
today <- format(Sys.Date(), "%Y%m%d")

# Determine output paths
if (is.null(opt$`output-dir`)) {
  base_output <- paste0('data/processed/', today, '_degs_analysis')
} else {
  base_output <- opt$`output-dir`
}

save_path <- paste0(base_output, '_limma/', opt$description)
save_path_wc <- paste0(base_output, '_wilcox/', opt$description)

group_to_analyze <- c(opt$group1, opt$group2)
group_valid_names <- sort(group_to_analyze)

alpha <- opt$alpha
path_data <- opt$`data-path`
path_clinical <- opt$`clinical-path`
path_equivalences <- opt$equivalences
description <- opt$description


# ==============================================================================
# Load Data
# ==============================================================================

log_info("=== Differential Expression Analysis Pipeline ===")
log_info(sprintf("Started at: %s", Sys.time()))

# Load RNA-Seq data
log_info(sprintf("Reading RNA-Seq data from: %s", path_data))
data <- tryCatch(
  {
    df <- read.csv(path_data, header = TRUE, row.names = 1)
    log_success(sprintf("Loaded RNA-Seq data: %d genes x %d samples",
                       nrow(df), ncol(df)))
    df
  },
  error = function(e) {
    log_error(sprintf("Failed to read RNA-Seq data: %s", e$message))
    stop(e)
  }
)

# Standardize sample names (change "." to "-")
colnames(data) <- gsub("\\.", "-", colnames(data))

# Load clinical data
log_info(sprintf("Reading clinical data from: %s", path_clinical))
clinical_df <- tryCatch(
  {
    df <- read.csv(path_clinical, row.names = 1)
    log_success(sprintf("Loaded clinical data: %d samples", nrow(df)))
    df
  },
  error = function(e) {
    log_error(sprintf("Failed to read clinical data: %s", e$message))
    stop(e)
  }
)

# ==============================================================================
# Process Clinical Data and Map Stages
# ==============================================================================

log_info("Processing clinical stage information...")

# Define mapping from TCGA stages to early/late groups
dict_early_late <- c(
  "Stage I" = "early",
  "Stage II" = "early",
  "Stage III" = "late",
  "Stage IV" = "late"
)

# Check if clinical data has expected stage values and map accordingly
if (any(c("Stage I", "Stage II", "Stage III", "Stage IV") %in%
        clinical_df$ajcc_pathologic_tumor_stage)) {
  log_info("Mapping TCGA stages (Stage I-IV) to early/late groups")
  clinical <- setNames(
    dict_early_late[clinical_df$ajcc_pathologic_tumor_stage],
    rownames(clinical_df)
  )
} else if (all(c("early", "late") %in% clinical_df$ajcc_pathologic_tumor_stage)) {
  log_info("Clinical data already has early/late stage mapping")
  clinical <- setNames(
    clinical_df$ajcc_pathologic_tumor_stage,
    rownames(clinical_df)
  )
} else {
  log_error("Clinical data does not have expected stage values (Stage I-IV or early/late)")
  stop("Clinical data does not have expected stage values.")
}

# Convert to data frame format for downstream analysis
clinical_df <- data.frame(
  Sample = names(clinical),
  Subgroup = as.vector(clinical),
  row.names = NULL,
  stringsAsFactors = FALSE
)

# Log stage distribution
stage_counts <- table(clinical_df$Subgroup)
log_info(sprintf("Stage distribution: %s",
                 paste(names(stage_counts), "=", stage_counts, collapse = ", ")))


# ==============================================================================
# Data Dimension Checking and Sample Selection
# ==============================================================================

# Check if rows match between data and clinical
log_info("Checking data dimensions...")
log_info(sprintf("Data dimensions: %d genes x %d samples", nrow(data), ncol(data)))
log_info(sprintf("Clinical samples: %d", nrow(clinical_df)))

if (nrow(data) != nrow(clinical_df)) {
  if (ncol(data) == nrow(clinical_df)) {
    log_warning("Data needs transposition (samples are columns in data)")
    log_info("Transposing data so samples are rows and genes are columns...")
    data <- t(data)
  } else {
    log_error(sprintf(
      "Data/clinical dimension mismatch: data has %d rows, clinical has %d samples",
      nrow(data), nrow(clinical_df)
    ))
    stop("Cannot match data dimensions with clinical samples")
  }
}

stopifnot(nrow(data) == nrow(clinical_df))
log_success(sprintf("Data dimensions verified: %d samples x %d genes",
                   nrow(data), ncol(data)))

# Subset data to match clinical samples
data_select <- data[rownames(data) %in% clinical_df$Sample, ]
log_info(sprintf("Selected samples: %d", nrow(data_select)))

# ==============================================================================
# Create Sample Groups for Analysis
# ==============================================================================

log_info(sprintf("Filtering for groups: %s", paste(group_to_analyze, collapse = ", ")))

patients_early <- clinical_df$Sample[clinical_df$Subgroup == opt$group1]
patients_late <- clinical_df$Sample[clinical_df$Subgroup == opt$group2]

log_info(sprintf("Group '%s': %d samples", opt$group1, length(patients_early)))
log_info(sprintf("Group '%s': %d samples", opt$group2, length(patients_late)))

# Extract group-specific expression matrices
early_df <- data_select[rownames(data_select) %in% patients_early, , drop = FALSE]
late_df <- data_select[rownames(data_select) %in% patients_late, , drop = FALSE]

log_success("Sample grouping completed")

# ==============================================================================
# Wilcoxon Rank-Sum Test Functions
# ==============================================================================

#' Perform Wilcoxon Rank-Sum Test on Gene Expression
#'
#' @description
#' Performs a two-sided Wilcoxon rank-sum test (Mann-Whitney U test) comparing
#' gene expression between two groups. This non-parametric test is robust to
#' outliers and does not assume normally distributed data.
#'
#' @param gene_a Numeric vector of expression values from first group
#' @param gene_b Numeric vector of expression values from second group
#' @param alpha Numeric significance threshold for confidence intervals
#'   (default: 0.01)
#'
#' @return A list with components:
#'   - statistic: The Wilcoxon test statistic (W value)
#'   - p.value: Two-sided p-value from the test
#'   - alternative: Character string "two.sided" (fixed)
#'   - method: Character string describing the test
#'
#' @details
#' Uses exact test when possible with continuity correction for tied ranks.
#' Paired=FALSE for independent samples (Mann-Whitney U test).
#'
#' @examples
#' group1 <- rnorm(20, mean = 5, sd = 1)
#' group2 <- rnorm(20, mean = 6, sd = 1)
#' result <- wilcox_function(group1, group2, alpha = 0.05)
#'
#' @keywords internal
wilcox_function <- function(gene_a, gene_b, alpha) {
  # Setting paired to FALSE to perform Mann-Whitney U test
  results <- wilcox.test(
    x = gene_a,
    y = gene_b,
    alternative = "two.sided",
    paired = FALSE,
    exact = NULL,
    correct = TRUE,
    conf.int = FALSE,
    conf.level = 1 - alpha
  )

  return(results)
}

#' Perform Wilcoxon Test on All Genes Between Two Groups
#'
#' @description
#' Applies the Wilcoxon rank-sum test to all genes comparing two experimental
#' groups. Includes calculation of effect sizes (Cohen's d) with confidence
#' intervals and multiple testing correction.
#'
#' @param data_early Numeric matrix or data frame with expression values for
#'   first group (samples x genes)
#' @param data_late Numeric matrix or data frame with expression values for
#'   second group (samples x genes)
#' @param alpha Numeric significance threshold (default: 0.01)
#'
#' @return A data frame with columns:
#'   - gene: Gene names (row names of result)
#'   - W: Wilcoxon test statistic
#'   - p.value: Unadjusted p-value
#'   - alternative: Test direction ("two.sided")
#'   - effect_size: Cohen's d estimate
#'   - cohens_d_lower: Lower 95% CI bound
#'   - cohens_d_upper: Upper 95% CI bound
#'   - adj.p.value: BH-adjusted p-value
#'
#' @details
#' Performs robust non-parametric testing with effect size estimation.
#' Automatically handles and logs errors for individual genes.
#' Multiple testing correction using Benjamini-Hochberg FDR method.
#'
#' @seealso
#' \code{\link{wilcox_function}} for individual gene tests
#'
#' @keywords internal
wilcox_all_genes <- function(data_early, data_late, alpha = 0.01) {
  # Get common genes between both dataframes
  genes <- colnames(data_late)

  # Initialize results dataframe
  results_df <- data.frame(
    gene = character(),
    W = numeric(),
    p.value = numeric(),
    statistic = numeric(),
    alternative = character(),
    method = character(),
    stringsAsFactors = FALSE
  )

  log_info(sprintf("Running Wilcoxon tests on %d genes...", length(genes)))

  # Loop through each gene
  for (i in seq_along(genes)) {
    gene <- genes[i]
    tryCatch(
      {
        # Perform wilcox test for current gene
        test_result <- wilcox_function(
          data_early[, gene],
          data_late[, gene],
          alpha
        )
        # Calculate Cohen's d effect size
        cohens_result <- cohen.d(
          data_early[, gene],
          data_late[, gene],
          conf.level = 1 - alpha
        )
        # Add results to dataframe
        results_df <- rbind(results_df, data.frame(
          gene = gene,
          W = test_result$statistic,
          p.value = test_result$p.value,
          alternative = test_result$alternative,
          effect_size = cohens_result$estimate,
          cohens_d_lower = cohens_result$conf.int[1],
          cohens_d_upper = cohens_result$conf.int[2],
          stringsAsFactors = FALSE
        ))
      },
      error = function(e) {
        log_warning(sprintf("Error with gene %s: %s", gene, e$message))
      }
    )
  }

  # Set gene names as row names and remove gene column
  rownames(results_df) <- results_df$gene
  results_df$gene <- NULL

  log_success(sprintf("Wilcoxon tests completed for %d genes", nrow(results_df)))
  return(results_df)
}

# ==============================================================================
# Execute Wilcoxon Tests
# ==============================================================================

log_info("========================================================================")
log_info("STEP 1: Running Wilcoxon Rank-Sum Tests")
log_info("========================================================================")

wilcox_results <- wilcox_all_genes(early_df, late_df, alpha)

# Add adjusted p-values
wilcox_results$adj.p.value <- p.adjust(wilcox_results$p.value, method = "BH")

# Filter for significant genes
significant_genes_wc <- wilcox_results[wilcox_results$adj.p.value < alpha, ]
log_success(sprintf("Wilcoxon: Found %d significant genes", nrow(significant_genes_wc)))

# Save the results
log_info("Saving Wilcoxon test results...")
if (!dir.exists(save_path_wc)) {
  dir.create(save_path_wc, recursive = TRUE)
  log_info(sprintf("Created output directory: %s", save_path_wc))
}

tryCatch(
  {
    write.table(
      wilcox_results,
      file = file.path(save_path_wc, 'wilcox_results.csv'),
      row.names = TRUE,
      sep = ",",
      quote = FALSE
    )
    log_success("Saved: wilcox_results.csv")

    write.table(
      significant_genes_wc,
      file = file.path(save_path_wc, 'significant_genes_wilcox.csv'),
      row.names = TRUE,
      sep = ",",
      quote = FALSE
    )
    log_success("Saved: significant_genes_wilcox.csv")
  },
  error = function(e) {
    log_error(sprintf("Failed to save Wilcoxon results: %s", e$message))
  }
)

# ==============================================================================
# Limma Differential Expression Analysis
# ==============================================================================

#' Fit Linear Models for Differential Expression (Limma)
#'
#' @description
#' Uses the Limma package to fit linear models to gene expression data and
#' perform empirical Bayes moderated t-tests. This approach is particularly
#' powerful for small sample sizes due to shrinkage of standard errors.
#'
#' @details
#' The analysis proceeds in three steps:
#' 1. Fit linear model with design matrix containing group assignments
#' 2. Apply contrast matrix to extract comparisons between groups
#' 3. Use empirical Bayes moderation to improve estimates of standard error
#'
#' Outputs include log-fold change (logFC), moderated t-statistics, and
#' adjusted p-values using the Benjamini-Hochberg FDR control method.
#'
#' @seealso
#' \code{\link[limma]{lmFit}} for model fitting
#' \code{\link[limma]{makeContrasts}} for contrast specification
#' \code{\link[limma]{eBayes}} for empirical Bayes moderation
#'
#' @keywords internal
#'
#' @return No explicit return; saves results to output files

log_info("========================================================================")
log_info("STEP 2: Running Limma Differential Expression Analysis")
log_info("========================================================================")

log_info("Creating design matrix...")

# Define model matrix with group factors
mm <- model.matrix(~0 + clinical_df$Subgroup)
# Adjust the column names
colnames(mm) <- gsub("clinical_df$Subgroup", "", colnames(mm), fixed = TRUE)

log_info(sprintf("Design matrix created with groups: %s",
                 paste(colnames(mm), collapse = ", ")))

# Fit linear models
log_info("Fitting linear models to expression data...")
fit <- lmFit(t(data_select), mm)

log_info(sprintf("Model fit completed. Coefficients shape: %d genes x %d groups",
                 nrow(fit$coefficients), ncol(fit$coefficients)))

# Define contrasts between groups
contr <- makeContrasts(
  comparisons = paste(opt$group1, opt$group2, sep = "-"),
  levels = group_valid_names
)

log_info(sprintf("Applying contrast: %s - %s", opt$group1, opt$group2))

# Estimate contrasts for each gene
tmp <- contrasts.fit(fit, contr)
log_success("Contrasts fitted for all genes")

# Bayes smoothing for standard errors
log_info("Applying empirical Bayes moderation of standard errors...")
tmp <- eBayes(tmp)
log_success("Empirical Bayes moderation completed")

# Get the sorted differentially expressed genes
log_info("Extracting differentially expressed genes...")
sorted.table <- topTable(tmp, sort.by = "logFC", n = Inf)

# Filter for significant genes
significant_genes_limma <- sorted.table[sorted.table$adj.P.Val < alpha, ]
log_success(sprintf("Limma: Found %d significant genes", nrow(significant_genes_limma)))

# Save the results
log_info("Saving Limma results...")
if (!dir.exists(save_path)) {
  dir.create(save_path, recursive = TRUE)
  log_info(sprintf("Created output directory: %s", save_path))
}

tryCatch(
  {
    write.table(
      sorted.table,
      file = file.path(save_path, 'sorted_table.csv'),
      row.names = TRUE,
      sep = ",",
      quote = FALSE
    )
    log_success("Saved: sorted_table.csv")

    write.table(
      significant_genes_limma,
      file = file.path(save_path, 'significant_genes.csv'),
      row.names = TRUE,
      sep = ",",
      quote = FALSE
    )
    log_success("Saved: significant_genes.csv")
  },
  error = function(e) {
    log_error(sprintf("Failed to save Limma results: %s", e$message))
  }
)

# ==============================================================================
# Combine Results from Both Methods
# ==============================================================================

log_info("========================================================================")
log_info("STEP 3: Combining Limma and Wilcoxon Results")
log_info("========================================================================")

log_info("Preparing results for combination...")

# Prepare limma results with renamed columns to avoid conflicts
limma_combined <- sorted.table
colnames(limma_combined) <- paste0("limma_", colnames(limma_combined))

# Prepare wilcox results with renamed columns
wilcox_combined <- wilcox_results
colnames(wilcox_combined) <- paste0("wilcox_", colnames(wilcox_combined))

# Combine results by gene name (row names)
# Create a full outer join
all_genes <- sort(unique(c(rownames(limma_combined), rownames(wilcox_combined))))
combined_results <- data.frame(row.names = all_genes)

log_info(sprintf("Merging results for %d total genes", length(all_genes)))

# Add limma results
for (col in colnames(limma_combined)) {
  combined_results[[col]] <- NA_real_
  matching_genes <- intersect(rownames(combined_results), rownames(limma_combined))
  combined_results[matching_genes, col] <- limma_combined[matching_genes, col]
}

# Add wilcox results
for (col in colnames(wilcox_combined)) {
  combined_results[[col]] <- NA_real_
  matching_genes <- intersect(rownames(combined_results), rownames(wilcox_combined))
  combined_results[matching_genes, col] <- wilcox_combined[matching_genes, col]
}

log_success(sprintf("Combined results created for %d genes", nrow(combined_results)))

# Save combined results (all genes)
log_info("Saving combined results (all genes)...")
tryCatch(
  {
    write.table(
      combined_results,
      file = file.path(save_path_wc, 'combined_results_all_genes.csv'),
      row.names = TRUE,
      sep = ",",
      quote = FALSE
    )
    log_success("Saved: combined_results_all_genes.csv (wilcox dir)")

    # Also save in limma folder for convenience
    write.table(
      combined_results,
      file = file.path(save_path, 'combined_results_all_genes.csv'),
      row.names = TRUE,
      sep = ",",
      quote = FALSE
    )
    log_success("Saved: combined_results_all_genes.csv (limma dir)")
  },
  error = function(e) {
    log_error(sprintf("Failed to save combined results: %s", e$message))
  }
)

log_info(sprintf("Combined results (all genes) contain %d genes", nrow(combined_results)))

# ==============================================================================
# Identify Consensus Genes (Significant in Both Methods)
# ==============================================================================

log_info("Identifying genes significant in BOTH methods...")

limma_sig <- rownames(sorted.table)[sorted.table$adj.P.Val < alpha]
wilcox_sig <- rownames(wilcox_results)[wilcox_results$adj.p.value < alpha]
both_significant <- intersect(limma_sig, wilcox_sig)

log_info(sprintf("Significant in Limma only: %d genes",
                 length(setdiff(limma_sig, wilcox_sig))))
log_info(sprintf("Significant in Wilcoxon only: %d genes",
                 length(setdiff(wilcox_sig, limma_sig))))
log_success(sprintf("Significant in BOTH methods: %d genes", length(both_significant)))

# Create combined results for genes significant in both methods
combined_significant <- combined_results[both_significant, ]

# Save combined results for genes significant in both methods
log_info("Saving consensus results (genes significant in both methods)...")
tryCatch(
  {
    write.table(
      combined_significant,
      file = file.path(save_path_wc, 'combined_results_both_significant.csv'),
      row.names = TRUE,
      sep = ",",
      quote = FALSE
    )
    log_success("Saved: combined_results_both_significant.csv (wilcox dir)")

    # Also save in limma folder for convenience
    write.table(
      combined_significant,
      file = file.path(save_path, 'combined_results_both_significant.csv'),
      row.names = TRUE,
      sep = ",",
      quote = FALSE
    )
    log_success("Saved: combined_results_both_significant.csv (limma dir)")
  },
  error = function(e) {
    log_error(sprintf("Failed to save consensus results: %s", e$message))
  }
)

# ==============================================================================
# Generate Summary Statistics
# ==============================================================================

log_info("========================================================================")
log_info("STEP 4: Summary Statistics")
log_info("========================================================================")

summary_stats <- data.frame(
  Method = c("Limma", "Wilcoxon", "Both methods"),
  Significant_genes = c(length(limma_sig), length(wilcox_sig), length(both_significant)),
  stringsAsFactors = FALSE
)

log_info("Summary of significant genes per method:")
for (i in seq_len(nrow(summary_stats))) {
  log_info(sprintf("  - %s: %d genes",
                   summary_stats$Method[i],
                   summary_stats$Significant_genes[i]))
}

# Save summary statistics
log_info("Saving summary statistics...")
tryCatch(
  {
    write.table(
      summary_stats,
      file = file.path(save_path_wc, 'summary_statistics.csv'),
      row.names = FALSE,
      sep = ",",
      quote = FALSE
    )
    log_success("Saved: summary_statistics.csv (wilcox dir)")

    write.table(
      summary_stats,
      file = file.path(save_path, 'summary_statistics.csv'),
      row.names = FALSE,
      sep = ",",
      quote = FALSE
    )
    log_success("Saved: summary_statistics.csv (limma dir)")
  },
  error = function(e) {
    log_error(sprintf("Failed to save summary statistics: %s", e$message))
  }
)

# ==============================================================================
# Execution Summary and Finalization
# ==============================================================================

end_time <- Sys.time()
elapsed <- end_time - start_time

log_info("========================================================================")
log_info("ANALYSIS COMPLETED")
log_info("========================================================================")
log_success(sprintf("Total execution time: %.2f seconds (%.2f minutes)",
                   as.numeric(difftime(end_time, start_time, units = "secs")),
                   as.numeric(difftime(end_time, start_time, units = "mins"))))

log_info("Output directories:")
log_info(sprintf("  Limma results: %s", save_path))
log_info(sprintf("  Wilcoxon results: %s", save_path_wc))

log_info("Key outputs:")
log_info(sprintf("  - %d total genes analyzed", nrow(combined_results)))
log_info(sprintf("  - %d genes significant (Limma, FDR < %s)", length(limma_sig), alpha))
log_info(sprintf("  - %d genes significant (Wilcoxon, FDR < %s)", length(wilcox_sig), alpha))
log_info(sprintf("  - %d genes significant in BOTH methods (consensus set)", length(both_significant)))

log_success("Differential expression analysis completed successfully!")
log_info(sprintf("Finished at: %s", Sys.time()))
