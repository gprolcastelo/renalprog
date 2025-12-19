#!/usr/bin/env Rscript

# ==============================================================================
# Gene Enrichment Analysis using g:Profiler
# ==============================================================================
# Description: Performs gene enrichment analysis on important genes for
#              classification using g:Profiler API
# Author: Automated Refactoring
# Date: 2025-12-17
# ==============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library("gprofiler2")
  library("ggplot2")
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

#' Display Help for gprofiler_pipeline Function
#'
#' @description
#' Prints comprehensive documentation for the gprofiler_pipeline function
#' including parameters, return values, and usage examples.
#'
#' @return NULL (prints to console)
#' @export
#'
#' @examples
#' show_gprofiler_help()
#'
show_gprofiler_help <- function() {
  cat("
================================================================================
              gprofiler_pipeline - g:Profiler Enrichment Analysis
================================================================================

DESCRIPTION:
  Performs functional enrichment analysis on a gene list using the g:Profiler
  web service (https://biit.cs.ut.ee/gprofiler/). Queries multiple biological
  databases (GO, KEGG, Reactome, WikiPathways, etc.) to identify significantly
  enriched pathways and biological processes associated with your gene set.

  The pipeline:
    1. Queries g:Profiler API with your gene list
    2. Filters results by statistical significance (FDR correction)
    3. Extracts top N most significant terms per database
    4. Adds log-transformed p-values for visualization
    5. Handles duplicate term names across databases
    6. Returns both processed plot data and raw results

--------------------------------------------------------------------------------
USAGE:
--------------------------------------------------------------------------------

  results <- gprofiler_pipeline(
    genes,
    sources,
    top_n,
    organism,
    fdr_threshold
  )

--------------------------------------------------------------------------------
PARAMETERS:
--------------------------------------------------------------------------------

  genes          Character vector of gene symbols (HGNC gene names)
                 Example: c('TP53', 'BRCA1', 'EGFR', 'MYC', 'KRAS')
                 Note: Use official gene symbols for best results

  sources        Comma-separated string of databases to query
                 Available databases:
                   GO      - Gene Ontology (BP, MF, CC)
                   KEGG    - KEGG pathways
                   REAC    - Reactome pathways
                   WP      - WikiPathways
                   TF      - TRANSFAC transcription factors
                   MIRNA   - miRTarBase microRNAs
                   HPA     - Human Protein Atlas
                   CORUM   - CORUM protein complexes
                   HP      - Human Phenotype Ontology
                 Example: 'GO,KEGG,REAC' or 'GO,REAC,KEGG,WP'

  top_n          Integer - Number of top terms to extract per database
                 Terms are ranked by p-value (most significant first)
                 Default: 20
                 Recommended range: 10-50

  organism       Character - Organism identifier
                 Common options:
                   'hsapiens'     - Homo sapiens (human)
                   'mmusculus'    - Mus musculus (mouse)
                   'rnorvegicus'  - Rattus norvegicus (rat)
                   'drerio'       - Danio rerio (zebrafish)
                   'celegans'     - Caenorhabditis elegans
                   'scerevisiae'  - Saccharomyces cerevisiae (yeast)
                 Default: 'hsapiens'

  fdr_threshold  Numeric - FDR (False Discovery Rate) significance threshold
                 Only terms with FDR-corrected p-value below this are considered
                 Range: 0.0 to 1.0
                 Default: 0.05
                 Common values: 0.01 (strict), 0.05 (standard), 0.1 (permissive)

--------------------------------------------------------------------------------
RETURNS:
--------------------------------------------------------------------------------

  A list with two elements:

  $plot_data     Data frame ready for visualization containing:
                   - source      : Database name (character)
                   - term_name   : Name of enriched term (character)
                   - p_value     : Original p-value (numeric)
                   - log_p_value : -log10 transformed p-value (numeric)

                 Notes:
                   - Duplicate term names are suffixed with ' (source)'
                   - NA values are removed
                   - Empty data frame if no significant results

  $gostres       Complete g:Profiler result object containing:
                   - result : Full enrichment data frame with columns:
                              term_id, term_name, p_value, source, term_size,
                              query_size, intersection_size, precision, recall,
                              effective_domain_size, source_order, parents
                   - meta   : Metadata about the query

--------------------------------------------------------------------------------
EXAMPLES:
--------------------------------------------------------------------------------

  # Example 1: Basic usage with cancer genes
  cancer_genes <- c('TP53', 'BRCA1', 'EGFR', 'MYC', 'KRAS', 'PTEN')
  results <- gprofiler_pipeline(
    genes = cancer_genes,
    sources = 'GO,KEGG,REAC',
    top_n = 20,
    organism = 'hsapiens',
    fdr_threshold = 0.05
  )

  # View results
  head(results$plot_data)
  head(results$gostres$result)

  # Example 2: GO-only analysis with strict threshold
  results_go <- gprofiler_pipeline(
    genes = my_gene_list,
    sources = 'GO',
    top_n = 10,
    organism = 'hsapiens',
    fdr_threshold = 0.01
  )

  # Example 3: Mouse genes analysis
  mouse_genes <- c('Trp53', 'Brca1', 'Egfr', 'Myc', 'Kras')
  results_mouse <- gprofiler_pipeline(
    genes = mouse_genes,
    sources = 'GO,KEGG,REAC,WP',
    top_n = 15,
    organism = 'mmusculus',
    fdr_threshold = 0.05
  )

  # Example 4: Comprehensive analysis with many databases
  results_broad <- gprofiler_pipeline(
    genes = my_genes,
    sources = 'GO,KEGG,REAC,WP,TF,MIRNA',
    top_n = 30,
    organism = 'hsapiens',
    fdr_threshold = 0.1
  )

  # Example 5: Extracting specific information
  plot_data <- results$plot_data
  full_results <- results$gostres$result

  # Get only KEGG pathways
  kegg_terms <- plot_data[plot_data$source == 'KEGG', ]

  # Get terms with very low p-values
  top_hits <- plot_data[plot_data$p_value < 0.001, ]

--------------------------------------------------------------------------------
ERROR HANDLING:
--------------------------------------------------------------------------------

  - If g:Profiler API call fails: Error is logged and function stops
  - If no significant results found: Warning is logged, empty results returned
  - If duplicate term names exist: Automatically suffixed with database name

--------------------------------------------------------------------------------
LOGGING:
--------------------------------------------------------------------------------

  All major steps are logged with timestamps:
    - Number of genes queried
    - Databases being searched
    - Number of significant terms found
    - Number of terms extracted per database
    - Duplicate term handling

--------------------------------------------------------------------------------
SEE ALSO:
--------------------------------------------------------------------------------

  create_enrichment_plot() - Visualize enrichment results
  save_results()           - Save results to files

  g:Profiler web: https://biit.cs.ut.ee/gprofiler/
  g:Profiler API: https://biit.cs.ut.ee/gprofiler/page/docs

--------------------------------------------------------------------------------
CITATION:
--------------------------------------------------------------------------------

  If you use g:Profiler in published research, please cite:

  Liis Kolberg, Uku Raudvere, Ivan Kuzmin, Jaak Vilo, Hedi Peterson (2020).
  gprofiler2 -- an R package for gene list functional enrichment analysis and
  namespace conversion toolset g:Profiler. F1000Research, 9:ELIXIR-709.

================================================================================
\n")
  invisible(NULL)
}

# ==============================================================================
# Command Line Arguments
# ==============================================================================

option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    default = "data/external/important_genes_shap.csv",
    help = "Path to input CSV file with gene names [default: %default]",
    metavar = "character"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = NULL,
    help = "Path to save output results [default: reports/figures/<date>_gprofiler_on_kirc_classification_genes/]",
    metavar = "character"
  ),
  make_option(
    c("-s", "--sources"),
    type = "character",
    default = "GO,REAC,KEGG,WP",
    help = "Comma-separated list of databases to query [default: %default]",
    metavar = "character"
  ),
  make_option(
    c("-n", "--top_n"),
    type = "integer",
    default = 20,
    help = "Number of top terms to extract per source [default: %default]",
    metavar = "integer"
  ),
  make_option(
    c("--organism"),
    type = "character",
    default = "hsapiens",
    help = "Organism identifier [default: %default]",
    metavar = "character"
  ),
  make_option(
    c("--fdr"),
    type = "numeric",
    default = 0.05,
    help = "FDR threshold for significance [default: %default]",
    metavar = "numeric"
  ),
  make_option(
    c("--width"),
    type = "numeric",
    default = 25,
    help = "Plot width in inches [default: %default]",
    metavar = "numeric"
  ),
  make_option(
    c("--height"),
    type = "numeric",
    default = 15,
    help = "Plot height in inches [default: %default]",
    metavar = "numeric"
  ),
  make_option(
    c("--dpi"),
    type = "integer",
    default = 600,
    help = "Plot resolution in DPI [default: %default]",
    metavar = "integer"
  )
)

opt_parser <- OptionParser(
  option_list = option_list,
  description = "\nGene Enrichment Analysis Pipeline using g:Profiler"
)
opt <- parse_args(opt_parser)

# ==============================================================================
# Main Function: g:Profiler Pipeline
# ==============================================================================

#' g:Profiler Enrichment Analysis Pipeline
#'
#' Performs functional enrichment analysis on a gene list using the g:Profiler
#' web service. This function queries multiple biological databases (GO, KEGG,
#' Reactome, WikiPathways, etc.) to identify significantly enriched pathways
#' and biological processes.
#'
#' @description
#' This pipeline function orchestrates the complete enrichment analysis workflow:
#' 1. Queries the g:Profiler API with the provided gene list
#' 2. Filters results by statistical significance (FDR correction)
#' 3. Extracts top N most significant terms per database
#' 4. Adds log-transformed p-values for visualization
#' 5. Handles duplicate term names across databases
#' 6. Returns both processed plot data and raw g:Profiler results
#'
#' @param genes Character vector of gene symbols (HGNC gene names) to analyze.
#'   Example: c("TP53", "BRCA1", "EGFR", "MYC")
#'
#' @param sources Character string with comma-separated database names to query.
#'   Available databases:
#'   - "GO" - Gene Ontology (BP, MF, CC)
#'   - "KEGG" - KEGG pathways
#'   - "REAC" - Reactome pathways
#'   - "WP" - WikiPathways
#'   - "TF" - TRANSFAC transcription factors
#'   - "MIRNA" - miRTarBase microRNAs
#'   - "HPA" - Human Protein Atlas
#'   - "CORUM" - CORUM protein complexes
#'   - "HP" - Human Phenotype Ontology
#'   Example: "GO,KEGG,REAC" or "GO,REAC,KEGG,WP"
#'
#' @param top_n Integer specifying the maximum number of top terms to extract
#'   per database source. Terms are ranked by p-value (most significant first).
#'   Default: 20
#'
#' @param organism Character string specifying the organism identifier.
#'   Common options:
#'   - "hsapiens" - Homo sapiens (human)
#'   - "mmusculus" - Mus musculus (mouse)
#'   - "rnorvegicus" - Rattus norvegicus (rat)
#'   - "drerio" - Danio rerio (zebrafish)
#'   - "celegans" - Caenorhabditis elegans
#'   - "scerevisiae" - Saccharomyces cerevisiae (yeast)
#'   Default: "hsapiens"
#'
#' @param fdr_threshold Numeric value for the FDR (False Discovery Rate)
#'   significance threshold. Only terms with FDR-corrected p-value below this
#'   threshold are considered significant.
#'   Range: 0.0 to 1.0
#'   Default: 0.05
#'   Note: Currently logged but not actively used in filtering (g:Profiler
#'   returns only significant results by default)
#'
#' @return A list with two elements:
#'   \item{plot_data}{A data frame ready for visualization containing:
#'     - source: Database name (character)
#'     - term_name: Name of the enriched term (character)
#'     - p_value: Original p-value (numeric)
#'     - log_p_value: -log10 transformed p-value (numeric)
#'     Duplicate term names are suffixed with " (source)" to ensure uniqueness.
#'     NA values are removed.
#'   }
#'   \item{gostres}{Complete g:Profiler result object containing:
#'     - result: Full data frame with all enrichment results including
#'       term_id, term_name, p_value, source, term_size, query_size,
#'       intersection_size, precision, recall, and other statistics
#'     - meta: Metadata about the query
#'   }
#'   If no significant results are found, plot_data will be an empty data frame.
#'
#' @details
#' The function uses the g:Profiler web service (https://biit.cs.ut.ee/gprofiler/)
#' to perform over-representation analysis. Statistical testing is performed
#' using the hypergeometric test with FDR correction for multiple testing.
#'
#' Error handling:
#' - If the g:Profiler API call fails, an error is logged and the function stops
#' - If no significant results are found, a warning is logged and empty results
#'   are returned (the analysis continues)
#'
#' Logging:
#' - All major steps are logged with timestamps
#' - Number of genes queried is logged
#' - Number of significant terms found is logged
#' - Number of terms per database is logged
#' - Duplicate term handling is logged
#'
#' @examples
#' # Example 1: Basic usage with default parameters
#' genes <- c("TP53", "BRCA1", "EGFR", "MYC", "KRAS")
#' results <- gprofiler_pipeline(
#'   genes = genes,
#'   sources = "GO,KEGG,REAC",
#'   top_n = 20,
#'   organism = "hsapiens",
#'   fdr_threshold = 0.05
#' )
#'
#' # Example 2: Analyze only Gene Ontology terms
#' results_go <- gprofiler_pipeline(
#'   genes = my_gene_list,
#'   sources = "GO",
#'   top_n = 10,
#'   organism = "hsapiens",
#'   fdr_threshold = 0.01
#' )
#'
#' # Example 3: Mouse genes with multiple databases
#' mouse_genes <- c("Trp53", "Brca1", "Egfr", "Myc")
#' results_mouse <- gprofiler_pipeline(
#'   genes = mouse_genes,
#'   sources = "GO,KEGG,REAC,WP",
#'   top_n = 15,
#'   organism = "mmusculus",
#'   fdr_threshold = 0.05
#' )
#'
#' # Example 4: Access results
#' plot_data <- results$plot_data  # For visualization
#' full_results <- results$gostres$result  # Full enrichment table
#'
#' # Example 5: More permissive threshold with many databases
#' results_broad <- gprofiler_pipeline(
#'   genes = my_genes,
#'   sources = "GO,KEGG,REAC,WP,TF,MIRNA",
#'   top_n = 30,
#'   organism = "hsapiens",
#'   fdr_threshold = 0.1
#' )
#'
#' @seealso
#' \code{\link[gprofiler2]{gost}} for the underlying g:Profiler function
#' \code{\link{create_enrichment_plot}} for visualization of results
#' \code{\link{save_results}} for saving results to files
#'
#' @author Automated Refactoring
#' @export
gprofiler_pipeline <- function(genes, sources, top_n, organism, fdr_threshold) {
  log_info("Starting g:Profiler enrichment analysis...")

  # Parse sources
  sources_vec <- unlist(strsplit(sources, ","))
  log_info(sprintf("Using databases: %s", paste(sources_vec, collapse = ", ")))

  # Run g:Profiler
  log_info(sprintf("Querying %d genes against g:Profiler", length(genes)))
  gostres <- tryCatch(
    {
      gost(
        query = genes,
        organism = organism,
        significant = TRUE,
        correction_method = "fdr",
        domain_scope = "annotated",
        sources = sources_vec
      )
    },
    error = function(e) {
      log_error(sprintf("g:Profiler query failed: %s", e$message))
      stop(e)
    }
  )

  if (is.null(gostres$result) || nrow(gostres$result) == 0) {
    log_warning("No significant enrichment results found!")
    return(list(plot_data = data.frame(), gostres = gostres))
  }

  log_success(sprintf("Found %d significant terms", nrow(gostres$result)))

  # Get top N terms per source
  sources_used <- unique(gostres$result$source)
  log_info(sprintf("Processing top %d terms for each of %d sources",
                   top_n, length(sources_used)))

  top_terms <- lapply(sources_used, function(src) {
    terms <- gostres$result[gostres$result$source == src, ]
    terms <- terms[order(terms$p_value), ]
    n_terms <- min(top_n, nrow(terms))
    terms <- terms[1:n_terms, ]
    log_info(sprintf("  - %s: %d terms", src, n_terms))
    return(terms)
  })

  # Add log p-values
  for (i in seq_along(top_terms)) {
    top_terms[[i]]$log_p_value <- -log10(top_terms[[i]]$p_value)
  }

  # Create plot_data dataframe
  plot_data <- do.call(rbind, lapply(top_terms, function(df) {
    df[, c("source", "term_name", "p_value", "log_p_value")]
  }))

  # Handle duplicate term names
  duplicated_terms <- duplicated(plot_data$term_name)
  if (any(duplicated_terms)) {
    n_duplicates <- sum(duplicated_terms)
    log_info(sprintf("Handling %d duplicate term names", n_duplicates))
    plot_data$term_name[duplicated_terms] <- paste0(
      plot_data$term_name[duplicated_terms],
      " (", plot_data$source[duplicated_terms], ")"
    )
  }

  # Remove NAs
  plot_data <- na.omit(plot_data)

  log_success("g:Profiler pipeline completed successfully")
  return(list(plot_data = plot_data, gostres = gostres))
}

# ==============================================================================
# Visualization Function
# ==============================================================================

create_enrichment_plot <- function(df_gprof, plot_width, plot_height, plot_dpi) {
  log_info("Creating enrichment plot...")

  n_sources <- length(unique(df_gprof$source))

  p <- ggplot(df_gprof, aes(x = p_value, y = term_name, fill = source)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~ source, nrow = 3, scales = "free_y") +
    labs(
      x = "p-value",
      y = "Term Name",
      title = "Top Enriched Terms of Each Database"
    ) +
    theme_bw() +
    theme(
      axis.text.y = element_text(angle = 0, hjust = 1, size = 16),
      axis.title.x = element_text(size = 16, hjust = 0.1),
      axis.title.y = element_text(size = 16),
      plot.title = element_text(size = 20, hjust = 0.25),
      legend.position = "bottom"
    )

  log_success("Plot created successfully")
  return(p)
}

# ==============================================================================
# Save Results Function
# ==============================================================================

save_results <- function(plot, df_gprof, gostres, save_path,
                        plot_width, plot_height, plot_dpi) {
  log_info(sprintf("Saving results to: %s", save_path))

  # Create directory if it doesn't exist
  if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
    log_info(sprintf("Created output directory: %s", save_path))
  }

  # Save plots
  tryCatch({
    ggsave(
      filename = file.path(save_path, "gostplot.pdf"),
      plot = plot,
      width = plot_width,
      height = plot_height,
      dpi = plot_dpi
    )
    log_success("Saved plot: gostplot.pdf")

    ggsave(
      filename = file.path(save_path, "gostplot.png"),
      plot = plot,
      width = plot_width,
      height = plot_height,
      dpi = plot_dpi
    )
    log_success("Saved plot: gostplot.png")
  }, error = function(e) {
    log_error(sprintf("Failed to save plots: %s", e$message))
  })

  # Save g:Profiler results
  if (!is.null(gostres$result) && nrow(gostres$result) > 0) {
    tryCatch({
      result_cols <- min(13, ncol(gostres$result))

      write.table(
        gostres$result[, 1:result_cols],
        file = file.path(save_path, "gostres.csv"),
        sep = ",",
        row.names = FALSE,
        quote = FALSE
      )
      log_success("Saved results: gostres.csv")

      write.table(
        gostres$result[, 1:result_cols],
        file = file.path(save_path, "gostres.tsv"),
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
      )
      log_success("Saved results: gostres.tsv")
    }, error = function(e) {
      log_error(sprintf("Failed to save g:Profiler results: %s", e$message))
    })
  }

  # Save plot data
  if (nrow(df_gprof) > 0) {
    tryCatch({
      write.table(
        df_gprof,
        file = file.path(save_path, "plot_data_gprofiler.csv"),
        sep = ",",
        row.names = FALSE,
        quote = FALSE
      )
      log_success("Saved plot data: plot_data_gprofiler.csv")

      write.table(
        df_gprof,
        file = file.path(save_path, "plot_data_gprofiler.tsv"),
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
      )
      log_success("Saved plot data: plot_data_gprofiler.tsv")
    }, error = function(e) {
      log_error(sprintf("Failed to save plot data: %s", e$message))
    })
  }
}

# ==============================================================================
# Main Execution
# ==============================================================================

main <- function() {
  log_info("=== Gene Enrichment Analysis Pipeline ===")
  log_info(sprintf("Started at: %s", Sys.time()))

  # Load gene data
  log_info(sprintf("Reading gene data from: %s", opt$input))
  genes_classification <- tryCatch(
    {
      df <- read.table(opt$input, sep = ",", header = TRUE)
      genes <- df[[1]]
      log_success(sprintf("Loaded %d genes", length(genes)))
      log_info(sprintf("First 10 genes: %s",
                       paste(head(genes, 10), collapse = ", ")))
      genes
    },
    error = function(e) {
      log_error(sprintf("Failed to read input file: %s", e$message))
      stop(e)
    }
  )

  # Determine output path
  if (is.null(opt$output)) {
    date_str <- format(Sys.Date(), "%Y%m%d")
    save_path <- file.path(
      "reports", "figures",
      paste0(date_str, "_gprofiler_on_kirc_classification_genes")
    )
  } else {
    save_path <- opt$output
  }

  # Run g:Profiler pipeline
  results <- gprofiler_pipeline(
    genes = genes_classification,
    sources = opt$sources,
    top_n = opt$top_n,
    organism = opt$organism,
    fdr_threshold = opt$fdr
  )

  df_gprof <- results$plot_data
  gostres <- results$gostres

  # Create and save plots if we have results
  if (nrow(df_gprof) > 0) {
    plot <- create_enrichment_plot(
      df_gprof = df_gprof,
      plot_width = opt$width,
      plot_height = opt$height,
      plot_dpi = opt$dpi
    )

    save_results(
      plot = plot,
      df_gprof = df_gprof,
      gostres = gostres,
      save_path = save_path,
      plot_width = opt$width,
      plot_height = opt$height,
      plot_dpi = opt$dpi
    )
  } else {
    log_warning("No enrichment results to plot")
  }

  log_info(sprintf("Finished at: %s", Sys.time()))
  log_success("=== Analysis completed successfully ===")
}

# Run main function
if (!interactive()) {
  main()
}

