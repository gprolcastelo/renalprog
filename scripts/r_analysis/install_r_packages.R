#!/usr/bin/env Rscript

# ==============================================================================
# Install R Dependencies for renalprog Package
# ==============================================================================
# Description: Installs required R packages for gene enrichment analysis
# Author: Automated Setup
# Date: 2025-12-17
# ==============================================================================

cat("================================================================================\n")
cat("Installing R dependencies for renalprog\n")
cat("================================================================================\n\n")

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org/"))

# Function to install package if not already installed
install_if_missing <- function(package, source = "CRAN", bioc_version = NULL) {
  if (!require(package, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("Installing %s from %s...\n", package, source))

    if (source == "CRAN") {
      install.packages(package, dependencies = TRUE, quiet = FALSE)
    } else if (source == "Bioconductor") {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      if (!is.null(bioc_version)) {
        BiocManager::install(version = bioc_version, ask = FALSE)
      }
      BiocManager::install(package, ask = FALSE, update = FALSE)
    }

    # Verify installation
    if (require(package, character.only = TRUE, quietly = TRUE)) {
      cat(sprintf("✓ Successfully installed %s\n\n", package))
      return(TRUE)
    } else {
      cat(sprintf("✗ Failed to install %s\n\n", package))
      return(FALSE)
    }
  } else {
    cat(sprintf("✓ %s is already installed\n\n", package))
    return(TRUE)
  }
}

# ==============================================================================
# Install Required Packages
# ==============================================================================

cat("Installing CRAN packages...\n")
cat("--------------------------------------------------------------------------------\n")

# CRAN packages
cran_packages <- c(
  "gprofiler2",  # For g:Profiler enrichment analysis
  "ggplot2",     # For plotting
  "optparse"     # For command-line argument parsing
)

cran_success <- sapply(cran_packages, function(pkg) {
  install_if_missing(pkg, source = "CRAN")
})

# ==============================================================================
# Summary
# ==============================================================================

cat("\n================================================================================\n")
cat("Installation Summary\n")
cat("================================================================================\n")

all_packages <- c(cran_packages)
all_success <- c(cran_success)

cat(sprintf("\nTotal packages: %d\n", length(all_packages)))
cat(sprintf("Successfully installed: %d\n", sum(all_success)))
cat(sprintf("Failed: %d\n", sum(!all_success)))

if (sum(!all_success) > 0) {
  cat("\nFailed packages:\n")
  failed_pkgs <- all_packages[!all_success]
  for (pkg in failed_pkgs) {
    cat(sprintf("  - %s\n", pkg))
  }
  cat("\nPlease install these packages manually.\n")
  quit(status = 1)
} else {
  cat("\n✓ All R dependencies installed successfully!\n")
  cat("\nYou can now run: Rscript scripts/r_analysis/gene_enrichment.R\n")
}

cat("================================================================================\n")

