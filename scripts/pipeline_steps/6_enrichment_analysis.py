"""
Step 6: Dynamic Enrichment Analysis Pipeline with PyDESeq2

This script orchestrates the dynamic enrichment analysis pipeline:

1. **PyDESeq2 Differential Expression Analysis**
   - Converts log2(RSEM+1) trajectory data back to RSEM integer counts
   - Runs PyDESeq2 for each trajectory timepoint vs healthy controls
   - Generates statistically valid log2FoldChange rankings

2. **GSEA Pathway Enrichment**
   - Creates ranked gene lists (.rnk files) from DESeq2 results
   - Executes GSEA in parallel using ReactomePathways gene sets
   - Computes pathway enrichment scores (NES, p-values, FDR q-values)

3. **Result Aggregation**
   - Combines GSEA results across all trajectories and timepoints
   - Generates final enrichment dataset

4. **Visualization**
   - Generates pathway enrichment heatmaps showing NES across pseudo-time

IMPORTANT: This pipeline uses PyDESeq2 (not simple fold-change) for proper
differential expression analysis. This ensures biologically valid and
publication-quality results.

Dependencies:
- PyDESeq2 (installed via requirements.txt)
- GSEA CLI tool (download from https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp)
- ReactomePathways.gmt file (in data/external/)

Author: Renalprog Team
Date: December 19, 2025
License: Apache 2.0

Example usage:
    python scripts/pipeline_steps/6_enrichment_analysis.py \\
        --trajectory_dir data/interim/20251216_synthetic_data/kirc/recnet/early_to_late/train_to_train/early_to_late \\
        --output_dir data/processed/20251219_enrichment \\
        --cancer_type kirc \\
        --n_threads 8 \\
        --gsea_path ./GSEA_4.3.2/gsea-cli.sh \\
        --pathways_file data/external/ReactomePathways.gmt
"""
import os
os.environ["NUMEXPR_MAX_THREADS"] = "112"
import argparse
import logging
import sys
from pathlib import Path
from datetime import datetime

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from renalprog.enrichment import (
    EnrichmentPipeline,
    process_trajectories_for_deseq,
    run_gsea_parallel,
    combine_gsea_results
)
from renalprog.config import PATHS
from renalprog.utils import set_seed

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Run dynamic enrichment analysis on synthetic trajectories',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    parser.add_argument(
        '--trajectory_dir',
        type=str,
        required=True,
        help='Directory containing trajectory CSV files'
    )

    parser.add_argument(
        '--output_dir',
        type=str,
        default=None,
        help='Output directory for enrichment results (default: data/processed/YYYYMMDD_enrichment)'
    )

    parser.add_argument(
        '--cancer_type',
        type=str,
        choices=['kirc', 'lobular', 'ductal'],
        default='kirc',
        help='Cancer type (default: kirc)'
    )

    parser.add_argument(
        '--data_dir',
        type=str,
        default=None,
        help='Path to preprocessed RNA-seq data (default: latest preprocessed KIRC data)'
    )

    parser.add_argument(
        '--metadata_dir',
        type=str,
        default=None,
        help='Path to clinical metadata (default: latest preprocessed KIRC data)'
    )

    parser.add_argument(
        '--control_data_dir',
        type=str,
        default=None,
        help='Path to control RNA-seq data (default: data/processed/controls/KIRC/rnaseq_control.csv)'
    )

    parser.add_argument(
        '--control_metadata_dir',
        type=str,
        default=None,
        help='Path to control metadata (default: data/processed/controls/KIRC/clinical_control.csv)'
    )

    parser.add_argument(
        '--gsea_path',
        type=str,
        default='./GSEA_4.3.2/gsea-cli.sh',
        help='Path to GSEA CLI tool (default: ./GSEA_4.3.2/gsea-cli.sh)'
    )

    parser.add_argument(
        '--pathways_file',
        type=str,
        default='data/external/ReactomePathways.gmt',
        help='Path to pathways GMT file (default: data/external/ReactomePathways.gmt)'
    )

    parser.add_argument(
        '--n_threads',
        type=int,
        default=4,
        help='Number of parallel threads for GSEA (default: 4)'
    )

    parser.add_argument(
        '--seed',
        type=int,
        default=42,
        help='Random seed for reproducibility (default: 42)'
    )

    parser.add_argument(
        '--skip_deseq',
        action='store_true',
        help='Skip DESeq processing (use if already completed)'
    )

    parser.add_argument(
        '--skip_gsea',
        action='store_true',
        help='Skip GSEA analysis (use if already completed)'
    )

    parser.add_argument(
        '--cleanup',
        action='store_true',
        help='Remove intermediate files after processing'
    )

    return parser.parse_args()


def main():
    """Main execution function."""
    args = parse_args()

    # Set random seed
    set_seed(args.seed)

    # Set default output directory
    if args.output_dir is None:
        today = datetime.now().strftime("%Y%m%d")
        args.output_dir = str(PATHS['processed'] / f"{today}_enrichment")

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 80)
    logger.info("ENRICHMENT ANALYSIS PIPELINE")
    logger.info("=" * 80)
    logger.info(f"Trajectory directory: {args.trajectory_dir}")
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Cancer type: {args.cancer_type}")
    logger.info(f"Number of threads: {args.n_threads}")
    logger.info(f"GSEA path: {args.gsea_path}")
    logger.info(f"Pathways file: {args.pathways_file}")
    logger.info("=" * 80)

    # Validate GSEA installation
    gsea_path = Path(args.gsea_path)
    if not gsea_path.exists():
        logger.error(f"GSEA CLI not found at {gsea_path}")
        logger.error("Please download GSEA from: https://www.gsea-msigdb.org/gsea/downloads.jsp")
        logger.error("And extract to the project root or specify path with --gsea_path")
        sys.exit(1)

    # Validate pathways file
    pathways_path = Path(args.pathways_file)
    if not pathways_path.exists():
        logger.error(f"Pathways file not found at {pathways_path}")
        logger.error("Please ensure ReactomePathways.gmt is in data/external/")
        sys.exit(1)

    # Initialize pipeline
    pipeline = EnrichmentPipeline(
        trajectory_dir=args.trajectory_dir,
        output_dir=output_dir,
        cancer_type=args.cancer_type,
        data_dir=args.data_dir,
        metadata_dir=args.metadata_dir,
        control_data_dir=args.control_data_dir,
        control_metadata_dir=args.control_metadata_dir,
        gsea_path=args.gsea_path,
        pathways_file=args.pathways_file,
        n_threads=args.n_threads
    )

    # Run pipeline
    try:
        pipeline.run(
            skip_deseq=args.skip_deseq,
            skip_gsea=args.skip_gsea,
            cleanup=args.cleanup
        )

        logger.info("=" * 80)
        logger.info("ENRICHMENT ANALYSIS COMPLETE")
        logger.info(f"Results saved to: {output_dir}")
        logger.info("=" * 80)

    except Exception as e:
        logger.error(f"Pipeline failed: {e}", exc_info=True)

        # Provide helpful guidance for common errors
        if "No GSEA results found" in str(e):
            logger.error("")
            logger.error("=" * 80)
            logger.error("TROUBLESHOOTING GSEA RESULTS ERROR")
            logger.error("=" * 80)
            logger.error("The enrichment pipeline could not find GSEA output files.")
            logger.error("")
            logger.error("To diagnose the issue, run:")
            logger.error(f"  python scripts/debug_gsea_results.py --deseq_dir {output_dir / 'deseq'}")
            logger.error("")
            logger.error("For detailed troubleshooting steps, see:")
            logger.error("  docs/GSEA_TROUBLESHOOTING.md")
            logger.error("")
            logger.error("Common causes:")
            logger.error("  1. GSEA commands failed to run (check failed_gsea_commands.txt)")
            logger.error("  2. GSEA not installed or wrong path")
            logger.error("  3. Java version incompatible")
            logger.error("  4. GSEA output structure different than expected")
            logger.error("=" * 80)

        sys.exit(1)


if __name__ == "__main__":
    main()

