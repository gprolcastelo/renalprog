"""
Step 6b: Generate Pathway Enrichment Heatmap

This script generates multiple pathway enrichment heatmaps from existing GSEA results.
It creates five heatmaps showing the sum of NES for significant pathways (FDR < 0.05):

1. Top 50 most changing pathways (comparing first and last timepoint)
2. Top 50 most upregulated pathways (with YlGn colormap)
3. Top 50 most downregulated pathways (with YlOrBr colormap)
4. Selected high-level pathways (29 pathways from Reactome highest level)
5. Selected literature pathways (33 pathways from literature review)

This can be run independently after Step 6 enrichment analysis, or to regenerate
the heatmaps with different parameters.

Author: Renalprog Team
Date: December 19, 2025
License: Apache 2.0

Example usage:
    # Generate from existing enrichment results
    python scripts/pipeline_steps/6b_generate_pathway_heatmap.py \\
        --enrichment_file data/processed/20251217_enrichment/trajectory_enrichment.csv \\
        --output_dir data/processed/20251217_enrichment \\
        --fdr_threshold 0.05

    # Use different FDR threshold
    python scripts/pipeline_steps/6b_generate_pathway_heatmap.py \\
        --enrichment_file data/processed/20251217_enrichment/trajectory_enrichment.csv \\
        --output_dir data/processed/20251217_enrichment \\
        --fdr_threshold 0.01 \\
        --yticks_fontsize 14
"""

import argparse
import logging
import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from renalprog.enrichment import generate_pathway_heatmap
import pandas as pd

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Generate multiple pathway enrichment heatmaps from GSEA results',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    parser.add_argument(
        '--enrichment_file',
        type=str,
        required=True,
        help='Path to trajectory_enrichment.csv file from Step 6'
    )

    parser.add_argument(
        '--output_dir',
        type=str,
        default=None,
        help='Output directory for heatmaps (default: same as enrichment file)'
    )

    parser.add_argument(
        '--fdr_threshold',
        type=float,
        default=0.05,
        help='FDR q-value threshold for significance (default: 0.05)'
    )

    parser.add_argument(
        '--colorbar',
        action='store_true',
        default=True,
        help='Show colorbar (default: True)'
    )

    parser.add_argument(
        '--legend',
        action='store_true',
        default=False,
        help='Show legend (default: False)'
    )

    parser.add_argument(
        '--yticks_fontsize',
        type=int,
        default=12,
        help='Font size for y-axis tick labels (default: 12)'
    )

    parser.add_argument(
        '--show',
        action='store_true',
        default=False,
        help='Display the plots interactively (default: False)'
    )

    return parser.parse_args()


def main():
    """Main execution function."""
    args = parse_args()

    # Validate input file
    enrichment_file = Path(args.enrichment_file)
    if not enrichment_file.exists():
        logger.error(f"Enrichment file not found: {enrichment_file}")
        sys.exit(1)

    # Set output directory
    if args.output_dir is None:
        output_dir = enrichment_file.parent
    else:
        output_dir = Path(args.output_dir)

    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 80)
    logger.info("PATHWAY ENRICHMENT HEATMAP GENERATION")
    logger.info("=" * 80)
    logger.info(f"Input file: {enrichment_file}")
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"FDR threshold: {args.fdr_threshold}")
    logger.info(f"Y-axis fontsize: {args.yticks_fontsize}")
    logger.info("=" * 80)

    try:
        # Load enrichment data
        logger.info("Loading enrichment data...")
        enrichment_df = pd.read_csv(enrichment_file)
        logger.info(f"Loaded {enrichment_df.shape[0]} enrichment results")

        # Validate required columns
        required_cols = ['Patient', 'Idx', 'Transition', 'NAME', 'NES', 'FDR q-val']
        missing_cols = [col for col in required_cols if col not in enrichment_df.columns]
        if missing_cols:
            logger.error(f"Missing required columns: {missing_cols}")
            logger.error(f"Available columns: {enrichment_df.columns.tolist()}")
            sys.exit(1)

        # Generate heatmaps
        logger.info("Generating pathway enrichment heatmaps...")
        heatmap_data, figures_dict = generate_pathway_heatmap(
            enrichment_df=enrichment_df,
            output_dir=output_dir,
            fdr_threshold=args.fdr_threshold,
            colorbar=args.colorbar,
            legend=args.legend,
            yticks_fontsize=args.yticks_fontsize,
            show=args.show
        )

        logger.info("=" * 80)
        logger.info("HEATMAP GENERATION COMPLETE")
        logger.info(f"Significant pathways: {heatmap_data.shape[0]}")
        logger.info(f"Timepoints: {heatmap_data.shape[1]}")
        logger.info(f"Heatmaps generated: {len(figures_dict)}")
        logger.info(f"Results saved to: {output_dir}")
        logger.info("  - pathway_nes_summary.csv: Full summary data (pathways × timepoints)")
        logger.info("  - top50_most_changing_pathways.{pdf,png,svg}: Top 50 pathways with largest change")
        logger.info("  - top50_most_upregulated_pathways.{pdf,png,svg}: Top 50 upregulated pathways")
        logger.info("  - top50_most_downregulated_pathways.{pdf,png,svg}: Top 50 downregulated pathways")
        logger.info("  - selected_pathways_highest_level.{pdf,png,svg}: 29 high-level Reactome pathways")
        logger.info("  - selected_pathways_literature.{pdf,png,svg}: 33 literature pathways")
        logger.info("=" * 80)

    except Exception as e:
        logger.error(f"Failed to generate heatmap: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()

