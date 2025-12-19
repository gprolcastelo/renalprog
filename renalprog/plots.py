"""
Visualization functions for renalprog using Plotly.

All plots are saved in multiple formats: HTML (interactive), PNG, PDF, and SVG.

Contains plotting functions for:
- VAE latent space visualization
- Training history plots
- Trajectory plots
- Gene expression heatmaps
- Enrichment analysis results
- Paper figures
"""

import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from pathlib import Path
from typing import Optional, List, Tuple, Union, Dict
import logging, umap, os
from datetime import datetime

logger = logging.getLogger(__name__)

# Default Plotly settings
DEFAULT_TEMPLATE = "plotly_white"
DEFAULT_COLORSCALE = "Viridis"
DEFAULT_WIDTH = 800
DEFAULT_HEIGHT = 600


def save_plot(
    fig: go.Figure,
    save_path: Union[str, Path],
    formats: List[str] = ["html", "png", "pdf", "svg"],
    width: int = DEFAULT_WIDTH,
    height: int = DEFAULT_HEIGHT,
) -> None:
    """
    Save plotly figure in multiple formats.

    Args:
        fig: Plotly figure object
        save_path: Base path for saving (without extension)
        formats: List of formats to save ['html', 'png', 'pdf', 'svg']
        width: Width in pixels for static formats
        height: Height in pixels for static formats
    """
    save_path = Path(save_path)
    save_path.parent.mkdir(parents=True, exist_ok=True)

    # Remove extension if present
    base_path = save_path.with_suffix('')

    for fmt in formats:
        output_path = base_path.with_suffix(f'.{fmt}')
        try:
            if fmt == 'html':
                fig.write_html(str(output_path))
            elif fmt in ['png', 'pdf', 'svg']:
                fig.write_image(str(output_path), width=width, height=height)
            logger.info(f"Saved plot to {output_path}")
        except Exception as e:
            logger.warning(f"Failed to save {fmt} format: {e}")
            if fmt in ['png', 'pdf', 'svg']:
                logger.warning("Note: Static image export requires kaleido package: pip install kaleido")


def plot_latent_space(
    latent_coords: np.ndarray,
    labels: np.ndarray,
    method: str = "UMAP",
    save_path: Optional[Path] = None,
    title: Optional[str] = None,
    color_discrete_map: Optional[Dict] = None,
    **kwargs
) -> go.Figure:
    """
    Plot 2D latent space representation colored by labels.

    Args:
        latent_coords: 2D array of latent coordinates (n_samples x 2)
        labels: Array of labels for coloring points
        method: Dimensionality reduction method name for title
        save_path: Optional path to save figure
        title: Optional custom title
        color_discrete_map: Optional color mapping for labels
        **kwargs: Additional arguments for scatter plot

    Returns:
        Plotly Figure object
    """
    if latent_coords.shape[1] != 2:
        raise ValueError(f"Expected 2D coordinates, got shape {latent_coords.shape}")

    df = pd.DataFrame({
        f'{method}_1': latent_coords[:, 0],
        f'{method}_2': latent_coords[:, 1],
        'Label': labels
    })

    fig = px.scatter(
        df,
        x=f'{method}_1',
        y=f'{method}_2',
        color='Label',
        color_discrete_map=color_discrete_map,
        title=title or f'{method} Latent Space Visualization',
        template=DEFAULT_TEMPLATE,
        **kwargs
    )

    fig.update_layout(
        width=DEFAULT_WIDTH,
        height=DEFAULT_HEIGHT,
        font=dict(size=12)
    )

    if save_path:
        save_plot(fig, save_path)

    return fig


def plot_training_history(
    history: Dict[str, List[float]],
    save_path: Optional[Path] = None,
    title: str = "Training History",
    log_scale: bool = False,
) -> go.Figure:
    """
    Plot training and validation losses over epochs.

    Args:
        history: Dictionary with 'train_loss' and 'val_loss' keys
        save_path: Optional path to save figure
        title: Plot title
        log_scale: Whether to use log scale for y-axis

    Returns:
        Plotly Figure object
    """
    epochs = list(range(1, len(history['train_loss']) + 1))

    fig = go.Figure()

    ## Total loss = KL + reconstruction

    # Train loss
    fig.add_trace(go.Scatter(
        x=epochs,
        y=history['train_loss'],
        mode='lines',
        name='Train Loss',
        line=dict(color='#1f77b4', width=2),
        marker=dict(size=4)
    ))

    # Validation loss
    if 'val_loss' in history:
        fig.add_trace(go.Scatter(
            x=epochs,
            y=history['val_loss'],
            mode='lines',
            name='Val Loss',
            line=dict(color='#ff7f0e', width=2),
            marker=dict(size=4)
        ))

    fig.update_layout(
        title=title,
        xaxis_title='Epoch',
        yaxis_title='Loss',
        yaxis_type='log' if log_scale else 'linear',
        template=DEFAULT_TEMPLATE,
        width=DEFAULT_WIDTH,
        height=DEFAULT_HEIGHT,
        hovermode='x unified'
    )



    # KL Divergence plot
    fig_kl = go.Figure()
    if 'train_kl_loss' in history and 'val_kl_loss' in history:
        fig_kl.add_trace(go.Scatter(
            x=epochs,
            y=history['train_kl_loss'],
            mode='lines',
            name='Train KL Divergence',
            line=dict(color='#1f77b4', width=2),
            marker=dict(size=4)
        ))
        fig_kl.add_trace(go.Scatter(
            x=epochs,
            y=history['val_kl_loss'],
            mode='lines',
            name='Val KL Divergence',
            line=dict(color='#ff7f0e', width=2),
            marker=dict(size=4)
        ))

    # Reconstruction Loss plot
    fig_rec = go.Figure()
    if 'train_recon_loss' in history and 'val_recon_loss' in history:
        fig_rec.add_trace(go.Scatter(
            x=epochs,
            y=history['train_recon_loss'],
            mode='lines',
            name='Train Reconstruction Loss',
            line=dict(color='#1f77b4', width=2),
            marker=dict(size=4)
        ))
        fig_rec.add_trace(go.Scatter(
            x=epochs,
            y=history['val_recon_loss'],
            mode='lines',
            name='Val Reconstruction Loss',
            line=dict(color='#ff7f0e', width=2),
            marker=dict(size=4)
        ))

    if save_path:
        save_plot(fig, save_path / "total_loss")
        save_plot(fig_kl, save_path / "kl_divergence")
        save_plot(fig_rec, save_path / "reconstruction_loss")

    return fig


def plot_reconstruction_losses(
    loss_train: List[float],
    loss_test: List[float],
    save_path: Optional[Path] = None,
    title: str = "Reconstruction Network Losses",
) -> go.Figure:
    """
    Plot training and test losses for reconstruction network.

    Args:
        loss_train: List of training losses
        loss_test: List of test losses
        save_path: Optional path to save figure
        title: Plot title

    Returns:
        Plotly Figure object
    """
    epochs = list(range(1, len(loss_train) + 1))

    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=epochs,
        y=loss_train,
        mode='lines',
        name='Train',
        line=dict(color='#1f77b4', width=2)
    ))

    fig.add_trace(go.Scatter(
        x=epochs,
        y=loss_test,
        mode='lines',
        name='Test',
        line=dict(color='#ff7f0e', width=2)
    ))

    fig.update_layout(
        title=title,
        xaxis_title='Epoch',
        yaxis_title='Loss',
        template=DEFAULT_TEMPLATE,
        width=DEFAULT_WIDTH,
        height=DEFAULT_HEIGHT,
        hovermode='x unified'
    )

    if save_path:
        save_plot(fig, save_path)

    return fig


def plot_gene_expression_heatmap(
    data: pd.DataFrame,
    save_path: Optional[Path] = None,
    title: str = "Gene Expression Heatmap",
    colorscale: str = DEFAULT_COLORSCALE,
    show_labels: bool = True,
    **kwargs
) -> go.Figure:
    """
    Create a heatmap of gene expression data.

    Args:
        data: DataFrame with genes as columns and samples as rows
        save_path: Optional path to save figure
        title: Plot title
        colorscale: Plotly colorscale name
        show_labels: Whether to show axis labels
        **kwargs: Additional arguments for heatmap

    Returns:
        Plotly Figure object
    """
    fig = go.Figure(data=go.Heatmap(
        z=data.values,
        x=data.columns if show_labels else None,
        y=data.index if show_labels else None,
        colorscale=colorscale,
        **kwargs
    ))

    fig.update_layout(
        title=title,
        template=DEFAULT_TEMPLATE,
        width=max(DEFAULT_WIDTH, len(data.columns) * 10),
        height=max(DEFAULT_HEIGHT, len(data) * 10),
        xaxis_title='Genes' if show_labels else None,
        yaxis_title='Samples' if show_labels else None
    )

    if save_path:
        save_plot(fig, save_path)

    return fig


def plot_trajectory(
    trajectory: np.ndarray,
    gene_names: Optional[List[str]] = None,
    save_path: Optional[Path] = None,
    title: str = "Gene Expression Trajectory",
    n_genes_to_show: int = 20,
) -> go.Figure:
    """
    Plot gene expression changes along a trajectory.

    Args:
        trajectory: Array of shape (n_timepoints, n_genes)
        gene_names: Optional list of gene names
        save_path: Optional path to save figure
        title: Plot title
        n_genes_to_show: Number of top varying genes to display

    Returns:
        Plotly Figure object
    """
    n_timepoints, n_genes = trajectory.shape

    # Calculate variance for each gene
    gene_variance = np.var(trajectory, axis=0)
    top_genes_idx = np.argsort(gene_variance)[-n_genes_to_show:]

    if gene_names is None:
        gene_names = [f'Gene_{i}' for i in range(n_genes)]

    fig = go.Figure()

    timepoints = list(range(n_timepoints))

    for idx in top_genes_idx:
        fig.add_trace(go.Scatter(
            x=timepoints,
            y=trajectory[:, idx],
            mode='lines',
            name=gene_names[idx],
            line=dict(width=1.5)
        ))

    fig.update_layout(
        title=title,
        xaxis_title='Timepoint',
        yaxis_title='Expression Level',
        template=DEFAULT_TEMPLATE,
        width=DEFAULT_WIDTH,
        height=DEFAULT_HEIGHT,
        hovermode='x unified'
    )

    if save_path:
        save_plot(fig, save_path)

    return fig


def plot_pca_variance(
    explained_variance_ratio: np.ndarray,
    save_path: Optional[Path] = None,
    title: str = "PCA Explained Variance",
    n_components: int = 20,
) -> go.Figure:
    """
    Plot explained variance ratio from PCA.

    Args:
        explained_variance_ratio: Array of explained variance ratios
        save_path: Optional path to save figure
        title: Plot title
        n_components: Number of components to show

    Returns:
        Plotly Figure object
    """
    n_show = min(n_components, len(explained_variance_ratio))
    components = list(range(1, n_show + 1))

    # Individual variance
    fig = make_subplots(
        rows=1, cols=2,
        subplot_titles=('Individual Explained Variance', 'Cumulative Explained Variance')
    )

    fig.add_trace(
        go.Bar(
            x=components,
            y=explained_variance_ratio[:n_show],
            name='Individual',
            marker_color='#1f77b4'
        ),
        row=1, col=1
    )

    fig.add_trace(
        go.Scatter(
            x=components,
            y=np.cumsum(explained_variance_ratio[:n_show]),
            mode='lines+markers',
            name='Cumulative',
            line=dict(color='#ff7f0e', width=2),
            marker=dict(size=6)
        ),
        row=1, col=2
    )

    fig.update_xaxes(title_text='Principal Component', row=1, col=1)
    fig.update_xaxes(title_text='Principal Component', row=1, col=2)
    fig.update_yaxes(title_text='Explained Variance Ratio', row=1, col=1)
    fig.update_yaxes(title_text='Cumulative Variance', row=1, col=2)

    fig.update_layout(
        title_text=title,
        template=DEFAULT_TEMPLATE,
        width=DEFAULT_WIDTH * 2,
        height=DEFAULT_HEIGHT,
        showlegend=False
    )

    if save_path:
        save_plot(fig, save_path, width=DEFAULT_WIDTH * 2)

    return fig


def plot_enrichment_results(
    enrichment_df: pd.DataFrame,
    save_path: Optional[Path] = None,
    title: str = "Pathway Enrichment Analysis",
    top_n: int = 20,
    p_value_col: str = 'p_value',
    pathway_col: str = 'pathway',
) -> go.Figure:
    """
    Plot enrichment analysis results.

    Args:
        enrichment_df: DataFrame with enrichment results
        save_path: Optional path to save figure
        title: Plot title
        top_n: Number of top pathways to show
        p_value_col: Column name for p-values
        pathway_col: Column name for pathway names

    Returns:
        Plotly Figure object
    """
    # Sort by p-value and take top N
    df_sorted = enrichment_df.sort_values(p_value_col).head(top_n)

    # Calculate -log10(p-value)
    df_sorted['neg_log_p'] = -np.log10(df_sorted[p_value_col])

    fig = go.Figure(go.Bar(
        x=df_sorted['neg_log_p'],
        y=df_sorted[pathway_col],
        orientation='h',
        marker=dict(
            color=df_sorted['neg_log_p'],
            colorscale='Reds',
            showscale=True,
            colorbar=dict(title='-log10(p-value)')
        )
    ))

    fig.update_layout(
        title=title,
        xaxis_title='-log10(p-value)',
        yaxis_title='Pathway',
        template=DEFAULT_TEMPLATE,
        width=DEFAULT_WIDTH,
        height=max(DEFAULT_HEIGHT, top_n * 25),
        yaxis={'categoryorder': 'total ascending'}
    )

    if save_path:
        save_plot(fig, save_path, height=max(DEFAULT_HEIGHT, top_n * 25))

    return fig


def plot_confusion_matrix(
    confusion_matrix: np.ndarray,
    class_names: List[str],
    save_path: Optional[Path] = None,
    title: str = "Confusion Matrix",
    normalize: bool = False,
) -> go.Figure:
    """
    Plot confusion matrix as heatmap.

    Args:
        confusion_matrix: Square confusion matrix
        class_names: List of class names
        save_path: Optional path to save figure
        title: Plot title
        normalize: Whether to normalize by row

    Returns:
        Plotly Figure object
    """
    if normalize:
        cm = confusion_matrix.astype('float') / confusion_matrix.sum(axis=1)[:, np.newaxis]
        text = [[f'{confusion_matrix[i, j]}<br>({cm[i, j]:.2%})'
                for j in range(len(class_names))]
               for i in range(len(class_names))]
    else:
        cm = confusion_matrix
        text = [[str(confusion_matrix[i, j])
                for j in range(len(class_names))]
               for i in range(len(class_names))]

    fig = go.Figure(data=go.Heatmap(
        z=cm,
        x=class_names,
        y=class_names,
        text=text,
        texttemplate='%{text}',
        colorscale='Blues',
        showscale=True
    ))

    fig.update_layout(
        title=title,
        xaxis_title='Predicted',
        yaxis_title='Actual',
        template=DEFAULT_TEMPLATE,
        width=DEFAULT_WIDTH,
        height=DEFAULT_HEIGHT,
        yaxis={'autorange': 'reversed'}
    )

    if save_path:
        save_plot(fig, save_path)

    return fig


def plot_loss_landscape(
    loss_values: np.ndarray,
    param1_values: np.ndarray,
    param2_values: np.ndarray,
    save_path: Optional[Path] = None,
    title: str = "Loss Landscape",
    param1_name: str = "Parameter 1",
    param2_name: str = "Parameter 2",
) -> go.Figure:
    """
    Plot 2D loss landscape for hyperparameter tuning.

    Args:
        loss_values: 2D array of loss values
        param1_values: Values for first parameter
        param2_values: Values for second parameter
        save_path: Optional path to save figure
        title: Plot title
        param1_name: Name of first parameter
        param2_name: Name of second parameter

    Returns:
        Plotly Figure object
    """
    fig = go.Figure(data=go.Contour(
        z=loss_values,
        x=param1_values,
        y=param2_values,
        colorscale='RdYlBu_r',
        contours=dict(
            showlabels=True,
            labelfont=dict(size=10, color='white')
        ),
        colorbar=dict(title='Loss')
    ))

    fig.update_layout(
        title=title,
        xaxis_title=param1_name,
        yaxis_title=param2_name,
        template=DEFAULT_TEMPLATE,
        width=DEFAULT_WIDTH,
        height=DEFAULT_HEIGHT
    )

    if save_path:
        save_plot(fig, save_path)

    return fig

def plot_umap_plotly(
    data,
    clinical,
    colors_dict,
    shapes_dict=None,
    n_components=2,
    save_fig=False,
    save_as=None,
    seed=None,
    title='UMAP',
    show=True,
    marker_size=8,
):
    """
    Plot UMAP of the data with Plotly using different colors for the different groups.

    Parameters
    ----------
    data : pandas.DataFrame
        Features as rows and samples as columns (same as in plot_umap).
    clinical : pandas.Series
        Category per sample (index must match data.columns).
    colors_dict : dict
        Mapping {group_name: color_hex_or_name}.
    shapes_dict: dict
        Mapping {group_name: shape}.
    n_components : int, optional
        2 or 3, by default 2.
    save_fig : bool, optional
        If True, save HTML/PNG/PDF/SVG, by default False.
    save_as : str or None, optional
        Base path (without extension) for saving, by default None.
    seed : int or None, optional
        Random state for UMAP, by default None.
    title : str, optional
        Plot title, by default 'UMAP'.
    show : bool, optional
        If True, display the plot, by default True.
    """

    # Check number of samples is the first dimension of data:
    if data.shape[0] != clinical.shape[0]:
        data = data.T
        if data.shape[0] != clinical.shape[0]:
            raise ValueError("Data and clinical metadata must have the same number of samples")


    if n_components not in (2, 3):
        raise ValueError("n_components must be 2 or 3 for plot_umap_plotly")

    today = datetime.now().strftime("%Y%m%d")
    if save_as is None:
        suffix = "UMAP" if n_components == 2 else "3D_UMAP"

    if seed is not None:
        umap_ = umap.UMAP(n_components=n_components, random_state=seed)
    else:
        umap_ = umap.UMAP(n_components=n_components)

    # data: samples x features
    X_umap = umap_.fit_transform(data)
    print("X_umap.shape", X_umap.shape)

    # Determine color and shape series from clinical
    if isinstance(clinical, pd.DataFrame):
        color_col = clinical.columns[0]
        color_series = clinical[color_col]
        # use second column for shapes if provided and shapes_dict is given
        if shapes_dict is not None and clinical.shape[1] >= 2:
            shape_col = clinical.columns[1]
            shape_series = clinical[shape_col]
        else:
            shape_series = None
    elif isinstance(clinical, pd.Series):
        color_series = clinical
        shape_series = None
    else:
        raise ValueError("clinical must be a pandas Series or DataFrame")
    print("color_series.shape", color_series.shape)

    # Build plotting DataFrame
    all_patients = data.index.tolist()
    print("len(all_patients)", len(all_patients))
    print("color_series.loc[all_patients].values.shape", color_series.loc[all_patients].values.shape)
    df_plot = pd.DataFrame(
        {
            "sample": all_patients,
            "group": color_series.loc[all_patients].values,
            "UMAP_1": X_umap[:, 0],
            "UMAP_2": X_umap[:, 1],
        }
    )
    if n_components == 3:
        df_plot["UMAP_3"] = X_umap[:, 2]

    # Attach shape column if available
    if shape_series is not None:
        df_plot["shape"] = shape_series.loc[all_patients].values

    # Build color sequence in the order of unique groups
    unique_groups = df_plot["group"].unique()
    color_sequence = [colors_dict[g] for g in unique_groups]

    # Prepare symbol mapping if shapes are used
    symbol_map = None
    if "shape" in df_plot.columns and shapes_dict is not None:
        # convert common Matplotlib markers to Plotly symbols if needed
        matplot_to_plotly = {
            'o': 'circle', 's': 'square', '^': 'triangle-up', 'v': 'triangle-down',
            'D': 'diamond', 'd': 'diamond-wide', 'X': 'x', 'x': 'x', '*': 'star',
            '+': 'cross', 'p': 'pentagon', 'h': 'hexagon', 'H': 'hexagon2'
        }
        unique_shapes = df_plot["shape"].unique()
        symbol_map = {}
        for sh in unique_shapes:
            # get marker definition from shapes_dict; fallback to the value itself
            marker = shapes_dict.get(sh, shapes_dict.get(str(sh), sh))
            # translate matplotlib marker codes to plotly symbol names when possible
            symbol = matplot_to_plotly.get(marker, marker)
            symbol_map[sh] = symbol

    # Create plotly figure with optional symbols
    if n_components == 2:
        fig = px.scatter(
            df_plot,
            x="UMAP_1",
            y="UMAP_2",
            color="group",
            color_discrete_sequence=color_sequence,
            hover_name="sample",
            template="simple_white",
            width=800,
            height=800,
            symbol="shape" if "shape" in df_plot.columns and symbol_map is not None else None,
            symbol_map=symbol_map if symbol_map is not None else None,
        )
        fig.update_layout(
            title=title,
            xaxis_title="UMAP 1",
            yaxis_title="UMAP 2",
        )
    else:
        fig = px.scatter_3d(
            df_plot,
            x="UMAP_1",
            y="UMAP_2",
            z="UMAP_3",
            color="group",
            color_discrete_sequence=color_sequence,
            hover_name="sample",
            template="simple_white",
            width=800,
            height=800,
            symbol="shape" if "shape" in df_plot.columns and symbol_map is not None else None,
            symbol_map=symbol_map if symbol_map is not None else None,
        )
        fig.update_layout(
            title=title,
            scene=dict(
                xaxis_title="UMAP 1",
                yaxis_title="UMAP 2",
                zaxis_title="UMAP 3",
            ),
        )
    fig.update_traces(marker=dict(size=marker_size))
    # Optional saving
    if save_fig:
        base_dir = os.path.dirname(save_as)
        if base_dir and not os.path.exists(base_dir):
            os.makedirs(base_dir, exist_ok=True)
        # Save as HTML
        fig.write_html(f"{save_as}.html")
        # Save static images
        for extension in ['png', 'pdf', 'svg']:
            print(f"Saved UMAP plotly figure to: {save_as}.{extension}")
            fig.write_image(f"{save_as}.{extension}", scale=2)
    if show:
        fig.show()