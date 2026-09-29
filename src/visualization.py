"""
Visualization module for generating publication-quality figures.

This module generates Figures 1-9 for the computational drug repurposing
workflow. Uses matplotlib and seaborn with a consistent scientific style.
"""

import os
import logging
from typing import List, Dict, Optional, Any, Tuple
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns
from scipy.stats import mannwhitneyu

# Setup consistent style
plt.style.use('seaborn-v0_8-whitegrid')
FIG_DPI = 300
FIG_FORMAT = 'png'
COLORS = {
    'positive_control': '#E74C3C', 
    'screening': '#3498DB', 
    'decoy': '#95A5A6', 
    'candidate': '#2ECC71'
}

logger = logging.getLogger(__name__)
if not logger.handlers:
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')


def _save_and_close(fig: plt.Figure, output_path: str) -> str:
    """Helper to save a figure, close it, and log the path."""
    try:
        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        fig.tight_layout()
        fig.savefig(output_path, dpi=FIG_DPI, format=FIG_FORMAT, bbox_inches='tight')
        logger.info(f"Figure saved to {output_path}")
    except Exception as e:
        logger.error(f"Failed to save figure to {output_path}: {e}")
    finally:
        plt.close(fig)
    return output_path


def plot_workflow_diagram(output_path: str) -> str:
    """
    Figure 1: Creates a flowchart-style diagram of the computational workflow.

    Args:
        output_path: Path to save the generated figure.

    Returns:
        Path where the figure was saved.
    """
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.axis('off')
    
    stages = [
        "Target Selection", "Structure Selection", "Binding Site ID",
        "Library Assembly", "Ligand Prep", "Protein Prep",
        "Docking", "Validation", "Screening",
        "Filtering", "Interaction Analysis", "Robustness", "Candidate Selection"
    ]
    
    box_w = 0.25
    box_h = 0.08
    y_starts = np.linspace(0.9, 0.1, len(stages))
    
    for i, stage in enumerate(stages):
        y = y_starts[i]
        rect = patches.Rectangle((0.375, y - box_h/2), box_w, box_h, 
                               linewidth=1, edgecolor='black', facecolor='#ECF0F1')
        ax.add_patch(rect)
        ax.text(0.5, y, stage, ha='center', va='center', fontsize=10, weight='bold')
        
        if i < len(stages) - 1:
            ax.arrow(0.5, y - box_h/2, 0, -(y_starts[i] - y_starts[i+1] - box_h),
                     head_width=0.02, head_length=0.02, fc='k', ec='k', 
                     length_includes_head=True)
            
    return _save_and_close(fig, output_path)


def plot_binding_pocket(receptor_pdb: str, ligand_pdb: str, center: Tuple[float, float, float],
                       box_size: Tuple[float, float, float], residue_list: List[str],
                       output_path: str) -> str:
    """
    Figure 2: Creates a schematic 2D representation of the binding pocket.

    Args:
        receptor_pdb: Path to receptor PDB file (unused in schematic).
        ligand_pdb: Path to ligand PDB file (unused in schematic).
        center: Coordinates of the docking box center.
        box_size: Dimensions of the docking box.
        residue_list: List of key residues.
        output_path: Path to save the generated figure.

    Returns:
        Path where the figure was saved.
    """
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.axis('equal')
    ax.axis('off')
    
    # Schematic central ligand
    circle = patches.Circle((0, 0), radius=2, edgecolor='black', facecolor='#E74C3C', alpha=0.6)
    ax.add_patch(circle)
    ax.text(0, 0, "Ligand Area\n" + f"Center: {center}", ha='center', va='center')
    
    # Schematic box
    rect = patches.Rectangle((-box_size[0]/2, -box_size[1]/2), box_size[0], box_size[1],
                             linewidth=2, edgecolor='blue', facecolor='none', linestyle='--')
    ax.add_patch(rect)
    
    # Plot residues in a circle around
    angles = np.linspace(0, 2 * np.pi, len(residue_list), endpoint=False)
    radius = max(box_size) / 2 + 2
    
    # Subsite colors S1=blue, S2=red, S2-ext=green, catalytic=gold
    # Simplified assigning just randomly for schematic purposes
    color_map = ['blue', 'red', 'green', 'gold']
    
    for i, (angle, res) in enumerate(zip(angles, residue_list)):
        x = radius * np.cos(angle)
        y = radius * np.sin(angle)
        color = color_map[i % len(color_map)]
        ax.plot(x, y, 'o', markersize=10, color=color)
        ax.text(x * 1.2, y * 1.2, res, ha='center', va='center', fontsize=9)

    ax.set_title("Schematic Binding Pocket", fontsize=14, weight='bold')
    return _save_and_close(fig, output_path)


def plot_score_distribution(docking_results_df: pd.DataFrame, output_path: str, score_col: str = 'best_score') -> str:
    """
    Figure 3: Histogram + KDE of all docking scores.

    Args:
        docking_results_df: DataFrame with docking results.
        output_path: Path to save the figure.
        score_col: Column name containing docking scores.

    Returns:
        Path where the figure was saved.
    """
    if score_col not in docking_results_df.columns:
        logger.error(f"Column '{score_col}' not found in dataframe.")
        return output_path
        
    fig, ax = plt.subplots(figsize=(8, 6))
    scores = docking_results_df[score_col].dropna()
    
    sns.histplot(scores, kde=True, ax=ax, color=COLORS['screening'], bins=30)
    
    mean_score = scores.mean()
    median_score = scores.median()
    
    ax.axvline(mean_score, color='red', linestyle='--', label=f'Mean: {mean_score:.2f}')
    ax.axvline(median_score, color='green', linestyle='-', label=f'Median: {median_score:.2f}')
    
    stats_text = (f"Count: {len(scores)}\n"
                  f"Min: {scores.min():.2f}\n"
                  f"Max: {scores.max():.2f}\n"
                  f"Std: {scores.std():.2f}")
    
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax.text(0.05, 0.95, stats_text, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', bbox=props)
            
    ax.set_title(f"Docking Score Distribution ({score_col})", fontsize=14)
    ax.set_xlabel("Docking Score (kcal/mol)")
    ax.set_ylabel("Frequency")
    ax.legend()
    
    return _save_and_close(fig, output_path)


def plot_controls_vs_library(docking_results_df: pd.DataFrame, positive_control_ids: List[str],
                             output_path: str, score_col: str = 'best_score') -> str:
    """
    Figure 4: Overlaid histograms comparing controls vs screening library.

    Args:
        docking_results_df: DataFrame with docking results.
        positive_control_ids: List of identifiers for positive controls.
        output_path: Path to save the figure.
        score_col: Column containing docking scores.

    Returns:
        Path where the figure was saved.
    """
    if score_col not in docking_results_df.columns:
        logger.error(f"Column '{score_col}' not found.")
        return output_path
        
    fig, ax = plt.subplots(figsize=(10, 6))
    
    controls = docking_results_df[docking_results_df.index.isin(positive_control_ids) | 
                                  docking_results_df.get('id', pd.Series()).isin(positive_control_ids)]
    
    library = docking_results_df[~docking_results_df.index.isin(positive_control_ids) & 
                                 ~docking_results_df.get('id', pd.Series()).isin(positive_control_ids)]
                                 
    c_scores = controls[score_col].dropna()
    l_scores = library[score_col].dropna()
    
    sns.kdeplot(l_scores, fill=True, ax=ax, color=COLORS['screening'], label='Screening Library')
    
    if len(c_scores) > 0 and len(l_scores) > 0:
        stat, pval = mannwhitneyu(c_scores, l_scores, alternative='two-sided')
        ax.set_title(f"Controls vs Library Scores\nMann-Whitney p-value: {pval:.2e}", fontsize=14)
    else:
        ax.set_title("Controls vs Library Scores", fontsize=14)

    # Plot individual controls
    if not c_scores.empty:
        # jitter y-axis for scatter
        y_pos = np.random.uniform(0.01, ax.get_ylim()[1]*0.1, size=len(c_scores))
        ax.scatter(c_scores, y_pos, color=COLORS['positive_control'], zorder=5, label='Controls', marker='X', s=100)
        
        for idx, (idx_name, row) in enumerate(controls.iterrows()):
            ax.annotate(str(idx_name)[:10], (c_scores.iloc[idx], y_pos[idx]), textcoords="offset points", xytext=(0,10), ha='center', fontsize=8)

    ax.set_xlabel("Docking Score (kcal/mol)")
    ax.set_ylabel("Density")
    ax.legend()
    
    return _save_and_close(fig, output_path)


def plot_top_candidate_scores(docking_results_df: pd.DataFrame, top_n: int = 20,
                              output_path: Optional[str] = None) -> str:
    """
    Figure 5: Horizontal bar chart of top N compounds by docking score.

    Args:
        docking_results_df: DataFrame with docking results.
        top_n: Number of top compounds to plot.
        output_path: Path to save the figure.

    Returns:
        Path where the figure was saved.
    """
    if 'best_score' not in docking_results_df.columns:
        return output_path if output_path else ""
        
    df_sorted = docking_results_df.sort_values(by='best_score', ascending=True).head(top_n)
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    y_pos = np.arange(len(df_sorted))
    scores = df_sorted['best_score'].values
    labels = df_sorted.index.astype(str)
    
    colors = [COLORS['candidate']] * len(df_sorted)
    
    bars = ax.barh(y_pos, scores, color=colors)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    
    for bar in bars:
        width = bar.get_width()
        ax.text(width + 0.1 if width > 0 else width - 0.1, bar.get_y() + bar.get_height()/2,
                f'{width:.2f}', ha='left' if width > 0 else 'right', va='center')
                
    ax.set_xlabel("Docking Score (kcal/mol)")
    ax.set_title(f"Top {top_n} Candidates by Docking Score", fontsize=14)
    
    if output_path is None:
        output_path = "top_candidates.png"
    return _save_and_close(fig, output_path)


def plot_interaction_heatmap(interactions_df: pd.DataFrame, output_path: str,
                             residue_cols: Optional[List[str]] = None) -> str:
    """
    Figure 6: Heatmap showing residues each top compound interacts with.

    Args:
        interactions_df: DataFrame where rows are compounds and cols are residues.
        output_path: Path to save the figure.
        residue_cols: List of column names representing residues.

    Returns:
        Path where the figure was saved.
    """
    if residue_cols:
        data = interactions_df[residue_cols].copy()
    else:
        # Try to infer residue columns by dropping non-numeric or known non-residue cols
        cols_to_drop = [c for c in interactions_df.columns if c.lower() in ('id', 'name', 'smiles', 'score')]
        data = interactions_df.drop(columns=cols_to_drop, errors='ignore').select_dtypes(include=[np.number, bool])
        
    if data.empty:
        logger.error("Interaction dataframe is empty or missing residue columns.")
        return output_path

    # Convert to numeric if boolean
    data = data.astype(float)
    
    fig, ax = plt.subplots(figsize=(12, 8))
    sns.heatmap(data, cmap="Blues", cbar_kws={'label': 'Interaction Presence/Strength'}, ax=ax)
    
    ax.set_title("Protein-Ligand Interaction Fingerprints", fontsize=14)
    ax.set_xlabel("Residue")
    ax.set_ylabel("Compound")
    
    return _save_and_close(fig, output_path)


def plot_chemical_properties(descriptors_df: pd.DataFrame, output_path: str) -> str:
    """
    Figure 7: Multi-panel figure showing distributions of chemical properties.

    Args:
        descriptors_df: DataFrame with property columns.
        output_path: Path to save the figure.

    Returns:
        Path where the figure was saved.
    """
    props = {
        'MW': ('Molecular Weight', 500),
        'LogP': ('LogP', 5),
        'HBD': ('H-Bond Donors', 5),
        'HBA': ('H-Bond Acceptors', 10),
        'RotBonds': ('Rotatable Bonds', 10),
        'TPSA': ('TPSA', 140)
    }
    
    available_props = [p for p in props.keys() if p in descriptors_df.columns]
    if not available_props:
        logger.error("No valid descriptor columns found in dataframe.")
        return output_path
        
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()
    
    for i, prop in enumerate(available_props):
        if i >= 6: break
        ax = axes[i]
        title, cutoff = props[prop]
        
        sns.histplot(descriptors_df[prop].dropna(), ax=ax, color=COLORS['screening'], kde=True)
        ax.axvline(cutoff, color='red', linestyle='--', label=f'Lipinski Cutoff ({cutoff})')
        ax.set_title(title)
        ax.set_xlabel(prop)
        ax.legend()
        
    for i in range(len(available_props), 6):
        axes[i].axis('off')
        
    fig.suptitle("Chemical Property Distributions (Ro5)", fontsize=16)
    return _save_and_close(fig, output_path)


def plot_robustness_scores(robustness_df: pd.DataFrame, output_path: str) -> str:
    """
    Figure 8: Multi-panel showing score variation across conditions.

    Args:
        robustness_df: DataFrame with robustness results.
        output_path: Path to save the figure.

    Returns:
        Path where the figure was saved.
    """
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    
    # Placeholder for actual robustness plotting logic
    # Assume dataframe has columns: 'compound_id', 'condition_type', 'condition_value', 'score'
    if all(c in robustness_df.columns for c in ['compound_id', 'condition_type', 'condition_value', 'score']):
        
        for i, cond in enumerate(['box_size', 'exhaustiveness', 'seed']):
            if i >= 3: break
            ax = axes[i]
            df_cond = robustness_df[robustness_df['condition_type'] == cond]
            if df_cond.empty: continue
            
            if cond == 'seed':
                sns.stripplot(data=df_cond, x='compound_id', y='score', ax=ax, jitter=True)
                ax.set_title("Effect of Random Seed")
            else:
                sns.lineplot(data=df_cond, x='condition_value', y='score', hue='compound_id', ax=ax, marker='o')
                ax.set_title(f"Effect of {cond.capitalize()}")
                ax.legend().set_visible(False)
                
            ax.set_ylabel("Docking Score")
    else:
        # Fallback empty plots
        for ax in axes: ax.axis('off')
        axes[1].text(0.5, 0.5, "Robustness Data Format Unknown", ha='center', va='center')
        
    fig.suptitle("Docking Robustness Analysis", fontsize=16)
    return _save_and_close(fig, output_path)


def plot_rank_correlation(rank_corr_df: pd.DataFrame, output_path: str) -> str:
    """
    Figure 9: Heatmap of Spearman rank correlations between docking conditions.

    Args:
        rank_corr_df: Correlation matrix dataframe.
        output_path: Path to save the figure.

    Returns:
        Path where the figure was saved.
    """
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.heatmap(rank_corr_df, annot=True, cmap="coolwarm", vmin=-1, vmax=1, center=0, ax=ax, square=True)
    ax.set_title("Spearman Rank Correlation Across Conditions", fontsize=14)
    
    return _save_and_close(fig, output_path)


def plot_enrichment_curve(docking_results_df: pd.DataFrame, positive_control_ids: List[str],
                          output_path: str) -> str:
    """
    Bonus: ROC-like enrichment curve.

    Args:
        docking_results_df: DataFrame with docking results.
        positive_control_ids: List of active identifiers.
        output_path: Path to save the figure.

    Returns:
        Path where the figure was saved.
    """
    if 'best_score' not in docking_results_df.columns:
        return output_path
        
    df_sorted = docking_results_df.sort_values(by='best_score', ascending=True)
    
    is_active = df_sorted.index.isin(positive_control_ids) | df_sorted.get('id', pd.Series()).isin(positive_control_ids)
    
    n_actives = sum(is_active)
    n_total = len(df_sorted)
    
    if n_actives == 0:
        logger.warning("No active compounds found for enrichment curve.")
        return output_path
        
    x = np.arange(1, n_total + 1) / n_total
    y = np.cumsum(is_active) / n_actives
    
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.plot(x, y, color=COLORS['positive_control'], lw=2, label='Screening Results')
    ax.plot([0, 1], [0, 1], color='gray', linestyle='--', label='Random')
    
    ax.set_xlabel('Fraction of Library Screened')
    ax.set_ylabel('Fraction of Actives Recovered')
    ax.set_title('Enrichment Curve', fontsize=14)
    ax.legend()
    ax.grid(True)
    
    return _save_and_close(fig, output_path)


def generate_all_figures(results_dir: str, figures_dir: str) -> None:
    """
    Master function that loads all result CSVs and generates all figures.

    Args:
        results_dir: Directory containing input CSV files.
        figures_dir: Directory to save generated figures.
    """
    os.makedirs(figures_dir, exist_ok=True)
    
    logger.info(f"Generating all figures from {results_dir} into {figures_dir}")
    
    plot_workflow_diagram(os.path.join(figures_dir, "Figure1_Workflow.png"))
    plot_binding_pocket("rec.pdb", "lig.pdb", (0, 0, 0), (20, 20, 20),
                        ["HIS100", "SER200", "ASP300"],
                        os.path.join(figures_dir, "Figure2_BindingPocket.png"))
                        
    results_path = os.path.join(results_dir, "docking_results.csv")
    if os.path.exists(results_path):
        try:
            df = pd.read_csv(results_path)
            if 'id' in df.columns:
                df.set_index('id', inplace=True)
                
            plot_score_distribution(df, os.path.join(figures_dir, "Figure3_ScoreDistribution.png"))
            
            controls = ['sitagliptin', 'vildagliptin'] # Example IDs
            plot_controls_vs_library(df, controls, os.path.join(figures_dir, "Figure4_ControlsVsLibrary.png"))
            plot_top_candidate_scores(df, top_n=20, output_path=os.path.join(figures_dir, "Figure5_TopCandidates.png"))
            plot_enrichment_curve(df, controls, os.path.join(figures_dir, "Figure10_EnrichmentCurve.png"))
        except Exception as e:
            logger.error(f"Error processing {results_path}: {e}")
            
    interactions_path = os.path.join(results_dir, "interactions.csv")
    if os.path.exists(interactions_path):
        try:
            int_df = pd.read_csv(interactions_path, index_col=0)
            plot_interaction_heatmap(int_df, os.path.join(figures_dir, "Figure6_Interactions.png"))
        except Exception as e:
            logger.error(f"Error processing {interactions_path}: {e}")
            
    descriptors_path = os.path.join(results_dir, "descriptors.csv")
    if os.path.exists(descriptors_path):
        try:
            desc_df = pd.read_csv(descriptors_path)
            plot_chemical_properties(desc_df, os.path.join(figures_dir, "Figure7_ChemicalProperties.png"))
        except Exception as e:
            logger.error(f"Error processing {descriptors_path}: {e}")
            
    robustness_path = os.path.join(results_dir, "robustness.csv")
    if os.path.exists(robustness_path):
        try:
            rob_df = pd.read_csv(robustness_path)
            plot_robustness_scores(rob_df, os.path.join(figures_dir, "Figure8_Robustness.png"))
        except Exception as e:
            logger.error(f"Error processing {robustness_path}: {e}")
            
    corr_path = os.path.join(results_dir, "rank_correlation.csv")
    if os.path.exists(corr_path):
        try:
            corr_df = pd.read_csv(corr_path, index_col=0)
            plot_rank_correlation(corr_df, os.path.join(figures_dir, "Figure9_RankCorrelation.png"))
        except Exception as e:
            logger.error(f"Error processing {corr_path}: {e}")
            
    logger.info("Figure generation complete.")
