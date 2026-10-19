#!/usr/bin/env python3
# Plots PCA scatters, probability boxplots, ancestry counts and the reference confusion matrix.
# Usage: python3 06_plot_results.py --pcs X --assignments Y --ref-predictions Z --ref-labels L --out-prefix P
import argparse
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path


POP_COLORS = {
    'AFR': '#E74C3C',
    'AMR': '#F39C12',
    'EAS': '#27AE60',
    'EUR': '#2980B9',
    'SAS': '#8E44AD',
    'Unassigned': '#95A5A6',
    'reference_1KGP3': '#BDC3C7',
}


def plot_pca_scatter(all_pcs: pd.DataFrame, study_assignments: pd.DataFrame,
                     ref_labels: pd.DataFrame, out_prefix: str, pc_x: int = 1, pc_y: int = 2):
    fig, axes = plt.subplots(1, 2, figsize=(18, 7))

    pcx, pcy = f'PC{pc_x}', f'PC{pc_y}'

    for ax_idx, (ax, title) in enumerate(zip(axes, ['Reference (1KGP3)', 'Study samples'])):
        ref = pd.merge(all_pcs[all_pcs['source'] == 'reference_1KGP3'],
                       ref_labels[['IID','SuperPop']], on='IID', how='left')
        ref['SuperPop'] = ref['SuperPop'].fillna('Unknown')
        for pop in sorted(ref['SuperPop'].unique()):
            sub = ref[ref['SuperPop'] == pop]
            ax.scatter(sub[pcx], sub[pcy],
                       c=POP_COLORS.get(pop, '#ccc'),
                       s=8, alpha=0.4, label=f'1KGP3 {pop}',
                       zorder=2)

        if ax_idx == 1:
            study = pd.merge(all_pcs[all_pcs['source'] == 'study'],
                             study_assignments[['IID','AncestryGEL']], on='IID', how='left')
            for pop in sorted(study['AncestryGEL'].dropna().unique()):
                sub = study[study['AncestryGEL'] == pop]
                ax.scatter(sub[pcx], sub[pcy],
                           c=POP_COLORS.get(pop, '#ccc'),
                           s=20, alpha=0.8, edgecolors='black', linewidths=0.3,
                           label=f'Study {pop}', marker='D', zorder=3)

        ax.set_xlabel(f'PC{pc_x}', fontsize=12)
        ax.set_ylabel(f'PC{pc_y}', fontsize=12)
        ax.set_title(title, fontsize=13)
        ax.legend(fontsize=8, markerscale=1.5, bbox_to_anchor=(1.01, 1), loc='upper left')

    plt.suptitle(f'Ancestry inference: PC{pc_x} vs PC{pc_y}', fontsize=14, y=1.01)
    plt.tight_layout()
    out = f"{out_prefix}_pca_PC{pc_x}_PC{pc_y}.png"
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {out}")


def plot_probability_distributions(study_assignments: pd.DataFrame, out_prefix: str):
    pops = sorted(study_assignments['AncestryGEL'].unique())
    fig, ax = plt.subplots(figsize=(10, 5))
    data = [study_assignments[study_assignments['AncestryGEL'] == p]['MaxProb'].values
            for p in pops]
    bp = ax.boxplot(data, patch_artist=True, notch=False)
    for patch, pop in zip(bp['boxes'], pops):
        patch.set_facecolor(POP_COLORS.get(pop, '#ccc'))
    ax.set_xticklabels(pops, fontsize=11)
    ax.set_ylabel('Max RF Probability', fontsize=12)
    ax.axhline(y=0.8, color='black', linestyle='--', linewidth=1, label='GEL threshold (0.8)')
    ax.legend(fontsize=10)
    ax.set_title('Classification Probability by Assigned Ancestry', fontsize=13)
    plt.tight_layout()
    out = f"{out_prefix}_probability_distribution.png"
    plt.savefig(out, dpi=150)
    plt.close()
    print(f"Saved: {out}")


def plot_ancestry_barplot(study_assignments: pd.DataFrame, out_prefix: str):
    counts = study_assignments['AncestryGEL'].value_counts()
    fig, ax = plt.subplots(figsize=(9, 5))
    colors = [POP_COLORS.get(p, '#ccc') for p in counts.index]
    bars = ax.bar(counts.index, counts.values, color=colors, edgecolor='black', linewidth=0.5)
    for bar, val in zip(bars, counts.values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1,
                str(val), ha='center', va='bottom', fontsize=10)
    ax.set_xlabel('Ancestry Group', fontsize=12)
    ax.set_ylabel('Number of Samples', fontsize=12)
    ax.set_title('Ancestry Assignment, Study Cohort', fontsize=13)
    total = counts.sum()
    ax.text(0.98, 0.98, f'Total = {total}', transform=ax.transAxes,
            ha='right', va='top', fontsize=10, color='#555')
    plt.tight_layout()
    out = f"{out_prefix}_ancestry_counts.png"
    plt.savefig(out, dpi=150)
    plt.close()
    print(f"Saved: {out}")


def plot_confusion_matrix(ref_predictions: pd.DataFrame, out_prefix: str):
    if 'TrueLabel' not in ref_predictions.columns:
        return
    from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay

    pops = sorted(ref_predictions['TrueLabel'].unique())
    cm = confusion_matrix(ref_predictions['TrueLabel'],
                          ref_predictions['PredLabel'], labels=pops)
    fig, ax = plt.subplots(figsize=(8, 7))
    disp = ConfusionMatrixDisplay(confusion_matrix=cm, display_labels=pops)
    disp.plot(ax=ax, colorbar=True, cmap='Blues', values_format='d')
    ax.set_title('Reference Sample Classification (QC)\nRandom Forest', fontsize=12)
    plt.tight_layout()
    out = f"{out_prefix}_confusion_matrix.png"
    plt.savefig(out, dpi=150)
    plt.close()
    print(f"Saved: {out}")


def main():
    parser = argparse.ArgumentParser(description='Plot ancestry inference results')
    parser.add_argument('--pcs',              required=True, help='All PCs TSV file')
    parser.add_argument('--assignments',      required=True, help='Study ancestry assignments TSV')
    parser.add_argument('--ref-predictions',  required=True, help='Reference predictions TSV (QC)')
    parser.add_argument('--ref-labels',       required=True, help='1KGP3 super-population labels TSV')
    parser.add_argument('--out-prefix',       required=True, help='Output prefix for plots')
    args = parser.parse_args()

    Path(args.out_prefix).parent.mkdir(parents=True, exist_ok=True)

    all_pcs = pd.read_csv(args.pcs, sep='\t')
    study_assignments = pd.read_csv(args.assignments, sep='\t')
    ref_predictions   = pd.read_csv(args.ref_predictions, sep='\t')
    ref_labels        = pd.read_csv(args.ref_labels, sep='\t')

    # Numeric-looking IIDs are read as int64 in one table and str in another,
    # which breaks the merges.
    for df in (all_pcs, study_assignments, ref_predictions):
        if 'IID' in df.columns:
            df['IID'] = df['IID'].astype(str)

    col_map = {}
    for c in ref_labels.columns:
        cl = c.lower()
        if cl in ('iid', 'sample', 'sample_id'): col_map[c] = 'IID'
        if 'super' in cl and 'pop' in cl:         col_map[c] = 'SuperPop'
    ref_labels.rename(columns=col_map, inplace=True)

    plot_pca_scatter(all_pcs, study_assignments, ref_labels, args.out_prefix, 1, 2)
    plot_pca_scatter(all_pcs, study_assignments, ref_labels, args.out_prefix, 1, 3)
    plot_pca_scatter(all_pcs, study_assignments, ref_labels, args.out_prefix, 2, 3)
    plot_probability_distributions(study_assignments, args.out_prefix)
    plot_ancestry_barplot(study_assignments, args.out_prefix)
    plot_confusion_matrix(ref_predictions, args.out_prefix)


if __name__ == '__main__':
    main()
