import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
from scipy.stats import pearsonr, linregress
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from constants import GENETIC_CODE_TABLES, CODON_GRID_ORDER, ALL_CODONS_SORTED, AA_CODON_MAPS

def generate_rscu_heatmap_and_table(all_data, output_folder, genetic_code_id, status_queue, palette='viridis'):
    base_filename = list(all_data.keys())[0]
    data = all_data[base_filename]
    codon_counts = data['counts']
    rscu_values = data['rscu']
    
    print(f"\n=== RSCU HEATMAP ANALYSIS ===")
    print(f"File: {base_filename}")
    
    status_queue.put(("progress", 20))
    
    print("  Saving codon count table...")
    df_counts = pd.DataFrame(codon_counts.items(), columns=['Codon', 'Count']).sort_values(by='Count', ascending=False)
    csv_path = os.path.join(output_folder, f"codon_counts_{base_filename}.csv")
    df_counts.to_csv(csv_path, index=False, sep=';')
    print(f"  ✅ Count table saved in: {csv_path}")
    
    status_queue.put(("progress", 50))
    
    print(f"  📊 Statistics:")
    print(f"     ENC (Effective Number of Codons): {data['enc']:.2f}")
    print(f"     GC3 (GC at third position): {data['gc3']:.2f}%")
    print(f"     CAI (Codon Adaptation Index): {data['cai']:.3f}")
    
    status_queue.put(("progress", 70))
    print("  Preparing data for charts...")
    rscu_data_grid = []
    codon_labels_grid = []
    codon_aa_map = GENETIC_CODE_TABLES.get(genetic_code_id, GENETIC_CODE_TABLES[1])
    
    for row_codons in CODON_GRID_ORDER:
        value_row, label_row = [], []
        for codon in row_codons:
            value = rscu_values.get(codon, 0.0)
            value_row.append(value)
            aa = codon_aa_map.get(codon, '?')
            label_row.append(f"{codon} ({aa})\n{value:.2f}")
        rscu_data_grid.append(value_row)
        codon_labels_grid.append(label_row)
        
    df_rscu = pd.DataFrame(rscu_data_grid)
    
    df_rscu_plot = pd.DataFrame(list(rscu_values.items()), columns=['Codon', 'RSCU'])
    df_rscu_plot['AminoAcid'] = df_rscu_plot['Codon'].map(codon_aa_map)
    df_rscu_plot = df_rscu_plot.sort_values(by=['AminoAcid', 'Codon'])
    csv_rscu_path = os.path.join(output_folder, f"rscu_values_{base_filename}.csv")
    df_rscu_plot.to_csv(csv_rscu_path, index=False, sep=';')
    print(f"  ✅ RSCU values table saved in: {csv_rscu_path}")

    print("  Generating Chart 1: RSCU Heatmap...")
    try:
        plt.figure(figsize=(20, 8))
        ax = sns.heatmap(df_rscu, annot=codon_labels_grid, fmt="", cmap=palette, linewidths=0.5, cbar_kws={'label': 'RSCU Value'})
        ax.set_yticklabels(['T', 'C', 'A', 'G'], rotation=0)
        ax.set_xticklabels([f"Pos {i+1}" for i in range(16)], rotation=0)
        ax.set_title(f"1. Relative Synonymous Codon Usage (RSCU) - {base_filename}\nENC={data['enc']:.2f}, GC3={data['gc3']:.2f}%, CAI={data['cai']:.3f}", fontsize=16)
        
        output_file1 = os.path.join(output_folder, f"rscu_1_heatmap_{base_filename}.png")
        plt.savefig(output_file1, dpi=150, bbox_inches="tight")
        plt.close()
        status_queue.put(("image_ready", (output_file1, f"1. RSCU Heatmap")))
    except Exception as e:
        print(f"\n❌ ERROR GENERATING CHART 1: {e}")

    print("  Generating Chart 2: RSCU Polar Plot...")
    try:
        codons_polar = df_rscu_plot['Codon'].values
        rscu_polar = df_rscu_plot['RSCU'].values
        angles = np.linspace(0, 2 * np.pi, len(codons_polar), endpoint=False)
        
        rscu_polar_closed = np.concatenate((rscu_polar, [rscu_polar[0]]))
        angles_closed = np.concatenate((angles, [angles[0]]))
        
        fig, ax = plt.subplots(figsize=(10, 10), subplot_kw={'projection': 'polar'})
        try:
            plot_color = sns.color_palette(palette)[0]
        except:
            plot_color = 'teal'
            
        ax.plot(angles_closed, rscu_polar_closed, linewidth=2, color=plot_color)
        ax.fill(angles_closed, rscu_polar_closed, alpha=0.3, color=plot_color)
        
        expected_circle = np.linspace(0, 2 * np.pi, 100)
        ax.plot(expected_circle, [1.0] * 100, color='red', linestyle='--', label='RSCU=1.0 (No Bias)')
        
        ax.set_xticks(angles)
        ax.set_xticklabels(codons_polar, fontsize=8)
        plt.title(f"2. Polar Plot of RSCU Distribution - {base_filename}", y=1.08, fontsize=16)
        plt.legend(loc='upper right', bbox_to_anchor=(1.3, 1.1))
        
        output_file2 = os.path.join(output_folder, f"rscu_2_polar_{base_filename}.png")
        plt.savefig(output_file2, dpi=150, bbox_inches="tight")
        plt.close()
        status_queue.put(("image_ready", (output_file2, f"2. RSCU Polar Plot")))
    except Exception as e:
        print(f"\n❌ ERROR GENERATING CHART 2: {e}")

    print("  Generating Chart 3: RSCU Line Plot...")
    try:
        aa_codon_map = AA_CODON_MAPS.get(genetic_code_id, AA_CODON_MAPS[1])
        ordered_codons = []
        ordered_aas = []
        for aa in sorted(aa_codon_map.keys()):
            if aa == '*': continue
            for codon in sorted(aa_codon_map[aa]):
                ordered_codons.append(codon)
                ordered_aas.append(aa)

        plt.figure(figsize=(26, 8))
        
        bg_colors = ['#ffffff', '#f0f4f8']
        current_aa = ordered_aas[0]
        start_idx = 0
        color_idx = 0
        
        for i, aa in enumerate(ordered_aas):
            if aa != current_aa:
                plt.axvspan(start_idx - 0.5, i - 0.5, facecolor=bg_colors[color_idx % 2], alpha=1.0, zorder=0)
                current_aa = aa
                start_idx = i
                color_idx += 1
        plt.axvspan(start_idx - 0.5, len(ordered_aas) - 0.5, facecolor=bg_colors[color_idx % 2], alpha=1.0, zorder=0)

        y_vals = [rscu_values.get(c, 0.0) for c in ordered_codons]
        try:
            line_color = sns.color_palette(palette)[0]
        except:
            line_color = 'teal'
            
        plt.plot(range(len(ordered_codons)), y_vals, label=base_filename, linewidth=2.5, marker='o', markersize=5, zorder=3, color=line_color)

        plt.axhline(1.5, color='darkred', linestyle='--', alpha=0.8, label="Optimal (>1.5)", zorder=2)
        plt.axhline(0.5, color='darkred', linestyle='--', alpha=0.8, label="Rare (<0.5)", zorder=2)
        
        labels = [f"{c}\n({a})" for c, a in zip(ordered_codons, ordered_aas)]
        plt.xticks(range(len(ordered_codons)), labels, rotation=90, fontsize=11)
        
        plt.xlabel("Codon / Amino Acid", fontsize=13, fontweight='bold')
        plt.ylabel("RSCU Value", fontsize=13, fontweight='bold')
        plt.title(f"3. RSCU Distribution Profile Across All Codons - {base_filename}", fontsize=18, fontweight='bold')
        
        plt.grid(axis='y', linestyle=':', alpha=0.6, zorder=1)
        plt.legend(bbox_to_anchor=(1.01, 1), loc='upper left')
        plt.margins(x=0.01)
        plt.tight_layout()

        output_file3 = os.path.join(output_folder, f"rscu_3_lineplot_{base_filename}.png")
        plt.savefig(output_file3, dpi=150, bbox_inches="tight")
        plt.close()
        status_queue.put(("image_ready", (output_file3, f"3. RSCU Line Plot")))
    except Exception as e:
        print(f"\n❌ ERROR GENERATING CHART 3: {e}")

def comparative_rscu_analysis(all_data, output_folder, status_queue, genetic_code_id=1, palette='viridis'):
    """Unified comparative RSCU module (analyses 4 + 6).

    Charts produced:
      1.  Clustermap        — hierarchical clustering of all species × 64 codons
      2A. PCA Scores        — species in PC1 × PC2 space (samples space)
      2B. PCA Loadings      — codon loading arrows in PC1 × PC2 space (biplot)
      3.  Line Plot         — per-codon RSCU profile for all species, grouped by AA
      4.  Box Plot by AA    — RSCU distribution per amino acid across all species
      5.  Variance Barplot  — top-30 codons by inter-species variance
      6.  Heatmap Top-30    — RSCU values for the 30 most divergent codons
    """
    print(f"\n=== COMPARATIVE RSCU ANALYSIS (unified) ===")
    all_rscu_data = {species: data['rscu'] for species, data in all_data.items()}

    # ── shared data ───────────────────────────────────────────────────────────
    status_queue.put(("progress", 40))
    df_rscu_matrix = pd.DataFrame.from_dict(
        all_rscu_data, orient='index', columns=ALL_CODONS_SORTED
    ).fillna(0.0)
    df_rscu_matrix.index.name = 'Species'

    csv_matrix_path = os.path.join(output_folder, 'comparative_rscu_matrix.csv')
    df_rscu_matrix.to_csv(csv_matrix_path, sep=';', decimal='.')
    print(f"  ✅ RSCU matrix saved: {csv_matrix_path}")

    try:
        colors = sns.color_palette(palette, max(len(all_data), 3))
    except Exception:
        colors = ['teal', 'coral', 'steelblue']

    # ── Chart 1: Clustermap ───────────────────────────────────────────────────
    print("  Generating Chart 1: Clustermap...")
    status_queue.put(("progress", 50))
    try:
        g = sns.clustermap(
            df_rscu_matrix,
            metric="euclidean", method="average",
            cmap=palette, annot=False, linewidths=0.3,
            figsize=(max(16, len(ALL_CODONS_SORTED) * 0.22), max(8, len(all_data) * 0.6)),
            cbar_kws={'label': 'RSCU'}
        )
        g.fig.suptitle("1. RSCU Hierarchical Clustermap", y=1.02, fontsize=16)
        plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize=9)
        clustermap_path = os.path.join(output_folder, 'comparative_1_clustermap.png')
        g.savefig(clustermap_path, dpi=150, bbox_inches="tight")
        plt.close(g.fig)
        status_queue.put(("image_ready", (clustermap_path, "1. RSCU Clustermap")))
    except Exception as e:
        print(f"  ❌ Chart 1 error: {e}")

    # ── Chart 2A: PCA Scores (species in samples space) ──────────────────────
    # ── Chart 2B: PCA Loadings biplot (codon arrows from origin) ─────────────
    print("  Generating Chart 2: PCA (scores + loadings, separate figures)...")
    status_queue.put(("progress", 60))
    try:
        X_scaled = StandardScaler().fit_transform(df_rscu_matrix.values)
        pca = PCA(n_components=2)
        scores = pca.fit_transform(X_scaled)
        loadings = pca.components_.T          # shape: (n_codons, 2)
        pc1_var, pc2_var = pca.explained_variance_ratio_ * 100

        score_pc1_max = max(np.max(np.abs(scores[:, 0])), 1e-9)
        score_pc2_max = max(np.max(np.abs(scores[:, 1])), 1e-9)

        # ---- 2A. Samples / scores plot ------------------------------------
        fig_a, ax_a = plt.subplots(figsize=(10, 8))
        sp_colors = sns.color_palette(palette, len(df_rscu_matrix))
        for i, species in enumerate(df_rscu_matrix.index):
            x, y = scores[i, 0], scores[i, 1]
            ax_a.scatter(x, y, s=180, color=sp_colors[i],
                         edgecolors='black', zorder=3)
            # Position-aware label placement so names near the edge don't get
            # clipped: anchor the label on the side facing the chart center.
            ha = 'right' if x > 0 else 'left'
            x_off = -8 if x > 0 else 8
            y_off = 6 if y >= 0 else -10
            va = 'bottom' if y >= 0 else 'top'
            ax_a.annotate(species, (x, y),
                          xytext=(x_off, y_off), textcoords='offset points',
                          fontsize=10, fontweight='bold',
                          ha=ha, va=va, clip_on=False, zorder=4)
        ax_a.axhline(0, color='grey', linestyle='--', linewidth=0.8)
        ax_a.axvline(0, color='grey', linestyle='--', linewidth=0.8)
        # Extra margin so even position-aware labels for long species names
        # at extreme scores still have room before bbox_inches='tight' kicks in.
        ax_a.set_xlim(-score_pc1_max * 1.35, score_pc1_max * 1.35)
        ax_a.set_ylim(-score_pc2_max * 1.35, score_pc2_max * 1.35)
        ax_a.set_xlabel(f'Principal Component 1 ({pc1_var:.2f}%)', fontsize=12)
        ax_a.set_ylabel(f'Principal Component 2 ({pc2_var:.2f}%)', fontsize=12)
        ax_a.set_title('2A. RSCU Principal Component Analysis (Samples Space)', fontsize=14)
        ax_a.grid(alpha=0.3)
        plt.tight_layout()
        path_2a = os.path.join(output_folder, 'comparative_2a_pca_samples.png')
        fig_a.savefig(path_2a, dpi=150, bbox_inches="tight")
        plt.close(fig_a)
        status_queue.put(("image_ready", (path_2a, "2A. PCA — Samples Space")))

        # ---- 2B. Loadings / codon-biplot ---------------------------------
        # Plot all 64 codons in loading space (PCA components_, not scaled).
        load_pc1_max = max(np.max(np.abs(loadings[:, 0])), 1e-9)
        load_pc2_max = max(np.max(np.abs(loadings[:, 1])), 1e-9)

        fig_b, ax_b = plt.subplots(figsize=(10, 8))
        for idx, codon in enumerate(ALL_CODONS_SORTED):
            xarr, yarr = loadings[idx, 0], loadings[idx, 1]
            ax_b.annotate("", xy=(xarr, yarr), xytext=(0, 0),
                          arrowprops=dict(arrowstyle='->', color='#3a3a82',
                                          lw=1.0, alpha=0.7))
            ax_b.text(xarr * 1.08, yarr * 1.08, codon,
                      fontsize=8, color='#3a3a82',
                      ha='center', va='center', clip_on=True)

        ax_b.axhline(0, color='grey', linestyle='--', linewidth=0.8)
        ax_b.axvline(0, color='grey', linestyle='--', linewidth=0.8)
        ax_b.set_xlim(-load_pc1_max * 1.25, load_pc1_max * 1.25)
        ax_b.set_ylim(-load_pc2_max * 1.25, load_pc2_max * 1.25)
        ax_b.set_xlabel(f'Component 1 ({pc1_var:.2f}%)', fontsize=12)
        ax_b.set_ylabel(f'Component 2 ({pc2_var:.2f}%)', fontsize=12)
        ax_b.set_title('2B. RSCU Principal Component Analysis (Codon Loadings / Biplot)', fontsize=14)
        ax_b.grid(alpha=0.3)
        plt.tight_layout()
        path_2b = os.path.join(output_folder, 'comparative_2b_pca_loadings.png')
        fig_b.savefig(path_2b, dpi=150, bbox_inches="tight")
        plt.close(fig_b)
        status_queue.put(("image_ready", (path_2b, "2B. PCA — Codon Loadings")))
    except Exception as e:
        print(f"  ❌ Chart 2 error: {e}")

    # ── Chart 3: Line Plot per codon, grouped by amino acid ───────────────────
    print("  Generating Chart 3: Comparative Line Plot...")
    status_queue.put(("progress", 70))
    try:
        # build codon order grouped by AA (using the selected genetic code)
        aa_codon_map = AA_CODON_MAPS.get(genetic_code_id, AA_CODON_MAPS[1])
        ordered_codons, ordered_aas = [], []
        for aa in sorted(aa_codon_map.keys()):
            if aa == '*':
                continue
            for codon in sorted(aa_codon_map[aa]):
                ordered_codons.append(codon)
                ordered_aas.append(aa)

        fig, ax = plt.subplots(figsize=(26, 8))

        # alternating AA background bands
        bg_colors = ['#ffffff', '#f0f4f8']
        current_aa, start_idx, col_idx = ordered_aas[0], 0, 0
        for i, aa in enumerate(ordered_aas):
            if aa != current_aa:
                ax.axvspan(start_idx - 0.5, i - 0.5,
                           facecolor=bg_colors[col_idx % 2], alpha=1.0, zorder=0)
                current_aa, start_idx, col_idx = aa, i, col_idx + 1
        ax.axvspan(start_idx - 0.5, len(ordered_aas) - 0.5,
                   facecolor=bg_colors[col_idx % 2], alpha=1.0, zorder=0)

        sp_colors = sns.color_palette(palette, len(df_rscu_matrix))
        for sp_name, col in zip(df_rscu_matrix.index, sp_colors):
            y_vals = [df_rscu_matrix.loc[sp_name, c] for c in ordered_codons]
            ax.plot(range(len(ordered_codons)), y_vals,
                    label=sp_name, color=col, linewidth=2, marker='o', markersize=4, zorder=3)

        ax.axhline(1.5, color='darkred', linestyle='--', alpha=0.7, label='Optimal (>1.5)', zorder=2)
        ax.axhline(0.5, color='darkred', linestyle=':', alpha=0.7, label='Rare (<0.5)', zorder=2)

        labels = [f"{c}\n({a})" for c, a in zip(ordered_codons, ordered_aas)]
        ax.set_xticks(range(len(ordered_codons)))
        ax.set_xticklabels(labels, rotation=90, fontsize=10)
        ax.set_xlabel("Codon (Amino Acid)", fontsize=13)
        ax.set_ylabel("RSCU", fontsize=13)
        ax.set_title("3. RSCU Comparative Profile — All Species × All Codons (grouped by AA)", fontsize=16)
        ax.legend(bbox_to_anchor=(1.01, 1), loc='upper left', fontsize=9)
        ax.grid(axis='y', linestyle=':', alpha=0.5, zorder=1)
        ax.margins(x=0.01)
        plt.tight_layout()
        line_path = os.path.join(output_folder, 'comparative_3_lineplot.png')
        plt.savefig(line_path, dpi=150, bbox_inches="tight")
        plt.close()
        status_queue.put(("image_ready", (line_path, "3. Comparative Line Plot")))
    except Exception as e:
        print(f"  ❌ Chart 3 error: {e}")

    # ── Chart 4: Box Plot by Amino Acid ───────────────────────────────────────
    print("  Generating Chart 4: Box Plot by Amino Acid...")
    status_queue.put(("progress", 80))
    try:
        codon_aa_map = GENETIC_CODE_TABLES.get(genetic_code_id, GENETIC_CODE_TABLES[1])
        df_long = (
            df_rscu_matrix.reset_index()
            .rename(columns={'Species': 'Species'})
            .melt(id_vars='Species', var_name='Codon', value_name='RSCU')
        )
        df_long['AminoAcid'] = df_long['Codon'].map(codon_aa_map)
        df_long = df_long[df_long['AminoAcid'].notna() & (df_long['AminoAcid'] != '*')]
        aa_order = sorted(df_long['AminoAcid'].unique())

        plt.figure(figsize=(20, 8))
        sns.boxplot(x='AminoAcid', y='RSCU', data=df_long,
                    order=aa_order, palette='Set3', showfliers=False)
        plt.axhline(1.0, color='red', linestyle='--', alpha=0.6, label='RSCU = 1.0 (uniform)')
        plt.xlabel('Amino Acid', fontsize=12)
        plt.ylabel('RSCU', fontsize=12)
        plt.title('4. RSCU Distribution by Amino Acid (All Species)', fontsize=16)
        plt.legend()
        plt.grid(axis='y', alpha=0.3)
        plt.tight_layout()
        box_path = os.path.join(output_folder, 'comparative_4_boxplot_by_aa.png')
        plt.savefig(box_path, dpi=150, bbox_inches="tight")
        plt.close()
        status_queue.put(("image_ready", (box_path, "4. Box Plot by Amino Acid")))
    except Exception as e:
        print(f"  ❌ Chart 4 error: {e}")

    # ── Chart 5: Variance barplot — top 30 most variable codons ──────────────
    print("  Generating Chart 5: Top 30 Most Variable Codons (variance barplot)...")
    status_queue.put(("progress", 90))
    codon_variance = None
    top30 = None
    try:
        codon_variance = df_rscu_matrix.var(axis=0).sort_values(ascending=False)
        top30 = codon_variance.head(30)

        fig5, ax5 = plt.subplots(figsize=(14, 6))
        sp_col = sns.color_palette(palette)[0]
        ax5.bar(range(len(top30)), top30.values, color=sp_col, edgecolor='black', alpha=0.85)
        ax5.set_xticks(range(len(top30)))
        ax5.set_xticklabels(top30.index, rotation=90, fontsize=10)
        ax5.set_xlabel('Codon', fontsize=12)
        ax5.set_ylabel('Variance across species', fontsize=12)
        ax5.set_title('5. Top 30 Most Variable Codons (inter-species variance)', fontsize=14)
        ax5.grid(axis='y', alpha=0.3)
        plt.tight_layout()
        path_5 = os.path.join(output_folder, 'comparative_5_variance.png')
        fig5.savefig(path_5, dpi=150, bbox_inches="tight")
        plt.close(fig5)
        status_queue.put(("image_ready", (path_5, "5. Top 30 Most Variable Codons")))
    except Exception as e:
        print(f"  ❌ Chart 5 error: {e}")

    # ── Chart 6: Heatmap — RSCU of top-30 divergent codons ───────────────────
    print("  Generating Chart 6: RSCU Heatmap — Top 30 Divergent Codons...")
    status_queue.put(("progress", 96))
    try:
        if top30 is None:
            codon_variance = df_rscu_matrix.var(axis=0).sort_values(ascending=False)
            top30 = codon_variance.head(30)
        df_top = df_rscu_matrix[top30.index]
        fig6, ax6 = plt.subplots(figsize=(16, max(4, len(df_rscu_matrix) * 0.8 + 2)))
        sns.heatmap(df_top, cmap=palette, annot=True, fmt=".2f",
                    linewidths=0.4, ax=ax6, cbar_kws={'label': 'RSCU'})
        ax6.set_title('6. RSCU Heatmap — Top 30 Divergent Codons', fontsize=14)
        ax6.set_xticklabels(ax6.get_xticklabels(), rotation=90, fontsize=9)
        plt.tight_layout()
        path_6 = os.path.join(output_folder, 'comparative_6_heatmap_top30.png')
        fig6.savefig(path_6, dpi=150, bbox_inches="tight")
        plt.close(fig6)
        status_queue.put(("image_ready", (path_6, "6. RSCU Heatmap — Top 30 Divergent")))
    except Exception as e:
        print(f"  ❌ Chart 6 error: {e}")

    status_queue.put(("progress", 100))
    print("  ✅ Comparative RSCU analysis complete.")


def rscu_correlation_analysis(all_data, output_folder, status_queue, palette='viridis'):
    print(f"\n=== RSCU CORRELATION ANALYSIS ===")
    all_rscu_data = {species: data['rscu'] for species, data in all_data.items()}
    
    status_queue.put(("progress", 60))
    df_rscu_matrix = pd.DataFrame.from_dict(all_rscu_data, orient='index', columns=ALL_CODONS_SORTED).fillna(0.0)
    
    if len(df_rscu_matrix) < 2:
        print("  ❌ Need at least 2 species for correlation.")
        return

    species_x, species_y = df_rscu_matrix.index[0], df_rscu_matrix.index[1]
    rscu_x, rscu_y = df_rscu_matrix.iloc[0], df_rscu_matrix.iloc[1]

    print(f"  Calculating correlation between '{species_x}' and '{species_y}'...")
    r, p_value = pearsonr(rscu_x, rscu_y)
    
    print(f"    📊 Pearson Coefficient (R): {r:.4f}")
    print(f"    📊 P-value: {p_value:.2e}")
    
    delta_rscu = np.abs(rscu_x - rscu_y)
    
    df_corr_details = pd.DataFrame({
        'Codon': ALL_CODONS_SORTED,
        f'RSCU_{species_x}': rscu_x.values,
        f'RSCU_{species_y}': rscu_y.values,
        'Delta_RSCU (Absolute Diff)': delta_rscu.values,
        'Global_Pearson_R': [r] * len(ALL_CODONS_SORTED),
        'Global_P_value': [p_value] * len(ALL_CODONS_SORTED)
    })
    
    csv_corr_path = os.path.join(output_folder, 'rscu_correlation_details.csv')
    df_corr_details.to_csv(csv_corr_path, index=False, sep=';')
    print(f"  ✅ Detailed correlation table saved in: {csv_corr_path}")

    try:
        plot_color_scatter = sns.color_palette(palette)[0]
        plot_color_line = sns.color_palette(palette)[-1]
    except:
        plot_color_scatter = 'blue'
        plot_color_line = 'darkred'

    print("  Generating Chart 1: Regression Scatter...")
    try:
        plt.figure(figsize=(10, 8))
        ax = sns.regplot(x=rscu_x, y=rscu_y, scatter_kws={'alpha': 0.6, 's': 50, 'color': plot_color_scatter}, line_kws={'color': plot_color_line, 'linewidth': 2})
        ax.set_xlabel(f"RSCU - {species_x}", fontsize=12)
        ax.set_ylabel(f"RSCU - {species_y}", fontsize=12)
        ax.set_title(f"1. RSCU Correlation Regression\n(R={r:.3f}, p-value={p_value:.2e})", fontsize=16)
        plt.grid(alpha=0.3)
        
        out1 = os.path.join(output_folder, 'corr_1_regression_scatter.png')
        plt.savefig(out1, dpi=150, bbox_inches="tight")
        plt.close()
        status_queue.put(("image_ready", (out1, "1. Correlation Scatter")))
    except Exception as e:
        print(f"  ❌ Error: {e}")

    print("  Generating Chart 2: Delta RSCU Barplot...")
    try:
        delta_sorted = delta_rscu.sort_values(ascending=False)
        plt.figure(figsize=(18, 6))
        sns.barplot(x=list(delta_sorted.index), y=list(delta_sorted.values),
                    hue=list(delta_sorted.index), palette=palette, legend=False)
        plt.title(f"2. Absolute Difference (Delta RSCU) per Codon: {species_x} vs {species_y}", fontsize=16)
        plt.xticks(rotation=90, fontsize=10)
        plt.ylabel("Delta |RSCU|")
        plt.grid(axis='y', alpha=0.3)
        plt.tight_layout()
        
        out2 = os.path.join(output_folder, 'corr_2_delta_rscu.png')
        plt.savefig(out2, dpi=150, bbox_inches="tight")
        plt.close()
        status_queue.put(("image_ready", (out2, "2. Delta Differences Barplot")))
    except Exception as e:
        print(f"  ❌ Error: {e}")

    print("  Generating Chart 3: Top Divergent Heatmap...")
    try:
        top_20_codons = delta_sorted.head(20).index
        df_top_divergent = df_rscu_matrix[top_20_codons].T
        
        plt.figure(figsize=(10, 8))
        sns.heatmap(df_top_divergent, annot=True, cmap=palette, fmt=".2f", linewidths=0.5)
        plt.title("3. Top 20 Most Divergent Codons Heatmap", fontsize=16)
        plt.tight_layout()
        
        out3 = os.path.join(output_folder, 'corr_3_divergent_heatmap.png')
        plt.savefig(out3, dpi=150, bbox_inches="tight")
        plt.close()
        status_queue.put(("image_ready", (out3, "3. Top Divergent Heatmap")))
    except Exception as e:
        print(f"  ❌ Error: {e}")


def generate_rscu_histograms(all_data, output_folder, genetic_code_id, status_queue, palette='viridis'):
    """Analysis 6 — redirects to the unified comparative_rscu_analysis module."""
    comparative_rscu_analysis(all_data, output_folder, status_queue, palette=palette)


def enc_gc3_analysis(all_data, output_folder, status_queue, file_list=None, genetic_code_id=1, gene_list=None, palette='viridis'):
    print(f"\n=== ENC vs GC3 ANALYSIS ===")
    status_queue.put(("progress", 60))

    plot_data = []
    for species, data in all_data.items():
        plot_data.append({'Species': species, 'ENC': data['enc'], 'GC3': data['gc3']})
    df_genome = pd.DataFrame(plot_data)

    # Build expected Wright (1990) curve: Nc = 2 + s + 29 / (s² + (1-s)²)
    s_values = np.linspace(0.01, 0.99, 200)
    enc_expected = []
    for s in s_values:
        f_s = s**2 + (1-s)**2
        enc_val = 2 + s + 29 / f_s   # Wright (1990) correct formula
        enc_expected.append(min(enc_val, 61.0))

    # --- Per-gene data (required for the scientifically correct Wright plot) ---
    df_genes = pd.DataFrame()
    if file_list is not None:
        from core_utils import get_w_reference_table, calculate_metrics_per_gene
        gene_rows = []
        for file_path in file_list:
            base_name = os.path.basename(file_path).split('.')[0]
            if base_name not in all_data:
                continue
            w_ref = get_w_reference_table(all_data[base_name]['rscu'], genetic_code_id)
            gene_results = calculate_metrics_per_gene(file_path, w_ref, genetic_code_id, gene_list)
            for res in gene_results:
                gene_rows.append({'Species': base_name, 'Gene': res['gene'],
                                   'ENC': res['enc'], 'GC3': res['gc3']})
        if gene_rows:
            df_genes = pd.DataFrame(gene_rows)
            print(f"  Per-gene Wright plot: {len(df_genes)} genes from {len(file_list)} file(s).")

    # KDE plots use per-gene data when available (richer distribution);
    # barplot and CSV use genome-level summary (one value per species).
    df_plot_density = df_genes if not df_genes.empty else df_genome

    print("  Generating Chart 1: Wright Plot (per-gene)...")
    try:
        plt.figure(figsize=(12, 8))
        plt.plot(s_values * 100, enc_expected, 'r--',
                 label='Expected Curve — Wright (1990)', linewidth=2, zorder=3)

        if not df_genes.empty:
            sns.scatterplot(x='GC3', y='ENC', hue='Species', data=df_genes,
                            palette=palette, s=20, alpha=0.4, edgecolor='none', legend=True)
            title_str = '1. Wright Plot — ENC vs GC3 (per gene)'
        else:
            # Fallback: genome-level points
            sns.scatterplot(x='GC3', y='ENC', hue='Species', data=df_genome,
                            palette=palette, s=150, alpha=0.8, edgecolor='black', legend=False)
            for _, row in df_genome.iterrows():
                plt.annotate(row['Species'], (row['GC3'], row['ENC']),
                             xytext=(5, 5), textcoords='offset points', fontsize=10)
            title_str = '1. Wright Plot — ENC vs GC3 (per genome — pass file_list for per-gene)'

        plt.xlabel('GC3 (%)', fontsize=12)
        plt.ylabel('ENC (Effective Number of Codons)', fontsize=12)
        plt.title(title_str, fontsize=16)
        plt.legend(loc='upper right')
        plt.grid(True, alpha=0.3)
        plt.xlim(0, 100)
        plt.ylim(20, 61)
        plt.tight_layout()

        out1 = os.path.join(output_folder, 'enc_1_wright_plot.png')
        plt.savefig(out1, dpi=150, bbox_inches='tight')
        plt.close()
        status_queue.put(("image_ready", (out1, "1. Wright Plot")))
    except Exception as e:
        print(f"  ❌ Error: {e}")

    try:
        plot_color1 = sns.color_palette(palette)[0]
        plot_color2 = sns.color_palette(palette)[-1]
    except:
        plot_color1 = 'purple'
        plot_color2 = 'green'

    print("  Generating Chart 2: ENC KDE Density...")
    try:
        plt.figure(figsize=(8, 6))
        sns.kdeplot(data=df_plot_density, x='ENC', fill=True, color=plot_color1, alpha=0.5, linewidth=2)
        plt.title("2. Density Distribution of ENC Values", fontsize=14)
        plt.xlabel("Effective Number of Codons (ENC)")
        plt.axvline(61, color='red', linestyle='--', label='Max ENC (No Bias)')
        plt.legend()
        plt.grid(alpha=0.3)

        out2 = os.path.join(output_folder, 'enc_2_kde_density.png')
        plt.savefig(out2, dpi=150, bbox_inches='tight')
        plt.close()
        status_queue.put(("image_ready", (out2, "2. ENC Density Plot")))
    except Exception as e:
        print(f"  ❌ Error: {e}")

    print("  Generating Chart 3: GC3 KDE Density...")
    try:
        plt.figure(figsize=(8, 6))
        sns.kdeplot(data=df_plot_density, x='GC3', fill=True, color=plot_color2, alpha=0.5, linewidth=2)
        plt.title("3. Density Distribution of GC3 Values", fontsize=14)
        plt.xlabel("GC at 3rd Position (%)")
        plt.xlim(0, 100)
        plt.grid(alpha=0.3)

        out3 = os.path.join(output_folder, 'enc_3_gc3_density.png')
        plt.savefig(out3, dpi=150, bbox_inches='tight')
        plt.close()
        status_queue.put(("image_ready", (out3, "3. GC3 Density Plot")))
    except Exception as e:
        print(f"  ❌ Error: {e}")

    print("  Generating Chart 4: ENC Barplot by Species...")
    try:
        plt.figure(figsize=(14, 6))
        df_sorted = df_genome.sort_values('ENC')
        sns.barplot(x='Species', y='ENC', data=df_sorted, hue='Species', palette=palette, legend=False)
        plt.axhline(61, color='red', linestyle='--', label='Theoretical Max (61)')
        plt.axhline(df_genome['ENC'].mean(), color='blue', linestyle='-.', label='Population Mean')
        plt.title("4. Effective Number of Codons (ENC) by Species", fontsize=16)
        plt.xticks(rotation=45, ha='right')
        plt.legend()
        plt.ylim(20, 65)
        plt.grid(axis='y', alpha=0.3)
        plt.tight_layout()

        out4 = os.path.join(output_folder, 'enc_4_barplot_species.png')
        plt.savefig(out4, dpi=150, bbox_inches='tight')
        plt.close()
        status_queue.put(("image_ready", (out4, "4. ENC Barplot by Species")))
    except Exception as e:
        print(f"  ❌ Error: {e}")

    df_genome.to_csv(os.path.join(output_folder, 'enc_gc3_results.csv'), sep=';', index=False)
    if not df_genes.empty:
        df_genes.to_csv(os.path.join(output_folder, 'enc_gc3_per_gene.csv'), sep=';', index=False)
        print(f"  📄 Per-gene CSV saved: enc_gc3_per_gene.csv ({len(df_genes)} rows)")


def optimal_rare_codons_analysis(all_data, output_folder, status_queue, palette='viridis'):
    print(f"\n=== OPTIMAL AND RARE CODONS ANALYSIS ===")
    status_queue.put(("progress", 60))
    
    df_data = []
    for species, data in all_data.items():
        df_data.append({
            'Species': species,
            'CAI': data['cai'],
            'Optimal_Count': len(data['optimal']),
            'Rare_Count': len(data['rare'])
        })
    df_stats = pd.DataFrame(df_data)

    try:
        col1 = sns.color_palette(palette)[0]
        col2 = sns.color_palette(palette)[-1]
    except:
        col1, col2 = 'green', 'red'

    print("  Generating Chart 1: Barplot Opt/Rare...")
    try:
        df_melt = df_stats.melt(id_vars='Species', value_vars=['Optimal_Count', 'Rare_Count'], 
                                var_name='Type', value_name='Count')
        plt.figure(figsize=(14, 6))
        sns.barplot(x='Species', y='Count', hue='Type', data=df_melt, palette=[col1, col2])
        plt.title('1. Number of Optimal (>1.2) vs Rare (<0.8) Codons', fontsize=16, fontweight='bold')
        plt.xticks(rotation=45, ha="right")
        plt.grid(axis='y', alpha=0.3)
        plt.tight_layout()
        
        out1 = os.path.join(output_folder, 'opt_1_barplot.png')
        plt.savefig(out1, dpi=150, bbox_inches='tight')
        plt.close()
        status_queue.put(("image_ready", (out1, "1. Optimal/Rare Barplot")))
    except Exception as e:
        print(f"  ❌ Error: {e}")

    print("  Generating Chart 2: CAI Barplot...")
    try:
        plt.figure(figsize=(14, 6))
        df_sorted_cai = df_stats.sort_values('CAI')
        bars = sns.barplot(x='Species', y='CAI', data=df_sorted_cai, hue='Species', palette=palette, legend=False)
        plt.title('2. Codon Adaptation Index (CAI) by Species', fontsize=16, fontweight='bold')
        plt.xticks(rotation=45, ha="right")
        plt.grid(axis='y', alpha=0.3)
        
        for bar in bars.patches:
            bars.annotate(format(bar.get_height(), '.3f'), 
                          (bar.get_x() + bar.get_width() / 2., bar.get_height()), 
                          ha = 'center', va = 'center', xytext = (0, 9), textcoords = 'offset points')
        
        plt.tight_layout()
        out2 = os.path.join(output_folder, 'opt_2_cai_barplot.png')
        plt.savefig(out2, dpi=150, bbox_inches='tight')
        plt.close()
        status_queue.put(("image_ready", (out2, "2. CAI Barplot")))
    except Exception as e:
        print(f"  ❌ Error: {e}")

    print("  Generating Chart 3: CAI vs Optimal Count Scatter...")
    try:
        plt.figure(figsize=(10, 8))
        sns.scatterplot(x='Optimal_Count', y='CAI', hue='Species', data=df_stats, palette=palette, s=150, edgecolor='black')
        plt.title("3. CAI vs Total Optimal Codons Count", fontsize=16)
        plt.xlabel("Number of Optimal Codons")
        plt.ylabel("Codon Adaptation Index (CAI)")
        plt.grid(alpha=0.3)
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        
        out3 = os.path.join(output_folder, 'opt_3_cai_scatter.png')
        plt.savefig(out3, dpi=150, bbox_inches='tight')
        plt.close()
        status_queue.put(("image_ready", (out3, "3. CAI vs Optimal Scatter")))
    except Exception as e:
        print(f"  ❌ Error: {e}")

    print("  Generating Chart 4: Optimal Codons RSCU Boxplot...")
    try:
        opt_rscu_data = []
        for sp, data in all_data.items():
            for aa, cod in data['optimal'].items():
                opt_rscu_data.append({'Species': sp, 'RSCU': data['rscu'][cod]})
                
        df_opt_rscu = pd.DataFrame(opt_rscu_data)
        
        plt.figure(figsize=(14, 6))
        sns.boxplot(x='Species', y='RSCU', data=df_opt_rscu, hue='Species', palette=palette, legend=False, showfliers=False)
        plt.title("4. Distribution of RSCU Values for Optimal Codons Only", fontsize=16)
        plt.axhline(1.5, color='red', linestyle='--', label='Optimal Threshold')
        plt.xticks(rotation=45, ha='right')
        plt.grid(axis='y', alpha=0.3)
        plt.legend()
        plt.tight_layout()
        
        out4 = os.path.join(output_folder, 'opt_4_optimal_rscu_boxplot.png')
        plt.savefig(out4, dpi=150, bbox_inches='tight')
        plt.close()
        status_queue.put(("image_ready", (out4, "4. Optimal RSCU Boxplot")))
    except Exception as e:
        print(f"  ❌ Error: {e}")

    for species_name, data in all_data.items():
        df_optimal = pd.DataFrame(list(data['optimal'].items()), columns=['Amino_Acid', 'Optimal_Codon'])
        df_rare = pd.DataFrame(list(data['rare'].items()), columns=['Amino_Acid', 'Rare_Codon'])
        df_optimal.to_csv(os.path.join(output_folder, f'{species_name}_optimal_codons.csv'), sep=';', index=False)
        df_rare.to_csv(os.path.join(output_folder, f'{species_name}_rare_codons.csv'), sep=';', index=False)


def neutrality_plot_analysis(all_data, output_folder, status_queue, file_list=None, genetic_code_id=1, gene_list=None, palette='viridis'):
    print(f"\n=== NEUTRALITY PLOT ANALYSIS ===")
    status_queue.put(("progress", 60))

    data_list = []
    for species, data in all_data.items():
        if 'gc12' in data and 'gc3' in data:
            counts = data['counts']
            total_codons = sum(counts.values())
            if total_codons == 0:
                continue
            gc1_count = sum(v for k, v in counts.items() if len(k) == 3 and k[0] in 'GC')
            gc2_count = sum(v for k, v in counts.items() if len(k) == 3 and k[1] in 'GC')
            data_list.append({
                'Species': species,
                'GC12': data['gc12'],
                'GC3': data['gc3'],
                'GC1': (gc1_count / total_codons) * 100,
                'GC2': (gc2_count / total_codons) * 100
            })

    if not data_list:
        print("  ❌ Error: GC12 or GC3 data not found.")
        return

    df_genome = pd.DataFrame(data_list).replace([np.inf, -np.inf], np.nan).dropna()

    # --- Per-gene data for the scientifically correct Neutrality Plot ---
    df_genes = pd.DataFrame()
    if file_list is not None:
        from core_utils import get_w_reference_table, calculate_metrics_per_gene
        gene_rows = []
        for file_path in file_list:
            base_name = os.path.basename(file_path).split('.')[0]
            if base_name not in all_data:
                continue
            w_ref = get_w_reference_table(all_data[base_name]['rscu'], genetic_code_id)
            gene_results = calculate_metrics_per_gene(file_path, w_ref, genetic_code_id, gene_list)
            for res in gene_results:
                gene_rows.append({'Species': base_name, 'Gene': res['gene'],
                                   'GC12': res['gc12'], 'GC3': res['gc3']})
        if gene_rows:
            df_genes = pd.DataFrame(gene_rows).replace([np.inf, -np.inf], np.nan).dropna()
            print(f"  Per-gene neutrality plot: {len(df_genes)} genes from {len(file_list)} file(s).")

    print("  Generating Chart 1: Neutrality Plot (GC12 vs GC3)...")
    try:
        plt.figure(figsize=(10, 8))

        if not df_genes.empty:
            # Per-gene scatter + per-species regression lines.
            # Single legend entry per species (scatter has no label; the regression
            # line carries species name + slope + R²).
            # Reference line y = x (slope=1): pure mutational neutrality (Sueoka).
            plt.plot([0, 100], [0, 100], color='black', linestyle='--',
                     linewidth=1.5, alpha=0.6,
                     label='y = x  (slope=1, pure mutational neutrality)', zorder=2)

            species_list = df_genes['Species'].unique()
            colors = sns.color_palette(palette, len(species_list))
            for sp, col in zip(species_list, colors):
                sub = df_genes[df_genes['Species'] == sp]
                plt.scatter(sub['GC3'], sub['GC12'], color=col, s=15, alpha=0.35,
                            edgecolor='none')
                if len(sub) >= 3:
                    sl, intercept, r_val, p_val, _ = linregress(sub['GC3'], sub['GC12'])
                    x_vals = np.linspace(sub['GC3'].min(), sub['GC3'].max(), 100)
                    plt.plot(x_vals, intercept + sl * x_vals, color=col, linewidth=2,
                             label=f'{sp}  slope={sl:.3f}  R²={r_val**2:.3f}')
                    print(f"  [{sp}] slope={sl:.3f}, R²={r_val**2:.3f}, p={p_val:.4f}")
                else:
                    # Still include the species in the legend even without regression
                    plt.plot([], [], color=col, linewidth=2, label=f'{sp}  (N<3)')
            title_str = '1. Neutrality Plot — GC12 vs GC3 (per gene, Sueoka 1992)'
        else:
            # Fallback: genome-level
            if len(df_genome) >= 3:
                slope, intercept, r_value, p_value, _ = linregress(df_genome['GC3'], df_genome['GC12'])
                x_vals = np.array([0, 100])
                plt.plot(x_vals, intercept + slope * x_vals, 'r--',
                         label=f'Regression Fit (Slope = {slope:.3f})')
                print(f"  Genome-level regression: slope={slope:.3f}, R²={r_value**2:.3f}")
            sns.scatterplot(x='GC3', y='GC12', data=df_genome, palette=palette,
                            s=150, alpha=0.8, hue='Species', edgecolor='black')
            for _, row in df_genome.iterrows():
                plt.annotate(row['Species'], (row['GC3'], row['GC12']),
                             xytext=(5, 5), textcoords='offset points', fontsize=10)
            title_str = '1. Neutrality Plot — GC12 vs GC3 (per genome — pass file_list for per-gene)'

        plt.xlabel('GC3 (%)', fontsize=12)
        plt.ylabel('GC12 (%)', fontsize=12)
        plt.title(title_str, fontsize=15)
        # Legend outside the plot area on the right so it never covers the data
        plt.legend(loc='upper left', bbox_to_anchor=(1.02, 1.0),
                   fontsize=9, frameon=True, borderaxespad=0.0)
        plt.grid(True, alpha=0.3)
        plt.xlim(0, 100)
        plt.ylim(0, 100)
        plt.tight_layout()

        out1 = os.path.join(output_folder, 'neu_1_gc12_vs_gc3.png')
        plt.savefig(out1, dpi=150, bbox_inches='tight')
        plt.close()
        status_queue.put(("image_ready", (out1, "1. Neutrality Plot (GC12 vs GC3)")))
    except Exception as e:
        print(f"  ❌ Error: {e}")

    try:
        col1 = sns.color_palette(palette)[0]
        col2 = sns.color_palette(palette)[-1]
    except:
        col1, col2 = 'blue', 'green'

    print("  Generating Chart 2: GC1 and GC2 Regressions...")
    try:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
        
        sns.regplot(x='GC3', y='GC1', data=df_genome, ax=ax1, color=col1, scatter_kws={'s': 80})
        ax1.set_title("2A. GC1 vs GC3", fontsize=14)
        ax1.set_xlim(0, 100); ax1.set_ylim(0, 100)
        ax1.grid(alpha=0.3)
        
        sns.regplot(x='GC3', y='GC2', data=df_genome, ax=ax2, color=col2, scatter_kws={'s': 80})
        ax2.set_title("2B. GC2 vs GC3", fontsize=14)
        ax2.set_xlim(0, 100); ax2.set_ylim(0, 100)
        ax2.grid(alpha=0.3)
        
        plt.suptitle("2. GC Frequency at 1st and 2nd Positions relative to 3rd Position", fontsize=16)
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        
        out2 = os.path.join(output_folder, 'neu_2_gc1_gc2_vs_gc3.png')
        plt.savefig(out2, dpi=150, bbox_inches='tight')
        plt.close()
        status_queue.put(("image_ready", (out2, "2. GC1 and GC2 Regressions")))
    except Exception as e:
        print(f"  ❌ Error: {e}")

    print("  Generating Chart 3: Position Specific GC Barplot...")
    try:
        df_melt = df_genome.melt(id_vars='Species', value_vars=['GC1', 'GC2', 'GC3'],
                                var_name='Position', value_name='GC%')
        plt.figure(figsize=(16, 6))
        sns.barplot(x='Species', y='GC%', hue='Position', data=df_melt, palette=palette)
        plt.title("3. GC Content by Codon Position", fontsize=16)
        plt.xticks(rotation=45, ha='right')
        plt.grid(axis='y', alpha=0.3)
        plt.tight_layout()
        
        out3 = os.path.join(output_folder, 'neu_3_gc_positions_barplot.png')
        plt.savefig(out3, dpi=150, bbox_inches='tight')
        plt.close()
        status_queue.put(("image_ready", (out3, "3. Position GC Barplot")))
    except Exception as e:
        print(f"  ❌ Error: {e}")

    print("  Generating Chart 4: KDE Densities...")
    try:
        plt.figure(figsize=(10, 6))
        sns.kdeplot(data=df_genome, x='GC12', fill=True, color=col1, label='GC12 Density', alpha=0.5)
        sns.kdeplot(data=df_genome, x='GC3', fill=True, color=col2, label='GC3 Density', alpha=0.5)
        plt.title("4. Density Distribution of GC12 and GC3 Values", fontsize=16)
        plt.xlabel("GC Content (%)")
        plt.xlim(0, 100)
        plt.legend()
        plt.grid(alpha=0.3)
        plt.tight_layout()
        
        out4 = os.path.join(output_folder, 'neu_4_gc_density_kde.png')
        plt.savefig(out4, dpi=150, bbox_inches='tight')
        plt.close()
        status_queue.put(("image_ready", (out4, "4. GC12/GC3 KDE Densities")))
    except Exception as e:
        print(f"  ❌ Error: {e}")

    df_genome.to_csv(os.path.join(output_folder, 'neutrality_plot_results.csv'), sep=';', index=False)
    if not df_genes.empty:
        df_genes.to_csv(os.path.join(output_folder, 'neutrality_plot_per_gene.csv'), sep=';', index=False)
        print(f"  📄 Per-gene CSV saved: neutrality_plot_per_gene.csv ({len(df_genes)} rows)")
