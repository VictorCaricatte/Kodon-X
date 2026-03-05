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

def comparative_rscu_analysis(all_data, output_folder, status_queue, palette='viridis'):
    print(f"\n=== COMPARATIVE RSCU ANALYSIS ===")
    all_rscu_data = {species: data['rscu'] for species, data in all_data.items()}

    print("\n  Creating comparative RSCU matrix...")
    status_queue.put(("progress", 55))
    df_rscu_matrix = pd.DataFrame.from_dict(all_rscu_data, orient='index', columns=ALL_CODONS_SORTED).fillna(0.0)
    csv_matrix_path = os.path.join(output_folder, 'comparative_rscu_matrix.csv')
    df_rscu_matrix.to_csv(csv_matrix_path, sep=';', decimal='.')
    print(f"  ✅ RSCU Matrix saved in: {csv_matrix_path}")
    
    print("  Generating Chart 1: Clustermap...")
    status_queue.put(("progress", 70))
    try:
        g = sns.clustermap(
            df_rscu_matrix, metric="euclidean", method="average", cmap=palette, annot=False, linewidths=0.5,
            figsize=(max(15, len(ALL_CODONS_SORTED) * 0.2), max(8, len(all_data) * 0.5))
        )
        g.fig.suptitle("1. Comparative RSCU Analysis (Clustermap)", y=1.02, fontsize=16)
        plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90)
        
        clustermap_path = os.path.join(output_folder, 'comparative_1_clustermap.png')
        g.savefig(clustermap_path, dpi=150, bbox_inches="tight")
        plt.close()
        status_queue.put(("image_ready", (clustermap_path, "1. RSCU Clustermap")))
    except Exception as e:
        print(f"\n❌ ERROR GENERATING CLUSTERMAP: {e}")

    print("  Generating PCA Data...")
    status_queue.put(("progress", 80))
    try:
        X_scaled = StandardScaler().fit_transform(df_rscu_matrix.values)
        pca = PCA(n_components=2)
        principal_components = pca.fit_transform(X_scaled)
        df_pca = pd.DataFrame(data=principal_components, columns=['PC1', 'PC2'], index=df_rscu_matrix.index)
        pc1_var, pc2_var = pca.explained_variance_ratio_ * 100
        
        print("  Generating Chart 2: PCA Samples Space...")
        plt.figure(figsize=(12, 10))
        try:
            plot_color = sns.color_palette(palette)[0]
        except:
            plot_color = 'teal'
            
        sns.scatterplot(x='PC1', y='PC2', data=df_pca, s=150, alpha=0.7, edgecolor='k', color=plot_color)
        
        for i, sample in enumerate(df_pca.index):
            plt.text(df_pca.iloc[i]['PC1'] + 0.05, df_pca.iloc[i]['PC2'], sample, fontsize=10)
            
        plt.xlabel(f'Principal Component 1 ({pc1_var:.2f}%)', fontsize=12)
        plt.ylabel(f'Principal Component 2 ({pc2_var:.2f}%)', fontsize=12)
        plt.title('2. RSCU Principal Component Analysis (Samples Space)', fontsize=16)
        plt.axhline(0, color='grey', linestyle='--', linewidth=1)
        plt.axvline(0, color='grey', linestyle='--', linewidth=1)
        
        pca_samples_path = os.path.join(output_folder, 'comparative_2_pca_samples.png')
        plt.savefig(pca_samples_path, dpi=150, bbox_inches="tight")
        plt.close()
        status_queue.put(("image_ready", (pca_samples_path, "2. PCA Samples Space")))

        print("  Generating Chart 3: PCA Biplot (Codons)...")
        plt.figure(figsize=(16, 16))
        loadings = pca.components_.T
        max_loading = np.max(np.abs(loadings))
        
        for i, codon in enumerate(df_rscu_matrix.columns):
            x_val = loadings[i, 0]
            y_val = loadings[i, 1]
            plt.arrow(0, 0, x_val, y_val, color=plot_color, alpha=0.6, 
                      head_width=max_loading*0.015, head_length=max_loading*0.015, linewidth=1.5)
            
            plt.text(x_val * 1.08, y_val * 1.08, codon, color=plot_color, 
                     ha='center', va='center', fontsize=10, fontweight='bold')
            
        plt.xlabel(f'Component 1 ({pc1_var:.2f}%)', fontsize=12)
        plt.ylabel(f'Component 2 ({pc2_var:.2f}%)', fontsize=12)
        plt.title('3. RSCU Principal Component Analysis (Codon Loadings / Biplot)', fontsize=16)
        plt.axhline(0, color='grey', linestyle='--', linewidth=1)
        plt.axvline(0, color='grey', linestyle='--', linewidth=1)
        
        limit = max_loading * 1.25
        plt.xlim(-limit, limit)
        plt.ylim(-limit, limit)
        plt.grid(alpha=0.3)
        
        pca_biplot_path = os.path.join(output_folder, 'comparative_3_pca_biplot.png')
        plt.savefig(pca_biplot_path, dpi=150, bbox_inches="tight")
        plt.close()
        status_queue.put(("image_ready", (pca_biplot_path, "3. PCA Biplot (Loadings)")))

    except Exception as e:
        print(f"\n❌ ERROR GENERATING PCA: {e}")

    print("  Generating Chart 4: Variance Boxplot...")
    status_queue.put(("progress", 95))
    try:
        plt.figure(figsize=(24, 8))
        df_long_variance = df_rscu_matrix.melt(var_name='Codon', value_name='RSCU')
        sns.boxplot(x='Codon', y='RSCU', data=df_long_variance, palette=palette, showfliers=False)
        plt.title("4. RSCU Variance Across Species per Codon", fontsize=16)
        plt.xticks(rotation=90, fontsize=10)
        plt.axhline(1.0, color='red', linestyle='--', alpha=0.5)
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        
        variance_path = os.path.join(output_folder, 'comparative_4_variance_boxplot.png')
        plt.savefig(variance_path, dpi=150, bbox_inches="tight")
        plt.close()
        status_queue.put(("image_ready", (variance_path, "4. RSCU Variance Boxplot")))
    except Exception as e:
        print(f"\n❌ ERROR GENERATING VARIANCE PLOT: {e}")


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
        sns.barplot(x=delta_sorted.index, y=delta_sorted.values, palette=palette)
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
    print(f"\n=== RSCU HISTOGRAMS ANALYSIS ===")
    all_rscu_data = {species: data['rscu'] for species, data in all_data.items()}
    
    status_queue.put(("progress", 40))
    df_rscu_matrix = pd.DataFrame.from_dict(all_rscu_data, orient='index', columns=ALL_CODONS_SORTED).fillna(0.0)
    
    df_long = df_rscu_matrix.reset_index().rename(columns={'index': 'Species'}).melt(id_vars='Species', var_name='Codon', value_name='RSCU')
    codon_aa_map = GENETIC_CODE_TABLES.get(genetic_code_id, GENETIC_CODE_TABLES[1])
    df_long['AminoAcid'] = df_long['Codon'].map(codon_aa_map)
    df_long = df_long.sort_values(by=['AminoAcid', 'Codon'])

    aa_codon_map = AA_CODON_MAPS.get(genetic_code_id, AA_CODON_MAPS[1])
    
    csv_hist_path = os.path.join(output_folder, 'rscu_histograms_data.csv')
    df_long.to_csv(csv_hist_path, sep=';', index=False)
    print(f"  ✅ RSCU Histogram data saved in: {csv_hist_path}")
    
    print(f"  📊 Processing {len(all_data)} species for histogram analysis")

    print("  Generating Chart 1: RSCU Box Plot by Amino Acid...")
    status_queue.put(("progress", 50))
    try:
        plt.figure(figsize=(20, 10))
        df_filtered = df_long[df_long['AminoAcid'].isin(aa_codon_map.keys())]
        df_filtered = df_filtered[df_filtered['AminoAcid'] != '*']
        
        aa_order = sorted(df_filtered['AminoAcid'].unique())
        
        sns.boxplot(x='AminoAcid', y='RSCU', data=df_filtered, order=aa_order, palette="Set3", showfliers=False)
        
        plt.xlabel('Amino Acids', fontsize=12)
        plt.ylabel('RSCU Value', fontsize=12)
        plt.title('1. RSCU Distribution by Amino Acid (All Species)', fontsize=16, fontweight='bold')
        plt.grid(axis='y', alpha=0.3)
        plt.axhline(y=1.0, color='red', linestyle='--', alpha=0.5)
        plt.tight_layout()
        
        out1 = os.path.join(output_folder, 'hist_1_boxplot_by_aminoacid.png')
        plt.savefig(out1, dpi=150, bbox_inches="tight")
        plt.close()
        status_queue.put(("image_ready", (out1, "1. RSCU Box Plot by AA")))
    except Exception as e:
        print(f"\n❌ ERROR GENERATING CHART 1: {e}")

    print("\n  Generating Chart 2: Stacked Vertical Bars (High Contrast Colors)...")
    status_queue.put(("progress", 70))
    try:
        aa_list = sorted([aa for aa, codons in aa_codon_map.items() if aa != '*' and len(codons) > 1]) 
        
        all_used_codons = []
        for aa in aa_list:
            all_used_codons.extend(sorted(aa_codon_map[aa]))
            
        cmap = plt.get_cmap('gist_ncar')
        distinct_colors = [mcolors.to_hex(cmap(i/len(all_used_codons))) for i in range(len(all_used_codons))]
        custom_cmap = mcolors.ListedColormap(distinct_colors)
        
        for species_name, rscu_series in df_rscu_matrix.iterrows():
            fig, ax = plt.subplots(figsize=(18, 10))
            plot_data = []
            for aa in aa_list:
                syn_codons = sorted(aa_codon_map[aa])
                for codon in syn_codons:
                    value = rscu_series.get(codon, 0)
                    plot_data.append({'AminoAcid': aa, 'Codon': codon, 'RSCU': value})
            
            df_plot = pd.DataFrame(plot_data)

            df_plot.pivot(index='AminoAcid', columns='Codon', values='RSCU').loc[aa_list].plot(
                kind='bar', stacked=True, ax=ax,
                colormap=custom_cmap, edgecolor='black', linewidth=0.5
            )

            ax.set_xlabel("Amino Acid", fontsize=12, fontweight='bold')
            ax.set_ylabel("RSCU Value (Stacked)", fontsize=12)
            ax.set_title(f"2. Codon Usage Profile (RSCU) - {species_name}", fontsize=16, fontweight='bold')
            ax.legend(title='Codons', loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=10, fontsize='small')
            ax.grid(axis='y', linestyle='--', alpha=0.7)
            plt.xticks(rotation=0)
            plt.tight_layout() 
            
            out2 = os.path.join(output_folder, f'hist_2_stacked_vertical_{species_name}.png')
            plt.savefig(out2, dpi=150, bbox_inches="tight")
            plt.close()
        status_queue.put(("image_ready", (out2, "2. High-Contrast Stacked Bars")))
    except Exception as e:
        print(f"\n❌ ERROR GENERATING CHART 2: {e}")

    print("\n  Generating Chart 3: RSCU Line Plot with Grouping...")
    status_queue.put(("progress", 85))
    try:
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

        for species_name in df_rscu_matrix.index:
            y_vals = [df_rscu_matrix.loc[species_name, c] for c in ordered_codons]
            plt.plot(range(len(ordered_codons)), y_vals, label=species_name, linewidth=2.5, marker='o', markersize=5, zorder=3)

        plt.axhline(1.5, color='darkred', linestyle='--', alpha=0.8, label="Optimal (>1.5)", zorder=2)
        plt.axhline(0.5, color='darkred', linestyle='--', alpha=0.8, label="Rare (<0.5)", zorder=2)
        
        labels = [f"{c}\n({a})" for c, a in zip(ordered_codons, ordered_aas)]
        plt.xticks(range(len(ordered_codons)), labels, rotation=90, fontsize=11)
        
        plt.xlabel("Codon / Amino Acid", fontsize=13, fontweight='bold')
        plt.ylabel("RSCU Value", fontsize=13, fontweight='bold')
        plt.title("3. RSCU Comparative Distribution Profile Across All Codons", fontsize=18, fontweight='bold')
        
        plt.grid(axis='y', linestyle=':', alpha=0.6, zorder=1)
        plt.legend(bbox_to_anchor=(1.01, 1), loc='upper left')
        plt.margins(x=0.01)
        plt.tight_layout()

        out3 = os.path.join(output_folder, 'hist_3_comparative_lineplot.png')
        plt.savefig(out3, dpi=150, bbox_inches="tight")
        plt.close()
        status_queue.put(("image_ready", (out3, "3. RSCU Comparative Line Plot")))

    except Exception as e:
        print(f"\n❌ ERROR GENERATING CHART 3: {e}")

    print("\n  Generating Chart 4: Matrix Heatmap...")
    status_queue.put(("progress", 95))
    try:
        plt.figure(figsize=(22, max(6, len(df_rscu_matrix)*0.5)))
        sns.heatmap(df_rscu_matrix[ordered_codons], cmap=palette, linewidths=0.5, cbar_kws={'label': 'RSCU'})
        plt.title("4. Global Comparative RSCU Heatmap (Grouped by Amino Acid)", fontsize=16)
        plt.xticks(rotation=90, fontsize=10)
        plt.tight_layout()
        
        out4 = os.path.join(output_folder, 'hist_4_matrix_heatmap.png')
        plt.savefig(out4, dpi=150, bbox_inches="tight")
        plt.close()
        status_queue.put(("image_ready", (out4, "4. Global RSCU Matrix")))
    except Exception as e:
        print(f"\n❌ ERROR GENERATING CHART 4: {e}")


def enc_gc3_analysis(all_data, output_folder, status_queue, palette='viridis'):
    print(f"\n=== ENC vs GC3 ANALYSIS ===")
    status_queue.put(("progress", 60))
    
    plot_data = []
    for species, data in all_data.items():
        plot_data.append({'Species': species, 'ENC': data['enc'], 'GC3': data['gc3']})
    df_plot = pd.DataFrame(plot_data)

    print("  Generating Chart 1: Wright Plot...")
    try:
        plt.figure(figsize=(12, 8))
        s_values = np.linspace(0.01, 0.99, 200)
        enc_expected = []
        for s in s_values:
            f_s = s**2 + (1-s)**2
            enc_val = 2 + (9/f_s) + (1/f_s) + (5/f_s) + (3/f_s)
            enc_expected.append(enc_val)
            
        plt.plot(s_values * 100, enc_expected, 'r--', label='Expected Curve (Mutational Bias)', linewidth=2)
        
        sns.scatterplot(x='GC3', y='ENC', hue='Species', data=df_plot, palette=palette, s=150, alpha=0.8, edgecolor='black', legend=False)
        
        for i, row in df_plot.iterrows():
            plt.annotate(row['Species'], (row['GC3'], row['ENC']), 
                        xytext=(5, 5), textcoords='offset points', fontsize=10)
        
        plt.xlabel('GC3 (%)', fontsize=12)
        plt.ylabel('ENC (Effective Number of Codons)', fontsize=12)
        plt.title('1. ENC vs GC3 Analysis (Wright Plot)', fontsize=16)
        plt.legend()
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
        sns.kdeplot(data=df_plot, x='ENC', fill=True, color=plot_color1, alpha=0.5, linewidth=2)
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
        sns.kdeplot(data=df_plot, x='GC3', fill=True, color=plot_color2, alpha=0.5, linewidth=2)
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

    print("  Generating Chart 4: ENC Barplot...")
    try:
        plt.figure(figsize=(14, 6))
        df_sorted = df_plot.sort_values('ENC')
        sns.barplot(x='Species', y='ENC', data=df_sorted, palette=palette)
        plt.axhline(61, color='red', linestyle='--', label='Theoretical Max (61)')
        plt.axhline(df_plot['ENC'].mean(), color='blue', linestyle='-.', label='Population Mean')
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

    df_plot.to_csv(os.path.join(output_folder, 'enc_gc3_results.csv'), sep=';', index=False)


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
        bars = sns.barplot(x='Species', y='CAI', data=df_sorted_cai, palette=palette)
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
        sns.boxplot(x='Species', y='RSCU', data=df_opt_rscu, palette=palette, showfliers=False)
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


def neutrality_plot_analysis(all_data, output_folder, status_queue, palette='viridis'):
    print(f"\n=== NEUTRALITY PLOT ANALYSIS ===")
    status_queue.put(("progress", 60))
    
    data_list = []
    for species, data in all_data.items():
        if 'gc12' in data and 'gc3' in data:
            counts = data['counts']
            total_codons = sum(counts.values())
            if total_codons == 0: continue
            
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
        
    df_plot = pd.DataFrame(data_list)
    df_plot = df_plot.replace([np.inf, -np.inf], np.nan).dropna()
    
    if df_plot.empty or len(df_plot) < 2:
        print("  ❌ Error: Insufficient data for linear regression.")
        return

    print("  Generating Chart 1: Neutrality Plot (GC12 vs GC3)...")
    try:
        slope, intercept, r_value, p_value, std_err = linregress(df_plot['GC3'], df_plot['GC12'])
        print(f"  📊 Regression: Slope = {slope:.3f}, R² = {r_value**2:.3f}")

        plt.figure(figsize=(10, 8))
        sns.scatterplot(x='GC3', y='GC12', data=df_plot, palette=palette, s=150, alpha=0.8, hue='Species', legend=False, edgecolor='black')
        
        x_vals = np.array([0, 100])
        y_vals = intercept + slope * x_vals
        plt.plot(x_vals, y_vals, 'r--', label=f'Regression Fit (Slope = {slope:.3f})')
        
        for i, row in df_plot.iterrows():
            plt.annotate(row['Species'], (row['GC3'], row['GC12']), 
                        xytext=(5, 5), textcoords='offset points', fontsize=10)
        
        plt.xlabel('GC3 (%)', fontsize=12)
        plt.ylabel('GC12 (%)', fontsize=12)
        plt.title(f"1. Neutrality Plot (GC12 vs GC3)\nSlope = {slope:.3f}, R² = {r_value**2:.3f}", fontsize=16)
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.xlim(0, 100); plt.ylim(0, 100)
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
        
        sns.regplot(x='GC3', y='GC1', data=df_plot, ax=ax1, color=col1, scatter_kws={'s': 80})
        ax1.set_title("2A. GC1 vs GC3", fontsize=14)
        ax1.set_xlim(0, 100); ax1.set_ylim(0, 100)
        ax1.grid(alpha=0.3)
        
        sns.regplot(x='GC3', y='GC2', data=df_plot, ax=ax2, color=col2, scatter_kws={'s': 80})
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
        df_melt = df_plot.melt(id_vars='Species', value_vars=['GC1', 'GC2', 'GC3'], 
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
        sns.kdeplot(data=df_plot, x='GC12', fill=True, color=col1, label='GC12 Density', alpha=0.5)
        sns.kdeplot(data=df_plot, x='GC3', fill=True, color=col2, label='GC3 Density', alpha=0.5)
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

    df_plot.to_csv(os.path.join(output_folder, 'neutrality_plot_results.csv'), sep=';', index=False)
