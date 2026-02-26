import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr, linregress
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from constants import GENETIC_CODE_TABLES, CODON_GRID_ORDER, ALL_CODONS_SORTED, AA_CODON_MAPS

def generate_rscu_heatmap_and_table(all_data, output_folder, genetic_code_id, status_queue):
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
    print("  Preparing data for heatmap...")
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
    
    print("  Generating RSCU Heatmap...")
    plt.figure(figsize=(20, 8))
    try:
        ax = sns.heatmap(df_rscu, annot=codon_labels_grid, fmt="", cmap="viridis", linewidths=0.5, cbar_kws={'label': 'RSCU Value'})
        ax.set_yticklabels(['T', 'C', 'A', 'G'], rotation=0)
        ax.set_xticklabels([f"Pos {i+1}" for i in range(16)], rotation=0)
        ax.set_title(f"Relative Synonymous Codon Usage (RSCU) - {base_filename}\nENC={data['enc']:.2f}, GC3={data['gc3']:.2f}%, CAI={data['cai']:.3f}", fontsize=16)
        
        output_file = os.path.join(output_folder, f"detailed_rscu_heatmap_{base_filename}.png")
        plt.savefig(output_file, dpi=150, bbox_inches="tight")
        plt.close()
        
        print(f"\n✅ RSCU Heatmap saved in: {output_file}")
        status_queue.put(("image_ready", (output_file, f"RSCU - {base_filename}")))
        
    except Exception as e:
        print(f"\n❌ ERROR GENERATING RSCU HEATMAP (Analysis 3): {e}")

def comparative_rscu_analysis(all_data, output_folder, status_queue):
    print(f"\n=== COMPARATIVE RSCU ANALYSIS ===")
    all_rscu_data = {species: data['rscu'] for species, data in all_data.items()}

    print("\n  Creating comparative RSCU matrix...")
    status_queue.put(("progress", 55))
    df_rscu_matrix = pd.DataFrame.from_dict(all_rscu_data, orient='index', columns=ALL_CODONS_SORTED).fillna(0.0)
    csv_matrix_path = os.path.join(output_folder, 'comparative_rscu_matrix.csv')
    df_rscu_matrix.to_csv(csv_matrix_path, sep=';', decimal='.')
    print(f"  ✅ RSCU Matrix saved in: {csv_matrix_path}")
    
    print("  Generating Clustermap...")
    status_queue.put(("progress", 70))
    try:
        g = sns.clustermap(
            df_rscu_matrix, metric="euclidean", method="average", cmap="viridis", annot=False, linewidths=0.5,
            figsize=(max(15, len(ALL_CODONS_SORTED) * 0.2), max(8, len(all_data) * 0.5))
        )
        g.fig.suptitle("Comparative RSCU Analysis (Clustermap)", y=1.02, fontsize=16)
        plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90)
        
        clustermap_path = os.path.join(output_folder, 'comparative_rscu_clustermap.png')
        g.savefig(clustermap_path, dpi=150, bbox_inches="tight")
        plt.close()
        
        print(f"  ✅ Clustermap saved in: {clustermap_path}")
        status_queue.put(("image_ready", (clustermap_path, "Comparative RSCU Clustermap")))
        
    except Exception as e:
        print(f"\n❌ ERROR GENERATING CLUSTERMAP: {e}")

    print("  Generating Principal Component Analysis (PCA)...")
    status_queue.put(("progress", 90))
    try:
        X_scaled = StandardScaler().fit_transform(df_rscu_matrix.values)
        pca = PCA(n_components=2)
        principal_components = pca.fit_transform(X_scaled)
        df_pca = pd.DataFrame(data=principal_components, columns=['PC1', 'PC2'], index=df_rscu_matrix.index)
        
        plt.figure(figsize=(12, 10))
        sns.scatterplot(x='PC1', y='PC2', data=df_pca, s=150, alpha=0.7)
        
        for i, sample in enumerate(df_pca.index):
            plt.text(df_pca.iloc[i]['PC1'] + 0.05, df_pca.iloc[i]['PC2'], sample, fontsize=9)
            
        pc1_var, pc2_var = pca.explained_variance_ratio_ * 100
        plt.xlabel(f'Principal Component 1 ({pc1_var:.2f}%)', fontsize=12)
        plt.ylabel(f'Principal Component 2 ({pc2_var:.2f}%)', fontsize=12)
        plt.title('RSCU Principal Component Analysis (PCA)', fontsize=16)
        plt.axhline(0, color='grey', linestyle='--', linewidth=0.5)
        plt.axvline(0, color='grey', linestyle='--', linewidth=0.5)
        
        pca_path = os.path.join(output_folder, 'comparative_rscu_pca.png')
        plt.savefig(pca_path, dpi=150, bbox_inches="tight")
        plt.close()
        
        print(f"  ✅ PCA Chart saved in: {pca_path}")
        print(f"  📈 Explained variance: PC1={pc1_var:.2f}%, PC2={pc2_var:.2f}%")
        status_queue.put(("image_ready", (pca_path, "RSCU PCA Analysis")))
        
    except Exception as e:
        print(f"\n❌ ERROR GENERATING PCA: {e}")

def rscu_correlation_analysis(all_data, output_folder, status_queue):
    print(f"\n=== RSCU CORRELATION ANALYSIS ===")
    all_rscu_data = {species: data['rscu'] for species, data in all_data.items()}
    
    status_queue.put(("progress", 60))
    df_rscu_matrix = pd.DataFrame.from_dict(all_rscu_data, orient='index', columns=ALL_CODONS_SORTED).fillna(0.0)
    
    species_x, species_y = df_rscu_matrix.index[0], df_rscu_matrix.index[1]
    rscu_x, rscu_y = df_rscu_matrix.iloc[0], df_rscu_matrix.iloc[1]

    print(f"  Calculating correlation between '{species_x}' and '{species_y}'...")
    r, p_value = pearsonr(rscu_x, rscu_y)
    
    print(f"    📊 Pearson Coefficient (R): {r:.4f}")
    print(f"    📊 P-value: {p_value:.2e}")

    status_queue.put(("progress", 80))
    print("  Generating correlation chart...")
    
    plt.figure(figsize=(10, 8))
    ax = sns.regplot(x=rscu_x, y=rscu_y,
                     scatter_kws={'alpha': 0.5},
                     line_kws={'color': 'darkred', 'linewidth': 2})
                     
    ax.set_xlabel(f"{species_x}", fontsize=12)
    ax.set_ylabel(f"{species_y}", fontsize=12)
    
    title_str = f"RSCU Correlation: {species_x} vs {species_y}\n(R={r:.3f}, p-value={p_value:.2e})"
    ax.set_title(title_str, fontsize=16)

    correlation_path = os.path.join(output_folder, f'rscu_correlation_{species_x}_vs_{species_y}.png')
    plt.savefig(correlation_path, dpi=150, bbox_inches="tight")
    plt.close()
    
    print(f"  ✅ Correlation chart saved in: {correlation_path}")
    status_queue.put(("image_ready", (correlation_path, f"RSCU Correlation")))

def generate_rscu_histograms(all_data, output_folder, genetic_code_id, status_queue):
    print(f"\n=== RSCU HISTOGRAMS ANALYSIS ===")
    all_rscu_data = {species: data['rscu'] for species, data in all_data.items()}
    
    status_queue.put(("progress", 50))
    df_rscu_matrix = pd.DataFrame.from_dict(all_rscu_data, orient='index', columns=ALL_CODONS_SORTED).fillna(0.0)
    
    df_long = df_rscu_matrix.reset_index().rename(columns={'index': 'Species'}).melt(id_vars='Species', var_name='Codon', value_name='RSCU')
    codon_aa_map = GENETIC_CODE_TABLES.get(genetic_code_id, GENETIC_CODE_TABLES[1])
    df_long['AminoAcid'] = df_long['Codon'].map(codon_aa_map)
    df_long = df_long.sort_values(by=['AminoAcid', 'Codon'])

    print(f"  📊 Processing {len(all_data)} species for histogram analysis")

    print("  Generating Chart 1: RSCU Box Plot by Amino Acid...")
    status_queue.put(("progress", 60))
    try:
        plt.figure(figsize=(20, 10))
        df_filtered = df_long[df_long['AminoAcid'].isin(AA_CODON_MAPS[genetic_code_id].keys())]
        df_filtered = df_filtered[df_filtered['AminoAcid'] != '*']
        
        aa_order = sorted(df_filtered['AminoAcid'].unique())
        
        sns.boxplot(x='AminoAcid', y='RSCU', data=df_filtered, order=aa_order, palette="Set3")
        sns.stripplot(x='AminoAcid', y='RSCU', data=df_filtered, order=aa_order, color=".25", size=3, alpha=0.5)
        
        plt.xlabel('Amino Acids', fontsize=12)
        plt.ylabel('RSCU Value', fontsize=12)
        plt.title('RSCU Distribution by Amino Acid (All Species)', fontsize=16, fontweight='bold')
        plt.grid(axis='y', alpha=0.3)
        plt.axhline(y=1.0, color='red', linestyle='--', alpha=0.5, label='RSCU = 1.0')
        plt.legend()
        plt.tight_layout()
        
        box_path = os.path.join(output_folder, 'rscu_boxplot_by_aminoacid.png')
        plt.savefig(box_path, dpi=150, bbox_inches="tight")
        plt.close()
        print(f"  ✅ Box plot saved in: {box_path}")
        status_queue.put(("image_ready", (box_path, "RSCU Box Plot by AA")))
    except Exception as e:
        print(f"\n❌ ERROR GENERATING BOX PLOT: {e}")

    print("\n  Generating Chart 2: Stacked Vertical Bars by Species...")
    status_queue.put(("progress", 85))
    aa_codon_map = AA_CODON_MAPS.get(genetic_code_id, AA_CODON_MAPS[1])
    aa_list = sorted([aa for aa, codons in aa_codon_map.items() if aa != '*' and len(codons) > 1]) 
    
    for species_name, rscu_series in df_rscu_matrix.iterrows():
        try:
            fig, ax = plt.subplots(figsize=(16, 10))
            plot_data = []
            for aa in aa_list:
                syn_codons = sorted(aa_codon_map[aa])
                for codon in syn_codons:
                    value = rscu_series.get(codon, 0)
                    plot_data.append({'AminoAcid': aa, 'Codon': codon, 'RSCU': value})
            
            df_plot = pd.DataFrame(plot_data)

            df_plot.pivot(index='AminoAcid', columns='Codon', values='RSCU').loc[aa_list].plot(
                kind='bar', stacked=True, ax=ax,
                colormap='tab20',
                edgecolor='black', linewidth=0.5
            )

            ax.set_xlabel("Amino Acid", fontsize=12, fontweight='bold')
            ax.set_ylabel("RSCU Value (Stacked)", fontsize=12)
            ax.set_title(f"Codon Usage Profile (RSCU) - {species_name}", fontsize=16, fontweight='bold')
            ax.legend(title='Codons',loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=8, fontsize='small')
            ax.grid(axis='y', linestyle='--', alpha=0.7)
            plt.xticks(rotation=0)
            plt.tight_layout() 
            
            hist_v_path = os.path.join(output_folder, f'rscu_vertical_histogram_{species_name}.png')
            plt.savefig(hist_v_path, dpi=150, bbox_inches="tight")
            plt.close()
            print(f"  ✅ Vertical bar chart saved for {species_name}")
            status_queue.put(("image_ready", (hist_v_path, f"Vertical RSCU - {species_name}")))
        except Exception as e:
            print(f"\n❌ ERROR GENERATING VERTICAL HISTOGRAM for {species_name}: {e}")

def enc_gc3_analysis(all_data, output_folder, status_queue):
    print(f"\n=== ENC vs GC3 ANALYSIS ===")
    status_queue.put(("progress", 60))
    
    plot_data = []
    for species, data in all_data.items():
        plot_data.append({'species': species, 'ENC': data['enc'], 'GC3': data['gc3']})
    df_plot = pd.DataFrame(plot_data)

    plt.figure(figsize=(12, 8))
    s_values = np.linspace(0.01, 0.99, 200)
    enc_expected = []
    for s in s_values:
        f_s = s**2 + (1-s)**2
        enc_val = 2 + (9/f_s) + (1/f_s) + (5/f_s) + (3/f_s)
        enc_expected.append(enc_val)
        
    plt.plot(s_values * 100, enc_expected, 'r--', label='Expected Curve (Mutational Bias)', linewidth=2)
    
    scatter = plt.scatter(df_plot['GC3'], df_plot['ENC'], c=df_plot.index, 
                         cmap='viridis', s=100, alpha=0.7, edgecolors='black')
    
    for i, row in df_plot.iterrows():
        plt.annotate(row['species'], (row['GC3'], row['ENC']), 
                    xytext=(5, 5), textcoords='offset points', fontsize=9)
    
    plt.xlabel('GC3 (%)', fontsize=12)
    plt.ylabel('ENC (Effective Number of Codons)', fontsize=12)
    plt.title('ENC vs GC3 Analysis (Wright Plot)', fontsize=16)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.xlim(0, 100)
    plt.ylim(20, 61)
    plt.tight_layout()
    
    output_path = os.path.join(output_folder, 'enc_gc3_analysis.png')
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"\n✅ ENC vs GC3 analysis saved in: {output_path}")
    status_queue.put(("image_ready", (output_path, "ENC vs GC3 Analysis")))
    df_plot.to_csv(os.path.join(output_folder, 'enc_gc3_results.csv'), sep=';', index=False)

def optimal_rare_codons_analysis(all_data, output_folder, status_queue):
    print(f"\n=== OPTIMAL AND RARE CODONS ANALYSIS ===")
    status_queue.put(("progress", 60))
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    species = list(all_data.keys())
    
    optimal_counts = [len(all_data[s]['optimal']) for s in species]
    rare_counts = [len(all_data[s]['rare']) for s in species]
    
    x = np.arange(len(species))
    width = 0.35
    
    ax1.bar(x - width/2, optimal_counts, width, label='Optimal (>1.2)', color='green', alpha=0.7)
    ax1.bar(x + width/2, rare_counts, width, label='Rare (<0.8)', color='red', alpha=0.7)
    
    ax1.set_title('Optimal vs Rare Codons by Species', fontweight='bold')
    ax1.set_ylabel('Number of Codons')
    ax1.set_xticks(x)
    ax1.set_xticklabels(species, rotation=45, ha="right")
    ax1.legend()
    ax1.grid(axis='y', alpha=0.3)
    
    cai_values = [all_data[s]['cai'] for s in species]
    bars = ax2.bar(species, cai_values, color='purple', alpha=0.7)
    ax2.set_title('Codon Adaptation Index (CAI)', fontweight='bold')
    ax2.set_ylabel('CAI')
    ax2.grid(axis='y', alpha=0.3)
    plt.setp(ax2.get_xticklabels(), rotation=45, ha="right")
    
    for bar, value in zip(bars, cai_values):
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01, 
                f'{value:.3f}', ha='center', va='bottom')
    
    plt.suptitle('Optimal, Rare Codons and CAI Analysis', fontsize=16, fontweight='bold')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    
    output_path = os.path.join(output_folder, 'optimal_rare_codons.png')
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"\n✅ Optimal and rare codons analysis saved in: {output_path}")
    status_queue.put(("image_ready", (output_path, "Optimal and Rare Codons")))
    
    for species_name, data in all_data.items():
        df_optimal = pd.DataFrame(list(data['optimal'].items()), 
                                 columns=['Amino_Acid', 'Optimal_Codon'])
        df_rare = pd.DataFrame(list(data['rare'].items()), 
                              columns=['Amino_Acid', 'Rare_Codon'])
        
        df_optimal.to_csv(os.path.join(output_folder, f'{species_name}_optimal_codons.csv'), 
                         sep=';', index=False)
        df_rare.to_csv(os.path.join(output_folder, f'{species_name}_rare_codons.csv'), 
                      sep=';', index=False)

def neutrality_plot_analysis(all_data, output_folder, status_queue):
    print(f"\n=== NEUTRALITY PLOT ANALYSIS ===")
    status_queue.put(("progress", 60))
    
    plot_data = []
    for species, data in all_data.items():
        if 'gc12' in data and 'gc3' in data:
            plot_data.append({
                'species': species, 
                'GC12': data['gc12'], 
                'GC3': data['gc3']
            })
    
    if not plot_data:
        print("  ❌ Error: GC12 or GC3 data not found.")
        return
        
    df_plot = pd.DataFrame(plot_data)
    df_plot = df_plot.replace([np.inf, -np.inf], np.nan).dropna()
    
    if df_plot.empty or len(df_plot) < 2:
        print("  ❌ Error: Insufficient data for linear regression.")
        return

    slope, intercept, r_value, p_value, std_err = linregress(df_plot['GC3'], df_plot['GC12'])
    
    print(f"  📊 Regression: Slope = {slope:.3f}")
    print(f"  📊 Regression: R² = {r_value**2:.3f}")

    plt.figure(figsize=(12, 8))
    sns.scatterplot(x='GC3', y='GC12', data=df_plot, s=100, alpha=0.7, hue='species', legend=False)
    
    x_vals = np.array(plt.xlim())
    y_vals = intercept + slope * x_vals
    plt.plot(x_vals, y_vals, 'r--', label=f'Regression (Slope = {slope:.3f})')
    
    for i, row in df_plot.iterrows():
        plt.annotate(row['species'], (row['GC3'], row['GC12']), 
                    xytext=(5, 5), textcoords='offset points', fontsize=9)
    
    plt.xlabel('GC3 (%)', fontsize=12)
    plt.ylabel('GC12 (%)', fontsize=12)
    plt.title(f"Neutrality Plot (GC12 vs GC3)\nSlope = {slope:.3f}, R² = {r_value**2:.3f}", fontsize=16)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.xlim(0, 100)
    plt.ylim(0, 100)
    plt.tight_layout()
    
    output_path = os.path.join(output_folder, 'neutrality_plot_gc12_vs_gc3.png')
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"\n✅ Neutrality Plot saved in: {output_path}")
    status_queue.put(("image_ready", (output_path, "Neutrality Plot")))
    df_plot.to_csv(os.path.join(output_folder, 'neutrality_plot_results.csv'), sep=';', index=False)