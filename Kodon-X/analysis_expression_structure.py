import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu, chi2_contingency, spearmanr
from Bio import SeqIO
from collections import Counter

try:
    import RNA
    VIENNARNA_AVAILABLE = True
except ImportError:
    VIENNARNA_AVAILABLE = False

from core_utils import extract_cds_sequences, process_genomes_for_bias_analysis, get_w_reference_table, calculate_metrics_per_gene

def _calculate_mfe_for_sequences(sequences, mfe_region_length):
    if not VIENNARNA_AVAILABLE:
        raise ImportError("ViennaRNA (import RNA) not found. Cannot calculate MFE.")
        
    mfe_results = []
    
    for seq in sequences:
        if len(seq) < mfe_region_length:
            continue
            
        region_5prime = seq[:mfe_region_length]
        
        try:
            (structure, mfe) = RNA.fold(region_5prime)
            mfe_results.append({'mfe': mfe, 'length': len(region_5prime)})
        except Exception as e:
            print(f"  Warning: Failed to calculate MFE for sequence {region_5prime[:10]}...: {e}")
            
    return pd.DataFrame(mfe_results)

def initiation_mfe_analysis(file_list, output_folder, genetic_code_id, status_queue, gene_list=None, mfe_region_length=50, palette='viridis'):
    print(f"\n=== INITIATION MFE ANALYSIS (5' UTR) ===")
    print(f"  Parameters: Region = First {mfe_region_length}bp")
    
    if not VIENNARNA_AVAILABLE:
        print("  ❌ CRITICAL ERROR: ViennaRNA (import RNA) not found.")
        print("     Analysis 17 cannot be executed. Install 'viennarna' via Conda.")
        status_queue.put(("message", "Error: ViennaRNA not found."))
        return

    all_seqs_by_species = extract_cds_sequences(file_list, status_queue, gene_list)
    if not all_seqs_by_species:
        print("  ❌ Error: No valid CDS sequence was extracted.")
        return

    status_queue.put(("progress", 40))
    all_results_long = []

    for species_name, sequences in all_seqs_by_species.items():
        print(f"  Calculating MFE for {species_name}...")
        status_queue.put(("message", f"Calculating MFE for {species_name}..."))
        
        try:
            df_results = _calculate_mfe_for_sequences(sequences, mfe_region_length)
        except ImportError as e:
            print(f"  ❌ Error: {e}")
            return
            
        if df_results.empty:
            print(f"  Warning: No MFE data calculated for {species_name}.")
            continue
            
        df_results['species'] = species_name
        all_results_long.append(df_results)
        
    if not all_results_long:
        print("  ❌ Error: No MFE data could be calculated.")
        return
        
    df_plot_long = pd.concat(all_results_long).dropna()
    df_plot_long.to_csv(os.path.join(output_folder, 'mfe_data_per_gene.csv'), sep=';', index=False)
    
    status_queue.put(("progress", 80))

    try:
        print("  Generating Chart 1: MFE Box Plot...")
        plt.figure(figsize=(max(12, len(all_seqs_by_species)*1.5), 8))
        sns.boxplot(x='species', y='mfe', data=df_plot_long, palette=palette, showfliers=False)
        plt.title(f"1. MFE Distribution (Minimum Free Energy) - 5' Region ({mfe_region_length}bp)", fontsize=16)
        plt.ylabel('MFE (kcal/mol)\n(More negative = more stable/structured)')
        plt.xlabel('Species')
        plt.xticks(rotation=45, ha='right')
        plt.grid(axis='y', linestyle='--', alpha=0.5)
        plt.tight_layout()
        
        plot_path1 = os.path.join(output_folder, 'mfe_1_comparative_boxplot.png')
        plt.savefig(plot_path1, dpi=150, bbox_inches='tight')
        plt.close()
        status_queue.put(("image_ready", (plot_path1, "1. MFE Boxplot (5' Region)")))
    except Exception as e:
        print(f"  ❌ Error generating MFE Boxplot: {e}")

    try:
        print("  Generating Chart 2: MFE KDE Density...")
        plt.figure(figsize=(12, 8))
        sns.kdeplot(x='mfe', hue='species', data=df_plot_long, palette=palette, fill=True, common_norm=False, alpha=0.4)
        plt.title(f"2. Density Distribution of MFE Values ({mfe_region_length}bp)", fontsize=16)
        plt.xlabel('Minimum Free Energy - MFE (kcal/mol)')
        plt.ylabel('Density')
        plt.grid(alpha=0.3)
        plt.tight_layout()
        
        plot_path2 = os.path.join(output_folder, 'mfe_2_density_kde.png')
        plt.savefig(plot_path2, dpi=150, bbox_inches='tight')
        plt.close()
        status_queue.put(("image_ready", (plot_path2, "2. MFE Density (KDE)")))
    except Exception as e:
        print(f"  ❌ Error generating MFE KDE: {e}")

    try:
        print("  Generating Chart 3: MFE Violin Plot...")
        plt.figure(figsize=(max(12, len(all_seqs_by_species)*1.5), 8))
        sns.violinplot(x='species', y='mfe', data=df_plot_long, palette=palette, inner='quartile')
        plt.title(f"3. Violin Plot of MFE Distributions ({mfe_region_length}bp)", fontsize=16)
        plt.xlabel('Species')
        plt.ylabel('MFE (kcal/mol)')
        plt.xticks(rotation=45, ha='right')
        plt.grid(axis='y', alpha=0.3)
        plt.tight_layout()
        
        plot_path3 = os.path.join(output_folder, 'mfe_3_violinplot.png')
        plt.savefig(plot_path3, dpi=150, bbox_inches='tight')
        plt.close()
        status_queue.put(("image_ready", (plot_path3, "3. MFE Violin Plot")))
    except Exception as e:
        print(f"  ❌ Error generating MFE Violin Plot: {e}")

def two_groups_comparative_analysis(file_list, output_folder, genetic_code_id, status_queue, gene_list_1, gene_list_2, palette='viridis'):
    print(f"\n=== COMPARATIVE ANALYSIS OF TWO GROUPS ===")
    print(f"  Group 1: {len(gene_list_1)} genes")
    print(f"  Group 2: {len(gene_list_2)} genes")
    
    if not gene_list_1 or not gene_list_2:
        print("  ❌ Error: Both gene groups must be provided.")
        return

    print("  Calculating reference table W (full genome)...")
    status_queue.put(("progress", 10))
    status_queue.put(("message", "Calculating reference (genome)..."))
    all_bias_data_full_genome = process_genomes_for_bias_analysis(file_list, genetic_code_id, status_queue, gene_list=None)
    
    if not all_bias_data_full_genome:
        print("  ❌ Error: Failed to process reference genomes.")
        return
        
    all_results_g1 = []
    all_results_g2 = []
    all_counts_g1 = Counter()
    all_counts_g2 = Counter()
    
    status_queue.put(("progress", 40))
    
    for i, full_path in enumerate(file_list):
        base_name = os.path.basename(full_path).split('.')[0]
        print(f"\n  Processing {base_name}...")
        
        if base_name not in all_bias_data_full_genome:
            print(f"  Warning: {base_name} not found in reference data. Skipping.")
            continue
            
        w_reference_table = get_w_reference_table(
            all_bias_data_full_genome[base_name]['rscu'], 
            genetic_code_id
        )
        
        status_queue.put(("message", f"Calculating Group 1 in {base_name}..."))
        results_g1 = calculate_metrics_per_gene(full_path, w_reference_table, genetic_code_id, gene_list_1)
        for res in results_g1:
            all_results_g1.append(res)
            all_counts_g1.update(res['counts'])
            
        status_queue.put(("message", f"Calculating Group 2 in {base_name}..."))
        results_g2 = calculate_metrics_per_gene(full_path, w_reference_table, genetic_code_id, gene_list_2)
        for res in results_g2:
            all_results_g2.append(res)
            all_counts_g2.update(res['counts'])

    if not all_results_g1 or not all_results_g2:
        print("  ❌ Error: No matching gene found in one or both groups.")
        return
        
    df_g1 = pd.DataFrame(all_results_g1)
    df_g2 = pd.DataFrame(all_results_g2)
    
    print("\n--- Statistical Results (Mann-Whitney U) ---")
    metrics_to_test = ['enc', 'gc3', 'cai']
    p_values = {}
    for metric in metrics_to_test:
        data1 = df_g1[metric].dropna()
        data2 = df_g2[metric].dropna()
        if len(data1) > 0 and len(data2) > 0:
            stat, p = mannwhitneyu(data1, data2, alternative='two-sided')
            p_values[metric] = p
            print(f"  {metric.upper()}: p-value = {p:.4e} (Median G1: {data1.median():.3f}, Median G2: {data2.median():.3f})")
        else:
            p_values[metric] = np.nan
            print(f"  {metric.upper()}: Insufficient data for test.")

    print("\n--- Statistical Results (Chi-Square on Codon Counts) ---")
    codons = sorted(list(set(all_counts_g1.keys()) | set(all_counts_g2.keys())))
    table = [
        [all_counts_g1.get(c, 0) for c in codons],
        [all_counts_g2.get(c, 0) for c in codons]
    ]
    
    try:
        chi2, p_chi, dof, expected = chi2_contingency(table)
        print(f"  Chi-Square (Total Counts): X²={chi2:.2f}, p-value = {p_chi:.4e}, DoF={dof}")
    except ValueError as e:
        print(f"  Chi-Square Error: {e} (probably low counts)")
        p_chi = np.nan
        
    status_queue.put(("progress", 80))
    df_g1['group'] = 'Group 1'
    df_g2['group'] = 'Group 2'
    df_plot = pd.concat([df_g1, df_g2])

    try:
        print("  Generating Chart 1: Group Boxplots...")
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(20, 8))
        
        sns.boxplot(x='group', y='enc', data=df_plot, ax=ax1, palette=palette, showfliers=False)
        ax1.set_title(f"ENC (p = {p_values.get('enc', np.nan):.2e})", fontsize=14)
        ax1.set_xlabel("")
        
        sns.boxplot(x='group', y='gc3', data=df_plot, ax=ax2, palette=palette, showfliers=False)
        ax2.set_title(f"GC3 (p = {p_values.get('gc3', np.nan):.2e})", fontsize=14)
        ax2.set_xlabel("")
        
        sns.boxplot(x='group', y='cai', data=df_plot, ax=ax3, palette=palette, showfliers=False)
        ax3.set_title(f"CAI (p = {p_values.get('cai', np.nan):.2e})", fontsize=14)
        ax3.set_xlabel("")
        
        plt.suptitle(f"1. Group Comparison Metrics (G1: {len(df_g1)} genes, G2: {len(df_g2)} genes)", fontsize=18)
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        
        plot_path1 = os.path.join(output_folder, 'group_1_boxplots.png')
        plt.savefig(plot_path1, dpi=150, bbox_inches='tight')
        plt.close()
        status_queue.put(("image_ready", (plot_path1, "1. Group Metrics Boxplots")))
    except Exception as e:
        print(f"  ❌ Error: {e}")

    try:
        print("  Generating Chart 2: KDE Densities...")
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(20, 6))
        
        sns.kdeplot(x='enc', hue='group', data=df_plot, ax=ax1, palette=palette, fill=True, common_norm=False)
        ax1.set_title("ENC Density", fontsize=14)
        
        sns.kdeplot(x='gc3', hue='group', data=df_plot, ax=ax2, palette=palette, fill=True, common_norm=False)
        ax2.set_title("GC3 Density", fontsize=14)
        
        sns.kdeplot(x='cai', hue='group', data=df_plot, ax=ax3, palette=palette, fill=True, common_norm=False)
        ax3.set_title("CAI Density", fontsize=14)
        
        plt.suptitle("2. Density Distributions (KDE) - Group 1 vs Group 2", fontsize=18)
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        
        plot_path2 = os.path.join(output_folder, 'group_2_kde_densities.png')
        plt.savefig(plot_path2, dpi=150, bbox_inches='tight')
        plt.close()
        status_queue.put(("image_ready", (plot_path2, "2. Group Metrics KDE Densities")))
    except Exception as e:
        print(f"  ❌ Error: {e}")

    try:
        print("  Generating Chart 3: GC3 vs ENC Scatter...")
        plt.figure(figsize=(10, 8))
        sns.scatterplot(x='gc3', y='enc', hue='group', data=df_plot, palette=palette, alpha=0.7, s=50, edgecolor='k')
        plt.title("3. GC3 vs ENC Space by Gene Group", fontsize=16)
        plt.xlabel("GC3 (%)")
        plt.ylabel("Effective Number of Codons (ENC)")
        plt.grid(alpha=0.3)
        plt.legend(title='Groups')
        plt.tight_layout()
        
        plot_path3 = os.path.join(output_folder, 'group_3_gc3_vs_enc_scatter.png')
        plt.savefig(plot_path3, dpi=150, bbox_inches='tight')
        plt.close()
        status_queue.put(("image_ready", (plot_path3, "3. GC3 vs ENC Scatter by Group")))
    except Exception as e:
        print(f"  ❌ Error: {e}")

    try:
        print("  Generating Chart 4: Group Means Barplot...")
        df_melt_metrics = df_plot.melt(id_vars='group', value_vars=['enc', 'gc3', 'cai'], 
                                       var_name='Metric', value_name='Value')
        
        fig, axes = plt.subplots(1, 3, figsize=(18, 6))
        for i, metric in enumerate(['enc', 'gc3', 'cai']):
            sns.barplot(x='group', y='Value', data=df_melt_metrics[df_melt_metrics['Metric'] == metric], 
                        ax=axes[i], capsize=.1, errorbar='sd', palette=palette)
            axes[i].set_title(f"Average {metric.upper()}", fontsize=14)
            axes[i].set_xlabel("")
            axes[i].grid(axis='y', alpha=0.3)
            
        plt.suptitle("4. Average Values (with Standard Deviation) per Metric", fontsize=18)
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        
        plot_path4 = os.path.join(output_folder, 'group_4_means_barplot.png')
        plt.savefig(plot_path4, dpi=150, bbox_inches='tight')
        plt.close()
        status_queue.put(("image_ready", (plot_path4, "4. Group Means Barplot")))
    except Exception as e:
        print(f"  ❌ Error: {e}")

def expression_correlation_analysis(file_list, output_folder, genetic_code_id, status_queue, gene_list, expression_data, gene_col, expr_col, palette='viridis'):
    print(f"\n=== CORRELATION WITH EXPRESSION ANALYSIS ===")
    print(f"  Expression File: {len(expression_data)} rows")
    print(f"  Gene Column: '{gene_col}', Expression Column: '{expr_col}'")

    try:
        expression_data[expr_col] = pd.to_numeric(expression_data[expr_col])
    except ValueError:
        print(f"  ❌ Error: The expression column '{expr_col}' is not numeric.")
        return
        
    min_expr = expression_data[expression_data[expr_col] > 0][expr_col].min()
    if pd.isna(min_expr): min_expr = 1e-3
    expression_data['expr_log'] = np.log10(expression_data[expr_col] + min_expr / 10)
    
    print("  Calculating reference table W (full genome)...")
    status_queue.put(("progress", 10))
    status_queue.put(("message", "Calculating reference (genome)..."))
    all_bias_data_full_genome = process_genomes_for_bias_analysis(file_list, genetic_code_id, status_queue, gene_list=None)
    
    if not all_bias_data_full_genome:
        print("  ❌ Error: Failed to process reference genomes.")
        return
        
    all_merged_data = []
    
    status_queue.put(("progress", 40))
    all_genes_in_gbk = set(gene_list) if gene_list else set()
    if not gene_list:
        print("  Collecting all genes from GenBank files...")
        for full_path in file_list:
            try:
                for record in SeqIO.parse(full_path, "genbank"):
                    for feature in record.features:
                        if feature.type != "CDS": continue
                        if "locus_tag" in feature.qualifiers: 
                            all_genes_in_gbk.add(feature.qualifiers["locus_tag"][0])
                        if "gene" in feature.qualifiers:
                            all_genes_in_gbk.add(feature.qualifiers["gene"][0])
            except Exception as e:
                print(f"  Warning: Error reading {full_path} for gene list: {e}")
    
    print(f"  Total of {len(all_genes_in_gbk)} genes for analysis.")

    for i, full_path in enumerate(file_list):
        base_name = os.path.basename(full_path).split('.')[0]
        print(f"\n  Processing {base_name}...")
        
        if base_name not in all_bias_data_full_genome:
            continue
            
        w_reference_table = get_w_reference_table(
            all_bias_data_full_genome[base_name]['rscu'], 
            genetic_code_id
        )
        
        status_queue.put(("message", f"Calculating metrics in {base_name}..."))
        results_genes = calculate_metrics_per_gene(
            full_path, w_reference_table, genetic_code_id, all_genes_in_gbk
        )
        
        if results_genes:
            df_genes = pd.DataFrame(results_genes)
            df_genes['species'] = base_name
            all_merged_data.append(df_genes)

    if not all_merged_data:
        print("  ❌ Error: No gene with metrics was calculated.")
        return

    df_plot_all_genes = pd.concat(all_merged_data)
    
    print("  Merging metrics and expression data...")
    df_merged = pd.merge(
        df_plot_all_genes, 
        expression_data, 
        left_on='gene', 
        right_on=gene_col,
        how='inner' 
    )
    
    if df_merged.empty:
        print(f"  ❌ Error: No common genes found between GenBank files and the expression file.")
        print(f"     Verify if identifiers in column '{gene_col}' match GenBank 'locus_tag' or 'gene'.")
        return
        
    print(f"  ✅ {len(df_merged)} genes with metric and expression data found.")
    
    print("\n--- Statistical Results (Spearman Correlation) ---")
    corr_cai, p_cai = spearmanr(df_merged['cai'].dropna(), df_merged['expr_log'].dropna())
    print(f"  CAI vs Expression (log10): Rho = {corr_cai:.3f}, p-value = {p_cai:.4e}")
    
    corr_enc, p_enc = spearmanr(df_merged['enc'].dropna(), df_merged['expr_log'].dropna())
    print(f"  ENC vs Expression (log10): Rho = {corr_enc:.3f}, p-value = {p_enc:.4e}")
    
    status_queue.put(("progress", 80))

    try:
        plot_col1 = sns.color_palette(palette)[-1]
        plot_col2 = sns.color_palette(palette)[0]
    except:
        plot_col1 = 'red'
        plot_col2 = 'blue'

    try:
        print("  Generating Chart 1: CAI vs Expression Scatter...")
        plt.figure(figsize=(10, 8))
        sns.regplot(x='expr_log', y='cai', data=df_merged, 
                    scatter_kws={'alpha': 0.3, 's': 30, 'edgecolor': 'w'}, 
                    line_kws={'color': plot_col1, 'linewidth': 2})
        plt.title(f"1. CAI vs Gene Expression (log10)\nSpearman Rho = {corr_cai:.3f} (p = {p_cai:.2e})", fontsize=16)
        plt.xlabel(f"Expression (log10 {expr_col})")
        plt.ylabel("Codon Adaptation Index (CAI)")
        plt.grid(alpha=0.3)
        plt.tight_layout()
        
        plot_path1 = os.path.join(output_folder, 'expr_1_cai_scatter_regression.png')
        plt.savefig(plot_path1, dpi=150, bbox_inches='tight')
        plt.close()
        status_queue.put(("image_ready", (plot_path1, "1. CAI vs Expression Regression")))
    except Exception as e:
        print(f"  ❌ Error: {e}")

    try:
        print("  Generating Chart 2: ENC vs Expression Scatter...")
        plt.figure(figsize=(10, 8))
        sns.regplot(x='expr_log', y='enc', data=df_merged, 
                    scatter_kws={'alpha': 0.3, 's': 30, 'edgecolor': 'w'}, 
                    line_kws={'color': plot_col2, 'linewidth': 2})
        plt.title(f"2. ENC vs Gene Expression (log10)\nSpearman Rho = {corr_enc:.3f} (p = {p_enc:.2e})", fontsize=16)
        plt.xlabel(f"Expression (log10 {expr_col})")
        plt.ylabel("Effective Number of Codons (ENC)")
        plt.grid(alpha=0.3)
        plt.tight_layout()
        
        plot_path2 = os.path.join(output_folder, 'expr_2_enc_scatter_regression.png')
        plt.savefig(plot_path2, dpi=150, bbox_inches='tight')
        plt.close()
        status_queue.put(("image_ready", (plot_path2, "2. ENC vs Expression Regression")))
    except Exception as e:
        print(f"  ❌ Error: {e}")

    try:
        print("  Generating Chart 3: CAI vs Expression Hexbin Density...")
        plt.figure(figsize=(10, 8))
        plt.hexbin(df_merged['expr_log'], df_merged['cai'], gridsize=40, cmap=palette, mincnt=1)
        cb = plt.colorbar(label='Number of Genes (Density)')
        plt.title("3. 2D Density Map (Hexbin) of CAI vs Expression", fontsize=16)
        plt.xlabel(f"Expression (log10 {expr_col})")
        plt.ylabel("Codon Adaptation Index (CAI)")
        plt.grid(alpha=0.2)
        plt.tight_layout()
        
        plot_path3 = os.path.join(output_folder, 'expr_3_cai_hexbin_density.png')
        plt.savefig(plot_path3, dpi=150, bbox_inches='tight')
        plt.close()
        status_queue.put(("image_ready", (plot_path3, "3. Expression Hexbin Density")))
    except Exception as e:
        print(f"  ❌ Error: {e}")

    try:
        print("  Generating Chart 4: Expression by CAI Quartiles...")
        df_merged['CAI_Quartile'] = pd.qcut(df_merged['cai'], 4, labels=['Q1 (Lowest Bias)', 'Q2', 'Q3', 'Q4 (Highest Bias)'])
        
        plt.figure(figsize=(12, 8))
        sns.boxplot(x='CAI_Quartile', y='expr_log', data=df_merged, palette=palette, showfliers=False)
        plt.title("4. Gene Expression Levels Distributed by CAI Quartiles", fontsize=16)
        plt.xlabel("CAI Quartiles")
        plt.ylabel(f"Expression (log10 {expr_col})")
        plt.grid(axis='y', alpha=0.3)
        plt.tight_layout()
        
        plot_path4 = os.path.join(output_folder, 'expr_4_quartiles_boxplot.png')
        plt.savefig(plot_path4, dpi=150, bbox_inches='tight')
        plt.close()
        status_queue.put(("image_ready", (plot_path4, "4. Expression by CAI Quartiles")))
    except Exception as e:
        print(f"  ❌ Error: {e}")
