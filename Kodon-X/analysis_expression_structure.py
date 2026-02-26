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

def initiation_mfe_analysis(file_list, output_folder, genetic_code_id, status_queue, gene_list=None, mfe_region_length=50):
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
    
    print("  Generating MFE Box Plot...")
    status_queue.put(("progress", 90))
    plt.figure(figsize=(max(10, len(all_seqs_by_species)*1.5), 8))
    
    sns.boxplot(x='species', y='mfe', data=df_plot_long, palette="coolwarm")
    sns.stripplot(x='species', y='mfe', data=df_plot_long, color=".25", size=2, alpha=0.2)
    
    plt.title(f"MFE Distribution (Minimum Free Energy) - 5' Region ({mfe_region_length}bp)", fontsize=16)
    plt.ylabel('MFE (kcal/mol) - (More negative = more stable/closed)')
    plt.xlabel('Species')
    plt.xticks(rotation=45, ha='right')
    plt.grid(axis='y', linestyle='--', alpha=0.5)
    plt.tight_layout()
    
    output_path = os.path.join(output_folder, 'mfe_5prime_comparative_boxplot.png')
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"\n✅ MFE Chart saved in: {output_path}")
    status_queue.put(("image_ready", (output_path, "MFE Analysis (5' UTR)")))

def two_groups_comparative_analysis(file_list, output_folder, genetic_code_id, status_queue, gene_list_1, gene_list_2):
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
        chi2, p, dof, expected = chi2_contingency(table)
        print(f"  Chi-Square (Total Counts): X²={chi2:.2f}, p-value = {p:.4e}, DoF={dof}")
    except ValueError as e:
        print(f"  Chi-Square Error: {e} (probably low counts)")
        p_values['chi2'] = np.nan
        
    print("\n  Generating comparative box plots...")
    status_queue.put(("progress", 80))
    
    df_g1['group'] = 'Group 1'
    df_g2['group'] = 'Group 2'
    df_plot = pd.concat([df_g1, df_g2])
    
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(20, 8))
    
    sns.boxplot(x='group', y='enc', data=df_plot, ax=ax1, palette="Set2")
    sns.stripplot(x='group', y='enc', data=df_plot, ax=ax1, color=".25", size=3, alpha=0.3)
    ax1.set_title(f"ENC (p = {p_values.get('enc', 'N/A'):.2e})", fontsize=14)
    ax1.set_xlabel("")
    
    sns.boxplot(x='group', y='gc3', data=df_plot, ax=ax2, palette="Set2")
    sns.stripplot(x='group', y='gc3', data=df_plot, ax=ax2, color=".25", size=3, alpha=0.3)
    ax2.set_title(f"GC3 (p = {p_values.get('gc3', 'N/A'):.2e})", fontsize=14)
    ax2.set_xlabel("")
    
    sns.boxplot(x='group', y='cai', data=df_plot, ax=ax3, palette="Set2")
    sns.stripplot(x='group', y='cai', data=df_plot, ax=ax3, color=".25", size=3, alpha=0.3)
    ax3.set_title(f"CAI (p = {p_values.get('cai', 'N/A'):.2e})", fontsize=14)
    ax3.set_xlabel("")
    
    plt.suptitle(f"Group Comparison (G1: {len(df_g1)} genes, G2: {len(df_g2)} genes)", fontsize=18)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    
    output_path = os.path.join(output_folder, 'group_comparison_boxplot.png')
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"\n✅ Group Comparison Chart saved in: {output_path}")
    status_queue.put(("image_ready", (output_path, "Gene Group Comparison")))

def expression_correlation_analysis(file_list, output_folder, genetic_code_id, status_queue, gene_list, expression_data, gene_col, expr_col):
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
    
    print("\n  Generating correlation plots...")
    status_queue.put(("progress", 80))
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 8))
    
    sns.regplot(x='expr_log', y='cai', data=df_merged, ax=ax1,
                scatter_kws={'alpha': 0.2, 's': 10}, 
                line_kws={'color': 'red'})
    ax1.set_title(f"CAI vs Expression (log10)\nSpearman Rho = {corr_cai:.3f} (p = {p_cai:.2e})", fontsize=14)
    ax1.set_xlabel(f"Expression (log10 {expr_col})")
    ax1.set_ylabel("CAI")
    
    sns.regplot(x='expr_log', y='enc', data=df_merged, ax=ax2,
                scatter_kws={'alpha': 0.2, 's': 10}, 
                line_kws={'color': 'blue'})
    ax2.set_title(f"ENC vs Expression (log10)\nSpearman Rho = {corr_enc:.3f} (p = {p_enc:.2e})", fontsize=14)
    ax2.set_xlabel(f"Expression (log10 {expr_col})")
    ax2.set_ylabel("ENC")
    
    plt.suptitle(f"Correlation with Expression (N = {len(df_merged)} genes)", fontsize=18)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    
    output_path = os.path.join(output_folder, 'expression_correlation.png')
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"\n✅ Expression Correlation Chart saved in: {output_path}")
    status_queue.put(("image_ready", (output_path, "Correlation with Expression")))