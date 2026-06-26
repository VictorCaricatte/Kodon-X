import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import re
from collections import Counter
from itertools import product
from Bio import SeqIO
from Bio.Seq import Seq
from constants import (
    GENETIC_CODE_TABLES, ALL_CODONS_SORTED, AA_CODON_MAPS, 
    KYTE_DOOLITTLE_HYDROPATHY, AROMATICITY, WOBBLE_MATRIX_BACTERIA, 
    WOBBLE_MATRIX_EUKARYA, ANTICODON_MODIFICATION_MAP
)
from core_utils import extract_cds_sequences, _apply_gene_filter

def calculate_codon_pair_bias(sequences, genetic_code_id):
    """Codon Pair Score — Coleman et al. (2008).

    CPS(Ci, Cj) = ln( F_obs(Ci,Cj) / F_exp(Ci,Cj) )

    F_exp(Ci,Cj) = F(Ci) * F(Cj) * F(Aa_i,Aa_j) / ( F(Aa_i) * F(Aa_j) )

    Amino acid pair frequencies use pseudocount 0.5 to avoid zeros.
    Stop codons are excluded from all counts.
    Returns 0 for undefined pairs (F_obs=0 or F_exp=0).
    """
    genetic_code = GENETIC_CODE_TABLES.get(genetic_code_id, GENETIC_CODE_TABLES[1])
    stop_codons = {codon for codon, aa in genetic_code.items() if aa == '*'}

    non_stop_codons = sorted([c for c in ALL_CODONS_SORTED if c not in stop_codons])
    codon_to_aa = {c: genetic_code[c] for c in non_stop_codons}
    non_stop_aas = sorted(set(codon_to_aa.values()))

    codon_index = {c: i for i, c in enumerate(non_stop_codons)}
    aa_index    = {a: i for i, a in enumerate(non_stop_aas)}

    n_c = len(non_stop_codons)
    n_a = len(non_stop_aas)

    pair_counts  = np.zeros((n_c, n_c))   # observed codon pair counts
    cnt_c1       = np.zeros(n_c)           # codon counts in position 1
    cnt_c2       = np.zeros(n_c)           # codon counts in position 2
    aa_pair_cnt  = np.zeros((n_a, n_a))   # observed amino-acid pair counts
    total_pairs  = 0

    for seq in sequences:
        for i in range(0, len(seq) - 5, 3):
            c1 = seq[i:i+3]
            c2 = seq[i+3:i+6]
            if c1 in codon_index and c2 in codon_index:
                i1, i2   = codon_index[c1], codon_index[c2]
                ai1, ai2 = aa_index[codon_to_aa[c1]], aa_index[codon_to_aa[c2]]
                pair_counts[i1, i2]    += 1
                cnt_c1[i1]             += 1
                cnt_c2[i2]             += 1
                aa_pair_cnt[ai1, ai2]  += 1
                total_pairs            += 1

    if total_pairs == 0:
        print("  Warning: No valid codon pair found.")
        return pd.DataFrame()

    # ── Codon frequencies (position-specific) ─────────────────────────────────
    F_c1 = cnt_c1 / total_pairs   # shape (n_c,)
    F_c2 = cnt_c2 / total_pairs   # shape (n_c,)

    # ── Amino-acid pair frequencies with pseudocount 0.5 ──────────────────────
    aa_pair_smooth  = aa_pair_cnt + 0.5
    total_aa_smooth = aa_pair_smooth.sum()
    F_aa_pair = aa_pair_smooth / total_aa_smooth           # shape (n_a, n_a)
    F_aa1     = F_aa_pair.sum(axis=1)                     # marginal P(Aa_i) in pos-1
    F_aa2     = F_aa_pair.sum(axis=0)                     # marginal P(Aa_j) in pos-2

    # ── AA ratio matrix: F(Aa_i,Aa_j) / (F(Aa_i) * F(Aa_j)) ─────────────────
    aa_outer = np.outer(F_aa1, F_aa2)
    aa_ratio = np.where(aa_outer > 0, F_aa_pair / aa_outer, 0.0)  # (n_a, n_a)

    # ── Map codon indices to their AA indices for broadcasting ─────────────────
    c1_aa = np.array([aa_index[codon_to_aa[c]] for c in non_stop_codons])
    c2_aa = np.array([aa_index[codon_to_aa[c]] for c in non_stop_codons])
    aa_ratio_codons = aa_ratio[np.ix_(c1_aa, c2_aa)]     # (n_c, n_c)

    # ── F_exp and F_obs ────────────────────────────────────────────────────────
    F_exp = np.outer(F_c1, F_c2) * aa_ratio_codons        # (n_c, n_c)
    F_obs = pair_counts / total_pairs                      # (n_c, n_c)

    # ── CPS = ln(F_obs / F_exp); 0 where either term is zero ──────────────────
    with np.errstate(divide='ignore', invalid='ignore'):
        cps_matrix = np.where(
            (F_obs > 0) & (F_exp > 0),
            np.log(F_obs / F_exp),
            0.0
        )

    return pd.DataFrame(cps_matrix, index=non_stop_codons, columns=non_stop_codons)

def codon_pair_bias_analysis(file_list, output_folder, genetic_code_id, status_queue, gene_list=None, palette='viridis'):
    print(f"\n=== CODON PAIR SCORE ANALYSIS (CPS — Coleman et al. 2008) ===")
    all_seqs_by_species = extract_cds_sequences(file_list, status_queue, gene_list)

    if not all_seqs_by_species:
        print("  ❌ Error: No valid CDS sequence was extracted.")
        return

    status_queue.put(("progress", 40))

    for i, (species_name, sequences) in enumerate(all_seqs_by_species.items()):
        print(f"  Calculating CPS for {species_name}...")
        status_queue.put(("message", f"Calculating CPS for {species_name}..."))

        df_cps = calculate_codon_pair_bias(sequences, genetic_code_id)

        if df_cps.empty:
            continue

        csv_path = os.path.join(output_folder, f"cps_matrix_{species_name}.csv")
        df_cps.to_csv(csv_path, sep=';', decimal='.')
        print(f"  ✅ CPS matrix saved: {csv_path}")

        try:
            print("  Generating Chart 1: CPS Heatmap...")
            plt.figure(figsize=(24, 20))
            sns.heatmap(df_cps, cmap='RdBu_r', center=0, annot=False,
                        cbar_kws={'label': 'Codon Pair Score — CPS (ln Obs/Exp)'})
            plt.title(f"1. Codon Pair Score Heatmap (CPS) — {species_name}", fontsize=18)
            plt.xlabel("Second Codon (Cj)", fontsize=12)
            plt.ylabel("First Codon (Ci)", fontsize=12)
            plt.tight_layout()

            plot_path1 = os.path.join(output_folder, f"cps_1_heatmap_{species_name}.png")
            plt.savefig(plot_path1, dpi=100, bbox_inches="tight")
            plt.close()
            status_queue.put(("image_ready", (plot_path1, f"1. CPS Heatmap — {species_name}")))
        except Exception as e:
            print(f"  ❌ Error generating CPS Heatmap: {e}")

        df_cps_long = df_cps.unstack().reset_index()
        df_cps_long.columns = ['Codon_1', 'Codon_2', 'CPS']
        df_cps_long['Pair'] = df_cps_long['Codon_1'] + "-" + df_cps_long['Codon_2']
        df_cps_long_nonzero = df_cps_long[df_cps_long['CPS'] != 0.0]

        try:
            plot_col = sns.color_palette(palette)[0]
        except Exception:
            plot_col = 'purple'

        try:
            print("  Generating Chart 2: CPS Distribution Histogram...")
            plt.figure(figsize=(10, 6))
            sns.histplot(df_cps_long_nonzero['CPS'], bins=60, kde=True,
                         color=plot_col, edgecolor='black')
            plt.axvline(0, color='red', linestyle='--', label='Neutral (CPS = 0)')
            plt.title(f"2. Distribution of Codon Pair Scores (CPS) — {species_name}", fontsize=14)
            plt.xlabel("Codon Pair Score (CPS = ln Obs/Exp)")
            plt.ylabel("Number of Codon Pairs")
            plt.legend()
            plt.grid(alpha=0.3)
            plt.tight_layout()

            plot_path2 = os.path.join(output_folder, f"cps_2_histogram_{species_name}.png")
            plt.savefig(plot_path2, dpi=100, bbox_inches="tight")
            plt.close()
            status_queue.put(("image_ready", (plot_path2, f"2. CPS Histogram — {species_name}")))
        except Exception as e:
            print(f"  ❌ Error generating CPS Histogram: {e}")

        try:
            print("  Generating Chart 3: Top 20 Over-represented Pairs...")
            top_over = df_cps_long_nonzero.sort_values('CPS', ascending=False).head(20)
            plt.figure(figsize=(12, 8))
            sns.barplot(x='CPS', y='Pair', data=top_over, hue='Pair', palette=palette, legend=False)
            plt.axvline(0, color='grey', linestyle='--', linewidth=0.8)
            plt.title(f"3. Top 20 Over-represented Codon Pairs (CPS) — {species_name}", fontsize=14)
            plt.xlabel("Codon Pair Score (CPS)")
            plt.ylabel("Codon Pair (Ci–Cj)")
            plt.grid(axis='x', alpha=0.3)
            plt.tight_layout()

            plot_path3 = os.path.join(output_folder, f"cps_3_over_{species_name}.png")
            plt.savefig(plot_path3, dpi=100, bbox_inches="tight")
            plt.close()
            status_queue.put(("image_ready", (plot_path3, f"3. Top Over-represented — {species_name}")))
        except Exception as e:
            print(f"  ❌ Error generating over-represented chart: {e}")

        try:
            print("  Generating Chart 4: Top 20 Under-represented Pairs...")
            top_under = df_cps_long_nonzero.sort_values('CPS', ascending=True).head(20)
            plt.figure(figsize=(12, 8))
            sns.barplot(x='CPS', y='Pair', data=top_under, hue='Pair', palette=palette, legend=False)
            plt.axvline(0, color='grey', linestyle='--', linewidth=0.8)
            plt.title(f"4. Top 20 Under-represented Codon Pairs (CPS) — {species_name}", fontsize=14)
            plt.xlabel("Codon Pair Score (CPS)")
            plt.ylabel("Codon Pair (Ci–Cj)")
            plt.grid(axis='x', alpha=0.3)
            plt.tight_layout()

            plot_path4 = os.path.join(output_folder, f"cps_4_under_{species_name}.png")
            plt.savefig(plot_path4, dpi=100, bbox_inches="tight")
            plt.close()
            status_queue.put(("image_ready", (plot_path4, f"4. Top Under-represented — {species_name}")))
        except Exception as e:
            print(f"  ❌ Error generating under-represented chart: {e}")

        status_queue.put(("progress", int(40 + ((i + 1) / len(all_seqs_by_species)) * 50)))

def calculate_gravy_aromo(sequences, genetic_code_id):
    codon_map = GENETIC_CODE_TABLES.get(genetic_code_id, GENETIC_CODE_TABLES[1])
    results = []
    
    for seq in sequences:
        aa_seq = []
        for i in range(0, len(seq), 3):
            codon = seq[i:i+3]
            aa = codon_map.get(codon, '*')
            if aa == '*':
                break 
            aa_seq.append(aa)
            
        if not aa_seq:
            continue
        
        total_len = len(aa_seq)
        
        try:
            gravy_score = sum(KYTE_DOOLITTLE_HYDROPATHY.get(aa, 0) for aa in aa_seq) / total_len
        except ZeroDivisionError:
            gravy_score = 0
            
        try:
            aromo_score = sum(AROMATICITY.get(aa, 0) for aa in aa_seq) / total_len
        except ZeroDivisionError:
            aromo_score = 0
            
        results.append({
            'gene_length_aa': total_len,
            'gravy': gravy_score,
            'aromo': aromo_score
        })
        
    return pd.DataFrame(results)

def gravy_aromo_analysis(file_list, output_folder, genetic_code_id, status_queue, gene_list=None, palette='viridis'):
    print(f"\n=== PHYSICOCHEMICAL ANALYSIS (GRAVY & AROMO) ===")
    all_seqs_by_species = extract_cds_sequences(file_list, status_queue, gene_list)
    
    if not all_seqs_by_species:
        print("  ❌ Error: No valid CDS sequence was extracted.")
        return

    status_queue.put(("progress", 40))
    all_results_long = []
    
    for species_name, sequences in all_seqs_by_species.items():
        print(f"  Calculating GRAVY/Aromo for {species_name}...")
        status_queue.put(("message", f"Calculating GRAVY/Aromo for {species_name}..."))
        
        df_results = calculate_gravy_aromo(sequences, genetic_code_id)
        if df_results.empty: continue
            
        csv_path = os.path.join(output_folder, f"gravy_aromo_per_gene_{species_name}.csv")
        df_results.to_csv(csv_path, sep=';', decimal='.', index=False)
        print(f"  ✅ GRAVY/Aromo table (per gene) saved in: {csv_path}")
        
        df_results['species'] = species_name
        all_results_long.append(df_results)

    if not all_results_long:
        print("  ❌ Error: No GRAVY/Aromo data could be calculated.")
        return
        
    df_plot_long = pd.concat(all_results_long)
    status_queue.put(("progress", 80))
    
    try:
        print("  Generating Chart 1: GRAVY vs Aromo Scatter Plot...")
        plt.figure(figsize=(14, 10))
        sns.scatterplot(data=df_plot_long, x='gravy', y='aromo', hue='species', 
                        palette=palette, alpha=0.6, s=30, edgecolor='w', linewidth=0.2)
        
        plt.title("1. Physicochemical Space (GRAVY vs Aromaticity)", fontsize=16, fontweight='bold')
        plt.xlabel("GRAVY Score (Hydropathicity)", fontsize=12)
        plt.ylabel("Aromaticity Score (% F, Y, W)", fontsize=12)
        
        plt.axhline(df_plot_long['aromo'].mean(), color='grey', linestyle='--', linewidth=1, alpha=0.5)
        plt.axvline(df_plot_long['gravy'].mean(), color='grey', linestyle='--', linewidth=1, alpha=0.5)
        
        plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', title="Species", borderaxespad=0)
        plt.grid(alpha=0.3)
        plt.tight_layout()
        
        plot_path1 = os.path.join(output_folder, "gravy_aromo_1_scatter.png")
        plt.savefig(plot_path1, dpi=150, bbox_inches="tight")
        plt.close()
        status_queue.put(("image_ready", (plot_path1, "1. GRAVY & Aromo Scatter Space"))) 
    except Exception as e:
        print(f"  ❌ Error generating Scatter Plot: {e}")

    try:
        print("  Generating Chart 2: GRAVY Boxplot...")
        plt.figure(figsize=(14, 6))
        sns.boxplot(x='species', y='gravy', data=df_plot_long, hue='species', palette=palette, legend=False, showfliers=False)
        plt.title("2. GRAVY Score Distribution by Species", fontsize=16)
        plt.ylabel("GRAVY Score")
        plt.xlabel("Species")
        plt.xticks(rotation=45, ha='right')
        plt.grid(axis='y', alpha=0.3)
        plt.tight_layout()
        
        plot_path2 = os.path.join(output_folder, "gravy_aromo_2_boxplot_gravy.png")
        plt.savefig(plot_path2, dpi=150, bbox_inches="tight")
        plt.close()
        status_queue.put(("image_ready", (plot_path2, "2. GRAVY Boxplot by Species")))
    except Exception as e:
        print(f"  ❌ Error generating GRAVY Boxplot: {e}")

    try:
        print("  Generating Chart 3: Aromaticity Boxplot...")
        plt.figure(figsize=(14, 6))
        sns.boxplot(x='species', y='aromo', data=df_plot_long, hue='species', palette=palette, legend=False, showfliers=False)
        plt.title("3. Aromaticity Score Distribution by Species", fontsize=16)
        plt.ylabel("Aromaticity Score (% F, Y, W)")
        plt.xlabel("Species")
        plt.xticks(rotation=45, ha='right')
        plt.grid(axis='y', alpha=0.3)
        plt.tight_layout()
        
        plot_path3 = os.path.join(output_folder, "gravy_aromo_3_boxplot_aromo.png")
        plt.savefig(plot_path3, dpi=150, bbox_inches="tight")
        plt.close()
        status_queue.put(("image_ready", (plot_path3, "3. Aromo Boxplot by Species")))
    except Exception as e:
        print(f"  ❌ Error generating Aromo Boxplot: {e}")

    try:
        print("  Generating Chart 4: Hexbin Density Map...")
        plt.figure(figsize=(12, 10))
        plt.hexbin(df_plot_long['gravy'], df_plot_long['aromo'], gridsize=40, cmap=palette, mincnt=1)
        cb = plt.colorbar(label='Number of Genes (Density)')
        plt.title("4. Density Map (Hexbin) of GRAVY vs Aromaticity (All Species)", fontsize=16)
        plt.xlabel("GRAVY Score")
        plt.ylabel("Aromaticity Score")
        plt.grid(alpha=0.3)
        plt.tight_layout()
        
        plot_path4 = os.path.join(output_folder, "gravy_aromo_4_hexbin_density.png")
        plt.savefig(plot_path4, dpi=150, bbox_inches="tight")
        plt.close()
        status_queue.put(("image_ready", (plot_path4, "4. GRAVY & Aromo Hexbin Density")))
    except Exception as e:
        print(f"  ❌ Error generating Hexbin Plot: {e}")

def dinucleotide_composition_analysis(file_list, output_folder, status_queue, palette='viridis'):
    print(f"\n=== DINUCLEOTIDE COMPOSITION ANALYSIS ===")
    all_results = {}
    total_files = len(file_list)
    dinu_order = [''.join(p) for p in product('ATGC', repeat=2)]
    
    for i, full_path in enumerate(file_list):
        base_name = os.path.basename(full_path).split('.')[0]
        print(f"  Analyzing dinucleotides of {base_name}...")
        status_queue.put(("message", f"Analyzing dinucleotides of {base_name}..."))
        
        dinu_counts = Counter()
        mono_counts = Counter()
        total_pairs = 0

        try:
            for record in SeqIO.parse(full_path, "genbank"):
                seq = str(record.seq).upper()
                for j in range(len(seq) - 1):
                    n1, n2 = seq[j], seq[j + 1]
                    if n1 in 'ATGC' and n2 in 'ATGC':
                        dinu_counts[n1 + n2] += 1
                        mono_counts[n1] += 1
                        total_pairs += 1
                # last nucleotide only contributes to mono frequency
                if seq and seq[-1] in 'ATGC':
                    mono_counts[seq[-1]] += 1

            if total_pairs > 0:
                total_nucs = sum(mono_counts.values())
                mono_freq = {n: mono_counts[n] / total_nucs for n in 'ATGC'}

                # Karlin-Burge (1995) odds ratio: ρ_XY = f(XY) / (f(X) × f(Y))
                # ρ > 1 = over-represented; ρ < 1 = under-represented (e.g., CpG suppression)
                rho = {}
                for dinu in dinu_order:
                    f_xy = dinu_counts.get(dinu, 0) / total_pairs
                    f_x = mono_freq.get(dinu[0], 0)
                    f_y = mono_freq.get(dinu[1], 0)
                    rho[dinu] = f_xy / (f_x * f_y) if (f_x * f_y) > 0 else 0.0

                all_results[base_name] = rho
                cg_rho = rho.get('CG', 0.0)
                print(f"  [{base_name}] ρ_CG = {cg_rho:.3f} {'(CpG suppressed)' if cg_rho < 1 else ''}")
            else:
                print(f"  Warning: No dinucleotide pair found in {base_name}")

        except Exception as e:
            print(f"  ❌ Error analyzing dinucleotides in {base_name}: {e}")
            
        status_queue.put(("progress", int(20 + (i / total_files) * 70)))
        
    if not all_results:
        print("  ❌ Error: No dinucleotide data processed.")
        return

    df_dinu = pd.DataFrame.from_dict(all_results, orient='index', columns=dinu_order)
    df_dinu.index.name = 'Species'
    csv_path = os.path.join(output_folder, 'dinucleotide_rho_karlin_burge.csv')
    df_dinu.to_csv(csv_path, sep=';', decimal='.')
    print(f"  ✅ Dinucleotide ρ_XY table (Karlin-Burge) saved in: {csv_path}")

    try:
        print("  Generating Chart 1: Dinucleotide Odds Ratio Heatmap (Karlin-Burge)...")
        status_queue.put(("progress", 90))
        plt.figure(figsize=(16, max(8, len(df_dinu) * 0.5)))
        # Diverging colormap centred at ρ=1.0 (neutral expectation)
        sns.heatmap(df_dinu, annot=True, fmt=".3f", cmap='RdBu_r', center=1.0,
                    linewidths=0.5, cbar_kws={'label': 'ρ_XY (Karlin-Burge odds ratio)'})
        plt.title("1. Dinucleotide Relative Abundance ρ_XY — Karlin & Burge (1995)\n"
                  "ρ < 1 = under-represented   |   ρ > 1 = over-represented   |   ρ = 1 = neutral",
                  fontsize=14)
        plt.xlabel("Dinucleotide")
        plt.ylabel("Species")
        plt.tight_layout()

        plot_path1 = os.path.join(output_folder, 'dinucleotide_1_heatmap.png')
        plt.savefig(plot_path1, dpi=150, bbox_inches="tight")
        plt.close()
        status_queue.put(("image_ready", (plot_path1, "1. Dinucleotide ρ_XY Heatmap")))
    except Exception as e:
        print(f"  ❌ Error generating Heatmap: {e}")

    try:
        print("  Generating Chart 2: Dinucleotide ρ_XY Variance Boxplot...")
        plt.figure(figsize=(16, 6))
        sns.boxplot(data=df_dinu, color=sns.color_palette(palette, 1)[0], showfliers=False)
        plt.axhline(y=1.0, color='red', linestyle='--', linewidth=1.5,
                    label='ρ = 1.0 (neutral expectation)')
        plt.title("2. Variation of Dinucleotide ρ_XY Across All Species", fontsize=16)
        plt.xlabel("Dinucleotide Pair")
        plt.ylabel("ρ_XY (Karlin-Burge odds ratio)")
        plt.legend()
        plt.grid(axis='y', alpha=0.3)
        plt.tight_layout()

        plot_path2 = os.path.join(output_folder, 'dinucleotide_2_variance_boxplot.png')
        plt.savefig(plot_path2, dpi=150, bbox_inches="tight")
        plt.close()
        status_queue.put(("image_ready", (plot_path2, "2. Dinucleotide ρ_XY Boxplot")))
    except Exception as e:
        print(f"  ❌ Error generating Boxplot: {e}")

    try:
        plot_col = sns.color_palette(palette)[0]
    except:
        plot_col = 'teal'

    try:
        print("  Generating Chart 3: Average Dinucleotide ρ_XY Barplot...")
        plt.figure(figsize=(16, 6))
        mean_rho = df_dinu.mean()
        std_rho = df_dinu.std()

        mean_rho.plot(kind='bar', yerr=std_rho, capsize=4, color=plot_col,
                      edgecolor='black', alpha=0.8)
        plt.axhline(y=1.0, color='red', linestyle='--', linewidth=1.5,
                    label='ρ = 1.0 (neutral)')
        plt.title("3. Mean Dinucleotide ρ_XY Across All Species (± StdDev)", fontsize=16)
        plt.xlabel("Dinucleotide Pair")
        plt.ylabel("Mean ρ_XY (Karlin-Burge odds ratio)")
        plt.legend()
        plt.grid(axis='y', alpha=0.3)
        plt.xticks(rotation=45)
        plt.tight_layout()

        plot_path3 = os.path.join(output_folder, 'dinucleotide_3_mean_barplot.png')
        plt.savefig(plot_path3, dpi=150, bbox_inches="tight")
        plt.close()
        status_queue.put(("image_ready", (plot_path3, "3. Mean Dinucleotide ρ_XY Barplot")))
    except Exception as e:
        print(f"  ❌ Error generating Barplot: {e}")

def calculate_pr2_per_gene(sequences, genetic_code_id):
    genetic_code = GENETIC_CODE_TABLES.get(genetic_code_id, GENETIC_CODE_TABLES[1])
    stop_codons = {codon for codon, aa in genetic_code.items() if aa == '*'}
    results = []

    for seq in sequences:
        counts = {'A3': 0, 'T3': 0, 'G3': 0, 'C3': 0}
        
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i+3]
            if codon in stop_codons:
                continue 
                
            nuc3 = codon[2]
            if nuc3 in ['A', 'T', 'G', 'C']:
                counts[f"{nuc3}3"] += 1
        
        a3_plus_t3 = counts['A3'] + counts['T3']
        g3_plus_c3 = counts['G3'] + counts['C3']
        
        a3_frac = (counts['A3'] / a3_plus_t3) if a3_plus_t3 > 0 else np.nan
        g3_frac = (counts['G3'] / g3_plus_c3) if g3_plus_c3 > 0 else np.nan
        
        at_skew = (counts['A3'] - counts['T3']) / a3_plus_t3 if a3_plus_t3 > 0 else np.nan
        gc_skew = (counts['G3'] - counts['C3']) / g3_plus_c3 if g3_plus_c3 > 0 else np.nan
        
        if not (np.isnan(a3_frac) or np.isnan(g3_frac)):
            results.append({
                'A3_frac': a3_frac, 
                'G3_frac': g3_frac,
                'AT_skew': at_skew,
                'GC_skew': gc_skew
            })

    return pd.DataFrame(results)

def pr2_plot_analysis(file_list, output_folder, genetic_code_id, status_queue, gene_list=None, palette='viridis'):
    print(f"\n=== PR2 PARITY PLOT ANALYSIS ===")
    all_seqs_by_species = extract_cds_sequences(file_list, status_queue, gene_list)
    
    if not all_seqs_by_species:
        print("  ❌ Error: No valid CDS sequence was extracted.")
        return

    status_queue.put(("progress", 40))
    all_results_long = []

    for species_name, sequences in all_seqs_by_species.items():
        print(f"  Calculating PR2 for {species_name}...")
        status_queue.put(("message", f"Calculating PR2 for {species_name}..."))
        
        df_results = calculate_pr2_per_gene(sequences, genetic_code_id)
        if df_results.empty:
            continue
            
        df_results['species'] = species_name
        all_results_long.append(df_results)
    
    if not all_results_long:
        print("  ❌ Error: No PR2 data could be calculated.")
        return
        
    df_plot_long = pd.concat(all_results_long).dropna()
    df_plot_long.to_csv(os.path.join(output_folder, 'pr2_plot_data_per_gene.csv'), sep=';', index=False)
    
    status_queue.put(("progress", 80))

    try:
        print("  Generating Chart 1: PR2 Scatter Plot...")
        plt.figure(figsize=(10, 10))
        ax = sns.scatterplot(x='G3_frac', y='A3_frac', data=df_plot_long, hue='species', palette=palette, alpha=0.5, s=20)
        
        ax.axhline(0.5, color='black', linestyle='--', linewidth=1)
        ax.axvline(0.5, color='black', linestyle='--', linewidth=1)
        
        ax.set_xlabel('G3 / (G3 + C3)', fontsize=12)
        ax.set_ylabel('A3 / (A3 + T3)', fontsize=12)
        ax.set_title('1. PR2 Parity Plot (3rd Position Bias)', fontsize=16)
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1))
        
        plt.tight_layout()
        plot_path1 = os.path.join(output_folder, 'pr2_1_scatter.png')
        plt.savefig(plot_path1, dpi=150, bbox_inches='tight')
        plt.close()
        status_queue.put(("image_ready", (plot_path1, "1. PR2 Parity Scatter")))
    except Exception as e:
        print(f"  ❌ Error generating Scatter Plot: {e}")

    try:
        print("  Generating Chart 2: PR2 Density KDE...")
        plt.figure(figsize=(10, 10))
        sns.kdeplot(x=df_plot_long['G3_frac'], y=df_plot_long['A3_frac'], cmap=palette, fill=True, thresh=0.05)
        plt.axhline(0.5, color='black', linestyle='--', linewidth=1)
        plt.axvline(0.5, color='black', linestyle='--', linewidth=1)
        plt.xlabel('G3 / (G3 + C3)', fontsize=12)
        plt.ylabel('A3 / (A3 + T3)', fontsize=12)
        plt.title('2. PR2 Parity Plot Density (KDE)', fontsize=16)
        plt.xlim(0, 1); plt.ylim(0, 1)
        plt.tight_layout()
        
        plot_path2 = os.path.join(output_folder, 'pr2_2_kde_density.png')
        plt.savefig(plot_path2, dpi=150, bbox_inches='tight')
        plt.close()
        status_queue.put(("image_ready", (plot_path2, "2. PR2 KDE Density")))
    except Exception as e:
        print(f"  ❌ Error generating KDE Plot: {e}")

    try:
        print("  Generating Chart 3: AT and GC Skews Boxplot...")
        df_skew = df_plot_long.melt(id_vars='species', value_vars=['AT_skew', 'GC_skew'], 
                                    var_name='Skew_Type', value_name='Skew_Value')
        plt.figure(figsize=(14, 6))
        sns.boxplot(x='species', y='Skew_Value', hue='Skew_Type', data=df_skew, palette=palette, showfliers=False)
        plt.axhline(0, color='red', linestyle='--', linewidth=1.5, label='Zero Skew')
        plt.title("3. Distribution of AT and GC Skews at 3rd Position", fontsize=16)
        plt.xlabel("Species")
        plt.ylabel("Skew Value")
        plt.xticks(rotation=45, ha='right')
        plt.legend(loc='upper left', bbox_to_anchor=(1.02, 1))
        plt.grid(axis='y', alpha=0.3)
        plt.tight_layout()
        
        plot_path3 = os.path.join(output_folder, 'pr2_3_skews_boxplot.png')
        plt.savefig(plot_path3, dpi=150, bbox_inches='tight')
        plt.close()
        status_queue.put(("image_ready", (plot_path3, "3. AT/GC Skews Boxplot")))
    except Exception as e:
        print(f"  ❌ Error generating Skews Boxplot: {e}")

    try:
        col1 = sns.color_palette(palette)[0]
        col2 = sns.color_palette(palette)[-1]
    except:
        col1, col2 = 'blue', 'orange'

    try:
        print("  Generating Chart 4: G3 and A3 Histograms...")
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
        
        sns.histplot(df_plot_long['G3_frac'], bins=50, kde=True, ax=ax1, color=col1, edgecolor='black')
        ax1.set_title("G3 / (G3+C3) Distribution", fontsize=14)
        ax1.axvline(0.5, color='red', linestyle='--', linewidth=2)
        ax1.set_xlabel("G3 Fraction")
        
        sns.histplot(df_plot_long['A3_frac'], bins=50, kde=True, ax=ax2, color=col2, edgecolor='black')
        ax2.set_title("A3 / (A3+T3) Distribution", fontsize=14)
        ax2.axvline(0.5, color='red', linestyle='--', linewidth=2)
        ax2.set_xlabel("A3 Fraction")
        
        plt.suptitle("4. Independent Distributions of 3rd Position Purine Bias", fontsize=16)
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        
        plot_path4 = os.path.join(output_folder, 'pr2_4_histograms.png')
        plt.savefig(plot_path4, dpi=150, bbox_inches='tight')
        plt.close()
        status_queue.put(("image_ready", (plot_path4, "4. G3 and A3 Histograms")))
    except Exception as e:
        print(f"  ❌ Error generating Histograms: {e}")

def count_tRNAs(file_list, status_queue):
    print("  Analyzing tRNA genes (anticodons)...")
    anticodon_counts = Counter()
    total_tRNAs = 0
    
    qualifier_regex = re.compile(r'seq:\s*([ATGCU]{3})', re.IGNORECASE)
    note_regex = re.compile(r'anticodon:\s*([ATGCU]{3})', re.IGNORECASE)
    
    for i, full_path in enumerate(file_list):
        base_name = os.path.basename(full_path).split('.')[0]
        status_queue.put(("message", f"Counting tRNAs in {base_name}..."))
        
        try:
            for record in SeqIO.parse(full_path, "genbank"):
                for feature in record.features:
                    if feature.type == "tRNA":
                        anticodon_str = ""
                        
                        if "anticodon" in feature.qualifiers:
                            qual_value = feature.qualifiers.get("anticodon")[0]
                            match = qualifier_regex.search(qual_value)
                            if match:
                                anticodon_str = match.group(1)
                        
                        if not anticodon_str:
                            note = str(feature.qualifiers.get("note", ""))
                            match = note_regex.search(note)
                            if match:
                                anticodon_str = match.group(1)
                        
                        if anticodon_str:
                            anticodon_norm = anticodon_str.upper().replace('T', 'U')
                            if len(anticodon_norm) == 3 and all(b in 'ATGCU' for b in anticodon_norm):
                                anticodon_counts[anticodon_norm] += 1
                                total_tRNAs += 1
                                    
        except Exception as e:
            print(f"  ❌ Error counting tRNAs in {base_name}: {e}")
            
    if total_tRNAs == 0:
         print("  ⚠️ Warning: No tRNA gene with recognizable anticodon was found.")
    else:
         print(f"  ✅ Found {total_tRNAs} tRNAs with defined anticodons.")
         
    return anticodon_counts

def calculate_wobble_W_weights(anticodon_counts, genetic_code_id, wobble_matrix):
    print("  Calculating W weights weighted by wobble (dos Reis)...")
    genetic_code = GENETIC_CODE_TABLES.get(genetic_code_id, GENETIC_CODE_TABLES[1])
    W_values = {}
    
    total_tRNAs = sum(anticodon_counts.values())
    if total_tRNAs == 0:
        print("  Warning: No tRNA count. Using pseudocounts.")
        for aa, codons in AA_CODON_MAPS[genetic_code_id].items():
            if aa == '*': continue
            for codon in codons:
                rc_anticodon = str(Seq(codon).reverse_complement()).replace('T', 'U')
                anticodon_counts[rc_anticodon] = 1
        total_tRNAs = sum(anticodon_counts.values())

    relative_abundance = {ac: count / total_tRNAs for ac, count in anticodon_counts.items()}
    
    for codon, aa in genetic_code.items():
        if aa == '*': continue
        
        codon_base_3 = codon[2].replace('T', 'U')  # wobble matrix keys use RNA 'U', not DNA 'T'
        codon_bases_1_2 = codon[:2]
        sum_wi = 0.0

        for anticodon_5_3, n_j in relative_abundance.items():
            # anticodon 5'→3' = [A34, A35, A36]; codon N1↔A36, N2↔A35 (Watson-Crick, antiparallel)
            if codon_bases_1_2 != str(Seq(anticodon_5_3[2:0:-1]).complement()):
                continue 
            
            anticodon_base_1_raw = anticodon_5_3[0].upper().replace('T', 'U')
            
            anticodon_base_1_mod = ANTICODON_MODIFICATION_MAP.get(anticodon_base_1_raw, anticodon_base_1_raw)
            s_ij = wobble_matrix.get((anticodon_base_1_mod, codon_base_3), 1.0)
            sum_wi += (1.0 - s_ij) * n_j
            
        W_values[codon] = sum_wi

    max_w = max(W_values.values())
    if max_w > 0:
        for codon in W_values:
            W_values[codon] /= max_w
    else:
        print("  Warning: Max W is 0. All tAI weights will be 0.")

    min_val = 1e-9
    for codon in W_values:
        if W_values[codon] < min_val:
            W_values[codon] = min_val
            
    print("  ✅ W weight calculation completed.")
    return W_values

def calculate_tai_per_gene(sequences, W_values):
    results = []
    for seq in sequences:
        log_w_sum = 0
        num_codons_validos = 0
        
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i+3]
            if codon in W_values:
                log_w_sum += np.log(W_values[codon])
                num_codons_validos += 1
        
        if num_codons_validos > 0:
            tai = np.exp(log_w_sum / num_codons_validos)
            results.append({'tAI': tai, 'gene_length_codons': num_codons_validos})
        
    return pd.DataFrame(results)

def tai_analysis(file_list, output_folder, genetic_code_id, status_queue, gene_list=None, super_kingdom="Bacteria", palette='viridis'):
    print(f"\n=== tRNA ADAPTATION INDEX ANALYSIS (Wobble-Weighted tAI) ===")
    wobble_rules_map = {"Bacteria": WOBBLE_MATRIX_BACTERIA, "Eukaryote": WOBBLE_MATRIX_EUKARYA}
    wobble_rules = wobble_rules_map.get(super_kingdom, WOBBLE_MATRIX_BACTERIA)
    print(f"  Using Wobble rules for: {super_kingdom}")

    status_queue.put(("progress", 20))
    anticodon_counts = count_tRNAs(file_list, status_queue)
    df_counts = pd.DataFrame(anticodon_counts.items(), columns=['Anticodon', 'Count']).sort_values(by='Count', ascending=False)
    df_counts.to_csv(os.path.join(output_folder, 'tRNA_anticodon_counts.csv'), sep=';', index=False)
    print(f"  ✅ Anticodon count saved.")

    status_queue.put(("progress", 30))
    W_values = calculate_wobble_W_weights(anticodon_counts, genetic_code_id, wobble_rules)
    df_W = pd.DataFrame(W_values.items(), columns=['Codon', 'Weight_W']).sort_values(by='Weight_W', ascending=False)
    df_W.to_csv(os.path.join(output_folder, f'tAI_codon_weights_wobble_{super_kingdom}.csv'), sep=';', index=False)
    
    all_seqs_by_species = extract_cds_sequences(file_list, status_queue, gene_list)
    if not all_seqs_by_species:
        print("  ❌ Error: No valid CDS sequence was extracted.")
        return

    status_queue.put(("progress", 60))
    all_results_long = []
    for species_name, sequences in all_seqs_by_species.items():
        print(f"  Calculating tAI (Wobble) for {species_name}...")
        status_queue.put(("message", f"Calculating tAI for {species_name}..."))
        
        df_results = calculate_tai_per_gene(sequences, W_values)
        if df_results.empty: continue
            
        df_results['species'] = species_name
        all_results_long.append(df_results)
        
    if not all_results_long:
        print("  ❌ Error: No tAI data could be calculated.")
        return
        
    df_plot_long = pd.concat(all_results_long).dropna()
    df_plot_long.to_csv(os.path.join(output_folder, 'tAI_wobble_data_per_gene.csv'), sep=';', index=False)
    
    status_queue.put(("progress", 85))
    tAI_floor = 1.1e-9 
    df_filtered_plot = df_plot_long[df_plot_long['tAI'] > tAI_floor]
    
    df_to_plot = df_filtered_plot if not df_filtered_plot.empty else df_plot_long

    try:
        print("  Generating Chart 1: tAI Boxplot...")
        plt.figure(figsize=(14, 8))
        plot_title_str = f'1. Wobble-Weighted tAI Distribution ({super_kingdom})'
        if df_filtered_plot.empty: plot_title_str += "\n(Warning: All values at minimum floor)"
        
        sns.boxplot(x='species', y='tAI', data=df_to_plot, hue='species', palette=palette, legend=False, showfliers=False)
        plt.title(plot_title_str, fontsize=16)
        plt.ylabel('tAI (Geometric Mean of W Weights)')
        plt.xlabel('Species')
        plt.yscale('log')
        plt.xticks(rotation=45, ha='right')
        plt.grid(axis='y', alpha=0.3)
        plt.tight_layout()
        
        plot_path1 = os.path.join(output_folder, 'tai_1_comparative_boxplot.png')
        plt.savefig(plot_path1, dpi=150, bbox_inches='tight')
        plt.close()
        status_queue.put(("image_ready", (plot_path1, "1. tAI Boxplot")))
    except Exception as e:
        print(f"  ❌ Error: {e}")

    try:
        print("  Generating Chart 2: tAI vs Length Scatter Plot...")
        plt.figure(figsize=(12, 8))
        sns.scatterplot(data=df_to_plot, x='gene_length_codons', y='tAI', hue='species', palette=palette, alpha=0.5, s=20)
        plt.title(f'2. tAI vs Gene Length ({super_kingdom})', fontsize=16)
        plt.xlabel('Gene Length (Number of Codons)')
        plt.ylabel('tAI')
        plt.xscale('log'); plt.yscale('log')
        plt.grid(alpha=0.3)
        plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left')
        plt.tight_layout()
        
        plot_path2 = os.path.join(output_folder, 'tai_2_vs_length_scatter.png')
        plt.savefig(plot_path2, dpi=150, bbox_inches='tight')
        plt.close()
        status_queue.put(("image_ready", (plot_path2, "2. tAI vs Gene Length")))
    except Exception as e:
        print(f"  ❌ Error: {e}")

    try:
        print("  Generating Chart 3: tAI KDE Density...")
        plt.figure(figsize=(12, 8))
        sns.kdeplot(data=df_to_plot, x='tAI', hue='species', palette=palette, fill=True, common_norm=False, log_scale=True, alpha=0.4)
        plt.title(f'3. Density Plot of tAI Values ({super_kingdom})', fontsize=16)
        plt.xlabel('tAI (Log Scale)')
        plt.ylabel('Density')
        plt.grid(alpha=0.3)
        plt.tight_layout()
        
        plot_path3 = os.path.join(output_folder, 'tai_3_kde_density.png')
        plt.savefig(plot_path3, dpi=150, bbox_inches='tight')
        plt.close()
        status_queue.put(("image_ready", (plot_path3, "3. tAI Density (KDE)")))
    except Exception as e:
        print(f"  ❌ Error: {e}")

    try:
        print("  Generating Chart 4: Calculated W-Weights Barplot...")
        plt.figure(figsize=(18, 6))
        sns.barplot(x='Codon', y='Weight_W', data=df_W, hue='Codon', palette=palette, legend=False)
        plt.title(f"4. Calculated Codon Adaptation Weights (W) based on Wobble Rules", fontsize=16)
        plt.xticks(rotation=90, fontsize=10)
        plt.ylabel("W-Weight (Relative Adaptability)")
        plt.xlabel("Codon")
        plt.grid(axis='y', alpha=0.3)
        plt.tight_layout()
        
        plot_path4 = os.path.join(output_folder, 'tai_4_w_weights_barplot.png')
        plt.savefig(plot_path4, dpi=150, bbox_inches='tight')
        plt.close()
        status_queue.put(("image_ready", (plot_path4, "4. W-Weights Barplot")))
    except Exception as e:
        print(f"  ❌ Error: {e}")

def upstream_motifs_analysis(file_list, output_folder, status_queue, gene_list=None, upstream_dist=200, kmer_size=6, palette='viridis'):
    print(f"\n=== UPSTREAM MOTIFS ANALYSIS ===")
    print(f"  Parameters: Distance={upstream_dist}bp, K-mer={kmer_size}bp")
    
    total_files = len(file_list)
    
    for i, full_path in enumerate(file_list):
        base_name = os.path.basename(full_path).split('.')[0]
        print(f"\n  Analyzing motifs in {base_name}...")
        status_queue.put(("message", f"Analyzing motifs in {base_name}..."))
        status_queue.put(("progress", int(10 + (i / total_files) * 80)))
        
        all_upstream_seqs = []
        processed_genes = 0
        
        try:
            for record in SeqIO.parse(full_path, "genbank"):
                record_seq_str = str(record.seq).upper()
                record_len = len(record_seq_str)
                
                for feature in record.features:
                    if feature.type != "CDS": continue
                    
                    if not _apply_gene_filter(feature, gene_list):
                        continue
                        
                    seq_upstream = ""
                    try:
                        if feature.location.strand == 1:
                            cds_start = feature.location.start
                            upstream_start = max(0, cds_start - upstream_dist)
                            upstream_end = cds_start
                            if upstream_start < upstream_end:
                                seq_upstream = record_seq_str[upstream_start:upstream_end]
                                
                        elif feature.location.strand == -1:
                            cds_end = feature.location.end
                            upstream_start = cds_end
                            upstream_end = min(record_len, cds_end + upstream_dist)
                            if upstream_start < upstream_end:
                                seq_fwd = record_seq_str[upstream_start:upstream_end]
                                seq_upstream = str(Seq(seq_fwd).reverse_complement())
                        
                        if seq_upstream:
                            all_upstream_seqs.append(seq_upstream)
                            processed_genes += 1

                    except Exception as e_feature:
                        print(f"    Warning: Skipping gene (location error): {e_feature}")

            print(f"  Extracted {len(all_upstream_seqs)} upstream regions from {processed_genes} genes.")
            if not all_upstream_seqs:
                print("  ❌ No upstream region found.")
                continue

            kmer_counts = Counter()
            for seq in all_upstream_seqs:
                for j in range(len(seq) - kmer_size + 1):
                    kmer = seq[j:j+kmer_size]
                    if 'N' not in kmer and all(b in 'ATGC' for b in kmer):
                        kmer_counts[kmer] += 1
            
            if not kmer_counts:
                print("  ❌ No valid k-mer found.")
                continue

            df_kmers = pd.DataFrame(kmer_counts.most_common(), columns=['K-mer', 'Count'])
            df_top_kmers = df_kmers.head(25)
            df_kmers.to_csv(os.path.join(output_folder, f'motif_counts_{kmer_size}mer_{base_name}.csv'), sep=';', index=False)
            
            try:
                print("  Generating Chart 1: Top 25 K-mers Barplot...")
                plt.figure(figsize=(10, 12))
                sns.barplot(x='Count', y='K-mer', data=df_top_kmers, hue='K-mer', palette=palette, legend=False)
                plt.title(f'1. Top 25 Most Frequent {kmer_size}-mers (Upstream {upstream_dist}bp)\n{base_name}', fontsize=14)
                plt.xlabel('Absolute Count', fontsize=12)
                plt.ylabel(f'{kmer_size}-mer', fontsize=12)
                plt.grid(axis='x', alpha=0.3)
                plt.tight_layout()
                
                plot_path1 = os.path.join(output_folder, f'motif_1_barplot_{kmer_size}mer_{base_name}.png')
                plt.savefig(plot_path1, dpi=150, bbox_inches='tight')
                plt.close()
                status_queue.put(("image_ready", (plot_path1, f"1. Top 25 Motifs - {base_name}")))
            except Exception as e:
                print(f"  ❌ Error: {e}")

            try:
                print("  Generating Chart 2: Rank-Frequency Curve...")
                plt.figure(figsize=(10, 6))
                
                try:
                    plot_col = sns.color_palette(palette)[-1]
                except:
                    plot_col = 'red'

                plt.plot(range(1, len(df_kmers) + 1), df_kmers['Count'], color=plot_col, linewidth=2)
                plt.xscale('log')
                plt.yscale('log')
                plt.title(f"2. Motif Frequency Rank Curve (Zipf's Law) - {base_name}", fontsize=14)
                plt.xlabel('Rank (Log scale)')
                plt.ylabel('Frequency Count (Log scale)')
                plt.grid(alpha=0.3)
                plt.tight_layout()
                
                plot_path2 = os.path.join(output_folder, f'motif_2_rank_curve_{base_name}.png')
                plt.savefig(plot_path2, dpi=150, bbox_inches='tight')
                plt.close()
                status_queue.put(("image_ready", (plot_path2, f"2. Motif Rank Curve - {base_name}")))
            except Exception as e:
                print(f"  ❌ Error: {e}")

            try:
                print("  Generating Chart 3: GC Content vs Frequency Boxplot...")
                df_kmers['GC_Content'] = df_kmers['K-mer'].apply(lambda x: (x.count('G') + x.count('C')) / len(x) * 100)
                df_kmers['GC_Content'] = df_kmers['GC_Content'].round(1)
                
                plt.figure(figsize=(12, 8))
                sns.boxplot(x='GC_Content', y='Count', data=df_kmers, hue='GC_Content', palette=palette, legend=False, showfliers=False)
                plt.title(f"3. K-mer GC Content vs Motif Frequency - {base_name}", fontsize=14)
                plt.xlabel('GC Content (%)')
                plt.ylabel('Absolute Frequency Count')
                plt.grid(axis='y', alpha=0.3)
                plt.tight_layout()
                
                plot_path3 = os.path.join(output_folder, f'motif_3_gc_vs_freq_boxplot_{base_name}.png')
                plt.savefig(plot_path3, dpi=150, bbox_inches='tight')
                plt.close()
                status_queue.put(("image_ready", (plot_path3, f"3. GC vs Motif Freq Boxplot - {base_name}")))
            except Exception as e:
                print(f"  ❌ Error: {e}")

        except Exception as e_file:
            print(f"  ❌ General error processing motifs in {base_name}: {e_file}")
