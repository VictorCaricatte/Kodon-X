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
    genetic_code = GENETIC_CODE_TABLES.get(genetic_code_id, GENETIC_CODE_TABLES[1])
    stop_codons = {codon for codon, aa in genetic_code.items() if aa == '*'}
    
    non_stop_codons = sorted([c for c in ALL_CODONS_SORTED if c not in stop_codons])
    codon_index = {codon: i for i, codon in enumerate(non_stop_codons)}
    
    num_codons = len(non_stop_codons)
    pair_counts_obs = np.zeros((num_codons, num_codons))
    codon_counts_c1 = np.zeros(num_codons) 
    codon_counts_c2 = np.zeros(num_codons) 
    total_pairs = 0
    
    for seq in sequences:
        for i in range(0, len(seq) - 5, 3): 
            c1 = seq[i:i+3]
            c2 = seq[i+3:i+6]
            
            if c1 in codon_index and c2 in codon_index:
                idx1 = codon_index[c1]
                idx2 = codon_index[c2]
                
                pair_counts_obs[idx1, idx2] += 1
                codon_counts_c1[idx1] += 1
                codon_counts_c2[idx2] += 1
                total_pairs += 1
    
    if total_pairs == 0:
        print("  Warning: No valid codon pair found.")
        return pd.DataFrame()
        
    pair_counts_exp = np.outer(codon_counts_c1, codon_counts_c2) / total_pairs
    epsilon = 1e-9
    cps_matrix = np.log((pair_counts_obs + epsilon) / (pair_counts_exp + epsilon))
    
    df_cps = pd.DataFrame(cps_matrix, index=non_stop_codons, columns=non_stop_codons)
    return df_cps

def codon_pair_bias_analysis(file_list, output_folder, genetic_code_id, status_queue, gene_list=None):
    print(f"\n=== CODON PAIR BIAS ANALYSIS (CPB) ===")
    all_seqs_by_species = extract_cds_sequences(file_list, status_queue, gene_list)
    
    if not all_seqs_by_species:
        print("  ❌ Error: No valid CDS sequence was extracted.")
        return

    status_queue.put(("progress", 40))
    
    for i, (species_name, sequences) in enumerate(all_seqs_by_species.items()):
        print(f"  Calculating CPB for {species_name}...")
        status_queue.put(("message", f"Calculating CPB for {species_name}..."))
        
        df_cps = calculate_codon_pair_bias(sequences, genetic_code_id)
        
        if df_cps.empty:
            continue
            
        csv_path = os.path.join(output_folder, f"cpb_matrix_{species_name}.csv")
        df_cps.to_csv(csv_path, sep=';', decimal='.')
        print(f"  ✅ CPB Matrix saved in: {csv_path}")
        
        print("  Generating CPB Heatmap...")
        plt.figure(figsize=(24, 20))
        sns.heatmap(df_cps, cmap="coolwarm", center=0, annot=False, 
                    cbar_kws={'label': 'Codon Pair Score (log(Obs/Exp))'})
        plt.title(f"Codon Pair Bias (CPB) - {species_name}", fontsize=18)
        plt.xlabel("Second Codon", fontsize=12)
        plt.ylabel("First Codon", fontsize=12)
        plt.tight_layout()
        
        plot_path = os.path.join(output_folder, f"cpb_heatmap_{species_name}.png")
        plt.savefig(plot_path, dpi=100, bbox_inches="tight") 
        plt.close()
        
        print(f"  ✅ CPB Heatmap saved in: {plot_path}")
        status_queue.put(("image_ready", (plot_path, f"CPB - {species_name}")))
        status_queue.put(("progress", int(40 + (i / len(all_seqs_by_species)) * 50)))

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

def gravy_aromo_analysis(file_list, output_folder, genetic_code_id, status_queue, gene_list=None):
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
    print("  Generating comparative Box Plots...")
    status_queue.put(("progress", 80))
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    
    sns.violinplot(x='species', y='gravy', data=df_plot_long, ax=ax1, palette="coolwarm", inner="quartile")
    ax1.set_title("GRAVY Distribution (Hydropathicity)", fontweight='bold')
    ax1.set_ylabel("GRAVY Score")
    ax1.set_xlabel("Species")
    ax1.tick_params(axis='x', rotation=45)
    
    sns.violinplot(x='species', y='aromo', data=df_plot_long, ax=ax2, palette="viridis", inner="quartile")
    ax2.set_title("Aromaticity Distribution", fontweight='bold')
    ax2.set_ylabel("Aromo Score (% F, Y, W)")
    ax2.set_xlabel("Species")
    ax2.tick_params(axis='x', rotation=45)
    
    plt.suptitle("Comparative Physicochemical Analysis", fontsize=16, fontweight='bold')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    
    plot_path = os.path.join(output_folder, "gravy_aromo_comparative.png")
    plt.savefig(plot_path, dpi=150, bbox_inches="tight")
    plt.close()
    
    print(f"  ✅ Comparative chart saved in: {plot_path}")
    status_queue.put(("image_ready", (plot_path, "GRAVY & Aromo Analysis"))) 

def dinucleotide_composition_analysis(file_list, output_folder, status_queue):
    print(f"\n=== DINUCLEOTIDE COMPOSITION ANALYSIS ===")
    all_results = {}
    total_files = len(file_list)
    dinu_order = [''.join(p) for p in product('ATGC', repeat=2)]
    
    for i, full_path in enumerate(file_list):
        base_name = os.path.basename(full_path).split('.')[0]
        print(f"  Analyzing dinucleotides of {base_name}...")
        status_queue.put(("message", f"Analyzing dinucleotides of {base_name}..."))
        
        counts = Counter()
        total_pairs = 0
        
        try:
            for record in SeqIO.parse(full_path, "genbank"):
                seq = record.seq.upper()
                for j in range(len(seq) - 1):
                    dinu = seq[j:j+2]
                    if dinu in dinu_order: 
                        counts[dinu] += 1
                        total_pairs += 1
                        
            if total_pairs > 0:
                freqs = {dinu: (counts.get(dinu, 0) / total_pairs) * 100 for dinu in dinu_order}
                all_results[base_name] = freqs
            else:
                print(f"  Warning: No dinucleotide pair found in {base_name}")
                
        except Exception as e:
            print(f"  ❌ Error analyzing dinucleotides in {base_name}: {e}")
            
        status_queue.put(("progress", int(20 + (i / total_files) * 70)))
        
    if not all_results:
        print("  ❌ Error: No dinucleotide data processed.")
        return

    df_dinu = pd.DataFrame.from_dict(all_results, orient='index', columns=dinu_order)
    csv_path = os.path.join(output_folder, 'dinucleotide_composition.csv')
    df_dinu.to_csv(csv_path, sep=';', decimal='.')
    print(f"  ✅ Dinucleotide composition table saved in: {csv_path}")

    print("  Generating Dinucleotide Heatmap...")
    status_queue.put(("progress", 90))
    plt.figure(figsize=(14, max(8, len(df_dinu) * 0.5)))
    sns.heatmap(df_dinu, annot=True, fmt=".2f", cmap="viridis", linewidths=0.5,
                cbar_kws={'label': 'Frequency (%)'})
    plt.title("Dinucleotide Composition (%)", fontsize=16)
    plt.xlabel("Dinucleotide")
    plt.ylabel("Species")
    plt.tight_layout()
    
    plot_path = os.path.join(output_folder, 'dinucleotide_composition_heatmap.png')
    plt.savefig(plot_path, dpi=150, bbox_inches="tight")
    plt.close()
    
    print(f"  ✅ Heatmap saved in: {plot_path}")
    status_queue.put(("image_ready", (plot_path, "Dinucleotide Composition")))

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
        
        if not (np.isnan(a3_frac) or np.isnan(g3_frac)):
            results.append({'A3_frac': a3_frac, 'G3_frac': g3_frac})

    return pd.DataFrame(results)

def pr2_plot_analysis(file_list, output_folder, genetic_code_id, status_queue, gene_list=None):
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
    
    print("  Generating PR2 Plot...")
    status_queue.put(("progress", 80))
    
    plt.figure(figsize=(12, 12))
    ax = sns.scatterplot(x='G3_frac', y='A3_frac', data=df_plot_long, hue='species', alpha=0.5, s=20)
    
    ax.axhline(0.5, color='black', linestyle='--', linewidth=1)
    ax.axvline(0.5, color='black', linestyle='--', linewidth=1)
    
    ax.set_xlabel('G3 / (G3 + C3)', fontsize=12)
    ax.set_ylabel('A3 / (A3 + T3)', fontsize=12)
    ax.set_title(f'PR2 Parity Plot (3rd Position Bias)', fontsize=16)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
    
    plt.tight_layout()
    output_path = os.path.join(output_folder, 'pr2_plot.png')
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"\n✅ PR2 Plot saved in: {output_path}")
    status_queue.put(("image_ready", (output_path, "PR2 Parity Plot")))

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
        
        codon_base_3 = codon[2]
        codon_bases_1_2 = codon[:2]
        sum_wi = 0.0
        
        for anticodon_5_3, n_j in relative_abundance.items():
            if codon_bases_1_2 != str(Seq(anticodon_5_3[1::-1]).complement()):
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

def tai_analysis(file_list, output_folder, genetic_code_id, status_queue, gene_list=None, super_kingdom="Bacteria"):
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
    
    print("  Generating tAI (Wobble) Box Plot...")
    status_queue.put(("progress", 90))
    plt.figure(figsize=(max(10, len(all_seqs_by_species)*1.5), 8))
    
    tAI_floor = 1.1e-9 
    df_filtered_plot = df_plot_long[df_plot_long['tAI'] > tAI_floor]
    
    plot_title_str = f'Wobble-Weighted tAI Distribution ({super_kingdom})'
    
    if df_filtered_plot.empty:
        df_to_plot = df_plot_long
        plot_title_str += "\n(Warning: All values at minimum floor)"
    else:
        df_to_plot = df_filtered_plot
        num_filtrados = len(df_plot_long) - len(df_filtered_plot)
        if num_filtrados > 0:
            print(f"  Info: Filtered {num_filtrados} genes from plot (tAI < {tAI_floor:.1e}).")
            plot_title_str += f"\n(Values < {tAI_floor:.1e} filtered for clarity)"
    
    sns.boxplot(x='species', y='tAI', data=df_to_plot, palette="Set3")
    sns.stripplot(x='species', y='tAI', data=df_to_plot, color=".25", size=2, alpha=0.2)
    
    plt.title(plot_title_str, fontsize=16)
    plt.ylabel('tAI (Geometric Mean of W Weights)')
    plt.xlabel('Species')
    plt.yscale('log')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    
    output_path = os.path.join(output_folder, 'tAI_wobble_comparative_boxplot.png')
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"\n✅ tAI (Wobble) Chart saved in: {output_path}")
    status_queue.put(("image_ready", (output_path, "Weighted tAI Analysis (Wobble)")))

def upstream_motifs_analysis(file_list, output_folder, status_queue, gene_list=None, upstream_dist=200, kmer_size=6):
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

            df_top_kmers = pd.DataFrame(kmer_counts.most_common(25), columns=['K-mer', 'Count'])
            df_top_kmers.to_csv(os.path.join(output_folder, f'motif_counts_{kmer_size}mer_{base_name}.csv'), sep=';', index=False)
            
            print("  Generating K-mers chart...")
            plt.figure(figsize=(10, 12))
            sns.barplot(x='Count', y='K-mer', data=df_top_kmers, palette='viridis')
            plt.title(f'Top 25 most frequent {kmer_size}-mers (Upstream {upstream_dist}bp)\n{base_name}', fontsize=14)
            plt.xlabel('Absolute Count', fontsize=12)
            plt.ylabel(f'{kmer_size}-mer', fontsize=12)
            plt.tight_layout()
            
            output_path = os.path.join(output_folder, f'motif_plot_{kmer_size}mer_{base_name}.png')
            plt.savefig(output_path, dpi=150, bbox_inches='tight')
            plt.close()
            
            print(f"  ✅ Motifs Chart saved in: {output_path}")
            status_queue.put(("image_ready", (output_path, f"Motifs {kmer_size}-mer - {base_name}")))

        except Exception as e_file:
            print(f"  ❌ General error processing motifs in {base_name}: {e_file}")