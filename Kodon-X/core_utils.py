import os
import numpy as np
import pandas as pd
from collections import Counter, defaultdict
from Bio import SeqIO
from constants import AA_CODON_MAPS, GENETIC_CODE_TABLES

def _apply_gene_filter(feature, gene_list):
    if not gene_list:
        return True 
        
    locus_tag = feature.qualifiers.get("locus_tag", [""])[0]
    gene_name = feature.qualifiers.get("gene", [""])[0]
    
    if locus_tag in gene_list or gene_name in gene_list:
        return True
        
    return False

def count_file_codons(file_path, status_queue=None, gene_list=None):
    print(f"  Reading all CDS (genes) from {os.path.basename(file_path)}...")
    codon_counts = Counter()
    total_cds_read = 0
    
    try:
        for record in SeqIO.parse(file_path, "genbank"):
            for feature in record.features:
                if feature.type != "CDS": continue
                
                if not _apply_gene_filter(feature, gene_list):
                    continue
                
                total_cds_read += 1
                seq_cds = feature.extract(record.seq)
                
                if feature.location.strand == -1:
                    seq_cds = seq_cds.reverse_complement()
                    
                codon_start = int(feature.qualifiers.get("codon_start", ["1"])[0])
                seq_cds = seq_cds[codon_start - 1:]
                
                for i in range(0, len(seq_cds) - 2, 3):
                    codon = str(seq_cds[i:i+3]).upper()
                    if len(codon) == 3 and 'N' not in codon and all(b in 'ATGC' for b in codon):
                        codon_counts[codon] += 1
                        
        print(f"  Total of {total_cds_read} CDS analyzed.")
        print(f"  Total of {sum(codon_counts.values())} codons counted.")
        
        if status_queue:
            status_queue.put(("message", f"Counted {sum(codon_counts.values())} codons in {os.path.basename(file_path)}"))
            
        return codon_counts
        
    except Exception as e:
        print(f"  ❌ Error counting codons in file {file_path}: {e}")
        if status_queue:
            status_queue.put(("message", f"Error reading {os.path.basename(file_path)}"))
        return None

def calculate_rscu(codon_counts, genetic_code_id=1):
    rscu_values = {}
    aa_codon_map = AA_CODON_MAPS.get(genetic_code_id, AA_CODON_MAPS[1]) 
    
    for aa, synonymous_codons in aa_codon_map.items():
        if len(synonymous_codons) == 1:
            if codon_counts[synonymous_codons[0]] > 0:
                rscu_values[synonymous_codons[0]] = 1.0
            else:
                rscu_values[synonymous_codons[0]] = 0.0
            continue
            
        total_for_aa = sum(codon_counts[c] for c in synonymous_codons)
        n_i = len(synonymous_codons)
        
        for codon in synonymous_codons:
            if total_for_aa == 0:
                rscu_values[codon] = 0.0
            else:
                rscu_values[codon] = (codon_counts[codon] * n_i) / total_for_aa
                
    return rscu_values

def calculate_gc12(codon_counts, genetic_code_id=1): 
    genetic_code = GENETIC_CODE_TABLES.get(genetic_code_id, GENETIC_CODE_TABLES[1])
    total_gc12 = 0
    total_codons = 0
    
    for codon, count in codon_counts.items():
        if len(codon) == 3 and codon in genetic_code:
            if genetic_code[codon] != '*':  
                total_codons += count
                if codon[0] in ['G', 'C']:
                    total_gc12 += count
                if codon[1] in ['G', 'C']:
                    total_gc12 += count
    
    return (total_gc12 / (total_codons * 2) * 100) if total_codons > 0 else 0

def calculate_enc(codon_counts, genetic_code_id=1):
    aa_codon_map = AA_CODON_MAPS.get(genetic_code_id, AA_CODON_MAPS[1])
    
    fold_groups = {2: [], 3: [], 4: [], 6: []}
    for aa, codons in aa_codon_map.items():
        if aa == '*' or len(codons) == 1:
            continue
        n_syn = len(codons)
        if n_syn in fold_groups:
            fold_groups[n_syn].append(aa)
            
    F_values = {2: 0, 3: 0, 4: 0, 6: 0}
    counts_per_group = {2: 0, 3: 0, 4: 0, 6: 0}
    
    for n_fold, aa_list in fold_groups.items():
        if not aa_list:
            continue
            
        sum_F_group = 0
        total_codons_group = 0
        
        for aa in aa_list:
            syn_codons = aa_codon_map[aa]
            total_for_aa = sum(codon_counts.get(c, 0) for c in syn_codons)
            
            if total_for_aa > 0:
                total_codons_group += total_for_aa
                sum_sq_p = sum((codon_counts.get(c, 0) / total_for_aa) ** 2 for c in syn_codons)
                F_aa = (total_for_aa * sum_sq_p - 1) / (total_for_aa - 1) if total_for_aa > 1 else 1
                sum_F_group += F_aa
                
        if len(aa_list) > 0 and total_codons_group > 0:
            F_values[n_fold] = sum_F_group / len(aa_list)
            counts_per_group[n_fold] = total_codons_group
            
    enc = 2.0  
    
    if F_values[2] > 0 and counts_per_group[2] > 0:
        enc += 9.0 / F_values[2]
    else:
        enc += 9.0 
        
    if F_values[3] > 0 and counts_per_group[3] > 0:
        enc += 1.0 / F_values[3]
    else:
        enc += 1.0
        
    if F_values[4] > 0 and counts_per_group[4] > 0:
        enc += 5.0 / F_values[4]
    else:
        enc += 5.0
        
    if F_values[6] > 0 and counts_per_group[6] > 0:
        enc += 3.0 / F_values[6]
    else:
        enc += 3.0
        
    return min(enc, 61.0) 

def calculate_gc3(codon_counts, genetic_code_id=1):
    genetic_code = GENETIC_CODE_TABLES.get(genetic_code_id, GENETIC_CODE_TABLES[1])
    total_gc3 = 0
    total_codons = 0
    
    for codon, count in codon_counts.items():
        if len(codon) == 3 and codon in genetic_code:
            if genetic_code[codon] != '*':  
                total_codons += count
                if codon[2] in ['G', 'C']:
                    total_gc3 += count
    
    return (total_gc3 / total_codons * 100) if total_codons > 0 else 0

def calculate_cai(codon_counts, rscu_values, genetic_code_id=1):
    aa_codon_map = AA_CODON_MAPS.get(genetic_code_id, AA_CODON_MAPS[1])
    
    w_values = {}
    for aa, codons in aa_codon_map.items():
        if aa == '*' or len(codons) == 1:
            continue
            
        max_rscu = 0
        for codon in codons:
            max_rscu = max(max_rscu, rscu_values.get(codon, 0))
            
        if max_rscu > 0:
            for codon in codons:
                w_values[codon] = rscu_values.get(codon, 0) / max_rscu
        else:
            for codon in codons:
                w_values[codon] = 1.0 
    
    log_w_sum = 0
    total_codons = 0
    
    for codon, count in codon_counts.items():
        if codon in w_values:
            if w_values[codon] > 0:
                log_w_sum += count * np.log(w_values[codon])
            total_codons += count
    
    if total_codons > 0:
        cai = np.exp(log_w_sum / total_codons)
        return cai
    else:
        return 0

def detect_optimal_rare_codons(rscu_values, genetic_code_id=1):
    aa_codon_map = AA_CODON_MAPS.get(genetic_code_id, AA_CODON_MAPS[1])
    optimal_codons = {}
    rare_codons = {}
    
    for aa, codons in aa_codon_map.items():
        if aa == '*' or len(codons) <= 1:
            continue
        
        codon_rscu = [(codon, rscu_values.get(codon, 0)) for codon in codons]
        
        codon_rscu.sort(key=lambda x: x[1], reverse=True)
        if codon_rscu[0][1] > 1.2:
             optimal_codons[aa] = codon_rscu[0][0]
        
        codon_rscu.sort(key=lambda x: x[1])
        if codon_rscu[0][1] < 0.8:
            rare_codons[aa] = codon_rscu[0][0]
    
    return optimal_codons, rare_codons

def process_genomes_for_bias_analysis(file_list, genetic_code_id, status_queue, gene_list=None):
    all_data = {}
    total_files = len(file_list)
    
    print(f"Starting codon bias processing for {total_files} files...")
    
    for i, file_path in enumerate(file_list):
        base_name = os.path.basename(file_path).split('.')[0]
        print(f"\nProcessing file {i+1}/{total_files}: {base_name}")
        status_queue.put(("progress", int(10 + (i / total_files) * 80)))
        
        counts = count_file_codons(file_path, status_queue, gene_list)
        if not counts or sum(counts.values()) == 0:
            print(f"  Warning: Skipping {base_name} (no codons or error).")
            continue
            
        rscu = calculate_rscu(counts, genetic_code_id)
        enc = calculate_enc(counts, genetic_code_id)
        gc3 = calculate_gc3(counts, genetic_code_id)
        cai = calculate_cai(counts, rscu, genetic_code_id)
        optimal, rare = detect_optimal_rare_codons(rscu, genetic_code_id)
        gc12 = calculate_gc12(counts, genetic_code_id)
        
        all_data[base_name] = {
            'counts': counts, 'rscu': rscu, 'enc': enc, 'gc3': gc3,
            'cai': cai, 'optimal': optimal, 'rare': rare, 'gc12': gc12
        }
        
        print(f"  📊 {base_name}: ENC={enc:.2f}, GC3={gc3:.2f}%, GC12={gc12:.2f}%, CAI={cai:.3f}")

    status_queue.put(("progress", 90))
    
    if len(all_data) < 1:
        print("  ❌ No file could be processed for bias analysis.")
        return None

    print("\n✅ Bias processing completed for all files.")
    return all_data

def extract_cds_sequences(file_list, status_queue, gene_list=None):
    all_seqs_by_species = defaultdict(list)
    total_files = len(file_list)
    
    print("  Extracting CDS sequences...")
    
    for i, full_path in enumerate(file_list):
        base_name = os.path.basename(full_path).split('.')[0]
        print(f"  Reading CDS of {base_name}...")
        status_queue.put(("message", f"Reading CDS of {base_name}..."))
        
        try:
            for record in SeqIO.parse(full_path, "genbank"):
                for feature in record.features:
                    if feature.type != "CDS": continue
                    
                    if not _apply_gene_filter(feature, gene_list):
                        continue
                    
                    seq_cds = feature.extract(record.seq).upper()
                    
                    if feature.location.strand == -1:
                        seq_cds = seq_cds.reverse_complement()
                        
                    codon_start = int(feature.qualifiers.get("codon_start", ["1"])[0])
                    seq_cds = seq_cds[codon_start - 1:]
                    
                    real_size = (len(seq_cds) // 3) * 3
                    seq_cds = seq_cds[:real_size]
                    
                    if 'N' not in seq_cds and len(seq_cds) >= 6: 
                        all_seqs_by_species[base_name].append(str(seq_cds))
                        
        except Exception as e:
            print(f"  ❌ Error extracting CDS from file {base_name}: {e}")
            
        status_queue.put(("progress", int(10 + (i / total_files) * 20)))
        
    print(f"  ✅ CDS extraction completed. {len(all_seqs_by_species)} species processed.")
    return all_seqs_by_species

def _count_sequence_codons(seq_str, genetic_code_id=1):
    genetic_code = GENETIC_CODE_TABLES.get(genetic_code_id, GENETIC_CODE_TABLES[1])
    codon_counts = Counter()
    seq_str = seq_str.upper()
    
    for i in range(0, len(seq_str) - 2, 3):
        codon = seq_str[i:i+3]
        if len(codon) == 3 and codon in genetic_code:
            codon_counts[codon] += 1
            
    return codon_counts

def _get_optimal_codon_map(rscu_values, genetic_code_id=1):
    aa_codon_map = AA_CODON_MAPS.get(genetic_code_id, AA_CODON_MAPS[1])
    optimal_map = {}
    
    for aa, codons in aa_codon_map.items():
        if not codons: continue
        
        best_codon = codons[0]
        max_rscu = -1.0
        
        for codon in codons:
            rscu = rscu_values.get(codon, 0.0)
            if rscu > max_rscu:
                max_rscu = rscu
                best_codon = codon
                
        optimal_map[aa] = best_codon
        
    return optimal_map

def _get_codon_frequency_rank(codon_counts, genetic_code_id=1):
    aa_codon_map = AA_CODON_MAPS.get(genetic_code_id, AA_CODON_MAPS[1])
    rank_map = {}
    
    for aa, codons in aa_codon_map.items():
        if not codons or len(codons) == 1:
            rank_map[aa] = codons 
            continue
            
        codon_pairs = [(codon, codon_counts.get(codon, 0)) for codon in codons]
        codon_pairs.sort(key=lambda x: x[1], reverse=True)
        rank_map[aa] = [codon for codon, count in codon_pairs]
        
    return rank_map

def get_w_reference_table(rscu_values, genetic_code_id=1):
    aa_codon_map = AA_CODON_MAPS.get(genetic_code_id, AA_CODON_MAPS[1])
    w_values = {}
    
    for aa, codons in aa_codon_map.items():
        if aa == '*' or len(codons) == 1:
            if codons:
                w_values[codons[0]] = 1.0
            continue
            
        max_rscu = 0
        for codon in codons:
            max_rscu = max(max_rscu, rscu_values.get(codon, 0))
            
        if max_rscu > 0:
            for codon in codons:
                w_values[codon] = rscu_values.get(codon, 0) / max_rscu
        else:
            for codon in codons:
                w_values[codon] = 1.0 
    
    return w_values

def calculate_metrics_per_gene(file_path, w_reference_table, genetic_code_id, gene_list):
    print(f"  Calculating metrics per gene for {len(gene_list) if gene_list else 'all'} genes...")
    results = []
    
    try:
        for record in SeqIO.parse(file_path, "genbank"):
            for feature in record.features:
                if feature.type != "CDS": continue
                
                if not _apply_gene_filter(feature, gene_list):
                    continue
                
                identifier = "N/A"
                if "locus_tag" in feature.qualifiers: identifier = feature.qualifiers["locus_tag"][0]
                elif "gene" in feature.qualifiers: identifier = feature.qualifiers["gene"][0]
                elif "protein_id" in feature.qualifiers: identifier = feature.qualifiers["protein_id"][0]
                
                seq_cds_str = ""
                try:
                    seq_cds = feature.extract(record.seq).upper()
                    if feature.location.strand == -1:
                        seq_cds = seq_cds.reverse_complement()
                    codon_start = int(feature.qualifiers.get("codon_start", ["1"])[0])
                    seq_cds_str = str(seq_cds[codon_start - 1:])
                except Exception as e_seq:
                    print(f"    Warning: Skipping gene {identifier} (extraction error): {e_seq}")
                    continue
                
                counts_gene = _count_sequence_codons(seq_cds_str, genetic_code_id)
                if sum(counts_gene.values()) == 0:
                    continue
                
                enc_gene = calculate_enc(counts_gene, genetic_code_id)
                gc3_gene = calculate_gc3(counts_gene, genetic_code_id)
                
                log_w_sum = 0
                total_codons = 0
                for codon, count in counts_gene.items():
                    if codon in w_reference_table:
                        if w_reference_table[codon] > 0:
                            log_w_sum += count * np.log(w_reference_table[codon])
                        total_codons += count
                
                cai_gene = np.exp(log_w_sum / total_codons) if total_codons > 0 else 0
                
                results.append({
                    'gene': identifier,
                    'enc': enc_gene,
                    'gc3': gc3_gene,
                    'cai': cai_gene,
                    'counts': counts_gene 
                })
                
    except Exception as e:
        print(f"  ❌ Error calculating metrics per gene: {e}")
        
    return results
