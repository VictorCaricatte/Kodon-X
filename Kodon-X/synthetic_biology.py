from constants import GENETIC_CODE_TABLES
from core_utils import _get_optimal_codon_map, _get_codon_frequency_rank, _count_sequence_codons

def optimize_codon_sequence(input_seq, host_rscu, genetic_code_id=1):
    print("  Starting Codon Optimization (Maximization)...")
    genetic_code = GENETIC_CODE_TABLES.get(genetic_code_id, GENETIC_CODE_TABLES[1])
    
    optimal_map = _get_optimal_codon_map(host_rscu, genetic_code_id)
    
    input_seq = input_seq.upper()
    aa_sequence = []
    stop_codon = ""
    for i in range(0, len(input_seq) - 2, 3):
        codon = input_seq[i:i+3]
        if codon not in genetic_code:
            aa_sequence.append('?') 
            continue
            
        aa = genetic_code[codon]
        if aa == '*':
            stop_codon = optimal_map.get('*', 'TAA') 
            break
        aa_sequence.append(aa)

    optimized_seq = []
    for aa in aa_sequence:
        if aa == '?':
            optimized_seq.append('NNN') 
        else:
            optimized_seq.append(optimal_map.get(aa, 'NNN'))
            
    optimized_seq.append(stop_codon)
    
    print("  ✅ Optimization Completed.")
    return "".join(optimized_seq)

def harmonize_codon_sequence(input_seq, host_counts, genetic_code_id=1):
    print("  Starting Codon Harmonization (Ranking)...")
    genetic_code = GENETIC_CODE_TABLES.get(genetic_code_id, GENETIC_CODE_TABLES[1])
    
    host_rank_map = _get_codon_frequency_rank(host_counts, genetic_code_id)
    
    input_counts = _count_sequence_codons(input_seq, genetic_code_id)
    input_rank_map = _get_codon_frequency_rank(input_counts, genetic_code_id)
    
    harmonization_map = {}
    for aa, input_ranked_codons in input_rank_map.items():
        if aa not in host_rank_map: continue
        
        host_ranked_codons = host_rank_map[aa]
        
        for i, input_codon in enumerate(input_ranked_codons):
            host_index = min(i, len(host_ranked_codons) - 1)
            harmonization_map[input_codon] = host_ranked_codons[host_index]
            
    input_seq = input_seq.upper()
    harmonized_seq = []
    
    for i in range(0, len(input_seq) - 2, 3):
        codon = input_seq[i:i+3]
        if codon in harmonization_map:
            harmonized_seq.append(harmonization_map[codon])
        else:
            harmonized_seq.append(codon)
            
    print("  ✅ Harmonization Completed.")
    return "".join(harmonized_seq)
