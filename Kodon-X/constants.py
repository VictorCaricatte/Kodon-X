import matplotlib.pyplot as plt
from collections import defaultdict

# Reference (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?chapter=tgencodes#SG16)

GENETIC_CODE_TABLES = {
    1: { 
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*', 'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    },
    2: { 
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*', 'TGT': 'C', 'TGC': 'C', 'TGA': 'W', 'TGG': 'W',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'M', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': '*', 'AGG': '*',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    },
    4: { 
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*', 'TGT': 'C', 'TGC': 'C', 'TGA': 'W', 'TGG': 'W',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    },
    11: { 
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*', 'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    },
}

AA_CODON_MAPS = {
    table_id: defaultdict(list) for table_id in GENETIC_CODE_TABLES
}
for table_id, codon_map in GENETIC_CODE_TABLES.items():
    for codon, aa in codon_map.items():
        AA_CODON_MAPS[table_id][aa].append(codon)

CODON_GRID_ORDER = [
    ['TTT', 'TTC', 'TTA', 'TTG', 'TCT', 'TCC', 'TCA', 'TCG', 'TAT', 'TAC', 'TAA', 'TAG', 'TGT', 'TGC', 'TGA', 'TGG'],
    ['CTT', 'CTC', 'CTA', 'CTG', 'CCT', 'CCC', 'CCA', 'CCG', 'CAT', 'CAC', 'CAA', 'CAG', 'CGT', 'CGC', 'CGA', 'CGG'],
    ['ATT', 'ATC', 'ATA', 'ATG', 'ACT', 'ACC', 'ACA', 'ACG', 'AAT', 'AAC', 'AAA', 'AAG', 'AGT', 'AGC', 'AGA', 'AGG'],
    ['GTT', 'GTC', 'GTA', 'GTG', 'GCT', 'GCC', 'GCA', 'GCG', 'GAT', 'GAC', 'GAA', 'GAG', 'GGT', 'GGC', 'GGA', 'GGG']
]
ALL_CODONS_SORTED = sorted(GENETIC_CODE_TABLES[1].keys())

CODON_COLORS = {codon: color for codon, color in zip(ALL_CODONS_SORTED, plt.cm.get_cmap('tab20').colors * 4)}

KYTE_DOOLITTLE_HYDROPATHY = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
    'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
    'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2,
    '*': 0.0 
}

AROMATICITY = {
    'A': 0, 'R': 0, 'N': 0, 'D': 0, 'C': 0,
    'Q': 0, 'E': 0, 'G': 0, 'H': 0, 'I': 0,
    'L': 0, 'K': 0, 'M': 0, 'F': 1, 'P': 0,
    'S': 0, 'T': 0, 'W': 1, 'Y': 1, 'V': 0,
    '*': 0 
}

WOBBLE_S_VALUES_BACTERIA = {
    's_GU': 0.69,  
    's_AU': 1.0,   
    's_UG': 0.69,  
    's_UU': 0.51,  
    's_CG': 1.0,   
    's_IA': 0.89,  
    's_IC': 0.45,  
    's_IU': 0.0,   
    's_UA': 1.0    
}

WOBBLE_S_VALUES_EUKARYA = {
    's_GU': 0.79,
    's_AU': 1.0,
    's_UG': 0.63,
    's_UU': 0.68,
    's_CG': 1.0,
    's_IA': 0.91,
    's_IC': 0.47,
    's_IU': 0.0,
    's_UA': 1.0
}

def _build_wobble_matrix(s_values):
    s_GU = s_values['s_GU']
    s_AU = s_values['s_AU']
    s_UG = s_values['s_UG']
    s_UU = s_values['s_UU']
    s_CG = s_values['s_CG']
    s_IA = s_values['s_IA']
    s_IC = s_values['s_IC']
    s_IU = s_values['s_IU']
    s_UA = s_values['s_UA']
    
    matrix = {
        ('G', 'C'): 0.0,   
        ('G', 'U'): s_GU,
        ('G', 'A'): 1.0,   
        ('G', 'G'): 1.0,   
        ('A', 'U'): 0.0,
        ('A', 'A'): 1.0,
        ('A', 'C'): 1.0,
        ('A', 'G'): s_AU,  
        ('U', 'A'): 0.0,   
        ('U', 'G'): s_UG,
        ('U', 'U'): s_UU,
        ('U', 'C'): 1.0,   
        ('C', 'G'): 0.0,   
        ('C', 'A'): 1.0,
        ('C', 'U'): 1.0,
        ('C', 'C'): 1.0,
        ('I', 'U'): s_IU,  
        ('I', 'C'): s_IC,
        ('I', 'A'): s_IA,
        ('I', 'G'): 1.0,   
    }
    return matrix

WOBBLE_MATRIX_BACTERIA = _build_wobble_matrix(WOBBLE_S_VALUES_BACTERIA)
WOBBLE_MATRIX_EUKARYA = _build_wobble_matrix(WOBBLE_S_VALUES_EUKARYA)

ANTICODON_MODIFICATION_MAP = {
    'A': 'I', 
    'G': 'G',
    'C': 'C',
    'U': 'U',
    'I': 'I'  
}