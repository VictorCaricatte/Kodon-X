# constants.py
import matplotlib.pyplot as plt
from collections import defaultdict

# ######################################################################
# --- DEFINIÇÕES DE CÓDONS E TABELAS GENÉTICAS (Constantes) ---
# ######################################################################

# Dicionário de dicionários para as tabelas genéticas
GENETIC_CODE_TABLES = {
    1: { # Tabela 1: Padrão
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*', 'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    },
    2: { # Tabela 2: Mitocondrial de Vertebrados
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*', 'TGT': 'C', 'TGC': 'C', 'TGA': 'W', 'TGG': 'W',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'M', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': '*', 'AGG': '*',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    },
    4: { # Tabela 4: Mitocondrial de Mofos, Protozoários e Celenterados
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*', 'TGT': 'C', 'TGC': 'C', 'TGA': 'W', 'TGG': 'W',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    },
    11: { # Tabela 11: Plasto de Bactérias e Plantas
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

# Gera mapas reversos (AA -> [Códons]) para cada tabela
AA_CODON_MAPS = {
    table_id: defaultdict(list) for table_id in GENETIC_CODE_TABLES
}
for table_id, codon_map in GENETIC_CODE_TABLES.items():
    for codon, aa in codon_map.items():
        AA_CODON_MAPS[table_id][aa].append(codon)

# Ordem para o heatmap da Opção 3
CODON_GRID_ORDER = [
    ['TTT', 'TTC', 'TTA', 'TTG', 'TCT', 'TCC', 'TCA', 'TCG', 'TAT', 'TAC', 'TAA', 'TAG', 'TGT', 'TGC', 'TGA', 'TGG'],
    ['CTT', 'CTC', 'CTA', 'CTG', 'CCT', 'CCC', 'CCA', 'CCG', 'CAT', 'CAC', 'CAA', 'CAG', 'CGT', 'CGC', 'CGA', 'CGG'],
    ['ATT', 'ATC', 'ATA', 'ATG', 'ACT', 'ACC', 'ACA', 'ACG', 'AAT', 'AAC', 'AAA', 'AAG', 'AGT', 'AGC', 'AGA', 'AGG'],
    ['GTT', 'GTC', 'GTA', 'GTG', 'GCT', 'GCC', 'GCA', 'GCG', 'GAT', 'GAC', 'GAA', 'GAG', 'GGT', 'GGC', 'GGA', 'GGG']
]
ALL_CODONS_SORTED = sorted(GENETIC_CODE_TABLES[1].keys())

# Paleta de cores para os histogramas (para consistência)
CODON_COLORS = {codon: color for codon, color in zip(ALL_CODONS_SORTED, plt.cm.get_cmap('tab20').colors * 4)}

# ######################################################################
# --- CONSTANTES FÍSICO-QUÍMICAS (Para Novas Análises) ---
# ######################################################################

# Escala de Hidropaticidade de Kyte & Doolittle (para GRAVY)
KYTE_DOOLITTLE_HYDROPATHY = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
    'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
    'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2,
    '*': 0.0 # Stop codon
}

# Aromaticidade (para Aromo)
AROMATICITY = {
    'A': 0, 'R': 0, 'N': 0, 'D': 0, 'C': 0,
    'Q': 0, 'E': 0, 'G': 0, 'H': 0, 'I': 0,
    'L': 0, 'K': 0, 'M': 0, 'F': 1, 'P': 0,
    'S': 0, 'T': 0, 'W': 1, 'Y': 1, 'V': 0,
    '*': 0 # Stop codon
}

# ######################################################################
# --- CONSTANTES DE PAREAMENTO tAI (WOBBLE) ---
# Baseado em dos Reis et al. (2004) Nucleic Acids Research
# e implementações canônicas (ex: pacote tAI em R)
# ######################################################################

# Estes são os 9 valores 's' (penalidades de pareamento) otimizados
# para os diferentes super-reinos. Um valor 's' de 0.0 é um pareamento
# perfeito (Watson-Crick), e 1.0 é nenhum pareamento.

# Valores 's' otimizados para Bactérias
# (Ref: Otimizado para E. coli no pacote tAI de dos Reis)
WOBBLE_S_VALUES_BACTERIA = {
    's_GU': 0.69,  # Anticódon G lendo C(0.0) ou U(s_GU)
    's_AU': 1.0,   # Anticódon A (raro) lendo U(0.0) - geralmente modificado para I
    's_UG': 0.69,  # Anticódon U lendo A(0.0) ou G(s_UG)
    's_UU': 0.51,  # Anticódon U lendo A(0.0) ou U(s_UU)
    's_CG': 1.0,   # Anticódon C lendo G(0.0) - sem wobble
    's_IA': 0.89,  # Anticódon I (de A) lendo A(s_IA)
    's_IC': 0.45,  # Anticódon I (de A) lendo C(s_IC)
    's_IU': 0.0,   # Anticódon I (de A) lendo U(s_IU) - pareamento mais forte
    's_UA': 1.0    # Anticódon U lendo A(0.0) - sem wobble com A (no modelo de 9 parâmetros)
}

# Valores 's' otimizados para Eucariotos
# (Ref: Otimizado para S. cerevisiae no pacote tAI de dos Reis)
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
# Nota: Archaea também tem seu próprio conjunto, mas Bactéria e Eukarya são os mais comuns.

def _build_wobble_matrix(s_values):
    """
    Função auxiliar para construir a matriz de pareamento S_ij completa
    a partir de um vetor de 9 parâmetros 's'.
    """
    s_GU = s_values['s_GU']
    s_AU = s_values['s_AU']
    s_UG = s_values['s_UG']
    s_UU = s_values['s_UU']
    s_CG = s_values['s_CG']
    s_IA = s_values['s_IA']
    s_IC = s_values['s_IC']
    s_IU = s_values['s_IU']
    s_UA = s_values['s_UA']
    
    # Chave: (Base do Anticódon, Base do Códon) -> Valor S_ij
    matrix = {
        # 1. Anticódon G (lê C ou U)
        ('G', 'C'): 0.0,   # Perfeito (Watson-Crick)
        ('G', 'U'): s_GU,
        ('G', 'A'): 1.0,   # Sem pareamento
        ('G', 'G'): 1.0,   # Sem pareamento
        
        # 2. Anticódon A (lê U) - Raro, quase sempre modificado para Inosina (I)
        ('A', 'U'): 0.0,
        ('A', 'A'): 1.0,
        ('A', 'C'): 1.0,
        ('A', 'G'): s_AU,  # No modelo de 9 parâmetros, A-G é permitido (embora 1.0 em E.coli)
        
        # 3. Anticódon U (lê A, G, ou U)
        ('U', 'A'): 0.0,   # Perfeito (Watson-Crick)
        ('U', 'G'): s_UG,
        ('U', 'U'): s_UU,
        ('U', 'C'): 1.0,   # Sem pareamento
        
        # 4. Anticódon C (lê G)
        ('C', 'G'): 0.0,   # Perfeito (Watson-Crick)
        ('C', 'A'): 1.0,
        ('C', 'U'): 1.0,
        ('C', 'C'): 1.0,
        
        # 5. Anticódon I (Inosina, modificado de A) (lê A, C, ou U)
        ('I', 'U'): s_IU,  # Mais forte
        ('I', 'C'): s_IC,
        ('I', 'A'): s_IA,
        ('I', 'G'): 1.0,   # Sem pareamento
    }
    return matrix

# Gerar as matrizes canônicas para exportação
WOBBLE_MATRIX_BACTERIA = _build_wobble_matrix(WOBBLE_S_VALUES_BACTERIA)
WOBBLE_MATRIX_EUKARYA = _build_wobble_matrix(WOBBLE_S_VALUES_EUKARYA)

# Mapeamento de modificação da 1ª base do anticódon (posição 34)
# Isso informa à nossa função qual base funcional usar (ex: A -> I)
ANTICODON_MODIFICATION_MAP = {
    # Chave: base 5' (posição 1) do anticódon
    # Valor: base modificada funcional para pareamento
    'A': 'I', # Adenina na pos 1 (wobble) quase sempre vira Inosina (I)
    'G': 'G',
    'C': 'C',
    'U': 'U',
    'I': 'I'  # Se já for Inosina
}
