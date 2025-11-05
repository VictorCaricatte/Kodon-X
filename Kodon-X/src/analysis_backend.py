# analysis_backend.py
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from collections import Counter
from scipy.cluster import hierarchy as sch
from scipy.stats import pearsonr
from scipy.spatial.distance import pdist, squareform
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.stats import linregress
from collections import defaultdict
from itertools import product
import re
from Bio.Seq import Seq
# --- NOVAS IMPORTAÇÕES ---
from scipy.stats import mannwhitneyu, chi2_contingency, spearmanr

# Tenta importar a biblioteca de MFE (ViennaRNA)
try:
    import RNA
    VIENNARNA_DISPONIVEL = True
except ImportError:
    VIENNARNA_DISPONIVEL = False
    print("Aviso de Backend: ViennaRNA não encontrado. Análise de MFE (17) desabilitada.")

# --- IMPORTAÇÕES DE CONSTANTES ATUALIZADAS ---
from constants import (
    KYTE_DOOLITTLE_HYDROPATHY, 
    AROMATICITY,
    WOBBLE_MATRIX_BACTERIA,  # <--- NOVA
    WOBBLE_MATRIX_EUKARYA, # <--- NOVA
    ANTICODON_MODIFICATION_MAP # <--- NOVA
)
# Importa as constantes do nosso arquivo local
from constants import (
    AA_CODON_MAPS, 
    GENETIC_CODE_TABLES, 
    CODON_GRID_ORDER, 
    ALL_CODONS_SORTED
)
# --- FIM DAS IMPORTAÇÕES ATUALIZADAS ---

def _apply_gene_filter(feature, gene_list):
    """Função auxiliar para verificar se uma feature passa no filtro de gene_list."""
    if not gene_list:
        return True # Se não há lista, passa
        
    locus_tag = feature.qualifiers.get("locus_tag", [""])[0]
    gene_name = feature.qualifiers.get("gene", [""])[0]
    
    if locus_tag in gene_list or gene_name in gene_list:
        return True
        
    return False

# ######################################################################
# --- BLOCO DE ANÁLISE (BACKEND) 
# ######################################################################

# --- Funções da Análise 1 ---

def processar_gbk_agregado(lista_arquivos, status_queue):
    """(Análise 1) Processa .gbk, calculando estatísticas AGREGADAS."""
    if not lista_arquivos:
        print("  ❌ Nenhum arquivo .gbk ou .gb fornecido.")
        return pd.DataFrame()
    
    status_queue.put(("progress", 20))
    resultados_gerais = []
    
    for i, caminho_completo in enumerate(lista_arquivos):
        arquivo = os.path.basename(caminho_completo)
        print(f"  Processando estatísticas de {arquivo}...")
        try:
            registros = list(SeqIO.parse(caminho_completo, "genbank"))
            num_records = len(registros)
            if num_records == 0:
                print(f"  ❌ Erro: Arquivo {arquivo} está vazio ou corrompido.")
                continue
            
            comprimento_total = 0
            contagem_total_gc = 0
            for record in registros:
                seq = record.seq.upper()
                comprimento_total += len(seq)
                contagem_total_gc += seq.count('G') + seq.count('C')
                
            if comprimento_total == 0:
                conteudo_gc_total = 0.0
            else:
                conteudo_gc_total = (contagem_total_gc / comprimento_total) * 100
                
            resultados_gerais.append({
                'Arquivo': arquivo, 'Num_Contigs_Total': num_records,
                'Comprimento_Total': comprimento_total, 'Conteudo_GC_%_Total': f'{conteudo_gc_total:.2f}'
            })
        except Exception as e:
            print(f"  ❌ Erro ao processar o arquivo {arquivo}: {e}")
            
    status_queue.put(("progress", 50))
    return pd.DataFrame(resultados_gerais)

def analisar_gbk_cds(lista_arquivos, filtro, status_queue, gene_list=None):
    """(Análise 1) Analisa CDS (Início/Tamanho) dos arquivos .gbk."""
    if not lista_arquivos:
        print("  ❌ Nenhum arquivo .gbk ou .gb fornecido.")
        return pd.DataFrame()
        
    resultados_finais = []
    total_files = len(lista_arquivos)
    
    for i, caminho_completo in enumerate(lista_arquivos):
        arquivo = os.path.basename(caminho_completo)
        print(f"  Analisando CDS de {arquivo}...")
        try:
            registros = list(SeqIO.parse(caminho_completo, "genbank"))
            if not registros:
                print(f"  ❌ Erro: Arquivo {arquivo} está vazio ou corrompido.")
                continue

            for record in registros:
                for feature in record.features:
                    if feature.type != "CDS": continue
                    
                    # --- APLICAR FILTRO ---
                    if not _apply_gene_filter(feature, gene_list):
                        continue
                    # --- FIM DO FILTRO ---
                    
                    identificador = "N/A"
                    if "locus_tag" in feature.qualifiers: identificador = feature.qualifiers["locus_tag"][0]
                    elif "protein_id" in feature.qualifiers: identificador = feature.qualifiers["protein_id"][0]
                    elif "gene" in feature.qualifiers: identificador = feature.qualifiers["gene"][0]
                    
                    seq_cds = feature.extract(record.seq)
                    if feature.location.strand == -1: seq_cds = seq_cds.reverse_complement()
                    
                    codon_start = int(feature.qualifiers.get("codon_start", ["1"])[0])
                    seq_cds = seq_cds[codon_start - 1:]
                    tamanho_real = (len(seq_cds) // 3) * 3
                    seq_cds = seq_cds[:tamanho_real]
                    
                    if len(seq_cds) < 3: continue
                    
                    codon_inicio = str(seq_cds[:3]).upper()
                    comeca_com_filtro = codon_inicio == filtro.upper()
                    
                    resultados_finais.append({
                        'Arquivo': arquivo, 'Registro': record.id, 'Identificador': identificador,
                        'Codon_Inicio': codon_inicio, 'Comeca_Com_Filtro': comeca_com_filtro,
                        'Codon_Start_Qualifier': codon_start, 'Tamanho_Real_ORF_nt': tamanho_real
                    })
        except Exception as e:
            print(f"  ❌ Erro ao analisar CDS do arquivo {arquivo}: {e}")
        
        status_queue.put(("progress", int(50 + (i / total_files) * 50)))
        
    return pd.DataFrame(resultados_finais)

# --- Função da Análise 2 ---

def listar_genes_do_arquivo(caminho_arquivo, status_queue, gene_list=None):
    """(Análise 2) Lista todos os genes, produtos e locus tags de um arquivo."""
    nome_arquivo = os.path.basename(caminho_arquivo)
    print(f"\nListando todos os genes/produtos de {nome_arquivo}...")
    resultados = []
    
    status_queue.put(("progress", 20))
    
    try:
        for i, record in enumerate(SeqIO.parse(caminho_arquivo, "genbank")):
            if i % 10 == 0: # Atualiza progresso a cada 10 contigs
                 status_queue.put(("message", f"Processando contig {record.id}..."))
                 
            for feature in record.features:
                if feature.type not in ["CDS", "tRNA", "rRNA", "gene"]:
                    continue
                
                # --- APLICAR FILTRO ---
                if not _apply_gene_filter(feature, gene_list):
                    continue
                # --- FIM DO FILTRO ---
                
                nome_gene = feature.qualifiers.get("gene", [""])[0]
                produto = feature.qualifiers.get("product", [""])[0]
                locus_tag = feature.qualifiers.get("locus_tag", ["N/A"])[0]
                
                if nome_gene or produto or locus_tag != "N/A":
                    resultados.append({
                        "Arquivo": nome_arquivo, "Contig": record.id,
                        "Locus_Tag": locus_tag, "Gene": nome_gene,
                        "Produto": produto, "Tipo": feature.type
                    })
                    
        status_queue.put(("progress", 80))
        
        if not resultados:
            print(f"  ❌ Nenhum gene ou produto encontrado em {nome_arquivo}.")
            return pd.DataFrame()
            
        print(f"  ✅ Encontrados {len(resultados)} genes/produtos no total.")
        return pd.DataFrame(resultados)
        
    except Exception as e:
        print(f"  ❌ Erro ao processar o arquivo {caminho_arquivo} para listar genes: {e}")
        return pd.DataFrame()

# --- Funções de Análise de Viés (Base para Análises 3 a 9) ---

def contar_codons_de_arquivo(caminho_arquivo, status_queue=None, gene_list=None):
    """Conta todos os códons de todos os CDS em um único arquivo .gbk."""
    print(f"  Lendo todos os CDS (genes) de {os.path.basename(caminho_arquivo)}...")
    codon_counts = Counter()
    total_cds_lidos = 0
    
    try:
        for record in SeqIO.parse(caminho_arquivo, "genbank"):
            for feature in record.features:
                if feature.type != "CDS": continue
                
                # --- APLICAR FILTRO ---
                if not _apply_gene_filter(feature, gene_list):
                    continue
                # --- FIM DO FILTRO ---
                
                total_cds_lidos += 1
                seq_cds = feature.extract(record.seq)
                
                if feature.location.strand == -1:
                    seq_cds = seq_cds.reverse_complement()
                    
                codon_start = int(feature.qualifiers.get("codon_start", ["1"])[0])
                seq_cds = seq_cds[codon_start - 1:]
                
                for i in range(0, len(seq_cds) - 2, 3):
                    codon = str(seq_cds[i:i+3]).upper()
                    if len(codon) == 3 and 'N' not in codon and all(b in 'ATGC' for b in codon):
                        codon_counts[codon] += 1
                        
        print(f"  Total de {total_cds_lidos} CDS analisados.")
        print(f"  Total de {sum(codon_counts.values())} códons contados.")
        
        if status_queue:
            status_queue.put(("message", f"Contados {sum(codon_counts.values())} códons em {os.path.basename(caminho_arquivo)}"))
            
        return codon_counts
        
    except Exception as e:
        print(f"  ❌ Erro ao contar códons no arquivo {caminho_arquivo}: {e}")
        if status_queue:
            status_queue.put(("message", f"Erro ao ler {os.path.basename(caminho_arquivo)}"))
        return None

def calcular_rscu(codon_counts, genetic_code_id=1):
    """Calcula o RSCU a partir de um dicionário de contagens de códons."""
    rscu_values = {}
    aa_codon_map = AA_CODON_MAPS.get(genetic_code_id, AA_CODON_MAPS[1]) # Fallback para tabela 1
    
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

# --- NOVAS FUNÇÕES DE ANÁLISE AVANÇADAS (Base) ---
def calcular_gc12(codon_counts, genetic_code_id=1): 
    """Calcula o conteúdo GC nas posições 1 e 2"""
    genetic_code = GENETIC_CODE_TABLES.get(genetic_code_id, GENETIC_CODE_TABLES[1])
    total_gc12 = 0
    total_codons = 0
    
    for codon, count in codon_counts.items():
        if len(codon) == 3 and codon in genetic_code:
            if genetic_code[codon] != '*':  # Excluir codons de parada
                total_codons += count
                if codon[0] in ['G', 'C']:
                    total_gc12 += count
                if codon[1] in ['G', 'C']:
                    total_gc12 += count
    
    # Cada códon tem 2 posições (1 e 2), então o total é total_codons * 2
    return (total_gc12 / (total_codons * 2) * 100) if total_codons > 0 else 0

def calcular_enc(codon_counts, genetic_code_id=1):
    """Calcula o Effective Number of Codons (ENC)"""
    aa_codon_map = AA_CODON_MAPS.get(genetic_code_id, AA_CODON_MAPS[1])
    
    # Agrupar AAs por número de sinônimos
    fold_groups = {2: [], 3: [], 4: [], 6: []}
    for aa, codons in aa_codon_map.items():
        if aa == '*' or len(codons) == 1:
            continue
        n_syn = len(codons)
        if n_syn in fold_groups:
            fold_groups[n_syn].append(aa)
            
    # Calcular F (homozigosidade) para cada grupo
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
            
    # Calcular ENC (Fórmula de Wright, 1990)
    enc = 2.0  # M (1 códon) + W (1 códon)
    
    if F_values[2] > 0 and counts_per_group[2] > 0:
        enc += 9.0 / F_values[2]
    else:
        enc += 9.0 # Valor máximo se não houver dados
        
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
        
    return min(enc, 61.0) # ENC máximo é 61

def calcular_gc3(codon_counts, genetic_code_id=1):
    """Calcula o conteúdo GC na terceira posição"""
    genetic_code = GENETIC_CODE_TABLES.get(genetic_code_id, GENETIC_CODE_TABLES[1])
    total_gc3 = 0
    total_codons = 0
    
    for codon, count in codon_counts.items():
        if len(codon) == 3 and codon in genetic_code:
            if genetic_code[codon] != '*':  # Excluir codons de parada
                total_codons += count
                if codon[2] in ['G', 'C']:
                    total_gc3 += count
    
    return (total_gc3 / total_codons * 100) if total_codons > 0 else 0

def calcular_cai(codon_counts, rscu_values, genetic_code_id=1):
    """Calcula o Codon Adaptation Index (CAI)"""
    aa_codon_map = AA_CODON_MAPS.get(genetic_code_id, AA_CODON_MAPS[1])
    
    # Calcular pesos relativos (w)
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
                w_values[codon] = 1.0 # Se não houver dados, peso é 1
    
    # Calcular CAI
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

def detectar_codons_otimos_raros(rscu_values, genetic_code_id=1):
    """Detecta códons ótimos e raros baseados no RSCU"""
    aa_codon_map = AA_CODON_MAPS.get(genetic_code_id, AA_CODON_MAPS[1])
    optimal_codons = {}
    rare_codons = {}
    
    for aa, codons in aa_codon_map.items():
        if aa == '*' or len(codons) <= 1:
            continue
        
        codon_rscu = [(codon, rscu_values.get(codon, 0)) for codon in codons]
        
        # Códons ótimos (RSCU > 1.2 e o maior do grupo)
        codon_rscu.sort(key=lambda x: x[1], reverse=True)
        if codon_rscu[0][1] > 1.2:
             optimal_codons[aa] = codon_rscu[0][0]
        
        # Códons raros (RSCU < 0.8)
        codon_rscu.sort(key=lambda x: x[1])
        if codon_rscu[0][1] < 0.8:
            rare_codons[aa] = codon_rscu[0][0]
    
    return optimal_codons, rare_codons

def analisar_composicao_nucleotideos(caminho_arquivo):
    """(Análise 8) Analisa a composição de nucleotídeos do genoma"""
    print(f"  Analisando composição de nucleotídeos de {os.path.basename(caminho_arquivo)}...")
    
    nt_counts = Counter()
    total_length = 0
    
    try:
        for record in SeqIO.parse(caminho_arquivo, "genbank"):
            seq = record.seq.upper()
            total_length += len(seq)
            for nt in seq:
                if nt in 'ATGC':
                    nt_counts[nt] += 1
        
        if total_length == 0:
            return None
        
        # Calcular percentagens
        composition = {}
        for nt in 'ATGC':
            composition[nt] = (nt_counts[nt] / total_length) * 100
        
        gc_content = (composition['G'] + composition['C'])
        
        print(f"  ✅ Composição analisada: GC={gc_content:.2f}%")
        
        return {
            'composition': composition,
            'gc_content': gc_content,
            'total_length': total_length
        }
    except Exception as e:
        print(f"  ❌ Erro ao analisar composição: {e}")
        return None


# --- FUNÇÃO MESTRA DE PROCESSAMENTO DE VIÉS (Refatorada) ---

def processar_genomas_para_analise_vies(lista_arquivos, genetic_code_id, status_queue, gene_list=None):
    """
    Função refatorada para processar múltiplos arquivos *uma única vez*.
    Calcula contagens, RSCU, ENC, GC3, GC12 e CAI e retorna um dicionário estruturado.
    """
    all_data = {}
    total_files = len(lista_arquivos)
    
    print(f"Iniciando processamento de viés de códons para {total_files} arquivos...")
    
    for i, caminho_arquivo in enumerate(lista_arquivos):
        nome_base = os.path.basename(caminho_arquivo).split('.')[0]
        print(f"\nProcessando arquivo {i+1}/{total_files}: {nome_base}")
        status_queue.put(("progress", int(10 + (i / total_files) * 80)))
        
        # 1. Contar Códons
        counts = contar_codons_de_arquivo(caminho_arquivo, status_queue, gene_list)
        if not counts or sum(counts.values()) == 0:
            print(f"  Aviso: Pulando {nome_base} (sem códons ou erro).")
            continue
            
        # 2. Calcular RSCU
        rscu = calcular_rscu(counts, genetic_code_id)
        
        # 3. Calcular ENC
        enc = calcular_enc(counts, genetic_code_id)
        
        # 4. Calcular GC3
        gc3 = calcular_gc3(counts, genetic_code_id)
        
        # 5. Calcular CAI
        cai = calcular_cai(counts, rscu, genetic_code_id)
        
        # 6. Detectar Códons Ótimos/Raros
        optimal, rare = detectar_codons_otimos_raros(rscu, genetic_code_id)
        
        # 7. Calcular GC12
        gc12 = calcular_gc12(counts, genetic_code_id)
        
        # Armazenar tudo
        all_data[nome_base] = {
            'counts': counts,
            'rscu': rscu,
            'enc': enc,
            'gc3': gc3,
            'cai': cai,
            'optimal': optimal,
            'rare': rare,
            'gc12': gc12
        }
        
        print(f"  📊 {nome_base}: ENC={enc:.2f}, GC3={gc3:.2f}%, GC12={gc12:.2f}%, CAI={cai:.3f}")

    status_queue.put(("progress", 90))
    
    if len(all_data) < 1:
        print("  ❌ Nenhum arquivo pôde ser processado para análise de viés.")
        return None

    print("\n✅ Processamento de viés concluído para todos os arquivos.")
    return all_data


# --- Funções Específicas das Análises (Frontend) ---

def gerar_rscu_heatmap_e_tabela(all_data, pasta_saida, genetic_code_id, status_queue):
    """(Análise 3) Análise de RSCU para UM arquivo, com heatmap detalhado."""
    
    nome_arquivo_base = list(all_data.keys())[0]
    data = all_data[nome_arquivo_base]
    codon_counts = data['counts']
    rscu_values = data['rscu']
    
    print(f"\n=== ANÁLISE RSCU HEATMAP ===")
    print(f"Arquivo: {nome_arquivo_base}")
    
    status_queue.put(("progress", 20))
    
    print("  Salvando tabela de contagem de códons...")
    df_counts = pd.DataFrame(codon_counts.items(), columns=['Cédon', 'Contagem']).sort_values(by='Contagem', ascending=False)
    csv_path = os.path.join(pasta_saida, f"contagem_codons_{nome_arquivo_base}.csv")
    df_counts.to_csv(csv_path, index=False, sep=';')
    print(f"  ✅ Tabela de contagem salva em: {csv_path}")
    
    status_queue.put(("progress", 50))
    
    print(f"  📊 Estatísticas:")
    print(f"     ENC (Effective Number of Codons): {data['enc']:.2f}")
    print(f"     GC3 (GC na terceira posição): {data['gc3']:.2f}%")
    print(f"     CAI (Codon Adaptation Index): {data['cai']:.3f}")
    
    status_queue.put(("progress", 70))
    print("  Preparando dados para o heatmap...")
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
    
    print("  Gerando Heatmap RSCU...")
    plt.figure(figsize=(20, 8))
    try:
        ax = sns.heatmap(df_rscu, annot=codon_labels_grid, fmt="", cmap="viridis", linewidths=0.5, cbar_kws={'label': 'Valor RSCU'})
        ax.set_yticklabels(['T', 'C', 'A', 'G'], rotation=0)
        ax.set_xticklabels([f"Pos {i+1}" for i in range(16)], rotation=0)
        ax.set_title(f"Viés de Uso de Códons (RSCU) - {nome_arquivo_base}\nENC={data['enc']:.2f}, GC3={data['gc3']:.2f}%, CAI={data['cai']:.3f}", fontsize=16)
        
        arquivo_saida = os.path.join(pasta_saida, f"rscu_heatmap_detalhado_{nome_arquivo_base}.png")
        plt.savefig(arquivo_saida, dpi=150, bbox_inches="tight")
        plt.close()
        
        print(f"\n✅ Heatmap RSCU salvo em: {arquivo_saida}")
        status_queue.put(("image_ready", (arquivo_saida, f"RSCU - {nome_arquivo_base}")))
        
    except Exception as e:
        print(f"\n❌ ERRO AO GERAR O HEATMAP RSCU (Análise 3): {e}")

def analise_rscu_comparativa(all_data, pasta_saida, status_queue):
    """(Análise 4) Análise comparativa (Clustermap & PCA) para MÚLTIPLOS arquivos."""
    print(f"\n=== ANÁLISE RSCU COMPARATIVA ===")
    
    # Extrair RSCU de all_data
    all_rscu_data = {species: data['rscu'] for species, data in all_data.items()}

    print("\n  Criando matriz RSCU comparativa...")
    status_queue.put(("progress", 55))
    df_rscu_matrix = pd.DataFrame.from_dict(all_rscu_data, orient='index', columns=ALL_CODONS_SORTED).fillna(0.0)
    csv_matrix_path = os.path.join(pasta_saida, 'rscu_matriz_comparativa.csv')
    df_rscu_matrix.to_csv(csv_matrix_path, sep=';', decimal='.')
    print(f"  ✅ Matriz RSCU salva em: {csv_matrix_path}")
    
    print("  Gerando Clustermap...")
    status_queue.put(("progress", 70))
    try:
        g = sns.clustermap(
            df_rscu_matrix, metric="euclidean", method="average", cmap="viridis", annot=False, linewidths=0.5,
            figsize=(max(15, len(ALL_CODONS_SORTED) * 0.2), max(8, len(all_data) * 0.5))
        )
        g.fig.suptitle("Análise Comparativa RSCU (Clustermap)", y=1.02, fontsize=16)
        plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90)
        
        clustermap_path = os.path.join(pasta_saida, 'rscu_comparativo_clustermap.png')
        g.savefig(clustermap_path, dpi=150, bbox_inches="tight")
        plt.close()
        
        print(f"  ✅ Clustermap salvo em: {clustermap_path}")
        status_queue.put(("image_ready", (clustermap_path, "RSCU Clustermap Comparativo")))
        
    except Exception as e:
        print(f"\n❌ ERRO AO GERAR O CLUSTERMAP: {e}")

    print("  Gerando Análise de Componentes Principais (PCA)...")
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
        plt.xlabel(f'Componente Principal 1 ({pc1_var:.2f}%)', fontsize=12)
        plt.ylabel(f'Componente Principal 2 ({pc2_var:.2f}%)', fontsize=12)
        plt.title('Análise de Componentes Principais (PCA) do RSCU', fontsize=16)
        plt.axhline(0, color='grey', linestyle='--', linewidth=0.5)
        plt.axvline(0, color='grey', linestyle='--', linewidth=0.5)
        
        pca_path = os.path.join(pasta_saida, 'rscu_comparativo_pca.png')
        plt.savefig(pca_path, dpi=150, bbox_inches="tight")
        plt.close()
        
        print(f"  ✅ Gráfico PCA salvo em: {pca_path}")
        print(f"  📈 Variância explicada: PC1={pc1_var:.2f}%, PC2={pc2_var:.2f}%")
        status_queue.put(("image_ready", (pca_path, "RSCU Análise PCA")))
        
    except Exception as e:
        print(f"\n❌ ERRO AO GERAR O PCA: {e}")

def analise_correlacao_rscu(all_data, pasta_saida, status_queue):
    """(Análise 5) Gera um gráfico de correlação entre DOIS arquivos."""
    print(f"\n=== ANÁLISE DE CORRELAÇÃO RSCU ===")

    # Extrair RSCU de all_data
    all_rscu_data = {species: data['rscu'] for species, data in all_data.items()}
    
    status_queue.put(("progress", 60))
    df_rscu_matrix = pd.DataFrame.from_dict(all_rscu_data, orient='index', columns=ALL_CODONS_SORTED).fillna(0.0)
    
    species_x, species_y = df_rscu_matrix.index[0], df_rscu_matrix.index[1]
    rscu_x, rscu_y = df_rscu_matrix.iloc[0], df_rscu_matrix.iloc[1]

    print(f"  Calculando correlação entre '{species_x}' e '{species_y}'...")
    r, p_value = pearsonr(rscu_x, rscu_y)
    
    print(f"    📊 Coeficiente de Pearson (R): {r:.4f}")
    print(f"    📊 P-valor: {p_value:.2e}")

    status_queue.put(("progress", 80))
    print("  Gerando gráfico de correlação...")
    
    plt.figure(figsize=(10, 8))
    ax = sns.regplot(x=rscu_x, y=rscu_y,
                     scatter_kws={'alpha': 0.5},
                     line_kws={'color': 'darkred', 'linewidth': 2})
                     
    ax.set_xlabel(f"{species_x}", fontsize=12)
    ax.set_ylabel(f"{species_y}", fontsize=12)
    
    title_str = f"Correlação RSCU: {species_x} vs {species_y}\n(R={r:.3f}, p-value={p_value:.2e})"
    ax.set_title(title_str, fontsize=16)

    correlation_path = os.path.join(pasta_saida, f'rscu_correlacao_{species_x}_vs_{species_y}.png')
    plt.savefig(correlation_path, dpi=150, bbox_inches="tight")
    plt.close()
    
    print(f"  ✅ Gráfico de correlação salvo em: {correlation_path}")
    status_queue.put(("image_ready", (correlation_path, f"Correlação RSCU")))

def gerar_histogramas_rscu(all_data, pasta_saida, genetic_code_id, status_queue):
    """(Análise 6) Gera histogramas comparativos de RSCU melhorados."""
    print(f"\n=== ANÁLISE DE HISTOGRAMAS RSCU ===")
    
    # Extrair RSCU de all_data
    all_rscu_data = {species: data['rscu'] for species, data in all_data.items()}
    
    status_queue.put(("progress", 50))
    df_rscu_matrix = pd.DataFrame.from_dict(all_rscu_data, orient='index', columns=ALL_CODONS_SORTED).fillna(0.0)
    
    # --- Preparação dos dados para formato "longo" ---
    df_long = df_rscu_matrix.reset_index().rename(columns={'index': 'Species'}).melt(id_vars='Species', var_name='Codon', value_name='RSCU')
    codon_aa_map = GENETIC_CODE_TABLES.get(genetic_code_id, GENETIC_CODE_TABLES[1])
    df_long['AminoAcid'] = df_long['Codon'].map(codon_aa_map)
    df_long = df_long.sort_values(by=['AminoAcid', 'Codon'])

    print(f"  📊 Processando {len(all_data)} espécies para análise de histogramas")

    # --- Gráfico 1: Box Plot por Aminoácido (Melhor que barras horizontais) ---
    print("  Gerando Gráfico 1: Box Plot de RSCU por Aminoácido...")
    status_queue.put(("progress", 60))
    try:
        plt.figure(figsize=(20, 10))
        
        # Filtrar AAs sinônimos
        df_filtered = df_long[df_long['AminoAcid'].isin(AA_CODON_MAPS[genetic_code_id].keys())]
        df_filtered = df_filtered[df_filtered['AminoAcid'] != '*']
        
        aa_order = sorted(df_filtered['AminoAcid'].unique())
        
        sns.boxplot(x='AminoAcid', y='RSCU', data=df_filtered, order=aa_order, palette="Set3")
        sns.stripplot(x='AminoAcid', y='RSCU', data=df_filtered, order=aa_order, color=".25", size=3, alpha=0.5)
        
        plt.xlabel('Aminoácidos', fontsize=12)
        plt.ylabel('Valor RSCU', fontsize=12)
        plt.title('Distribuição de RSCU por Aminoácido (Todas as Espécies)', fontsize=16, fontweight='bold')
        plt.grid(axis='y', alpha=0.3)
        plt.axhline(y=1.0, color='red', linestyle='--', alpha=0.5, label='RSCU = 1.0')
        plt.legend()
        
        plt.tight_layout()
        
        box_path = os.path.join(pasta_saida, 'rscu_boxplot_por_aminoacido.png')
        plt.savefig(box_path, dpi=150, bbox_inches="tight")
        plt.close()
        print(f"  ✅ Box plot salvo em: {box_path}")
        status_queue.put(("image_ready", (box_path, "RSCU Box Plot por AA")))
    except Exception as e:
        print(f"\n❌ ERRO AO GERAR BOX PLOT: {e}")

    # --- Gráfico 2: Barras Verticais Empilhadas (por Espécie) ---
    print("\n  Gerando Gráfico 2: Barras Verticais Empilhadas por Espécie...")
    status_queue.put(("progress", 85))
    aa_codon_map = AA_CODON_MAPS.get(genetic_code_id, AA_CODON_MAPS[1])
    aa_list = sorted([aa for aa, codons in aa_codon_map.items() if aa != '*' and len(codons) > 1]) # Apenas sinônimos
    
    for species_name, rscu_series in df_rscu_matrix.iterrows():
        try:
            fig, ax = plt.subplots(figsize=(16, 10))
            
            # Preparar dados para plot
            plot_data = []
            for aa in aa_list:
                syn_codons = sorted(aa_codon_map[aa])
                for codon in syn_codons:
                    value = rscu_series.get(codon, 0)
                    plot_data.append({'AminoAcid': aa, 'Codon': codon, 'RSCU': value})
            
            df_plot = pd.DataFrame(plot_data)

            # Criar o gráfico de barras empilhadas
            df_plot.pivot(index='AminoAcid', columns='Codon', values='RSCU').loc[aa_list].plot(
                kind='bar', stacked=True, ax=ax,
                colormap='tab20',
                edgecolor='black', linewidth=0.5
            )

            ax.set_xlabel("Aminoácido", fontsize=12, fontweight='bold')
            ax.set_ylabel("Valor RSCU (Empilhado)", fontsize=12)
            ax.set_title(f"Perfil de Uso de Códons (RSCU) - {species_name}", fontsize=16, fontweight='bold')
            ax.legend(title='Códons',loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=8, fontsize='small')
            ax.grid(axis='y', linestyle='--', alpha=0.7)
            plt.xticks(rotation=0)
            
            plt.tight_layout() # Ajustar para a legenda
            
            hist_v_path = os.path.join(pasta_saida, f'rscu_histograma_vertical_{species_name}.png')
            plt.savefig(hist_v_path, dpi=150, bbox_inches="tight")
            plt.close()
            print(f"  ✅ Gráfico de barras verticais salvo para {species_name}")
            status_queue.put(("image_ready", (hist_v_path, f"RSCU Vertical - {species_name}")))
        except Exception as e:
            print(f"\n❌ ERRO AO GERAR HISTOGRAMA VERTICAL para {species_name}: {e}")

def analise_enc_gc3(all_data, pasta_saida, status_queue):
    """(Análise 7) Análise de ENC vs GC3"""
    print(f"\n=== ANÁLISE ENC vs GC3 ===")
    
    status_queue.put(("progress", 60))
    
    # Extrair dados de ENC e GC3
    plot_data = []
    for species, data in all_data.items():
        plot_data.append({'species': species, 'ENC': data['enc'], 'GC3': data['gc3']})
    df_plot = pd.DataFrame(plot_data)

    # Gerar gráfico
    plt.figure(figsize=(12, 8))
    
    # Curva teórica esperada (Wright, 1990)
    # Correto (Novembre, 2002):
    s_values = np.linspace(0.01, 0.99, 200)
    enc_expected = []
    for s in s_values:
        f_s = s**2 + (1-s)**2
        enc_val = 2 + (9/f_s) + (1/f_s) + (5/f_s) + (3/f_s)
        enc_expected.append(enc_val)
        
    # Plotar a curva com s * 100 (para %)
    plt.plot(s_values * 100, enc_expected, 'r--', label='Curva Esperada (Viés Mutacional)', linewidth=2)
    
    # Pontos de dados
    scatter = plt.scatter(df_plot['GC3'], df_plot['ENC'], c=df_plot.index, 
                         cmap='viridis', s=100, alpha=0.7, edgecolors='black')
    
    # Adicionar labels
    for i, row in df_plot.iterrows():
        plt.annotate(row['species'], (row['GC3'], row['ENC']), 
                    xytext=(5, 5), textcoords='offset points', fontsize=9)
    
    plt.xlabel('GC3 (%)', fontsize=12)
    plt.ylabel('ENC (Effective Number of Codons)', fontsize=12)
    plt.title('Análise ENC vs GC3 (Gráfico de Wright)', fontsize=16)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.xlim(0, 100)
    plt.ylim(20, 61)
    
    plt.tight_layout()
    
    output_path = os.path.join(pasta_saida, 'enc_gc3_analysis.png')
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"\n✅ Análise ENC vs GC3 salva em: {output_path}")
    status_queue.put(("image_ready", (output_path, "ENC vs GC3 Analysis")))
    
    # Salvar dados
    df_plot.to_csv(os.path.join(pasta_saida, 'enc_gc3_results.csv'), sep=';', index=False)

def analise_composicao_genomica(lista_arquivos, pasta_saida, status_queue):
    """(Análise 8) Análise da composição genômica"""
    print(f"\n=== ANÁLISE DE COMPOSIÇÃO GENÔMICA ===")
    
    all_data = {}
    
    for i, caminho_completo in enumerate(lista_arquivos):
        nome_base = os.path.basename(caminho_completo).split('.')[0]
        print(f"\nProcessando arquivo {i+1}/{len(lista_arquivos)}: {nome_base}")
        status_queue.put(("progress", int(20 + (i / len(lista_arquivos)) * 70)))
        
        comp_data = analisar_composicao_nucleotideos(caminho_arquivo)
        if comp_data:
            all_data[nome_base] = comp_data
            print(f"  📊 {nome_base}: GC={comp_data['gc_content']:.2f}%, Tamanho={comp_data['total_length']} pb")
    
    if not all_data:
        print("  ❌ Erro: Nenhum dado processado.")
        return
    
    status_queue.put(("progress", 90))
    
    # Gerar gráficos
    fig, ((ax1, ax2)) = plt.subplots(1, 2, figsize=(16, 8))
    
    species = list(all_data.keys())
    
    # Gráfico 1: Composição de nucleotídeos (Barras agrupadas horizontais)
    nt_compositions = [all_data[s]['composition'] for s in species]
    df_nt = pd.DataFrame(nt_compositions, index=species)
    
    df_nt.plot(kind='barh', stacked=False, ax=ax1)
    
    ax1.set_title('Composição de Nucleotídeos (%)', fontweight='bold')
    ax1.set_xlabel('Percentagem (%)') 
    ax1.legend(title='Nucleotídeos', bbox_to_anchor=(1.0, 1), loc='upper left')
    ax1.grid(axis='x', alpha=0.3)
    
    # Gráfico 2: Tamanho do genoma vs Conteúdo de GC (Scatter plot)
    gc_contents = [all_data[s]['gc_content'] for s in species]
    genome_sizes = [all_data[s]['total_length'] for s in species]

    scatter = ax2.scatter(gc_contents, genome_sizes, c=range(len(species)), 
                          cmap='viridis', s=100, alpha=0.7, edgecolors='black')
                          
    ax2.set_title('Tamanho do Genoma vs Conteúdo de GC', fontweight='bold')
    ax2.set_xlabel('GC (%)')
    ax2.set_ylabel('Tamanho (pb)')
    ax2.grid(True, alpha=0.3)
    
    # Adicionar labels
    for i, species_name in enumerate(species):
        ax2.annotate(species_name, (gc_contents[i], genome_sizes[i]), 
                    xytext=(5, 5), textcoords='offset points', fontsize=9)
    
    ax2.ticklabel_format(style='scientific', axis='y', scilimits=(0,0))
    
    plt.suptitle('Análise de Composição Genômica', fontsize=16, fontweight='bold')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    
    output_path = os.path.join(pasta_saida, 'composicao_genomica.png')
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"\n✅ Análise de composição genômica salva em: {output_path}")
    status_queue.put(("image_ready", (output_path, "Composição Genômica")))
    
    # Salvar dados
    df_results = pd.DataFrame({
        'Especie': species,
        'Tamanho_pb': genome_sizes,
        'GC_percent': gc_contents,
        'A_percent': [all_data[s]['composition']['A'] for s in species],
        'T_percent': [all_data[s]['composition']['T'] for s in species],
        'G_percent': [all_data[s]['composition']['G'] for s in species],
        'C_percent': [all_data[s]['composition']['C'] for s in species],
    })
    df_results.to_csv(os.path.join(pasta_saida, 'composicao_genomica_results.csv'), sep=';', index=False)

def analise_codons_otimos_raros(all_data, pasta_saida, status_queue):
    """(Análise 9) Análise de códons ótimos e raros"""
    print(f"\n=== ANÁLISE DE CÓDONS ÓTIMOS E RAROS ===")
    
    status_queue.put(("progress", 60))
    
    # Gerar visualização
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    
    species = list(all_data.keys())
    
    # Gráfico 1: Contagem de códons ótimos e raros
    optimal_counts = [len(all_data[s]['optimal']) for s in species]
    rare_counts = [len(all_data[s]['rare']) for s in species]
    
    x = np.arange(len(species))
    width = 0.35
    
    ax1.bar(x - width/2, optimal_counts, width, label='Ótimos (>1.2)', color='green', alpha=0.7)
    ax1.bar(x + width/2, rare_counts, width, label='Raros (<0.8)', color='red', alpha=0.7)
    
    ax1.set_title('Códons Ótimos vs Raros por Espécie', fontweight='bold')
    ax1.set_ylabel('Número de Códons')
    ax1.set_xticks(x)
    ax1.set_xticklabels(species, rotation=45, ha="right")
    ax1.legend()
    ax1.grid(axis='y', alpha=0.3)
    
    # Gráfico 2: CAI values
    cai_values = [all_data[s]['cai'] for s in species]
    bars = ax2.bar(species, cai_values, color='purple', alpha=0.7)
    ax2.set_title('Codon Adaptation Index (CAI)', fontweight='bold')
    ax2.set_ylabel('CAI')
    ax2.grid(axis='y', alpha=0.3)
    plt.setp(ax2.get_xticklabels(), rotation=45, ha="right")
    
    # Adicionar valores nas barras
    for bar, value in zip(bars, cai_values):
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01, 
                f'{value:.3f}', ha='center', va='bottom')
    
    plt.suptitle('Análise de Códons Ótimos, Raros e CAI', fontsize=16, fontweight='bold')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    
    output_path = os.path.join(pasta_saida, 'codons_otimos_raros.png')
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"\n✅ Análise de códons ótimos e raros salva em: {output_path}")
    status_queue.put(("image_ready", (output_path, "Códons Ótimos e Raros")))
    
    # Salvar dados detalhados
    for species_name, data in all_data.items():
        df_optimal = pd.DataFrame(list(data['optimal'].items()), 
                                 columns=['Aminoácido', 'Códon Ótimo'])
        df_rare = pd.DataFrame(list(data['rare'].items()), 
                              columns=['Aminoácido', 'Códon Raro'])
        
        df_optimal.to_csv(os.path.join(pasta_saida, f'{species_name}_codons_otimos.csv'), 
                         sep=';', index=False)
        df_rare.to_csv(os.path.join(pasta_saida, f'{species_name}_codons_raros.csv'), 
                      sep=';', index=False)

# ######################################################################
# --- NOVAS FUNÇÕES DE ANÁLISE (10, 11, 12, 13) ---
# ######################################################################

# --- Funções Auxiliares Comuns ---

def extrair_sequencias_cds(lista_arquivos, status_queue, gene_list=None):
    """Extrai todas as sequências de CDS de uma lista de arquivos .gbk."""
    all_seqs_by_species = defaultdict(list)
    total_files = len(lista_arquivos)
    
    print("  Extraindo sequências de CDS...")
    
    for i, caminho_completo in enumerate(lista_arquivos):
        nome_base = os.path.basename(caminho_completo).split('.')[0]
        print(f"  Lendo CDS de {nome_base}...")
        status_queue.put(("message", f"Lendo CDS de {nome_base}..."))
        
        try:
            for record in SeqIO.parse(caminho_completo, "genbank"):
                for feature in record.features:
                    if feature.type != "CDS": continue
                    
                    # --- APLICAR FILTRO ---
                    if not _apply_gene_filter(feature, gene_list):
                        continue
                    # --- FIM DO FILTRO ---
                    
                    seq_cds = feature.extract(record.seq).upper()
                    
                    if feature.location.strand == -1:
                        seq_cds = seq_cds.reverse_complement()
                        
                    codon_start = int(feature.qualifiers.get("codon_start", ["1"])[0])
                    seq_cds = seq_cds[codon_start - 1:]
                    
                    tamanho_real = (len(seq_cds) // 3) * 3
                    seq_cds = seq_cds[:tamanho_real]
                    
                    if 'N' not in seq_cds and len(seq_cds) >= 6: # Precisa de pelo menos 2 códons
                        all_seqs_by_species[nome_base].append(str(seq_cds))
                        
        except Exception as e:
            print(f"  ❌ Erro ao extrair CDS do arquivo {nome_base}: {e}")
            
        status_queue.put(("progress", int(10 + (i / total_files) * 20)))
        
    print(f"  ✅ Extração de CDS concluída. {len(all_seqs_by_species)} espécies processadas.")
    return all_seqs_by_species

# --- INÍCIO: NOVAS FUNÇÕES AUXILIARES (BIOLOGIA SINTÉTICA) ---

def _contar_codons_de_sequencia(seq_str, genetic_code_id=1):
    """Conta códons de uma única string de DNA."""
    genetic_code = GENETIC_CODE_TABLES.get(genetic_code_id, GENETIC_CODE_TABLES[1])
    codon_counts = Counter()
    
    seq_str = seq_str.upper()
    
    for i in range(0, len(seq_str) - 2, 3):
        codon = seq_str[i:i+3]
        if len(codon) == 3 and codon in genetic_code:
            codon_counts[codon] += 1
            
    return codon_counts

def _get_optimal_codon_map(rscu_values, genetic_code_id=1):
    """
    Cria um mapa de AA -> Códon Ótimo (baseado no maior RSCU).
    """
    aa_codon_map = AA_CODON_MAPS.get(genetic_code_id, AA_CODON_MAPS[1])
    optimal_map = {}
    
    for aa, codons in aa_codon_map.items():
        if not codons: continue
        
        # Encontra o códon com o maior RSCU para este AA
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
    """
    Para cada AA, retorna uma lista ordenada de seus códons sinônimos,
    do mais frequente para o menos frequente.
    """
    aa_codon_map = AA_CODON_MAPS.get(genetic_code_id, AA_CODON_MAPS[1])
    rank_map = {}
    
    for aa, codons in aa_codon_map.items():
        if not codons or len(codons) == 1:
            rank_map[aa] = codons # Apenas um, sem ranking
            continue
            
        # Criar uma lista de (códon, contagem)
        codon_pairs = [(codon, codon_counts.get(codon, 0)) for codon in codons]
        
        # Ordenar pela contagem (decrescente)
        codon_pairs.sort(key=lambda x: x[1], reverse=True)
        
        # Armazenar apenas a lista ordenada de códons
        rank_map[aa] = [codon for codon, count in codon_pairs]
        
    return rank_map

def _calcular_mfe_para_sequencias(sequences, mfe_region_length):
    """
    Calcula o MFE para a região 5' de uma lista de sequências.
    Requer VIENNARNA_DISPONIVEL = True.
    """
    if not VIENNARNA_DISPONIVEL:
        raise ImportError("ViennaRNA (import RNA) não encontrado. Não é possível calcular MFE.")
        
    mfe_results = []
    
    for seq in sequences:
        if len(seq) < mfe_region_length:
            continue
            
        # Extrai a região 5' (início)
        region_5prime = seq[:mfe_region_length]
        
        # Calcula MFE
        try:
            # RNA.fold() retorna (estrutura_em_pontos, mfe_kcal_mol)
            (structure, mfe) = RNA.fold(region_5prime)
            mfe_results.append({'mfe': mfe, 'length': len(region_5prime)})
        except Exception as e:
            print(f"  Aviso: Falha ao calcular MFE para sequência {region_5prime[:10]}...: {e}")
            
    return pd.DataFrame(mfe_results)

# --- FIM: NOVAS FUNÇÕES AUXILIARES ---

# --- NOVAS FUNÇÕES AUXILIARES (Para Análises 18 e 19) ---

def get_w_reference_table(rscu_values, genetic_code_id=1):
    """
    Calcula a tabela de pesos (w) de adaptabilidade relativa, usada para o CAI.
    (Extraído da função calcular_cai)
    """
    aa_codon_map = AA_CODON_MAPS.get(genetic_code_id, AA_CODON_MAPS[1])
    w_values = {}
    
    for aa, codons in aa_codon_map.items():
        if aa == '*' or len(codons) == 1:
            # Atribuir peso 1 a códons únicos (Met, Trp)
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
                w_values[codon] = 1.0 # Se não houver dados, peso é 1
    
    return w_values

def calcular_metricas_por_gene(caminho_arquivo, w_reference_table, genetic_code_id, gene_list):
    """
    Calcula métricas (ENC, GC3, CAI) para genes INDIVIDUAIS em um arquivo.
    Usa um 'w_reference_table' pré-calculado (do genoma inteiro) para o CAI.
    Retorna uma lista de dicionários, um por gene.
    """
    print(f"  Calculando métricas por gene para {len(gene_list)} genes...")
    results = []
    
    try:
        for record in SeqIO.parse(caminho_arquivo, "genbank"):
            for feature in record.features:
                if feature.type != "CDS": continue
                
                # --- APLICAR FILTRO DE GENE ---
                if not _apply_gene_filter(feature, gene_list):
                    continue
                
                # Obter identificador
                identificador = "N/A"
                if "locus_tag" in feature.qualifiers: identificador = feature.qualifiers["locus_tag"][0]
                elif "gene" in feature.qualifiers: identificador = feature.qualifiers["gene"][0]
                elif "protein_id" in feature.qualifiers: identificador = feature.qualifiers["protein_id"][0]
                
                # Extrair sequência
                seq_cds_str = ""
                try:
                    seq_cds = feature.extract(record.seq).upper()
                    if feature.location.strand == -1:
                        seq_cds = seq_cds.reverse_complement()
                    codon_start = int(feature.qualifiers.get("codon_start", ["1"])[0])
                    seq_cds_str = str(seq_cds[codon_start - 1:])
                except Exception as e_seq:
                    print(f"    Aviso: Pulando gene {identificador} (erro de extração): {e_seq}")
                    continue
                
                # Contar códons para ESTE gene
                counts_gene = _contar_codons_de_sequencia(seq_cds_str, genetic_code_id)
                if sum(counts_gene.values()) == 0:
                    continue
                
                # 1. Calcular ENC do gene
                enc_gene = calcular_enc(counts_gene, genetic_code_id)
                
                # 2. Calcular GC3 do gene
                gc3_gene = calcular_gc3(counts_gene, genetic_code_id)
                
                # 3. Calcular CAI do gene (usando tabela de referência W)
                log_w_sum = 0
                total_codons = 0
                for codon, count in counts_gene.items():
                    if codon in w_reference_table:
                        if w_reference_table[codon] > 0:
                            log_w_sum += count * np.log(w_reference_table[codon])
                        total_codons += count
                
                cai_gene = np.exp(log_w_sum / total_codons) if total_codons > 0 else 0
                
                results.append({
                    'gene': identificador,
                    'enc': enc_gene,
                    'gc3': gc3_gene,
                    'cai': cai_gene,
                    'counts': counts_gene # Retorna contagens para Qui-Quadrado
                })
                
    except Exception as e:
        print(f"  ❌ Erro ao calcular métricas por gene: {e}")
        
    return results

# --- FIM DAS NOVAS FUNÇÕES AUXILIARES ---


# --- Análise 10: Viés de Pares de Códons (CPB) ---

def calcular_codon_pair_bias(sequences, genetic_code_id):
    """Calcula a matriz de Codon Pair Bias (CPS) para um conjunto de sequências."""
    
    genetic_code = GENETIC_CODE_TABLES.get(genetic_code_id, GENETIC_CODE_TABLES[1])
    stop_codons = {codon for codon, aa in genetic_code.items() if aa == '*'}
    
    # Usar ALL_CODONS_SORTED para garantir a ordem, mas excluir stops
    non_stop_codons = sorted([c for c in ALL_CODONS_SORTED if c not in stop_codons])
    codon_index = {codon: i for i, codon in enumerate(non_stop_codons)}
    
    num_codons = len(non_stop_codons)
    pair_counts_obs = np.zeros((num_codons, num_codons))
    codon_counts_c1 = np.zeros(num_codons) # Contagem de códons na Posição 1
    codon_counts_c2 = np.zeros(num_codons) # Contagem de códons na Posição 2
    
    total_pairs = 0
    
    for seq in sequences:
        for i in range(0, len(seq) - 5, 3): # Itera em pares de códons
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
        print("  Aviso: Nenhum par de códons válido encontrado.")
        return pd.DataFrame()
        
    # Calcular frequências esperadas
    # E(ij) = (Obs(Ci) * Obs(Cj)) / N_total
    pair_counts_exp = np.outer(codon_counts_c1, codon_counts_c2) / total_pairs
    
    # Calcular CPS = log(Obs / Exp)
    # Adicionar um pequeno epsilon para evitar divisão por zero
    epsilon = 1e-9
    cps_matrix = np.log((pair_counts_obs + epsilon) / (pair_counts_exp + epsilon))
    
    df_cps = pd.DataFrame(cps_matrix, index=non_stop_codons, columns=non_stop_codons)
    return df_cps

def analise_codon_pair_bias(lista_arquivos, pasta_saida, genetic_code_id, status_queue, gene_list=None):
    """(Análise 10) Gera heatmap e CSV para Codon Pair Bias."""
    print(f"\n=== ANÁLISE DE VIÉS DE PARES DE CÓDONS (CPB) ===")
    
    all_seqs_by_species = extrair_sequencias_cds(lista_arquivos, status_queue, gene_list)
    
    if not all_seqs_by_species:
        print("  ❌ Erro: Nenhuma sequência de CDS válida foi extraída.")
        return

    status_queue.put(("progress", 40))
    
    for i, (species_name, sequences) in enumerate(all_seqs_by_species.items()):
        print(f"  Calculando CPB para {species_name}...")
        status_queue.put(("message", f"Calculando CPB para {species_name}..."))
        
        df_cps = calcular_codon_pair_bias(sequences, genetic_code_id)
        
        if df_cps.empty:
            print(f"  Aviso: CPB vazio para {species_name}.")
            continue
            
        # Salvar CSV
        csv_path = os.path.join(pasta_saida, f"cpb_matrix_{species_name}.csv")
        df_cps.to_csv(csv_path, sep=';', decimal='.')
        print(f"  ✅ Matriz CPB salva em: {csv_path}")
        
        # Gerar Heatmap
        print("  Gerando Heatmap CPB...")
        plt.figure(figsize=(24, 20))
        sns.heatmap(df_cps, cmap="coolwarm", center=0, annot=False, 
                    cbar_kws={'label': 'Codon Pair Score (log(Obs/Exp))'})
        plt.title(f"Viés de Pares de Códons (CPB) - {species_name}", fontsize=18)
        plt.xlabel("Segundo Códon", fontsize=12)
        plt.ylabel("Primeiro Códon", fontsize=12)
        plt.tight_layout()
        
        plot_path = os.path.join(pasta_saida, f"cpb_heatmap_{species_name}.png")
        plt.savefig(plot_path, dpi=100, bbox_inches="tight") # DPI 100 para não ficar gigante
        plt.close()
        
        print(f"  ✅ Heatmap CPB salvo em: {plot_path}")
        status_queue.put(("image_ready", (plot_path, f"CPB - {species_name}")))
        
        status_queue.put(("progress", int(40 + (i / len(all_seqs_by_species)) * 50)))

# --- Análise 11: Análise Físico-Química (GRAVY & Aromo) ---

def calcular_gravy_aromo(sequences, genetic_code_id):
    """Calcula GRAVY e Aromo para uma lista de sequências de CDS."""
    
    codon_map = GENETIC_CODE_TABLES.get(genetic_code_id, GENETIC_CODE_TABLES[1])
    results = []
    
    for seq in sequences:
        aa_seq = []
        for i in range(0, len(seq), 3):
            codon = seq[i:i+3]
            aa = codon_map.get(codon, '*')
            if aa == '*':
                break # Parar tradução no stop codon
            aa_seq.append(aa)
            
        if not aa_seq:
            continue
        
        total_len = len(aa_seq)
        
        # Calcular GRAVY
        try:
            gravy_score = sum(KYTE_DOOLITTLE_HYDROPATHY.get(aa, 0) for aa in aa_seq) / total_len
        except ZeroDivisionError:
            gravy_score = 0
            
        # Calcular Aromaticidade
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

def analise_gravy_aromo(lista_arquivos, pasta_saida, genetic_code_id, status_queue, gene_list=None):
    """(Análise 11) Gera boxplots e CSV para GRAVY e Aromo."""
    print(f"\n=== ANÁLISE FÍSICO-QUÍMICA (GRAVY & AROMO) ===")
    
    all_seqs_by_species = extrair_sequencias_cds(lista_arquivos, status_queue, gene_list)
    
    if not all_seqs_by_species:
        print("  ❌ Erro: Nenhuma sequência de CDS válida foi extraída.")
        return

    status_queue.put(("progress", 40))
    all_results_long = []
    
    for species_name, sequences in all_seqs_by_species.items():
        print(f"  Calculando GRAVY/Aromo para {species_name}...")
        status_queue.put(("message", f"Calculando GRAVY/Aromo para {species_name}..."))
        
        df_results = calcular_gravy_aromo(sequences, genetic_code_id)
        
        if df_results.empty:
            continue
            
        # Salvar CSV por gene
        csv_path = os.path.join(pasta_saida, f"gravy_aromo_per_gene_{species_name}.csv")
        df_results.to_csv(csv_path, sep=';', decimal='.', index=False)
        print(f"  ✅ Tabela GRAVY/Aromo (por gene) salva em: {csv_path}")
        
        # Preparar dados para plot agregado
        df_results['species'] = species_name
        all_results_long.append(df_results)

    if not all_results_long:
        print("  ❌ Erro: Nenhum dado de GRAVY/Aromo pôde ser calculado.")
        return
        
    df_plot_long = pd.concat(all_results_long)
    
    # Gerar Plots
    print("  Gerando Box Plots comparativos...")
    status_queue.put(("progress", 80))
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    
    # --- INÍCIO DA MODIFICAÇÃO (Violin Plot) ---
    
    # Boxplot para GRAVY
    sns.violinplot(x='species', y='gravy', data=df_plot_long, ax=ax1, palette="coolwarm",
                   inner="quartile") # inner="quartile" mostra os quartis (como o boxplot)
    ax1.set_title("Distribuição GRAVY (Hidropaticidade)", fontweight='bold')
    ax1.set_ylabel("Escore GRAVY")
    ax1.set_xlabel("Espécie")
    ax1.tick_params(axis='x', rotation=45)
    
    # Boxplot para Aromo
    sns.violinplot(x='species', y='aromo', data=df_plot_long, ax=ax2, palette="viridis",
                   inner="quartile") # inner="quartile" mostra os quartis (como o boxplot)
    ax2.set_title("Distribuição de Aromaticidade", fontweight='bold')
    ax2.set_ylabel("Escore Aromo (% F, Y, W)")
    ax2.set_xlabel("Espécie")
    ax2.tick_params(axis='x', rotation=45)
    
    # --- FIM DA MODIFICAÇÃO ---
    
    plt.suptitle("Análise Físico-Química Comparativa", fontsize=16, fontweight='bold')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    
    plot_path = os.path.join(pasta_saida, "gravy_aromo_comparativo.png")
    plt.savefig(plot_path, dpi=150, bbox_inches="tight")
    plt.close()
    
    print(f"  ✅ Gráfico comparativo salvo em: {plot_path}")
    status_queue.put(("image_ready", (plot_path, "Análise GRAVY & Aromo"))) # <-- LINHA CORRIGIDA

# --- Análise 12: Gráfico de Neutralidade (GC12 vs GC3) ---

def analise_neutrality_plot(all_data, pasta_saida, status_queue):
    """(Análise 12) Gera o Gráfico de Neutralidade (GC12 vs GC3)."""
    print(f"\n=== ANÁLISE GRÁFICO DE NEUTRALIDADE ===")
    
    status_queue.put(("progress", 60))
    
    # Extrair dados de GC12 e GC3 (calculados em processar_genomas_para_analise_vies)
    plot_data = []
    for species, data in all_data.items():
        if 'gc12' in data and 'gc3' in data:
            plot_data.append({
                'species': species, 
                'GC12': data['gc12'], 
                'GC3': data['gc3']
            })
    
    if not plot_data:
        print("  ❌ Erro: Dados de GC12 ou GC3 não encontrados. (A Análise 12 depende do processamento de viés)")
        return
        
    df_plot = pd.DataFrame(plot_data)

    # Calcular Regressão Linear
    # Remover NaNs/Infs se houver
    df_plot = df_plot.replace([np.inf, -np.inf], np.nan).dropna()
    
    if df_plot.empty or len(df_plot) < 2:
        print("  ❌ Erro: Dados insuficientes para regressão linear.")
        return

    slope, intercept, r_value, p_value, std_err = linregress(df_plot['GC3'], df_plot['GC12'])
    
    print(f"  📊 Regressão: Inclinação (Slope) = {slope:.3f}")
    print(f"  📊 Regressão: R² = {r_value**2:.3f}")

    # Gerar gráfico
    plt.figure(figsize=(12, 8))
    
    # Pontos de dados
    sns.scatterplot(x='GC3', y='GC12', data=df_plot, s=100, alpha=0.7, hue='species', legend=False)
    
    # Linha de regressão
    x_vals = np.array(plt.xlim())
    y_vals = intercept + slope * x_vals
    plt.plot(x_vals, y_vals, 'r--', label=f'Regressão (Inclinação = {slope:.3f})')
    
    # Adicionar labels
    for i, row in df_plot.iterrows():
        plt.annotate(row['species'], (row['GC3'], row['GC12']), 
                    xytext=(5, 5), textcoords='offset points', fontsize=9)
    
    plt.xlabel('GC3 (%)', fontsize=12)
    plt.ylabel('GC12 (%)', fontsize=12)
    plt.title(f"Gráfico de Neutralidade (GC12 vs GC3)\nSlope = {slope:.3f}, R² = {r_value**2:.3f}", fontsize=16)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.xlim(0, 100)
    plt.ylim(0, 100)
    
    plt.tight_layout()
    
    output_path = os.path.join(pasta_saida, 'neutrality_plot_gc12_vs_gc3.png')
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"\n✅ Gráfico de Neutralidade salvo em: {output_path}")
    status_queue.put(("image_ready", (output_path, "Gráfico de Neutralidade")))
    
    # Salvar dados
    df_plot.to_csv(os.path.join(pasta_saida, 'neutrality_plot_results.csv'), sep=';', index=False)

# --- Análise 13: Composição de Dinucleotídeos ---

def analise_dinucleotide_composition(lista_arquivos, pasta_saida, status_queue):
    """(Análise 13) Calcula a frequência de dinucleotídeos no genoma completo."""
    print(f"\n=== ANÁLISE DE COMPOSIÇÃO DE DINUCLEOTÍDEOS ===")
    
    all_results = {}
    total_files = len(lista_arquivos)
    
    # Criar todos os 16 pares
    dinu_order = [''.join(p) for p in product('ATGC', repeat=2)]
    
    for i, caminho_completo in enumerate(lista_arquivos):
        nome_base = os.path.basename(caminho_completo).split('.')[0]
        print(f"  Analisando dinucleotídeos de {nome_base}...")
        status_queue.put(("message", f"Analisando dinucleotídeos de {nome_base}..."))
        
        counts = Counter()
        total_pairs = 0
        
        try:
            for record in SeqIO.parse(caminho_completo, "genbank"):
                seq = record.seq.upper()
                for j in range(len(seq) - 1):
                    dinu = seq[j:j+2]
                    if dinu in dinu_order: # Apenas conta pares ATGC válidos
                        counts[dinu] += 1
                        total_pairs += 1
                        
            if total_pairs > 0:
                freqs = {dinu: (counts.get(dinu, 0) / total_pairs) * 100 for dinu in dinu_order}
                all_results[nome_base] = freqs
            else:
                print(f"  Aviso: Nenhum par de dinucleotídeos encontrado em {nome_base}")
                
        except Exception as e:
            print(f"  ❌ Erro ao analisar dinucleotídeos em {nome_base}: {e}")
            
        status_queue.put(("progress", int(20 + (i / total_files) * 70)))
        
    if not all_results:
        print("  ❌ Erro: Nenhum dado de dinucleotídeos processado.")
        return

    # Criar DataFrame e salvar CSV
    df_dinu = pd.DataFrame.from_dict(all_results, orient='index', columns=dinu_order)
    csv_path = os.path.join(pasta_saida, 'dinucleotide_composition.csv')
    df_dinu.to_csv(csv_path, sep=';', decimal='.')
    print(f"  ✅ Tabela de composição de dinucleotídeos salva em: {csv_path}")

    # Gerar Heatmap
    print("  Gerando Heatmap de Dinucleotídeos...")
    status_queue.put(("progress", 90))
    plt.figure(figsize=(14, max(8, len(df_dinu) * 0.5)))
    sns.heatmap(df_dinu, annot=True, fmt=".2f", cmap="viridis", linewidths=0.5,
                cbar_kws={'label': 'Frequência (%)'})
    plt.title("Composição de Dinucleotídeos (%)", fontsize=16)
    plt.xlabel("Dinucleotídeo")
    plt.ylabel("Espécie")
    plt.tight_layout()
    
    plot_path = os.path.join(pasta_saida, 'dinucleotide_composition_heatmap.png')
    plt.savefig(plot_path, dpi=150, bbox_inches="tight")
    plt.close()
    
    print(f"  ✅ Heatmap salvo em: {plot_path}")
    status_queue.put(("image_ready", (plot_path, "Composição de Dinucleotídeos")))


# ######################################################################
# --- NOVAS FUNÇÕES DE ANÁLISE (14, 15, 16) ---
# ######################################################################

# --- Análise 14: Gráfico de Paridade PR2 (Feature c) ---

def calcular_pr2_por_gene(sequences, genetic_code_id):
    """
    Calcula os valores A3/(A3+T3) e G3/(G3+C3) para uma lista de genes (sequências).
    """
    genetic_code = GENETIC_CODE_TABLES.get(genetic_code_id, GENETIC_CODE_TABLES[1])
    stop_codons = {codon for codon, aa in genetic_code.items() if aa == '*'}
    results = []

    for seq in sequences:
        counts = {'A3': 0, 'T3': 0, 'G3': 0, 'C3': 0}
        
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i+3]
            if codon in stop_codons:
                continue # Ignora stop codons
                
            nuc3 = codon[2]
            if nuc3 in ['A', 'T', 'G', 'C']:
                counts[f"{nuc3}3"] += 1
        
        # Calcular frações
        a3_plus_t3 = counts['A3'] + counts['T3']
        g3_plus_c3 = counts['G3'] + counts['C3']
        
        a3_frac = (counts['A3'] / a3_plus_t3) if a3_plus_t3 > 0 else np.nan
        g3_frac = (counts['G3'] / g3_plus_c3) if g3_plus_c3 > 0 else np.nan
        
        if not (np.isnan(a3_frac) or np.isnan(g3_frac)):
            results.append({'A3_frac': a3_frac, 'G3_frac': g3_frac})

    return pd.DataFrame(results)

def analise_pr2_plot(lista_arquivos, pasta_saida, genetic_code_id, status_queue, gene_list=None):
    """(Análise 14) Gera o Gráfico de Paridade PR2 (A3/T3 vs G3/C3)."""
    print(f"\n=== ANÁLISE GRÁFICO DE PARIDADE PR2 ===")
    
    all_seqs_by_species = extrair_sequencias_cds(lista_arquivos, status_queue, gene_list)
    
    if not all_seqs_by_species:
        print("  ❌ Erro: Nenhuma sequência de CDS válida foi extraída.")
        return

    status_queue.put(("progress", 40))
    all_results_long = []

    for species_name, sequences in all_seqs_by_species.items():
        print(f"  Calculando PR2 para {species_name}...")
        status_queue.put(("message", f"Calculando PR2 para {species_name}..."))
        
        df_results = calcular_pr2_por_gene(sequences, genetic_code_id)
        if df_results.empty:
            continue
            
        df_results['species'] = species_name
        all_results_long.append(df_results)
    
    if not all_results_long:
        print("  ❌ Erro: Nenhum dado de PR2 pôde ser calculado.")
        return
        
    df_plot_long = pd.concat(all_results_long).dropna()
    df_plot_long.to_csv(os.path.join(pasta_saida, 'pr2_plot_data_per_gene.csv'), sep=';', index=False)
    
    print("  Gerando Gráfico PR2...")
    status_queue.put(("progress", 80))
    
    plt.figure(figsize=(12, 12))
    ax = sns.scatterplot(x='G3_frac', y='A3_frac', data=df_plot_long, hue='species', alpha=0.5, s=20)
    
    # Centro e linhas
    ax.axhline(0.5, color='black', linestyle='--', linewidth=1)
    ax.axvline(0.5, color='black', linestyle='--', linewidth=1)
    
    ax.set_xlabel('G3 / (G3 + C3)', fontsize=12)
    ax.set_ylabel('A3 / (A3 + T3)', fontsize=12)
    ax.set_title(f'Gráfico de Paridade PR2 (Viés da 3ª Posição)', fontsize=16)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
    
    plt.tight_layout()
    output_path = os.path.join(pasta_saida, 'pr2_plot.png')
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"\n✅ Gráfico PR2 salvo em: {output_path}")
    status_queue.put(("image_ready", (output_path, "Gráfico de Paridade PR2")))

# --- INÍCIO DO BLOCO DE ANÁLISE 15 (tAI) ATUALIZADO ---

def contar_tRNAs(lista_arquivos, status_queue):
    """
    Conta a ocorrência de cada anticódon nos arquivos GenBank.
    """
    print("  Analisando genes de tRNA (anticódons)...")
    anticodon_counts = Counter()
    total_tRNAs = 0
    
    # Regex para o qualifier /anticodon, ex: (pos:34..36,aa:Leu,seq:UAA)
    qualifier_regex = re.compile(r'seq:\s*([ATGCU]{3})', re.IGNORECASE)
    
    # Regex de fallback para o campo /note, ex: (anticodon: GGC)
    note_regex = re.compile(r'anticodon:\s*([ATGCU]{3})', re.IGNORECASE)
    
    for i, caminho_completo in enumerate(lista_arquivos):
        nome_base = os.path.basename(caminho_completo).split('.')[0]
        status_queue.put(("message", f"Contando tRNAs em {nome_base}..."))
        
        try:
            for record in SeqIO.parse(caminho_completo, "genbank"):
                for feature in record.features:
                    if feature.type == "tRNA":
                        anticodon_str = ""
                        
                        # --- MÉTODO 1 (Preferencial): Ler o qualifier /anticodon ---
                        if "anticodon" in feature.qualifiers:
                            qual_value = feature.qualifiers.get("anticodon")[0]
                            match = qualifier_regex.search(qual_value)
                            if match:
                                anticodon_str = match.group(1)
                        
                        # --- MÉTODO 2 (Fallback): Ler o qualifier /note ---
                        if not anticodon_str:
                            note = str(feature.qualifiers.get("note", ""))
                            match = note_regex.search(note)
                            if match:
                                anticodon_str = match.group(1)
                        
                        if anticodon_str:
                            anticodon_norm = anticodon_str.upper().replace('T', 'U')
                            if len(anticodon_norm) == 3 and all(b in 'ATGCU' for b in anticodon_norm):
                                # Armazena o anticódon 5'-3'
                                anticodon_counts[anticodon_norm] += 1
                                total_tRNAs += 1
                                    
        except Exception as e:
            print(f"  ❌ Erro ao contar tRNAs em {nome_base}: {e}")
            
    if total_tRNAs == 0:
         print("  ⚠️ Aviso: Nenhum gene de tRNA com anticódon reconhecível foi encontrado.")
         print("     A análise tAI pode não ser significativa (todos os pesos W serão 1.0).")
    else:
         print(f"  ✅ Encontrados {total_tRNAs} tRNAs com anticódons definidos.")
         
    return anticodon_counts

def calcular_pesos_W_wobble(anticodon_counts, genetic_code_id, wobble_matrix):
    """
    Calcula os pesos de adaptação (Wi) para cada códon, usando regras de
    pareamento wobble (dos Reis et al., 2004) FORNECIDAS.
    """
    print("  Calculando pesos (W) ponderados por wobble (dos Reis)...")
    genetic_code = GENETIC_CODE_TABLES.get(genetic_code_id, GENETIC_CODE_TABLES[1])
    W_values = {}
    
    # 1. Normalizar contagens de tRNA (abundância relativa)
    total_tRNAs = sum(anticodon_counts.values())
    if total_tRNAs == 0:
        print("  Aviso: Nenhuma contagem de tRNA. Usando pseudocontagens.")
        for aa, codons in AA_CODON_MAPS[genetic_code_id].items():
            if aa == '*': continue
            for codon in codons:
                rc_anticodon = str(Seq(codon).reverse_complement()).replace('T', 'U')
                anticodon_counts[rc_anticodon] = 1
        total_tRNAs = sum(anticodon_counts.values())

    relative_abundance = {ac: count / total_tRNAs for ac, count in anticodon_counts.items()}
    
    # 2. Iterar sobre todos os 61 códons codificantes
    for codon, aa in genetic_code.items():
        if aa == '*': continue
        
        codon_base_3 = codon[2]
        codon_bases_1_2 = codon[:2]
        
        sum_wi = 0.0
        
        # 3. Iterar sobre todos os tRNAs (anticódons) disponíveis
        for anticodon_5_3, n_j in relative_abundance.items():
            
            # 3a. Pareamento das posições 1 e 2 (Watson-Crick reverso)
            if codon_bases_1_2 != str(Seq(anticodon_5_3[1::-1]).complement()):
                continue 
            
            # 3b. Pareamento da posição 3 (Wobble)
            anticodon_base_1_raw = anticodon_5_3[0].upper().replace('T', 'U')
            
            # Aplicar modificação (ex: A -> I)
            anticodon_base_1_mod = ANTICODON_MODIFICATION_MAP.get(
                anticodon_base_1_raw, anticodon_base_1_raw
            )
            
            # Obter o valor S_ij da matriz de wobble FORNECIDA
            s_ij = wobble_matrix.get(
                (anticodon_base_1_mod, codon_base_3), 1.0 # 1.0 = sem pareamento (padrão)
            )
            
            # 3c. Adicionar à soma (Fórmula: Wi = Σ [ (1 - Sij) * nj ] )
            sum_wi += (1.0 - s_ij) * n_j
            
        W_values[codon] = sum_wi

    # 4. Normalizar os valores de W (dividir pelo W máximo)
    max_w = max(W_values.values())
    
    if max_w > 0:
        for codon in W_values:
            W_values[codon] /= max_w
    else:
        print("  Aviso: W máximo é 0. Todos os pesos de tAI serão 0.")

    min_val = 1e-9
    for codon in W_values:
        if W_values[codon] < min_val:
            W_values[codon] = min_val
            
    print("  ✅ Cálculo de pesos W concluído.")
    return W_values

def calcular_tAI_por_gene(sequences, W_values):
    """Calcula o tAI (média geométrica dos pesos W) para cada gene."""
    results = []
    
    for seq in sequences:
        log_w_sum = 0
        num_codons_validos = 0
        
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i+3]
            if codon in W_values:
                # W_values[codon] nunca será 0 por causa do piso
                log_w_sum += np.log(W_values[codon])
                num_codons_validos += 1
        
        if num_codons_validos > 0:
            tai = np.exp(log_w_sum / num_codons_validos)
            results.append({'tAI': tai, 'gene_length_codons': num_codons_validos})
        
    return pd.DataFrame(results)

def analise_tAI(lista_arquivos, pasta_saida, genetic_code_id, status_queue, gene_list=None, super_kingdom="Bactéria"):
    """
    (Análise 15 - NOVA) Calcula e plota o tAI Ponderado por Wobble.
    """
    print(f"\n=== ANÁLISE DE ADAPTAÇÃO AO tRNA (tAI Ponderado por Wobble) ===")
    
    # ####################################################################
    # ## AQUI ESTÁ A ESCOLHA (AGORA DINÂMICA) ##
    wobble_rules_map = {
        "Bactéria": WOBBLE_MATRIX_BACTERIA,
        "Eucarioto": WOBBLE_MATRIX_EUKARYA
    }
    
    # Seleciona a matriz correta, com fallback para Bactéria
    wobble_rules = wobble_rules_map.get(super_kingdom, WOBBLE_MATRIX_BACTERIA)
    print(f"  Usando regras de Wobble para: {super_kingdom}")
    # ####################################################################

    # 1. Contar tRNAs
    status_queue.put(("progress", 20))
    anticodon_counts = contar_tRNAs(lista_arquivos, status_queue)
    df_counts = pd.DataFrame(anticodon_counts.items(), columns=['Anticodon', 'Count']).sort_values(by='Count', ascending=False)
    df_counts.to_csv(os.path.join(pasta_saida, 'tRNA_anticodon_counts.csv'), sep=';', index=False)
    print(f"  ✅ Contagem de anticódons salva.")

    # 2. Calcular Pesos W (Wobble)
    status_queue.put(("progress", 30))
    W_values = calcular_pesos_W_wobble(anticodon_counts, genetic_code_id, wobble_rules)
    df_W = pd.DataFrame(W_values.items(), columns=['Codon', 'Weight_W']).sort_values(by='Weight_W', ascending=False)
    df_W.to_csv(os.path.join(pasta_saida, f'tAI_codon_weights_wobble_{super_kingdom}.csv'), sep=';', index=False)
    
    # 3. Extrair sequências (com filtro, se houver)
    all_seqs_by_species = extrair_sequencias_cds(lista_arquivos, status_queue, gene_list)
    if not all_seqs_by_species:
        print("  ❌ Erro: Nenhuma sequência de CDS válida foi extraída.")
        return

    # 4. Calcular tAI para cada gene em cada espécie
    status_queue.put(("progress", 60))
    all_results_long = []
    for species_name, sequences in all_seqs_by_species.items():
        print(f"  Calculando tAI (Wobble) para {species_name}...")
        status_queue.put(("message", f"Calculando tAI para {species_name}..."))
        
        df_results = calcular_tAI_por_gene(sequences, W_values)
        if df_results.empty:
            continue
            
        df_results['species'] = species_name
        all_results_long.append(df_results)
        
    if not all_results_long:
        print("  ❌ Erro: Nenhum dado de tAI pôde ser calculado.")
        return
        
    df_plot_long = pd.concat(all_results_long).dropna()
    df_plot_long.to_csv(os.path.join(pasta_saida, 'tAI_wobble_data_per_gene.csv'), sep=';', index=False)
    
    # 5. Gerar Boxplot
    print("  Gerando Box Plot de tAI (Wobble)...")
    status_queue.put(("progress", 90))
    plt.figure(figsize=(max(10, len(all_seqs_by_species)*1.5), 8))
    
    # --- INÍCIO DA MODIFICAÇÃO (Filtro de Piso) ---
    # Definir o "piso" de tAI (ligeiramente acima do valor mínimo de 1e-9)
    tAI_floor = 1.1e-9 
    
    # Filtrar os dados para plotagem
    df_filtered_plot = df_plot_long[df_plot_long['tAI'] > tAI_floor]
    
    plot_title_str = f'Distribuição tAI Ponderado por Wobble ({super_kingdom})'
    
    # Checar se o filtro removeu tudo
    if df_filtered_plot.empty:
        print("  ⚠️ Aviso: Todos os valores de tAI estão no piso mínimo (provavelmente 1e-9).")
        print("     Isso sugere que não há dados de tRNA nos arquivos GenBank.")
        print("     O gráfico será gerado com dados originais, mas pode parecer vazio/distorcido.")
        # Se estiver vazio, plotar os dados originais para não falhar
        df_to_plot = df_plot_long
        plot_title_str += "\n(Aviso: Todos os valores no piso mínimo)"
    else:
        df_to_plot = df_filtered_plot
        num_filtrados = len(df_plot_long) - len(df_filtered_plot)
        if num_filtrados > 0:
            print(f"  Info: Filtrados {num_filtrados} genes do gráfico (tAI < {tAI_floor:.1e}).")
            plot_title_str += f"\n(Valores < {tAI_floor:.1e} filtrados para clareza)"
    
    # Revertendo para boxplot + stripplot
    sns.boxplot(x='species', y='tAI', data=df_to_plot, palette="Set3")
    sns.stripplot(x='species', y='tAI', data=df_to_plot, color=".25", size=2, alpha=0.2)
    
    plt.title(plot_title_str, fontsize=16)
    plt.ylabel('tAI (Média Geométrica dos Pesos W)')
    plt.xlabel('Espécie')
    
    # Adicionar a escala logarítmica
    plt.yscale('log')
    
    # --- FIM DA MODIFICAÇÃO ---
    
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    
    output_path = os.path.join(pasta_saida, 'tAI_wobble_comparative_boxplot.png')
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"\n✅ Gráfico de tAI (Wobble) salvo em: {output_path}")
    status_queue.put(("image_ready", (output_path, "Análise tAI Ponderada (Wobble)")))

# --- FIM DO BLOCO DE ANÁLISE 15 (tAI) ATUALIZADO ---


# --- Análise 16: Análise de Motifs Upstream (Feature b) ---

def analise_motifs_upstream(lista_arquivos, pasta_saida, status_queue, gene_list=None, upstream_dist=200, kmer_size=6):
    """(Análise 16) Extrai regiões upstream, conta k-mers e plota os mais frequentes."""
    print(f"\n=== ANÁLISE DE MOTIFS UPSTREAM ===")
    print(f"  Parâmetros: Distância={upstream_dist}pb, K-mer={kmer_size}pb")
    
    total_files = len(lista_arquivos)
    
    for i, caminho_completo in enumerate(lista_arquivos):
        nome_base = os.path.basename(caminho_completo).split('.')[0]
        print(f"\n  Analisando motifs em {nome_base}...")
        status_queue.put(("message", f"Analisando motifs em {nome_base}..."))
        status_queue.put(("progress", int(10 + (i / total_files) * 80)))
        
        all_upstream_seqs = []
        genes_processados = 0
        
        try:
            # Iterar pelos contigs/registros no arquivo
            for record in SeqIO.parse(caminho_completo, "genbank"):
                record_seq_str = str(record.seq).upper()
                record_len = len(record_seq_str)
                
                for feature in record.features:
                    if feature.type != "CDS": continue
                    
                    # Aplicar filtro de gene
                    if not _apply_gene_filter(feature, gene_list):
                        continue
                        
                    seq_upstream = ""
                    try:
                        if feature.location.strand == 1:
                            # Gene na fita (+)
                            cds_start = feature.location.start
                            # Evitar ler antes do início do contig
                            upstream_start = max(0, cds_start - upstream_dist)
                            upstream_end = cds_start
                            if upstream_start < upstream_end:
                                seq_upstream = record_seq_str[upstream_start:upstream_end]
                                
                        elif feature.location.strand == -1:
                            # Gene na fita (-)
                            cds_end = feature.location.end
                            # Evitar ler além do fim do contig
                            upstream_start = cds_end
                            upstream_end = min(record_len, cds_end + upstream_dist)
                            if upstream_start < upstream_end:
                                # Extrair a fita (+) e fazer reverso complementar
                                seq_fwd = record_seq_str[upstream_start:upstream_end]
                                seq_upstream = str(Seq(seq_fwd).reverse_complement())
                        
                        if seq_upstream:
                            all_upstream_seqs.append(seq_upstream)
                            genes_processados += 1

                    except Exception as e_feature:
                        print(f"    Aviso: Pulando gene (erro de localização): {e_feature}")

            print(f"  Extraídas {len(all_upstream_seqs)} regiões upstream de {genes_processados} genes.")
            if not all_upstream_seqs:
                print("  ❌ Nenhuma região upstream encontrada.")
                continue

            # Contar k-mers
            kmer_counts = Counter()
            for seq in all_upstream_seqs:
                for j in range(len(seq) - kmer_size + 1):
                    kmer = seq[j:j+kmer_size]
                    if 'N' not in kmer and all(b in 'ATGC' for b in kmer):
                        kmer_counts[kmer] += 1
            
            if not kmer_counts:
                print("  ❌ Nenhum k-mer válido encontrado.")
                continue

            # Preparar dados para plot e CSV
            df_top_kmers = pd.DataFrame(kmer_counts.most_common(25), columns=['K-mer', 'Count'])
            df_top_kmers.to_csv(os.path.join(pasta_saida, f'motif_counts_{kmer_size}mer_{nome_base}.csv'), sep=';', index=False)
            
            # Gerar gráfico
            print("  Gerando gráfico de K-mers...")
            plt.figure(figsize=(10, 12))
            sns.barplot(x='Count', y='K-mer', data=df_top_kmers, palette='viridis')
            plt.title(f'Top 25 {kmer_size}-mers mais Frequentes (Upstream de {upstream_dist}pb)\n{nome_base}', fontsize=14)
            plt.xlabel('Contagem Absoluta', fontsize=12)
            plt.ylabel(f'{kmer_size}-mer', fontsize=12)
            plt.tight_layout()
            
            output_path = os.path.join(pasta_saida, f'motif_plot_{kmer_size}mer_{nome_base}.png')
            plt.savefig(output_path, dpi=150, bbox_inches='tight')
            plt.close()
            
            print(f"  ✅ Gráfico de Motifs salvo em: {output_path}")
            status_queue.put(("image_ready", (output_path, f"Motifs {kmer_size}-mer - {nome_base}")))

        except Exception as e_file:
            print(f"  ❌ Erro geral ao processar motifs em {nome_base}: {e_file}")

# ######################################################################
# --- NOVAS FERRAMENTAS DE BIOLOGIA SINTÉTICA (1, 2) ---
# ######################################################################

def otimizar_sequencia_codons(input_seq, host_rscu, genetic_code_id=1):
    """
    Otimização de Códons (Feature 1):
    Reescreve a sequência de entrada usando os códons ótimos do hospedeiro.
    """
    print("  Iniciando Otimização de Códons (Maximização)...")
    genetic_code = GENETIC_CODE_TABLES.get(genetic_code_id, GENETIC_CODE_TABLES[1])
    
    # 1. Obter o mapa de códons ótimos do hospedeiro
    optimal_map = _get_optimal_codon_map(host_rscu, genetic_code_id)
    
    # 2. Traduzir a sequência de entrada para AAs
    input_seq = input_seq.upper()
    aa_sequence = []
    stop_codon = ""
    for i in range(0, len(input_seq) - 2, 3):
        codon = input_seq[i:i+3]
        if codon not in genetic_code:
            aa_sequence.append('?') # Códon desconhecido
            continue
            
        aa = genetic_code[codon]
        if aa == '*':
            stop_codon = optimal_map.get('*', 'TAA') # Usa o "ótimo" stop
            break
        aa_sequence.append(aa)

    # 3. Re-traduzir AAs para DNA usando o mapa ótimo
    optimized_seq = []
    for aa in aa_sequence:
        if aa == '?':
            optimized_seq.append('NNN') # Preserva códons desconhecidos
        else:
            optimized_seq.append(optimal_map.get(aa, 'NNN'))
            
    # Adicionar o stop codon
    optimized_seq.append(stop_codon)
    
    print("  ✅ Otimização Concluída.")
    return "".join(optimized_seq)

def harmonizar_sequencia_codons(input_seq, host_counts, genetic_code_id=1):
    """
    Harmonização de Códons (Feature 2):
    Reescreve a sequência de entrada igualando o *ranking* de frequência
    dos códons do hospedeiro.
    """
    print("  Iniciando Harmonização de Códons (Ranking)...")
    genetic_code = GENETIC_CODE_TABLES.get(genetic_code_id, GENETIC_CODE_TABLES[1])
    
    # 1. Obter rankings de frequência do HOSPEDEIRO
    host_rank_map = _get_codon_frequency_rank(host_counts, genetic_code_id)
    
    # 2. Obter rankings de frequência da sequência de ENTRADA
    input_counts = _contar_codons_de_sequencia(input_seq, genetic_code_id)
    input_rank_map = _get_codon_frequency_rank(input_counts, genetic_code_id)
    
    # 3. Criar mapa de tradução (Input Rank -> Host Rank)
    # Ex: O 2º códon mais raro de Leucina no input será o 2º mais raro no host
    harmonization_map = {}
    for aa, input_ranked_codons in input_rank_map.items():
        if aa not in host_rank_map: continue
        
        host_ranked_codons = host_rank_map[aa]
        
        for i, input_codon in enumerate(input_ranked_codons):
            # Se o host tiver menos códons sinônimos (raro), usa o último
            host_index = min(i, len(host_ranked_codons) - 1)
            harmonization_map[input_codon] = host_ranked_codons[host_index]
            
    # 4. Iterar pela sequência de entrada e traduzir
    input_seq = input_seq.upper()
    harmonized_seq = []
    
    for i in range(0, len(input_seq) - 2, 3):
        codon = input_seq[i:i+3]
        if codon in harmonization_map:
            harmonized_seq.append(harmonization_map[codon])
        else:
            # Preserva códons não-mapeados (ex: Stop, NNN)
            harmonized_seq.append(codon)
            
    print("  ✅ Harmonização Concluída.")
    return "".join(harmonized_seq)

# ######################################################################
# --- NOVA ANÁLISE (17) ---
# ######################################################################

# --- Análise 17: Análise de MFE (Estrutura Secundária) ---

def analise_mfe_iniciacao(lista_arquivos, pasta_saida, genetic_code_id, status_queue, gene_list=None, mfe_region_length=50):
    """(Análise 17) Calcula o MFE da região 5' dos CDS."""
    print(f"\n=== ANÁLISE DE MFE DE INICIAÇÃO (5' UTR) ===")
    print(f"  Parâmetros: Região = Primeiros {mfe_region_length}pb")
    
    if not VIENNARNA_DISPONIVEL:
        print("  ❌ ERRO CRÍTICO: ViennaRNA (import RNA) não encontrado.")
        print("     A Análise 17 não pode ser executada. Instale 'viennarna' via Conda.")
        status_queue.put(("message", "Erro: ViennaRNA não encontrado."))
        return

    # 1. Extrair sequências de CDS
    all_seqs_by_species = extrair_sequencias_cds(lista_arquivos, status_queue, gene_list)
    if not all_seqs_by_species:
        print("  ❌ Erro: Nenhuma sequência de CDS válida foi extraída.")
        return

    status_queue.put(("progress", 40))
    all_results_long = []

    # 2. Calcular MFE para cada gene em cada espécie
    for species_name, sequences in all_seqs_by_species.items():
        print(f"  Calculando MFE para {species_name}...")
        status_queue.put(("message", f"Calculando MFE para {species_name}..."))
        
        try:
            df_results = _calcular_mfe_para_sequencias(sequences, mfe_region_length)
        except ImportError as e:
            # Captura o erro se o ViennaRNA falhar no meio
            print(f"  ❌ Erro: {e}")
            return
            
        if df_results.empty:
            print(f"  Aviso: Nenhum dado de MFE calculado para {species_name}.")
            continue
            
        df_results['species'] = species_name
        all_results_long.append(df_results)
        
    if not all_results_long:
        print("  ❌ Erro: Nenhum dado de MFE pôde ser calculado.")
        return
        
    df_plot_long = pd.concat(all_results_long).dropna()
    df_plot_long.to_csv(os.path.join(pasta_saida, 'mfe_data_per_gene.csv'), sep=';', index=False)
    
    # 3. Gerar Boxplot
    print("  Gerando Box Plot de MFE...")
    status_queue.put(("progress", 90))
    plt.figure(figsize=(max(10, len(all_seqs_by_species)*1.5), 8))
    
    sns.boxplot(x='species', y='mfe', data=df_plot_long, palette="coolwarm")
    sns.stripplot(x='species', y='mfe', data=df_plot_long, color=".25", size=2, alpha=0.2)
    
    plt.title(f'Distribuição de MFE (Energia Livre Mínima) - Região 5\' ({mfe_region_length}pb)', fontsize=16)
    plt.ylabel('MFE (kcal/mol) - (Mais negativo = mais estável/fechado)')
    plt.xlabel('Espécie')
    plt.xticks(rotation=45, ha='right')
    plt.grid(axis='y', linestyle='--', alpha=0.5)
    plt.tight_layout()
    
    output_path = os.path.join(pasta_saida, 'mfe_5prime_comparative_boxplot.png')
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"\n✅ Gráfico de MFE salvo em: {output_path}")
    status_queue.put(("image_ready", (output_path, "Análise MFE (5' UTR)")))

# ######################################################################
# --- NOVA ANÁLISE (18) ---
# ######################################################################

# --- Análise 18: Análise Comparativa de Conjuntos de Genes ---

def analise_comparativa_dois_grupos(lista_arquivos, pasta_saida, genetic_code_id, status_queue, gene_list_1, gene_list_2):
    """(Análise 18) Compara CUB (ENC, GC3, CAI) entre dois grupos de genes."""
    print(f"\n=== ANÁLISE COMPARATIVA DE DOIS GRUPOS ===")
    print(f"  Grupo 1: {len(gene_list_1)} genes")
    print(f"  Grupo 2: {len(gene_list_2)} genes")
    
    if not gene_list_1 or not gene_list_2:
        print("  ❌ Erro: Ambos os grupos de genes devem ser fornecidos.")
        return

    # 1. Calcular dados de viés do GENOMA INTEIRO para obter a tabela de referência 'W' para CAI
    print("  Calculando tabela de referência W (genoma completo)...")
    status_queue.put(("progress", 10))
    status_queue.put(("message", "Calculando referência (genoma)..."))
    all_bias_data_full_genome = processar_genomas_para_analise_vies(lista_arquivos, genetic_code_id, status_queue, gene_list=None)
    
    if not all_bias_data_full_genome:
        print("  ❌ Erro: Falha ao processar genomas de referência.")
        return
        
    all_results_g1 = []
    all_results_g2 = []
    all_counts_g1 = Counter()
    all_counts_g2 = Counter()
    
    status_queue.put(("progress", 40))
    
    # 2. Iterar sobre os arquivos e calcular métricas por gene para cada grupo
    for i, caminho_completo in enumerate(lista_arquivos):
        nome_base = os.path.basename(caminho_completo).split('.')[0]
        print(f"\n  Processando {nome_base}...")
        
        # Obter a tabela W específica da espécie
        if nome_base not in all_bias_data_full_genome:
            print(f"  Aviso: {nome_base} não encontrado nos dados de referência. Pulando.")
            continue
            
        w_reference_table = get_w_reference_table(
            all_bias_data_full_genome[nome_base]['rscu'], 
            genetic_code_id
        )
        
        # Calcular para Grupo 1
        status_queue.put(("message", f"Calculando Grupo 1 em {nome_base}..."))
        results_g1 = calcular_metricas_por_gene(caminho_completo, w_reference_table, genetic_code_id, gene_list_1)
        for res in results_g1:
            all_results_g1.append(res)
            all_counts_g1.update(res['counts'])
            
        # Calcular para Grupo 2
        status_queue.put(("message", f"Calculando Grupo 2 em {nome_base}..."))
        results_g2 = calcular_metricas_por_gene(caminho_completo, w_reference_table, genetic_code_id, gene_list_2)
        for res in results_g2:
            all_results_g2.append(res)
            all_counts_g2.update(res['counts'])

    if not all_results_g1 or not all_results_g2:
        print("  ❌ Erro: Nenhum gene correspondente encontrado em um ou ambos os grupos.")
        return
        
    df_g1 = pd.DataFrame(all_results_g1)
    df_g2 = pd.DataFrame(all_results_g2)
    
    # 3. Executar Testes Estatísticos
    print("\n--- Resultados Estatísticos (Mann-Whitney U) ---")
    
    metrics_to_test = ['enc', 'gc3', 'cai']
    p_values = {}
    for metric in metrics_to_test:
        data1 = df_g1[metric].dropna()
        data2 = df_g2[metric].dropna()
        if len(data1) > 0 and len(data2) > 0:
            stat, p = mannwhitneyu(data1, data2, alternative='two-sided')
            p_values[metric] = p
            print(f"  {metric.upper()}: p-valor = {p:.4e} (Mediana G1: {data1.median():.3f}, Mediana G2: {data2.median():.3f})")
        else:
            p_values[metric] = np.nan
            print(f"  {metric.upper()}: Dados insuficientes para o teste.")

    # 4. Teste Qui-Quadrado
    print("\n--- Resultados Estatísticos (Qui-Quadrado nas Contagens de Códons) ---")
    codons = sorted(list(set(all_counts_g1.keys()) | set(all_counts_g2.keys())))
    table = [
        [all_counts_g1.get(c, 0) for c in codons],
        [all_counts_g2.get(c, 0) for c in codons]
    ]
    
    try:
        chi2, p, dof, expected = chi2_contingency(table)
        print(f"  Qui-Quadrado (Contagens Totais): X²={chi2:.2f}, p-valor = {p:.4e}, DoF={dof}")
    except ValueError as e:
        print(f"  Erro no Qui-Quadrado: {e} (provavelmente contagens baixas)")
        p_values['chi2'] = np.nan
        
    # 5. Gerar Box Plots Comparativos
    print("\n  Gerando box plots comparativos...")
    status_queue.put(("progress", 80))
    
    df_g1['group'] = 'Grupo 1'
    df_g2['group'] = 'Grupo 2'
    df_plot = pd.concat([df_g1, df_g2])
    
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(20, 8))
    
    # Plot ENC
    sns.boxplot(x='group', y='enc', data=df_plot, ax=ax1, palette="Set2")
    sns.stripplot(x='group', y='enc', data=df_plot, ax=ax1, color=".25", size=3, alpha=0.3)
    ax1.set_title(f"ENC (p = {p_values.get('enc', 'N/A'):.2e})", fontsize=14)
    ax1.set_xlabel("")
    
    # Plot GC3
    sns.boxplot(x='group', y='gc3', data=df_plot, ax=ax2, palette="Set2")
    sns.stripplot(x='group', y='gc3', data=df_plot, ax=ax2, color=".25", size=3, alpha=0.3)
    ax2.set_title(f"GC3 (p = {p_values.get('gc3', 'N/A'):.2e})", fontsize=14)
    ax2.set_xlabel("")
    
    # Plot CAI
    sns.boxplot(x='group', y='cai', data=df_plot, ax=ax3, palette="Set2")
    sns.stripplot(x='group', y='cai', data=df_plot, ax=ax3, color=".25", size=3, alpha=0.3)
    ax3.set_title(f"CAI (p = {p_values.get('cai', 'N/A'):.2e})", fontsize=14)
    ax3.set_xlabel("")
    
    plt.suptitle(f"Comparação de Grupos (G1: {len(df_g1)} genes, G2: {len(df_g2)} genes)", fontsize=18)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    
    output_path = os.path.join(pasta_saida, 'comparacao_grupos_boxplot.png')
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"\n✅ Gráfico de Comparação de Grupos salvo em: {output_path}")
    status_queue.put(("image_ready", (output_path, "Comparação de Grupos de Genes")))

# ######################################################################
# --- NOVA ANÁLISE (19) ---
# ######################################################################

# --- Análise 19: Correlação com Dados de Expressão ---

def analise_correlacao_expressao(lista_arquivos, pasta_saida, genetic_code_id, status_queue, gene_list, expression_data, gene_col, expr_col):
    """(Análise 19) Correlaciona CUB (ENC, CAI) com dados de expressão."""
    print(f"\n=== ANÁLISE DE CORRELAÇÃO COM EXPRESSÃO ===")
    print(f"  Arquivo de Expressão: {len(expression_data)} linhas")
    print(f"  Coluna de Gene: '{gene_col}', Coluna de Expressão: '{expr_col}'")

    # 0. Preparar dados de expressão
    # Assegurar que a coluna de expressão seja numérica e > 0 para log
    try:
        expression_data[expr_col] = pd.to_numeric(expression_data[expr_col])
    except ValueError:
        print(f"  ❌ Erro: A coluna de expressão '{expr_col}' não é numérica.")
        return
        
    # Adicionar um pequeno "piso" para evitar log(0)
    min_expr = expression_data[expression_data[expr_col] > 0][expr_col].min()
    if pd.isna(min_expr): min_expr = 1e-3
    expression_data['expr_log'] = np.log10(expression_data[expr_col] + min_expr / 10)
    
    # 1. Calcular dados de viés do GENOMA INTEIRO para obter a tabela 'W'
    print("  Calculando tabela de referência W (genoma completo)...")
    status_queue.put(("progress", 10))
    status_queue.put(("message", "Calculando referência (genoma)..."))
    all_bias_data_full_genome = processar_genomas_para_analise_vies(lista_arquivos, genetic_code_id, status_queue, gene_list=None)
    
    if not all_bias_data_full_genome:
        print("  ❌ Erro: Falha ao processar genomas de referência.")
        return
        
    all_merged_data = []
    
    # 2. Iterar sobre os arquivos e calcular métricas por gene
    status_queue.put(("progress", 40))
    all_genes_in_gbk = set(gene_list) if gene_list else set()
    if not gene_list:
        # Se nenhum filtro global for aplicado, pegamos TODOS os genes
        print("  Coletando todos os genes dos arquivos GenBank...")
        for caminho_completo in lista_arquivos:
            try:
                for record in SeqIO.parse(caminho_completo, "genbank"):
                    for feature in record.features:
                        if feature.type != "CDS": continue
                        if "locus_tag" in feature.qualifiers: 
                            all_genes_in_gbk.add(feature.qualifiers["locus_tag"][0])
                        if "gene" in feature.qualifiers:
                            all_genes_in_gbk.add(feature.qualifiers["gene"][0])
            except Exception as e:
                print(f"  Aviso: Erro ao ler {caminho_completo} para lista de genes: {e}")
    
    print(f"  Total de {len(all_genes_in_gbk)} genes para análise.")

    for i, caminho_completo in enumerate(lista_arquivos):
        nome_base = os.path.basename(caminho_completo).split('.')[0]
        print(f"\n  Processando {nome_base}...")
        
        if nome_base not in all_bias_data_full_genome:
            continue
            
        w_reference_table = get_w_reference_table(
            all_bias_data_full_genome[nome_base]['rscu'], 
            genetic_code_id
        )
        
        # Calcular métricas para todos os genes (ou os filtrados)
        status_queue.put(("message", f"Calculando métricas em {nome_base}..."))
        results_genes = calcular_metricas_por_gene(
            caminho_completo, w_reference_table, genetic_code_id, all_genes_in_gbk
        )
        
        if results_genes:
            df_genes = pd.DataFrame(results_genes)
            df_genes['species'] = nome_base
            all_merged_data.append(df_genes)

    if not all_merged_data:
        print("  ❌ Erro: Nenhum gene com métricas foi calculado.")
        return

    df_plot_all_genes = pd.concat(all_merged_data)
    
    # 3. Mesclar dados de métricas com dados de expressão
    print("  Mesclando dados de métricas e expressão...")
    df_merged = pd.merge(
        df_plot_all_genes, 
        expression_data, 
        left_on='gene', 
        right_on=gene_col,
        how='inner' # Apenas genes presentes em AMBOS os conjuntos
    )
    
    if df_merged.empty:
        print(f"  ❌ Erro: Nenhum gene em comum encontrado entre os arquivos GenBank e o arquivo de expressão.")
        print(f"     Verifique se os identificadores na coluna '{gene_col}' correspondem aos 'locus_tag' ou 'gene' do GenBank.")
        return
        
    print(f"  ✅ {len(df_merged)} genes com dados de métrica e expressão encontrados.")
    
    # 4. Executar Testes de Correlação (Spearman Rank)
    print("\n--- Resultados Estatísticos (Correlação de Spearman) ---")
    
    corr_cai, p_cai = spearmanr(df_merged['cai'].dropna(), df_merged['expr_log'].dropna())
    print(f"  CAI vs Expressão (log10): Rho = {corr_cai:.3f}, p-valor = {p_cai:.4e}")
    
    corr_enc, p_enc = spearmanr(df_merged['enc'].dropna(), df_merged['expr_log'].dropna())
    print(f"  ENC vs Expressão (log10): Rho = {corr_enc:.3f}, p-valor = {p_enc:.4e}")
    
    # 5. Gerar Gráficos de Correlação
    print("\n  Gerando gráficos de correlação...")
    status_queue.put(("progress", 80))
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 8))
    
    # Plot CAI vs Expressão
    sns.regplot(x='expr_log', y='cai', data=df_merged, ax=ax1,
                scatter_kws={'alpha': 0.2, 's': 10}, 
                line_kws={'color': 'red'})
    ax1.set_title(f"CAI vs Expressão (log10)\nSpearman Rho = {corr_cai:.3f} (p = {p_cai:.2e})", fontsize=14)
    ax1.set_xlabel(f"Expressão (log10 {expr_col})")
    ax1.set_ylabel("CAI")
    
    # Plot ENC vs Expressão
    sns.regplot(x='expr_log', y='enc', data=df_merged, ax=ax2,
                scatter_kws={'alpha': 0.2, 's': 10}, 
                line_kws={'color': 'blue'})
    ax2.set_title(f"ENC vs Expressão (log10)\nSpearman Rho = {corr_enc:.3f} (p = {p_enc:.2e})", fontsize=14)
    ax2.set_xlabel(f"Expressão (log10 {expr_col})")
    ax2.set_ylabel("ENC")
    
    plt.suptitle(f"Correlação com Expressão (N = {len(df_merged)} genes)", fontsize=18)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    
    output_path = os.path.join(pasta_saida, 'correlacao_expressao.png')
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"\n✅ Gráfico de Correlação com Expressão salvo em: {output_path}")
    status_queue.put(("image_ready", (output_path, "Correlação com Expressão")))
