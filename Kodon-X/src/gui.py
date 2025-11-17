# gui.py
import os
import sys
import tkinter as tk # Ainda necessário para as variáveis (StringVar) e para o matplotlib
from tkinter import filedialog, messagebox # Manter os pop-ups de diálogo
import customtkinter as ctk # A nova biblioteca de widgets!
import threading
import queue
import glob
from datetime import datetime
import webbrowser

# Importa as bibliotecas de visualização
from PIL import Image, ImageTk
import matplotlib
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

# Essencial para threads com Tkinter.
matplotlib.use('Agg')

# Importa TODAS as funções de análise do nosso arquivo
# Isso permite que a GUI as chame pelo nome, como no script original
from analysis_backend import *


# ######################################################################
# --- BLOCO DA INTERFACE GRÁFICA (FRONTEND) - Versão CustomTkinter
# ######################################################################

# Configuração global do CustomTkinter
ctk.set_appearance_mode("Dark")  # "Dark", "Light", ou "System"
ctk.set_default_color_theme("blue")  # <--- MUDANÇA DE TEMA

class Redirector:
    """Redireciona stdout para a GUI"""
    def __init__(self, queue):
        self.queue = queue
    
    def write(self, string):
        self.queue.put(string)
    
    def flush(self):
        pass

class KodonE_GUI(ctk.CTk):
    """Interface gráfica profissional para análise de códons (Versão CustomTkinter)"""
    
    def __init__(self):
        super().__init__()
        
        self.title("Kodon-X") # <--- MUDANÇA DE NOME
        self.geometry("1440x850") # Um pouco maior para a nova interface
        
        # --- Paleta de Cores (Baseada no Tema) ---
        self.SUCCESS_COLOR = "#388E3C"
        self.ERROR_COLOR = "#D32F2F"
        self.WARNING_COLOR = "#F57C00"
        self.INFO_COLOR = "#0288D1"
        self.H3_COLOR = "#3498DB" # Cor azul estável para títulos de ajuda <--- MUDANÇA DE TEMA
        
        # --- Fontes Redondas (Tamanhos Ajustados) ---
        self.default_font = ctk.CTkFont(family="Calibri", size=14)
        self.default_bold_font = ctk.CTkFont(family="Calibri", size=14, weight="bold")
        self.console_font = ctk.CTkFont(family="Courier New", size=12)
        self.h1_font = ctk.CTkFont(family="Calibri", size=24, weight="bold")
        self.h2_font = ctk.CTkFont(family="Calibri", size=20, weight="bold")
        self.h3_font = ctk.CTkFont(family="Calibri", size=16, weight="bold")
        
        # --- Variáveis de Estado ---
        self.stdout_queue = queue.Queue()
        self.status_queue = queue.Queue()
        self.file_paths_map = {}
        self.file_checkboxes = {} # Novo: para substituir o Listbox
        self.current_images = []
        self.analysis_data = self._get_analysis_definitions()
        
        # Variáveis Tkinter (ctk usa as mesmas)
        self.input_folder_var = tk.StringVar()
        self.output_folder_var = tk.StringVar()
        self.genetic_code_var = tk.StringVar(value="1: Padrão (Universal)")
        self.filter_cds_var = tk.StringVar(value="ATG")
        
        self.super_kingdom_var = tk.StringVar(value="Bactéria")
        
        self.expression_dataframe = None
        self.expression_file_var = tk.StringVar(value="Nenhum arquivo carregado.")
        self.expr_gene_col_var = tk.StringVar(value="locus_tag")
        self.expr_value_col_var = tk.StringVar(value="TPM")
        
        # --- Configuração da Interface ---
        self._create_widgets()
        self._bind_events()
        
        # --- Redirecionar stdout ---
        self.stdout_original = sys.stdout
        sys.stdout = Redirector(self.stdout_queue)
        
        # Iniciar processamento de filas
        self._process_queues()
        
        # Configurar fechamento
        self.protocol("WM_DELETE_WINDOW", self._on_closing)
        
        # Mensagem de boas-vindas
        self._write_to_console("Bem-vindo ao Kodon-X v1.0!\n", "info") # <--- MUDANÇA DE NOME E VERSÃO
        self._write_to_console("Por favor, selecione sua pasta de entrada e clique em 'Carregar Arquivos'.\n\n")

    # --- INÍCIO DA GRANDE MUDANÇA: DESCRIÇÕES DETALHADAS ---
    def _get_analysis_definitions(self):
            """Retorna descrições detalhadas e requisitos das análises (sem alteração)"""
            return {
                '1: Estatística e CDS': {
                    'id': '1',
                    'description': ('Faz uma varredura inicial em todos os genomas. \n'
                                  'Aplicação em Pesquisa: Essencial para o *controle de qualidade (QC)* dos seus dados. Permite validar se os genomas baixados estão completos (Tamanho, Contigs) e se a anotação de genes (CDS) é consistente (ex: a maioria começa com "ATG").'),
                    'files_required': '1+', # 1 ou mais
                    'function': None 
                },
                '2: Listagem de Genes': {
                    'id': '2',
                    'description': ('Lista todos os genes, produtos e locus tags presentes em um único arquivo GenBank. \n'
                                  'Aplicação em Pesquisa: Útil para criar um "inventário" de genes. Permite verificar rapidamente se um gene de interesse (ex: "dnaA") está presente no genoma ou qual é o seu "locus_tag" (identificador único).'),
                    'files_required': '1', # Exatamente 1
                    'function': listar_genes_do_arquivo
                },
                '3: RSCU Heatmap Individual': {
                    'id': '3',
                    'description': ('Calcula e visualiza o RSCU (Relative Synonymous Codon Usage) para um único genoma. \n'
                                  'Aplicação em Pesquisa: É a visualização fundamental do viés de códons. Permite identificar quais códons são "preferenciais" (RSCU > 1.0) e quais são "raros" (RSCU < 1.0) para cada aminoácido. Também calcula métricas globais (ENC, GC3, CAI).'),
                    'files_required': '1', # Exatamente 1
                    'function': gerar_rscu_heatmap_e_tabela
                },
                '4: RSCU Comparativo': {
                    'id': '4',
                    'description': ('Agrupa genomas com base na similaridade de seu viés de códons (RSCU). \n'
                                  'Aplicação em Pesquisa: Usado para construir *árvores filogenéticas* baseadas em uso de códons (filogenômica). O PCA (Análise de Componentes Principais) revela quais códons mais contribuem para a variação entre grupos (ex: separar bactérias de arqueias).'),
                    'files_required': '2+', # 2 ou mais
                    'function': analise_rscu_comparativa
                },
                '5: Correlação RSCU (2 Genomas)': {
                    'id': '5',
                    'description': ('Calcula a correlação de Pearson (R) entre os padrões de RSCU de exatamente dois genomas. \n'
                                  'Aplicação em Pesquisa: Mede quantitativamente o quão similar é o viés entre dois organismos. Um R² alto (ex: 0.9) sugere que eles usam os códons de forma muito parecida, indicando possível proximidade evolutiva ou pressões de seleção similares.'),
                    'files_required': '2', # Exatamente 2
                    'function': analise_correlacao_rscu
                },
                '6: Histogramas RSCU Comparativos': {
                    'id': '6',
                    'description': ('Gera visualizações detalhadas do uso de códons, incluindo Box Plots por aminoácido e gráficos de barras empilhadas por espécie. \n'
                                  'Aplicação em Pesquisa: Permite uma análise "micro" do viés. Ajuda a ver quais aminoácidos (ex: Leucina, Arginina) sofrem mais viés e como cada espécie "escolhe" seus códons preferidos de forma diferente.'),
                    'files_required': '2+', # 2 ou mais
                    'function': gerar_histogramas_rscu
                },
                '7: Análise ENC vs GC3 (Wright Plot)': {
                    'id': '7',
                    'description': ("Gera o 'Gráfico de Wright' (ENC vs. GC3). \n"
                                  "Aplicação em Pesquisa: É a principal ferramenta para responder à pergunta: *'O viés de códons é causado por seleção natural (para eficiência) ou por pressão mutacional (viés de GC)?'*. Pontos na curva = mutação; pontos abaixo da curva = seleção."),
                    'files_required': '1+', # 1 ou mais
                    'function': analise_enc_gc3
                },
                '8: Composição Genômica': {
                    'id': '8',
                    'description': ('Análise detalhada da composição nucleotídica (A, T, G, C), conteúdo de GC total e tamanho do genoma para múltiplos arquivos. \n'
                                  'Aplicação em Pesquisa: Permite comparar características genômicas fundamentais. Útil para identificar outliers ou para correlacionar o conteúdo de GC geral com o viés de códons (GC3) e outras métricas.'),
                    'files_required': '1+', # 1 ou mais
                    'function': analise_composicao_genomica
                },
                '9: Códons Ótimos, Raros e CAI': {
                    'id': '9',
                    'description': ('Identifica códons preferenciais/raros e compara o CAI (Codon Adaptation Index) entre espécies. \n'
                                  'Aplicação em Pesquisa: O CAI é uma métrica de "quão otimizado" um genoma é. Genomas com CAI alto são geralmente de crescimento rápido ou alta expressão. A lista de códons ótimos é crucial para a *biologia sintética* (ver Aba 2).'),
                    'files_required': '1+', # 1 ou mais
                    'function': analise_codons_otimos_raros
                },
                '10: Viés de Pares de Códons (CPB)': {
                    'id': '10',
                    'description': ('Analisa a frequência de pares de códons adjacentes (ex: ATG-CGC). \n'
                                  'Aplicação em Pesquisa: A eficiência da tradução também depende da *interação* entre códons vizinhos (disponibilidade de tRNAs). Esta análise revela pares "rápidos" (favorecidos) e "lentos" (evitados), que são gargalos de tradução.'),
                    'files_required': '1+', # 1 ou mais
                    'function': analise_codon_pair_bias
                },
                '11: Análise Físico-Química (GRAVY & Aromo)': {
                    'id': '11',
                    'description': ('Calcula o escore GRAVY (hidropaticidade) e Aromo (aromaticidade) para todos os genes. \n'
                                  'Aplicação em Pesquisa: Compara o "perfil" químico geral dos proteomas. Permite testar hipóteses como: *'
                                  "'Organismos extremófilos (termofílicos) usam proteínas mais estáveis (ex: mais hidrofóbicas)?'*"),
                    'files_required': '1+', # 1 ou mais
                    'function': analise_gravy_aromo
                },
                '12: Gráfico de Neutralidade (GC12 vs GC3)': {
                    'id': '12',
                    'description': ('Plota o conteúdo de GC das posições 1+2 (GC12) contra o GC3. \n'
                                  'Aplicação em Pesquisa: Similar ao Gráfico de Wright (Análise 7), é outro método para separar a pressão de mutação da seleção. Se a seleção domina, a inclinação da reta (slope) se aproxima de 0. Se a mutação domina, a inclinação se aproxima de 1.'),
                    'files_required': '2+', # 2 ou mais para regressão
                    'function': analise_neutrality_plot
                },
                '13: Composição de Dinucleotídeos': {
                    'id': '13',
                    'description': ('Calcula e compara a frequência de todos os 16 dinucleotídeos (AA, AT, GC, etc.) no genoma completo. \n'
                                  'Aplicação em Pesquisa: É uma "impressão digital" do genoma, muito usada em virologia. Permite ver se um vírus (ex: Coronavírus) está se "camuflando" para ter uma assinatura de dinucleotídeos similar à do seu hospedeiro (humano).'),
                    'files_required': '1+', # 1 ou mais
                    'function': analise_dinucleotide_composition
                },
                '14: Gráfico de Paridade PR2 (A3/T3 vs G3/C3)': {
                    'id': '14',
                    'description': ('Plota A3/(A3+T3) vs G3/(G3+C3) para cada gene individualmente. \n'
                                  'Aplicação em Pesquisa: Analisa vieses de mutação e seleção que atuam de forma diferente nas bases (A/T vs G/C). Desvios do centro (0.5, 0.5) indicam que a seleção ou a mutação não está tratando todas as bases da 3ª posição de forma igual.'),
                    'files_required': '1+', # 1 ou mais
                    'function': analise_pr2_plot
                },
                '15: Índice de Adaptação ao tRNA (tAI)': {
                    'id': '15',
                    'description': ('Calcula o \'Índice de Adaptação ao tRNA\' (tAI) usando regras de *wobble* (padrão-ouro de dos Reis et al.). \n'
                                  'Aplicação em Pesquisa: Mede a *eficiência de tradução* teórica. Permite correlacionar o viés de códons com a abundância real de tRNAs no genoma, sendo um indicador mais forte de seleção translacional do que o CAI.'),
                    'files_required': '1+', # 1 ou mais
                    'function': analise_tAI 
                },
                '16: Análise de Motifs Upstream': {
                    'id': '16',
                    'description': ('Busca por sequências curtas (k-mers) sobrerrepresentadas nas regiões *antes* do início de cada gene (upstream). \n'
                                  'Aplicação em Pesquisa: Usado para descobrir *motifs* regulatórios, como sítios de ligação de fatores de transcrição (promotores, operadores) ou o sítio de ligação do ribossomo (Shine-Dalgarno em bactérias).'),
                    'files_required': '1+', # 1 ou mais
                    'function': analise_motifs_upstream
                },
                '17: Análise de MFE (Estrutura 5\')': {
                    'id': '17',
                    'description': ("Calcula a Energia Livre Mínima (MFE) da região 5' (início) dos genes. \n"
                                  "Aplicação em Pesquisa: Testa a hipótese da 'rampa traducional'. Estruturas de RNA muito fortes (MFE negativo) no início de um gene podem *bloquear* o ribossomo, diminuindo a eficiência da tradução."),
                    'files_required': '1+', # 1 ou mais
                    'function': analise_mfe_iniciacao
                },
                '18: Comparação de Grupos de Genes': {
                    'id': '18',
                    'description': ('Compara métricas de CUB (ENC, CAI) entre dois conjuntos de genes definidos pelo usuário (ex: genes de alta expressão vs. baixa expressão). \n'
                                  'Aplicação em Pesquisa: Permite testar hipóteses específicas. Ex: *"Genes ribossomais (Grupo 1) têm um viés (CAI) significativamente maior do que genes de membrana (Grupo 2)?"* O teste estatístico (Mann-Whitney U) responde a essa pergunta.'),
                    'files_required': '1+', # 1 ou mais
                    'function': analise_comparativa_dois_grupos
                },
                '19: Correlação com Expressão (RNA-Seq)': {
                    'id': '19',
                    'description': ('Correlaciona o viés (CAI, ENC) com dados de expressão (RNA-Seq, proteômica) fornecidos pelo usuário. \n'
                                  'Aplicação em Pesquisa: É o *teste de hipótese final*. Se a seleção translacional for a força motriz, você *espera* ver uma correlação forte: genes altamente expressos (ex: alto TPM) *devem* ter um viés maior (ex: alto CAI).'),
                    'files_required': '1+', # 1 ou mais
                    'function': analise_correlacao_expressao
                }
            }
    # --- FIM DA GRANDE MUDANÇA ---

    # --- Criação de Widgets ---

    def _create_widgets(self):
        """Cria todos os widgets da interface"""
        
        # --- Configuração do Grid Principal ---
        self.grid_columnconfigure(1, weight=3) # Painel principal (abas) é 3x maior
        self.grid_columnconfigure(0, weight=1) # Painel esquerdo
        self.grid_rowconfigure(0, weight=1)    # Linha principal
        self.grid_rowconfigure(1, weight=0)    # Barra de status
        
        # --- Painel Esquerdo (Configuração) ---
        left_frame = self._create_left_panel()
        left_frame.grid(row=0, column=0, padx=(10, 5), pady=10, sticky="nsew")
        
        # --- Painel Principal (Análise e Resultados) ---
        main_frame = self._create_main_panel()
        main_frame.grid(row=0, column=1, padx=(5, 10), pady=10, sticky="nsew")
        
        # --- Barra de Status ---
        self._create_status_bar()

    def _create_left_panel(self):
        """Cria o painel de configuração da esquerda"""
        left_panel = ctk.CTkFrame(self, corner_radius=10)
        
        # --- Passo 1: Arquivos de Entrada ---
        input_frame = ctk.CTkFrame(left_panel, corner_radius=10)
        input_frame.pack(fill=tk.X, padx=10, pady=(10, 15)) # <--- Mais espaçamento
        
        # <--- MUDANÇA DE EMOJI E FONTE ---
        ctk.CTkLabel(input_frame, text="📂 Passo 1: Arquivos de Entrada", font=ctk.CTkFont(size=16, weight="bold")).pack(anchor="w", padx=10, pady=(10, 5))
        
        # Frame para entrada e botão
        input_entry_frame = ctk.CTkFrame(input_frame, fg_color="transparent")
        input_entry_frame.pack(fill=tk.X, padx=10, pady=(0, 5))
        input_entry_frame.grid_columnconfigure(0, weight=1)
        
        input_entry = ctk.CTkEntry(input_entry_frame, textvariable=self.input_folder_var, corner_radius=8, font=self.default_font)
        input_entry.grid(row=0, column=0, sticky="ew", padx=(0, 5), pady=2)
        
        browse_input_btn = ctk.CTkButton(input_entry_frame, text="Procurar...", command=self._on_browse_input, width=80, corner_radius=8, font=self.default_bold_font)
        browse_input_btn.grid(row=0, column=1, padx=(5, 0), pady=2)
        
        load_files_btn = ctk.CTkButton(input_frame, text="Carregar Arquivos", command=self._on_load_files, corner_radius=8, font=self.default_bold_font)
        load_files_btn.pack(fill=tk.X, padx=10, pady=5)
        
        # --- SUBSTITUIÇÃO DO LISTBOX ---
        self.file_list_frame = ctk.CTkScrollableFrame(input_frame, corner_radius=8, height=250)
        self.file_list_frame.pack(fill=tk.X, expand=True, padx=10, pady=5)
        ctk.CTkLabel(self.file_list_frame, text="Nenhum arquivo carregado.", font=self.default_font).pack(padx=5, pady=5)
        # --- FIM DA SUBSTITUIÇÃO ---
        
        # Botões de seleção
        file_buttons = ctk.CTkFrame(input_frame, fg_color="transparent")
        file_buttons.pack(fill=tk.X, padx=10, pady=(5, 10))
        
        ctk.CTkButton(file_buttons, text="Todos", command=self._on_select_all, corner_radius=8, font=self.default_bold_font).pack(side=tk.LEFT, padx=2, expand=True)
        ctk.CTkButton(file_buttons, text="Nenhum", command=self._on_select_none, corner_radius=8, font=self.default_bold_font).pack(side=tk.LEFT, padx=2, expand=True)
        self.file_status_label = ctk.CTkLabel(file_buttons, text="0 arquivos", font=self.default_font)
        self.file_status_label.pack(side=tk.RIGHT, padx=5)

        # --- Passo 2: Pasta de Saída ---
        output_frame = ctk.CTkFrame(left_panel, corner_radius=10)
        output_frame.pack(fill=tk.X, padx=10, pady=15) # <--- Mais espaçamento
        
        # <--- MUDANÇA DE EMOJI E FONTE ---
        ctk.CTkLabel(output_frame, text="💾 Passo 2: Pasta de Saída", font=ctk.CTkFont(size=16, weight="bold")).pack(anchor="w", padx=10, pady=(10, 5))

        output_entry_frame = ctk.CTkFrame(output_frame, fg_color="transparent")
        output_entry_frame.pack(fill=tk.X, padx=10, pady=(0, 10))
        output_entry_frame.grid_columnconfigure(0, weight=1)
        
        output_entry = ctk.CTkEntry(output_entry_frame, textvariable=self.output_folder_var, corner_radius=8, font=self.default_font)
        output_entry.grid(row=0, column=0, sticky="ew", padx=(0, 5), pady=2)
        
        browse_output_btn = ctk.CTkButton(output_entry_frame, text="Procurar...", command=self._on_browse_output, width=80, corner_radius=8, font=self.default_bold_font)
        browse_output_btn.grid(row=0, column=1, padx=(5, 0), pady=2)

        # --- Passo 3: Tabela Genética ---
        genetic_frame = ctk.CTkFrame(left_panel, corner_radius=10)
        genetic_frame.pack(fill=tk.X, padx=10, pady=(15, 10)) # <--- Mais espaçamento
        
        # <--- MUDANÇA DE EMOJI E FONTE ---
        ctk.CTkLabel(genetic_frame, text="⚙️ Passo 3: Tabela Genética", font=ctk.CTkFont(size=16, weight="bold")).pack(anchor="w", padx=10, pady=(10, 5))
        
        genetic_combo = ctk.CTkComboBox(genetic_frame, variable=self.genetic_code_var, state="readonly", corner_radius=8,
                                     font=self.default_font,
                                     values=[
                                         "1: Padrão (Universal)", 
                                         "2: Mitocondrial de Vertebrados",
                                         "4: Mitocondrial de Mofos/Protozoários", 
                                         "11: Plasto de Bactérias/Plantas"
                                     ])
        genetic_combo.pack(fill=tk.X, padx=10, pady=(0, 10))

        return left_panel
    
    def _create_main_panel(self):
        """Cria o painel principal com abas (CTkTabview)"""
        main_panel = ctk.CTkFrame(self, fg_color="transparent")
        
        self.notebook = ctk.CTkTabview(main_panel, corner_radius=10, border_width=1)
        self.notebook.pack(fill=tk.BOTH, expand=True)
        
        # Adicionar abas
        analysis_tab = self.notebook.add("▶️ Executar Análise")
        synth_bio_tab = self.notebook.add("🛠️ Biologia Sintética")
        viz_tab = self.notebook.add("📈 Visualizador de Gráficos")
        console_tab = self.notebook.add("🖥️ Console de Saída")
        help_tab = self.notebook.add("❓ Ajuda")

        # --- Aba 1: Executar Análise ---
        self._create_analysis_tab(analysis_tab)
        
        # --- Aba 2: Biologia Sintética ---
        self._create_synth_bio_tab(synth_bio_tab)

        # --- Aba 3: Visualizador de Gráficos ---
        self._create_visualization_tab(viz_tab)

        # --- Aba 4: Console de Saída ---
        self._create_console_tab(console_tab)
        
        # --- Aba 5: Ajuda ---
        self._create_help_tab(help_tab)
        
        self.notebook.set("▶️ Executar Análise") # Definir aba inicial
        
        return main_panel

    def _create_analysis_tab(self, tab):
        """Cria a aba de seleção e execução de análise"""
        tab.grid_columnconfigure(0, weight=1)
        
        # Topo: Seleção da Análise
        ctk.CTkLabel(tab, text="Selecione a Análise:", font=ctk.CTkFont(size=16, weight="bold")).grid(row=0, column=0, sticky="w", padx=10, pady=(10, 5))
        
        self.analysis_combo = ctk.CTkComboBox(tab, state="readonly", font=self.default_font,
                                              values=list(self.analysis_data.keys()),
                                              command=self._on_analysis_selected, # Comando em vez de bind
                                              corner_radius=8)
        self.analysis_combo.grid(row=1, column=0, sticky="ew", padx=10, pady=5)
        
        # Meio: Descrição e Opções
        desc_frame = ctk.CTkFrame(tab, corner_radius=8)
        desc_frame.grid(row=2, column=0, sticky="ew", padx=10, pady=10)
        
        self.analysis_desc_text = ctk.CTkTextbox(desc_frame, height=120, wrap=tk.WORD, # <--- Mais altura
                                                 font=self.default_font, corner_radius=8)
        self.analysis_desc_text.pack(fill=tk.X, expand=True, padx=5, pady=5)
        self.analysis_desc_text.insert(1.0, "Selecione uma análise acima para ver sua descrição e requisitos.")
        self.analysis_desc_text.configure(state="disabled")

        self.analysis_options_frame = ctk.CTkFrame(tab, fg_color="transparent")
        self.analysis_options_frame.grid(row=3, column=0, sticky="ew", padx=10, pady=5)
        
        # --- Filtro de Subconjunto de Genes (Feature d) ---
        gene_filter_frame = ctk.CTkFrame(tab, corner_radius=8)
        gene_filter_frame.grid(row=4, column=0, sticky="nsew", padx=10, pady=10)
        tab.grid_rowconfigure(4, weight=1) # Fazer este frame crescer
        
        ctk.CTkLabel(gene_filter_frame, text="Filtro de Genes Global (Opcional)", font=ctk.CTkFont(size=14, weight="bold")).pack(anchor="w", padx=10, pady=5)
        ctk.CTkLabel(gene_filter_frame, 
                  text="Cole uma lista de locus_tags ou nomes de genes (um por linha).\nSe preenchido, TODAS as análises usarão apenas este subconjunto de genes.",
                  font=self.default_font).pack(anchor="w", padx=10, pady=(0, 5))
        
        self.gene_filter_text = ctk.CTkTextbox(gene_filter_frame, wrap=tk.WORD, 
                                               font=self.console_font, corner_radius=8) # Usar fonte monoespaçada
        self.gene_filter_text.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Baixo: Botão de Execução
        self.run_button = ctk.CTkButton(tab, text="🚀 Executar Análise Selecionada", 
                                     font=ctk.CTkFont(size=16, weight="bold"), # <--- Fonte maior
                                     command=self._on_run_analysis,
                                     corner_radius=8)
        self.run_button.grid(row=5, column=0, sticky="ew", padx=10, pady=(10, 10), ipady=15) # <--- Mais altura
        
        return tab

    def _create_synth_bio_tab(self, tab):
        """Cria a aba de ferramentas de Biologia Sintética"""
        
        # --- Variáveis de estado para esta aba ---
        self.synth_host_var = tk.StringVar()
        
        # Frame principal dividido (Grid de 2 colunas)
        tab.grid_columnconfigure(0, weight=1)
        tab.grid_columnconfigure(1, weight=1)
        tab.grid_rowconfigure(0, weight=1)

        # --- Painel Esquerdo (Configuração) ---
        left_frame = ctk.CTkFrame(tab, fg_color="transparent")
        left_frame.grid(row=0, column=0, sticky="nsew", padx=(10, 5), pady=10)
        left_frame.grid_rowconfigure(1, weight=1) # Fazer o input crescer
        left_frame.grid_columnconfigure(0, weight=1)
        
        # 1. Seleção do Hospedeiro
        host_frame = ctk.CTkFrame(left_frame, corner_radius=8)
        host_frame.grid(row=0, column=0, sticky="ew", pady=(0, 10))
        
        ctk.CTkLabel(host_frame, text="1. Selecione o Hospedeiro", font=ctk.CTkFont(size=14, weight="bold")).pack(anchor="w", padx=10, pady=5)
        ctk.CTkLabel(host_frame, text="Selecione um genoma carregado para ser o 'Hospedeiro':", font=self.default_font).pack(anchor="w", padx=10)
        
        self.synth_host_combo = ctk.CTkComboBox(host_frame, variable=self.synth_host_var, state="readonly", font=self.default_font, corner_radius=8)
        self.synth_host_combo.pack(fill=tk.X, padx=10, pady=(5, 10))
        
        # 2. Sequência de Entrada
        input_frame = ctk.CTkFrame(left_frame, corner_radius=8)
        input_frame.grid(row=1, column=0, sticky="nsew", pady=5)
        input_frame.grid_rowconfigure(1, weight=1)
        input_frame.grid_columnconfigure(0, weight=1)

        ctk.CTkLabel(input_frame, text="2. Cole a Sequência de DNA (CDS)", font=ctk.CTkFont(size=14, weight="bold")).grid(row=0, column=0, sticky="w", padx=10, pady=5)
        self.synth_input_seq = ctk.CTkTextbox(input_frame, wrap=tk.WORD, 
                                              font=self.console_font, corner_radius=8)
        self.synth_input_seq.grid(row=1, column=0, sticky="nsew", padx=5, pady=5)
        
        # 3. Botões de Ação
        button_frame = ctk.CTkFrame(left_frame, fg_color="transparent")
        button_frame.grid(row=2, column=0, sticky="ew", pady=(10, 0))
        button_frame.grid_columnconfigure((0,1), weight=1)
        
        self.run_optimize_btn = ctk.CTkButton(button_frame, text="Otimizar\n(Maximizar RSCU)", 
                                           command=self._on_run_optimization, corner_radius=8, font=self.default_bold_font)
        self.run_optimize_btn.grid(row=0, column=0, sticky="ew", padx=(0, 5), ipady=5)
        
        self.run_harmonize_btn = ctk.CTkButton(button_frame, text="Harmonizar\n(Manter Ranking)", 
                                            command=self._on_run_harmonization, corner_radius=8, font=self.default_bold_font)
        self.run_harmonize_btn.grid(row=0, column=1, sticky="ew", padx=(5, 0), ipady=5)


        # --- Painel Direito (Resultados) ---
        right_frame = ctk.CTkFrame(tab, fg_color="transparent")
        right_frame.grid(row=0, column=1, sticky="nsew", padx=(5, 10), pady=10)
        right_frame.grid_rowconfigure(0, weight=1)
        right_frame.grid_columnconfigure(0, weight=1)
        
        output_frame = ctk.CTkFrame(right_frame, corner_radius=8)
        output_frame.grid(row=0, column=0, sticky="nsew")
        output_frame.grid_rowconfigure(1, weight=1)
        output_frame.grid_columnconfigure(0, weight=1)
        
        ctk.CTkLabel(output_frame, text="Sequência Resultante", font=ctk.CTkFont(size=14, weight="bold")).grid(row=0, column=0, sticky="w", padx=10, pady=5)
        self.synth_output_seq = ctk.CTkTextbox(output_frame, wrap=tk.WORD, 
                                               font=self.console_font, corner_radius=8, state='disabled')
        self.synth_output_seq.grid(row=1, column=0, sticky="nsew", padx=5, pady=5)
        
        # Conectar o combobox ao evento de carregamento de arquivos
        self.bind("<<FilesLoaded>>", self._update_synth_host_list)
        return tab

    def _create_visualization_tab(self, tab):
        """Cria a aba de visualização de gráficos"""
        tab.grid_rowconfigure(0, weight=1)
        tab.grid_columnconfigure(0, weight=1)

        self.viz_canvas_frame = ctk.CTkFrame(tab, corner_radius=8)
        self.viz_canvas_frame.grid(row=0, column=0, sticky="nsew", padx=10, pady=(10, 5))
        
        # Figura para visualização
        self.fig = Figure(figsize=(10, 8), dpi=100)
        self.fig.patch.set_facecolor("#2B2B2B") # Cor de fundo escura
        self.fig.set_tight_layout(True)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.viz_canvas_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        toolbar_frame = ctk.CTkFrame(tab, corner_radius=8)
        toolbar_frame.grid(row=1, column=0, sticky="ew", padx=10, pady=(5, 10))
        self.toolbar = NavigationToolbar2Tk(self.canvas, toolbar_frame)
        self.toolbar.config(background="#2B2B2B") 
        self.toolbar._message_label.config(background="#2B2B2B", foreground="#DCE4EE")
        

        for widget in self.toolbar.winfo_children():
            try:
                # Tenta configurar botões, checkbuttons, e labels
                widget.config(background="#2B2B2B", fg="#DCE4EE")
            except tk.TclError:
                try:
                    # Tenta configurar apenas o fundo (para o frame do subplots)
                    widget.config(background="#2B2B2B")
                except tk.TclError:
                    # Ignora widgets que não podem ser configurados (ex: Separadores)
                    pass
            
        self.toolbar.update()
        
        self.image_status_label = ctk.CTkLabel(toolbar_frame, text="Nenhum gráfico carregado.", font=self.default_font)
        self.image_status_label.pack(side=tk.LEFT, padx=10)
        
        return tab

    def _create_console_tab(self, tab):
        """Cria a aba do console de saída"""
        tab.grid_rowconfigure(0, weight=1)
        tab.grid_columnconfigure(0, weight=1)

        self.console_text = ctk.CTkTextbox(tab, wrap=tk.WORD, 
                                           font=self.console_font, corner_radius=8)
        self.console_text.grid(row=0, column=0, sticky="nsew", padx=10, pady=10)
        self.console_text.configure(state='disabled')
        
        # Configurar tags para coloração (REMOVIDO o parâmetro 'font')
        self.console_text.tag_config("error", foreground=self.ERROR_COLOR)
        self.console_text.tag_config("success", foreground=self.SUCCESS_COLOR)
        self.console_text.tag_config("info", foreground=self.INFO_COLOR)
        self.console_text.tag_config("warning", foreground=self.WARNING_COLOR)
        
        return tab

    
    # Mover 'add_text' para ser um método da classe
    def add_text(self, help_frame, text, style="normal", wrap_length=900):
        """Função auxiliar para adicionar texto formatado à aba de Ajuda."""
        font = self.default_font
        text_color = None  
        pady = (0, 2)

        if style == "h1":
            font = self.h1_font
            pady = (5, 10)
        elif style == "h2":
            font = self.h2_font
            pady = (15, 8)
        elif style == "h3":
            font = self.h3_font
            text_color = self.H3_COLOR
            pady = (10, 5)
        elif style == "bold":
            font = self.default_bold_font
            pady = (0, 2)
        
        ctk.CTkLabel(help_frame, text=text, font=font, text_color=text_color, 
                     anchor="w", justify="left", wraplength=wrap_length).pack(fill="x", padx=10, pady=pady)
    
    def _create_help_tab(self, tab):
        """Cria a aba de ajuda com texto detalhado (Reconstruída com Labels)"""
        tab.grid_rowconfigure(0, weight=1)
        tab.grid_columnconfigure(0, weight=1)

        help_frame = ctk.CTkScrollableFrame(tab, corner_radius=8)
        help_frame.grid(row=0, column=0, sticky="nsew", padx=10, pady=10)

        # --- Conteúdo do Manual (agora usando a função add_text) ---
        
        self.add_text(help_frame, "Kodon-X v1.0 - GUIA DO USUÁRIO", 'h1')
        self.add_text(help_frame, "Bem-vindo à interface de análise Kodon-X. Siga este guia para um fluxo de trabalho eficiente.", 'normal')
        
        self.add_text(help_frame, "FLUXO DE TRABALHO RÁPIDO", 'h2')
        
        self.add_text(help_frame, "1. Carregar Arquivos", 'h3')
        self.add_text(help_frame, "Use o painel 'Passo 1' para 'Procurar...' sua pasta com arquivos .gbk ou .gb. Clique em 'Carregar Arquivos'.", 'normal')
        
        self.add_text(help_frame, "2. Selecionar Arquivos", 'h3')
        self.add_text(help_frame, "Na lista, marque os arquivos que deseja analisar ou use 'Todos'/'Nenhum'.", 'normal')
        
        self.add_text(help_frame, "3. Definir Saída", 'h3')
        self.add_text(help_frame, "Use o 'Passo 2' para escolher onde salvará os gráficos e tabelas .csv.", 'normal')
        
        self.add_text(help_frame, "4. Selecionar Análise", 'h3')
        self.add_text(help_frame, "Vá para a aba '▶️ Executar Análise' e escolha uma análise no menu. A descrição e os requisitos de arquivos aparecerão abaixo.", 'normal')
        
        self.add_text(help_frame, "5. Executar", 'h3')
        self.add_text(help_frame, "Clique no botão '🚀 Executar Análise Selecionada'.", 'normal')
        
        self.add_text(help_frame, "6. Ver Resultados", 'h3')
        self.add_text(help_frame, "Acompanhe o progresso no '🖥️ Console' e veja os gráficos gerados no '📈 Visualizador'.", 'normal')
        

        self.add_text(help_frame, "DESCRIÇÃO DAS ANÁLISES", 'h2')

        for analysis_name, data in self.analysis_data.items():
            self.add_text(help_frame, analysis_name, 'h3')
            self.add_text(help_frame, f"Requisito de Arquivos: {data['files_required']}", 'bold')
            self.add_text(help_frame, data['description'], 'normal')

        
        return tab
    # --- FIM DA RECONSTRUÇÃO DA ABA DE AJUDA ---

    def _create_status_bar(self):
        """Cria a barra de status inferior"""
        status_frame = ctk.CTkFrame(self, height=30, corner_radius=0)
        status_frame.grid(row=1, column=0, columnspan=2, sticky="ew", padx=0, pady=0)
        
        self.status_label = ctk.CTkLabel(status_frame, text="Pronto", anchor="w", font=self.default_font)
        self.status_label.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=10, pady=5)
        
        self.progress_bar = ctk.CTkProgressBar(status_frame, orientation='horizontal', mode='determinate')
        self.progress_bar.set(0)
        # O progress bar será adicionado com .pack() quando for usado

    # --- Bindings e Eventos ---

    def _bind_events(self):
        """Vincula eventos de widgets a funções"""
        # Atalhos
        self.bind('<Control-o>', lambda e: self._on_browse_input())
        self.bind('<Control-s>', lambda e: self._on_browse_output())
        self.bind('<Control-q>', lambda e: self._on_closing())
        self.bind('<F5>', lambda e: self._on_run_analysis())

    # --- Callbacks (Ações do Usuário) ---

    def _on_browse_input(self):
        folder = filedialog.askdirectory(title="Selecione a pasta de entrada com arquivos GenBank")
        if folder:
            self.input_folder_var.set(folder)
            self._on_load_files() # Carregar automaticamente
    
    def _on_browse_output(self):
        folder = filedialog.askdirectory(title="Selecione a pasta de saída para os resultados")
        if folder:
            self.output_folder_var.set(folder)
            
    def _on_load_files(self):
        folder = self.input_folder_var.get()
        if not folder:
            messagebox.showwarning("Aviso", "Selecione uma pasta de entrada primeiro.", parent=self)
            return
        
        # Limpar widgets antigos
        for widget in self.file_list_frame.winfo_children():
            widget.destroy()
        self.file_paths_map.clear()
        self.file_checkboxes.clear()
        
        patterns = ["*.gbk", "*.gb", "*.gbff"]
        files = []
        for pattern in patterns:
            files.extend(glob.glob(os.path.join(folder, pattern)))
        
        if not files:
            self._write_to_console(f"Nenhum arquivo .gbk ou .gb encontrado em '{folder}'\n", "warning")
            self.file_status_label.configure(text="0 arquivos")
            ctk.CTkLabel(self.file_list_frame, text="Nenhum arquivo .gbk/.gb encontrado.", font=self.default_font).pack(padx=5, pady=5)
            return
        
        self._write_to_console(f"Encontrados {len(files)} arquivos .gbk/.gb.\n", "success")
        
        for file_path in sorted(files):
            filename = os.path.basename(file_path)
            self.file_paths_map[filename] = file_path
            
            # Criar checkbox para cada arquivo
            var = tk.StringVar(value="on") # Ligar por padrão
            cb = ctk.CTkCheckBox(self.file_list_frame, text=filename, variable=var, onvalue="on", offvalue="off",
                                 font=ctk.CTkFont(size=11))
            cb.pack(fill=tk.X, padx=10, pady=2)
            self.file_checkboxes[filename] = (cb, var)
        
        self.file_status_label.configure(text=f"{len(files)} arquivos")
        
        self.event_generate("<<FilesLoaded>>")
    
    def _on_select_all(self):
        for cb, var in self.file_checkboxes.values():
            var.set("on")
    
    def _on_select_none(self):
        for cb, var in self.file_checkboxes.values():
            var.set("off")

    def _on_load_expression_file(self):
        """Callback para o botão 'Carregar Arquivo de Expressão'."""
        filepath = filedialog.askopenfilename(
            title="Selecione o arquivo de expressão (CSV ou TSV)",
            filetypes=[("Arquivos de Texto", "*.csv;*.tsv;*.txt"), ("Todos os Arquivos", "*.*")],
            parent=self
        )
        if not filepath: return
        try:
            with open(filepath, 'r') as f:
                first_line = f.readline()
                if '\t' in first_line: sep = '\t'
                elif ',' in first_line: sep = ','
                elif ';' in first_line: sep = ';'
                else: sep = r'\s+'
            
            self.expression_dataframe = pd.read_csv(filepath, sep=sep)
            filename = os.path.basename(filepath)
            self.expression_file_var.set(f"Carregado: {filename} ({len(self.expression_dataframe)} linhas)")
            self._write_to_console(f"Arquivo de expressão '{filename}' carregado com sucesso.\n", "success")
            self._write_to_console(f"Colunas detectadas: {list(self.expression_dataframe.columns)}\n")
            self.status_queue.put(("message", f"Arquivo de expressão '{filename}' carregado."))
        except Exception as e:
            self.expression_dataframe = None
            self.expression_file_var.set("Falha ao carregar o arquivo.")
            messagebox.showerror("Erro ao Carregar", f"Não foi possível ler o arquivo:\n{e}", parent=self)
            self._write_to_console(f"Falha ao carregar arquivo de expressão: {e}\n", "error")

    def _on_analysis_selected(self, selection):
        """Atualiza a descrição e as opções da análise"""
        if not selection: return
        
        data = self.analysis_data[selection]
        
        # Atualizar descrição
        self.analysis_desc_text.configure(state="normal")
        self.analysis_desc_text.delete(1.0, tk.END)
        
        # --- INÍCIO DA MUDANÇA: Usar a nova descrição detalhada ---
        desc = f"Descrição:\n{data['description']}\n\n"
        req = f"Requisito de Arquivos: {data['files_required']}"
        
        self.analysis_desc_text.insert(1.0, desc)
        self.analysis_desc_text.insert(tk.END, req, 'bold_req') # Tag para o requisito
        self.analysis_desc_text.configure(state="disabled")
        
        # --- INÍCIO DA CORREÇÃO ---
        self.analysis_desc_text.tag_config('bold_req', foreground=self.WARNING_COLOR) # Removido 'font'
        # --- FIM DA CORREÇÃO ---
        
        # Limpar opções antigas
        for widget in self.analysis_options_frame.winfo_children():
            widget.destroy()
        
        # Adicionar opções específicas (se houver)
        if data['id'] == '1':
            self.gene_filter_text.configure(state='normal')
            
            filter_frame = ctk.CTkFrame(self.analysis_options_frame, fg_color="transparent")
            filter_frame.pack(fill=tk.X, pady=5)
            ctk.CTkLabel(filter_frame, text="Filtro de Códon de Início:", font=self.default_font).pack(side=tk.LEFT, padx=(5,10))
            ctk.CTkEntry(filter_frame, textvariable=self.filter_cds_var, width=60, corner_radius=8, font=self.default_font).pack(side=tk.LEFT, padx=5)

        elif data['id'] == '15':
            self.gene_filter_text.configure(state='normal')
            options_frame = ctk.CTkFrame(self.analysis_options_frame, fg_color="transparent")
            options_frame.pack(fill=tk.X, pady=5)
            ctk.CTkLabel(options_frame, text="Super-reino (Regras de Wobble):", font=self.default_font).pack(side=tk.LEFT, padx=(5,10))
            kingdom_combo = ctk.CTkComboBox(options_frame, variable=self.super_kingdom_var, state="readonly", values=["Bactéria", "Eucarioto"], corner_radius=8, font=self.default_font)
            kingdom_combo.pack(side=tk.LEFT, padx=5)

        elif data['id'] == '16':
            self.gene_filter_text.configure(state='normal')
            self.upstream_dist_var = tk.StringVar(value="200")
            self.kmer_size_var = tk.StringVar(value="6")
            options_frame = ctk.CTkFrame(self.analysis_options_frame, fg_color="transparent")
            options_frame.pack(fill=tk.X, pady=5)
            ctk.CTkLabel(options_frame, text="Distância Upstream (pb):", font=self.default_font).pack(side=tk.LEFT, padx=(5,10))
            ctk.CTkEntry(options_frame, textvariable=self.upstream_dist_var, width=60, corner_radius=8, font=self.default_font).pack(side=tk.LEFT, padx=5)
            ctk.CTkLabel(options_frame, text="Tamanho do K-mer:", font=self.default_font).pack(side=tk.LEFT, padx=(15,10))
            ctk.CTkEntry(options_frame, textvariable=self.kmer_size_var, width=40, corner_radius=8, font=self.default_font).pack(side=tk.LEFT, padx=5)
        
        elif data['id'] == '17':
            self.gene_filter_text.configure(state='normal')
            self.mfe_region_var = tk.StringVar(value="50")
            options_frame = ctk.CTkFrame(self.analysis_options_frame, fg_color="transparent")
            options_frame.pack(fill=tk.X, pady=5)
            ctk.CTkLabel(options_frame, text="Região 5' (pb):", font=self.default_font).pack(side=tk.LEFT, padx=(5,10))
            ctk.CTkEntry(options_frame, textvariable=self.mfe_region_var, width=60, corner_radius=8, font=self.default_font).pack(side=tk.LEFT, padx=5)
        
        elif data['id'] == '18':
            self.gene_filter_text.configure(state='disabled')
            self._write_to_console("Filtro global desabilitado. Use as caixas 'Grupo 1' e 'Grupo 2'.\n", "info")
            
            # Frame para Grupo 1
            group1_frame = ctk.CTkFrame(self.analysis_options_frame, corner_radius=8)
            group1_frame.pack(fill=tk.X, pady=(5,10), expand=True)
            ctk.CTkLabel(group1_frame, text="Grupo de Genes 1 (um por linha):", font=self.default_bold_font).pack(anchor="w", padx=10, pady=5)
            self.gene_list_1_text = ctk.CTkTextbox(group1_frame, height=80, wrap=tk.WORD, 
                                                   font=self.console_font, corner_radius=8)
            self.gene_list_1_text.pack(fill=tk.X, expand=True, padx=5, pady=(0, 5))
            
            # Frame para Grupo 2
            group2_frame = ctk.CTkFrame(self.analysis_options_frame, corner_radius=8)
            group2_frame.pack(fill=tk.X, pady=5, expand=True)
            ctk.CTkLabel(group2_frame, text="Grupo de Genes 2 (um por linha):", font=self.default_bold_font).pack(anchor="w", padx=10, pady=5)
            self.gene_list_2_text = ctk.CTkTextbox(group2_frame, height=80, wrap=tk.WORD, 
                                                   font=self.console_font, corner_radius=8)
            self.gene_list_2_text.pack(fill=tk.X, expand=True, padx=5, pady=(0, 5))

        elif data['id'] == '19':
            self.gene_filter_text.configure(state='normal')
            expr_frame = ctk.CTkFrame(self.analysis_options_frame, corner_radius=8)
            expr_frame.pack(fill=tk.X, pady=5, expand=True)
            
            ctk.CTkLabel(expr_frame, text="Configuração de Expressão", font=self.default_bold_font).pack(anchor="w", padx=10, pady=5)

            btn_frame = ctk.CTkFrame(expr_frame, fg_color="transparent")
            btn_frame.pack(fill=tk.X, padx=5, pady=5)
            ctk.CTkButton(btn_frame, text="Carregar Arquivo de Expressão (.csv/.tsv)...", 
                       command=self._on_load_expression_file, corner_radius=8, font=self.default_bold_font).pack(side=tk.LEFT, padx=5)
            ctk.CTkLabel(btn_frame, textvariable=self.expression_file_var, 
                         font=ctk.CTkFont(size=11, slant="italic")).pack(side=tk.LEFT, padx=10)
            
            cols_frame = ctk.CTkFrame(expr_frame, fg_color="transparent")
            cols_frame.pack(fill=tk.X, padx=5, pady=(5,10))
            ctk.CTkLabel(cols_frame, text="Coluna de Gene:", font=self.default_font).pack(side=tk.LEFT, padx=5)
            ctk.CTkEntry(cols_frame, textvariable=self.expr_gene_col_var, width=120, corner_radius=8, font=self.default_font).pack(side=tk.LEFT, padx=5)
            ctk.CTkLabel(cols_frame, text="Coluna de Expressão:", font=self.default_font).pack(side=tk.LEFT, padx=15)
            ctk.CTkEntry(cols_frame, textvariable=self.expr_value_col_var, width=100, corner_radius=8, font=self.default_font).pack(side=tk.LEFT, padx=5)

        else:
            self.gene_filter_text.configure(state='normal')

    def _on_run_analysis(self):
        """Inicia a validação e a thread de análise"""
        
        files = self._get_selected_files()
        output_folder = self.output_folder_var.get()
        analysis_name = self.analysis_combo.get()
        analysis_data = self.analysis_data.get(analysis_name, {})
        
        gene_list = None
        extra_args = {}
        
        if analysis_data.get('id') != '18':
            gene_list_raw = self.gene_filter_text.get(1.0, tk.END)
            gene_list = set(line.strip() for line in gene_list_raw.splitlines() if line.strip())
            if not gene_list:
                gene_list = None
            else:
                self._write_to_console(f"Aplicando filtro global de {len(gene_list)} genes.\n", "info")

        if not analysis_name or analysis_name == "":
            messagebox.showwarning("Aviso", "Selecione uma análise para executar.", parent=self)
            return
        if not files:
            messagebox.showwarning("Aviso", "Selecione pelo menos um arquivo na lista.", parent=self)
            return
        if not output_folder or not os.path.isdir(output_folder):
            messagebox.showwarning("Aviso", "Selecione uma pasta de saída válida.", parent=self)
            return
            
        req = analysis_data.get('files_required', '1+')
        num_files = len(files)
        
        if req == '1' and num_files != 1:
            messagebox.showerror("Erro", f"A análise '{analysis_name}' requer EXATAMENTE 1 arquivo. Você selecionou {num_files}.", parent=self)
            return
        if req == '2' and num_files != 2:
            messagebox.showerror("Erro", f"A análise '{analysis_name}' requer EXATAMENTE 2 arquivos. Você selecionou {num_files}.", parent=self)
            return
        if req == '2+' and num_files < 2:
            messagebox.showerror("Erro", f"A análise '{analysis_name}' requer 2 OU MAIS arquivos. Você selecionou {num_files}.", parent=self)
            return
            
        try:
            if analysis_data.get('id') == '15':
                extra_args['super_kingdom'] = self.super_kingdom_var.get()
            elif analysis_data.get('id') == '16':
                extra_args['upstream_dist'] = int(self.upstream_dist_var.get())
                extra_args['kmer_size'] = int(self.kmer_size_var.get())
            elif analysis_data.get('id') == '17':
                extra_args['mfe_region_length'] = int(self.mfe_region_var.get())
            elif analysis_data.get('id') == '18':
                gene_list_1 = set(line.strip() for line in self.gene_list_1_text.get(1.0, tk.END).splitlines() if line.strip())
                gene_list_2 = set(line.strip() for line in self.gene_list_2_text.get(1.0, tk.END).splitlines() if line.strip())
                if not gene_list_1 or not gene_list_2:
                    messagebox.showerror("Erro", "A Análise 18 requer que AMBOS os 'Grupo 1' e 'Grupo 2' sejam preenchidos.", parent=self)
                    return
                extra_args['gene_list_1'] = gene_list_1
                extra_args['gene_list_2'] = gene_list_2
                gene_list = None
            elif analysis_data.get('id') == '19':
                gene_col = self.expr_gene_col_var.get()
                expr_col = self.expr_value_col_var.get()
                if self.expression_dataframe is None:
                    messagebox.showerror("Erro", "Carregue um arquivo de expressão para a Análise 19.", parent=self)
                    return
                if not gene_col or not expr_col:
                    messagebox.showerror("Erro", "Especifique os nomes das colunas de 'Gene' e 'Expressão'.", parent=self)
                    return
                if gene_col not in self.expression_dataframe.columns:
                    messagebox.showerror("Erro", f"Coluna de Gene '{gene_col}' não encontrada no arquivo de expressão.", parent=self)
                    return
                if expr_col not in self.expression_dataframe.columns:
                    messagebox.showerror("Erro", f"Coluna de Expressão '{expr_col}' não encontrada no arquivo de expressão.", parent=self)
                    return
                extra_args['expression_data'] = self.expression_dataframe
                extra_args['gene_col'] = gene_col
                extra_args['expr_col'] = expr_col
        except ValueError:
            messagebox.showerror("Erro", "Opções de análise inválidas. Verifique se os números são inteiros.", parent=self)
            return
        
        self._start_analysis_thread(files, output_folder, analysis_name, gene_list, extra_args)

    def _update_synth_host_list(self, event=None):
        """Atualiza a lista do combobox de hospedeiros."""
        hosts = sorted(self.file_paths_map.keys())
        self.synth_host_combo.configure(values=hosts)
        if hosts:
            self.synth_host_var.set(hosts[0])

    def _get_host_data(self, host_filename):
        """Função auxiliar para obter dados de viés (RSCU, Contagens)"""
        if host_filename not in self.file_paths_map:
            print(f"❌ Erro: Arquivo hospedeiro '{host_filename}' não encontrado no mapa.")
            return None
        host_filepath = self.file_paths_map[host_filename]
        genetic_code_id = self._get_genetic_code_id()
        print(f"  Calculando dados de viés para o hospedeiro: {host_filename}...")
        self.status_queue.put(("message", f"Calculando dados de {host_filename}..."))
        all_data = processar_genomas_para_analise_vies(
            [host_filepath], 
            genetic_code_id, 
            self.status_queue, 
            gene_list=None
        )
        if not all_data:
            print(f"❌ Falha ao processar dados do hospedeiro: {host_filename}")
            return None
        return all_data[list(all_data.keys())[0]]

    def _on_run_optimization(self):
        self._start_optimization_thread(mode="optimize")

    def _on_run_harmonization(self):
        self._start_optimization_thread(mode="harmonize")

    def _start_optimization_thread(self, mode="optimize"):
        """Valida e inicia a thread da ferramenta de Biologia Sintética."""
        host_filename = self.synth_host_var.get()
        input_seq = self.synth_input_seq.get(1.0, tk.END).strip().replace("\n", "")
        
        if not host_filename:
            messagebox.showwarning("Aviso", "Selecione um genoma hospedeiro.", parent=self)
            return
        if len(input_seq) < 10:
            messagebox.showwarning("Aviso", "Insira uma sequência de DNA válida (mínimo 10bp).", parent=self)
            return
            
        self.notebook.set("🖥️ Console de Saída")
        self.run_optimize_btn.configure(state='disabled')
        self.run_harmonize_btn.configure(state='disabled')
        self.progress_bar.pack(side=tk.RIGHT, padx=10, pady=5)
        self.progress_bar.set(0)
        self.progress_bar.start()
        self.status_label.configure(text=f"Iniciando {mode}...")
        
        genetic_code_id = self._get_genetic_code_id()

        thread = threading.Thread(
            target=self._optimization_thread_target,
            args=(input_seq, host_filename, genetic_code_id, mode),
            daemon=True
        )
        thread.start()

    def _optimization_thread_target(self, input_seq, host_filename, genetic_code_id, mode):
        """A função que executa na thread (backend) para otimização/harmonização."""
        try:
            print(f"\n{'='*60}")
            print(f"🚀 INICIANDO FERRAMENTA DE BIOLOGIA SINTÉTICA")
            print(f"   Modo: {mode}")
            print(f"   Hospedeiro: {host_filename}")
            print(f"   Tamanho Input: {len(input_seq)}bp")
            print(f"{'='*60}")
            
            host_data = self._get_host_data(host_filename)
            if not host_data:
                raise Exception("Não foi possível obter os dados do hospedeiro.")
            
            self.status_queue.put(("progress", 50))
            
            result_seq = ""
            if mode == "optimize":
                self.status_queue.put(("message", "Otimizando sequência..."))
                result_seq = otimizar_sequencia_codons(input_seq, host_data['rscu'], genetic_code_id)
            elif mode == "harmonize":
                self.status_queue.put(("message", "Harmonizando sequência..."))
                result_seq = harmonizar_sequencia_codons(input_seq, host_data['counts'], genetic_code_id)
            
            self.status_queue.put(("optimization_complete", result_seq))
            print("✅ Ferramenta concluída com sucesso.")
        except Exception as e:
            print(f"\n❌ ERRO DURANTE A OPERAÇÃO: {e}", "error")
            import traceback
            print(traceback.format_exc())
            self.status_queue.put(("done", None)) 
        finally:
            self.status_queue.put(("tool_done", None))

    # --- Lógica de Thread e Análise ---

    def _start_analysis_thread(self, files, output_folder, analysis_name, gene_list, extra_args):
        """Prepara e inicia a thread de análise de backend"""
        
        self.console_text.configure(state="normal")
        self.console_text.delete(1.0, tk.END)
        self.console_text.configure(state="disabled")
        self.notebook.set("🖥️ Console de Saída")
        
        self.run_button.configure(state='disabled')
        self.progress_bar.pack(side=tk.RIGHT, padx=10, pady=5)
        self.progress_bar.set(0)
        self.progress_bar.configure(mode='indeterminate') # Mudar para modo indeterminado
        self.progress_bar.start() 
        self.status_label.configure(text=f"Iniciando: {analysis_name}...")
        
        genetic_code_id = self._get_genetic_code_id()
        analysis_data = self.analysis_data[analysis_name]
        
        thread = threading.Thread(
            target=self._analysis_thread_target,
            args=(files, output_folder, analysis_name, analysis_data, genetic_code_id, gene_list, extra_args),
            daemon=True
        )
        thread.start()

    def _analysis_thread_target(self, files, output_folder, analysis_name, analysis_data, genetic_code_id, gene_list, extra_args):
            """A função que executa na thread (backend) (Lógica idêntica ao original)"""
            try:
                print(f"🚀 INICIANDO ANÁLISE: {analysis_name}")
                if gene_list:
                    print(f"FILTRO ATIVO: Analisando apenas {len(gene_list)} genes especificados.")
                print(f"{'='*60}")
                print(f"Arquivos: {len(files)}")
                print(f"Pasta de Saída: {output_folder}")
                print(f"Tabela Genética: {genetic_code_id}")
                print(f"{'='*60}")
                
                # --- Roteamento da Análise (Sem alterações) ---
                if analysis_data['id'] in ['3', '4', '5', '6', '7', '9', '12']:
                    all_bias_data = processar_genomas_para_analise_vies(files, genetic_code_id, self.status_queue, gene_list)
                    if not all_bias_data:
                        raise Exception("Falha ao processar dados de viés. Verifique os arquivos de entrada.")
                    analysis_function = analysis_data['function']
                    if analysis_function:
                        if analysis_data['id'] in ['3', '6']: 
                            analysis_function(all_bias_data, output_folder, genetic_code_id, self.status_queue)
                        else:
                            analysis_function(all_bias_data, output_folder, self.status_queue)
                elif analysis_data['id'] == '1':
                    print("\n--- [Análise 1] Parte 1: Estatísticas Agregadas ---")
                    df_stats = processar_gbk_agregado(files, self.status_queue)
                    if not df_stats.empty:
                        self._print_dataframe_limitado(df_stats, "Estatísticas Agregadas")
                        df_stats.to_csv(os.path.join(output_folder, 'estatisticas_genoma.csv'), index=False, sep=';')
                    print("\n--- [Análise 1] Parte 2: Análise de CDS ---")
                    filtro = self.filter_cds_var.get()
                    df_cds = analisar_gbk_cds(files, filtro, self.status_queue, gene_list)
                    if not df_cds.empty:
                        self._print_dataframe_limitado(df_cds, "Análise de CDS")
                        df_cds.to_csv(os.path.join(output_folder, 'analise_cds.csv'), index=False, sep=';')
                elif analysis_data['id'] == '2':
                    df_genes = listar_genes_do_arquivo(files[0], self.status_queue, gene_list)
                    if not df_genes.empty:
                        self._print_dataframe_limitado(df_genes, "Lista de Genes")
                        df_genes.to_csv(os.path.join(output_folder, f"lista_genes_{os.path.basename(files[0])}.csv"), index=False, sep=';')
                elif analysis_data['id'] == '8':
                    analise_composicao_genomica(files, output_folder, self.status_queue)
                elif analysis_data['id'] == '10':
                    analise_codon_pair_bias(files, output_folder, genetic_code_id, self.status_queue, gene_list)
                elif analysis_data['id'] == '11':
                    analise_gravy_aromo(files, output_folder, genetic_code_id, self.status_queue, gene_list)
                elif analysis_data['id'] == '13':
                    analise_dinucleotide_composition(files, output_folder, self.status_queue)
                elif analysis_data['id'] == '14':
                    analise_pr2_plot(files, output_folder, genetic_code_id, self.status_queue, gene_list)
                elif analysis_data['id'] == '15':
                    kingdom_name = extra_args.get('super_kingdom', 'Bactéria')
                    analise_tAI(files, output_folder, genetic_code_id, self.status_queue, gene_list, 
                                super_kingdom=kingdom_name)
                elif analysis_data['id'] == '16':
                    analise_motifs_upstream(files, output_folder, self.status_queue, gene_list, 
                                            extra_args['upstream_dist'], extra_args['kmer_size'])
                elif analysis_data['id'] == '17':
                    analise_mfe_iniciacao(files, output_folder, genetic_code_id, self.status_queue, gene_list, 
                                          extra_args['mfe_region_length'])
                elif analysis_data['id'] == '18':
                    analise_comparativa_dois_grupos(files, output_folder, genetic_code_id, self.status_queue, 
                                                    extra_args['gene_list_1'], extra_args['gene_list_2'])
                elif analysis_data['id'] == '19':
                    analise_correlacao_expressao(files, output_folder, genetic_code_id, self.status_queue, 
                                                 gene_list, extra_args['expression_data'], 
                                                 extra_args['gene_col'], extra_args['expr_col'])
                else:
                    print(f"Aviso: Análise '{analysis_name}' não implementada.", "warning")
                print(f"\n{'='*60}")
                print(f"✅ ANÁLISE CONCLUÍDA: {analysis_name}")
                print(f"📁 Resultados salvos em: {output_folder}")
                print(f"{'='*60}")
            except Exception as e:
                print(f"\n❌ ERRO DURANTE A ANÁLISE: {e}", "error")
                import traceback
                print(traceback.format_exc())
            finally:
                self.status_queue.put(("done", None))

    # --- Funções de Fila e Atualização da GUI ---

    def _process_queues(self):
        """Processa as filas de stdout e status (loop principal da GUI)"""
        try:
            while True:
                text = self.stdout_queue.get_nowait()
                if "❌" in text or "Erro" in text.capitalize() or "Falha" in text:
                    self._write_to_console(text, "error")
                elif "✅" in text or "Sucesso" in text or "concluída" in text.lower():
                    self._write_to_console(text, "success")
                elif "📊" in text or "📈" in text or "Iniciando" in text or "🚀" in text:
                    self._write_to_console(text, "info")
                elif "⚠️" in text or "Aviso" in text:
                    self._write_to_console(text, "warning")
                else:
                    self._write_to_console(text)
        except queue.Empty:
            pass
        
        try:
            command, data = self.status_queue.get_nowait()
            
            if command == "done":
                self.run_button.configure(state='normal')
                self.progress_bar.stop()
                self.progress_bar.set(1.0)
                self.progress_bar.pack_forget() # Esconder a barra
                self.status_label.configure(text="Análise concluída com sucesso!")
            
            elif command == "tool_done":
                self.run_optimize_btn.configure(state='normal')
                self.run_harmonize_btn.configure(state='normal')
                if "Processando" in self.status_label.cget("text"):
                    self.status_label.configure(text="Operação finalizada.")
            
            elif command == "optimization_complete":
                self.synth_output_seq.configure(state='normal')
                self.synth_output_seq.delete(1.0, tk.END)
                self.synth_output_seq.insert(1.0, data)
                self.synth_output_seq.configure(state='disabled')
                self.notebook.set("🛠️ Biologia Sintética")
                self.status_label.configure(text="Sequência processada com sucesso!")
                self.progress_bar.stop()
                self.progress_bar.set(1.0)
                self.progress_bar.pack_forget()
                
            elif command == "progress":
                self.progress_bar.configure(mode='determinate') # Mudar para modo determinado
                self.progress_bar.set(data / 100.0)
                self.status_label.configure(text=f"Processando... {data}%")
                
            elif command == "image_ready":
                image_path, title = data
                self._display_image(image_path, title)
                
            elif command == "message":
                self.status_label.configure(text=data)
                
        except queue.Empty:
            pass
        
        self.after(100, self._process_queues)

    def _write_to_console(self, text, tag=None):
        """Escreve texto no console da GUI de forma segura"""
        self.console_text.configure(state='normal')
        if tag:
            self.console_text.insert(tk.END, text, tag)
        else:
            self.console_text.insert(tk.END, text)
        self.console_text.see(tk.END)
        self.console_text.configure(state='disabled')
    
    def _display_image(self, image_path, title):
        """Exibe uma imagem no canvas"""
        try:
            self.fig.clear()
            img = Image.open(image_path)
            ax = self.fig.add_subplot(111)
            ax.imshow(img)
            ax.axis('off')
            ax.set_title(title, fontsize=14, fontweight='bold', color="#DCE4EE") # Cor do texto clara
            ax.set_facecolor("#2B2B2B")
            self.fig.patch.set_facecolor("#2B2B2B")
            self.fig.set_tight_layout(True)
            self.canvas.draw()
            
            self.current_images.append((image_path, title))
            self.image_status_label.configure(text=f"Exibindo: {title}")
            
            self.notebook.set("📈 Visualizador de Gráficos")
            
        except Exception as e:
            self._write_to_console(f"❌ Erro ao exibir imagem {image_path}: {e}\n", "error")

    # --- Funções Utilitárias ---
    
    def _get_selected_files(self):
        """Retorna a lista de caminhos completos dos arquivos selecionados (dos checkboxes)"""
        selected_files = []
        for filename, (cb, var) in self.file_checkboxes.items():
            if var.get() == "on":
                selected_files.append(self.file_paths_map[filename])
        return selected_files
    
    def _get_genetic_code_id(self):
        """Retorna o ID numérico da tabela genética selecionada"""
        try:
            return int(self.genetic_code_var.get().split(':')[0])
        except:
            return 1 # Fallback
            
    def _print_dataframe_limitado(self, df, titulo):
        """Imprime DataFrame no console com limite de linhas"""
        limite = 50
        print(f"\n--- {titulo} (Primeiras {min(len(df), limite)} linhas) ---")
        if len(df) > limite:
            print(df.head(limite).to_string())
            print(f"\n... e mais {len(df) - limite} linhas (veja o .csv completo na pasta de saída).")
        else:
            print(df.to_string())
        print("--------------------------------------------------\n")

    def _on_closing(self):
        """Evento de fechamento da janela"""
        # A importação do messagebox deve ser do tkinter padrão, não ctk
        from tkinter import messagebox
        if messagebox.askokcancel("Sair", "Deseja sair do Kodon-X?", parent=self): # <--- MUDANÇA DE NOME
            sys.stdout = self.stdout_original
            self.destroy()
