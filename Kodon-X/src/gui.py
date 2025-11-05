# gui.py
import os
import sys
import tkinter as tk
from tkinter import ttk, filedialog, scrolledtext, messagebox
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
# --- BLOCO DA INTERFACE GRÁFICA (FRONTEND)
# ######################################################################

class Redirector:
    """Redireciona stdout para a GUI"""
    def __init__(self, queue):
        self.queue = queue
    
    def write(self, string):
        self.queue.put(string)
    
    def flush(self):
        pass

class KodonE_GUI:
    """Interface gráfica profissional para análise de códons"""
    
    def __init__(self, root):
        self.root = root
        self.root.title("Kódon-E - Análise Visual de Genomas")
        self.root.geometry("1440x850") # Um pouco maior para a nova interface
        
        # --- Paleta de Cores Profissional (Clara) ---
        self.BG_COLOR = "#F0F0F0"
        self.LEFT_PANEL_BG = "#FDFDFD"
        self.FRAME_BG = "#FFFFFF"
        self.ACCENT_COLOR = "#0078D4"
        self.ACCENT_ACTIVE = "#005A9E"
        self.TEXT_COLOR = "#222222"
        self.HEADER_COLOR = "#004578" # Azul escuro para títulos
        self.SUCCESS_COLOR = "#388E3C"
        self.ERROR_COLOR = "#D32F2F"
        self.WARNING_COLOR = "#F57C00"
        self.INFO_COLOR = "#0288D1"
        
        # Configurar tema
        self.style = ttk.Style()
        self.style.theme_use('clam')
        
        # --- Configuração de Estilos Personalizados ---
        
        # Geral
        self.style.configure('.', font=('Calibri', 10), background=self.BG_COLOR, foreground=self.TEXT_COLOR)
        
        # Frames
        self.style.configure('TFrame', background=self.BG_COLOR)
        self.style.configure('Left.TFrame', background=self.LEFT_PANEL_BG)
        self.style.configure('Content.TFrame', background=self.FRAME_BG)
        
        # Labels
        self.style.configure('TLabel', background=self.BG_COLOR, foreground=self.TEXT_COLOR, font=('Calibri', 10))
        self.style.configure('Left.TLabel', background=self.LEFT_PANEL_BG, font=('Calibri', 10))
        self.style.configure('Header.TLabel', font=('Arial', 14, 'bold'), foreground=self.HEADER_COLOR, background=self.FRAME_BG)
        self.style.configure('Content.TLabel', background=self.FRAME_BG, foreground=self.TEXT_COLOR, font=('Calibri', 10))
        
        # Botões
        self.style.configure('TButton', font=('Arial', 10, 'bold'), background='#E1E1E1', relief='flat', padding=(10, 5))
        self.style.map('TButton', background=[('active', '#CCCCCC')])
        
        self.style.configure('Accent.TButton', font=('Arial', 12, 'bold'), background=self.ACCENT_COLOR, foreground='white')
        self.style.map('Accent.TButton', background=[('active', self.ACCENT_ACTIVE)])
        
        # LabelFrames (Painéis com título)
        self.style.configure('TLabelframe', background=self.FRAME_BG, relief='solid', borderwidth=1, bordercolor='#DDD')
        self.style.configure('TLabelframe.Label', font=('Arial', 12, 'bold'), foreground=self.HEADER_COLOR, background=self.FRAME_BG, padding=(10, 5))
        
        self.style.configure('Left.TLabelframe', background=self.LEFT_PANEL_BG, relief='flat')
        self.style.configure('Left.TLabelframe.Label', font=('Arial', 11, 'bold'), foreground=self.HEADER_COLOR, background=self.LEFT_PANEL_BG)
        
        # Notebook (Abas)
        self.style.configure('TNotebook', background=self.BG_COLOR, tabposition='n')
        self.style.configure('TNotebook.Tab', font=('Arial', 10, 'bold'), padding=(15, 8), foreground='#555555')
        self.style.map('TNotebook.Tab', 
                       foreground=[('selected', self.ACCENT_COLOR)],
                       background=[('selected', self.FRAME_BG)],
                       expand=[('selected', (1, 1, 1, 0))])
        
        self.style.configure('TProgressbar', thickness=25, background=self.ACCENT_COLOR, troughcolor='#E0E0E0')

        self.style.configure('Status.TFrame', background='#EAEAEA')
        self.style.configure('Status.TLabel', background='#EAEAEA', font=('Calibri', 10))
        
        # Entradas e Combobox
        self.style.configure('TEntry', fieldbackground='white', foreground=self.TEXT_COLOR, relief='solid')
        self.style.configure('TCombobox', fieldbackground='white', foreground=self.TEXT_COLOR, relief='solid')
        self.style.map('TCombobox', fieldbackground=[('readonly', 'white')])

        # --- Variáveis de Estado ---
        self.stdout_queue = queue.Queue()
        self.status_queue = queue.Queue()
        self.file_paths_map = {}
        self.current_images = []
        self.analysis_data = self._get_analysis_definitions()
        
        # Variáveis Tkinter
        self.input_folder_var = tk.StringVar()
        self.output_folder_var = tk.StringVar()
        self.genetic_code_var = tk.StringVar()
        self.filter_cds_var = tk.StringVar(value="ATG")
        
        # --- NOVO: Variável para Análise tAI ---
        self.super_kingdom_var = tk.StringVar(value="Bactéria")
        
        # Variáveis para Análise 19 (Expressão)
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
        self.root.protocol("WM_DELETE_WINDOW", self._on_closing)
        
        # Mensagem de boas-vindas
        self._write_to_console("Bem-vindo ao Kódon-E v1.0!\n", "info")
        self._write_to_console("Interface redesenhada para uma experiência mais limpa e moderna.\n")
        self._write_to_console("Por favor, selecione sua pasta de entrada e clique em 'Carregar Arquivos'.\n\n")

    def _get_analysis_definitions(self):
            """Retorna descrições detalhadas e requisitos das análises"""
            return {
                '1: Estatística e CDS': {
                    'id': '1',
                    'description': ('Realiza análise estatística básica dos genomas (Tamanho, GC, Contigs) '
                                'e uma análise detalhada dos CDS (Codon de Início, Tamanho). '
                                'Útil para caracterização geral e validação.'),
                    'files_required': '1+', # 1 ou mais
                    'function': None 
                },
                '2: Listagem de Genes': {
                    'id': '2',
                    'description': ('Lista todos os genes, produtos e locus tags presentes em um único '
                                'arquivo GenBank. Essencial para inventário genético.'),
                    'files_required': '1', # Exatamente 1
                    'function': listar_genes_do_arquivo
                },
                '3: RSCU Heatmap Individual': {
                    'id': '3',
                    'description': ('Calcula e visualiza o RSCU (Relative Synonymous Codon Usage) para um '
                                'único genoma. Inclui métricas: ENC, GC3 e CAI.'),
                    'files_required': '1', # Exatamente 1
                    'function': gerar_rscu_heatmap_e_tabela
                },
                '4: RSCU Comparativo': {
                    'id': '4',
                    'description': ('Compara padrões de RSCU entre múltiplos genomas usando '
                                'Clustermap (Agrupamento) e PCA (Análise de Componentes Principais).'),
                    'files_required': '2+', # 2 ou mais
                    'function': analise_rscu_comparativa
                },
                '5: Correlação RSCU (2 Genomas)': {
                    'id': '5',
                    'description': ('Calcula a correlação de Pearson (R) entre os padrões de RSCU de '
                                'exatamente dois genomas. Indica similaridade no viés.'),
                    'files_required': '2', # Exatamente 2
                    'function': analise_correlacao_rscu
                },
                '6: Histogramas RSCU Comparativos': {
                    'id': '6',
                    'description': ('Gera visualizações detalhadas do uso de códons, incluindo Box Plots '
                                'por aminoácido e gráficos de barras empilhadas por espécie.'),
                    'files_required': '2+', # 2 ou mais
                    'function': gerar_histogramas_rscu
                },
                '7: Análise ENC vs GC3 (Wright Plot)': {
                    'id': '7',
                    'description': ('Plota o "Effective Number of Codons" (ENC) contra o conteúdo de GC '
                                'na 3ª posição (GC3) para determinar pressões evolutivas (Mutação vs Seleção).'),
                    'files_required': '1+', # 1 ou mais
                    'function': analise_enc_gc3
                },
                '8: Composição Genômica': {
                    'id': '8',
                    'description': ('Análise detalhada da composição nucleotídica (A, T, G, C), conteúdo de GC '
                                'total e tamanho do genoma para múltiplos arquivos.'),
                    'files_required': '1+', # 1 ou mais
                    'function': analise_composicao_genomica
                },
                '9: Códons Ótimos, Raros e CAI': {
                    'id': '9',
                    'description': ('Identifica códons preferenciais (RSCU > 1.2) e raros (RSCU < 0.8) e '
                                'compara o "Codon Adaptation Index" (CAI) entre espécies.'),
                    'files_required': '1+', # 1 ou mais
                    'function': analise_codons_otimos_raros
                },
                '10: Viés de Pares de Códons (CPB)': {
                    'id': '10',
                    'description': ('Analisa a frequência de pares de códons adjacentes (ex: '
                                    'ATG-CGC). Gera uma matriz 61x61 (heatmap) mostrando pares '
                                    'favorecidos (positivos) ou evitados (negativos) na tradução.'),
                    'files_required': '1+', # 1 ou mais
                    'function': analise_codon_pair_bias
                },
                '11: Análise Físico-Química (GRAVY & Aromo)': {
                    'id': '11',
                    'description': ('Calcula o escore GRAVY (hidropaticidade média) e Aromo '
                                    '(aromaticidade média) para todos os genes. '
                                    'Gera boxplots comparando as distribuições entre espécies.'),
                    'files_required': '1+', # 1 ou mais
                    'function': analise_gravy_aromo
                },
                '12: Gráfico de Neutralidade (GC12 vs GC3)': {
                    'id': '12',
                    'description': ('Plota o conteúdo de GC das posições 1+2 (GC12) contra o GC3. '
                                    'A inclinação (slope) da regressão indica o balanço entre '
                                    'pressão mutacional (slope ~1) e seleção (slope ~0).'),
                    'files_required': '2+', # 2 ou mais para regressão
                    'function': analise_neutrality_plot
                },
                '13: Composição de Dinucleotídeos': {
                    'id': '13',
                    'description': ('Calcula e compara a frequência de todos os 16 dinucleotídeos '
                                    '(AA, AT, GC, etc.) no genoma completo. '
                                    'Gera um heatmap comparativo.'),
                    'files_required': '1+', # 1 ou mais
                    'function': analise_dinucleotide_composition
                },
                '14: Gráfico de Paridade PR2 (A3/T3 vs G3/C3)': {
                    'id': '14',
                    'description': ('(Feature c) Plota A3/(A3+T3) vs G3/(G3+C3) para cada gene. '
                                    'Analisa viés mutacional e de seleção específico da fita. '
                                    'Desvios do centro (0.5, 0.5) indicam viés.'),
                    'files_required': '1+', # 1 ou mais
                    'function': analise_pr2_plot
                },
                
                # --- ANÁLISE 15 ATUALIZADA ---
                '15: Índice de Adaptação ao tRNA (tAI)': {
                    'id': '15',
                    'description': ('(Feature a - ATUALIZADA) Calcula o tAI (tRNA Adaptation Index) '
                                    'usando a metodologia de **pareamento wobble (dos Reis et al.)**. '
                                    'Este é o padrão-ouro para medir a adaptação translacional '
                                    'com base na contagem de genes de tRNA do genoma.'),
                    'files_required': '1+', # 1 ou mais
                    'function': analise_tAI 
                },
                '16: Análise de Motifs Upstream': {
                    'id': '16',
                    'description': ('(Feature b) Busca por k-mers (motifs) sobrerrepresentados nas regiões '
                                    'N-pb upstream (antes do início) de cada gene. '
                                    'Útil para encontrar sítios de ligação (promotores).'),
                    'files_required': '1+', # 1 ou mais
                    'function': analise_motifs_upstream
                },
                '17: Análise de MFE (Estrutura 5\')': {
                    'id': '17',
                    'description': ('Calcula a Energia Livre Mínima (MFE) da região 5\' (início) de '
                                    'todos os CDS. Valores muito negativos indicam estruturas '
                                    'secundárias (grampos) que podem bloquear a iniciação da tradução.'),
                    'files_required': '1+', # 1 ou mais
                    'function': analise_mfe_iniciacao
                },
                
                # --- INÍCIO DAS NOVAS ANÁLISES (18 e 19) ---
                '18: Comparação de Grupos de Genes': {
                    'id': '18',
                    'description': ('(Feature 5) Compara métricas de CUB (ENC, CAI, GC3) entre '
                                    'dois conjuntos de genes (Grupo 1 vs Grupo 2). '
                                    'Executa testes estatísticos (Mann-Whitney U, Qui-Quadrado) '
                                    'e gera boxplots comparativos.'),
                    'files_required': '1+', # 1 ou mais
                    'function': analise_comparativa_dois_grupos
                },
                '19: Correlação com Expressão (RNA-Seq)': {
                    'id': '19',
                    'description': ('(Feature 7) Correlaciona métricas de CUB (ENC, CAI) com dados '
                                    'de expressão (ex: TPM, RPKM) de um arquivo CSV/TSV. '
                                    'Gera gráficos de dispersão com correlação de Spearman.'),
                    'files_required': '1+', # 1 ou mais
                    'function': analise_correlacao_expressao
                }
                # --- FIM DAS NOVAS ANÁLISES ---
            }

    # --- Criação de Widgets ---

    def _create_widgets(self):
        """Cria todos os widgets da interface"""
        
        # Menu principal
        self._create_menu()
        
        # Container principal (Divide em Esquerda e Direita)
        main_paned_window = ttk.PanedWindow(self.root, orient=tk.HORIZONTAL)
        main_paned_window.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # --- Painel Esquerdo (Configuração) ---
        left_frame = self._create_left_panel()
        main_paned_window.add(left_frame, weight=1)
        
        # --- Painel Principal (Análise e Resultados) ---
        main_frame = self._create_main_panel()
        main_paned_window.add(main_frame, weight=3)
        
        # --- Barra de Status ---
        self._create_status_bar()

    def _create_menu(self):
        """Cria o menu principal"""
        menubar = tk.Menu(self.root, font=('Calibri', 10))
        self.root.config(menu=menubar)
        
        # Menu Arquivo
        file_menu = tk.Menu(menubar, tearoff=0, font=('Calibri', 10))
        menubar.add_cascade(label="Arquivo", menu=file_menu)
        file_menu.add_command(label="Carregar Pasta...", command=self._on_browse_input, accelerator="Ctrl+O")
        file_menu.add_command(label="Definir Saída...", command=self._on_browse_output, accelerator="Ctrl+S")
        file_menu.add_separator()
        file_menu.add_command(label="Sair", command=self._on_closing, accelerator="Ctrl+Q")
        
        # Menu Ajuda
        help_menu = tk.Menu(menubar, tearoff=0, font=('Calibri', 10))
        menubar.add_cascade(label="Ajuda", menu=help_menu)
        help_menu.add_command(label="Sobre Kódon-E", command=self._show_about)

    def _create_left_panel(self):
        """Cria o painel de configuração da esquerda"""
        left_panel = ttk.Frame(self.root, style='Left.TFrame', width=450)
        left_panel.pack(fill=tk.Y, side=tk.LEFT, padx=(0, 5))
        
        # --- Passo 1: Arquivos de Entrada ---
        input_frame = ttk.LabelFrame(left_panel, text=" 📁 Passo 1: Arquivos de Entrada", style='Left.TLabelframe', padding="15")
        input_frame.pack(fill=tk.X, padx=10, pady=10)
        
        ttk.Label(input_frame, text="Pasta de Entrada:", style='Left.TLabel').grid(row=0, column=0, sticky="w", padx=5, pady=(0, 5))
        input_entry = ttk.Entry(input_frame, textvariable=self.input_folder_var, width=35)
        input_entry.grid(row=1, column=0, sticky="ew", padx=5, pady=2)
        
        browse_input_btn = ttk.Button(input_frame, text="Procurar...", command=self._on_browse_input)
        browse_input_btn.grid(row=1, column=1, padx=5, pady=2)
        
        load_files_btn = ttk.Button(input_frame, text="Carregar Arquivos", command=self._on_load_files)
        load_files_btn.grid(row=2, column=0, columnspan=2, sticky="ew", padx=5, pady=10)
        
        # Listbox com scrollbar
        list_container = ttk.Frame(input_frame, style='Left.TFrame')
        list_container.grid(row=3, column=0, columnspan=2, sticky="nsew", padx=5, pady=5)
        input_frame.grid_rowconfigure(3, weight=1)
        input_frame.grid_columnconfigure(0, weight=1)
        
        self.file_listbox = tk.Listbox(list_container, selectmode=tk.MULTIPLE, height=15, font=('Courier', 10),
                                        bg='white', fg=self.TEXT_COLOR, relief='solid', bd=1,
                                        selectbackground=self.ACCENT_COLOR, selectforeground='white',
                                        highlightthickness=0)
        list_scrollbar = ttk.Scrollbar(list_container, orient=tk.VERTICAL, command=self.file_listbox.yview)
        self.file_listbox.config(yscrollcommand=list_scrollbar.set)
        
        self.file_listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        list_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
        # Botões de seleção
        file_buttons = ttk.Frame(input_frame, style='Left.TFrame')
        file_buttons.grid(row=4, column=0, columnspan=2, sticky="ew", padx=5, pady=(10, 5))
        
        ttk.Button(file_buttons, text="Todos", command=self._on_select_all).pack(side=tk.LEFT, padx=2, expand=True, fill=tk.X)
        ttk.Button(file_buttons, text="Nenhum", command=self._on_select_none).pack(side=tk.LEFT, padx=2, expand=True, fill=tk.X)
        self.file_status_label = ttk.Label(file_buttons, text="0 arquivos", style='Left.TLabel', anchor="e")
        self.file_status_label.pack(side=tk.RIGHT, padx=5)

        # --- Passo 2: Pasta de Saída ---
        output_frame = ttk.LabelFrame(left_panel, text=" 💾 Passo 2: Pasta de Saída", style='Left.TLabelframe', padding="15")
        output_frame.pack(fill=tk.X, padx=10, pady=10)
        
        output_entry = ttk.Entry(output_frame, textvariable=self.output_folder_var, width=35)
        output_entry.grid(row=0, column=0, sticky="ew", padx=5, pady=2)
        output_frame.grid_columnconfigure(0, weight=1)
        
        browse_output_btn = ttk.Button(output_frame, text="Procurar...", command=self._on_browse_output)
        browse_output_btn.grid(row=0, column=1, padx=5, pady=2)

        # --- Passo 3: Tabela Genética ---
        genetic_frame = ttk.LabelFrame(left_panel, text=" ⚙️ Passo 3: Tabela Genética", style='Left.TLabelframe', padding="15")
        genetic_frame.pack(fill=tk.X, padx=10, pady=10)
        
        genetic_combo = ttk.Combobox(genetic_frame, textvariable=self.genetic_code_var, state="readonly", width=35,
                                     font=('Calibri', 10))
        genetic_combo['values'] = [
            "1: Padrão (Universal)", 
            "2: Mitocondrial de Vertebrados",
            "4: Mitocondrial de Mofos/Protozoários", 
            "11: Plasto de Bactérias/Plantas"
        ]
        genetic_combo.set("1: Padrão (Universal)")
        genetic_combo.pack(fill=tk.X, padx=5, pady=5)

        return left_panel
    
    def _create_main_panel(self):
        """Cria o painel principal com abas (Notebook)"""
        main_panel = ttk.Frame(self.root, style='TFrame')
        main_panel.pack(fill=tk.BOTH, expand=True, side=tk.RIGHT, padx=(5, 0))
        
        self.notebook = ttk.Notebook(main_panel)
        self.notebook.pack(fill=tk.BOTH, expand=True)
        
        analysis_tab = ttk.Frame(self.notebook, style='Content.TFrame', padding="10")
        synth_bio_tab = ttk.Frame(self.notebook, style='Content.TFrame', padding="10")
        viz_tab = ttk.Frame(self.notebook, style='Content.TFrame', padding="10")
        console_tab = ttk.Frame(self.notebook, style='Content.TFrame', padding="10")
        help_tab = ttk.Frame(self.notebook, style='Content.TFrame', padding="10")

        # --- Aba 1: Executar Análise ---
        self._create_analysis_tab(analysis_tab)
        self.notebook.add(analysis_tab, text=" ▶️ Executar Análise ")
        
        # --- Aba 2: Biologia Sintética ---
        self._create_synth_bio_tab(synth_bio_tab)
        self.notebook.add(synth_bio_tab, text=" 🛠️ Biologia Sintética ")

        # --- Aba 3: Visualizador de Gráficos ---
        self._create_visualization_tab(viz_tab)
        self.notebook.add(viz_tab, text=" 📈 Visualizador de Gráficos ")

        # --- Aba 4: Console de Saída ---
        self._create_console_tab(console_tab)
        self.notebook.add(console_tab, text=" 🖥️ Console de Saída ")
        
        # --- Aba 5: Ajuda ---
        self._create_help_tab(help_tab)
        self.notebook.add(help_tab, text=" ❓ Ajuda ")
        
        return main_panel

    def _create_analysis_tab(self, tab):
        """Cria a aba de seleção e execução de análise"""
        
        # Topo: Seleção da Análise
        ttk.Label(tab, text="Selecione a Análise:", style='Header.TLabel').pack(anchor="w", padx=5, pady=5)
        
        self.analysis_combo = ttk.Combobox(tab, state="readonly", font=('Calibri', 11),
                                           values=list(self.analysis_data.keys()))
        self.analysis_combo.pack(fill=tk.X, padx=5, pady=(5, 10))
        
        # Meio: Descrição e Opções
        desc_frame = ttk.LabelFrame(tab, text="Descrição da Análise", padding="10")
        desc_frame.pack(fill=tk.X, padx=5, pady=10)
        
        self.analysis_desc_text = tk.Text(desc_frame, height=6, wrap=tk.WORD, font=('Calibri', 10), 
                                          bg='#FDFDFD', relief="flat", fg=self.TEXT_COLOR,
                                          highlightthickness=1, highlightbackground='#DDD')
        self.analysis_desc_text.pack(fill=tk.X, expand=True, padx=5, pady=5)
        self.analysis_desc_text.insert(1.0, "Selecione uma análise acima para ver sua descrição e requisitos.")
        self.analysis_desc_text.config(state="disabled")

        self.analysis_options_frame = ttk.Frame(tab, style='Content.TFrame')
        self.analysis_options_frame.pack(fill=tk.X, padx=5, pady=5)
        
        # --- Filtro de Subconjunto de Genes (Feature d) ---
        gene_filter_frame = ttk.LabelFrame(tab, text="Filtro de Genes Global (Opcional)", padding="10")
        gene_filter_frame.pack(fill=tk.X, padx=5, pady=(15, 5))
        
        ttk.Label(gene_filter_frame, 
                  text="Cole uma lista de locus_tags ou nomes de genes (um por linha).\nSe preenchido, TODAS as análises usarão apenas este subconjunto de genes.", 
                  style='Content.TLabel').pack(anchor="w", padx=5)
        
        self.gene_filter_text = scrolledtext.ScrolledText(gene_filter_frame, height=5, wrap=tk.WORD, 
                                                          font=('Courier', 10), relief="solid", borderwidth=1,
                                                          bg="#FDFDFD", fg=self.TEXT_COLOR,
                                                          highlightthickness=0)
        self.gene_filter_text.pack(fill=tk.X, expand=True, padx=5, pady=5)
        
        # Baixo: Botão de Execução
        self.run_button = ttk.Button(tab, text="🚀 Executar Análise Selecionada", 
                                     style="Accent.TButton", command=self._on_run_analysis)
        self.run_button.pack(fill=tk.X, padx=5, pady=20, ipady=10)
        
        return tab

    def _create_synth_bio_tab(self, tab):
        """Cria a aba de ferramentas de Biologia Sintética"""
        
        # --- Variáveis de estado para esta aba ---
        self.synth_host_var = tk.StringVar()
        self.synth_input_seq = None # Widget de texto
        self.synth_output_seq = None # Widget de texto
        
        # Frame principal dividido (Esquerda: Input/Host, Direita: Output)
        main_pane = ttk.PanedWindow(tab, orient=tk.HORIZONTAL)
        main_pane.pack(fill=tk.BOTH, expand=True)

        # --- Painel Esquerdo (Configuração) ---
        left_frame = ttk.Frame(main_pane, style='Content.TFrame', padding="5")
        main_pane.add(left_frame, weight=1)
        
        # 1. Seleção do Hospedeiro
        host_frame = ttk.LabelFrame(left_frame, text="1. Selecione o Hospedeiro", padding="10")
        host_frame.pack(fill=tk.X, padx=5, pady=5)
        
        ttk.Label(host_frame, text="Selecione um genoma carregado para ser o 'Hospedeiro':", style='Content.TLabel').pack(anchor="w")
        
        self.synth_host_combo = ttk.Combobox(host_frame, textvariable=self.synth_host_var, state="readonly", font=('Calibri', 11))
        self.synth_host_combo.pack(fill=tk.X, padx=5, pady=5)
        
        # 2. Sequência de Entrada
        input_frame = ttk.LabelFrame(left_frame, text="2. Cole a Sequência de DNA (CDS)", padding="10")
        input_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        self.synth_input_seq = scrolledtext.ScrolledText(input_frame, height=15, wrap=tk.WORD, 
                                                          font=('Courier', 10), relief="solid", borderwidth=1,
                                                          bg="#FDFDFD", fg=self.TEXT_COLOR,
                                                          highlightthickness=0)
        self.synth_input_seq.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # 3. Botões de Ação
        button_frame = ttk.Frame(left_frame, style='Content.TFrame')
        button_frame.pack(fill=tk.X, padx=5, pady=10)
        
        self.run_optimize_btn = ttk.Button(button_frame, text="Otimizar\n(Maximizar RSCU)", 
                                           style="TButton", command=self._on_run_optimization)
        self.run_optimize_btn.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=5, ipady=5)
        
        self.run_harmonize_btn = ttk.Button(button_frame, text="Harmonizar\n(Manter Ranking)", 
                                            style="TButton", command=self._on_run_harmonization)
        self.run_harmonize_btn.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=5, ipady=5)


        # --- Painel Direito (Resultados) ---
        right_frame = ttk.Frame(main_pane, style='Content.TFrame', padding="5")
        main_pane.add(right_frame, weight=1)
        
        output_frame = ttk.LabelFrame(right_frame, text="Sequência Resultante", padding="10")
        output_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        self.synth_output_seq = scrolledtext.ScrolledText(output_frame, height=20, wrap=tk.WORD, 
                                                          font=('Courier', 10), relief="solid", borderwidth=1,
                                                          bg="#F0F0F0", fg=self.TEXT_COLOR,
                                                          highlightthickness=0, state='disabled')
        self.synth_output_seq.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Conectar o combobox ao evento de carregamento de arquivos
        self.root.bind("<<FilesLoaded>>", self._update_synth_host_list)
        return tab

    def _create_visualization_tab(self, tab):
        """Cria a aba de visualização de gráficos"""
        
        self.viz_canvas_frame = ttk.Frame(tab, style='Content.TFrame')
        self.viz_canvas_frame.pack(fill=tk.BOTH, expand=True)
        
        # Figura para visualização
        self.fig = Figure(figsize=(10, 8), dpi=100)
        self.fig.patch.set_facecolor(self.FRAME_BG) # Fundo da figura
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.viz_canvas_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        toolbar_frame = ttk.Frame(tab, style='Content.TFrame')
        toolbar_frame.pack(fill=tk.X, pady=(5,0))
        self.toolbar = NavigationToolbar2Tk(self.canvas, toolbar_frame)
        self.toolbar.config(background=self.FRAME_BG) 
        for button in self.toolbar.winfo_children():
            button.config(background=self.FRAME_BG)
        self.toolbar.update()
        
        controls_frame = ttk.Frame(tab, style='Content.TFrame')
        controls_frame.pack(fill=tk.X, pady=5)
        
        self.image_status_label = ttk.Label(controls_frame, text="Nenhum gráfico carregado.", style='Content.TLabel')
        self.image_status_label.pack(side=tk.LEFT, padx=5)
        
        return tab

    def _create_console_tab(self, tab):
        """Cria a aba do console de saída"""
        
        self.console_text = scrolledtext.ScrolledText(tab, height=15, wrap=tk.WORD, 
                                                      font=('Courier', 10), relief="solid", borderwidth=1,
                                                      bg="#FDFDFD", fg=self.TEXT_COLOR,
                                                      highlightthickness=0)
        
        self.console_text.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        self.console_text.configure(state='disabled')
        
        # Configurar tags para coloração
        self.console_text.tag_configure("error", foreground=self.ERROR_COLOR, font=('Courier', 10, 'bold'))
        self.console_text.tag_configure("success", foreground=self.SUCCESS_COLOR, font=('Courier', 10, 'bold'))
        self.console_text.tag_configure("info", foreground=self.INFO_COLOR, font=('Courier', 10, 'bold'))
        self.console_text.tag_configure("warning", foreground=self.WARNING_COLOR, font=('Courier', 10, 'bold'))
        
        return tab

    def _create_help_tab(self, tab):
        """Cria a aba de ajuda com texto detalhado (NOVO)"""
        
        help_text_widget = scrolledtext.ScrolledText(tab, wrap=tk.WORD, font=('Calibri', 11), 
                                                 bg="#FDFDFD", relief='flat', padx=10)
        help_text_widget.pack(fill=tk.BOTH, expand=True)
        
        # --- Tags de Estilo para o Texto de Ajuda ---
        help_text_widget.tag_configure('h1', font=('Arial', 16, 'bold'), foreground=self.HEADER_COLOR, spacing3=10)
        help_text_widget.tag_configure('h2', font=('Arial', 14, 'bold'), foreground=self.TEXT_COLOR, spacing3=8, spacing1=10)
        help_text_widget.tag_configure('h3', font=('Arial', 12, 'bold'), foreground=self.ACCENT_COLOR, spacing3=5, spacing1=8)
        help_text_widget.tag_configure('bold', font=('Calibri', 11, 'bold'))
        help_text_widget.tag_configure('code', font=('Courier', 10), background='#EAEAEA')
        
        # --- Conteúdo do Manual ---
        
        help_text_widget.insert(tk.END, "KÓDON-E v1.0 - GUIA DO USUÁRIO\n", 'h1')
        
        help_text_widget.insert(tk.END, "Bem-vindo à interface de análise Kódon-E. Siga este guia para um fluxo de trabalho eficiente.\n", 'normal')
        
        help_text_widget.insert(tk.END, "\nFLUXO DE TRABALHO RÁPIDO\n", 'h2')
        help_text_widget.insert(tk.END, "1. ", 'h3')
        help_text_widget.insert(tk.END, "Carregar Arquivos: ", 'bold')
        help_text_widget.insert(tk.END, "Use o painel 'Passo 1' para 'Procurar...' sua pasta com arquivos .gbk ou .gb. Clique em 'Carregar Arquivos'.\n", 'normal')
        
        help_text_widget.insert(tk.END, "2. ", 'h3')
        help_text_widget.insert(tk.END, "Selecionar Arquivos: ", 'bold')
        help_text_widget.insert(tk.END, "Na lista, selecione os arquivos que deseja analisar (Ctrl+Click para múltiplos ou use 'Todos').\n", 'normal')
        
        help_text_widget.insert(tk.END, "3. ", 'h3')
        help_text_widget.insert(tk.END, "Definir Saída: ", 'bold')
        help_text_widget.insert(tk.END, "Use o 'Passo 2' para escolher onde salvará os gráficos e tabelas .csv.\n", 'normal')
        
        help_text_widget.insert(tk.END, "4. ", 'h3')
        help_text_widget.insert(tk.END, "Selecionar Análise: ", 'bold')
        help_text_widget.insert(tk.END, "Vá para a aba '▶️ Executar Análise' e escolha uma análise no menu. A descrição e os requisitos de arquivos aparecerão abaixo.\n", 'normal')
        
        help_text_widget.insert(tk.END, "5. ", 'h3')
        help_text_widget.insert(tk.END, "Executar: ", 'bold')
        help_text_widget.insert(tk.END, "Clique no botão '🚀 Executar Análise Selecionada'.\n", 'normal')
        
        help_text_widget.insert(tk.END, "6. ", 'h3')
        help_text_widget.insert(tk.END, "Ver Resultados: ", 'bold')
        help_text_widget.insert(tk.END, "Acompanhe o progresso no '🖥️ Console' e veja os gráficos gerados no '📈 Visualizador'.\n\n", 'normal')
        
        
        help_text_widget.insert(tk.END, "\nGUIA DE ANÁLISE (FUNÇÃO POR FUNÇÃO)\n", 'h2')
        
        # --- Análise 1 ---
        help_text_widget.insert(tk.END, "1: Estatística e CDS\n", 'h3')
        help_text_widget.insert(tk.END, "O que faz: ", 'bold')
        help_text_widget.insert(tk.END, "Gera duas tabelas. A primeira ('estatisticas_genoma.csv') mostra dados gerais: número de contigs, tamanho total do genoma (pb) e conteúdo GC total (%). A segunda ('analise_cds.csv') lista todas as sequências codificantes (CDS), seu códon de início, e se ele bate com o filtro (padrão 'ATG').\n", 'normal')
        help_text_widget.insert(tk.END, "Para que serve: ", 'bold')
        help_text_widget.insert(tk.END, "Excelente para uma visão geral e controle de qualidade. Permite validar rapidamente a montagem (tamanho, GC) e verificar a anotação dos genes (ex: quantos genes começam com ATG vs GTG).\n", 'normal')
        help_text_widget.insert(tk.END, "Como usar: ", 'bold')
        help_text_widget.insert(tk.END, "Selecione 1 ou mais arquivos. Você pode alterar o códon de início no filtro que aparece abaixo da descrição.\n", 'normal')

        # --- Análise 2 ---
        help_text_widget.insert(tk.END, "2: Listagem de Genes\n", 'h3')
        help_text_widget.insert(tk.END, "O que faz: ", 'bold')
        help_text_widget.insert(tk.END, "Varre o arquivo GenBank e extrai cada gene anotado (CDS, tRNA, rRNA), listando seu Locus Tag, Nome do Gene (se houver) e Produto (descrição da proteína).\n", 'normal')
        help_text_widget.insert(tk.END, "Para que serve: ", 'bold')
        help_text_widget.insert(tk.END, "Ideal para criar um inventário completo de todos os genes em um genoma, facilitando buscas e comparações.\n", 'normal')
        help_text_widget.insert(tk.END, "Como usar: ", 'bold')
        help_text_widget.insert(tk.END, "Selecione EXATAMENTE 1 arquivo.\n", 'normal')

        # --- Análise 3 ---
        help_text_widget.insert(tk.END, "3: RSCU Heatmap Individual\n", 'h3')
        help_text_widget.insert(tk.END, "O que faz: ", 'bold')
        help_text_widget.insert(tk.END, "Calcula o 'Relative Synonymous Codon Usage' (RSCU) para todos os códons. O RSCU mede o viés de uso (se um códon é usado mais ou menos do que o esperado). Gera um heatmap visual (imagem) e uma tabela de contagem ('contagem_codons_...csv'). Também calcula métricas avançadas: ENC, GC3 e CAI.\n", 'normal')
        help_text_widget.insert(tk.END, "Para que serve: ", 'bold')
        help_text_widget.insert(tk.END, "Permite uma análise profunda do viés de uso de códons de um único organismo. O heatmap mostra visualmente quais códons são preferidos (RSCU > 1.0) ou evitados (RSCU < 1.0).\n", 'normal')
        help_text_widget.insert(tk.END, "Como usar: ", 'bold')
        help_text_widget.insert(tk.END, "Selecione EXATAMENTE 1 arquivo.\n", 'normal')

        # --- Análise 4 ---
        help_text_widget.insert(tk.END, "4: RSCU Comparativo\n", 'h3')
        help_text_widget.insert(tk.END, "O que faz: ", 'bold')
        help_text_widget.insert(tk.END, "Calcula o RSCU para múltiplos genomas e os compara usando duas técnicas estatísticas: Agrupamento Hierárquico (Clustermap) e Análise de Componentes Principais (PCA).\n", 'normal')
        help_text_widget.insert(tk.END, "Para que serve: ", 'bold')
        help_text_widget.insert(tk.END, "Determina a 'distância' evolutiva ou adaptativa entre diferentes espécies com base em seu viés de códons. O Clustermap agrupa espécies com padrões similares, e o PCA mostra essa relação em um gráfico 2D.\n", 'normal')
        help_text_widget.insert(tk.END, "Como usar: ", 'bold')
        help_text_widget.insert(tk.END, "Selecione 2 ou mais arquivos.\n", 'normal')
        
        # --- Análise 5 ---
        help_text_widget.insert(tk.END, "5: Correlação RSCU (2 Genomas)\n", 'h3')
        help_text_widget.insert(tk.END, "O que faz: ", 'bold')
        help_text_widget.insert(tk.END, "Compara o RSCU de cada um dos 64 códons entre dois genomas e calcula a Correlação de Pearson (R). Gera um gráfico de dispersão.\n", 'normal')
        help_text_widget.insert(tk.END, "Para que serve: ", 'bold')
        help_text_widget.insert(tk.END, "Mede diretamente o quão similar é o viés de uso de códons entre duas espécies. Um R próximo de 1.0 indica um padrão de uso quase idêntico.\n", 'normal')
        help_text_widget.insert(tk.END, "Como usar: ", 'bold')
        help_text_widget.insert(tk.END, "Selecione EXATAMENTE 2 arquivos.\n", 'normal')
        
        # --- Análise 6 ---
        help_text_widget.insert(tk.END, "6: Histogramas RSCU Comparativos\n", 'h3')
        help_text_widget.insert(tk.END, "O que faz: ", 'bold')
        help_text_widget.insert(tk.END, "Gera dois tipos de gráficos: (1) Um Box Plot comparando a distribuição de RSCU para cada aminoácido entre todas as espécies; e (2) um gráfico de barras empilhadas para cada espécie, mostrando a contribuição de cada códon sinônimo.\n", 'normal')
        help_text_widget.insert(tk.END, "Para que serve: ", 'bold')
        help_text_widget.insert(tk.END, "Uma visualização detalhada e 'zoom-in' do RSCU. Ótimo para ver quais aminoácidos sofrem mais viés (ex: Leucina) e como cada espécie 'resolve' esse viés.\n", 'normal')
        help_text_widget.insert(tk.END, "Como usar: ", 'bold')
        help_text_widget.insert(tk.END, "Selecione 2 ou mais arquivos.\n", 'normal')
        
        # --- Análise 7 ---
        help_text_widget.insert(tk.END, "7: Análise ENC vs GC3 (Wright Plot)\n", 'h3')
        help_text_widget.insert(tk.END, "O que faz: ", 'bold')
        help_text_widget.insert(tk.END, "Plota o 'Effective Number of Codons' (ENC) contra o conteúdo de GC na 3ª posição (GC3). O ENC mede o quão extremo é o viés (61 = sem viés, 20 = viés extremo). O gráfico inclui a 'Curva Esperada' de Wright.\n", 'normal')
        help_text_widget.insert(tk.END, "Para que serve: ", 'bold')
        help_text_widget.insert(tk.END, "É a análise clássica para determinar a principal força evolutiva no genoma. Pontos na curva esperada sugerem que o viés é causado por pressão mutacional (ex: GC-bias). Pontos abaixo da curva sugerem seleção natural (tradução eficiente).\n", 'normal')
        help_text_widget.insert(tk.END, "Como usar: ", 'bold')
        help_text_widget.insert(tk.END, "Selecione 1 ou mais arquivos.\n", 'normal')
        
        # --- Análise 8 ---
        help_text_widget.insert(tk.END, "8: Composição Genômica\n", 'h3')
        help_text_widget.insert(tk.END, "O que faz: ", 'bold')
        help_text_widget.insert(tk.END, "Calcula a composição de A, T, G, C (%), o GC total (%) e o tamanho total do genoma (pb) para todos os arquivos. Gera gráficos comparativos.\n", 'normal')
        help_text_widget.insert(tk.END, "Para que serve: ", 'bold')
        help_text_widget.insert(tk.END, "Semelhante à Análise 1, mas focada na comparação visual entre múltiplos genomas. Útil para identificar outliers ou grupos com base em composição.\n", 'normal')
        help_text_widget.insert(tk.END, "Como usar: ", 'bold')
        help_text_widget.insert(tk.END, "Selecione 1 ou mais arquivos.\n", 'normal')
        
        # --- Análise 9 ---
        help_text_widget.insert(tk.END, "9: Códons Ótimos, Raros e CAI\n", 'h3')
        help_text_widget.insert(tk.END, "O que faz: ", 'bold')
        help_text_widget.insert(tk.END, "Identifica códons 'ótimos' (RSCU > 1.2) e 'raros' (RSCU < 0.8). Também calcula o 'Codon Adaptation Index' (CAI), que mede o quão 'otimizado' um genoma é para alta expressão. Gera gráficos comparativos.\n", 'normal')
        help_text_widget.insert(tk.END, "Para que serve: ", 'bold')
        help_text_widget.insert(tk.END, "Focado em engenharia genética e expressão. O CAI é uma métrica chave para prever níveis de expressão gênica. A lista de códons ótimos/raros é fundamental para otimizar sequências gênicas.\n", 'normal')
        help_text_widget.insert(tk.END, "Como usar: ", 'bold')
        help_text_widget.insert(tk.END, "Selecione 1 ou mais arquivos.\n", 'normal')

        # --- Análise 10 ---
        help_text_widget.insert(tk.END, "10: Viés de Pares de Códons (CPB)\n", 'h3')
        help_text_widget.insert(tk.END, "O que faz: ", 'bold')
        help_text_widget.insert(tk.END, "Calcula a frequência observada de pares de códons (ex: ATG seguido de CGC) em relação à frequência esperada. Gera um 'heatmap' 61x61 que mostra quais pares são preferidos (valor > 0) ou evitados (valor < 0).\n", 'normal')
        help_text_widget.insert(tk.END, "Para que serve: ", 'bold')
        help_text_widget.insert(tk.END, "A velocidade e eficiência da tradução (produção de proteína) não depende apenas de códons individuais, mas também de pares. Esta análise revela 'gargalos' traducionais ou otimizações.\n", 'normal')
        help_text_widget.insert(tk.END, "Como usar: ", 'bold')
        help_text_widget.insert(tk.END, "Selecione 1 ou mais arquivos. Gerará um heatmap e CSV por arquivo.\n", 'normal')
        
        # --- Análise 11 ---
        help_text_widget.insert(tk.END, "11: Análise Físico-Química (GRAVY & Aromo)\n", 'h3')
        help_text_widget.insert(tk.END, "O que faz: ", 'bold')
        help_text_widget.insert(tk.END, "Traduz todos os genes e calcula duas métricas para cada proteína: (1) GRAVY (hidropaticidade média), que indica se a proteína é hidrofílica (negativa) ou hidrofóbica (positiva); e (2) Aromo (aromaticidade média).\n", 'normal')
        help_text_widget.insert(tk.END, "Para que serve: ", 'bold')
        help_text_widget.insert(tk.END, "Compara o perfil físico-químico geral dos proteomas. Permite ver se uma espécie tem, em média, mais proteínas de membrana (GRAVY > 0) ou mais proteínas aromáticas.\n", 'normal')
        help_text_widget.insert(tk.END, "Como usar: ", 'bold')
        help_text_widget.insert(tk.END, "Selecione 1 ou mais arquivos. Gera boxplots comparativos.\n", 'normal')
        
        # --- Análise 12 ---
        help_text_widget.insert(tk.END, "12: Gráfico de Neutralidade (GC12 vs GC3)\n", 'h3')
        help_text_widget.insert(tk.END, "O que faz: ", 'bold')
        help_text_widget.insert(tk.END, "Plota o conteúdo GC das posições 1 e 2 (GC12) contra o conteúdo GC da posição 3 (GC3). Calcula uma linha de regressão sobre os pontos.\n", 'normal')
        help_text_widget.insert(tk.END, "Para que serve: ", 'bold')
        help_text_widget.insert(tk.END, "Semelhante ao Gráfico de Wright (Análise 7), mas foca no balanço entre mutação e seleção. Uma inclinação (slope) próxima de 1 sugere que o GC3 é ditado pela pressão mutacional (neutralidade). Uma inclinação próxima de 0 sugere forte seleção natural nas posições 1 e 2, quebrando a neutralidade.\n", 'normal')
        help_text_widget.insert(tk.END, "Como usar: ", 'bold')
        help_text_widget.insert(tk.END, "Selecione 2 ou mais arquivos (necessário para a regressão).\n", 'normal')

        # --- Análise 13 ---
        help_text_widget.insert(tk.END, "13: Composição de Dinucleotídeos\n", 'h3')
        help_text_widget.insert(tk.END, "O que faz: ", 'bold')
        help_text_widget.insert(tk.END, "Calcula a frequência de todos os 16 pares de nucleotídeos (AA, AT, GC, CG, etc.) no genoma completo (não apenas CDS).\n", 'normal')
        help_text_widget.insert(tk.END, "Para que serve: ", 'bold')
        help_text_widget.insert(tk.END, "Revela 'assinaturas genômicas'. Por exemplo, genomas de vertebrados tendem a ter uma sub-representação de 'CG' (supressão de CpG). É muito usado em virologia para comparar o vírus ao seu hospedeiro.\n", 'normal')
        help_text_widget.insert(tk.END, "Como usar: ", 'bold')
        help_text_widget.insert(tk.END, "Selecione 1 ou mais arquivos. Gera um heatmap comparativo.\n", 'normal')
        
        # --- Análise 14 ---
        help_text_widget.insert(tk.END, "14: Gráfico de Paridade PR2 (A3/T3 vs G3/C3)\n", 'h3')
        help_text_widget.insert(tk.END, "O que faz: ", 'bold')
        help_text_widget.insert(tk.END, "Plota A3/(A3+T3) contra G3/(G3+C3) para cada gene individualmente. O centro do gráfico (0.5, 0.5) representa ausência de viés.\n", 'normal')
        help_text_widget.insert(tk.END, "Para que serve: ", 'bold')
        help_text_widget.insert(tk.END, "Analisa o viés mutacional e de seleção que atua de forma diferente nas posições A/T e G/C da terceira base do códon. Desvios do centro indicam viés.\n", 'normal')
        help_text_widget.insert(tk.END, "Como usar: ", 'bold')
        help_text_widget.insert(tk.END, "Selecione 1 ou mais arquivos. Gera um gráfico de dispersão com todos os genes.\n", 'normal')

        # --- Análise 15 (ATUALIZADA) ---
        help_text_widget.insert(tk.END, "15: Índice de Adaptação ao tRNA (tAI)\n", 'h3')
        help_text_widget.insert(tk.END, "O que faz: ", 'bold')
        help_text_widget.insert(tk.END, "Calcula o 'tRNA Adaptation Index' (tAI) usando as regras de pareamento **wobble** (dos Reis et al.), que é o padrão-ouro. Mede a eficiência traducional com base na contagem de genes de tRNA.\n", 'normal')
        help_text_widget.insert(tk.END, "Para que serve: ", 'bold')
        help_text_widget.insert(tk.END, "Mede o quão 'otimizado' um gene é para ser traduzido pela maquinaria de tRNAs disponível. Genes com tAI alto são (teoricamente) traduzidos mais eficientemente.\n", 'normal')
        help_text_widget.insert(tk.END, "Como usar: ", 'bold')
        help_text_widget.insert(tk.END, "Selecione 1 ou mais arquivos. A contagem de tRNAs é agregada. **Escolha o Super-reino (Bactéria/Eucarioto)** nas opções para usar as regras de wobble corretas.\n", 'normal')

        # --- Análise 16 ---
        help_text_widget.insert(tk.END, "16: Análise de Motifs Upstream\n", 'h3')
        help_text_widget.insert(tk.END, "O que faz: ", 'bold')
        help_text_widget.insert(tk.END, "Busca por k-mers (pequenas sequências, ex: 'TATAAT') que aparecem com alta frequência nas regiões 'upstream' (antes do início) de todos os genes.\n", 'normal')
        help_text_widget.insert(tk.END, "Para que serve: ", 'bold')
        help_text_widget.insert(tk.END, "Útil para encontrar possíveis sítios de ligação de fatores de transcrição, como promotores (ex: TATA box) ou sítios Shine-Dalgarno.\n", 'normal')
        help_text_widget.insert(tk.END, "Como usar: ", 'bold')
        help_text_widget.insert(tk.END, "Selecione 1 ou mais arquivos. Opções para definir a distância (ex: 200pb) e o tamanho do k-mer (ex: 6pb) aparecerão na aba de análise.\n", 'normal')

        # --- Análise 17 ---
        help_text_widget.insert(tk.END, "17: Análise de MFE (Estrutura 5')\n", 'h3')
        help_text_widget.insert(tk.END, "O que faz: ", 'bold')
        help_text_widget.insert(tk.END, "Calcula a Energia Livre Mínima (MFE) da região 5' (primeiros N nucleotídeos) de todos os CDS. Requer a instalação do 'ViennaRNA'.\n", 'normal')
        help_text_widget.insert(tk.END, "Para que serve: ", 'bold')
        help_text_widget.insert(tk.END, "Estruturas de RNA muito estáveis (MFE muito negativo) no início de um gene podem bloquear fisicamente o ribossomo, impedindo a tradução. Esta análise identifica genomas ou genes com potencial 'bloqueio' traducional.\n", 'normal')
        help_text_widget.insert(tk.END, "Como usar: ", 'bold')
        help_text_widget.insert(tk.END, "Selecione 1 ou mais arquivos. Opções para definir o tamanho da região (ex: 50pb) aparecerão na aba de análise.\n", 'normal')

        # --- ANÁLISES NOVAS (18, 19) ---
        help_text_widget.insert(tk.END, "18: Comparação de Grupos de Genes\n", 'h3')
        help_text_widget.insert(tk.END, "O que faz: ", 'bold')
        help_text_widget.insert(tk.END, "Compara as métricas de viés (ENC, GC3, CAI) entre dois grupos de genes definidos pelo usuário. Realiza testes estatísticos (Mann-Whitney U, Qui-Quadrado).\n", 'normal')
        help_text_widget.insert(tk.END, "Para que serve: ", 'bold')
        help_text_widget.insert(tk.END, "Permite testar hipóteses, ex: 'Genes ribossomais (Grupo 1) têm viés significativamente maior que genes de membrana (Grupo 2)?'.\n", 'normal')
        help_text_widget.insert(tk.END, "Como usar: ", 'bold')
        help_text_widget.insert(tk.END, "Selecione 1 ou mais arquivos. Nas opções, cole a lista de genes (locus_tags ou nomes) para o Grupo 1 e Grupo 2. O filtro de gene global é ignorado.\n", 'normal')

        help_text_widget.insert(tk.END, "19: Correlação com Expressão (RNA-Seq)\n", 'h3')
        help_text_widget.insert(tk.END, "O que faz: ", 'bold')
        help_text_widget.insert(tk.END, "Carrega um arquivo CSV/TSV de dados de expressão (ex: TPM, RPKM) e correlaciona (Spearman) os valores de expressão com as métricas de CUB (ENC, CAI).\n", 'normal')
        help_text_widget.insert(tk.END, "Para que serve: ", 'bold')
        help_text_widget.insert(tk.END, "É o teste clássico para validar se o viés de códons (medido pelo CAI/ENC) é realmente impulsionado pela seleção para alta expressão.\n", 'normal')
        help_text_widget.insert(tk.END, "Como usar: ", 'bold')
        help_text_widget.insert(tk.END, "Selecione 1 ou mais arquivos. Nas opções, carregue seu arquivo de expressão, especifique a coluna com o identificador do gene (deve bater com o locus_tag/gene) e a coluna com o valor de expressão.\n", 'normal')


        help_text_widget.configure(state='disabled')
        return tab

    def _create_status_bar(self):
        """Cria a barra de status inferior"""
        status_frame = ttk.Frame(self.root, relief="sunken", style='Status.TFrame')
        status_frame.pack(fill=tk.X, side=tk.BOTTOM, padx=0, pady=0, ipady=3)
        
        self.status_label = ttk.Label(status_frame, text="Pronto", anchor="w", style='Status.TLabel')
        self.status_label.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=10, pady=2)
        
        self.progress_bar = ttk.Progressbar(status_frame, mode='determinate', style='TProgressbar', length=200)
        self.progress_bar.pack(side=tk.RIGHT, padx=10, pady=2)

    # --- Bindings e Eventos ---

    def _bind_events(self):
        """Vincula eventos de widgets a funções"""
        self.analysis_combo.bind('<<ComboboxSelected>>', self._on_analysis_selected)
        
        # Atalhos
        self.root.bind('<Control-o>', lambda e: self._on_browse_input())
        self.root.bind('<Control-s>', lambda e: self._on_browse_output())
        self.root.bind('<Control-q>', lambda e: self._on_closing())
        self.root.bind('<F5>', lambda e: self._on_run_analysis())

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
            messagebox.showwarning("Aviso", "Selecione uma pasta de entrada primeiro.", parent=self.root)
            return
        
        self.file_listbox.delete(0, tk.END)
        self.file_paths_map.clear()
        
        patterns = ["*.gbk", "*.gb", "*.gbff"]
        files = []
        for pattern in patterns:
            files.extend(glob.glob(os.path.join(folder, pattern)))
        
        if not files:
            self._write_to_console(f"Nenhum arquivo .gbk ou .gb encontrado em '{folder}'\n", "warning")
            self.file_status_label.config(text="0 arquivos")
            return
        
        self._write_to_console(f"Encontrados {len(files)} arquivos .gbk/.gb.\n", "success")
        
        for file_path in sorted(files):
            filename = os.path.basename(file_path)
            self.file_listbox.insert(tk.END, filename)
            self.file_paths_map[filename] = file_path
        
        self.file_status_label.config(text=f"{len(files)} arquivos")
        self._on_select_all() # Selecionar todos por padrão
        
        self.root.event_generate("<<FilesLoaded>>")
    
    def _on_select_all(self):
        self.file_listbox.selection_set(0, tk.END)
    
    def _on_select_none(self):
        self.file_listbox.selection_clear(0, tk.END)

    # --- NOVO CALLBACK ---
    def _on_load_expression_file(self):
        """Callback para o botão 'Carregar Arquivo de Expressão'."""
        
        filepath = filedialog.askopenfilename(
            title="Selecione o arquivo de expressão (CSV ou TSV)",
            filetypes=[("Arquivos de Texto", "*.csv;*.tsv;*.txt"), ("Todos os Arquivos", "*.*")],
            parent=self.root
        )
        
        if not filepath:
            return
            
        try:
            # Tentar detectar separador
            with open(filepath, 'r') as f:
                first_line = f.readline()
                if '\t' in first_line:
                    sep = '\t'
                elif ',' in first_line:
                    sep = ','
                elif ';' in first_line:
                    sep = ';'
                else:
                    sep = r'\s+' # Fallback para whitespace
            
            self.expression_dataframe = pd.read_csv(filepath, sep=sep)
            
            filename = os.path.basename(filepath)
            self.expression_file_var.set(f"Carregado: {filename} ({len(self.expression_dataframe)} linhas)")
            
            self._write_to_console(f"Arquivo de expressão '{filename}' carregado com sucesso.\n", "success")
            self._write_to_console(f"Colunas detectadas: {list(self.expression_dataframe.columns)}\n")
            self.status_queue.put(("message", f"Arquivo de expressão '{filename}' carregado."))
            
        except Exception as e:
            self.expression_dataframe = None
            self.expression_file_var.set("Falha ao carregar o arquivo.")
            messagebox.showerror("Erro ao Carregar", 
                                 f"Não foi possível ler o arquivo:\n{e}",
                                 parent=self.root)
            self._write_to_console(f"Falha ao carregar arquivo de expressão: {e}\n", "error")
    # --- FIM DO NOVO CALLBACK ---


    def _on_analysis_selected(self, event=None):
        """Atualiza a descrição e as opções da análise"""
        selection = self.analysis_combo.get()
        if not selection:
            return
            
        data = self.analysis_data[selection]
        
        # Atualizar descrição
        self.analysis_desc_text.config(state="normal")
        self.analysis_desc_text.delete(1.0, tk.END)
        desc = f"Descrição:\n{data['description']}\n\n"
        req = f"Requisito de Arquivos: {data['files_required']}"
        self.analysis_desc_text.insert(1.0, desc)
        self.analysis_desc_text.insert(tk.END, req, ('bold')) # Destaca o requisito
        self.analysis_desc_text.config(state="disabled")
        self.analysis_desc_text.tag_configure('bold', font=('Calibri', 10, 'bold'), foreground=self.ERROR_COLOR)
        
        # Limpar opções antigas
        for widget in self.analysis_options_frame.winfo_children():
            widget.destroy()
            
        # Adicionar opções específicas (se houver)
        if data['id'] == '1':
            # Garantir que o filtro global esteja habilitado
            self.gene_filter_text.config(state='normal', bg='#FDFDFD')
            
            filter_frame = ttk.Frame(self.analysis_options_frame, style='Content.TFrame')
            filter_frame.pack(fill=tk.X, pady=5)
            ttk.Label(filter_frame, text="Filtro de Códon de Início (para Análise CDS):", style='Content.TLabel').pack(side=tk.LEFT, padx=(5,10))
            ttk.Entry(filter_frame, textvariable=self.filter_cds_var, width=10).pack(side=tk.LEFT, padx=5)

        # --- INÍCIO DO NOVO BLOCO 15 (Opção tAI) ---
        elif data['id'] == '15':
            # Garantir que o filtro global esteja habilitado
            self.gene_filter_text.config(state='normal', bg='#FDFDFD')
            
            options_frame = ttk.Frame(self.analysis_options_frame, style='Content.TFrame')
            options_frame.pack(fill=tk.X, pady=5)
            
            ttk.Label(options_frame, text="Super-reino (Regras de Wobble):", style='Content.TLabel').pack(side=tk.LEFT, padx=(5,10))
            
            kingdom_combo = ttk.Combobox(options_frame, textvariable=self.super_kingdom_var, state="readonly", width=15)
            kingdom_combo['values'] = ["Bactéria", "Eucarioto"]
            kingdom_combo.pack(side=tk.LEFT, padx=5)
        # --- FIM DO NOVO BLOCO 15 ---

        elif data['id'] == '16':
            # Garantir que o filtro global esteja habilitado
            self.gene_filter_text.config(state='normal', bg='#FDFDFD')
            
            # Variáveis específicas para esta análise
            self.upstream_dist_var = tk.StringVar(value="200")
            self.kmer_size_var = tk.StringVar(value="6")
            
            options_frame = ttk.Frame(self.analysis_options_frame, style='Content.TFrame')
            options_frame.pack(fill=tk.X, pady=5)
            
            # Opção de Distância Upstream
            ttk.Label(options_frame, text="Distância Upstream (pb):", style='Content.TLabel').pack(side=tk.LEFT, padx=(5,10))
            ttk.Entry(options_frame, textvariable=self.upstream_dist_var, width=6).pack(side=tk.LEFT, padx=5)
            
            # Opção de Tamanho do K-mer
            ttk.Label(options_frame, text="Tamanho do K-mer:", style='Content.TLabel').pack(side=tk.LEFT, padx=(15,10))
            ttk.Entry(options_frame, textvariable=self.kmer_size_var, width=4).pack(side=tk.LEFT, padx=5)
        
        elif data['id'] == '17':
            # Garantir que o filtro global esteja habilitado
            self.gene_filter_text.config(state='normal', bg='#FDFDFD')
            
            # Variáveis específicas para esta análise
            self.mfe_region_var = tk.StringVar(value="50")
            
            options_frame = ttk.Frame(self.analysis_options_frame, style='Content.TFrame')
            options_frame.pack(fill=tk.X, pady=5)
            
            # Opção de Distância Upstream
            ttk.Label(options_frame, text="Região 5' (pb):", style='Content.TLabel').pack(side=tk.LEFT, padx=(5,10))
            ttk.Entry(options_frame, textvariable=self.mfe_region_var, width=6).pack(side=tk.LEFT, padx=5)
        
        # --- INÍCIO DO NOVO BLOCO 18 (Comparação de Grupos) ---
        elif data['id'] == '18':
            # Esta análise usa widgets de texto próprios, desabilitar o global
            self.gene_filter_text.config(state='disabled', bg='#E0E0E0')
            self._write_to_console("Filtro global desabilitado. Use as caixas 'Grupo 1' e 'Grupo 2'.\n", "info")

            # Frame para Grupo 1
            group1_frame = ttk.LabelFrame(self.analysis_options_frame, text="Grupo de Genes 1", padding="5")
            group1_frame.pack(fill=tk.X, pady=(5,10), expand=True)
            
            ttk.Label(group1_frame, text="Cole a lista de locus_tags ou nomes (um por linha):", style='Content.TLabel').pack(anchor="w", padx=5)
            self.gene_list_1_text = scrolledtext.ScrolledText(group1_frame, height=6, wrap=tk.WORD, 
                                                              font=('Courier', 10), relief="solid", borderwidth=1,
                                                              bg="#FFFFFF", fg=self.TEXT_COLOR)
            self.gene_list_1_text.pack(fill=tk.X, expand=True, padx=5, pady=5)
            
            # Frame para Grupo 2
            group2_frame = ttk.LabelFrame(self.analysis_options_frame, text="Grupo de Genes 2", padding="5")
            group2_frame.pack(fill=tk.X, pady=5, expand=True)
            
            ttk.Label(group2_frame, text="Cole a lista de locus_tags ou nomes (um por linha):", style='Content.TLabel').pack(anchor="w", padx=5)
            self.gene_list_2_text = scrolledtext.ScrolledText(group2_frame, height=6, wrap=tk.WORD, 
                                                              font=('Courier', 10), relief="solid", borderwidth=1,
                                                              bg="#FFFFFF", fg=self.TEXT_COLOR)
            self.gene_list_2_text.pack(fill=tk.X, expand=True, padx=5, pady=5)

        # --- INÍCIO DO NOVO BLOCO 19 (Correlação com Expressão) ---
        elif data['id'] == '19':
            # Esta análise pode usar o filtro global
            self.gene_filter_text.config(state='normal', bg='#FDFDFD')
            
            expr_frame = ttk.LabelFrame(self.analysis_options_frame, text="Configuração de Expressão", padding="10")
            expr_frame.pack(fill=tk.X, pady=5, expand=True)
            
            # Linha 1: Botão e Status do Arquivo
            btn_frame = ttk.Frame(expr_frame, style='Content.TFrame')
            btn_frame.pack(fill=tk.X)
            
            ttk.Button(btn_frame, text="Carregar Arquivo de Expressão (.csv/.tsv)...", 
                       command=self._on_load_expression_file).pack(side=tk.LEFT, padx=5, pady=5)
            
            ttk.Label(btn_frame, textvariable=self.expression_file_var, style='Content.TLabel', 
                      font=('Calibri', 9, 'italic')).pack(side=tk.LEFT, padx=10)
            
            # Linha 2: Configuração das Colunas
            cols_frame = ttk.Frame(expr_frame, style='Content.TFrame')
            cols_frame.pack(fill=tk.X, pady=(10,5))
            
            ttk.Label(cols_frame, text="Coluna de Gene (ex: locus_tag):", style='Content.TLabel').pack(side=tk.LEFT, padx=5)
            ttk.Entry(cols_frame, textvariable=self.expr_gene_col_var, width=20).pack(side=tk.LEFT, padx=5)
            
            ttk.Label(cols_frame, text="Coluna de Expressão (ex: TPM):", style='Content.TLabel').pack(side=tk.LEFT, padx=15)
            ttk.Entry(cols_frame, textvariable=self.expr_value_col_var, width=15).pack(side=tk.LEFT, padx=5)

        # --- FIM DOS NOVOS BLOCOS ---
        else:
            # Garantir que o filtro global seja reabilitado se não for a Análise 18
            self.gene_filter_text.config(state='normal', bg='#FDFDFD')


    def _on_run_analysis(self):
        """Inicia a validação e a thread de análise"""
        
        # 1. Obter seleções
        files = self._get_selected_files()
        output_folder = self.output_folder_var.get()
        analysis_name = self.analysis_combo.get()
        analysis_data = self.analysis_data.get(analysis_name, {}) # Obter dados da análise
        
        # --- INÍCIO DA MODIFICAÇÃO (Variáveis de Análise) ---
        gene_list = None
        extra_args = {} # Dicionário para argumentos especiais
        
        # Obter lista de genes (Filtro Global)
        # Só lemos isso se NÃO for a Análise 18
        if analysis_data.get('id') != '18':
            gene_list_raw = self.gene_filter_text.get(1.0, tk.END)
            gene_list = set(line.strip() for line in gene_list_raw.splitlines() if line.strip())
            if not gene_list:
                gene_list = None # Passar None se estiver vazio
            else:
                self._write_to_console(f"Aplicando filtro global de {len(gene_list)} genes.\n", "info")
        # --- FIM DA MODIFICAÇÃO ---

        # 2. Validar Entradas
        if not analysis_name:
            messagebox.showwarning("Aviso", "Selecione uma análise para executar.", parent=self.root)
            return
            
        if not files:
            messagebox.showwarning("Aviso", "Selecione pelo menos um arquivo na lista.", parent=self.root)
            return
            
        if not output_folder or not os.path.isdir(output_folder):
            messagebox.showwarning("Aviso", "Selecione uma pasta de saída válida.", parent=self.root)
            return
            
        # 3. Validar Requisitos de Arquivos
        req = analysis_data.get('files_required', '1+')
        num_files = len(files)
        
        if req == '1' and num_files != 1:
            messagebox.showerror("Erro", f"A análise '{analysis_name}' requer EXATAMENTE 1 arquivo. Você selecionou {num_files}.", parent=self.root)
            return
        if req == '2' and num_files != 2:
            messagebox.showerror("Erro", f"A análise '{analysis_name}' requer EXATAMENTE 2 arquivos. Você selecionou {num_files}.", parent=self.root)
            return
        if req == '2+' and num_files < 2:
            messagebox.showerror("Erro", f"A análise '{analysis_name}' requer 2 OU MAIS arquivos. Você selecionou {num_files}.", parent=self.root)
            return
            
        # 4. Validar e Obter Argumentos Específicos
        try:
            # --- INÍCIO DO NOVO BLOCO (Análise 15) ---
            if analysis_data.get('id') == '15':
                extra_args['super_kingdom'] = self.super_kingdom_var.get()
            # --- FIM DO NOVO BLOCO (Análise 15) ---
            
            elif analysis_data.get('id') == '16':
                extra_args['upstream_dist'] = int(self.upstream_dist_var.get())
                extra_args['kmer_size'] = int(self.kmer_size_var.get())
            elif analysis_data.get('id') == '17':
                extra_args['mfe_region_length'] = int(self.mfe_region_var.get())
            
            # --- INÍCIO DO NOVO BLOCO (Validação 18 e 19) ---
            elif analysis_data.get('id') == '18':
                # Obter listas de genes dos widgets específicos
                gene_list_1 = set(line.strip() for line in self.gene_list_1_text.get(1.0, tk.END).splitlines() if line.strip())
                gene_list_2 = set(line.strip() for line in self.gene_list_2_text.get(1.0, tk.END).splitlines() if line.strip())
                
                if not gene_list_1 or not gene_list_2:
                    messagebox.showerror("Erro", "A Análise 18 requer que AMBOS os 'Grupo 1' e 'Grupo 2' sejam preenchidos.", parent=self.root)
                    return
                
                extra_args['gene_list_1'] = gene_list_1
                extra_args['gene_list_2'] = gene_list_2
                gene_list = None # Ignora o filtro global
            
            elif analysis_data.get('id') == '19':
                gene_col = self.expr_gene_col_var.get()
                expr_col = self.expr_value_col_var.get()

                if self.expression_dataframe is None:
                    messagebox.showerror("Erro", "Carregue um arquivo de expressão para a Análise 19.", parent=self.root)
                    return
                if not gene_col or not expr_col:
                    messagebox.showerror("Erro", "Especifique os nomes das colunas de 'Gene' e 'Expressão'.", parent=self.root)
                    return
                if gene_col not in self.expression_dataframe.columns:
                    messagebox.showerror("Erro", f"Coluna de Gene '{gene_col}' não encontrada no arquivo de expressão.", parent=self.root)
                    return
                if expr_col not in self.expression_dataframe.columns:
                    messagebox.showerror("Erro", f"Coluna de Expressão '{expr_col}' não encontrada no arquivo de expressão.", parent=self.root)
                    return
                
                extra_args['expression_data'] = self.expression_dataframe
                extra_args['gene_col'] = gene_col
                extra_args['expr_col'] = expr_col
            # --- FIM DO NOVO BLOCO ---
                
        except ValueError:
            messagebox.showerror("Erro", "Opções de análise inválidas. Verifique se os números são inteiros.", parent=self.root)
            return
        
        # 5. Se tudo OK, iniciar thread
        self._start_analysis_thread(files, output_folder, analysis_name, gene_list, extra_args)

    def _update_synth_host_list(self, event=None):
        """Atualiza a lista do combobox de hospedeiros."""
        hosts = sorted(self.file_paths_map.keys())
        self.synth_host_combo['values'] = hosts
        if hosts:
            self.synth_host_var.set(hosts[0])

    # --- INÍCIO: NOVAS FUNÇÕES (BIOLOGIA SINTÉTICA) ---

    def _get_host_data(self, host_filename):
        """
        Função auxiliar para obter dados de viés (RSCU, Contagens)
        para um único arquivo hospedeiro, sob demanda.
        """
        if host_filename not in self.file_paths_map:
            print(f"❌ Erro: Arquivo hospedeiro '{host_filename}' não encontrado no mapa.")
            return None
            
        host_filepath = self.file_paths_map[host_filename]
        genetic_code_id = self._get_genetic_code_id()
        
        print(f"  Calculando dados de viés para o hospedeiro: {host_filename}...")
        self.status_queue.put(("message", f"Calculando dados de {host_filename}..."))
        
        # Usamos 'gene_list=None' aqui, pois as ferramentas de otimização
        # geralmente devem ser baseadas no genoma COMPLETO do hospedeiro.
        all_data = processar_genomas_para_analise_vies(
            [host_filepath], 
            genetic_code_id, 
            self.status_queue, 
            gene_list=None
        )
        
        if not all_data:
            print(f"❌ Falha ao processar dados do hospedeiro: {host_filename}")
            return None
            
        # Retorna o primeiro (e único) item
        return all_data[list(all_data.keys())[0]]

    def _on_run_optimization(self):
        """Inicia a thread de OTIMIZAÇÃO"""
        self._start_optimization_thread(mode="optimize")

    def _on_run_harmonization(self):
        """Inicia a thread de HARMONIZAÇÃO"""
        self._start_optimization_thread(mode="harmonize")

    def _start_optimization_thread(self, mode="optimize"):
        """Valida e inicia a thread da ferramenta de Biologia Sintética."""
        
        # 1. Obter seleções
        host_filename = self.synth_host_var.get()
        input_seq = self.synth_input_seq.get(1.0, tk.END).strip().replace("\n", "")
        
        # 2. Validar
        if not host_filename:
            messagebox.showwarning("Aviso", "Selecione um genoma hospedeiro.", parent=self.root)
            return
            
        if len(input_seq) < 10:
            messagebox.showwarning("Aviso", "Insira uma sequência de DNA válida (mínimo 10bp).", parent=self.root)
            return
            
        # 3. Mudar para o console, desabilitar botões
        self.notebook.select(3) # Mudar para a aba do console (agora é a 4ª, índice 3)
        self.run_optimize_btn.config(state='disabled')
        self.run_harmonize_btn.config(state='disabled')
        self.progress_bar['value'] = 0
        self.status_label.config(text=f"Iniciando {mode}...")
        
        genetic_code_id = self._get_genetic_code_id()

        # 4. Iniciar thread
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
            
            # 1. Obter dados do hospedeiro (RSCU, Contagens)
            # Esta é a parte demorada (leitura de arquivo)
            host_data = self._get_host_data(host_filename)
            if not host_data:
                raise Exception("Não foi possível obter os dados do hospedeiro.")
            
            self.status_queue.put(("progress", 50))
            
            # 2. Executar a função de backend
            result_seq = ""
            if mode == "optimize":
                self.status_queue.put(("message", "Otimizando sequência..."))
                result_seq = otimizar_sequencia_codons(input_seq, host_data['rscu'], genetic_code_id)
                
            elif mode == "harmonize":
                self.status_queue.put(("message", "Harmonizando sequência..."))
                result_seq = harmonizar_sequencia_codons(input_seq, host_data['counts'], genetic_code_id)
            
            # 3. Enviar resultado de volta para a GUI
            self.status_queue.put(("optimization_complete", result_seq))
            print("✅ Ferramenta concluída com sucesso.")

        except Exception as e:
            print(f"\n❌ ERRO DURANTE A OPERAÇÃO: {e}", "error")
            import traceback
            print(traceback.format_exc())
            # Certifique-se de reabilitar os botões em caso de falha
            self.status_queue.put(("done", None)) 
            
        finally:
            # 'done' re-habilita os botões principais, vamos criar um
            # comando específico para re-habilitar os botões da ferramenta
            self.status_queue.put(("tool_done", None))

    # --- FIM: NOVAS FUNÇÕES (BIOLOGIA SINTÉTICA) ---

    # --- Lógica de Thread e Análise ---

    def _start_analysis_thread(self, files, output_folder, analysis_name, gene_list, extra_args):
        """Prepara e inicia a thread de análise de backend"""
        
        # Limpar console e mudar para a aba
        self.console_text.config(state="normal")
        self.console_text.delete(1.0, tk.END)
        self.console_text.config(state="disabled")
        self.notebook.select(3) # Mudar para a aba do console (Índice 3)
        
        # Desabilitar botão e resetar progresso
        self.run_button.config(state='disabled')
        self.progress_bar['value'] = 0
        self.status_label.config(text=f"Iniciando: {analysis_name}...")
        
        # Obter dados de configuração
        genetic_code_id = self._get_genetic_code_id()
        analysis_data = self.analysis_data[analysis_name]
        
        # O bloco que preenchia 'extra_args' foi movido para _on_run_analysis
        # onde a validação acontece.

        # Iniciar thread
        thread = threading.Thread(
            target=self._analysis_thread_target,
            args=(files, output_folder, analysis_name, analysis_data, genetic_code_id, gene_list, extra_args),
            daemon=True
        )
        thread.start()

    def _analysis_thread_target(self, files, output_folder, analysis_name, analysis_data, genetic_code_id, gene_list, extra_args):
            """A função que executa na thread (backend)"""
            try:
                print(f"🚀 INICIANDO ANÁLISE: {analysis_name}")
                if gene_list:
                    print(f"FILTRO ATIVO: Analisando apenas {len(gene_list)} genes especificados.")
                print(f"{'='*60}")
                print(f"Arquivos: {len(files)}")
                print(f"Pasta de Saída: {output_folder}")
                print(f"Tabela Genética: {genetic_code_id}")
                print(f"{'='*60}")
                
                # --- Roteamento da Análise ---
                
                # Análises que usam a função de processamento unificado (3, 4, 5, 6, 7, 9, 12)
                if analysis_data['id'] in ['3', '4', '5', '6', '7', '9', '12']:
                    
                    # Processar os genomas (contagens, rscu, enc, etc.)
                    all_bias_data = processar_genomas_para_analise_vies(files, genetic_code_id, self.status_queue, gene_list)
                    
                    if not all_bias_data:
                        raise Exception("Falha ao processar dados de viés. Verifique os arquivos de entrada.")
                    
                    # Chamar a função de plotagem específica
                    analysis_function = analysis_data['function']
                    if analysis_function:
                        # Análises 3 e 6 precisam do genetic_code_id
                        if analysis_data['id'] in ['3', '6']: 
                            analysis_function(all_bias_data, output_folder, genetic_code_id, self.status_queue)
                        # Análise 12 e outras
                        else:
                            analysis_function(all_bias_data, output_folder, self.status_queue)
                        
                # Análise 1: Estatística e CDS (Função própria)
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
                
                # Análise 2: Listar Genes (Função própria)
                elif analysis_data['id'] == '2':
                    df_genes = listar_genes_do_arquivo(files[0], self.status_queue, gene_list) # files[0] é garantido pela validação
                    if not df_genes.empty:
                        self._print_dataframe_limitado(df_genes, "Lista de Genes")
                        df_genes.to_csv(os.path.join(output_folder, f"lista_genes_{os.path.basename(files[0])}.csv"), index=False, sep=';')
                
                # Análise 8: Composição Genômica (Função própria)
                elif analysis_data['id'] == '8':
                    analise_composicao_genomica(files, output_folder, self.status_queue)

                # Análise 10: Viés de Pares de Códons (Função própria)
                elif analysis_data['id'] == '10':
                    analise_codon_pair_bias(files, output_folder, genetic_code_id, self.status_queue, gene_list)

                # Análise 11: GRAVY & Aromo (Função própria)
                elif analysis_data['id'] == '11':
                    analise_gravy_aromo(files, output_folder, genetic_code_id, self.status_queue, gene_list)

                # Análise 13: Composição de Dinucleotídeos (Função própria)
                elif analysis_data['id'] == '13':
                    analise_dinucleotide_composition(files, output_folder, self.status_queue)

                # Análise 14: PR2 Plot (Feature c)
                elif analysis_data['id'] == '14':
                    analise_pr2_plot(files, output_folder, genetic_code_id, self.status_queue, gene_list)

                # Análise 15: tAI (Feature a) - ATUALIZADO
                elif analysis_data['id'] == '15':
                    kingdom_name = extra_args.get('super_kingdom', 'Bactéria') # Padrão 'Bactéria'
                    analise_tAI(files, output_folder, genetic_code_id, self.status_queue, gene_list, 
                                super_kingdom=kingdom_name)
                
                # Análise 16: Motifs Upstream (Feature b)
                elif analysis_data['id'] == '16':
                    analise_motifs_upstream(files, output_folder, self.status_queue, gene_list, 
                                            extra_args['upstream_dist'], extra_args['kmer_size'])

                # Análise 17: MFE (Feature c - Estrutura)
                elif analysis_data['id'] == '17':
                    analise_mfe_iniciacao(files, output_folder, genetic_code_id, self.status_queue, gene_list, 
                                          extra_args['mfe_region_length'])
                
                # --- INÍCIO DO NOVO BLOCO (18 e 19) ---
                # Análise 18: Comparação de Grupos de Genes
                elif analysis_data['id'] == '18':
                    analise_comparativa_dois_grupos(files, output_folder, genetic_code_id, self.status_queue, 
                                                    extra_args['gene_list_1'], extra_args['gene_list_2'])
                
                # Análise 19: Correlação com Expressão
                elif analysis_data['id'] == '19':
                    analise_correlacao_expressao(files, output_folder, genetic_code_id, self.status_queue, 
                                                 gene_list, # Passa o filtro global (pode ser None)
                                                 extra_args['expression_data'], 
                                                 extra_args['gene_col'], 
                                                 extra_args['expr_col'])
                # --- FIM DO NOVO BLOCO ---
                    
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
        # Processar stdout
        try:
            while True:
                text = self.stdout_queue.get_nowait()
                # Determinar a tag de cor
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
        
        # Processar status
        try:
            command, data = self.status_queue.get_nowait()
            
            if command == "done":
                self.run_button.config(state='normal')
                self.progress_bar['value'] = 100
                self.status_label.config(text="Análise concluída com sucesso!")
            
            elif command == "tool_done":
                # Comando separado para reabilitar os botões da ferramenta
                self.run_optimize_btn.config(state='normal')
                self.run_harmonize_btn.config(state='normal')
                # Não sobrescreva a mensagem de "concluído" se houver
                if "Processando" in self.status_label.cget("text"):
                    self.status_label.config(text="Operação finalizada.")
            
            elif command == "optimization_complete":
                # Recebe a sequência otimizada
                self.synth_output_seq.config(state='normal')
                self.synth_output_seq.delete(1.0, tk.END)
                self.synth_output_seq.insert(1.0, data)
                self.synth_output_seq.config(state='disabled')
                
                # Mudar para a aba da ferramenta
                self.notebook.select(1) # A aba "Biologia Sintética" é o índice 1
                self.status_label.config(text="Sequência processada com sucesso!")
                self.progress_bar['value'] = 100
                
            elif command == "progress":
                self.progress_bar['value'] = data
                self.status_label.config(text=f"Processando... {data}%")
                
            elif command == "image_ready":
                image_path, title = data
                self._display_image(image_path, title)
                
            elif command == "message":
                self.status_label.config(text=data)
                
        except queue.Empty:
            pass
        
        # Agendar próxima verificação
        self.root.after(100, self._process_queues)

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
            ax.set_title(title, fontsize=14, fontweight='bold', color=self.TEXT_COLOR)
            self.fig.tight_layout()
            self.canvas.draw()
            
            self.current_images.append((image_path, title))
            self.image_status_label.config(text=f"Exibindo: {title}")
            
            # Mudar para a aba do visualizador (Índice 2)
            self.notebook.select(2)
            
        except Exception as e:
            self._write_to_console(f"❌ Erro ao exibir imagem {image_path}: {e}\n", "error")

    # --- Funções Utilitárias ---
    
    def _get_selected_files(self):
        """Retorna a lista de caminhos completos dos arquivos selecionados"""
        indices = self.file_listbox.curselection()
        return [self.file_paths_map[self.file_listbox.get(i)] for i in indices]
    
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

    def _show_about(self):
        """Mostra a janela 'Sobre'"""
        messagebox.showinfo("Sobre Kódon-E v1.0",
            "Kódon-E v1.0 - Análise Visual de Genomas\n\n"
            "Interface redesenhada para uma experiência de usuário "
            "limpa, moderna e intuitiva.\n\n"
            "O backend utiliza BioPython, Pandas, Matplotlib/Seaborn, "
            "Scikit-learn e SciPy para análises robustas.",
            parent=self.root)
    
    def _on_closing(self):
        """Evento de fechamento da janela"""
        if messagebox.askokcancel("Sair", "Deseja sair do Kódon-E?", parent=self.root):
            sys.stdout = self.stdout_original
            self.root.destroy()
