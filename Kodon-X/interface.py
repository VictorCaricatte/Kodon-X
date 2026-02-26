import os
import sys
import threading
import queue
import glob
from datetime import datetime
import pandas as pd
from PIL import Image

import matplotlib
matplotlib.use('qtagg') 
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure

from PyQt6.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton, 
    QLineEdit, QComboBox, QCheckBox, QListWidget, QListWidgetItem, 
    QTabWidget, QTextEdit, QFileDialog, QMessageBox, QProgressBar, 
    QDockWidget, QSplitter, QFrame, QScrollArea, QGridLayout, QFormLayout
)
from PyQt6.QtCore import Qt, QTimer, pyqtSignal, QObject
from PyQt6.QtGui import QFont, QIcon, QAction, QColor, QTextCursor, QTextCharFormat


from core_utils import process_genomes_for_bias_analysis
from analysis_basic import process_aggregated_gbk, analyze_gbk_cds, list_genes_from_file, analyze_genomic_composition
from analysis_bias import generate_rscu_heatmap_and_table, comparative_rscu_analysis, rscu_correlation_analysis, generate_rscu_histograms, enc_gc3_analysis, optimal_rare_codons_analysis, neutrality_plot_analysis
from analysis_advanced import codon_pair_bias_analysis, gravy_aromo_analysis, dinucleotide_composition_analysis, pr2_plot_analysis, tai_analysis, upstream_motifs_analysis
from analysis_expression_structure import initiation_mfe_analysis, two_groups_comparative_analysis, expression_correlation_analysis
from synthetic_biology import optimize_codon_sequence, harmonize_codon_sequence

class Redirector(QObject):
    new_text = pyqtSignal(str)

    def __init__(self, queue):
        super().__init__()
        self.queue = queue
    
    def write(self, string):
        self.queue.put(string)
    
    def flush(self):
        pass

class KodonE_GUI(QMainWindow):
    
    files_loaded_signal = pyqtSignal()

    def __init__(self):
        super().__init__()
        
        self.setWindowTitle("Kodon-X")
        self.resize(1600, 900)
        
        
        self.SUCCESS_COLOR = "#43A047"
        self.ERROR_COLOR = "#E53935"
        self.WARNING_COLOR = "#FB8C00"
        self.INFO_COLOR = "#039BE5"
        
        
        self.stdout_queue = queue.Queue()
        self.status_queue = queue.Queue()
        
        
        self.file_paths_map = {}
        self.current_images = []
        self.expression_dataframe = None
        self.analysis_data = self._get_analysis_definitions()
        
        self.stdout_original = sys.stdout
        sys.stdout = Redirector(self.stdout_queue)
        
        self._apply_stylesheet()
        self._create_widgets()
        
        
        self.timer = QTimer(self)
        self.timer.timeout.connect(self._process_queues)
        self.timer.start(100) # 100ms
        
        self.files_loaded_signal.connect(self._update_synth_host_list)
        
        self._write_to_console("Welcome to Kodon-X!\n", "info") 
        self._write_to_console("Please select your input folder and click 'Load Files'.\n\n")

    def _get_analysis_definitions(self):
        return {
            '1: Statistics and CDS': {
                'id': '1',
                'description': ('Performs an initial scan on all genomes. \nResearch Application: Essential for data *quality control (QC)*.'),
                'files_required': '1+', 
                'function': None 
            },
            '2: Gene Listing': {
                'id': '2',
                'description': ('Lists all genes, products, and locus tags present in a single GenBank file.'),
                'files_required': '1', 
                'function': list_genes_from_file
            },
            '3: Individual RSCU Heatmap': {
                'id': '3',
                'description': ('Calculates and visualizes the RSCU (Relative Synonymous Codon Usage) for a single genome.'),
                'files_required': '1', 
                'function': generate_rscu_heatmap_and_table
            },
            '4: Comparative RSCU': {
                'id': '4',
                'description': ('Groups genomes based on the similarity of their codon bias (RSCU).'),
                'files_required': '2+', 
                'function': comparative_rscu_analysis
            },
            '5: RSCU Correlation (2 Genomes)': {
                'id': '5',
                'description': ('Calculates the Pearson correlation (R) between the RSCU patterns of exactly two genomes.'),
                'files_required': '2', 
                'function': rscu_correlation_analysis
            },
            '6: Comparative RSCU Histograms': {
                'id': '6',
                'description': ('Generates detailed visualizations of codon usage, including Box Plots by amino acid and stacked bar charts by species.'),
                'files_required': '2+', 
                'function': generate_rscu_histograms
            },
            '7: ENC vs GC3 Analysis (Wright Plot)': {
                'id': '7',
                'description': ("Generates the 'Wright Plot' (ENC vs. GC3)."),
                'files_required': '1+', 
                'function': enc_gc3_analysis
            },
            '8: Genomic Composition': {
                'id': '8',
                'description': ('Detailed analysis of nucleotide composition (A, T, G, C), total GC content, and genome size for multiple files.'),
                'files_required': '1+', 
                'function': analyze_genomic_composition
            },
            '9: Optimal, Rare Codons and CAI': {
                'id': '9',
                'description': ('Identifies preferred/rare codons and compares the CAI (Codon Adaptation Index) between species.'),
                'files_required': '1+', 
                'function': optimal_rare_codons_analysis
            },
            '10: Codon Pair Bias (CPB)': {
                'id': '10',
                'description': ('Analyzes the frequency of adjacent codon pairs (e.g., ATG-CGC).'),
                'files_required': '1+', 
                'function': codon_pair_bias_analysis
            },
            '11: Physicochemical Analysis (GRAVY & Aromo)': {
                'id': '11',
                'description': ('Calculates the GRAVY score (hydropathicity) and Aromo (aromaticity) for all genes.'),
                'files_required': '1+', 
                'function': gravy_aromo_analysis
            },
            '12: Neutrality Plot (GC12 vs GC3)': {
                'id': '12',
                'description': ('Plots the GC content of positions 1+2 (GC12) against GC3.'),
                'files_required': '2+', 
                'function': neutrality_plot_analysis
            },
            '13: Dinucleotide Composition': {
                'id': '13',
                'description': ('Calculates and compares the frequency of all 16 dinucleotides (AA, AT, GC, etc.) in the complete genome.'),
                'files_required': '1+', 
                'function': dinucleotide_composition_analysis
            },
            '14: PR2 Parity Plot (A3/T3 vs G3/C3)': {
                'id': '14',
                'description': ('Plots A3/(A3+T3) vs G3/(G3+C3) for each individual gene.'),
                'files_required': '1+', 
                'function': pr2_plot_analysis
            },
            '15: tRNA Adaptation Index (tAI)': {
                'id': '15',
                'description': ('Calculates the \'tRNA Adaptation Index\' (tAI) using wobble rules.'),
                'files_required': '1+', 
                'function': tai_analysis 
            },
            '16: Upstream Motifs Analysis': {
                'id': '16',
                'description': ('Searches for short sequences (k-mers) overrepresented in regions before the start of each gene (upstream).'),
                'files_required': '1+', 
                'function': upstream_motifs_analysis
            },
            '17: MFE Analysis (5\' Structure)': {
                'id': '17',
                'description': ("Calculates the Minimum Free Energy (MFE) of the 5' region (start) of genes."),
                'files_required': '1+', 
                'function': initiation_mfe_analysis
            },
            '18: Gene Group Comparison': {
                'id': '18',
                'description': ('Compares CUB metrics (ENC, CAI) between two user-defined gene sets.'),
                'files_required': '1+', 
                'function': two_groups_comparative_analysis
            },
            '19: Correlation with Expression (RNA-Seq)': {
                'id': '19',
                'description': ('Correlates bias (CAI, ENC) with expression data provided by the user.'),
                'files_required': '1+', 
                'function': expression_correlation_analysis
            }
        }

    def _apply_stylesheet(self):
        
        qss = """
        QMainWindow {
            background-color: #1a1b26;
        }
        QWidget {
            color: #a9b1d6;
            font-family: "Segoe UI", "Calibri", sans-serif;
            font-size: 13px;
        }
        QDockWidget {
            titlebar-close-icon: url(close.png);
            titlebar-normal-icon: url(float.png);
            background: #24283b;
            color: #7aa2f7;
            font-weight: bold;
        }
        QDockWidget::title {
            text-align: left;
            background: #1f2335;
            padding-left: 10px;
            padding-top: 5px;
            padding-bottom: 5px;
            border-bottom: 1px solid #16161e;
        }
        QTabWidget::pane {
            border: 1px solid #1f2335;
            background: #1a1b26;
        }
        QTabBar::tab {
            background: #1f2335;
            color: #565f89;
            padding: 8px 15px;
            border-top-left-radius: 4px;
            border-top-right-radius: 4px;
            margin-right: 2px;
        }
        QTabBar::tab:selected {
            background: #24283b;
            color: #7aa2f7;
            border-bottom: 2px solid #7aa2f7;
        }
        QTabBar::tab:hover {
            background: #292e42;
            color: #c0caf5;
        }
        QPushButton {
            background-color: #3d59a1;
            color: #ffffff;
            border: none;
            border-radius: 4px;
            padding: 6px 12px;
            font-weight: bold;
        }
        QPushButton:hover {
            background-color: #7aa2f7;
        }
        QPushButton:pressed {
            background-color: #2ac3de;
        }
        QPushButton:disabled {
            background-color: #3b4261;
            color: #565f89;
        }
        QLineEdit, QTextEdit, QComboBox, QListWidget {
            background-color: #16161e;
            color: #c0caf5;
            border: 1px solid #3b4261;
            border-radius: 4px;
            padding: 4px;
        }
        QLineEdit:focus, QTextEdit:focus, QComboBox:focus, QListWidget:focus {
            border: 1px solid #7aa2f7;
        }
        QListWidget::item {
            padding: 5px;
        }
        QListWidget::item:selected {
            background-color: #3d59a1;
            color: white;
        }
        QScrollBar:vertical {
            border: none;
            background: #16161e;
            width: 10px;
            margin: 0px 0px 0px 0px;
        }
        QScrollBar::handle:vertical {
            background: #3b4261;
            min-height: 20px;
            border-radius: 5px;
        }
        QScrollBar::handle:vertical:hover {
            background: #565f89;
        }
        QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {
            border: none;
            background: none;
            height: 0px;
        }
        QProgressBar {
            border: 1px solid #3b4261;
            border-radius: 4px;
            text-align: center;
            background-color: #16161e;
            color: white;
        }
        QProgressBar::chunk {
            background-color: #7aa2f7;
            width: 20px;
        }
        QStatusBar {
            background-color: #16161e;
            color: #a9b1d6;
        }
        """
        self.setStyleSheet(qss)

    def _create_widgets(self):
        
        self.setCentralWidget(QWidget())
        self.centralWidget().hide()

        
        self.setDockOptions(QMainWindow.DockOption.AllowNestedDocks | QMainWindow.DockOption.AnimatedDocks)

        
        self._create_left_dock()
        
        
        self._create_bottom_dock()
        
        
        self._create_main_dock()

        
        self.statusBar = self.statusBar()
        self.status_label = QLabel("Ready")
        self.statusBar.addWidget(self.status_label, 1)
        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 100)
        self.progress_bar.setValue(0)
        self.progress_bar.setFixedSize(200, 15)
        self.progress_bar.hide()
        self.statusBar.addPermanentWidget(self.progress_bar)

    def _create_left_dock(self):
        dock = QDockWidget("Project Workspace", self)
        dock.setAllowedAreas(Qt.DockWidgetArea.LeftDockWidgetArea | Qt.DockWidgetArea.RightDockWidgetArea)
        dock.setFeatures(QDockWidget.DockWidgetFeature.DockWidgetMovable | QDockWidget.DockWidgetFeature.DockWidgetFloatable)
        
        container = QWidget()
        layout = QVBoxLayout(container)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(15)

        
        input_group = QFrame()
        input_layout = QVBoxLayout(input_group)
        input_layout.setContentsMargins(0, 0, 0, 0)
        
        lbl_input = QLabel("📂 Step 1: Input Files")
        lbl_input.setFont(QFont("Calibri", 11, QFont.Weight.Bold))
        input_layout.addWidget(lbl_input)
        
        h_input = QHBoxLayout()
        self.input_folder_edit = QLineEdit()
        self.input_folder_edit.setReadOnly(True)
        h_input.addWidget(self.input_folder_edit)
        btn_browse_in = QPushButton("Browse...")
        btn_browse_in.clicked.connect(self._on_browse_input)
        h_input.addWidget(btn_browse_in)
        input_layout.addLayout(h_input)
        
        btn_load = QPushButton("Load Files")
        btn_load.clicked.connect(self._on_load_files)
        input_layout.addWidget(btn_load)
        
        self.file_list_widget = QListWidget()
        input_layout.addWidget(self.file_list_widget)
        
        h_sel = QHBoxLayout()
        btn_sel_all = QPushButton("All")
        btn_sel_all.clicked.connect(self._on_select_all)
        btn_sel_none = QPushButton("None")
        btn_sel_none.clicked.connect(self._on_select_none)
        self.lbl_file_count = QLabel("0 files")
        h_sel.addWidget(btn_sel_all)
        h_sel.addWidget(btn_sel_none)
        h_sel.addStretch()
        h_sel.addWidget(self.lbl_file_count)
        input_layout.addLayout(h_sel)
        
        layout.addWidget(input_group)

        
        output_group = QFrame()
        output_layout = QVBoxLayout(output_group)
        output_layout.setContentsMargins(0, 0, 0, 0)
        
        lbl_out = QLabel("💾 Step 2: Output Folder")
        lbl_out.setFont(QFont("Calibri", 11, QFont.Weight.Bold))
        output_layout.addWidget(lbl_out)
        
        h_out = QHBoxLayout()
        self.output_folder_edit = QLineEdit()
        self.output_folder_edit.setReadOnly(True)
        h_out.addWidget(self.output_folder_edit)
        btn_browse_out = QPushButton("Browse...")
        btn_browse_out.clicked.connect(self._on_browse_output)
        h_out.addWidget(btn_browse_out)
        output_layout.addLayout(h_out)
        
        layout.addWidget(output_group)

        
        gen_group = QFrame()
        gen_layout = QVBoxLayout(gen_group)
        gen_layout.setContentsMargins(0, 0, 0, 0)
        
        lbl_gen = QLabel("⚙️ Step 3: Genetic Table")
        lbl_gen.setFont(QFont("Calibri", 11, QFont.Weight.Bold))
        gen_layout.addWidget(lbl_gen)
        
        self.genetic_combo = QComboBox()
        self.genetic_combo.addItems([
            "1: Standard (Universal)", 
            "2: Vertebrate Mitochondrial",
            "4: Mold/Protozoan Mitochondrial", 
            "11: Bacterial/Plant Plastid"
        ])
        gen_layout.addWidget(self.genetic_combo)
        
        layout.addWidget(gen_group)
        layout.addStretch()

        dock.setWidget(container)
        self.addDockWidget(Qt.DockWidgetArea.LeftDockWidgetArea, dock)

    def _create_bottom_dock(self):
        dock = QDockWidget("Output Console", self)
        dock.setAllowedAreas(Qt.DockWidgetArea.BottomDockWidgetArea)
        
        container = QWidget()
        layout = QVBoxLayout(container)
        layout.setContentsMargins(5, 5, 5, 5)
        
        self.console_text = QTextEdit()
        self.console_text.setReadOnly(True)
        self.console_text.setFont(QFont("Consolas", 10))
        
        self.console_text.setStyleSheet("background-color: #0d0d12; color: #a9b1d6; border: none;")
        layout.addWidget(self.console_text)
        
        dock.setWidget(container)
        self.addDockWidget(Qt.DockWidgetArea.BottomDockWidgetArea, dock)

    def _create_main_dock(self):
        dock = QDockWidget("Main Canvas", self)
        dock.setFeatures(QDockWidget.DockWidgetFeature.NoDockWidgetFeatures) # Fixa no centro
        dock.setTitleBarWidget(QWidget()) # Esconde a barra de título
        
        self.notebook = QTabWidget()
        
        self.tab_analysis = QWidget()
        self.tab_synth = QWidget()
        self.tab_viz = QWidget()
        self.tab_help = QWidget()
        
        self.notebook.addTab(self.tab_analysis, " Run Analysis")
        self.notebook.addTab(self.tab_synth, " Synthetic Biology")
        self.notebook.addTab(self.tab_viz, " Chart Viewer")
        self.notebook.addTab(self.tab_help, "❓ Help")
        
        self._create_analysis_tab(self.tab_analysis)
        self._create_synth_bio_tab(self.tab_synth)
        self._create_visualization_tab(self.tab_viz)
        self._create_help_tab(self.tab_help)
        
        dock.setWidget(self.notebook)
        self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, dock)

    def _create_analysis_tab(self, tab):
        layout = QVBoxLayout(tab)
        
        lbl_sel = QLabel("Select Analysis:")
        lbl_sel.setFont(QFont("Calibri", 12, QFont.Weight.Bold))
        layout.addWidget(lbl_sel)
        
        self.analysis_combo = QComboBox()
        self.analysis_combo.addItems(list(self.analysis_data.keys()))
        self.analysis_combo.currentTextChanged.connect(self._on_analysis_selected)
        layout.addWidget(self.analysis_combo)
        
        self.analysis_desc_text = QTextEdit()
        self.analysis_desc_text.setReadOnly(True)
        self.analysis_desc_text.setMaximumHeight(100)
        layout.addWidget(self.analysis_desc_text)
        
        
        self.analysis_options_container = QWidget()
        self.analysis_options_layout = QVBoxLayout(self.analysis_options_container)
        self.analysis_options_layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.analysis_options_container)
        
        # Filtro Global
        self.gene_filter_frame = QFrame()
        gf_layout = QVBoxLayout(self.gene_filter_frame)
        gf_layout.setContentsMargins(0, 10, 0, 0)
        
        lbl_gf = QLabel("Global Gene Filter (Optional)")
        lbl_gf.setFont(QFont("Calibri", 11, QFont.Weight.Bold))
        gf_layout.addWidget(lbl_gf)
        
        lbl_gf_desc = QLabel("Paste a list of locus_tags or gene names (one per line).\nIf filled, ALL analyses will only use this subset of genes.")
        gf_layout.addWidget(lbl_gf_desc)
        
        self.gene_filter_text = QTextEdit()
        self.gene_filter_text.setFont(QFont("Consolas", 10))
        gf_layout.addWidget(self.gene_filter_text)
        layout.addWidget(self.gene_filter_frame)
        
        self.run_button = QPushButton(" Run Selected Analysis")
        self.run_button.setFont(QFont("Calibri", 12, QFont.Weight.Bold))
        self.run_button.setMinimumHeight(50)
        self.run_button.clicked.connect(self._on_run_analysis)
        layout.addWidget(self.run_button)
        
        # Iniciar com a primeira análise selecionada
        self._on_analysis_selected(self.analysis_combo.currentText())

    def _create_synth_bio_tab(self, tab):
        layout = QHBoxLayout(tab)
        
        left_panel = QWidget()
        lp_layout = QVBoxLayout(left_panel)
        
        lbl_host = QLabel("1. Select Host:")
        lbl_host.setFont(QFont("Calibri", 11, QFont.Weight.Bold))
        lp_layout.addWidget(lbl_host)
        
        self.synth_host_combo = QComboBox()
        lp_layout.addWidget(self.synth_host_combo)
        
        lbl_seq = QLabel("2. Paste DNA Sequence (CDS):")
        lbl_seq.setFont(QFont("Calibri", 11, QFont.Weight.Bold))
        lp_layout.addWidget(lbl_seq)
        
        self.synth_input_seq = QTextEdit()
        self.synth_input_seq.setFont(QFont("Consolas", 10))
        lp_layout.addWidget(self.synth_input_seq)
        
        btn_layout = QHBoxLayout()
        self.run_optimize_btn = QPushButton("Optimize\n(Maximize RSCU)")
        self.run_optimize_btn.clicked.connect(self._on_run_optimization)
        self.run_harmonize_btn = QPushButton("Harmonize\n(Keep Ranking)")
        self.run_harmonize_btn.clicked.connect(self._on_run_harmonization)
        btn_layout.addWidget(self.run_optimize_btn)
        btn_layout.addWidget(self.run_harmonize_btn)
        lp_layout.addLayout(btn_layout)
        
        right_panel = QWidget()
        rp_layout = QVBoxLayout(right_panel)
        
        lbl_res = QLabel("Resulting Sequence:")
        lbl_res.setFont(QFont("Calibri", 11, QFont.Weight.Bold))
        rp_layout.addWidget(lbl_res)
        
        self.synth_output_seq = QTextEdit()
        self.synth_output_seq.setFont(QFont("Consolas", 10))
        self.synth_output_seq.setReadOnly(True)
        rp_layout.addWidget(self.synth_output_seq)
        
        layout.addWidget(left_panel)
        layout.addWidget(right_panel)

    def _create_visualization_tab(self, tab):
        layout = QVBoxLayout(tab)
        
        self.fig = Figure(figsize=(10, 8), dpi=100)
        self.fig.patch.set_facecolor("#1a1b26") # Fundo escuro
        self.fig.set_tight_layout(True)
        
        self.canvas = FigureCanvasQTAgg(self.fig)
        layout.addWidget(self.canvas, stretch=1)
        
        self.toolbar = NavigationToolbar2QT(self.canvas, self)
        # Ajuste de estilo da toolbar (cores padrão no dark theme costumam falhar)
        self.toolbar.setStyleSheet("QToolButton { background-color: #24283b; color: white; }")
        
        h_tool = QHBoxLayout()
        h_tool.addWidget(self.toolbar)
        self.image_status_label = QLabel("No chart loaded.")
        h_tool.addStretch()
        h_tool.addWidget(self.image_status_label)
        
        layout.addLayout(h_tool)

    def _create_help_tab(self, tab):
        layout = QVBoxLayout(tab)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        
        content = QWidget()
        c_layout = QVBoxLayout(content)
        
        def add_text(text, style="normal"):
            lbl = QLabel(text)
            lbl.setWordWrap(True)
            if style == "h1":
                lbl.setFont(QFont("Calibri", 20, QFont.Weight.Bold))
                lbl.setStyleSheet("color: #7aa2f7;")
            elif style == "h2":
                lbl.setFont(QFont("Calibri", 16, QFont.Weight.Bold))
            elif style == "h3":
                lbl.setFont(QFont("Calibri", 14, QFont.Weight.Bold))
                lbl.setStyleSheet("color: #2ac3de;")
            elif style == "bold":
                lbl.setFont(QFont("Calibri", 11, QFont.Weight.Bold))
            c_layout.addWidget(lbl)

        add_text("Kodon-X USER GUIDE", "h1")
        add_text("Welcome to the Kodon-X analysis interface. Follow this guide for an efficient workflow.", "normal")
        add_text("QUICK WORKFLOW", "h2")
        add_text("1. Load Files", "h3")
        add_text("Use the 'Step 1' panel to 'Browse...' your folder with .gbk or .gb files. Click 'Load Files'.", "normal")
        
        add_text("ANALYSIS DESCRIPTIONS", "h2")
        for analysis_name, data in self.analysis_data.items():
            add_text(analysis_name, "h3")
            add_text(f"File Requirement: {data['files_required']}", "bold")
            add_text(data['description'], "normal")
            
        c_layout.addStretch()
        scroll.setWidget(content)
        layout.addWidget(scroll)

    def _clear_layout(self, layout):
        if layout is not None:
            while layout.count():
                child = layout.takeAt(0)
                if child.widget() is not None:
                    child.widget().deleteLater()
                elif child.layout() is not None:
                    self._clear_layout(child.layout())

    def _on_browse_input(self):
        folder = QFileDialog.getExistingDirectory(self, "Select input folder with GenBank files")
        if folder:
            self.input_folder_edit.setText(folder)
            self._on_load_files()
            
    def _on_browse_output(self):
        folder = QFileDialog.getExistingDirectory(self, "Select output folder for results")
        if folder:
            self.output_folder_edit.setText(folder)

    def _on_load_files(self):
        folder = self.input_folder_edit.text()
        if not folder:
            QMessageBox.warning(self, "Warning", "Select an input folder first.")
            return
            
        self.file_list_widget.clear()
        self.file_paths_map.clear()
        
        patterns = ["*.gbk", "*.gb", "*.gbff"]
        files = []
        for pattern in patterns:
            files.extend(glob.glob(os.path.join(folder, pattern)))
            
        if not files:
            self._write_to_console(f"No .gbk or .gb files found in '{folder}'\n", "warning")
            self.lbl_file_count.setText("0 files")
            return
            
        self._write_to_console(f"Found {len(files)} .gbk/.gb files.\n", "success")
        
        for file_path in sorted(files):
            filename = os.path.basename(file_path)
            self.file_paths_map[filename] = file_path
            
            item = QListWidgetItem(filename)
            item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
            item.setCheckState(Qt.CheckState.Checked)
            self.file_list_widget.addItem(item)
            
        self.lbl_file_count.setText(f"{len(files)} files")
        self.files_loaded_signal.emit()

    def _on_select_all(self):
        for i in range(self.file_list_widget.count()):
            self.file_list_widget.item(i).setCheckState(Qt.CheckState.Checked)

    def _on_select_none(self):
        for i in range(self.file_list_widget.count()):
            self.file_list_widget.item(i).setCheckState(Qt.CheckState.Unchecked)

    def _on_load_expression_file(self):
        filepath, _ = QFileDialog.getOpenFileName(
            self, "Select expression file (CSV or TSV)", "", "Text Files (*.csv *.tsv *.txt);;All Files (*.*)"
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
            self.lbl_expr_file.setText(f"Loaded: {filename} ({len(self.expression_dataframe)} rows)")
            self._write_to_console(f"Expression file '{filename}' loaded successfully.\n", "success")
            self.status_queue.put(("message", f"Expression file '{filename}' loaded."))
        except Exception as e:
            self.expression_dataframe = None
            if hasattr(self, 'lbl_expr_file'):
                self.lbl_expr_file.setText("Failed to load file.")
            QMessageBox.critical(self, "Load Error", f"Could not read the file:\n{e}")
            self._write_to_console(f"Failed to load expression file: {e}\n", "error")

    def _on_analysis_selected(self, selection):
        if not selection: return
        
        data = self.analysis_data[selection]
        
        desc = f"Description:\n{data['description']}\n\nFile Requirement: {data['files_required']}"
        self.analysis_desc_text.setText(desc)
        
        self._clear_layout(self.analysis_options_layout)
        
        # Variáveis dinâmicas da GUI baseadas na análise
        self.gene_filter_text.setEnabled(True)
        
        if data['id'] == '1':
            w = QWidget()
            lo = QHBoxLayout(w)
            lo.addWidget(QLabel("Start Codon Filter:"))
            self.filter_cds_edit = QLineEdit("ATG")
            self.filter_cds_edit.setFixedWidth(80)
            lo.addWidget(self.filter_cds_edit)
            lo.addStretch()
            self.analysis_options_layout.addWidget(w)

        elif data['id'] == '15':
            w = QWidget()
            lo = QHBoxLayout(w)
            lo.addWidget(QLabel("Superkingdom (Wobble Rules):"))
            self.super_kingdom_combo = QComboBox()
            self.super_kingdom_combo.addItems(["Bacteria", "Eukaryote"])
            lo.addWidget(self.super_kingdom_combo)
            lo.addStretch()
            self.analysis_options_layout.addWidget(w)

        elif data['id'] == '16':
            w = QWidget()
            lo = QHBoxLayout(w)
            lo.addWidget(QLabel("Upstream Dist (bp):"))
            self.upstream_dist_edit = QLineEdit("200")
            self.upstream_dist_edit.setFixedWidth(60)
            lo.addWidget(self.upstream_dist_edit)
            
            lo.addWidget(QLabel("K-mer Size:"))
            self.kmer_size_edit = QLineEdit("6")
            self.kmer_size_edit.setFixedWidth(60)
            lo.addWidget(self.kmer_size_edit)
            lo.addStretch()
            self.analysis_options_layout.addWidget(w)
        
        elif data['id'] == '17':
            w = QWidget()
            lo = QHBoxLayout(w)
            lo.addWidget(QLabel("5' Region (bp):"))
            self.mfe_region_edit = QLineEdit("50")
            self.mfe_region_edit.setFixedWidth(60)
            lo.addWidget(self.mfe_region_edit)
            lo.addStretch()
            self.analysis_options_layout.addWidget(w)
            
        elif data['id'] == '18':
            self.gene_filter_text.setEnabled(False)
            w = QWidget()
            lo = QVBoxLayout(w)
            
            lo.addWidget(QLabel("Gene Group 1 (one per line):"))
            self.gene_list_1_text = QTextEdit()
            self.gene_list_1_text.setMaximumHeight(80)
            lo.addWidget(self.gene_list_1_text)
            
            lo.addWidget(QLabel("Gene Group 2 (one per line):"))
            self.gene_list_2_text = QTextEdit()
            self.gene_list_2_text.setMaximumHeight(80)
            lo.addWidget(self.gene_list_2_text)
            self.analysis_options_layout.addWidget(w)

        elif data['id'] == '19':
            w = QWidget()
            lo = QVBoxLayout(w)
            
            h1 = QHBoxLayout()
            btn_load_expr = QPushButton("Load Expression File...")
            btn_load_expr.clicked.connect(self._on_load_expression_file)
            h1.addWidget(btn_load_expr)
            self.lbl_expr_file = QLabel("No file loaded.")
            h1.addWidget(self.lbl_expr_file)
            h1.addStretch()
            lo.addLayout(h1)
            
            h2 = QHBoxLayout()
            h2.addWidget(QLabel("Gene Column:"))
            self.expr_gene_col_edit = QLineEdit("locus_tag")
            h2.addWidget(self.expr_gene_col_edit)
            
            h2.addWidget(QLabel("Expression Column:"))
            self.expr_value_col_edit = QLineEdit("TPM")
            h2.addWidget(self.expr_value_col_edit)
            h2.addStretch()
            lo.addLayout(h2)
            
            self.analysis_options_layout.addWidget(w)

    def _get_selected_files(self):
        selected = []
        for i in range(self.file_list_widget.count()):
            item = self.file_list_widget.item(i)
            if item.checkState() == Qt.CheckState.Checked:
                selected.append(self.file_paths_map[item.text()])
        return selected
        
    def _get_genetic_code_id(self):
        try:
            return int(self.genetic_combo.currentText().split(':')[0])
        except:
            return 1 

    def _on_run_analysis(self):
        files = self._get_selected_files()
        output_folder = self.output_folder_edit.text()
        analysis_name = self.analysis_combo.currentText()
        analysis_data = self.analysis_data.get(analysis_name, {})
        
        gene_list = None
        extra_args = {}
        
        if analysis_data.get('id') != '18':
            raw_text = self.gene_filter_text.toPlainText()
            gene_list = set(line.strip() for line in raw_text.splitlines() if line.strip())
            if not gene_list:
                gene_list = None
            else:
                self._write_to_console(f"Applying global filter of {len(gene_list)} genes.\n", "info")

        if not analysis_name:
            QMessageBox.warning(self, "Warning", "Select an analysis to run.")
            return
        if not files:
            QMessageBox.warning(self, "Warning", "Select at least one file in the list.")
            return
        if not output_folder or not os.path.isdir(output_folder):
            QMessageBox.warning(self, "Warning", "Select a valid output folder.")
            return
            
        req = analysis_data.get('files_required', '1+')
        num_files = len(files)
        
        if req == '1' and num_files != 1:
            QMessageBox.critical(self, "Error", f"Analysis requires EXACTLY 1 file. You selected {num_files}.")
            return
        if req == '2' and num_files != 2:
            QMessageBox.critical(self, "Error", f"Analysis requires EXACTLY 2 files. You selected {num_files}.")
            return
        if req == '2+' and num_files < 2:
            QMessageBox.critical(self, "Error", f"Analysis requires 2 OR MORE files. You selected {num_files}.")
            return
            
        try:
            if analysis_data.get('id') == '1':
                extra_args['filter_cds'] = self.filter_cds_edit.text()
            elif analysis_data.get('id') == '15':
                extra_args['super_kingdom'] = self.super_kingdom_combo.currentText()
            elif analysis_data.get('id') == '16':
                extra_args['upstream_dist'] = int(self.upstream_dist_edit.text())
                extra_args['kmer_size'] = int(self.kmer_size_edit.text())
            elif analysis_data.get('id') == '17':
                extra_args['mfe_region_length'] = int(self.mfe_region_edit.text())
            elif analysis_data.get('id') == '18':
                g1_text = self.gene_list_1_text.toPlainText()
                g2_text = self.gene_list_2_text.toPlainText()
                gene_list_1 = set(line.strip() for line in g1_text.splitlines() if line.strip())
                gene_list_2 = set(line.strip() for line in g2_text.splitlines() if line.strip())
                if not gene_list_1 or not gene_list_2:
                    QMessageBox.critical(self, "Error", "Analysis 18 requires BOTH groups to be filled.")
                    return
                extra_args['gene_list_1'] = gene_list_1
                extra_args['gene_list_2'] = gene_list_2
                gene_list = None
            elif analysis_data.get('id') == '19':
                gene_col = self.expr_gene_col_edit.text()
                expr_col = self.expr_value_col_edit.text()
                if self.expression_dataframe is None:
                    QMessageBox.critical(self, "Error", "Load an expression file for Analysis 19.")
                    return
                if not gene_col or not expr_col:
                    QMessageBox.critical(self, "Error", "Specify 'Gene' and 'Expression' columns.")
                    return
                if gene_col not in self.expression_dataframe.columns:
                    QMessageBox.critical(self, "Error", f"Gene column '{gene_col}' not found.")
                    return
                if expr_col not in self.expression_dataframe.columns:
                    QMessageBox.critical(self, "Error", f"Expression column '{expr_col}' not found.")
                    return
                extra_args['expression_data'] = self.expression_dataframe
                extra_args['gene_col'] = gene_col
                extra_args['expr_col'] = expr_col
        except ValueError:
            QMessageBox.critical(self, "Error", "Invalid analysis options. Ensure numbers are integers.")
            return
            
        self._start_analysis_thread(files, output_folder, analysis_name, gene_list, extra_args)

    def _update_synth_host_list(self):
        self.synth_host_combo.clear()
        hosts = sorted(self.file_paths_map.keys())
        self.synth_host_combo.addItems(hosts)

    def _get_host_data(self, host_filename):
        if host_filename not in self.file_paths_map:
            print(f"❌ Error: Host file '{host_filename}' not found.")
            return None
        host_filepath = self.file_paths_map[host_filename]
        genetic_code_id = self._get_genetic_code_id()
        print(f"  Calculating bias data for host: {host_filename}...")
        self.status_queue.put(("message", f"Calculating data for {host_filename}..."))
        all_data = process_genomes_for_bias_analysis(
            [host_filepath], 
            genetic_code_id, 
            self.status_queue, 
            gene_list=None
        )
        if not all_data:
            return None
        return all_data[list(all_data.keys())[0]]

    def _on_run_optimization(self):
        self._start_optimization_thread(mode="optimize")

    def _on_run_harmonization(self):
        self._start_optimization_thread(mode="harmonize")

    def _start_optimization_thread(self, mode="optimize"):
        host_filename = self.synth_host_combo.currentText()
        input_seq = self.synth_input_seq.toPlainText().strip().replace("\n", "")
        
        if not host_filename:
            QMessageBox.warning(self, "Warning", "Select a host genome.")
            return
        if len(input_seq) < 10:
            QMessageBox.warning(self, "Warning", "Enter a valid DNA sequence (minimum 10bp).")
            return
            
        self.run_optimize_btn.setEnabled(False)
        self.run_harmonize_btn.setEnabled(False)
        self.progress_bar.show()
        self.progress_bar.setRange(0, 0) # Indeterminado
        self.status_label.setText(f"Starting {mode}...")
        
        genetic_code_id = self._get_genetic_code_id()

        thread = threading.Thread(
            target=self._optimization_thread_target,
            args=(input_seq, host_filename, genetic_code_id, mode),
            daemon=True
        )
        thread.start()

    def _optimization_thread_target(self, input_seq, host_filename, genetic_code_id, mode):
        try:
            print(f"\n{'='*60}")
            print(f"🚀 STARTING SYNTHETIC BIOLOGY TOOL")
            print(f"   Mode: {mode}")
            print(f"   Host: {host_filename}")
            print(f"   Input Size: {len(input_seq)}bp")
            print(f"{'='*60}")
            
            host_data = self._get_host_data(host_filename)
            if not host_data:
                raise Exception("Could not obtain host data.")
            
            self.status_queue.put(("progress", 50))
            
            result_seq = ""
            if mode == "optimize":
                self.status_queue.put(("message", "Optimizing sequence..."))
                result_seq = optimize_codon_sequence(input_seq, host_data['rscu'], genetic_code_id)
            elif mode == "harmonize":
                self.status_queue.put(("message", "Harmonizing sequence..."))
                result_seq = harmonize_codon_sequence(input_seq, host_data['counts'], genetic_code_id)
            
            self.status_queue.put(("optimization_complete", result_seq))
            print("✅ Tool completed successfully.")
        except Exception as e:
            print(f"\n❌ ERROR DURING OPERATION: {e}", "error")
            import traceback
            print(traceback.format_exc())
            self.status_queue.put(("done", None)) 
        finally:
            self.status_queue.put(("tool_done", None))

    def _start_analysis_thread(self, files, output_folder, analysis_name, gene_list, extra_args):
        self.console_text.clear()
        
        self.run_button.setEnabled(False)
        self.progress_bar.show()
        self.progress_bar.setRange(0, 0)
        self.status_label.setText(f"Starting: {analysis_name}...")
        
        genetic_code_id = self._get_genetic_code_id()
        analysis_data = self.analysis_data[analysis_name]
        
        thread = threading.Thread(
            target=self._analysis_thread_target,
            args=(files, output_folder, analysis_name, analysis_data, genetic_code_id, gene_list, extra_args),
            daemon=True
        )
        thread.start()

    def _print_dataframe_limitado(self, df, titulo):
        limite = 50
        print(f"\n--- {titulo} (First {min(len(df), limite)} rows) ---")
        if len(df) > limite:
            print(df.head(limite).to_string())
            print(f"\n... and {len(df) - limite} more rows (see full .csv in output folder).")
        else:
            print(df.to_string())
        print("--------------------------------------------------\n")

    def _analysis_thread_target(self, files, output_folder, analysis_name, analysis_data, genetic_code_id, gene_list, extra_args):
        try:
            print(f"🚀 STARTING ANALYSIS: {analysis_name}")
            if gene_list:
                print(f"FILTER ACTIVE: Analyzing only {len(gene_list)} specified genes.")
            print(f"{'='*60}")
            print(f"Files: {len(files)}")
            print(f"Output Folder: {output_folder}")
            print(f"Genetic Table: {genetic_code_id}")
            print(f"{'='*60}")
            
            if analysis_data['id'] in ['3', '4', '5', '6', '7', '9', '12']:
                all_bias_data = process_genomes_for_bias_analysis(files, genetic_code_id, self.status_queue, gene_list)
                if not all_bias_data:
                    raise Exception("Failed to process bias data. Check input files.")
                analysis_function = analysis_data['function']
                if analysis_function:
                    if analysis_data['id'] in ['3', '6']: 
                        analysis_function(all_bias_data, output_folder, genetic_code_id, self.status_queue)
                    else:
                        analysis_function(all_bias_data, output_folder, self.status_queue)
            elif analysis_data['id'] == '1':
                print("\n--- [Analysis 1] Part 1: Aggregated Statistics ---")
                df_stats = process_aggregated_gbk(files, self.status_queue)
                if not df_stats.empty:
                    self._print_dataframe_limitado(df_stats, "Aggregated Statistics")
                    df_stats.to_csv(os.path.join(output_folder, 'genome_statistics.csv'), index=False, sep=';')
                print("\n--- [Analysis 1] Part 2: CDS Analysis ---")
                filtro = extra_args.get('filter_cds', 'ATG')
                df_cds = analyze_gbk_cds(files, filtro, self.status_queue, gene_list)
                if not df_cds.empty:
                    self._print_dataframe_limitado(df_cds, "CDS Analysis")
                    df_cds.to_csv(os.path.join(output_folder, 'cds_analysis.csv'), index=False, sep=';')
            elif analysis_data['id'] == '2':
                df_genes = list_genes_from_file(files[0], self.status_queue, gene_list)
                if not df_genes.empty:
                    self._print_dataframe_limitado(df_genes, "Gene List")
                    df_genes.to_csv(os.path.join(output_folder, f"gene_list_{os.path.basename(files[0])}.csv"), index=False, sep=';')
            elif analysis_data['id'] == '8':
                analyze_genomic_composition(files, output_folder, self.status_queue)
            elif analysis_data['id'] == '10':
                codon_pair_bias_analysis(files, output_folder, genetic_code_id, self.status_queue, gene_list)
            elif analysis_data['id'] == '11':
                gravy_aromo_analysis(files, output_folder, genetic_code_id, self.status_queue, gene_list)
            elif analysis_data['id'] == '13':
                dinucleotide_composition_analysis(files, output_folder, self.status_queue)
            elif analysis_data['id'] == '14':
                pr2_plot_analysis(files, output_folder, genetic_code_id, self.status_queue, gene_list)
            elif analysis_data['id'] == '15':
                kingdom_name = extra_args.get('super_kingdom', 'Bacteria')
                tai_analysis(files, output_folder, genetic_code_id, self.status_queue, gene_list, 
                            super_kingdom=kingdom_name)
            elif analysis_data['id'] == '16':
                upstream_motifs_analysis(files, output_folder, self.status_queue, gene_list, 
                                        extra_args['upstream_dist'], extra_args['kmer_size'])
            elif analysis_data['id'] == '17':
                initiation_mfe_analysis(files, output_folder, genetic_code_id, self.status_queue, gene_list, 
                                      extra_args['mfe_region_length'])
            elif analysis_data['id'] == '18':
                two_groups_comparative_analysis(files, output_folder, genetic_code_id, self.status_queue, 
                                                extra_args['gene_list_1'], extra_args['gene_list_2'])
            elif analysis_data['id'] == '19':
                expression_correlation_analysis(files, output_folder, genetic_code_id, self.status_queue, 
                                             gene_list, extra_args['expression_data'], 
                                             extra_args['gene_col'], extra_args['expr_col'])
            else:
                print(f"Warning: Analysis '{analysis_name}' not implemented.", "warning")
            print(f"\n{'='*60}")
            print(f"✅ ANALYSIS COMPLETED: {analysis_name}")
            print(f"📁 Results saved in: {output_folder}")
            print(f"{'='*60}")
        except Exception as e:
            print(f"\n❌ ERROR DURING ANALYSIS: {e}", "error")
            import traceback
            print(traceback.format_exc())
        finally:
            self.status_queue.put(("done", None))

    def _process_queues(self):
        # Processa o stdout customizado
        try:
            while True:
                text = self.stdout_queue.get_nowait()
                if "❌" in text or "Error" in text.capitalize() or "Failed" in text:
                    self._write_to_console(text, "error")
                elif "✅" in text or "Success" in text or "completed" in text.lower():
                    self._write_to_console(text, "success")
                elif "📊" in text or "📈" in text or "Starting" in text or "🚀" in text:
                    self._write_to_console(text, "info")
                elif "⚠️" in text or "Warning" in text:
                    self._write_to_console(text, "warning")
                else:
                    self._write_to_console(text)
        except queue.Empty:
            pass
        
        
        try:
            while True:
                command, data = self.status_queue.get_nowait()
                
                if command == "done":
                    self.run_button.setEnabled(True)
                    self.progress_bar.setRange(0, 100)
                    self.progress_bar.setValue(100)
                    self.status_label.setText("Analysis completed successfully!")
                
                elif command == "tool_done":
                    self.run_optimize_btn.setEnabled(True)
                    self.run_harmonize_btn.setEnabled(True)
                    if "Processing" in self.status_label.text():
                        self.status_label.setText("Operation finished.")
                
                elif command == "optimization_complete":
                    self.synth_output_seq.setReadOnly(False)
                    self.synth_output_seq.clear()
                    self.synth_output_seq.setText(data)
                    self.synth_output_seq.setReadOnly(True)
                    self.notebook.setCurrentWidget(self.tab_synth)
                    self.status_label.setText("Sequence processed successfully!")
                    self.progress_bar.setRange(0, 100)
                    self.progress_bar.setValue(100)
                    
                elif command == "progress":
                    self.progress_bar.setRange(0, 100)
                    self.progress_bar.setValue(int(data))
                    self.status_label.setText(f"Processing... {data}%")
                    
                elif command == "image_ready":
                    image_path, title = data
                    self._display_image(image_path, title)
                    
                elif command == "message":
                    self.status_label.setText(data)
        except queue.Empty:
            pass

    def _write_to_console(self, text, tag=None):
        cursor = self.console_text.textCursor()
        cursor.movePosition(QTextCursor.MoveOperation.End)
        
        format = QTextCharFormat()
        if tag == "error":
            format.setForeground(QColor(self.ERROR_COLOR))
        elif tag == "success":
            format.setForeground(QColor(self.SUCCESS_COLOR))
        elif tag == "info":
            format.setForeground(QColor(self.INFO_COLOR))
        elif tag == "warning":
            format.setForeground(QColor(self.WARNING_COLOR))
        else:
            format.setForeground(QColor("#a9b1d6"))
            
        cursor.insertText(text, format)
        self.console_text.setTextCursor(cursor)
        self.console_text.ensureCursorVisible()

    def _display_image(self, image_path, title):
        try:
            self.fig.clear()
            img = Image.open(image_path)
            ax = self.fig.add_subplot(111)
            ax.imshow(img)
            ax.axis('off')
            ax.set_title(title, fontsize=14, fontweight='bold', color="#c0caf5") 
            ax.set_facecolor("#1a1b26")
            self.fig.patch.set_facecolor("#1a1b26")
            self.fig.set_tight_layout(True)
            self.canvas.draw()
            
            self.current_images.append((image_path, title))
            self.image_status_label.setText(f"Displaying: {title}")
            
            self.notebook.setCurrentWidget(self.tab_viz)
            
        except Exception as e:
            self._write_to_console(f"❌ Error displaying image {image_path}: {e}\n", "error")

    def closeEvent(self, event):
        reply = QMessageBox.question(
            self, 'Exit', 'Do you want to exit Kodon-X?',
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No, QMessageBox.StandardButton.No
        )
        if reply == QMessageBox.StandardButton.Yes:
            sys.stdout = self.stdout_original
            event.accept()
        else:
            event.ignore()