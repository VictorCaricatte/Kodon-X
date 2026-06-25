import os
import sys
import threading
import queue
import glob
import platform
import subprocess
from datetime import datetime
import pandas as pd

from PyQt6.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton, 
    QLineEdit, QComboBox, QListWidget, QListWidgetItem, 
    QTabWidget, QTextEdit, QFileDialog, QMessageBox, QProgressBar, 
    QDockWidget, QFrame, QScrollArea,
    QMenuBar, QMenu, QSizePolicy, QGraphicsView, QGraphicsScene
)
from PyQt6.QtCore import Qt, QTimer, pyqtSignal, QObject
from PyQt6.QtGui import QFont, QColor, QTextCursor, QTextCharFormat, QPixmap, QPainter

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

from core_utils import process_genomes_for_bias_analysis, calculate_rscu
from constants import GENETIC_CODE_TABLES, CODON_GRID_ORDER
from analysis_basic import process_aggregated_gbk, analyze_gbk_cds, list_genes_from_file, analyze_genomic_composition
from analysis_bias import generate_rscu_heatmap_and_table, comparative_rscu_analysis, rscu_correlation_analysis, enc_gc3_analysis, optimal_rare_codons_analysis, neutrality_plot_analysis
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
        self.image_history = []
        self.current_image_index = -1
        
        self.expression_dataframe = None
        self.analysis_data = self._get_analysis_definitions()
        
        self.stdout_original = sys.stdout
        sys.stdout = Redirector(self.stdout_queue)
        
        self._apply_stylesheet()
        self._create_widgets()
        
        self.timer = QTimer(self)
        self.timer.timeout.connect(self._process_queues)
        self.timer.start(100)
        
        self.files_loaded_signal.connect(self._update_synth_host_list)
        
        self._write_to_console("Welcome to Kodon-X!\n", "info") 
        self._write_to_console("Please select your input folder and click 'Load Files'.\n\n")

    def _get_analysis_definitions(self):
        return {
            '1: Statistics and CDS': {
                'id': '1',
                'description': ('Performs an initial scan on all genomes.\nResearch Application: Essential for data quality control (QC).'),
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
                'description': ('Calculates and visualizes the RSCU (Relative Synonymous Codon Usage) for a single genome.\nOutputs heatmap, RSCU table, and global metrics (ENC, GC3, CAI).'),
                'files_required': '1',
                'function': generate_rscu_heatmap_and_table
            },
            '4: Comparative RSCU': {
                'id': '4',
                'description': ('Unified comparative RSCU module (2+ genomes). Generates 6 figures:\n'
                                '  1. Hierarchical clustermap (64 codons × N species)\n'
                                '  2. PCA biplot (species scores + top-20 codon loading arrows)\n'
                                '  3. Comparative line plot per codon (grouped by amino acid)\n'
                                '  4. Box plot of RSCU distribution by amino acid\n'
                                '  5. Top-30 most variable codons — variance barplot\n'
                                '  6. RSCU heatmap for the 30 most divergent codons'),
                'files_required': '2+',
                'function': comparative_rscu_analysis
            },
            '5: RSCU Correlation (2 Genomes)': {
                'id': '5',
                'description': ('Calculates the Pearson correlation (R) between the RSCU patterns of exactly two genomes.\n'
                                'Outputs regression scatter, delta-RSCU barplot, and top-20 divergent codons heatmap.'),
                'files_required': '2',
                'function': rscu_correlation_analysis
            },
            '6: ENC vs GC3 Analysis (Wright Plot)': {
                'id': '6',
                'description': ('Generates the Wright Plot — ENC vs. GC3 per gene CDS.\n'
                                'Each point = one gene; red dashed curve = neutral expectation (Wright 1990).\n'
                                'Genes below the curve are under translational selection.'),
                'files_required': '1+',
                'function': enc_gc3_analysis
            },
            '7: Genomic Composition': {
                'id': '7',
                'description': ('Detailed analysis of nucleotide composition (A, T, G, C), total GC content, and genome size.'),
                'files_required': '1+',
                'function': analyze_genomic_composition
            },
            '8: Optimal, Rare Codons and CAI': {
                'id': '8',
                'description': ('Identifies preferred/rare codons and compares the CAI (Codon Adaptation Index) between species.'),
                'files_required': '1+',
                'function': optimal_rare_codons_analysis
            },
            '9: Codon Pair Bias (CPB)': {
                'id': '9',
                'description': ('Analyzes the frequency of adjacent codon pairs (e.g., ATG-CGC).\nIdentifies translation bottlenecks via codon pair scores (CPS).'),
                'files_required': '1+',
                'function': codon_pair_bias_analysis
            },
            '10: Physicochemical Analysis (GRAVY & Aromo)': {
                'id': '10',
                'description': ('Calculates the GRAVY score (Kyte-Doolittle hydropathicity) and Aromaticity for all genes.'),
                'files_required': '1+',
                'function': gravy_aromo_analysis
            },
            '11: Neutrality Plot (GC12 vs GC3)': {
                'id': '11',
                'description': ('Sueoka Neutrality Plot — GC12 vs GC3 per gene CDS, with one OLS regression per species.\n'
                                'Slope ≈ 1: mutation pressure dominant. Slope ≈ 0: translational selection dominant.'),
                'files_required': '2+',
                'function': neutrality_plot_analysis
            },
            '12: Dinucleotide Composition (ρ_XY)': {
                'id': '12',
                'description': ('Calculates the Karlin-Burge (1995) dinucleotide odds ratio ρ_XY = f(XY)/(f(X)×f(Y)).\n'
                                'ρ < 1 = under-represented (e.g., CpG suppression); ρ > 1 = over-represented.\n'
                                'Diverging heatmap centred at ρ = 1.0 (neutral).'),
                'files_required': '1+',
                'function': dinucleotide_composition_analysis
            },
            '13: PR2 Parity Plot (A3/T3 vs G3/C3)': {
                'id': '13',
                'description': ('Plots A3/(A3+T3) vs G3/(G3+C3) for each individual gene (Sueoka 1995).\nReveals strand-specific mutation biases.'),
                'files_required': '1+',
                'function': pr2_plot_analysis
            },
            '14: tRNA Adaptation Index (tAI)': {
                'id': '14',
                'description': ('Calculates the tRNA Adaptation Index (tAI) using wobble rules (dos Reis 2004).\n'
                                'Measures theoretical translation efficiency per gene.\n'
                                'Optional: set Superkingdom (Bacteria / Eukaryote) for wobble s-values.'),
                'files_required': '1+',
                'function': tai_analysis
            },
            '15: Upstream Motifs Analysis': {
                'id': '15',
                'description': ('Searches for k-mers over-represented in upstream regions before each gene start.\n'
                                'Optional: upstream distance (bp) and k-mer size.'),
                'files_required': '1+',
                'function': upstream_motifs_analysis
            },
            '16: MFE Analysis (5\' Structure)': {
                'id': '16',
                'description': ("Calculates the Minimum Free Energy (MFE) of the 5' region of each gene.\n"
                                "Tests the translational ramp hypothesis. Requires ViennaRNA installed."),
                'files_required': '1+',
                'function': initiation_mfe_analysis
            },
            '17: Gene Group Comparison': {
                'id': '17',
                'description': ('Compares CUB metrics (ENC, CAI) between two user-defined gene sets using\nMann-Whitney U test and Chi-squared test.'),
                'files_required': '1+',
                'function': two_groups_comparative_analysis
            },
            '18: Correlation with Expression (RNA-Seq)': {
                'id': '18',
                'description': ('Correlates bias metrics (CAI, ENC) with external expression data (e.g., RNA-Seq TPM)\nusing Spearman rank correlation.'),
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
        self.setDockOptions(QMainWindow.DockOption.AllowNestedDocks | QMainWindow.DockOption.AnimatedDocks)

        self._create_left_dock()
        self._create_bottom_dock()
        self._create_central_area() 
        
        self._create_menu_bar()

        self.statusBar = self.statusBar()
        self.status_label = QLabel("Ready")
        self.statusBar.addWidget(self.status_label, 1)
        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 100)
        self.progress_bar.setValue(0)
        self.progress_bar.setFixedSize(200, 15)
        self.progress_bar.hide()
        self.statusBar.addPermanentWidget(self.progress_bar)

    def _create_menu_bar(self):
        menu_bar = self.menuBar()
        view_menu = menu_bar.addMenu("View")
        
        if hasattr(self, 'dock_workspace'):
            view_menu.addAction(self.dock_workspace.toggleViewAction())
        if hasattr(self, 'dock_console'):
            view_menu.addAction(self.dock_console.toggleViewAction())

    def _create_left_dock(self):
        self.dock_workspace = QDockWidget("Project Workspace", self)
        self.dock_workspace.setAllowedAreas(Qt.DockWidgetArea.LeftDockWidgetArea | Qt.DockWidgetArea.RightDockWidgetArea)
        self.dock_workspace.setFeatures(QDockWidget.DockWidgetFeature.DockWidgetMovable | QDockWidget.DockWidgetFeature.DockWidgetFloatable | QDockWidget.DockWidgetFeature.DockWidgetClosable)
        
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

        self.dock_workspace.setWidget(container)
        self.addDockWidget(Qt.DockWidgetArea.LeftDockWidgetArea, self.dock_workspace)

    def _create_bottom_dock(self):
        self.dock_console = QDockWidget("Output Console", self)
        self.dock_console.setAllowedAreas(Qt.DockWidgetArea.BottomDockWidgetArea | Qt.DockWidgetArea.TopDockWidgetArea)
        self.dock_console.setFeatures(QDockWidget.DockWidgetFeature.DockWidgetMovable | QDockWidget.DockWidgetFeature.DockWidgetFloatable | QDockWidget.DockWidgetFeature.DockWidgetClosable)
        
        container = QWidget()
        layout = QVBoxLayout(container)
        layout.setContentsMargins(5, 5, 5, 5)
        
        self.console_text = QTextEdit()
        self.console_text.setReadOnly(True)
        self.console_text.setFont(QFont("Consolas", 11)) 
        self.console_text.setLineWrapMode(QTextEdit.LineWrapMode.NoWrap) 
        self.console_text.document().setMaximumBlockCount(2500) 
        self.console_text.setStyleSheet("background-color: #0d0d12; color: #a9b1d6; border: none;")
        
        self.console_text.setMinimumHeight(50)
        self.console_text.setMinimumWidth(50)
        self.console_text.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Ignored)
        
        layout.addWidget(self.console_text)
        
        self.dock_console.setWidget(container)
        self.addDockWidget(Qt.DockWidgetArea.BottomDockWidgetArea, self.dock_console)

    def _create_central_area(self):
        self.notebook = QTabWidget()
        
        self.tab_analysis = QWidget()
        self.tab_synth = QWidget()
        self.tab_viz = QWidget()
        self.tab_help = QWidget()
        
        self.notebook.addTab(self.tab_analysis, " Run Analysis")
        self.notebook.addTab(self.tab_synth, " Synthetic Biology")
        self.notebook.addTab(self.tab_viz, " Output")
        self.notebook.addTab(self.tab_help, "❓ Help")
        
        self._create_analysis_tab(self.tab_analysis)
        self._create_synth_bio_tab(self.tab_synth)
        self._create_visualization_tab(self.tab_viz)
        self._create_help_tab(self.tab_help)
        
        self.setCentralWidget(self.notebook)

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

        lbl_palette = QLabel("Chart Color Palette:")
        lbl_palette.setFont(QFont("Calibri", 11, QFont.Weight.Bold))
        layout.addWidget(lbl_palette)
        
        self.palette_combo = QComboBox()
        self.palette_combo.addItems([
            "viridis", "plasma", "inferno", "magma",
            "cividis", "colorblind", "Set2", "Pastel1",
            "Dark2", "Blues", "Reds", "Greens"
        ])
        layout.addWidget(self.palette_combo)
        
        self.gene_filter_frame = QFrame()
        gf_layout = QVBoxLayout(self.gene_filter_frame)
        gf_layout.setContentsMargins(0, 10, 0, 0)
        
        lbl_gf_layout = QHBoxLayout()
        lbl_gf = QLabel("Global Gene Filter (Optional)")
        lbl_gf.setFont(QFont("Calibri", 11, QFont.Weight.Bold))
        lbl_gf_layout.addWidget(lbl_gf)
        
        btn_load_gf = QPushButton("📂 Load File (TXT, FASTA, GBK)")
        btn_load_gf.clicked.connect(lambda: self._on_load_gene_list_file(self.gene_filter_text))
        lbl_gf_layout.addWidget(btn_load_gf)
        lbl_gf_layout.addStretch()
        gf_layout.addLayout(lbl_gf_layout)
        
        lbl_gf_desc = QLabel("Paste a list of locus_tags or gene names (one per line), or load a file.\nIf filled, ALL analyses will only use this subset of genes.")
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
        
        self._on_analysis_selected(self.analysis_combo.currentText())

    def _on_load_gene_list_file(self, target_text_edit):
        filepath, _ = QFileDialog.getOpenFileName(
            self, "Select File with Gene IDs", "", "All Supported (*.txt *.fasta *.fas *.fa *.gbk *.gb);;Text Files (*.txt);;FASTA Files (*.fasta *.fas *.fa);;GenBank Files (*.gbk *.gb)"
        )
        if not filepath:
            return
        
        extracted_ids = []
        ext = filepath.lower().split('.')[-1]
        
        try:
            if ext == 'txt':
                with open(filepath, 'r') as f:
                    for line in f:
                        val = line.strip()
                        if val: extracted_ids.append(val)
            elif ext in ['fasta', 'fas', 'fa']:
                for record in SeqIO.parse(filepath, "fasta"):
                    extracted_ids.append(record.id)
            elif ext in ['gbk', 'gb']:
                for record in SeqIO.parse(filepath, "genbank"):
                    for feature in record.features:
                        if feature.type == "CDS":
                            if "locus_tag" in feature.qualifiers:
                                extracted_ids.append(feature.qualifiers["locus_tag"][0])
                            elif "protein_id" in feature.qualifiers:
                                extracted_ids.append(feature.qualifiers["protein_id"][0])
                            elif "gene" in feature.qualifiers:
                                extracted_ids.append(feature.qualifiers["gene"][0])
            else:
                QMessageBox.warning(self, "Unsupported Format", f"Cannot extract IDs from .{ext} files.")
                return
            
            if not extracted_ids:
                QMessageBox.information(self, "Info", "No gene IDs or locus_tags found in the selected file.")
                return
                
            seen = set()
            unique_ids = []
            for x in extracted_ids:
                if x not in seen:
                    unique_ids.append(x)
                    seen.add(x)
                    
            current_text = target_text_edit.toPlainText().strip()
            if current_text:
                target_text_edit.setText(current_text + "\n" + "\n".join(unique_ids))
            else:
                target_text_edit.setText("\n".join(unique_ids))
                
            self._write_to_console(f"Loaded {len(unique_ids)} unique IDs from {os.path.basename(filepath)}.\n", "success")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Could not read file:\n{e}")

    def _create_synth_bio_tab(self, tab):
        layout = QHBoxLayout(tab)
        
        left_panel = QWidget()
        lp_layout = QVBoxLayout(left_panel)
        
        lbl_host = QLabel("1. Select Host:")
        lbl_host.setFont(QFont("Calibri", 11, QFont.Weight.Bold))
        lp_layout.addWidget(lbl_host)
        
        self.synth_host_combo = QComboBox()
        lp_layout.addWidget(self.synth_host_combo)
        
        lbl_seq = QLabel("2. Paste CDS Sequence or Load FASTA:")
        lbl_seq.setFont(QFont("Calibri", 11, QFont.Weight.Bold))
        lp_layout.addWidget(lbl_seq)
        
        btn_fasta = QPushButton("📂 Load FASTA")
        btn_fasta.clicked.connect(self._on_load_fasta_synth)
        lp_layout.addWidget(btn_fasta)
        
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

    def _on_load_fasta_synth(self):
        filepath, _ = QFileDialog.getOpenFileName(self, "Select FASTA file", "", "FASTA Files (*.fasta *.fas *.fa);;All Files (*.*)")
        if filepath:
            try:
                record = next(SeqIO.parse(filepath, "fasta"))
                self.synth_input_seq.setText(str(record.seq))
                self._write_to_console(f"FASTA loaded successfully: {filepath}\n", "success")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Could not read FASTA:\n{e}")

    def _create_visualization_tab(self, tab):
        layout = QVBoxLayout(tab)
        layout.setContentsMargins(0, 0, 0, 0)
        self._setup_gallery_tab(layout)

    def _setup_gallery_tab(self, layout):
        toolbar_layout = QHBoxLayout()
        
        self.btn_prev_chart = QPushButton("◀ Previous")
        self.btn_prev_chart.setFixedWidth(100)
        self.btn_prev_chart.clicked.connect(self._on_prev_chart)
        
        self.image_status_label = QLabel("No chart loaded.")
        self.image_status_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.image_status_label.setStyleSheet("color: #c0caf5; font-weight: bold; font-size: 14px;")
        
        self.btn_next_chart = QPushButton("Next ▶")
        self.btn_next_chart.setFixedWidth(100)
        self.btn_next_chart.clicked.connect(self._on_next_chart)

        self.btn_zoom_in = QPushButton("🔍 +")
        self.btn_zoom_in.setFixedWidth(50)
        self.btn_zoom_in.clicked.connect(lambda: self.view.scale(1.2, 1.2))

        self.btn_zoom_out = QPushButton("🔍 -")
        self.btn_zoom_out.setFixedWidth(50)
        self.btn_zoom_out.clicked.connect(lambda: self.view.scale(0.8, 0.8))

        self.btn_zoom_fit = QPushButton("Fit Window")
        self.btn_zoom_fit.clicked.connect(self._fit_image)
        
        self.btn_open_external = QPushButton("🖵 Open Fullscreen")
        self.btn_open_external.clicked.connect(self._on_open_external)
        
        toolbar_layout.addWidget(self.btn_prev_chart)
        toolbar_layout.addWidget(self.image_status_label, stretch=1)
        toolbar_layout.addWidget(self.btn_next_chart)
        
        line = QFrame()
        line.setFrameShape(QFrame.Shape.VLine)
        toolbar_layout.addWidget(line)
        
        toolbar_layout.addWidget(self.btn_zoom_in)
        toolbar_layout.addWidget(self.btn_zoom_out)
        toolbar_layout.addWidget(self.btn_zoom_fit)
        
        line2 = QFrame()
        line2.setFrameShape(QFrame.Shape.VLine)
        toolbar_layout.addWidget(line2)

        toolbar_layout.addWidget(self.btn_open_external)
        
        layout.addLayout(toolbar_layout)

        self.scene = QGraphicsScene()
        self.view = QGraphicsView(self.scene)
        self.view.setRenderHint(QPainter.RenderHint.Antialiasing)
        self.view.setRenderHint(QPainter.RenderHint.SmoothPixmapTransform)
        self.view.setStyleSheet("background-color: #1a1b26; border: 1px solid #3b4261;")
        layout.addWidget(self.view, stretch=1)

        self._update_chart_buttons()

    def _fit_image(self):
        if self.scene.items():
            self.view.fitInView(self.scene.sceneRect(), Qt.AspectRatioMode.KeepAspectRatio)

    def _on_prev_chart(self):
        if self.current_image_index > 0:
            self.current_image_index -= 1
            self._render_current_image()

    def _on_next_chart(self):
        if self.current_image_index < len(self.image_history) - 1:
            self.current_image_index += 1
            self._render_current_image()

    def _on_open_external(self):
        if 0 <= self.current_image_index < len(self.image_history):
            filepath, _ = self.image_history[self.current_image_index]
            try:
                if platform.system() == 'Darwin':
                    subprocess.call(('open', filepath))
                elif platform.system() == 'Windows':
                    os.startfile(filepath)
                else:
                    subprocess.call(('xdg-open', filepath))
            except Exception as e:
                self._write_to_console(f"❌ Error opening external viewer: {e}\n", "error")

    def _update_chart_buttons(self):
        self.btn_prev_chart.setEnabled(self.current_image_index > 0)
        self.btn_next_chart.setEnabled(self.current_image_index < len(self.image_history) - 1)
        self.btn_open_external.setEnabled(self.current_image_index >= 0)
        self.btn_zoom_in.setEnabled(self.current_image_index >= 0)
        self.btn_zoom_out.setEnabled(self.current_image_index >= 0)
        self.btn_zoom_fit.setEnabled(self.current_image_index >= 0)

    def _create_help_tab(self, tab):
        layout = QVBoxLayout(tab)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)

        content = QWidget()
        c_layout = QVBoxLayout(content)
        c_layout.setContentsMargins(24, 16, 24, 16)
        c_layout.setSpacing(2)

        def add_text(text, style="normal"):
            lbl = QLabel(text)
            lbl.setWordWrap(True)
            lbl.setAlignment(Qt.AlignmentFlag.AlignJustify)
            if style == "h1":
                lbl.setFont(QFont("Calibri", 22, QFont.Weight.Bold))
                lbl.setStyleSheet("color: #7aa2f7; margin-top: 4px; margin-bottom: 10px;")
                lbl.setAlignment(Qt.AlignmentFlag.AlignLeft)
            elif style == "h2":
                lbl.setFont(QFont("Calibri", 16, QFont.Weight.Bold))
                lbl.setStyleSheet("color: #9ece6a; margin-top: 18px; margin-bottom: 4px; "
                                  "border-bottom: 1px solid #3d4255; padding-bottom: 3px;")
                lbl.setAlignment(Qt.AlignmentFlag.AlignLeft)
            elif style == "h3":
                lbl.setFont(QFont("Calibri", 13, QFont.Weight.Bold))
                lbl.setStyleSheet("color: #2ac3de; margin-top: 10px; margin-bottom: 2px;")
                lbl.setAlignment(Qt.AlignmentFlag.AlignLeft)
            elif style == "tag":
                lbl.setFont(QFont("Calibri", 12))
                lbl.setStyleSheet("color: #e0af68; background-color: #1e2030; "
                                  "padding: 1px 6px; border-radius: 3px; margin-bottom: 2px;")
                lbl.setAlignment(Qt.AlignmentFlag.AlignLeft)
            elif style == "code":
                lbl.setFont(QFont("Consolas", 12))
                lbl.setStyleSheet("color: #a9ff68; background-color: #1a1b26; "
                                  "padding: 4px 8px; border-radius: 4px; margin: 4px 0;")
                lbl.setAlignment(Qt.AlignmentFlag.AlignLeft)
            else:
                lbl.setFont(QFont("Calibri", 13))
                lbl.setStyleSheet("color: #a9b1d6; margin-bottom: 4px;")
            c_layout.addWidget(lbl)

        def add_separator():
            line = QFrame()
            line.setFrameShape(QFrame.Shape.HLine)
            line.setStyleSheet("color: #3d4255; margin: 6px 0;")
            c_layout.addWidget(line)

        # ── Header ────────────────────────────────────────────────────────────
        add_text("Kodon-X  —  User Guide", "h1")
        add_text(
            "Kodon-X is an integrated platform for Codon Usage Bias (CUB) analysis of prokaryotic "
            "and eukaryotic genomes. All analyses accept GenBank files (.gbk / .gb / .gbff) as input "
            "and write figures (PNG, 150 dpi) and tables (CSV, semicolon-separated) to the chosen "
            "output folder.",
            "normal"
        )

        # ── Quick workflow ─────────────────────────────────────────────────────
        add_separator()
        add_text("QUICK WORKFLOW", "h2")
        add_text("Step 1 — Load files", "h3")
        add_text(
            "Click Browse… in the Step 1 panel and select the folder containing your .gbk / .gb files. "
            "Then click Load Files. Each file is listed with a checkbox — uncheck files you want to exclude "
            "from the current run.",
            "normal"
        )
        add_text("Step 2 — Configure options", "h3")
        add_text(
            "Select an analysis from the list. The panel below the list updates automatically to show "
            "analysis-specific options (genetic code, palette, optional parameters). "
            "A global Gene Filter (locus_tag, one per line) restricts all analyses to the listed genes.",
            "normal"
        )
        add_text("Step 3 — Run and navigate results", "h3")
        add_text(
            "Click Run Analysis. Progress is shown in the console. Figures appear in the Output tab "
            "as soon as each one is ready. Use the ← → buttons to navigate between figures, "
            "or Open External to view the file in your default image viewer.",
            "normal"
        )

        # ── Input format ──────────────────────────────────────────────────────
        add_separator()
        add_text("INPUT FORMAT", "h2")
        add_text(
            "All analyses require annotated GenBank files. Each CDS feature must have a /translation "
            "qualifier (protein sequence) and a /locus_tag qualifier (unique gene identifier). "
            "Files generated by Prokka, RAST, or downloaded from NCBI RefSeq are fully compatible. "
            "FASTA files are not supported for bias analyses.",
            "normal"
        )
        add_text("Supported extensions:", "tag")
        add_text("   .gbk   .gb   .gbff   .gbf", "code")

        # ── Synthetic biology ──────────────────────────────────────────────────
        add_separator()
        add_text("SYNTHETIC BIOLOGY MODULE", "h2")
        add_text(
            "The Synthetic Biology tab adapts a heterologous DNA sequence for optimal expression in a "
            "selected host organism, using the host's own codon usage statistics derived from its "
            "GenBank file. Two modes are available:",
            "normal"
        )
        add_text("Optimize", "h3")
        add_text(
            "Each codon is replaced by the synonymous codon with the highest RSCU in the host genome. "
            "This maximises the use of the most abundant tRNAs and typically increases recombinant "
            "protein yield (Gustafsson et al. 2004).",
            "normal"
        )
        add_text("Harmonize", "h3")
        add_text(
            "The codon usage profile of the source sequence is mapped onto the host: if a codon is "
            "in the n-th percentile of usage frequency in the source, it is replaced by the codon in "
            "the same percentile in the host. This preserves the original translation rhythm — "
            "ribosome pauses associated with rare codons that are important for correct co-translational "
            "protein folding are maintained (Angov et al. 2008).",
            "normal"
        )
        add_text(
            "Output: an optimized / harmonized FASTA file and a comparison table (original vs. new "
            "codons) are saved to the output folder.",
            "normal"
        )

        # ── Analysis descriptions ──────────────────────────────────────────────
        add_separator()
        add_text("ANALYSIS REFERENCE", "h2")
        add_text(
            "Each analysis is numbered. Select it by ID in the Run Analysis tab. "
            "Required file count is shown in brackets.",
            "normal"
        )

        for analysis_name, data in self.analysis_data.items():
            add_text(analysis_name, "h3")
            add_text(f"Files required: {data['files_required']}", "tag")
            add_text(data['description'], "normal")

        # ── Tips ──────────────────────────────────────────────────────────────
        add_separator()
        add_text("TIPS & TROUBLESHOOTING", "h2")
        add_text("Genetic code", "h3")
        add_text(
            "The Standard code (table 1) is the default. Change it via the Genetic Code selector "
            "if analysing mitochondrial, mycoplasma, or other non-standard genomes. "
            "The code affects RSCU, ENC, CAI, tAI, CPB, GRAVY, and group comparison analyses.",
            "normal"
        )
        add_text("Gene filter", "h3")
        add_text(
            "Paste locus_tags (one per line) in the Gene Filter box to restrict analyses to a subset "
            "of genes — e.g. ribosomal proteins for a reference set, or virulence genes of interest. "
            "Leave it empty to use all CDS features.",
            "normal"
        )
        add_text("Color palette", "h3")
        add_text(
            "All figures use the selected palette. Recommended options: viridis, plasma (sequential); "
            "RdBu_r, coolwarm (diverging); Set2, Dark2, tab10 (categorical). "
            "The dinucleotide ρ_XY heatmap always uses RdBu_r regardless of this setting.",
            "normal"
        )
        add_text("Large genomes", "h3")
        add_text(
            "Analyses 4, 6, 11, and 12 scale figure height automatically with the number of genomes. "
            "For more than 20 genomes, consider splitting runs into smaller batches for readability.",
            "normal"
        )

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
        
        patterns = ["*.gbk", "*.gb", "*.gbff", "*.gbf"]
        files = []
        for pattern in patterns:
            files.extend(glob.glob(os.path.join(folder, pattern)))
            
        if not files:
            self._write_to_console(f"No .gbk, .gb, .gbff or .gbf files found in '{folder}'\n", "warning")
            self.lbl_file_count.setText("0 files")
            return
            
        self._write_to_console(f"Found {len(files)} valid GenBank files.\n", "success")
        
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

        elif data['id'] == '14':
            w = QWidget()
            lo = QHBoxLayout(w)
            lo.addWidget(QLabel("Superkingdom (Wobble Rules):"))
            self.super_kingdom_combo = QComboBox()
            self.super_kingdom_combo.addItems(["Bacteria", "Eukaryote"])
            lo.addWidget(self.super_kingdom_combo)
            lo.addStretch()
            self.analysis_options_layout.addWidget(w)

        elif data['id'] == '15':
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

        elif data['id'] == '16':
            w = QWidget()
            lo = QHBoxLayout(w)
            lo.addWidget(QLabel("5' Region (bp):"))
            self.mfe_region_edit = QLineEdit("50")
            self.mfe_region_edit.setFixedWidth(60)
            lo.addWidget(self.mfe_region_edit)
            lo.addStretch()
            self.analysis_options_layout.addWidget(w)

        elif data['id'] == '17':
            self.gene_filter_text.setEnabled(False)
            w = QWidget()
            lo = QVBoxLayout(w)
            
            h1 = QHBoxLayout()
            h1.addWidget(QLabel("Gene Group 1 (one per line):"))
            btn_g1 = QPushButton("📂 Load File")
            btn_g1.clicked.connect(lambda: self._on_load_gene_list_file(self.gene_list_1_text))
            h1.addWidget(btn_g1)
            h1.addStretch()
            lo.addLayout(h1)
            self.gene_list_1_text = QTextEdit()
            self.gene_list_1_text.setMaximumHeight(80)
            lo.addWidget(self.gene_list_1_text)
            
            h2 = QHBoxLayout()
            h2.addWidget(QLabel("Gene Group 2 (one per line):"))
            btn_g2 = QPushButton("📂 Load File")
            btn_g2.clicked.connect(lambda: self._on_load_gene_list_file(self.gene_list_2_text))
            h2.addWidget(btn_g2)
            h2.addStretch()
            lo.addLayout(h2)
            self.gene_list_2_text = QTextEdit()
            self.gene_list_2_text.setMaximumHeight(80)
            lo.addWidget(self.gene_list_2_text)
            
            self.analysis_options_layout.addWidget(w)

        elif data['id'] == '18':
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
        
        extra_args['palette'] = self.palette_combo.currentText().split()[0]
        
        if analysis_data.get('id') != '17':
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
            elif analysis_data.get('id') == '14':
                extra_args['super_kingdom'] = self.super_kingdom_combo.currentText()
            elif analysis_data.get('id') == '15':
                extra_args['upstream_dist'] = int(self.upstream_dist_edit.text())
                extra_args['kmer_size'] = int(self.kmer_size_edit.text())
            elif analysis_data.get('id') == '16':
                extra_args['mfe_region_length'] = int(self.mfe_region_edit.text())
            elif analysis_data.get('id') == '17':
                g1_text = self.gene_list_1_text.toPlainText()
                g2_text = self.gene_list_2_text.toPlainText()
                gene_list_1 = set(line.strip() for line in g1_text.splitlines() if line.strip())
                gene_list_2 = set(line.strip() for line in g2_text.splitlines() if line.strip())
                if not gene_list_1 or not gene_list_2:
                    QMessageBox.critical(self, "Error", "Analysis 17 requires BOTH groups to be filled.")
                    return
                extra_args['gene_list_1'] = gene_list_1
                extra_args['gene_list_2'] = gene_list_2
                gene_list = None
            elif analysis_data.get('id') == '18':
                gene_col = self.expr_gene_col_edit.text()
                expr_col = self.expr_value_col_edit.text()
                if self.expression_dataframe is None:
                    QMessageBox.critical(self, "Error", "Load an expression file for Analysis 18.")
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
            print(f"Error: Host file '{host_filename}' not found.")
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
            QMessageBox.warning(self, "Warning", "Insert a valid DNA sequence (minimum 10bp).")
            return
            
        output_filepath, _ = QFileDialog.getSaveFileName(self, "Save Result as GenBank", f"result_{mode}.gbk", "GenBank Files (*.gbk *.gb)")
        if not output_filepath:
            return 
            
        self.dock_console.show() 
        self.run_optimize_btn.setEnabled(False)
        self.run_harmonize_btn.setEnabled(False)
        self.progress_bar.show()
        self.progress_bar.setRange(0, 0)
        self.status_label.setText(f"Starting {mode}...")
        
        genetic_code_id = self._get_genetic_code_id()

        thread = threading.Thread(
            target=self._optimization_thread_target,
            args=(input_seq, host_filename, genetic_code_id, mode, output_filepath),
            daemon=True
        )
        thread.start()

    def _optimization_thread_target(self, input_seq, host_filename, genetic_code_id, mode, output_filepath):
        try:
            print(f"\n{'='*60}")
            print(f"STARTING SYNTHETIC BIOLOGY TOOL")
            print(f"   Mode: {mode}")
            print(f"   Host: {host_filename}")
            print(f"   Input Size: {len(input_seq)}bp")
            print(f"{'='*60}")
            
            host_data = self._get_host_data(host_filename)
            if not host_data:
                raise Exception("Could not retrieve host data.")
            
            self.status_queue.put(("progress", 50))
            
            result_seq = ""
            if mode == "optimize":
                self.status_queue.put(("message", "Optimizing sequence..."))
                result_seq = optimize_codon_sequence(input_seq, host_data['rscu'], genetic_code_id)
            elif mode == "harmonize":
                self.status_queue.put(("message", "Harmonizing sequence..."))
                result_seq = harmonize_codon_sequence(input_seq, host_data['counts'], genetic_code_id)
            
            record = SeqRecord(Seq(result_seq), id=f"Synth_{mode}", name="SyntheticBio", description=f"Sequence {mode}d for {host_filename}")
            
            record.annotations["molecule_type"] = "DNA"
            
            feature = SeqFeature(FeatureLocation(start=0, end=len(result_seq)), type="CDS")
            record.features.append(feature)
            
            with open(output_filepath, "w") as out_handle:
                SeqIO.write(record, out_handle, "genbank")
                
            print(f"\n✅ GenBank file saved successfully at:\n{output_filepath}")
            self.status_queue.put(("message", f"GenBank saved at {os.path.basename(output_filepath)}"))
            self.status_queue.put(("optimization_complete", result_seq))
            print("✅ Process completed successfully.")
        except Exception as e:
            print(f"\n❌ ERROR DURING OPERATION: {e}", "error")
            import traceback
            print(traceback.format_exc())
            self.status_queue.put(("done", None)) 
        finally:
            self.status_queue.put(("tool_done", None))

    def _start_analysis_thread(self, files, output_folder, analysis_name, gene_list, extra_args):
        self.dock_console.show() 
        self.console_text.clear()
        
        self.run_button.setEnabled(False)
        self.progress_bar.show()
        self.progress_bar.setRange(0, 0)
        self.status_label.setText(f"Starting: {analysis_name}...")
        
        genetic_code_id = self._get_genetic_code_id()
        analysis_data = self.analysis_data[analysis_name]
        
        self.image_history.clear()
        self.current_image_index = -1
        self._update_chart_buttons()
        self.scene.clear()
        self.image_status_label.setText("No chart loaded.")
        
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
            print(f"STARTING ANALYSIS: {analysis_name}")
            if gene_list:
                print(f"FILTER ACTIVE: Analyzing only {len(gene_list)} specified genes.")
            print(f"{'='*60}")
            print(f"Files: {len(files)}")
            print(f"Output Folder: {output_folder}")
            print(f"Genetic Table: {genetic_code_id}")
            print(f"Palette: {extra_args.get('palette', 'viridis')}")
            print(f"{'='*60}")
            
            palette = extra_args.get('palette', 'viridis')
            
            if analysis_data['id'] in ['3', '4', '5', '6', '8', '11']:
                all_bias_data = process_genomes_for_bias_analysis(files, genetic_code_id, self.status_queue, gene_list)
                if not all_bias_data:
                    raise Exception("Failed to process bias data. Check input files.")
                analysis_function = analysis_data['function']
                if analysis_function:
                    if analysis_data['id'] == '3':
                        analysis_function(all_bias_data, output_folder, genetic_code_id, self.status_queue, palette=palette)
                    elif analysis_data['id'] in ['6', '11']:
                        analysis_function(all_bias_data, output_folder, self.status_queue,
                                          file_list=files, genetic_code_id=genetic_code_id,
                                          gene_list=gene_list, palette=palette)
                    else:
                        analysis_function(all_bias_data, output_folder, self.status_queue, palette=palette)
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
            elif analysis_data['id'] == '7':
                analyze_genomic_composition(files, output_folder, self.status_queue, palette=palette)
            elif analysis_data['id'] == '9':
                codon_pair_bias_analysis(files, output_folder, genetic_code_id, self.status_queue, gene_list, palette=palette)
            elif analysis_data['id'] == '10':
                gravy_aromo_analysis(files, output_folder, genetic_code_id, self.status_queue, gene_list, palette=palette)
            elif analysis_data['id'] == '12':
                dinucleotide_composition_analysis(files, output_folder, self.status_queue, palette=palette)
            elif analysis_data['id'] == '13':
                pr2_plot_analysis(files, output_folder, genetic_code_id, self.status_queue, gene_list, palette=palette)
            elif analysis_data['id'] == '14':
                kingdom_name = extra_args.get('super_kingdom', 'Bacteria')
                tai_analysis(files, output_folder, genetic_code_id, self.status_queue, gene_list,
                             super_kingdom=kingdom_name, palette=palette)
            elif analysis_data['id'] == '15':
                upstream_motifs_analysis(files, output_folder, self.status_queue, gene_list,
                                         extra_args['upstream_dist'], extra_args['kmer_size'], palette=palette)
            elif analysis_data['id'] == '16':
                initiation_mfe_analysis(files, output_folder, genetic_code_id, self.status_queue, gene_list,
                                        extra_args['mfe_region_length'], palette=palette)
            elif analysis_data['id'] == '17':
                two_groups_comparative_analysis(files, output_folder, genetic_code_id, self.status_queue,
                                                extra_args['gene_list_1'], extra_args['gene_list_2'], palette=palette)
            elif analysis_data['id'] == '18':
                expression_correlation_analysis(files, output_folder, genetic_code_id, self.status_queue,
                                                gene_list, extra_args['expression_data'],
                                                extra_args['gene_col'], extra_args['expr_col'], palette=palette)
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
        updates_made = False
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
                updates_made = True
        except queue.Empty:
            pass
            
        if updates_made:
            self.console_text.ensureCursorVisible()
        
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

    def _display_image(self, image_path, title):
        try:
            self.image_history.append((image_path, title))
            self.current_image_index = len(self.image_history) - 1
            self._render_current_image()
            self.notebook.setCurrentWidget(self.tab_viz)
        except Exception as e:
            self._write_to_console(f"❌ Error displaying image {image_path}: {e}\n", "error")

    def _render_current_image(self):
        if 0 <= self.current_image_index < len(self.image_history):
            image_path, title = self.image_history[self.current_image_index]
            self.scene.clear()
            pixmap = QPixmap(image_path)
            self.scene.addPixmap(pixmap)
            self._fit_image()
            self.image_status_label.setText(f"({self.current_image_index + 1}/{len(self.image_history)}) {title}")
        self._update_chart_buttons()

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
