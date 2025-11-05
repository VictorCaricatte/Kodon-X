# main.py
import sys
import tkinter as tk

# --- Verificação e Importação de Bibliotecas ---
# Isso garante que o ambiente seja verificado antes de tentar iniciar a GUI.
try:
    from PIL import Image, ImageTk
except ImportError:
    print("Erro: A biblioteca 'Pillow' é necessária. Instale com: pip install Pillow")
    sys.exit(1)

try:
    import pandas as pd
except ImportError:
    print("Erro: A biblioteca 'pandas' é necessária. Instale com: pip install pandas")
    sys.exit(1)

try:
    import matplotlib
    import matplotlib.pyplot as plt
    import seaborn as sns
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
    from matplotlib.figure import Figure
    matplotlib.use('Agg')
except ImportError:
    print("Erro: 'matplotlib' e 'seaborn' são necessários. Instale com: pip install matplotlib seaborn")
    sys.exit(1)

try:
    import numpy as np
except ImportError:
    print("Erro: 'numpy' é necessário. Instale com: pip install numpy")
    sys.exit(1)

try:
    from scipy.cluster import hierarchy as sch
    from scipy.stats import pearsonr, chi2_contingency
    from scipy.spatial.distance import pdist, squareform
except ImportError:
    print("Erro: 'scipy' é necessário. Instale com: pip install scipy")
    sys.exit(1)

try:
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
    from sklearn.manifold import MDS
except ImportError:
    print("Erro: 'scikit-learn' é necessário. Instale com: pip install scikit-learn")
    sys.exit(1)

try:
    from Bio import SeqIO
except ImportError:
    print("Erro: 'biopython' é necessário. Instale com: pip install biopython")
    sys.exit(1)

# Importa a classe da GUI do nosso arquivo
from gui import KodonE_GUI

try:
    import RNA
except ImportError:
    print("-------------------------------------------------------------------")
    print("Aviso: A biblioteca 'RNA' (ViennaRNA) não foi encontrada.")
    print("A Análise de MFE (Análise 17) estará desabilitada.")
    print("Para habilitar, instale via Conda: 'conda install -c bioconda viennarna'")
    print("Ou instale o pacote 'vienna-rna' e o wrapper 'rna-wrapper'.")
    print("-------------------------------------------------------------------")
# ######################################################################
# --- PONTO DE ENTRADA PRINCIPAL ---
# ######################################################################

def main():
    """Função principal para iniciar a aplicação"""
    root = tk.Tk()
    
    # Criar aplicação
    app = KodonE_GUI(root)
    
    # Iniciar loop principal
    root.mainloop()

if __name__ == "__main__":
    main()
