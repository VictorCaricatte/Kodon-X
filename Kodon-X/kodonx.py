import sys

try:
    from PyQt6.QtWidgets import QApplication
except ImportError:
    print("Error: The 'PyQt6' library is required. Install with: pip install PyQt6")
    sys.exit(1)

try:
    from PIL import Image
except ImportError:
    print("Error: The 'Pillow' library is required. Install with: pip install Pillow")
    sys.exit(1)

try:
    import pandas as pd
except ImportError:
    print("Error: The 'pandas' library is required. Install with: pip install pandas")
    sys.exit(1)

try:
    import matplotlib
    import matplotlib.pyplot as plt
    import seaborn as sns
    from matplotlib.figure import Figure
    
    matplotlib.use('Agg')
except ImportError:
    print("Error: 'matplotlib' and 'seaborn' are required. Install with: pip install matplotlib seaborn")
    sys.exit(1)

try:
    import numpy as np
except ImportError:
    print("Error: 'numpy' is required. Install with: pip install numpy")
    sys.exit(1)

try:
    from scipy.cluster import hierarchy as sch
    from scipy.stats import pearsonr, chi2_contingency
    from scipy.spatial.distance import pdist, squareform
except ImportError:
    print("Error: 'scipy' is required. Install with: pip install scipy")
    sys.exit(1)

try:
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
    from sklearn.manifold import MDS
except ImportError:
    print("Error: 'scikit-learn' is required. Install with: pip install scikit-learn")
    sys.exit(1)

try:
    from Bio import SeqIO
except ImportError:
    print("Error: 'biopython' is required. Install with: pip install biopython")
    sys.exit(1)

from interface import KodonE_GUI

try:
    import RNA
except ImportError:
    print("-------------------------------------------------------------------")
    print("Warning: The 'RNA' (ViennaRNA) library was not found.")
    print("MFE Analysis (Analysis 17) will be disabled.")
    print("To enable it, install via Conda: 'conda install -c bioconda viennarna'")
    print("Or install the 'vienna-rna' package and the 'rna-wrapper'.")
    print("-------------------------------------------------------------------")

def main():
    app = QApplication(sys.argv)
    app.setStyle("Fusion") 
    window = KodonE_GUI()
    window.show()
    sys.exit(app.exec())

if __name__ == "__main__":
    main()
