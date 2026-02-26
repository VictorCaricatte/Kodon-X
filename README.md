# KodonX — Comprehensive Codon Usage & Genomic Analysis Tool

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**KodonX** is a comprehensive bioinformatics platform for codon usage bias (CUB) analysis, comparative genomics, and synthetic biology. It provides a graphical user interface (PyQt6) and a fully functional command-line interface (CLI), enabling researchers to extract, calculate, and visualize a broad set of codon-level metrics from annotated GenBank (`.gbk`/`.gb`/`.gbff`) files.

---

## Table of Contents

- [Features Overview](#features-overview)
- [Installation](#installation)
- [Input Format](#input-format)
- [Graphical Interface (GUI)](#graphical-interface-gui)
- [Command-Line Interface (CLI)](#command-line-interface-cli)
- [Supported Genetic Codes](#supported-genetic-codes)
- [Analysis Modules](#analysis-modules)
  - [Analysis 0 — Synthetic Biology](#analysis-0--synthetic-biology)
  - [Analysis 1 — Statistics and CDS](#analysis-1--statistics-and-cds)
  - [Analysis 2 — Gene Listing](#analysis-2--gene-listing)
  - [Analysis 3 — Individual RSCU Heatmap](#analysis-3--individual-rscu-heatmap)
  - [Analysis 4 — Comparative RSCU (Clustermap + PCA)](#analysis-4--comparative-rscu-clustermap--pca)
  - [Analysis 5 — RSCU Correlation (2 Genomes)](#analysis-5--rscu-correlation-2-genomes)
  - [Analysis 6 — Comparative RSCU Histograms](#analysis-6--comparative-rscu-histograms)
  - [Analysis 7 — ENC vs GC3 (Wright Plot)](#analysis-7--enc-vs-gc3-wright-plot)
  - [Analysis 8 — Genomic Composition](#analysis-8--genomic-composition)
  - [Analysis 9 — Optimal, Rare Codons & CAI](#analysis-9--optimal-rare-codons--cai)
  - [Analysis 10 — Codon Pair Bias (CPB)](#analysis-10--codon-pair-bias-cpb)
  - [Analysis 11 — Physicochemical Analysis (GRAVY & Aromo)](#analysis-11--physicochemical-analysis-gravy--aromo)
  - [Analysis 12 — Neutrality Plot (GC12 vs GC3)](#analysis-12--neutrality-plot-gc12-vs-gc3)
  - [Analysis 13 — Dinucleotide Composition](#analysis-13--dinucleotide-composition)
  - [Analysis 14 — PR2 Parity Plot](#analysis-14--pr2-parity-plot)
  - [Analysis 15 — tRNA Adaptation Index (tAI)](#analysis-15--trna-adaptation-index-tai)
  - [Analysis 16 — Upstream Motifs Analysis](#analysis-16--upstream-motifs-analysis)
  - [Analysis 17 — MFE Analysis (5' Structure)](#analysis-17--mfe-analysis-5-structure)
  - [Analysis 18 — Gene Group Comparison](#analysis-18--gene-group-comparison)
  - [Analysis 19 — Correlation with Expression](#analysis-19--correlation-with-expression)
- [Project Structure](#project-structure)
- [References](#references)

---

## Features Overview

- **19 analysis modules** covering codon bias, comparative genomics, translation efficiency, and gene expression.
- **Synthetic Biology module** for sequence optimization and harmonization.
- Supports multiple **genetic code tables** (Standard, Vertebrate Mitochondrial, Bacterial/Plant Plastid, Mold/Protozoan).
- Optional **gene-level filtering** via locus tag lists.
- Dual interface: **PyQt6 GUI** and **CLI** for use in automated pipelines.
- All outputs saved as publication-ready **PNG figures** and **CSV tables**.

---

## Installation

### Prerequisites

```bash
pip install PyQt6 Pillow pandas matplotlib seaborn numpy scipy scikit-learn biopython
```

### Optional (for MFE Analysis — Analysis 17)

ViennaRNA must be installed via Conda:

```bash
conda install -c bioconda viennarna
```

### Running

**GUI:**
```bash
python kodonx.py
```

**CLI:**
```bash
python kodonargs.py -i /path/to/gbk_folder -o /path/to/output -a <analysis_id>
```

---

## Input Format

KodonX accepts annotated genome files in **GenBank flat-file format** (`.gbk`, `.gb`, `.gbff`), as provided by NCBI. Files must contain `CDS` features with at least a `locus_tag` or `gene` qualifier, and optionally `tRNA` features (for Analysis 15).

Multiple files can be provided in a single input directory, enabling comparative analyses across multiple genomes or strains.

### Gene Filter File

A plain-text file with one `locus_tag` or gene name per line can be passed with `-f` to restrict all analyses to a specific gene subset:

```
b0001
b0002
recA
rpoB
```

---

## Graphical Interface (GUI)

Launch with `python kodonx.py`. The GUI provides:

- File selection with multi-genome support.
- Dropdown menus for genetic code and analysis selection.
- Real-time progress tracking and status messages.
- Integrated image viewer for generated plots.
- Gene filter input via text field or file.

---

## Command-Line Interface (CLI)

```
python kodonargs.py -i <input_dir> -o <output_dir> -a <analysis_id> [options]

Options:
  -g, --genetic-code INT    Genetic code table ID (default: 1)
  -f, --filter-file FILE    Text file with locus_tags for filtering
  --start-codon STR         Start codon filter for Analysis 1 (default: ATG)
  --superkingdom STR        Bacteria or Eukaryote, for tAI wobble rules (Analysis 15)
  --upstream-dist INT       Upstream distance in bp (Analysis 16, default: 200)
  --kmer-size INT           K-mer size for motif search (Analysis 16, default: 6)
  --mfe-region INT          5' region size for MFE (Analysis 17, default: 50)
  --group1 FILE             Gene list for Group 1 (Analysis 18)
  --group2 FILE             Gene list for Group 2 (Analysis 18)
  --expr-file FILE          Expression data CSV/TSV (Analysis 19)
  --expr-gene STR           Gene column name in expression file (default: locus_tag)
  --expr-val STR            Expression value column name (default: TPM)
  --synth-mode STR          optimize or harmonize (Analysis 0)
  --synth-host STR          Host genome filename inside input dir (Analysis 0)
  --synth-seq STR           Input DNA sequence (Analysis 0)
```

---

## Supported Genetic Codes

| ID | Description |
|----|-------------|
| 1  | Standard (NCBI Table 1) |
| 2  | Vertebrate Mitochondrial (NCBI Table 2) |
| 4  | Mold, Protozoan, Coelenterate Mitochondrial (NCBI Table 4) |
| 11 | Bacterial, Archaeal & Plant Plastid (NCBI Table 11) |

Genetic code definitions follow NCBI: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi

---

## Analysis Modules

---

### Analysis 0 — Synthetic Biology

**Mode:** `optimize` or `harmonize`  
**Input:** Any single host `.gbk` file + a raw DNA sequence string  
**Output:** Recoded DNA sequence (printed to terminal / GUI)

This module recodes an input CDS for optimal expression in a target host organism. Two strategies are implemented:

#### Codon Optimization (Maximization)

The `optimize` mode replaces each codon in the input sequence with the **most frequent synonymous codon** in the host genome — i.e., the codon with the highest RSCU value for each amino acid. This maximizes the use of highly expressed codons and is the classical approach used in recombinant protein production.

**Methodology:** For each amino acid, the optimal codon is identified as `argmax(RSCU_i)` across all synonymous codons. RSCU is computed from the host genome as described in Analysis 3.

**References:**
- Sharp, P.M. & Li, W.H. (1987). The codon adaptation index — a measure of directional synonymous codon usage bias, and its potential applications. *Nucleic Acids Research*, 15(3), 1281–1295. https://doi.org/10.1093/nar/15.3.1281
- Gustafsson, C., et al. (2004). Codon bias and heterologous protein expression. *Trends in Biotechnology*, 22(7), 346–353. https://doi.org/10.1016/j.tibtech.2004.04.006
- Puigbò, P., et al. (2007). OPTIMIZER: A web server for optimizing the codon usage of DNA sequences. *Nucleic Acids Research*, 35(Web Server issue), W126–W131. https://doi.org/10.1093/nar/gkm219

#### Codon Harmonization (Rank-Based)

The `harmonize` mode preserves the **relative codon usage pattern** of the input sequence, remapping codons to their rank-equivalent in the host genome. Instead of always maximizing, this approach maintains rare codons where they appear in the original sequence, which has been shown to improve proper folding of complex proteins by preserving translational pausing.

**Methodology:** For each amino acid, synonymous codons in the input sequence are ranked by frequency; the corresponding rank positions in the host codon frequency table are used as replacements.

**References:**
- Angov, E., et al. (2008). Heterologous protein expression is enhanced by harmonizing the codon usage frequencies of the target gene with those of the expression host. *PLoS ONE*, 3(5), e2189. https://doi.org/10.1371/journal.pone.0002189
- Chaney, J.L. & Clark, P.L. (2015). Roles for synonymous codon usage in protein biogenesis. *Annual Review of Biophysics*, 44, 143–166. https://doi.org/10.1146/annurev-biophys-060414-034333

---

### Analysis 1 — Statistics and CDS

**Input:** One or more `.gbk` files  
**Output:** `genome_statistics.csv`, `cds_analysis.csv`

Performs an initial scan of all input genomes, reporting total number of contigs, total genome length, and global GC content (%). Additionally, all CDS features are extracted and inspected for start codon usage, ORF size, and consistency with the `codon_start` qualifier. This analysis is the recommended first step before any codon bias analysis.

**References:**
- Cock, P.J.A., et al. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics*, 25(11), 1422–1423. https://doi.org/10.1093/bioinformatics/btp163

---

### Analysis 2 — Gene Listing

**Input:** Exactly 1 `.gbk` file  
**Output:** `gene_list_<filename>.csv`

Enumerates all annotated genomic features (CDS, tRNA, rRNA, gene) in a GenBank file, reporting locus tags, gene names, products, and feature types. Useful for constructing gene filter lists for subsequent analyses.

---

### Analysis 3 — Individual RSCU Heatmap

**Input:** Exactly 1 `.gbk` file  
**Output:** RSCU heatmap (PNG), codon counts table (CSV)  
**Metrics reported:** RSCU, ENC, GC3, CAI

**Relative Synonymous Codon Usage (RSCU)** is the ratio of the observed frequency of a codon to the expected frequency if all synonymous codons were used equally:

```
RSCU_ij = X_ij / (1/n_i * Σ X_ij)
```

Where `X_ij` is the count of codon `j` for amino acid `i`, and `n_i` is the number of synonymous codons. RSCU = 1.0 indicates unbiased usage; values > 1.0 indicate preferred codons; values < 1.0 indicate underused codons.

Global metrics (ENC, GC3, CAI) are also calculated for the genome — see Analyses 7 and 9 for their detailed descriptions.

**References:**
- Sharp, P.M. & Li, W.H. (1987). The codon adaptation index — a measure of directional synonymous codon usage bias, and its potential applications. *Nucleic Acids Research*, 15(3), 1281–1295. https://doi.org/10.1093/nar/15.3.1281

---

### Analysis 4 — Comparative RSCU (Clustermap + PCA)

**Input:** 2 or more `.gbk` files  
**Output:** Comparative RSCU matrix (CSV), clustermap (PNG), PCA plot (PNG)

Compares codon usage profiles across multiple genomes by constructing a matrix of RSCU values (genomes × 64 codons), followed by hierarchical clustering (UPGMA with Euclidean distance) and Principal Component Analysis (PCA).

The **clustermap** reveals groups of organisms with similar codon preferences, while the **PCA** reduces the 64-dimensional RSCU space to two principal components for visualization.

**Methodology:** RSCU matrix is standardized (Z-scores) prior to PCA using `StandardScaler`. PCA is performed using scikit-learn with `n_components=2`. Hierarchical clustering uses the `average` linkage method.

**References:**
- Karlin, S. & Mrazek, J. (1996). What drives codon choices in human genes? *Journal of Molecular Biology*, 262(4), 459–472. https://doi.org/10.1006/jmbi.1996.0528
- Pedregosa, F., et al. (2011). Scikit-learn: Machine learning in Python. *Journal of Machine Learning Research*, 12, 2825–2830.
- Novembre, J. (2002). Accounting for background nucleotide composition when measuring codon usage bias. *Molecular Biology and Evolution*, 19(8), 1390–1394. https://doi.org/10.1093/oxfordjournals.molbev.a004198

---

### Analysis 5 — RSCU Correlation (2 Genomes)

**Input:** Exactly 2 `.gbk` files  
**Output:** Pearson correlation scatter plot (PNG)

Calculates the **Pearson correlation coefficient (R)** between the RSCU vectors of two organisms, treating each of the 64 codons as a variable. A high R indicates similar codon usage preferences between the two genomes, which may reflect evolutionary proximity, shared translational machinery, or similar mutational biases.

**References:**
- Sharp, P.M., et al. (1995). Codon usage in Escherichia coli: genes for proteins with differing cellular abundance. *EMBO Journal*, 5(6), 1233–1244.
- Ikemura, T. (1981). Correlation between the abundance of Escherichia coli transfer RNAs and the occurrence of the respective codons in its protein genes. *Journal of Molecular Biology*, 146(1), 1–21. https://doi.org/10.1016/0022-2836(81)90363-6

---

### Analysis 6 — Comparative RSCU Histograms

**Input:** 2 or more `.gbk` files  
**Output:** Box plots by amino acid (PNG), stacked bar charts by species (PNG)

Generates detailed visual summaries of codon usage profiles. For each species, a **stacked bar chart** shows the RSCU contribution of each codon per amino acid. **Box plots** display the distribution of RSCU values across synonymous codons for each amino acid, enabling visual comparison across genomes.

---

### Analysis 7 — ENC vs GC3 (Wright Plot)

**Input:** 2 or more `.gbk` files  
**Output:** Wright plot (PNG), ENC/GC3 results table (CSV)

The **Wright Plot** (ENC vs GC3) is a classic tool for distinguishing the driving forces behind codon bias. Points falling on or near the theoretical expected curve indicate codon bias primarily driven by **mutational pressure** (GC-biased mutation). Points falling well below the curve indicate that **natural selection** on translational efficiency is an additional force shaping codon usage.

**ENC (Effective Number of Codons)** ranges from 20 (extreme bias: one codon per amino acid) to 61 (no bias: all synonymous codons equally used). It is calculated following Wright (1990):

```
ENC = 2 + 9/F2 + 1/F3 + 5/F4 + 3/F6
```

Where `F_k` is the mean F-value (homozygosity) across all k-fold degenerate amino acid families.

**GC3** is the GC content at the third codon position (synonymous sites), expressed as a percentage.

The theoretical expected ENC curve under purely mutational bias is derived from:

```
ENC_expected = 2 + s² + [(1-s)² / (s² + (1-s)²)] × (remaining terms)
```

Where `s` = GC3/100.

**References:**
- Wright, F. (1990). The 'effective number of codons' used in a gene. *Gene*, 87(1), 23–29. https://doi.org/10.1016/0378-1119(90)90491-9
- Sueoka, N. (1988). Directional mutation pressure and neutral molecular evolution. *PNAS*, 85(8), 2653–2657. https://doi.org/10.1073/pnas.85.8.2653
- dos Reis, M., et al. (2004). Solving the riddle of codon usage preferences: a test for translational selection. *Nucleic Acids Research*, 32(17), 5036–5044. https://doi.org/10.1093/nar/gkh834

---

### Analysis 8 — Genomic Composition

**Input:** One or more `.gbk` files  
**Output:** Nucleotide composition plot (PNG), composition results table (CSV)

Calculates the percentage of each nucleotide (A, T, G, C) and total GC content across entire genomic sequences. Generates a stacked bar chart of nucleotide composition per genome and a scatter plot of genome size versus GC content, which can reveal phylogenetic groupings and evolutionary trends.

**References:**
- Sueoka, N. (1962). On the genetic basis of variation and heterogeneity of DNA base composition. *PNAS*, 48(4), 582–592. https://doi.org/10.1073/pnas.48.4.582

---

### Analysis 9 — Optimal, Rare Codons & CAI

**Input:** One or more `.gbk` files  
**Output:** Optimal/rare codon tables (CSV), CAI bar chart (PNG)

**Optimal codons** are defined as those with RSCU > 1.2 (preferred, above-average usage). **Rare codons** are those with RSCU < 0.8 (underused).

The **Codon Adaptation Index (CAI)** measures how closely a gene's codon usage matches the codon preference of a reference set of highly expressed genes. In KodonX, CAI is calculated genome-wide using RSCU as the reference:

```
w_ij = RSCU_ij / max(RSCU_i)
CAI = exp(1/L × Σ ln(w_ij))
```

Where `L` is the total number of synonymous codons in the sequence and `w_ij` is the relative adaptedness of codon `j` for amino acid `i`.

**References:**
- Sharp, P.M. & Li, W.H. (1987). The codon adaptation index — a measure of directional synonymous codon usage bias, and its potential applications. *Nucleic Acids Research*, 15(3), 1281–1295. https://doi.org/10.1093/nar/15.3.1281
- Ikemura, T. (1985). Codon usage and tRNA content in unicellular and multicellular organisms. *Molecular Biology and Evolution*, 2(1), 13–34. https://doi.org/10.1093/oxfordjournals.molbev.a040335

---

### Analysis 10 — Codon Pair Bias (CPB)

**Input:** One or more `.gbk` files  
**Output:** CPB score matrix (CSV), CPB heatmap (PNG)

**Codon Pair Bias (CPB)** measures the over- or underrepresentation of adjacent codon pairs in a genome. The **Codon Pair Score (CPS)** for each pair (C1, C2) is defined as:

```
CPS(C1, C2) = log[ f(C1,C2)_observed / f(C1,C2)_expected ]
```

Where the expected frequency is estimated as:

```
f(C1,C2)_expected = f(C1) × f(C2) / N_total_pairs
```

CPB analysis has been instrumental in the field of **synthetic attenuation** (creating live-attenuated vaccines by introducing underrepresented codon pairs to reduce translational efficiency), as pioneered by the Coleman et al. (2008) poliovirus study.

**References:**
- Coleman, J.R., et al. (2008). Virus attenuation by genome-scale changes in codon pair bias. *Science*, 320(5884), 1784–1787. https://doi.org/10.1126/science.1155761
- Gutman, G.A. & Hatfield, G.W. (1989). Nonrandom utilization of codon pairs in Escherichia coli. *PNAS*, 86(10), 3699–3703. https://doi.org/10.1073/pnas.86.10.3699
- Suzuki, H., et al. (2005). Codon usage patterns in Escherichia coli. *Journal of Molecular Biology*, 345(5), 1141–1147. https://doi.org/10.1016/j.jmb.2004.11.014

---

### Analysis 11 — Physicochemical Analysis (GRAVY & Aromo)

**Input:** One or more `.gbk` files  
**Output:** Violin plots (PNG), per-gene data table (CSV)

**GRAVY (Grand Average of Hydropathicity)** is calculated as the arithmetic mean of the Kyte-Doolittle hydropathy index across all amino acid residues in a protein. Positive GRAVY scores indicate hydrophobic proteins (membrane-associated); negative scores indicate hydrophilic proteins.

```
GRAVY = (Σ hydropathy_i) / protein_length
```

**Aromaticity** is the fraction of residues that are aromatic amino acids (Phe, Tyr, Trp). Both metrics are derived from the translated CDS sequences using the genetic code selected.

The Kyte-Doolittle hydropathy scale values used in KodonX are those from the original 1982 publication.

**References:**
- Kyte, J. & Doolittle, R.F. (1982). A simple method for displaying the hydropathic character of a protein. *Journal of Molecular Biology*, 157(1), 105–132. https://doi.org/10.1016/0022-2836(82)90515-0
- Lobry, J.R. & Gautier, C. (1994). Hydrophobicity, expressivity and aromaticity are the major trends of amino-acid usage in 999 Escherichia coli chromosome-encoded genes. *Nucleic Acids Research*, 22(15), 3174–3180. https://doi.org/10.1093/nar/22.15.3174

---

### Analysis 12 — Neutrality Plot (GC12 vs GC3)

**Input:** 2 or more `.gbk` files  
**Output:** Neutrality plot with linear regression (PNG), results table (CSV)

The **Neutrality Plot** (GC12 vs GC3) provides a framework to quantify the relative contribution of **mutational pressure** versus **natural selection** to codon usage bias. GC12 is the mean GC content at first and second codon positions, while GC3 is GC at the third position.

Under pure mutational drift, a slope of 1.0 is expected (GC12 ≈ GC3). Values of slope < 1 indicate that selection is constraining amino acid composition (GC12 is less free to vary than GC3). The proportion attributable to selection is estimated as `1 - slope`.

**References:**
- Sueoka, N. (1999). Translation-coupled violation of Parity Rule 2 in human genes is not the cause of heterogeneity of the DNA G+C content of third codon position. *Gene*, 238(1), 53–58. https://doi.org/10.1016/S0378-1119(99)00257-6
- Liu, Q. (2006). Comparative analysis of base compositions, codon usages, and codon pair usages in invertebrates and vertebrates. *Biochemistry (Moscow)*, 71(10), 1079–1086.
- Yang, Z. & Nielsen, R. (2008). Mutation-selection models of codon substitution and their use to estimate selective strengths on codon usage. *Molecular Biology and Evolution*, 25(3), 568–579. https://doi.org/10.1093/molbev/msm284

---

### Analysis 13 — Dinucleotide Composition

**Input:** One or more `.gbk` files  
**Output:** Dinucleotide frequency heatmap (PNG), frequency table (CSV)

Calculates the **relative dinucleotide frequency** (odds ratio) for all 16 possible dinucleotides across the complete genome sequence. The dinucleotide odds ratio (ρ) for dinucleotide XY is:

```
ρ_XY = f(XY) / (f(X) × f(Y))
```

Where `f(XY)` is the observed dinucleotide frequency and `f(X)`, `f(Y)` are the individual mononucleotide frequencies. Values > 1.0 indicate over-representation; < 1.0 indicates under-representation. The CpG dinucleotide suppression (ρ_CpG < 1) is a well-known signature of vertebrate genomes.

**References:**
- Karlin, S. & Burge, C. (1995). Dinucleotide relative abundance extremes: a genomic signature. *Trends in Genetics*, 11(7), 283–290. https://doi.org/10.1016/S0168-9525(00)89076-9
- Karlin, S., et al. (1997). Compositional differences within and between eukaryotic genomes. *PNAS*, 94(19), 10227–10232. https://doi.org/10.1073/pnas.94.19.10227

---

### Analysis 14 — PR2 Parity Plot

**Input:** One or more `.gbk` files  
**Output:** PR2 scatter plot per genome (PNG), per-gene data table (CSV)

The **PR2 (Parity Rule 2) Plot** (Sueoka, 1995) plots A3/(A3+T3) on the x-axis versus G3/(G3+C3) on the y-axis for each individual CDS, where the subscript denotes the third codon position. Under Parity Rule 2 (PR2), in the absence of selection, A ≈ T and G ≈ C within each strand, so points should cluster around (0.5, 0.5).

Deviations from this center reveal asymmetric mutational pressures: points displaced horizontally indicate A→T or T→A biases, while vertical displacements reveal G→C or C→G biases. This is particularly informative for distinguishing leading versus lagging strand mutational asymmetry.

**References:**
- Sueoka, N. (1995). Intrastrand parity rules of DNA base composition and usage biases of synonymous codons. *Journal of Molecular Evolution*, 40(3), 318–325. https://doi.org/10.1007/BF00163236
- Lobry, J.R. (1996). Asymmetric substitution patterns in the two DNA strands of bacteria. *Molecular Biology and Evolution*, 13(5), 660–665. https://doi.org/10.1093/oxfordjournals.molbev.a025626
- Mackiewicz, P., et al. (2004). Replication associated mutational pressure generating long-range correlation in DNA. *Physica A*, 273(1-2), 103–115. https://doi.org/10.1016/S0378-4371(99)00509-6

---

### Analysis 15 — tRNA Adaptation Index (tAI)

**Input:** One or more `.gbk` files (must contain `tRNA` features)  
**Output:** tAI box plots (PNG), anticodon counts (CSV), codon W-weight table (CSV)  
**Parameter:** `--superkingdom` (Bacteria or Eukaryote)

The **tRNA Adaptation Index (tAI)** measures the compatibility between a gene's codon usage and the tRNA pool of the organism, providing an estimate of **translational efficiency**. It is based on the geometric mean of codon-specific weights (`w_i`) calculated from tRNA gene copy numbers and wobble pairing efficiencies.

KodonX implements the **wobble-weighted tAI** methodology:

1. tRNA anticodons are identified from `tRNA` features in the GenBank files, and 5' anticodon base `A` is modified to `I` (inosine), following known anticodon modification rules.
2. For each codon, the wobble interaction weight `W` is computed by summing the weighted contributions of all compatible tRNA anticodons, using empirically derived wobble s-values:

```
W_codon = Σ (tRNA_count_j × s_wobble(anticodon_j, codon))
```

3. The tAI for each gene is the geometric mean of W values across all its codons:

```
tAI = exp(1/L × Σ ln(W_i))
```

Separate wobble s-value matrices are implemented for **Bacteria** and **Eukaryotes**, reflecting differences in wobble modification repertoires.

**References:**
- dos Reis, M., et al. (2004). Solving the riddle of codon usage preferences: a test for translational selection. *Nucleic Acids Research*, 32(17), 5036–5044. https://doi.org/10.1093/nar/gkh834
- Sharp, P.M., et al. (2010). Variation in the strength of selected codon usage bias among bacteria. *Nucleic Acids Research*, 38(15), 4941–4950. https://doi.org/10.1093/nar/gkq142
- Pechmann, S. & Frydman, J. (2013). Evolutionary conservation of codon optimality reveals hidden signatures of cotranslational folding. *Nature Structural & Molecular Biology*, 20(2), 237–243. https://doi.org/10.1038/nsmb.2508
- Crick, F.H.C. (1966). Codon–anticodon pairing: the wobble hypothesis. *Journal of Molecular Biology*, 19(2), 548–555. https://doi.org/10.1016/S0022-2836(66)80022-0

---

### Analysis 16 — Upstream Motifs Analysis

**Input:** One or more `.gbk` files  
**Output:** K-mer frequency bar chart (PNG), k-mer count table (CSV)  
**Parameters:** `--upstream-dist` (default: 200 bp), `--kmer-size` (default: 6)

Extracts intergenic regions immediately upstream of each CDS (up to `upstream_dist` bp) and counts all k-mers of the specified size. Top 25 most frequent k-mers are reported. This analysis can reveal conserved regulatory elements, including **Shine-Dalgarno sequences** (bacterial ribosome binding sites), **TATA boxes**, and other promoter-associated motifs.

**References:**
- Shine, J. & Dalgarno, L. (1974). The 3'-terminal sequence of Escherichia coli 16S ribosomal RNA: complementarity to nonsense triplets and ribosome binding sites. *PNAS*, 71(4), 1342–1346. https://doi.org/10.1073/pnas.71.4.1342
- Nakagawa, S., et al. (2010). The overrepresentation of regulatory sequences near the start of Shine-Dalgarno codons in Escherichia coli. *PLOS Genetics*, 6(12), e1001257. https://doi.org/10.1371/journal.pgen.1001257
- Stormo, G.D. (2000). DNA binding sites: representation and discovery. *Bioinformatics*, 16(1), 16–23. https://doi.org/10.1093/bioinformatics/16.1.16

---

### Analysis 17 — MFE Analysis (5' Structure)

**Input:** One or more `.gbk` files  
**Output:** MFE box plots (PNG), per-gene MFE table (CSV)  
**Parameter:** `--mfe-region` (default: 50 bp)  
**Requires:** ViennaRNA (`conda install -c bioconda viennarna`)

Calculates the **Minimum Free Energy (MFE)** of the predicted mRNA secondary structure of the 5' region (first N nucleotides) of each CDS, using the thermodynamic folding algorithm implemented in the ViennaRNA package (Lorenz et al., 2011).

The **translational ramp hypothesis** proposes that lowly structured (less negative MFE) 5' regions facilitate ribosome loading, while strong secondary structures can occlude the ribosome binding site and inhibit translation initiation. This analysis tests whether highly expressed genes tend to have less stable 5' structures.

Statistical comparisons across species are presented as box plots, with per-gene MFE values saved for downstream analysis.

**References:**
- Lorenz, R., et al. (2011). ViennaRNA Package 2.0. *Algorithms for Molecular Biology*, 6, 26. https://doi.org/10.1186/1748-7188-6-26
- Tuller, T., et al. (2010). An evolutionarily conserved mechanism for controlling the efficiency of protein translation. *Cell*, 141(2), 344–354. https://doi.org/10.1016/j.cell.2010.03.031
- Kudla, G., et al. (2009). Coding-sequence determinants of gene expression in Escherichia coli. *Science*, 324(5924), 255–258. https://doi.org/10.1126/science.1170160
- Gu, W., et al. (2010). Secondary structure of the 5' untranslated region of the hepatitis B virus mRNA influences translation. *Journal of Virology*, 80(22), 11460–11461. https://doi.org/10.1128/JVI.02289-05

---

### Analysis 18 — Gene Group Comparison

**Input:** One or more `.gbk` files + two gene list files (`--group1`, `--group2`)  
**Output:** Box plot comparison (PNG), summary statistics (printed)

Compares CUB metrics (ENC, GC3, CAI) between two user-defined gene groups using the **Mann-Whitney U test** (non-parametric, no assumption of normality). Additionally, a **Chi-squared test of independence** is applied to the aggregated codon count tables of both groups.

This analysis is useful for comparing, for example, **highly expressed vs. lowly expressed genes**, **horizontally transferred genes vs. core genome genes**, or **pathogenicity islands vs. housekeeping genes**.

**References:**
- Mann, H.B. & Whitney, D.R. (1947). On a test of whether one of two random variables is stochastically larger than the other. *The Annals of Mathematical Statistics*, 18(1), 50–60. https://doi.org/10.1214/aoms/1177730491
- Sharp, P.M., et al. (1988). Codon usage in Saccharomyces cerevisiae: evidence for translational selection. *Nucleic Acids Research*, 16(17), 8207–8211. https://doi.org/10.1093/nar/16.17.8207
- Karlin, S. & Mrazek, J. (2000). Predicted highly expressed genes of diverse prokaryotic genomes. *Journal of Bacteriology*, 182(18), 5238–5250. https://doi.org/10.1128/JB.182.18.5238-5250.2000

---

### Analysis 19 — Correlation with Expression

**Input:** One or more `.gbk` files + expression data file (CSV/TSV)  
**Output:** Spearman correlation scatter plots (PNG), merged data table (CSV)  
**Parameters:** `--expr-gene` (gene column), `--expr-val` (expression value column)

Integrates external gene expression data (e.g., RNA-Seq TPM/RPKM counts) with per-gene CUB metrics (CAI, ENC) and calculates **Spearman rank correlation** coefficients. Expression values are log10-transformed prior to correlation.

Under the translational selection hypothesis, highly expressed genes are expected to have **higher CAI** and **lower ENC** (more biased toward optimal codons). This analysis empirically tests this relationship for the input dataset.

**References:**
- Spearman, C. (1904). The proof and measurement of association between two things. *The American Journal of Psychology*, 15(1), 72–101. https://doi.org/10.2307/1412159
- Ikemura, T. (1982). Correlation between the abundance of yeast transfer RNAs and the occurrence of the respective codons in protein genes. *Journal of Molecular Biology*, 158(4), 573–597. https://doi.org/10.1016/0022-2836(82)90250-9
- Quax, T.E.F., et al. (2015). Codon bias as a means to fine-tune gene expression. *Molecular Cell*, 59(2), 149–161. https://doi.org/10.1016/j.molcel.2015.05.035
- Frumkin, I., et al. (2018). Gene architectures that minimize cost of gene expression. *Molecular Cell*, 60(1), 142–153. https://doi.org/10.1016/j.molcel.2015.09.008

---

## Project Structure

```
kodonx/
├── kodonx.py                        # Main entry point (GUI)
├── kodonxa.py                       # CLI entry point
├── interface.py                     # PyQt6 GUI implementation
├── constants.py                     # Genetic code tables, physicochemical constants, wobble matrices
├── core_utils.py                    # Core calculation functions (RSCU, ENC, GC3, CAI, tAI weights)
├── analysis_basic.py                # Analyses 1, 2, 8 (statistics, gene listing, composition)
├── analysis_bias.py                 # Analyses 3–7, 9, 12 (RSCU, ENC/GC3, CAI, neutrality plot)
├── analysis_advanced.py             # Analyses 10, 11, 13–16 (CPB, GRAVY, dinucleotides, PR2, tAI, motifs)
├── analysis_expression_structure.py # Analyses 17–19 (MFE, group comparison, expression correlation)
└── synthetic_biology.py             # Analysis 0 (codon optimization and harmonization)
```

---

## References

The following is a consolidated list of the primary references underlying the methodologies implemented in KodonX:

1. Angov, E., et al. (2008). Heterologous protein expression is enhanced by harmonizing the codon usage frequencies of the target gene with those of the expression host. *PLoS ONE*, 3(5), e2189. https://doi.org/10.1371/journal.pone.0002189
2. Chaney, J.L. & Clark, P.L. (2015). Roles for synonymous codon usage in protein biogenesis. *Annual Review of Biophysics*, 44, 143–166. https://doi.org/10.1146/annurev-biophys-060414-034333
3. Cock, P.J.A., et al. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics*, 25(11), 1422–1423. https://doi.org/10.1093/bioinformatics/btp163
4. Coleman, J.R., et al. (2008). Virus attenuation by genome-scale changes in codon pair bias. *Science*, 320(5884), 1784–1787. https://doi.org/10.1126/science.1155761
5. Crick, F.H.C. (1966). Codon–anticodon pairing: the wobble hypothesis. *Journal of Molecular Biology*, 19(2), 548–555. https://doi.org/10.1016/S0022-2836(66)80022-0
6. dos Reis, M., et al. (2004). Solving the riddle of codon usage preferences: a test for translational selection. *Nucleic Acids Research*, 32(17), 5036–5044. https://doi.org/10.1093/nar/gkh834
7. Frumkin, I., et al. (2018). Gene architectures that minimize cost of gene expression. *Molecular Cell*, 60(1), 142–153. https://doi.org/10.1016/j.molcel.2015.09.008
8. Gustafsson, C., et al. (2004). Codon bias and heterologous protein expression. *Trends in Biotechnology*, 22(7), 346–353. https://doi.org/10.1016/j.tibtech.2004.04.006
9. Gutman, G.A. & Hatfield, G.W. (1989). Nonrandom utilization of codon pairs in Escherichia coli. *PNAS*, 86(10), 3699–3703. https://doi.org/10.1073/pnas.86.10.3699
10. Ikemura, T. (1981). Correlation between the abundance of Escherichia coli transfer RNAs and the occurrence of the respective codons in its protein genes. *Journal of Molecular Biology*, 146(1), 1–21. https://doi.org/10.1016/0022-2836(81)90363-6
11. Ikemura, T. (1982). Correlation between the abundance of yeast transfer RNAs and the occurrence of the respective codons in protein genes. *Journal of Molecular Biology*, 158(4), 573–597. https://doi.org/10.1016/0022-2836(82)90250-9
12. Ikemura, T. (1985). Codon usage and tRNA content in unicellular and multicellular organisms. *Molecular Biology and Evolution*, 2(1), 13–34. https://doi.org/10.1093/oxfordjournals.molbev.a040335
13. Karlin, S. & Burge, C. (1995). Dinucleotide relative abundance extremes: a genomic signature. *Trends in Genetics*, 11(7), 283–290. https://doi.org/10.1016/S0168-9525(00)89076-9
14. Karlin, S. & Mrazek, J. (1996). What drives codon choices in human genes? *Journal of Molecular Biology*, 262(4), 459–472. https://doi.org/10.1006/jmbi.1996.0528
15. Karlin, S. & Mrazek, J. (2000). Predicted highly expressed genes of diverse prokaryotic genomes. *Journal of Bacteriology*, 182(18), 5238–5250. https://doi.org/10.1128/JB.182.18.5238-5250.2000
16. Karlin, S., et al. (1997). Compositional differences within and between eukaryotic genomes. *PNAS*, 94(19), 10227–10232. https://doi.org/10.1073/pnas.94.19.10227
17. Kudla, G., et al. (2009). Coding-sequence determinants of gene expression in Escherichia coli. *Science*, 324(5924), 255–258. https://doi.org/10.1126/science.1170160
18. Kyte, J. & Doolittle, R.F. (1982). A simple method for displaying the hydropathic character of a protein. *Journal of Molecular Biology*, 157(1), 105–132. https://doi.org/10.1016/0022-2836(82)90515-0
19. Liu, Q. (2006). Comparative analysis of base compositions, codon usages, and codon pair usages in invertebrates and vertebrates. *Biochemistry (Moscow)*, 71(10), 1079–1086.
20. Lobry, J.R. (1996). Asymmetric substitution patterns in the two DNA strands of bacteria. *Molecular Biology and Evolution*, 13(5), 660–665. https://doi.org/10.1093/oxfordjournals.molbev.a025626
21. Lobry, J.R. & Gautier, C. (1994). Hydrophobicity, expressivity and aromaticity are the major trends of amino-acid usage in 999 Escherichia coli chromosome-encoded genes. *Nucleic Acids Research*, 22(15), 3174–3180. https://doi.org/10.1093/nar/22.15.3174
22. Lorenz, R., et al. (2011). ViennaRNA Package 2.0. *Algorithms for Molecular Biology*, 6, 26. https://doi.org/10.1186/1748-7188-6-26
23. Mann, H.B. & Whitney, D.R. (1947). On a test of whether one of two random variables is stochastically larger than the other. *The Annals of Mathematical Statistics*, 18(1), 50–60. https://doi.org/10.1214/aoms/1177730481
24. Nakagawa, S., et al. (2010). The overrepresentation of regulatory sequences near the start of Shine-Dalgarno codons in Escherichia coli. *PLOS Genetics*, 6(12), e1001257. https://doi.org/10.1371/journal.pgen.1001257
25. Novembre, J. (2002). Accounting for background nucleotide composition when measuring codon usage bias. *Molecular Biology and Evolution*, 19(8), 1390–1394. https://doi.org/10.1093/oxfordjournals.molbev.a004198
26. Pedregosa, F., et al. (2011). Scikit-learn: Machine learning in Python. *Journal of Machine Learning Research*, 12, 2825–2830.
27. Pechmann, S. & Frydman, J. (2013). Evolutionary conservation of codon optimality reveals hidden signatures of cotranslational folding. *Nature Structural & Molecular Biology*, 20(2), 237–243. https://doi.org/10.1038/nsmb.2508
28. Puigbò, P., et al. (2007). OPTIMIZER: A web server for optimizing the codon usage of DNA sequences. *Nucleic Acids Research*, 35(Web Server issue), W126–W131. https://doi.org/10.1093/nar/gkm219
29. Quax, T.E.F., et al. (2015). Codon bias as a means to fine-tune gene expression. *Molecular Cell*, 59(2), 149–161. https://doi.org/10.1016/j.molcel.2015.05.035
30. Sharp, P.M. & Li, W.H. (1987). The codon adaptation index — a measure of directional synonymous codon usage bias, and its potential applications. *Nucleic Acids Research*, 15(3), 1281–1295. https://doi.org/10.1093/nar/15.3.1281
31. Sharp, P.M., et al. (1988). Codon usage in Saccharomyces cerevisiae: evidence for translational selection. *Nucleic Acids Research*, 16(17), 8207–8211. https://doi.org/10.1093/nar/16.17.8207
32. Sharp, P.M., et al. (2010). Variation in the strength of selected codon usage bias among bacteria. *Nucleic Acids Research*, 38(15), 4941–4950. https://doi.org/10.1093/nar/gkq142
33. Shine, J. & Dalgarno, L. (1974). The 3'-terminal sequence of Escherichia coli 16S ribosomal RNA: complementarity to nonsense triplets and ribosome binding sites. *PNAS*, 71(4), 1342–1346. https://doi.org/10.1073/pnas.71.4.1342
34. Spearman, C. (1904). The proof and measurement of association between two things. *The American Journal of Psychology*, 15(1), 72–101. https://doi.org/10.2307/1412159
35. Stormo, G.D. (2000). DNA binding sites: representation and discovery. *Bioinformatics*, 16(1), 16–23. https://doi.org/10.1093/bioinformatics/16.1.16
36. Sueoka, N. (1962). On the genetic basis of variation and heterogeneity of DNA base composition. *PNAS*, 48(4), 582–592. https://doi.org/10.1073/pnas.48.4.582
37. Sueoka, N. (1988). Directional mutation pressure and neutral molecular evolution. *PNAS*, 85(8), 2653–2657. https://doi.org/10.1073/pnas.85.8.2653
38. Sueoka, N. (1995). Intrastrand parity rules of DNA base composition and usage biases of synonymous codons. *Journal of Molecular Evolution*, 40(3), 318–325. https://doi.org/10.1007/BF00163236
39. Sueoka, N. (1999). Translation-coupled violation of Parity Rule 2 in human genes is not the cause of heterogeneity of the DNA G+C content of third codon position. *Gene*, 238(1), 53–58. https://doi.org/10.1016/S0378-1119(99)00257-6
40. Tuller, T., et al. (2010). An evolutionarily conserved mechanism for controlling the efficiency of protein translation. *Cell*, 141(2), 344–354. https://doi.org/10.1016/j.cell.2010.03.031
41. Wright, F. (1990). The 'effective number of codons' used in a gene. *Gene*, 87(1), 23–29. https://doi.org/10.1016/0378-1119(90)90491-9
42. Yang, Z. & Nielsen, R. (2008). Mutation-selection models of codon substitution and their use to estimate selective strengths on codon usage. *Molecular Biology and Evolution*, 25(3), 568–579. https://doi.org/10.1093/molbev/msm284

---

*KodonX — Developed for codon usage analysis in comparative genomics and synthetic biology research.*
