# Kódon-X — Integrated Platform for Codon Usage Bias Analysis

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Kódon-X** is an integrated bioinformatics platform for analyzing Codon Usage Bias (CUB) in prokaryotic and eukaryotic genomes. Developed in Python 3 with a modular architecture, the system consolidates nineteen analytical modules into a unified environment accessible via both a graphical user interface (GUI) and a command-line interface (CLI). It combines classical CUB metrics (RSCU, ENC, CAI) with advanced biophysical descriptors (tAI, GRAVY, Aromaticity, MFE), comparative genomics tools, and a synthetic biology module for rational gene design.

> **Citation:** Barbanti ACC, Araújo VSC, Rosário AEC, Carneiro DVD, Tavares GC, Aburjaile FF. *Kódon-X: An Integrated Platform for Analyzing Codon Usage Bias in Prokaryotic and Eukaryotic Genomes.* Laboratory of Integrative Bioinformatics (LBI), Federal University of Minas Gerais (UFMG), Brazil.

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
  - [Analysis 4 — Comparative RSCU (Clustermap + PCA + Histograms)](#analysis-4--comparative-rscu-clustermap--pca--histograms)
  - [Analysis 5 — RSCU Correlation (2 Genomes)](#analysis-5--rscu-correlation-2-genomes)
  - [Analysis 6 — ENC vs GC3 (Wright Plot)](#analysis-6--enc-vs-gc3-wright-plot)
  - [Analysis 7 — Genomic Composition](#analysis-7--genomic-composition)
  - [Analysis 8 — Optimal, Rare Codons & CAI](#analysis-8--optimal-rare-codons--cai)
  - [Analysis 9 — Codon Pair Bias (CPB)](#analysis-9--codon-pair-bias-cpb)
  - [Analysis 10 — Physicochemical Analysis (GRAVY & Aromo)](#analysis-10--physicochemical-analysis-gravy--aromo)
  - [Analysis 11 — Neutrality Plot (GC12 vs GC3)](#analysis-11--neutrality-plot-gc12-vs-gc3)
  - [Analysis 12 — Dinucleotide Composition](#analysis-12--dinucleotide-composition)
  - [Analysis 13 — PR2 Parity Plot](#analysis-13--pr2-parity-plot)
  - [Analysis 14 — tRNA Adaptation Index (tAI)](#analysis-14--trna-adaptation-index-tai)
  - [Analysis 15 — Upstream Motifs Analysis](#analysis-15--upstream-motifs-analysis)
  - [Analysis 16 — MFE Analysis (5' Structure)](#analysis-16--mfe-analysis-5-structure)
  - [Analysis 17 — Gene Group Comparison](#analysis-17--gene-group-comparison)
  - [Analysis 18 — Correlation with Expression](#analysis-18--correlation-with-expression)
- [Output Architecture](#output-architecture)
- [Case Studies](#case-studies)
- [Project Structure](#project-structure)
- [References](#references)

---

## Features Overview

- **18 analysis modules** covering codon bias, comparative genomics, translation efficiency, physicochemistry, and gene expression.
- **Synthetic Biology module** for sequence optimization and harmonization targeting heterologous expression systems.
- Supports 4 **NCBI genetic code tables** (Standard, Vertebrate Mitochondrial, Mold/Protozoan, Bacterial/Plant Plastid).
- Optional **gene-level filtering** via locus tag lists for targeted subset analysis.
- Dual interface: **PyQt6 GUI** with real-time progress monitoring and **CLI** for automated pipelines.
- All outputs are publication-ready **PNG figures** (150 DPI, Agg backend) and **CSV tables** (semicolon-delimited, decimal point notation) with genome-identifier-based dynamic filenames.
- Threading-based execution ensures a responsive interface during long analyses.

---

## Installation

### Prerequisites

```bash
pip install PyQt6 Pillow pandas matplotlib seaborn numpy scipy scikit-learn biopython
```

### Optional (for MFE Analysis — Analysis 16)

ViennaRNA must be installed via Conda:

```bash
conda install -c bioconda viennarna
```

### Running

**GUI:**
```bash
python KodonX.py
```

**CLI:**
```bash
python KodonXa.py -i /path/to/gbk_folder -o /path/to/output -a <analysis_id>
```

---

## Input Format

Kódon-X accepts annotated genome files in **GenBank flat-file format** (`.gbk`, `.gb`, `.gbff`), as provided by NCBI. Files may contain multiple contigs (typical of short-read assemblies), which are processed in aggregate. An automatic quality filter excludes sequences containing ambiguous nucleotides (N) and those shorter than three nucleotides.

Files must contain `CDS` features with at least a `locus_tag` or `gene` qualifier; `tRNA` features are additionally required for Analysis 14 (tAI).

Multiple files in a single input directory enable comparative analyses across multiple genomes or strains.

### Gene Filter File

A plain-text file with one `locus_tag` or gene name per line passed via `-f` restricts all analyses to a specific gene subset:

```
b0001
b0002
recA
rpoB
```

### Synthetic Biology Input

In the synthetic biology module, an exogenous coding sequence in plain text is entered alongside one or more host genome files already loaded in the workspace.

---

## Graphical Interface (GUI)

Launch with `python KodonX.py`. The GUI provides:

- Multi-genome file selection dialog.
- Dropdown menus for genetic code table and analysis selection.
- Gene filter input via text field or file.
- Real-time progress tracking and log console (Python stdout redirected to the interface).
- Integrated image viewer with zoom, pan, and individual export.
- Reports CDS counts, codon tallies, per-genome metric values, and annotation warnings in real time.

---

## Command-Line Interface (CLI)

```
python KodonXa.py -i <input_dir> -o <output_dir> -a <analysis_id> [options]

Options:
  -g, --genetic-code INT    Genetic code table ID (default: 1)
  -f, --filter-file FILE    Text file with locus_tags for filtering
  --start-codon STR         Start codon filter for Analysis 1 (default: ATG)
  --superkingdom STR        Bacteria or Eukaryote, for tAI wobble rules (Analysis 14)
  --upstream-dist INT       Upstream distance in bp (Analysis 15, default: 200)
  --kmer-size INT           K-mer size for motif search (Analysis 15, default: 6)
  --mfe-region INT          5' region size for MFE (Analysis 16, default: 50)
  --group1 FILE             Gene list for Group 1 (Analysis 17)
  --group2 FILE             Gene list for Group 2 (Analysis 17)
  --expr-file FILE          Expression data CSV/TSV (Analysis 18)
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
| 1  | Standard / Universal (NCBI Table 1) |
| 2  | Vertebrate Mitochondrial (NCBI Table 2) |
| 4  | Mold, Protozoan, Coelenterate Mitochondrial (NCBI Table 4) |
| 11 | Bacterial, Archaeal & Plant Plastid (NCBI Table 11) |

Genetic code definitions follow NCBI: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi

For each table, reverse amino acid maps and synonymous codon lists are generated automatically and used in all RSCU, ENC, and CAI calculations. Support for organisms with atypical translation rules requires only the insertion of a new matrix into the `GENETIC_CODE_TABLES` dictionary in `constants.py`.

---

## Analysis Modules

---

### Analysis 0 — Synthetic Biology

**Mode:** `optimize` or `harmonize`
**Input:** Host `.gbk` file(s) + raw DNA sequence string
**Output:** Recoded DNA sequence (printed to terminal / GUI) and exported GenBank file

This module recodes an input CDS for expression in a target host organism using two strategies. The output is formatted as a standard GenBank file with a CDS feature annotation, ensuring direct compatibility with in silico cloning software and commercial gene synthesis platforms.

#### Codon Optimization (Maximization)

`optimize` mode replaces each codon with the **most frequent synonymous codon** in the host genome (highest RSCU value). This maximizes translational elongation speed by ensuring all codons are drawn from the pool most abundant in the host tRNA repertoire.

**Methodology:** For each amino acid, the optimal codon is `argmax(RSCU_i)` across all synonymous codons, following Sharp & Li (1986).

**References:**
- Sharp PM, Li WH (1987). The codon adaptation index — a measure of directional synonymous codon usage bias, and its potential applications. *Nucleic Acids Research*, 15(3), 1281–1295. https://doi.org/10.1093/nar/15.3.1281
- Gustafsson C, et al. (2004). Codon bias and heterologous protein expression. *Trends in Biotechnology*, 22(7), 346–353. https://doi.org/10.1016/j.tibtech.2004.04.006
- Wu G, Bashir-Bello N, Freeland SJ (2006). The Synthetic Gene Designer: A flexible web platform to explore sequence manipulation for heterologous expression. *Protein Expression and Purification*, 47, 441–445. https://doi.org/10.1016/j.pep.2005.10.020

#### Codon Harmonization (Rank-Based)

`harmonize` mode preserves the **relative codon usage pattern** of the input sequence rather than maximizing speed. The algorithm ranks synonymous codons per amino acid by frequency, then maps each codon in the input sequence to the codon occupying the same rank position in the host genome. This mimics the translational kinetics of the source organism, preserving ribosomal pausing at positions associated with cotranslational protein folding, following Angov et al. (2008).

**References:**
- Angov E, et al. (2008). Heterologous protein expression is enhanced by harmonizing the codon usage frequencies of the target gene with those of the expression host. *PLoS ONE*, 3(5), e2189. https://doi.org/10.1371/journal.pone.0002189
- Chaney JL & Clark PL (2015). Roles for synonymous codon usage in protein biogenesis. *Annual Review of Biophysics*, 44, 143–166. https://doi.org/10.1146/annurev-biophys-060414-034333

---

### Analysis 1 — Statistics and CDS

**Input:** One or more `.gbk` files
**Output:** `genome_statistics.csv`, `cds_analysis.csv`

Quality-control step for genomic annotation. Processes multiple files in batch and reports, per genome: number of contigs, total length (bp), and overall GC content (%). For each annotated CDS: identifier (locus_tag, protein_id, or gene), observed start codon, compliance with the configured start codon, and ORF size (nt).

**References:**
- Cock PJA, et al. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics*, 25(11), 1422–1423. https://doi.org/10.1093/bioinformatics/btp163

---

### Analysis 2 — Gene Listing

**Input:** Exactly 1 `.gbk` file
**Output:** `gene_list_<filename>.csv`

Extracts and organizes all entries labeled as CDS, tRNA, rRNA, and gene from a GenBank file, reporting locus_tag, gene name, functional product, and feature type. The CSV output facilitates the construction of gene filter lists for subsequent targeted analyses.

---

### Analysis 3 — Individual RSCU Heatmap

**Input:** Exactly 1 `.gbk` file
**Output:** RSCU heatmap (PNG), polar plot (PNG), line plot (PNG), codon counts table (CSV)
**Metrics reported:** RSCU, ENC, GC3, CAI

**Relative Synonymous Codon Usage (RSCU)** is the ratio of the observed frequency of a codon to the frequency expected if all synonymous codons were used equally (Sharp & Li, 1986):

```
RSCU(c_i) = X_ij / (1/n_A × Σ X_ij)
```

Where `X_ij` is the count of codon `j` for amino acid `A`, and `n_A` is the number of synonymous codons. RSCU = 1.0 indicates unbiased usage; values > 1.0 indicate preferred codons; values < 1.0 indicate underused codons.

Three figures are generated: (i) a 4×16 heatmap on the canonical genetic code grid organized by positions 1 and 2; (ii) a polar graph representing codon usage in circular coordinates; (iii) a line plot mapping codons by amino acid. Global metrics (ENC, GC3, CAI) are reported in the figure header.

**References:**
- Sharp PM, Li WH (1986). Codon usage in regulatory genes in Escherichia coli does not reflect selection for 'rare' codons. *Nucleic Acids Research*, 14, 7737–7749. https://doi.org/10.1093/nar/14.19.7737

---

### Analysis 4 — Comparative RSCU (Clustermap + PCA + Histograms)

**Input:** 2 or more `.gbk` files
**Output:** Comparative RSCU matrix (CSV), 8 figures (PNG)

Unified comparative RSCU module. Compares codon usage profiles across multiple genomes via an integrated approach combining descriptive statistics, eight visualization methods, and multivariate analysis:

1. **Clustermap** — UPGMA hierarchical clustering (Euclidean distance; Sneath & Sokal 1973)
2. **PCA projection** — Z-score normalization (StandardScaler) prior to PCA with explained variance reported on axes
3. **PCA biplot (loadings)** — identifies the codons driving variance
4. **Variance boxplot** — RSCU variance per codon across genomes
5. **Aggregated boxplot** — RSCU distribution per amino acid for all samples, with reference line at RSCU = 1.0
6. **Stacked bar charts** — fractional synonymous codon contribution per amino acid per genome
7. **Comparative line plots** — RSCU values per codon organized by amino acid for each genome
8. **Global heatmap** — RSCU matrix grouped by amino acid; color intensity represents preference level

The consolidated RSCU matrix (genomes × 64 codons) is exported in CSV format.

**References:**
- Karlin S & Mrazek J (1996). What drives codon choices in human genes? *Journal of Molecular Biology*, 262(4), 459–472. https://doi.org/10.1006/jmbi.1996.0528
- Pedregosa F, et al. (2011). Scikit-learn: Machine learning in Python. *Journal of Machine Learning Research*, 12, 2825–2830.
- Sneath PH & Sokal RR (1973). *Numerical Taxonomy*. W. H. Freeman, San Francisco.
- Novembre J (2002). Accounting for background nucleotide composition when measuring codon usage bias. *Molecular Biology and Evolution*, 19(8), 1390–1394. https://doi.org/10.1093/oxfordjournals.molbev.a004201

---

### Analysis 5 — RSCU Correlation (2 Genomes)

**Input:** Exactly 2 `.gbk` files
**Output:** 3 figures (PNG): scatter plot, delta RSCU bar plot, divergent codons heatmap

Assesses evolutionary and mutational proximity between a pair of genomes via three plots: (i) a scatter plot of RSCU vectors with the **Pearson correlation coefficient** and associated p-value (Virtanen et al., 2020); (ii) a bar plot of absolute differences (ΔRSCU) per codon; (iii) a heatmap focused on the 20 most divergent codons. A high R indicates similar translational selection pressure and tRNA repertoire.

**References:**
- Sharp PM, et al. (1995). Codon usage in Escherichia coli: genes for proteins with differing cellular abundance. *EMBO Journal*, 5(6), 1233–1244.
- Ikemura T (1981). Correlation between the abundance of Escherichia coli transfer RNAs and the occurrence of the respective codons in its protein genes. *Journal of Molecular Biology*, 146(1), 1–21. https://doi.org/10.1016/0022-2836(81)90363-6
- Virtanen P, et al. (2020). SciPy 1.0: fundamental algorithms for scientific computing in Python. *Nature Methods*, 17, 261–272. https://doi.org/10.1038/s41592-019-0686-2

---

### Analysis 6 — ENC vs GC3 (Wright Plot)

**Input:** 2 or more `.gbk` files
**Output:** Wright plot (PNG), ENC/GC3 per-gene results table (CSV), KDE plots and comparative bar charts (PNG)

The **Wright Plot** (ENC vs GC3) distinguishes the driving forces behind codon bias. **Each point represents an individual CDS gene**, colored by species. Points on or near the theoretical expected curve indicate codon bias driven primarily by **mutational pressure** (GC-biased mutation). Points falling substantially below the curve indicate **natural selection** on translational efficiency as an additional force.

**ENC (Effective Number of Codons)** quantifies codon bias independently of gene length, ranging from 20 (extreme: one codon per amino acid) to 61 (no bias). Following Wright (1990), it is calculated from the mean homozygosity F_k across amino acid families:

```
F_k = (n × Σp²_i − 1) / (n − 1)
```

```
ENC = 2 + 9/F₂ + 1/F₃ + 5/F₄ + 3/F₆
```

**GC3** is the GC content at the third codon position (synonymous sites), expressed as a proportion (0–1).

The theoretical expected ENC curve under pure mutational bias follows the formula derived by Wright (1990):

```
Nc_exp = 2 + s + 29 / (s² + (1−s)²)
```

Where `s = GC3` (as a proportion). This curve reaches a maximum of ≈ 60.5 at GC3 = 0.5 and decreases symmetrically toward the extremes. KDE density plots for ENC and GC3 distributions and comparative bar charts of population means are also generated.

**References:**
- Wright F (1990). The 'effective number of codons' used in a gene. *Gene*, 87(1), 23–29. https://doi.org/10.1016/0378-1119(90)90491-9
- Sueoka N (1988). Directional mutation pressure and neutral molecular evolution. *PNAS*, 85(8), 2653–2657. https://doi.org/10.1073/pnas.85.8.2653
- Novembre J (2002). Accounting for background nucleotide composition when measuring codon usage bias. *Molecular Biology and Evolution*, 19(8), 1390–1394. https://doi.org/10.1093/oxfordjournals.molbev.a004201

---

### Analysis 7 — Genomic Composition

**Input:** One or more `.gbk` files
**Output:** Nucleotide composition grouped bar chart (PNG), composition results table (CSV)

Calculates the percentage of each nucleotide (A, T, G, C) and total GC content across complete genomic sequences. Generates a grouped horizontal bar chart comparing nucleotide composition across multiple species, enabling identification of GC-biased versus AT-biased genomes.

**References:**
- Sueoka N (1962). On the genetic basis of variation and heterogeneity of DNA base composition. *PNAS*, 48(4), 582–592. https://doi.org/10.1073/pnas.48.4.582

---

### Analysis 8 — Optimal, Rare Codons & CAI

**Input:** One or more `.gbk` files
**Output:** Optimal/rare codon tables (CSV), 4 figures (PNG)

**Optimal codons** are those with the highest RSCU in their synonymous group (RSCU > 1.2; threshold for preferred usage). **Rare codons** are the lowest RSCU in their group (RSCU < 0.8). Both are identified per amino acid per genome.

The **Codon Adaptation Index (CAI)** measures adaptation of a gene's codon usage to the host codon preference, following Sharp & Li (1987). For each codon, the relative adaptedness weight w_i is:

```
w_i = RSCU(c_i) / max(RSCU_A)
```

Where `max(RSCU_A)` is the maximum RSCU among synonymous codons for amino acid A. Stop codons and codons encoding single-codon amino acids (Met, Trp) are excluded from the calculation.

```
CAI = exp(1/L × Σ ln(w_i))
```

Where `L` is the total number of eligible codons in the sequence. Four figures are generated: (i) CAI comparison bar plot, (ii) optimal/rare codon frequency bar plot, (iii) CAI vs. optimal codon count scatter plot, (iv) RSCU distribution box plot by genome for optimal codons.

**References:**
- Sharp PM & Li WH (1987). The codon adaptation index — a measure of directional synonymous codon usage bias, and its potential applications. *Nucleic Acids Research*, 15(3), 1281–1295. https://doi.org/10.1093/nar/15.3.1281
- Ikemura T (1985). Codon usage and tRNA content in unicellular and multicellular organisms. *Molecular Biology and Evolution*, 2(1), 13–34. https://doi.org/10.1093/oxfordjournals.molbev.a040335

---

### Analysis 9 — Codon Pair Bias (CPB)

**Input:** One or more `.gbk` files
**Output:** CPS matrix (CSV), 64×64 CPB heatmap (PNG), CPS distribution histogram (PNG), top-20 over/underrepresented pairs bar charts (PNG)

**Codon Pair Bias (CPB)** quantifies the statistical preference or aversion for adjacent codon pairs in coding sequences. The **Codon Pair Score (CPS)** for each pair (C1, C2), following Coleman et al. (2008):

```
CPS(C1, C2) = log[ F_obs(C1,C2) / F_exp(C1,C2) ]
```

Where the expected frequency under independence is:

```
F_exp(C1,C2) = F(C1) × F(C2) × [F(AA1,AA2) / (F(AA1) × F(AA2))]
```

This formulation corrects for amino acid pair frequency, isolating the codon-level bias from the amino acid-level bias. CPB has been applied to synthetic attenuation of viruses (poliovirus, influenza) by introducing underrepresented codon pairs to reduce translational efficiency without altering the encoded protein.

**References:**
- Coleman JR, et al. (2008). Virus attenuation by genome-scale changes in codon pair bias. *Science*, 320(5884), 1784–1787. https://doi.org/10.1126/science.1155761
- Gutman GA & Hatfield GW (1989). Nonrandom utilization of codon pairs in Escherichia coli. *PNAS*, 86(10), 3699–3703. https://doi.org/10.1073/pnas.86.10.3699

---

### Analysis 10 — Physicochemical Analysis (GRAVY & Aromo)

**Input:** One or more `.gbk` files
**Output:** Scatter plot (PNG), hexbin density map (PNG), box plots (PNG), per-gene data table (CSV)

Each CDS is translated using the selected genetic code table, then two physicochemical descriptors are calculated per gene.

**GRAVY (Grand Average of Hydropathicity)** — arithmetic mean of the Kyte-Doolittle hydropathy index across all residues (Kyte & Doolittle, 1982):

```
GRAVY = (Σ H_i) / L
```

Positive values indicate hydrophobic proteins (membrane-associated); negative values indicate hydrophilic proteins. Isoleucine and valine have positive values; arginine and aspartic acid have negative values.

**Aromaticity** — fraction of aromatic residues (Phe [F], Tyr [Y], Trp [W]) following Lobry & Gautier (1994):

```
Aromo = (count_F + count_Y + count_W) / L
```

Results are displayed in a two-dimensional GRAVY vs. Aromaticity scatter, supplemented by hexbin density maps and box plots comparing proteomic means per species.

**References:**
- Kyte J & Doolittle RF (1982). A simple method for displaying the hydropathic character of a protein. *Journal of Molecular Biology*, 157(1), 105–132. https://doi.org/10.1016/0022-2836(82)90515-0
- Lobry JR & Gautier C (1994). Hydrophobicity, expressivity and aromaticity are the major trends of amino-acid usage in 999 Escherichia coli chromosome-encoded genes. *Nucleic Acids Research*, 22(15), 3174–3180. https://doi.org/10.1093/nar/22.15.3174

---

### Analysis 11 — Neutrality Plot (GC12 vs GC3)

**Input:** 2 or more `.gbk` files
**Output:** Neutrality plot with per-species linear regressions (PNG), GC by codon position bar charts (PNG), KDE density plots (PNG), regression statistics table (CSV)

The **Neutrality Plot** (GC12 vs GC3), following Sueoka (1988, 1999), quantifies the relative contribution of **mutational pressure** versus **natural selection** to codon usage. For each species, an OLS linear regression is fitted using **one point per CDS gene** (GC12 of the gene vs. GC3 of the gene). The slope, intercept, R², and p-value are reported per species.

- **Slope ≈ 1.0:** mutational pressure dominates (all codon positions shift equally with GC content)
- **Slope ≈ 0.0:** translational selection dominates (GC3 diverges from GC12 independently)
- **0 < slope < 1:** combined action of mutation and selection; proportion attributable to selection ≈ 1 − slope

GC12 is the mean GC content at first and second codon positions (functionally constrained, encoding amino acid identity). GC3 is GC at the third (synonymous) position. Under pure mutational drift, GC12 ≈ GC3 (slope ≈ 1).

**References:**
- Sueoka N (1988). Directional mutation pressure and neutral molecular evolution. *PNAS*, 85(8), 2653–2657. https://doi.org/10.1073/pnas.85.8.2653
- Sueoka N (1999). Two aspects of DNA base composition: G+C content and translation-coupled deviation from intra-strand rule of A=T and G=C. *Journal of Molecular Evolution*, 49, 49–62. https://doi.org/10.1007/PL00006534
- Yang Z & Nielsen R (2008). Mutation-selection models of codon substitution and their use to estimate selective strengths on codon usage. *Molecular Biology and Evolution*, 25(3), 568–579. https://doi.org/10.1093/molbev/msm284

---

### Analysis 12 — Dinucleotide Composition

**Input:** One or more `.gbk` files
**Output:** ρ_XY heatmap (PNG), variance boxplot (PNG), mean barplot (PNG), odds ratio table (CSV)

Constructs the **genomic dinucleotide signature** following Karlin & Burge (1995). For each of the 16 possible dinucleotides XY, the **odds ratio (ρ_XY)** is calculated:

```
ρ_XY = f(XY) / (f(X) × f(Y))
```

Where `f(XY)` is the observed dinucleotide frequency and `f(X)`, `f(Y)` are the individual mononucleotide frequencies across the complete genome sequence. Values ρ > 1.0 indicate over-representation; ρ < 1.0 indicates under-representation. The well-known CpG suppression (ρ_CG < 1) is a genomic signature of vertebrate genomes. Results are visualized as comparative heatmaps between species, variance boxplots, and mean barplots.

**References:**
- Karlin S & Burge C (1995). Dinucleotide relative abundance extremes: a genomic signature. *Trends in Genetics*, 11(7), 283–290. https://doi.org/10.1016/S0168-9525(00)89076-9
- Karlin S, et al. (1997). Compositional differences within and between eukaryotic genomes. *PNAS*, 94(19), 10227–10232. https://doi.org/10.1073/pnas.94.19.10227

---

### Analysis 13 — PR2 Parity Plot

**Input:** One or more `.gbk` files
**Output:** PR2 scatter plot per genome (PNG), KDE density map (PNG), skew boxplots (PNG), histograms (PNG), per-gene data table (CSV)

The **PR2 (Parity Rule 2) Plot** (Sueoka, 1995) plots `A3/(A3+T3)` on the x-axis versus `G3/(G3+C3)` on the y-axis for each individual CDS, where subscript 3 denotes the third codon position. Exclusively four-fold degenerate codon families are used; stop codons are excluded.

Under Parity Rule 2 (PR2), in the absence of mutation bias, A ≈ T and G ≈ C within each strand, so points should cluster around (0.5, 0.5). Deviations reveal asymmetric mutational pressures: horizontal displacement indicates A→T or T→A bias; vertical displacement indicates G→C or C→G bias. AT-skew and GC-skew are quantified and reported. This analysis is particularly informative for detecting leading-vs-lagging strand mutational asymmetry.

**References:**
- Sueoka N (1995). Intrastrand parity rules of DNA base composition and usage biases of synonymous codons. *Journal of Molecular Evolution*, 40(3), 318–325. https://doi.org/10.1007/BF00163236
- Lobry JR (1996). Asymmetric substitution patterns in the two DNA strands of bacteria. *Molecular Biology and Evolution*, 13(5), 660–665. https://doi.org/10.1093/oxfordjournals.molbev.a025626
- Parvathy ST, Udayasuriyan V & Bhadana V (2022). Codon usage bias. *Molecular Biology Reports*, 49, 539–565. https://doi.org/10.1007/s11033-021-06749-4

---

### Analysis 14 — tRNA Adaptation Index (tAI)

**Input:** One or more `.gbk` files (must contain `tRNA` features)
**Output:** tAI distribution box plots (PNG), tAI vs. gene length scatter (PNG), W-values bar plot (PNG), anticodon counts (CSV), codon W-weight table (CSV)
**Parameter:** `--superkingdom` (Bacteria or Eukaryote)

The **tRNA Adaptation Index (tAI)** measures the compatibility between a gene's codon usage and the organism's tRNA pool, estimating translational efficiency. Implementation follows dos Reis, Savva & Wernisch (2004) in three steps:

**Step 1 — Anticodon extraction:** tRNA anticodons are extracted from `tRNA` features in the GenBank file. The post-transcriptional modification of adenine (A) to inosine (I) at anticodon position 34 is applied automatically via `ANTICODON_MODIFICATION_MAP`.

**Step 2 — Codon weight calculation:** For each codon, the wobble interaction weight W is the sum of contributions from all compatible tRNA anticodons:

```
W_codon = Σ (tRNA_count_j × s_wobble(anticodon_j[34], codon[3]))
```

Where `s_wobble` is the selectivity parameter for the wobble base pair between the first anticodon base (position 34) and the third codon base (wobble position). Separate 9-parameter wobble vectors are implemented for **Bacteria** and **Eukarya**, following the values published in dos Reis et al. (2004).

**Step 3 — tAI per gene:** W values are normalized by the maximum observed W. The tAI of each gene is the geometric mean of normalized W values across all its codons:

```
tAI = exp(1/L × Σ ln(W_i))
```

**References:**
- dos Reis M, Savva R & Wernisch L (2004). Solving the riddle of codon usage preferences: a test for translational selection. *Nucleic Acids Research*, 32(17), 5036–5044. https://doi.org/10.1093/nar/gkh834
- Crick FHC (1966). Codon–anticodon pairing: the wobble hypothesis. *Journal of Molecular Biology*, 19(2), 548–555. https://doi.org/10.1016/S0022-2836(66)80022-0
- Pechmann S & Frydman J (2013). Evolutionary conservation of codon optimality reveals hidden signatures of cotranslational folding. *Nature Structural & Molecular Biology*, 20(2), 237–243. https://doi.org/10.1038/nsmb.2508

---

### Analysis 15 — Upstream Motifs Analysis

**Input:** One or more `.gbk` files
**Output:** K-mer frequency bar chart (PNG), Zipf rank-frequency curve (PNG), GC vs. frequency scatter (PNG), k-mer count table (CSV)
**Parameters:** `--upstream-dist` (default: 200 bp), `--kmer-size` (default: 6)

Extracts intergenic regions immediately upstream of each CDS (up to `upstream_dist` bp) and counts all k-mers of the specified size using sliding windows. Handles both DNA strands; reverse complementation is applied for negative-strand genes. Top-25 most frequent motifs are reported and ranked. Generates bar charts of the 25 most prevalent motifs, Zipf's Law ranking curves (rank-frequency on logarithmic scale), and GC content vs. absolute frequency scatter plots. This analysis can reveal **Shine-Dalgarno sequences** (bacterial ribosome binding sites), TATA boxes, and other conserved regulatory elements.

**References:**
- Shine J & Dalgarno L (1974). The 3'-terminal sequence of Escherichia coli 16S ribosomal RNA: complementarity to nonsense triplets and ribosome binding sites. *PNAS*, 71(4), 1342–1346. https://doi.org/10.1073/pnas.71.4.1342
- Stormo GD (2000). DNA binding sites: representation and discovery. *Bioinformatics*, 16(1), 16–23. https://doi.org/10.1093/bioinformatics/16.1.16

---

### Analysis 16 — MFE Analysis (5' Structure)

**Input:** One or more `.gbk` files
**Output:** MFE box plots (PNG), violin plots (PNG), interspecies energy density maps (PNG), per-gene MFE table (CSV)
**Parameter:** `--mfe-region` (default: 50 bp)
**Requires:** ViennaRNA (`conda install -c bioconda viennarna`)

Calculates the **Minimum Free Energy (MFE, kcal/mol)** of the predicted mRNA secondary structure of the 5' region (first N nucleotides) of each CDS, using the Zuker-Stiegler thermodynamic folding algorithm via the `RNA.fold()` function from ViennaRNA (Lorenz et al., 2011).

Highly stable secondary structures (very negative MFE) near the ribosome binding site can repress translation initiation by masking the start codon. The **translational ramp hypothesis** proposes that genes with less stable 5' structures (less negative MFE) facilitate ribosome loading and achieve higher expression. Results are compared across species using box plots, violin plots, and interspecies energy density maps.

**References:**
- Lorenz R, et al. (2011). ViennaRNA Package 2.0. *Algorithms for Molecular Biology*, 6, 26. https://doi.org/10.1186/1748-7188-6-26
- Zuker M & Stiegler P (1981). Optimal computer folding of large RNA sequences using thermodynamics and auxiliary information. *Nucleic Acids Research*, 9, 133–148. https://doi.org/10.1093/nar/9.1.133
- Tuller T, et al. (2010). An evolutionarily conserved mechanism for controlling the efficiency of protein translation. *Cell*, 141(2), 344–354. https://doi.org/10.1016/j.cell.2010.03.031
- Kudla G, et al. (2009). Coding-sequence determinants of gene expression in Escherichia coli. *Science*, 324(5924), 255–258. https://doi.org/10.1126/science.1170160

---

### Analysis 17 — Gene Group Comparison

**Input:** One or more `.gbk` files + two gene list files (`--group1`, `--group2`)
**Output:** Box plots with annotated p-values (PNG), KDE distributions (PNG), GC3 vs. ENC projection (PNG), summary statistics

Allows the user to define two subsets of genes and compares ENC, GC3, and CAI distributions between groups. The **Mann-Whitney U test** (non-parametric; Mann & Whitney 1947) is applied to each continuous metric. **Pearson's chi-squared test** (Pearson 1900) is applied to aggregated codon count tables. Results are visualized as box plots with automatically annotated p-values, KDE distributions, and GC3 vs. ENC projections for the subgroups.

Useful for comparing highly expressed vs. lowly expressed genes, horizontally transferred vs. core genome genes, or pathogenicity islands vs. housekeeping genes.

**References:**
- Mann HB & Whitney DR (1947). On a test of whether one of two random variables is stochastically larger than the other. *The Annals of Mathematical Statistics*, 18(1), 50–60.
- Pearson K (1900). Mathematical contributions to the theory of evolution. *Philosophical Transactions A*, 195, 1–47. https://doi.org/10.1098/rsta.1900.0022
- Karlin S & Mrazek J (2000). Predicted highly expressed genes of diverse prokaryotic genomes. *Journal of Bacteriology*, 182(18), 5238–5250. https://doi.org/10.1128/JB.182.18.5238-5250.2000

---

### Analysis 18 — Correlation with Expression

**Input:** One or more `.gbk` files + expression data file (CSV/TSV)
**Output:** Scatter plots with regression lines (PNG), hexbin density plots (PNG), empirical distributions by CUB quartile (PNG), merged data table (CSV)
**Parameters:** `--expr-gene` (gene column), `--expr-val` (expression value column)

Integrates external gene expression data (RNA-Seq TPM/RPKM or microarray intensities) with per-gene CUB metrics (CAI, ENC). Gene identifiers are cross-referenced with GenBank locus_tags. For overlapping genes, **Spearman's rank correlation** (Spearman 1904; Virtanen et al. 2020) is calculated between CAI vs. log10(expression) and ENC vs. log10(expression). ρ coefficients and p-values are reported. Results include seaborn.regplot scatter plots and empirical distributions stratified by codon bias quartiles (Q1 to Q4).

Under the translational selection hypothesis, highly expressed genes are expected to show **higher CAI** and **lower ENC** (stronger bias toward optimal codons).

**References:**
- Spearman C (1904). The proof and measurement of association between two things. *The American Journal of Psychology*, 15(1), 72–101. https://doi.org/10.2307/1412159
- Ikemura T (1982). Correlation between the abundance of yeast transfer RNAs and the occurrence of the respective codons in protein genes. *Journal of Molecular Biology*, 158(4), 573–597. https://doi.org/10.1016/0022-2836(82)90250-9
- Quax TEF, et al. (2015). Codon bias as a means to fine-tune gene expression. *Molecular Cell*, 59(2), 149–161. https://doi.org/10.1016/j.molcel.2015.05.035
- Plotkin JB & Kudla G (2011). Synonymous but not the same: the causes and consequences of codon bias. *Nature Reviews Genetics*, 12, 32–42. https://doi.org/10.1038/nrg2899

---

## Output Architecture

All output files follow a consistent structure designed for scientific traceability:

- **PNG figures:** rendered at 150 DPI using the Matplotlib Agg backend (hardware-independent, suitable for headless servers and HPC clusters).
- **CSV tables:** semicolon-delimited with decimal point notation, ensuring interoperability with R and pandas.
- **Dynamic filenames:** all output files incorporate the genome identifier, ensuring data provenance in automated pipelines.

The GUI provides an integrated image viewer with zoom, pan, and individual export, plus a log console that reports analysis progress, CDS and codon counts, per-genome metric values, and annotation warnings in real time.

---

## Case Studies

Kódon-X was validated through case studies on four bacterial genomes from two clinically relevant species:

| Genome | Species | Source |
|--------|---------|--------|
| SA22AQUAVET | *Streptococcus agalactiae* (strain 22) | AQUAVET Lab, UFMG |
| SA90AQUAVET | *Streptococcus agalactiae* (strain 90) | AQUAVET Lab, UFMG |
| AB_ACTC19606 | *Acinetobacter baumannii* ATCC 19606 | NCBI |
| AB_XH1056 | *Acinetobacter baumannii* XH1056 | NCBI |

**Key findings:**
- *A. baumannii* strains showed greater inter-strain similarity than *S. agalactiae* strains, as evidenced by RSCU correlation regression.
- The Wright Plot (per-gene scatter) revealed distinct patterns of mutational vs. selective pressure between the Gram-positive (*S. agalactiae*) and Gram-negative (*A. baumannii*) organisms.
- Dinucleotide odds ratios (ρ_XY) clearly distinguished genomic signatures between the two species.
- The combined four-genome analysis captured marked divergence in codon usage between phylogenetically distant organisms while maintaining analytical consistency.

---

## Project Structure

```
Kodon-X/
├── KodonX.py                        # Main entry point (GUI)
├── KodonXa.py                       # CLI entry point
├── interface.py                     # PyQt6 GUI implementation
├── constants.py                     # Genetic code tables, physicochemical constants, wobble matrices
├── core_utils.py                    # Core functions: RSCU, ENC, GC3, CAI, per-gene metrics
├── analysis_basic.py                # Analyses 1, 2, 8 (statistics, gene listing, composition)
├── analysis_bias.py                 # Analyses 3–7, 9, 12 (RSCU, ENC/GC3, CAI, neutrality plot)
├── analysis_advanced.py             # Analyses 10, 11, 13–16 (CPB, GRAVY, dinucleotides, PR2, tAI, motifs)
├── analysis_expression_structure.py # Analyses 17–19 (MFE, group comparison, expression correlation)
└── synthetic_biology.py             # Analysis 0 (codon optimization and harmonization)
```

---

## References

1. Angov E, Hillier CJ, Kincaid RL, Lyon JA (2008). Heterologous protein expression is enhanced by harmonizing the codon usage frequencies of the target gene with those of the expression host. *PLoS ONE*, 3(5), e2189. https://doi.org/10.1371/journal.pone.0002189
2. Chaney JL & Clark PL (2015). Roles for synonymous codon usage in protein biogenesis. *Annual Review of Biophysics*, 44, 143–166. https://doi.org/10.1146/annurev-biophys-060414-034333
3. Cock PJA, et al. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics*, 25(11), 1422–1423. https://doi.org/10.1093/bioinformatics/btp163
4. Coleman JR, et al. (2008). Virus attenuation by genome-scale changes in codon pair bias. *Science*, 320(5884), 1784–1787. https://doi.org/10.1126/science.1155761
5. Crick FHC (1966). Codon–anticodon pairing: the wobble hypothesis. *Journal of Molecular Biology*, 19(2), 548–555. https://doi.org/10.1016/S0022-2836(66)80022-0
6. dos Reis M, Savva R & Wernisch L (2004). Solving the riddle of codon usage preferences: a test for translational selection. *Nucleic Acids Research*, 32(17), 5036–5044. https://doi.org/10.1093/nar/gkh834
7. Gustafsson C, et al. (2004). Codon bias and heterologous protein expression. *Trends in Biotechnology*, 22(7), 346–353. https://doi.org/10.1016/j.tibtech.2004.04.006
8. Gutman GA & Hatfield GW (1989). Nonrandom utilization of codon pairs in Escherichia coli. *PNAS*, 86(10), 3699–3703. https://doi.org/10.1073/pnas.86.10.3699
9. Harris CR, et al. (2020). Array programming with NumPy. *Nature*, 585, 357–362. https://doi.org/10.1038/s41586-020-2649-2
10. Hunter J (2007). Matplotlib: A 2D graphics environment. *Computing in Science & Engineering*, 9(3), 90–95.
11. Ikemura T (1981). Correlation between the abundance of Escherichia coli transfer RNAs and the occurrence of the respective codons in its protein genes. *Journal of Molecular Biology*, 146(1), 1–21. https://doi.org/10.1016/0022-2836(81)90363-6
12. Ikemura T (1982). Correlation between the abundance of yeast transfer RNAs and the occurrence of the respective codons in protein genes. *Journal of Molecular Biology*, 158(4), 573–597. https://doi.org/10.1016/0022-2836(82)90250-9
13. Ikemura T (1985). Codon usage and tRNA content in unicellular and multicellular organisms. *Molecular Biology and Evolution*, 2(1), 13–34. https://doi.org/10.1093/oxfordjournals.molbev.a040335
14. Karlin S & Burge C (1995). Dinucleotide relative abundance extremes: a genomic signature. *Trends in Genetics*, 11(7), 283–290. https://doi.org/10.1016/S0168-9525(00)89076-9
15. Karlin S & Mrazek J (1996). What drives codon choices in human genes? *Journal of Molecular Biology*, 262(4), 459–472. https://doi.org/10.1006/jmbi.1996.0528
16. Karlin S & Mrazek J (2000). Predicted highly expressed genes of diverse prokaryotic genomes. *Journal of Bacteriology*, 182(18), 5238–5250. https://doi.org/10.1128/JB.182.18.5238-5250.2000
17. Karlin S, et al. (1997). Compositional differences within and between eukaryotic genomes. *PNAS*, 94(19), 10227–10232. https://doi.org/10.1073/pnas.94.19.10227
18. Kudla G, et al. (2009). Coding-sequence determinants of gene expression in Escherichia coli. *Science*, 324(5924), 255–258. https://doi.org/10.1126/science.1170160
19. Kyte J & Doolittle RF (1982). A simple method for displaying the hydropathic character of a protein. *Journal of Molecular Biology*, 157(1), 105–132. https://doi.org/10.1016/0022-2836(82)90515-0
20. Liu Y (2020). A code within the genetic code: codon usage regulates co-translational protein folding. *Cell Communication and Signaling*, 18, 145. https://doi.org/10.1186/s12964-020-00642-6
21. Lobry JR (1996). Asymmetric substitution patterns in the two DNA strands of bacteria. *Molecular Biology and Evolution*, 13(5), 660–665. https://doi.org/10.1093/oxfordjournals.molbev.a025626
22. Lobry JR & Gautier C (1994). Hydrophobicity, expressivity and aromaticity are the major trends of amino-acid usage in 999 Escherichia coli chromosome-encoded genes. *Nucleic Acids Research*, 22(15), 3174–3180. https://doi.org/10.1093/nar/22.15.3174
23. Lorenz R, et al. (2011). ViennaRNA Package 2.0. *Algorithms for Molecular Biology*, 6, 26. https://doi.org/10.1186/1748-7188-6-26
24. Mann HB & Whitney DR (1947). On a test of whether one of two random variables is stochastically larger than the other. *The Annals of Mathematical Statistics*, 18(1), 50–60.
25. McKinney W (2010). Data structures for statistical computing in Python. *SciPy 2010*. https://doi.org/10.25080/Majora-92bf1922-00a
26. Novembre J (2002). Accounting for background nucleotide composition when measuring codon usage bias. *Molecular Biology and Evolution*, 19(8), 1390–1394. https://doi.org/10.1093/oxfordjournals.molbev.a004201
27. Parvathy ST, Udayasuriyan V & Bhadana V (2022). Codon usage bias. *Molecular Biology Reports*, 49, 539–565. https://doi.org/10.1007/s11033-021-06749-4
28. Pearson K (1900). Mathematical contributions to the theory of evolution. *Philosophical Transactions A*, 195, 1–47. https://doi.org/10.1098/rsta.1900.0022
29. Pechmann S & Frydman J (2013). Evolutionary conservation of codon optimality reveals hidden signatures of cotranslational folding. *Nature Structural & Molecular Biology*, 20(2), 237–243. https://doi.org/10.1038/nsmb.2508
30. Pedregosa F, et al. (2011). Scikit-learn: Machine learning in Python. *Journal of Machine Learning Research*, 12, 2825–2830.
31. Plotkin JB & Kudla G (2011). Synonymous but not the same: the causes and consequences of codon bias. *Nature Reviews Genetics*, 12, 32–42. https://doi.org/10.1038/nrg2899
32. Puigbò P, Bravo IG & Garcia-Vallve S (2008). CAIcal: A combined set of tools to assess codon usage adaptation. *Biology Direct*, 3, 38. https://doi.org/10.1186/1745-6150-3-38
33. Quax TEF, et al. (2015). Codon bias as a means to fine-tune gene expression. *Molecular Cell*, 59(2), 149–161. https://doi.org/10.1016/j.molcel.2015.05.035
34. Sharp PM & Li WH (1986). Codon usage in regulatory genes in Escherichia coli does not reflect selection for 'rare' codons. *Nucleic Acids Research*, 14, 7737–7749. https://doi.org/10.1093/nar/14.19.7737
35. Sharp PM & Li WH (1987). The codon adaptation index — a measure of directional synonymous codon usage bias, and its potential applications. *Nucleic Acids Research*, 15(3), 1281–1295. https://doi.org/10.1093/nar/15.3.1281
36. Sharp PM, et al. (1988). Codon usage in Saccharomyces cerevisiae: evidence for translational selection. *Nucleic Acids Research*, 16(17), 8207–8211. https://doi.org/10.1093/nar/16.17.8207
37. Sharp PM, et al. (2010). Variation in the strength of selected codon usage bias among bacteria. *Nucleic Acids Research*, 38(15), 4941–4950. https://doi.org/10.1093/nar/gkq142
38. Shine J & Dalgarno L (1974). The 3'-terminal sequence of Escherichia coli 16S ribosomal RNA: complementarity to nonsense triplets and ribosome binding sites. *PNAS*, 71(4), 1342–1346. https://doi.org/10.1073/pnas.71.4.1342
39. Sneath PH & Sokal RR (1973). *Numerical Taxonomy: The Principles and Practice of Numerical Classification*. W. H. Freeman, San Francisco.
40. Spearman C (1904). The proof and measurement of association between two things. *The American Journal of Psychology*, 15(1), 72–101. https://doi.org/10.2307/1412159
41. Stormo GD (2000). DNA binding sites: representation and discovery. *Bioinformatics*, 16(1), 16–23. https://doi.org/10.1093/bioinformatics/16.1.16
42. Sueoka N (1962). On the genetic basis of variation and heterogeneity of DNA base composition. *PNAS*, 48(4), 582–592. https://doi.org/10.1073/pnas.48.4.582
43. Sueoka N (1988). Directional mutation pressure and neutral molecular evolution. *PNAS*, 85(8), 2653–2657. https://doi.org/10.1073/pnas.85.8.2653
44. Sueoka N (1995). Intrastrand parity rules of DNA base composition and usage biases of synonymous codons. *Journal of Molecular Evolution*, 40(3), 318–325. https://doi.org/10.1007/BF00163236
45. Sueoka N (1999). Two aspects of DNA base composition: G+C content and translation-coupled deviation from intra-strand rule of A=T and G=C. *Journal of Molecular Evolution*, 49, 49–62. https://doi.org/10.1007/PL00006534
46. Tuller T, et al. (2010). An evolutionarily conserved mechanism for controlling the efficiency of protein translation. *Cell*, 141(2), 344–354. https://doi.org/10.1016/j.cell.2010.03.031
47. Virtanen P, et al. (2020). SciPy 1.0: fundamental algorithms for scientific computing in Python. *Nature Methods*, 17, 261–272. https://doi.org/10.1038/s41592-019-0686-2
48. Waskom ML (2021). seaborn: statistical data visualization. *Journal of Open Source Software*, 6, 3021. https://doi.org/10.21105/joss.03021
49. Wright F (1990). The 'effective number of codons' used in a gene. *Gene*, 87(1), 23–29. https://doi.org/10.1016/0378-1119(90)90491-9
50. Wu G, Bashir-Bello N & Freeland SJ (2006). The Synthetic Gene Designer: A flexible web platform to explore sequence manipulation for heterologous expression. *Protein Expression and Purification*, 47, 441–445. https://doi.org/10.1016/j.pep.2005.10.020
51. Yang Z & Nielsen R (2008). Mutation-selection models of codon substitution and their use to estimate selective strengths on codon usage. *Molecular Biology and Evolution*, 25(3), 568–579. https://doi.org/10.1093/molbev/msm284
52. Zuker M & Stiegler P (1981). Optimal computer folding of large RNA sequences using thermodynamics and auxiliary information. *Nucleic Acids Research*, 9, 133–148. https://doi.org/10.1093/nar/9.1.133

---

*Kódon-X — Developed at the Laboratory of Integrative Bioinformatics (LBI), Federal University of Minas Gerais (UFMG), Brazil.*
