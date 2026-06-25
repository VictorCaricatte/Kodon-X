import os
import sys
import glob
import argparse
import pandas as pd

from core_utils import process_genomes_for_bias_analysis
from analysis_basic import process_aggregated_gbk, analyze_gbk_cds, list_genes_from_file, analyze_genomic_composition
from analysis_bias import generate_rscu_heatmap_and_table, comparative_rscu_analysis, rscu_correlation_analysis, enc_gc3_analysis, optimal_rare_codons_analysis, neutrality_plot_analysis
from analysis_advanced import codon_pair_bias_analysis, gravy_aromo_analysis, dinucleotide_composition_analysis, pr2_plot_analysis, tai_analysis, upstream_motifs_analysis
from analysis_expression_structure import initiation_mfe_analysis, two_groups_comparative_analysis, expression_correlation_analysis
from synthetic_biology import optimize_codon_sequence, harmonize_codon_sequence

class DummyQueue:
    def put(self, item):
        command, data = item
        if command == "message":
            print(f"[*] {data}")
        elif command == "progress":
            sys.stdout.write(f"\r[+] Progress: {data}% ")
            sys.stdout.flush()
            if data >= 100:
                print()
        elif command == "image_ready":
            print(f"\n[+] Image saved: {data[0]} ({data[1]})")
        elif command == "error":
            print(f"\n[!] ERROR: {data}")
        elif command == "optimization_complete":
            print(f"\n[+] Resulting Sequence:\n{data}")
        elif command in ["done", "tool_done"]:
            pass

def get_files_from_dir(directory):
    patterns = ["*.gbk", "*.gb", "*.gbff"]
    files = []
    for pattern in patterns:
        files.extend(glob.glob(os.path.join(directory, pattern)))
    return sorted(files)

def read_gene_list(filepath):
    if not filepath or not os.path.isfile(filepath):
        return None
    with open(filepath, 'r') as f:
        return set(line.strip() for line in f if line.strip())

def main():
    help_description = """
================================================================================
Kodon-X - Comprehensive Codon Usage Bias Analysis
================================================================================

Run all 18 analyses and the synthetic biology module without a graphical interface.

AVAILABLE ANALYSES (-a / --analysis):

[0]  Synthetic Biology
     Optimizes or harmonizes a DNA sequence for a specific host genome.
     Requires: --synth-mode, --synth-host, and (--synth-seq OR --synth-fasta)

[1]  Statistics and CDS
     Initial scan: genome completeness, CDS consistency, start codon distribution.

[2]  Gene Listing
     Lists all genes, products, and locus tags in a single GenBank file.
     Requires exactly 1 file.

[3]  Individual RSCU Heatmap
     RSCU heatmap + table for a single genome. Global metrics: ENC, GC3, CAI.
     Requires exactly 1 file.

[4]  Comparative RSCU  [2+ files]
     Unified comparative RSCU module — 6 figures:
       • (1) Hierarchical clustermap (64 codons × N species)
       • (2) PCA biplot (species scores + top-20 codon loadings)
       • (3) Comparative line plot per codon (grouped by amino acid)
       • (4) Box plot by amino acid
       • (5) Top-30 most variable codons — variance barplot
       • (6) RSCU heatmap — top-30 most divergent codons

[5]  RSCU Correlation  [exactly 2 files]
     Pearson correlation between RSCU vectors of two genomes.
     Outputs: regression scatter, delta-RSCU barplot, top-20 divergent codons heatmap.

[6]  ENC vs GC3 — Wright Plot  [1+ files]
     Per-gene scatter: ENC vs GC3 with Wright (1990) expected curve overlaid.
     Genes below the curve are under translational selection.

[7]  Genomic Composition  [1+ files]
     Nucleotide frequencies (A, T, G, C), GC content, and genome size.

[8]  Optimal, Rare Codons and CAI  [1+ files]
     Preferred/rare codon identification; CAI comparison across species.

[9]  Codon Pair Bias (CPB)  [1+ files]
     Adjacent codon pair frequency matrix; codon pair scores (CPS).

[10] Physicochemical Analysis — GRAVY & Aromo  [1+ files]
     GRAVY (Kyte-Doolittle hydropathicity) and Aromaticity per gene.

[11] Neutrality Plot — GC12 vs GC3  [2+ files]
     Per-gene OLS regression (Sueoka 1992) per species.
     Slope ≈ 1: mutation pressure; slope ≈ 0: translational selection.

[12] Dinucleotide Composition ρ_XY  [1+ files]
     Karlin-Burge (1995) odds ratio ρ_XY = f(XY)/(f(X)×f(Y)).
     Diverging heatmap centred at ρ = 1.0 (CpG suppression detectable).

[13] PR2 Parity Plot — A3/T3 vs G3/C3  [1+ files]
     Sueoka (1995) strand-bias plot per gene.

[14] tRNA Adaptation Index (tAI)  [1+ files]
     Wobble-weighted tAI (dos Reis 2004). Optional: --superkingdom

[15] Upstream Motifs Analysis  [1+ files]
     Over-represented k-mers in upstream regions.
     Optional: --upstream-dist, --kmer-size

[16] MFE Analysis — 5' Structure  [1+ files]
     Minimum Free Energy of 5' region (translational ramp). Requires ViennaRNA.
     Optional: --mfe-region

[17] Gene Group Comparison  [1+ files]
     Mann-Whitney U + Chi-squared between two user-defined gene sets.
     Requires: --group1, --group2

[18] Correlation with Expression (RNA-Seq)  [1+ files]
     Spearman correlation of CAI/ENC with external expression data (TPM/RPKM).
     Requires: --expr-file, --expr-gene, --expr-val
================================================================================
"""

    parser = argparse.ArgumentParser(
        description=help_description,
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    parser.add_argument("-i", "--input", help="Input directory containing .gbk/.gb files")
    parser.add_argument("-o", "--output", help="Output directory for results (tables and charts)")
    parser.add_argument("-a", "--analysis", type=int, help="Analysis ID to execute (0 to 18)")
    parser.add_argument("-g", "--genetic-code", type=int, default=1, help="Genetic Code ID (Default: 1 — Standard)")
    parser.add_argument("-f", "--filter-file", help="Text file with locus_tags for global gene filtering (one per line)")
    parser.add_argument("--palette", default="viridis", help="Chart color palette (e.g. viridis, plasma, Set2, Dark2)")
    parser.add_argument("--start-codon", default="ATG", help="Start codon filter for Analysis 1")
    parser.add_argument("--superkingdom", default="Bacteria", choices=["Bacteria", "Eukaryote"], help="Superkingdom for wobble s-values — Analysis 14")
    parser.add_argument("--upstream-dist", type=int, default=200, help="Upstream region size in bp — Analysis 15")
    parser.add_argument("--kmer-size", type=int, default=6, help="K-mer length for upstream motif search — Analysis 15")
    parser.add_argument("--mfe-region", type=int, default=50, help="5' region size in bp for MFE calculation — Analysis 16")
    parser.add_argument("--group1", help="Text file with gene locus_tags for Group 1 — Analysis 17")
    parser.add_argument("--group2", help="Text file with gene locus_tags for Group 2 — Analysis 17")
    parser.add_argument("--expr-file", help="CSV/TSV file with expression data — Analysis 18")
    parser.add_argument("--expr-gene", default="locus_tag", help="Gene identifier column name in expression file — Analysis 18")
    parser.add_argument("--expr-val", default="TPM", help="Expression value column name (TPM, RPKM, etc.) — Analysis 18")
    parser.add_argument("--synth-mode", choices=["optimize", "harmonize"], default="optimize", help="Mode for synthetic biology (Analysis 0)")
    parser.add_argument("--synth-host", help="Host genome filename inside input dir (Analysis 0)")
    parser.add_argument("--synth-seq", help="Input DNA sequence as string (Analysis 0)")
    parser.add_argument("--synth-fasta", help="Input FASTA file for Synthetic Biology (Analysis 0)")

    args = parser.parse_args()

    if args.analysis is None:
        parser.print_help()
        sys.exit(1)

    status_queue = DummyQueue()
    
    if args.analysis == 0:
        if not args.input or not args.synth_host:
            print("[!] For Synthetic Biology, provide: --input and --synth-host")
            sys.exit(1)
        if not args.synth_seq and not args.synth_fasta:
            print("[!] For Synthetic Biology, provide --synth-seq OR --synth-fasta")
            sys.exit(1)
            
        host_filepath = os.path.join(args.input, args.synth_host)
        if not os.path.isfile(host_filepath):
            print(f"[!] Host file not found: {host_filepath}")
            sys.exit(1)
            
        seq = ""
        if args.synth_fasta:
            from Bio import SeqIO
            try:
                record = next(SeqIO.parse(args.synth_fasta, "fasta"))
                seq = str(record.seq)
            except Exception as e:
                print(f"[!] Could not read FASTA file: {e}")
                sys.exit(1)
        else:
            seq = args.synth_seq

        print(f"STARTING SYNTHETIC BIOLOGY ({args.synth_mode.upper()})")
        print(f"Host: {args.synth_host}")
        
        all_data = process_genomes_for_bias_analysis([host_filepath], args.genetic_code, status_queue, gene_list=None)
        if not all_data:
            print("[!] Failed to process host data.")
            sys.exit(1)
            
        host_data = all_data[list(all_data.keys())[0]]
        
        if args.synth_mode == "optimize":
            result = optimize_codon_sequence(seq, host_data['rscu'], args.genetic_code)
        else:
            result = harmonize_codon_sequence(seq, host_data['counts'], args.genetic_code)
            
        print(f"\n[+] Resulting Sequence:\n{result}")

        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        from Bio import SeqIO
        
        record = SeqRecord(Seq(result), id=f"Synth_{args.synth_mode}", name="SyntheticBio", description=f"Sequence {args.synth_mode}d for {args.synth_host}")
        record.annotations["molecule_type"] = "DNA"
        
        feature = SeqFeature(FeatureLocation(start=0, end=len(result)), type="CDS")
        record.features.append(feature)
        
        if args.output and not os.path.isdir(args.output):
            os.makedirs(args.output, exist_ok=True)
            
        out_path = os.path.join(args.output, f"synth_{args.synth_mode}_result.gbk") if args.output else f"synth_{args.synth_mode}_result.gbk"
        SeqIO.write(record, out_path, "genbank")
        print(f"[+] Saved GenBank file to: {out_path}")
        sys.exit(0)

    if not args.input or not args.output:
        print("[!] Please provide input (-i) and output (-o) directories.")
        sys.exit(1)

    if not os.path.isdir(args.output):
        os.makedirs(args.output, exist_ok=True)

    files = get_files_from_dir(args.input)
    if not files:
        print(f"[!] No .gbk/.gb files found in: {args.input}")
        sys.exit(1)

    gene_list = read_gene_list(args.filter_file)
    genetic_code_id = args.genetic_code

    print(f"STARTING ANALYSIS {args.analysis}")
    print(f"Detected files: {len(files)}")
    print(f"Genetic Code: {genetic_code_id}")
    print(f"Palette: {args.palette}")
    if gene_list:
        print(f"Global Filter Active: {len(gene_list)} genes.")
    print("=" * 60)

    bias_dependent_analyses = [3, 4, 5, 6, 8, 11]
    all_bias_data = None
    if args.analysis in bias_dependent_analyses:
        all_bias_data = process_genomes_for_bias_analysis(files, genetic_code_id, status_queue, gene_list)
        if not all_bias_data:
            print("[!] Failed to process bias data. Check input files.")
            sys.exit(1)

    try:
        if args.analysis == 1:
            print("\n--- [Analysis 1] Part 1: Aggregated Statistics ---")
            df_stats = process_aggregated_gbk(files, status_queue)
            if not df_stats.empty:
                df_stats.to_csv(os.path.join(args.output, 'genome_statistics.csv'), index=False, sep=';')
            print("\n--- [Analysis 1] Part 2: CDS Analysis ---")
            df_cds = analyze_gbk_cds(files, args.start_codon, status_queue, gene_list)
            if not df_cds.empty:
                df_cds.to_csv(os.path.join(args.output, 'cds_analysis.csv'), index=False, sep=';')

        elif args.analysis == 2:
            if len(files) != 1:
                print("[!] Analysis 2 requires EXACTLY 1 file.")
                sys.exit(1)
            df_genes = list_genes_from_file(files[0], status_queue, gene_list)
            if not df_genes.empty:
                df_genes.to_csv(os.path.join(args.output, f"gene_list_{os.path.basename(files[0])}.csv"), index=False, sep=';')

        elif args.analysis == 3:
            if len(files) != 1:
                print("[!] Analysis 3 requires EXACTLY 1 file.")
                sys.exit(1)
            generate_rscu_heatmap_and_table(all_bias_data, args.output, genetic_code_id, status_queue, palette=args.palette)

        elif args.analysis == 4:
            if len(files) < 2:
                print("[!] Analysis 4 requires 2 or more files.")
                sys.exit(1)
            comparative_rscu_analysis(all_bias_data, args.output, status_queue, palette=args.palette)

        elif args.analysis == 5:
            if len(files) != 2:
                print("[!] Analysis 5 requires EXACTLY 2 files.")
                sys.exit(1)
            rscu_correlation_analysis(all_bias_data, args.output, status_queue, palette=args.palette)

        elif args.analysis == 6:
            enc_gc3_analysis(all_bias_data, args.output, status_queue,
                             file_list=files, genetic_code_id=genetic_code_id,
                             gene_list=gene_list, palette=args.palette)

        elif args.analysis == 7:
            analyze_genomic_composition(files, args.output, status_queue, palette=args.palette)

        elif args.analysis == 8:
            optimal_rare_codons_analysis(all_bias_data, args.output, status_queue, palette=args.palette)

        elif args.analysis == 9:
            codon_pair_bias_analysis(files, args.output, genetic_code_id, status_queue, gene_list, palette=args.palette)

        elif args.analysis == 10:
            gravy_aromo_analysis(files, args.output, genetic_code_id, status_queue, gene_list, palette=args.palette)

        elif args.analysis == 11:
            if len(files) < 2:
                print("[!] Analysis 11 requires 2 or more files.")
                sys.exit(1)
            neutrality_plot_analysis(all_bias_data, args.output, status_queue,
                                     file_list=files, genetic_code_id=genetic_code_id,
                                     gene_list=gene_list, palette=args.palette)

        elif args.analysis == 12:
            dinucleotide_composition_analysis(files, args.output, status_queue, palette=args.palette)

        elif args.analysis == 13:
            pr2_plot_analysis(files, args.output, genetic_code_id, status_queue, gene_list, palette=args.palette)

        elif args.analysis == 14:
            tai_analysis(files, args.output, genetic_code_id, status_queue, gene_list,
                         super_kingdom=args.superkingdom, palette=args.palette)

        elif args.analysis == 15:
            upstream_motifs_analysis(files, args.output, status_queue, gene_list,
                                     args.upstream_dist, args.kmer_size, palette=args.palette)

        elif args.analysis == 16:
            initiation_mfe_analysis(files, args.output, genetic_code_id, status_queue, gene_list,
                                    args.mfe_region, palette=args.palette)

        elif args.analysis == 17:
            list_1 = read_gene_list(args.group1)
            list_2 = read_gene_list(args.group2)
            if not list_1 or not list_2:
                print("[!] Analysis 17 requires --group1 and --group2.")
                sys.exit(1)
            two_groups_comparative_analysis(files, args.output, genetic_code_id, status_queue,
                                            list_1, list_2, palette=args.palette)

        elif args.analysis == 18:
            if not args.expr_file:
                print("[!] Analysis 18 requires --expr-file.")
                sys.exit(1)
            try:
                with open(args.expr_file, 'r') as f:
                    first_line = f.readline()
                    sep = '\t' if '\t' in first_line else (';' if ';' in first_line else ',')
                df_expr = pd.read_csv(args.expr_file, sep=sep)
            except Exception as e:
                print(f"[!] Error reading expression file: {e}")
                sys.exit(1)
            expression_correlation_analysis(
                files, args.output, genetic_code_id, status_queue, gene_list,
                expression_data=df_expr, gene_col=args.expr_gene, expr_col=args.expr_val,
                palette=args.palette
            )

        else:
            print(f"[!] Analysis ID '{args.analysis}' is invalid. Choose from 1 to 19 (or 0).")

        print("=" * 60)
        print(f"ANALYSIS {args.analysis} COMPLETED.")
        print(f"Results saved in: {os.path.abspath(args.output)}")
        print("=" * 60)

    except Exception as e:
        print(f"\n[!] ERROR DURING EXECUTION: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
