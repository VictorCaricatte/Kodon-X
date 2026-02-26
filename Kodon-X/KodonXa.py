import os
import sys
import glob
import argparse
import pandas as pd

from core_utils import process_genomes_for_bias_analysis
from analysis_basic import process_aggregated_gbk, analyze_gbk_cds, list_genes_from_file, analyze_genomic_composition
from analysis_bias import generate_rscu_heatmap_and_table, comparative_rscu_analysis, rscu_correlation_analysis, generate_rscu_histograms, enc_gc3_analysis, optimal_rare_codons_analysis, neutrality_plot_analysis
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
Kodon-X CLI - Comprehensive Codon Usage and Genomic Analysis Tool
================================================================================

This command-line interface allows you to run all 19 analyses and the synthetic 
biology module without a graphical interface. Perfect for servers and pipelines.

AVAILABLE ANALYSES (-a / --analysis):

[0] Synthetic Biology
    Optimizes or harmonizes a DNA sequence for a specific host genome.
    Requires: --synth-mode, --synth-host, --synth-seq

[1] Statistics and CDS
    Performs an initial scan on all genomes. Validates genome completeness 
    and checks CDS consistency (e.g., start codons).

[2] Gene Listing
    Lists all genes, products, and locus tags present in a single GenBank file.
    Requires exactly 1 input file.

[3] Individual RSCU Heatmap
    Calculates and visualizes the RSCU for a single genome. Identifies preferred 
    and rare codons, calculating global metrics (ENC, GC3, CAI).
    Requires exactly 1 input file.

[4] Comparative RSCU
    Groups genomes based on the similarity of their codon bias. Generates a 
    clustermap and PCA to separate organisms by codon usage.
    Requires 2 or more files.

[5] RSCU Correlation (2 Genomes)
    Calculates the Pearson correlation (R) between the RSCU patterns of exactly 
    two genomes to measure evolutionary proximity.
    Requires exactly 2 files.

[6] Comparative RSCU Histograms
    Generates detailed visualizations of codon usage, including Box Plots by 
    amino acid and stacked bar charts by species.
    Requires 2 or more files.

[7] ENC vs GC3 Analysis (Wright Plot)
    Generates the Wright Plot to determine if codon bias is driven by natural 
    selection or mutational pressure.

[8] Genomic Composition
    Detailed analysis of nucleotide composition (A, T, G, C), total GC content, 
    and genome size for multiple files.

[9] Optimal, Rare Codons and CAI
    Identifies preferred/rare codons and compares the Codon Adaptation Index (CAI) 
    between species.

[10] Codon Pair Bias (CPB)
     Analyzes the frequency of adjacent codon pairs to identify translation bottlenecks.

[11] Physicochemical Analysis (GRAVY & Aromo)
     Calculates the GRAVY score (hydropathicity) and Aromo (aromaticity) for all genes.

[12] Neutrality Plot (GC12 vs GC3)
     Plots the GC content of positions 1+2 against GC3 to separate mutation pressure 
     from selection.
     Requires 2 or more files.

[13] Dinucleotide Composition
     Calculates and compares the frequency of all 16 dinucleotides in the genome.

[14] PR2 Parity Plot
     Plots A3/(A3+T3) vs G3/(G3+C3) for each individual gene to analyze mutation biases.

[15] tRNA Adaptation Index (tAI)
     Calculates tAI using wobble rules to measure theoretical translation efficiency.
     Optional argument: --superkingdom

[16] Upstream Motifs Analysis
     Searches for overrepresented short sequences (k-mers) in regions before the 
     start of each gene.
     Optional arguments: --upstream-dist, --kmer-size

[17] MFE Analysis (5' Structure)
     Calculates the Minimum Free Energy (MFE) of the 5' region to test the 
     translational ramp hypothesis. (Requires ViennaRNA)
     Optional argument: --mfe-region

[18] Gene Group Comparison
     Compares CUB metrics (ENC, CAI) between two user-defined gene sets.
     Requires: --group1, --group2

[19] Correlation with Expression
     Correlates bias metrics with external expression data (e.g., RNA-Seq TPM).
     Requires: --expr-file, --expr-gene, --expr-val
================================================================================
"""

    parser = argparse.ArgumentParser(
        description=help_description,
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    parser.add_argument("-i", "--input", help="Input directory containing .gbk/.gb files")
    parser.add_argument("-o", "--output", help="Output directory for results (tables and charts)")
    parser.add_argument("-a", "--analysis", type=int, help="Analysis ID to execute (0 to 19)")
    parser.add_argument("-g", "--genetic-code", type=int, default=1, help="Genetic Code ID (Default: 1)")
    parser.add_argument("-f", "--filter-file", help="Text file with locus_tags for global filtering (one per line)")
    parser.add_argument("--start-codon", default="ATG", help="Start codon filter for Analysis 1")
    parser.add_argument("--superkingdom", default="Bacteria", choices=["Bacteria", "Eukaryote"], help="Superkingdom for Wobble rules (Analysis 15)")
    parser.add_argument("--upstream-dist", type=int, default=200, help="Upstream distance in bp (Analysis 16)")
    parser.add_argument("--kmer-size", type=int, default=6, help="K-mer size (Analysis 16)")
    parser.add_argument("--mfe-region", type=int, default=50, help="5' region size for MFE calculation (Analysis 17)")
    parser.add_argument("--group1", help="Text file with genes for Group 1 (Analysis 18)")
    parser.add_argument("--group2", help="Text file with genes for Group 2 (Analysis 18)")
    parser.add_argument("--expr-file", help="CSV/TSV file with expression data (Analysis 19)")
    parser.add_argument("--expr-gene", default="locus_tag", help="Gene column name in expression file (Analysis 19)")
    parser.add_argument("--expr-val", default="TPM", help="Value column name in expression file (Analysis 19)")
    parser.add_argument("--synth-mode", choices=["optimize", "harmonize"], default="optimize", help="Mode for synthetic biology (Analysis 0)")
    parser.add_argument("--synth-host", help="Host genome filename inside input dir (Analysis 0)")
    parser.add_argument("--synth-seq", help="Input DNA sequence (Analysis 0)")

    args = parser.parse_args()

    if args.analysis is None:
        parser.print_help()
        sys.exit(1)

    status_queue = DummyQueue()
    
    if args.analysis == 0:
        if not args.input or not args.synth_host or not args.synth_seq:
            print("[!] For Synthetic Biology, provide: --input, --synth-host and --synth-seq")
            sys.exit(1)
            
        host_filepath = os.path.join(args.input, args.synth_host)
        if not os.path.isfile(host_filepath):
            print(f"[!] Host file not found: {host_filepath}")
            sys.exit(1)
            
        print(f"STARTING SYNTHETIC BIOLOGY ({args.synth_mode.upper()})")
        print(f"Host: {args.synth_host}")
        
        all_data = process_genomes_for_bias_analysis([host_filepath], args.genetic_code, status_queue, gene_list=None)
        if not all_data:
            print("[!] Failed to process host data.")
            sys.exit(1)
            
        host_data = all_data[list(all_data.keys())[0]]
        
        if args.synth_mode == "optimize":
            result = optimize_codon_sequence(args.synth_seq, host_data['rscu'], args.genetic_code)
        else:
            result = harmonize_codon_sequence(args.synth_seq, host_data['counts'], args.genetic_code)
            
        print(f"\n[+] Resulting Sequence:\n{result}")
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
    if gene_list:
        print(f"Global Filter Active: {len(gene_list)} genes.")
    print("=" * 60)

    bias_dependent_analyses = [3, 4, 5, 6, 7, 9, 12]
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
            generate_rscu_heatmap_and_table(all_bias_data, args.output, genetic_code_id, status_queue)

        elif args.analysis == 4:
            if len(files) < 2:
                print("[!] Analysis 4 requires 2 or more files.")
                sys.exit(1)
            comparative_rscu_analysis(all_bias_data, args.output, status_queue)

        elif args.analysis == 5:
            if len(files) != 2:
                print("[!] Analysis 5 requires EXACTLY 2 files.")
                sys.exit(1)
            rscu_correlation_analysis(all_bias_data, args.output, status_queue)

        elif args.analysis == 6:
            if len(files) < 2:
                print("[!] Analysis 6 requires 2 or more files.")
                sys.exit(1)
            generate_rscu_histograms(all_bias_data, args.output, genetic_code_id, status_queue)

        elif args.analysis == 7:
            enc_gc3_analysis(all_bias_data, args.output, status_queue)

        elif args.analysis == 8:
            analyze_genomic_composition(files, args.output, status_queue)

        elif args.analysis == 9:
            optimal_rare_codons_analysis(all_bias_data, args.output, status_queue)

        elif args.analysis == 10:
            codon_pair_bias_analysis(files, args.output, genetic_code_id, status_queue, gene_list)

        elif args.analysis == 11:
            gravy_aromo_analysis(files, args.output, genetic_code_id, status_queue, gene_list)

        elif args.analysis == 12:
            if len(files) < 2:
                print("[!] Analysis 12 requires 2 or more files.")
                sys.exit(1)
            neutrality_plot_analysis(all_bias_data, args.output, status_queue)

        elif args.analysis == 13:
            dinucleotide_composition_analysis(files, args.output, status_queue)

        elif args.analysis == 14:
            pr2_plot_analysis(files, args.output, genetic_code_id, status_queue, gene_list)

        elif args.analysis == 15:
            tai_analysis(files, args.output, genetic_code_id, status_queue, gene_list, super_kingdom=args.superkingdom)

        elif args.analysis == 16:
            upstream_motifs_analysis(files, args.output, status_queue, gene_list, args.upstream_dist, args.kmer_size)

        elif args.analysis == 17:
            initiation_mfe_analysis(files, args.output, genetic_code_id, status_queue, gene_list, args.mfe_region)

        elif args.analysis == 18:
            list_1 = read_gene_list(args.group1)
            list_2 = read_gene_list(args.group2)
            if not list_1 or not list_2:
                print("[!] For Analysis 18, it is mandatory to provide --group1 and --group2.")
                sys.exit(1)
            two_groups_comparative_analysis(files, args.output, genetic_code_id, status_queue, list_1, list_2)

        elif args.analysis == 19:
            if not args.expr_file:
                print("[!] For Analysis 19, provide the CSV/TSV file using --expr-file.")
                sys.exit(1)
            try:
                with open(args.expr_file, 'r') as f:
                    first_line = f.readline()
                    if '\t' in first_line: sep = '\t'
                    elif ';' in first_line: sep = ';'
                    else: sep = ','
                df_expr = pd.read_csv(args.expr_file, sep=sep)
            except Exception as e:
                print(f"[!] Error reading expression file: {e}")
                sys.exit(1)
            
            expression_correlation_analysis(
                files, args.output, genetic_code_id, status_queue, gene_list, 
                expression_data=df_expr, gene_col=args.expr_gene, expr_col=args.expr_val
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