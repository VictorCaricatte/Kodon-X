import os
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
from Bio import SeqIO
from core_utils import _apply_gene_filter

def process_aggregated_gbk(file_list, status_queue):
    if not file_list:
        print("  ❌ No .gbk or .gb file provided.")
        return pd.DataFrame()
    
    status_queue.put(("progress", 20))
    general_results = []
    
    for i, full_path in enumerate(file_list):
        filename = os.path.basename(full_path)
        print(f"  Processing statistics of {filename}...")
        try:
            records = list(SeqIO.parse(full_path, "genbank"))
            num_records = len(records)
            if num_records == 0:
                print(f"  ❌ Error: File {filename} is empty or corrupted.")
                continue
            
            total_length = 0
            total_gc_count = 0
            for record in records:
                seq = record.seq.upper()
                total_length += len(seq)
                total_gc_count += seq.count('G') + seq.count('C')
                
            if total_length == 0:
                total_gc_content = 0.0
            else:
                total_gc_content = (total_gc_count / total_length) * 100
                
            general_results.append({
                'File': filename, 'Total_Contigs': num_records,
                'Total_Length': total_length, 'Total_GC_%': f'{total_gc_content:.2f}'
            })
        except Exception as e:
            print(f"  ❌ Error processing file {filename}: {e}")
            
    status_queue.put(("progress", 50))
    return pd.DataFrame(general_results)

def analyze_gbk_cds(file_list, filter_string, status_queue, gene_list=None):
    if not file_list:
        print("  ❌ No .gbk or .gb file provided.")
        return pd.DataFrame()
        
    final_results = []
    total_files = len(file_list)
    
    for i, full_path in enumerate(file_list):
        filename = os.path.basename(full_path)
        print(f"  Analyzing CDS of {filename}...")
        try:
            records = list(SeqIO.parse(full_path, "genbank"))
            if not records:
                print(f"  ❌ Error: File {filename} is empty or corrupted.")
                continue

            for record in records:
                for feature in record.features:
                    if feature.type != "CDS": continue
                    
                    if not _apply_gene_filter(feature, gene_list):
                        continue
                    
                    identifier = "N/A"
                    if "locus_tag" in feature.qualifiers: identifier = feature.qualifiers["locus_tag"][0]
                    elif "protein_id" in feature.qualifiers: identifier = feature.qualifiers["protein_id"][0]
                    elif "gene" in feature.qualifiers: identifier = feature.qualifiers["gene"][0]
                    
                    seq_cds = feature.extract(record.seq)
                    if feature.location.strand == -1: seq_cds = seq_cds.reverse_complement()
                    
                    codon_start = int(feature.qualifiers.get("codon_start", ["1"])[0])
                    seq_cds = seq_cds[codon_start - 1:]
                    real_size = (len(seq_cds) // 3) * 3
                    seq_cds = seq_cds[:real_size]
                    
                    if len(seq_cds) < 3: continue
                    
                    start_codon = str(seq_cds[:3]).upper()
                    starts_with_filter = start_codon == filter_string.upper()
                    
                    final_results.append({
                        'File': filename, 'Record': record.id, 'Identifier': identifier,
                        'Start_Codon': start_codon, 'Starts_With_Filter': starts_with_filter,
                        'Codon_Start_Qualifier': codon_start, 'Real_ORF_Size_nt': real_size
                    })
        except Exception as e:
            print(f"  ❌ Error analyzing CDS from file {filename}: {e}")
        
        status_queue.put(("progress", int(50 + (i / total_files) * 50)))
        
    return pd.DataFrame(final_results)

def list_genes_from_file(file_path, status_queue, gene_list=None):
    filename = os.path.basename(file_path)
    print(f"\nListing all genes/products of {filename}...")
    results = []
    
    status_queue.put(("progress", 20))
    
    try:
        for i, record in enumerate(SeqIO.parse(file_path, "genbank")):
            if i % 10 == 0: 
                 status_queue.put(("message", f"Processing contig {record.id}..."))
                 
            for feature in record.features:
                if feature.type not in ["CDS", "tRNA", "rRNA", "gene"]:
                    continue
                
                if not _apply_gene_filter(feature, gene_list):
                    continue
                
                gene_name = feature.qualifiers.get("gene", [""])[0]
                product = feature.qualifiers.get("product", [""])[0]
                locus_tag = feature.qualifiers.get("locus_tag", ["N/A"])[0]
                
                if gene_name or product or locus_tag != "N/A":
                    results.append({
                        "File": filename, "Contig": record.id,
                        "Locus_Tag": locus_tag, "Gene": gene_name,
                        "Product": product, "Type": feature.type
                    })
                    
        status_queue.put(("progress", 80))
        
        if not results:
            print(f"  ❌ No genes or products found in {filename}.")
            return pd.DataFrame()
            
        print(f"  ✅ Found {len(results)} total genes/products.")
        return pd.DataFrame(results)
        
    except Exception as e:
        print(f"  ❌ Error processing file {file_path} to list genes: {e}")
        return pd.DataFrame()

def analyze_nucleotide_composition(file_path):
    print(f"  Analyzing nucleotide composition of {os.path.basename(file_path)}...")
    nt_counts = Counter()
    total_length = 0
    
    try:
        for record in SeqIO.parse(file_path, "genbank"):
            seq = record.seq.upper()
            total_length += len(seq)
            for nt in seq:
                if nt in 'ATGC':
                    nt_counts[nt] += 1
        
        if total_length == 0:
            return None
        
        composition = {}
        for nt in 'ATGC':
            composition[nt] = (nt_counts[nt] / total_length) * 100
        
        gc_content = (composition['G'] + composition['C'])
        
        print(f"  ✅ Composition analyzed: GC={gc_content:.2f}%")
        
        return {
            'composition': composition,
            'gc_content': gc_content,
            'total_length': total_length
        }
    except Exception as e:
        print(f"  ❌ Error analyzing composition: {e}")
        return None

def analyze_genomic_composition(file_list, output_folder, status_queue, palette='viridis'):
    print(f"\n=== GENOMIC COMPOSITION ANALYSIS ===")
    all_data = {}
    
    for i, full_path in enumerate(file_list):
        base_name = os.path.basename(full_path).split('.')[0]
        print(f"\nProcessing file {i+1}/{len(file_list)}: {base_name}")
        status_queue.put(("progress", int(20 + (i / len(file_list)) * 70)))
        
        comp_data = analyze_nucleotide_composition(full_path)
        if comp_data:
            all_data[base_name] = comp_data
            print(f"  📊 {base_name}: GC={comp_data['gc_content']:.2f}%, Size={comp_data['total_length']} bp")
    
    if not all_data:
        print("  ❌ Error: No data processed.")
        return
    
    status_queue.put(("progress", 90))
    
    fig, ((ax1, ax2)) = plt.subplots(1, 2, figsize=(16, 8))
    species = list(all_data.keys())
    
    nt_compositions = [all_data[s]['composition'] for s in species]
    df_nt = pd.DataFrame(nt_compositions, index=species)
    df_nt.plot(kind='barh', stacked=False, ax=ax1, colormap=palette)
    
    ax1.set_title('Nucleotide Composition (%)', fontweight='bold')
    ax1.set_xlabel('Percentage (%)') 
    ax1.legend(title='Nucleotides', bbox_to_anchor=(1.0, 1), loc='upper left')
    ax1.grid(axis='x', alpha=0.3)
    
    gc_contents = [all_data[s]['gc_content'] for s in species]
    genome_sizes = [all_data[s]['total_length'] for s in species]

    scatter = ax2.scatter(gc_contents, genome_sizes, c=range(len(species)), 
                          cmap=palette, s=100, alpha=0.7, edgecolors='black')
                          
    ax2.set_title('Genome Size vs GC Content', fontweight='bold')
    ax2.set_xlabel('GC (%)')
    ax2.set_ylabel('Size (bp)')
    ax2.grid(True, alpha=0.3)
    
    for i, species_name in enumerate(species):
        ax2.annotate(species_name, (gc_contents[i], genome_sizes[i]), 
                    xytext=(5, 5), textcoords='offset points', fontsize=9)
    
    ax2.ticklabel_format(style='scientific', axis='y', scilimits=(0,0))
    
    plt.suptitle('Genomic Composition Analysis', fontsize=16, fontweight='bold')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    
    output_path = os.path.join(output_folder, 'genomic_composition.png')
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"\n✅ Genomic composition analysis saved in: {output_path}")
    status_queue.put(("image_ready", (output_path, "Genomic Composition")))
    
    df_results = pd.DataFrame({
        'Species': species,
        'Size_bp': genome_sizes,
        'GC_percent': gc_contents,
        'A_percent': [all_data[s]['composition']['A'] for s in species],
        'T_percent': [all_data[s]['composition']['T'] for s in species],
        'G_percent': [all_data[s]['composition']['G'] for s in species],
        'C_percent': [all_data[s]['composition']['C'] for s in species],
    })
    df_results.to_csv(os.path.join(output_folder, 'genomic_composition_results.csv'), sep=';', index=False)
