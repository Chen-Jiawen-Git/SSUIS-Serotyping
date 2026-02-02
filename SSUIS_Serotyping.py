import os
import subprocess
import csv
import glob
import tempfile
import shutil
import argparse
from Bio import SeqIO

# Configuration (Defaults, can be overridden by arguments)
DEFAULT_OUTPUT_FILE = "serotyping_results_v1.csv"
DEFAULT_MARKER_DB = "marker_genes.fasta" 
DEFAULT_MIN_COVERAGE = 95.0
DEFAULT_MIN_IDENTITY = 95.0

def prepare_marker_db(ref_dir, temp_dir):
    """
    Combines specific marker gene files from ref_dir into a single FASTA.
    Returns path to the combined file.
    """
    marker_files = [
        "Serotype_cps_2_cpsK.fna",
        "Serotype_cps_0.5_cpsK.fna",
        "Marker_1_14_cpsK.fna"
    ]
    
    combined_path = os.path.join(temp_dir, "combined_markers.fasta")
    
    with open(combined_path, "w") as outfile:
        found_any = False
        for fname in marker_files:
            fpath = os.path.join(ref_dir, fname)
            if os.path.exists(fpath):
                found_any = True
                with open(fpath, "r") as infile:
                    outfile.write(infile.read())
                    if not infile.read().endswith("\n"):
                        outfile.write("\n")
            else:
                pass
                # print(f"Warning: Marker file {fname} not found in {ref_dir}")
                
    if not found_any:
        return None
        
    return combined_path

def create_merged_ref(ref_dir, output_path):
    """Merges all .fna, .fa, .fasta files in ref_dir into a single FASTA file."""
    ref_files = []
    for ext in ["*.fna", "*.fa", "*.fasta"]:
        ref_files.extend(glob.glob(os.path.join(ref_dir, ext)))
    
    count = 0
    with open(output_path, "w") as outfile:
        for fname in ref_files:
            basename = os.path.basename(fname)
            
            # Skip marker gene files from the primary reference set
            if "_cpsK" in basename or "_pssD" in basename or "_unique" in basename or "Marker_" in basename:
                continue
                
            if "Serotype_cps_" in basename:
                serotype = os.path.splitext(basename)[0].replace("Serotype_cps_", "")
            else:
                serotype = os.path.splitext(basename)[0]
            
            for record in SeqIO.parse(fname, "fasta"):
                record.id = f"{serotype}|{record.id}"
                record.description = ""
                SeqIO.write(record, outfile, "fasta")
                count += 1
    return count

def run_blast(query_fasta, subject_db, output_txt, task="megablast"):
    """Runs blastn."""
    # Added qseq sseq for SNP analysis
    cmd = [
        "blastn",
        "-query", query_fasta,
        "-db", subject_db,
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qseq sseq",
        "-out", output_txt,
        "-task", task
    ]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

def calculate_coverage(blast_results_file):
    """
    Parses BLAST results and calculates coverage and weighted identity for each query (reference).
    Returns a dict: {ref_id: {'coverage': float, 'identity': float, 'bitscore': float}}
    """
    hits = {}
    q_lengths = {}

    if not os.path.exists(blast_results_file) or os.path.getsize(blast_results_file) == 0:
        return {}

    with open(blast_results_file, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if len(row) < 13: continue # Safety check
            qseqid = row[0]
            pident = float(row[2])
            length = int(row[3])
            qstart = int(row[6])
            qend = int(row[7])
            bitscore = float(row[11])
            qlen = int(row[12])

            q_lengths[qseqid] = qlen
            
            if qseqid not in hits:
                hits[qseqid] = []
            
            # Normalize start/end
            start = min(qstart, qend)
            end = max(qstart, qend)
            hits[qseqid].append((start, end, pident, bitscore, length))

    results = {}
    for ref_id, intervals in hits.items():
        # Calculate coverage (merge overlapping intervals)
        intervals.sort(key=lambda x: x[0])
        merged = []
        if intervals:
            curr_start, curr_end = intervals[0][0], intervals[0][1]
            for next_start, next_end, _, _, _ in intervals[1:]:
                if next_start <= curr_end + 1: # Overlap or adjacent
                    curr_end = max(curr_end, next_end)
                else:
                    merged.append((curr_start, curr_end))
                    curr_start, curr_end = next_start, next_end
            merged.append((curr_start, curr_end))
        
        covered_bases = sum(end - start + 1 for start, end in merged)
        total_length = q_lengths[ref_id]
        coverage = (covered_bases / total_length) * 100

        # Calculate weighted identity
        total_aligned_len = sum(x[4] for x in intervals)
        weighted_id_sum = sum(x[2] * x[4] for x in intervals)
        avg_identity = weighted_id_sum / total_aligned_len if total_aligned_len > 0 else 0
        
        max_bitscore = sum(x[3] for x in intervals)
        
        results[ref_id] = {
            'coverage': coverage,
            'identity': avg_identity,
            'bitscore': max_bitscore
        }

    return results

def check_snps_in_blast(blast_results_file, marker_id, snp_definitions):
    """
    Checks specific SNP positions in the BLAST alignment.
    snp_definitions: list of dicts {'pos': int, 'ref1': char, 'ref14': char}
    Returns: dict {'1': count, '14': count, 'total': count}
    """
    votes = {'1': 0, '14': 0, 'total': 0}
    
    if not os.path.exists(blast_results_file):
        return votes

    with open(blast_results_file, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if len(row) < 15: continue
            
            # Filter low identity hits to avoid paralogs/noise
            pident = float(row[2])
            if pident < 95.0: continue

            qseqid = row[0]
            if qseqid != marker_id: continue
            
            qstart = int(row[6])
            qend = int(row[7])
            qseq = row[13] # Query is Marker (Ref1)
            sseq = row[14] # Subject is Sample
            
            # Iterate through desired SNPs
            for snp in snp_definitions:
                target_pos = snp['pos']
                
                # Check if target_pos is within this HSP
                if qstart <= target_pos <= qend:
                    # Find index in qseq corresponding to target_pos
                    # qseq contains gaps, so we must count non-gap chars
                    
                    current_q_pos = qstart
                    found_idx = -1
                    
                    for i, char in enumerate(qseq):
                        if char != '-':
                            if current_q_pos == target_pos:
                                found_idx = i
                                break
                            current_q_pos += 1
                            
                    if found_idx != -1:
                        # Check subject base at this position
                        # sseq[found_idx] is the base in the sample
                        sample_base = sseq[found_idx].upper()
                        ref1_base = snp['ref1'].upper()
                        ref14_base = snp['ref14'].upper()
                        
                        # Note: If sample_base is '-', it's a deletion in sample.
                        if sample_base == ref1_base:
                            votes['1'] += 1
                        elif sample_base == ref14_base:
                            votes['14'] += 1
                        votes['total'] += 1
                        
    return votes

def check_marker_genes(genome_fasta, marker_db_path, temp_dir):
    """
    Checks marker genes against the genome.
    Returns detailed stats (coverage, identity, bitscore) for each marker.
    """
    if not os.path.exists(marker_db_path):
        return None, None

    # Reuse the genome database created in main
    genome_db = os.path.join(temp_dir, "genome_db")
    marker_out = os.path.join(temp_dir, "marker_results.txt")
    
    try:
        # Query=Markers, Subject=Genome
        # task blastn is better for markers than megablast if they are short/divergent, 
        # but megablast is fine for high identity.
        run_blast(marker_db_path, genome_db, marker_out, task="blastn")
    except Exception as e:
        print(f"Marker blast failed: {e}")
        return None, None
        
    # Use calculate_coverage to get full stats
    stats = calculate_coverage(marker_out)
    return stats, marker_out

def main():
    parser = argparse.ArgumentParser(description="Streptococcus suis CPS Serotyping Tool v3 (SNP-based)")
    parser.add_argument("-f", "--fasta_dir", required=True, help="Directory containing genome FASTA files")
    parser.add_argument("-r", "--ref_dir", required=True, help="Directory containing reference CPS sequences")
    parser.add_argument("--marker_db", default=DEFAULT_MARKER_DB, help=f"Path to marker genes FASTA (default: {DEFAULT_MARKER_DB})")
    parser.add_argument("-o", "--output", default=DEFAULT_OUTPUT_FILE, help=f"Output CSV file path (default: {DEFAULT_OUTPUT_FILE})")
    parser.add_argument("-min_co", "--min_co", dest="min_coverage", type=float, default=DEFAULT_MIN_COVERAGE, help=f"Minimum coverage percentage (default: {DEFAULT_MIN_COVERAGE})")
    parser.add_argument("-min_id", "--min_id", dest="min_identity", type=float, default=DEFAULT_MIN_IDENTITY, help=f"Minimum identity percentage (default: {DEFAULT_MIN_IDENTITY})")
    
    args = parser.parse_args()
    
    QUERY_DIR = args.fasta_dir
    REF_DIR = args.ref_dir
    MARKER_DB = args.marker_db
    OUTPUT_FILE = args.output
    MIN_COVERAGE = args.min_coverage
    MIN_IDENTITY = args.min_identity

    with tempfile.TemporaryDirectory() as temp_dir:
        print(f"Working in temporary directory: {temp_dir}")
        
        print("Step 1: Preparing reference sequences...")
        merged_ref_file = os.path.join(temp_dir, "all_refs.fasta")
        ref_count = create_merged_ref(REF_DIR, merged_ref_file)
        print(f"Merged {ref_count} reference sequences.")

        input_files = []
        for ext in ["*.fasta", "*.fa", "*.fna"]:
            input_files.extend(glob.glob(os.path.join(QUERY_DIR, ext)))
            
        if not input_files:
            print(f"No .fasta/.fa/.fna files found in {QUERY_DIR}")
            return

        print(f"Found {len(input_files)} genome files to process.")

        results_table = []

        for genome_file in input_files:
            sample_name = os.path.splitext(os.path.basename(genome_file))[0]
            print(f"Processing {sample_name}...")

            temp_genome_path = os.path.join(temp_dir, "current_genome.fasta")
            try:
                records = list(SeqIO.parse(genome_file, "fasta"))
                SeqIO.write(records, temp_genome_path, "fasta")
            except Exception as e:
                print(f"Error reading {genome_file}: {e}")
                continue

            # 1. Make BLAST DB for the genome (Subject)
            db_path = os.path.join(temp_dir, "genome_db")
            for ext in ['.ndb', '.nhr', '.nin', '.not', '.nsq', '.ntf', '.nto']:
                 fpath = db_path + ext
                 if os.path.exists(fpath):
                     os.remove(fpath)
                     
            cmd_db = ["makeblastdb", "-in", temp_genome_path, "-dbtype", "nucl", "-out", db_path]
            try:
                subprocess.run(cmd_db, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            except subprocess.CalledProcessError as e:
                print(f"Error running makeblastdb for {sample_name}: {e}")
                continue

            # 2. Run BLAST (Query=Refs, Subject=Genome)
            blast_out = os.path.join(temp_dir, "blast_results.txt")
            try:
                run_blast(merged_ref_file, db_path, blast_out)
            except subprocess.CalledProcessError as e:
                print(f"Error running blastn for {sample_name}: {e}")
                continue

            # 3. Analyze results
            stats = calculate_coverage(blast_out)
            
            effective_min_id = min(MIN_IDENTITY, 90.0)
            effective_min_cov = min(MIN_COVERAGE, 80.0)
            
            valid_refs = {
                k: v for k, v in stats.items() 
                if v['coverage'] >= effective_min_cov and v['identity'] >= effective_min_id
            }

            sorted_refs = sorted(valid_refs.items(), key=lambda x: (x[1]['coverage'], x[1]['identity']), reverse=True)
            all_sorted_refs = sorted(stats.items(), key=lambda x: (x[1]['coverage'], x[1]['identity']), reverse=True)

            best_match = None
            notes = ""
            top_cov = 0
            top_id = 0

            if not sorted_refs:
                best_match = "NA"
                if all_sorted_refs:
                    top_raw = all_sorted_refs[0]
                    raw_name = top_raw[0].split('|')[0]
                    top_cov = top_raw[1]['coverage']
                    top_id = top_raw[1]['identity']
                    notes = f"Best match: {raw_name} (Below threshold, Cov:{top_cov:.1f}%, ID:{top_id:.1f}%)"
                else:
                    notes = "No hits found"
            else:
                top_hit = sorted_refs[0]
                ref_name = top_hit[0].split('|')[0]
                top_cov = top_hit[1]['coverage']
                top_id = top_hit[1]['identity']
                
                is_putative = False
                note_parts = []
                if top_id > 90.0 and top_id < 95.0:
                    is_putative = True
                    note_parts.append("Identity between 90% and 95%")
                if top_cov > 80.0 and top_cov < 90.0:
                    is_putative = True
                    note_parts.append("Coverage between 80% and 90%")

                if is_putative:
                    best_match = f"Putative {ref_name}"
                    notes += ". ".join(note_parts) + ". " if note_parts else ""
                else:
                    best_match = ref_name

                # --- Marker Gene Tie-Breaker ---
                
                check_group = None
                if ref_name in ['2', '0.5']:
                    check_group = '2_0.5'
                elif ref_name in ['1', '14']:
                    check_group = '1_14'
                
                # Also check second hit
                if len(sorted_refs) > 1:
                    second_hit = sorted_refs[1]
                    second_name = second_hit[0].split('|')[0]
                    if not check_group:
                        if (ref_name in ['2', '0.5'] and second_name in ['2', '0.5']):
                            check_group = '2_0.5'
                        elif (ref_name in ['1', '14'] and second_name in ['1', '14']):
                            check_group = '1_14'

                if check_group:
                    print(f"  Ambiguous group {check_group} detected. Running marker gene check...")
                    
                    dynamic_marker_db = prepare_marker_db(REF_DIR, temp_dir)
                    if not dynamic_marker_db:
                        print("    Error: No marker genes found.")
                        marker_hits = None
                        marker_blast_file = None
                    else:
                        marker_hits, marker_blast_file = check_marker_genes(temp_genome_path, dynamic_marker_db, temp_dir)
                    
                    if marker_hits:
                        if check_group == '2_0.5':
                            # SNP-based check using Serotype_cps_2_cpsK
                            # SNPs: Pos 603 (A=2, G=0.5), Pos 714 (G=2, A=0.5)
                            snps_cpsK_2_05 = [
                                {'pos': 603, 'ref1': 'A', 'ref14': 'G'}, # ref1=2, ref14=0.5
                                {'pos': 714, 'ref1': 'G', 'ref14': 'A'}
                            ]
                            
                            # Using check_snps_in_blast. 
                            # '1' count -> ref1 (Serotype 2)
                            # '14' count -> ref14 (Serotype 0.5)
                            votes_2_05 = check_snps_in_blast(marker_blast_file, 'Serotype_cps_2_cpsK', snps_cpsK_2_05)
                            
                            v2 = votes_2_05['1']
                            v05 = votes_2_05['14']
                            total_votes = v2 + v05
                            
                            print(f"    SNP Voting (cpsK 2 vs 0.5): 2={v2}, 0.5={v05}, Total={total_votes}")
                            
                            if total_votes > 0:
                                if v2 > v05:
                                    best_match = "2" if not is_putative else "Putative 2"
                                    notes += f"; Re-assigned to 2 (SNP Vote: 2={v2}, 0.5={v05})"
                                elif v05 > v2:
                                    best_match = "0.5" if not is_putative else "Putative 0.5"
                                    notes += f"; Re-assigned to 0.5 (SNP Vote: 0.5={v05}, 2={v2})"
                                else:
                                    notes += f"; Ambiguous: Equal SNP votes (2={v2}, 0.5={v05})"
                            else:
                                notes += "; Ambiguous: No valid SNP votes found for 2 vs 0.5"
                        
                        elif check_group == '1_14':
                            # SNP-based check using Marker_1_14_cpsK (Serotype 1 cpsK gene)
                            # Key SNPs identified in cpsK homologs (S1 vs S14):
                            # Pos 492: T (S1) vs G (S14) - This is the primary determinant (AA 164).
                            # Pos 653/663 are variable in S14 samples (often S1-like), so we exclude them from voting.
                            snps_cpsK = [
                                {'pos': 492, 'ref1': 'T', 'ref14': 'G'}
                            ]
                            
                            # Note: The marker file ID is 'S1_cpsK' (as we copied S1_cpsK.fna to Marker_1_14_cpsK.fna)
                            # We need to make sure we match the correct ID in blast results.
                            # Usually blast ID is the first word.
                            votes_cpsK = check_snps_in_blast(marker_blast_file, 'S1_cpsK', snps_cpsK)

                            
                            v1 = votes_cpsK['1']
                            v14 = votes_cpsK['14']
                            total_votes = v1 + v14
                            
                            print(f"    SNP Voting (cpsK): 1={v1}, 14={v14}, Total={total_votes}")
                            
                            if total_votes > 0:
                                if v14 > v1:
                                    best_match = "14" if not is_putative else "Putative 14"
                                    notes += f"; Re-assigned to 14 (SNP Vote: 14={v14}, 1={v1})"
                                elif v1 > v14:
                                    best_match = "1" if not is_putative else "Putative 1"
                                    notes += f"; Re-assigned to 1 (SNP Vote: 1={v1}, 14={v14})"
                                else:
                                    notes += f"; Ambiguous: Equal SNP votes (1={v1}, 14={v14})"
                            else:
                                notes += "; Ambiguous: No valid SNP votes found"

            print(f"  => Final Prediction for {sample_name}: {best_match}")

            results_table.append({
                'Sample': sample_name,
                'Predicted_Serotype': best_match,
                'Coverage': f"{top_cov:.2f}",
                'Identity': f"{top_id:.2f}",
                'Notes': notes
            })

    try:
        with open(OUTPUT_FILE, "w", newline="") as csvfile:
            fieldnames = ['Sample', 'Predicted_Serotype', 'Coverage', 'Identity', 'Notes']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            for row in results_table:
                writer.writerow(row)
        print(f"Done. Results saved to {OUTPUT_FILE}")
    except Exception as e:
        print(f"Error writing output file: {e}")

if __name__ == "__main__":
    main()
