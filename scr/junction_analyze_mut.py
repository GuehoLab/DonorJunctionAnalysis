import pysam
import os
import sys
import matplotlib.pyplot as plt
import argparse
from typing import Dict, Tuple, Any

def rc(s: str) -> str:
    """Computes the reverse complement of a DNA sequence."""
    return s.translate(str.maketrans('ATGCatgc','TACGtacg'))[::-1]

def main():
    parser = argparse.ArgumentParser(
        description="Analyze mutations near a cleavage site using UMI clustering from a sorted BAM file."
    )

    # Required positional arguments (corresponding to sys.argv in the original code)
    parser.add_argument("input_file_prefix", help="Prefix for input BAM file (e.g., 'sample_plus' for 'sample_plus.sorted.bam').")
    parser.add_argument("donor_sequence", help="Full donor sequence used for HDR.")
    parser.add_argument("standard_chromosome", help="Chromosome name of the cleavage site (e.g., 'chrX').")
    parser.add_argument("standard_position", type=int, help="1-based position of the cleavage site.")

    # Optional arguments (corresponding to global variables with default values)
    filter_group = parser.add_argument_group("Filtering Options")
    filter_group.add_argument("--umi_size_min", type=int, default=1, 
                              help="Minimum required size of a UMI cluster (default: 0).")
    filter_group.add_argument("--read_count_min", type=int, default=1, 
                              help="Minimum read count for the consensus mutation within a UMI cluster (default: 0).")

    args = parser.parse_args()

    # --- Variables Initialization from Args ---
    fn = args.input_file_prefix
    donor = args.donor_sequence
    standard_chr = args.standard_chromosome
    standard_pos = args.standard_position
    umi_size_min = args.umi_size_min
    read_cnt_min = args.read_count_min
    
    BAM_FILE_PATH = f'{fn}.sorted.bam'
    
    if not os.path.exists(BAM_FILE_PATH):
        sys.exit(f"Error: BAM file not found at {BAM_FILE_PATH}")
        
    # --- Start Core Logic (Original code translated to function scope) ---
    
    with pysam.AlignmentFile(BAM_FILE_PATH, 'rb') as f:
        umi_d: Dict[str, Dict[Tuple[Any, ...], int]] = {}
        
        # Open file handle for low-level mutation IDs
        with open(f'{fn}_low_level_mut_ids.txt', 'w') as fw:
        
            if 'plus' in fn:
                align_strand = True
                donor_end = donor[-5:]
            else:
                align_strand = False
                donor_end = donor[-5:]
            
            # Fetch reads from the BAM file
            for line in f:
                if line.reference_name != standard_chr or line.mapq < 30:
                    continue
                
                # UMI parsing (assuming format: sid_umi_mut_donor)
                try:
                    sid, umi, mut_donor_str = line.query_name.split('_')
                    mut_donor = int(mut_donor_str)
                except ValueError:
                    continue # Skip if query_name format is wrong

                query_seq = line.query_sequence
                end_seq = ''
                ins_len = 0
                
                # --- Complex Alignment Logic (Identical to original) ---
                if align_strand: 
                    if line.is_reverse:
                        del_len = line.reference_end - standard_pos
                        insert_pos = line.reference_end
                        if line.cigar and line.cigar[-1][0] == 4:
                            ins_len = line.cigar[-1][1]
                            end_seq = rc(query_seq[-ins_len - 10:])
                        else:
                            end_seq = rc(query_seq[-10:])
                        inv = True
                    else:
                        del_len = standard_pos - line.reference_start
                        insert_pos = line.reference_start
                        if line.cigar and line.cigar[0][0] == 4:
                            ins_len = line.cigar[0][1]
                            end_seq = query_seq[:ins_len+10]
                        else:
                            end_seq = query_seq[:10]
                        inv = False
                else: 
                    if line.is_reverse:
                        del_len = line.reference_end - standard_pos
                        insert_pos = line.reference_end
                        if line.cigar and line.cigar[-1][0] == 4:
                            ins_len = line.cigar[-1][1]
                            end_seq = rc(query_seq[-ins_len-10:])
                        else:
                            end_seq = rc(query_seq[-10:])
                        inv = False
                    else:
                        del_len = standard_pos - line.reference_start
                        insert_pos = line.reference_start
                        if line.cigar and line.cigar[0][0] == 4:
                            ins_len = line.cigar[0][1]
                            end_seq = query_seq[:ins_len+10]
                        else:
                            end_seq = query_seq[:10]
                        inv = True
                # --- End Complex Alignment Logic ---

                # --- Junction Sequence and Filtering Logic (Identical to original) ---
                if ins_len >= 5:
                    if end_seq and end_seq[-15: -10] == donor_end:
                        ins_len = 0
                        end_seq = end_seq[-10:]
                        mut_donor = 0
                        
                donor_part = donor[:len(donor) + mut_donor]
                end_seq = donor_part[-10:] + '/' + end_seq
                info = (mut_donor, del_len, ins_len, end_seq, inv, insert_pos)
                
                # Filtering low level mutations
                if info[:5] not in [(0, 0, 0, False, 7111535), (0, 56, 0, False, 71111591)]:
                    fw.write(line.query_name + '\n')

                # UMI Clustering
                if umi not in umi_d:
                    umi_d[umi] = {info: 1}
                else:
                    umi_d[umi][info] = umi_d[umi].get(info, 0) + 1
        
    # --- UMI Filtering and Mutation Aggregation ---
    umi_hist: Dict[int, int] = {}
    filtered_mut: Dict[Tuple[Any, ...], int] = {}
    small_mut: Dict[Tuple[Any, ...], int] = {}
    all_cnt = 0
    small_all_cnt = 0
    
    for x, y in umi_d.items():
        n = sum(y.values())
        max_m = sorted(y.items(), key=lambda x: x[1], reverse=True)[0]
        
        # UMI size distribution
        umi_hist[n] = umi_hist.get(n, 0) + 1
        
        # Filtering based on min size and min read count
        if n > umi_size_min and max_m[1] > read_cnt_min:
            m = max_m[0]
            
            # Aggregate all filtered mutations
            filtered_mut[m] = filtered_mut.get(m, 0) + 1
            all_cnt += 1
            
            # Aggregate small mutations (< 1000bp deletion)
            if abs(m[1]) < 1000: 
                small_mut[m] = small_mut.get(m, 0) + 1
                small_all_cnt += 1
                
    # --- Output Writing ---

    # 1. UMI Size Distribution File
    with open(f'{fn}_umi_size.txt', 'w') as fw:
        for x, y in sorted(umi_hist.items(), key= lambda x: x[0]):
            fw.write(f'{x}\t{y}\n')

    # 2. Full Mutation Table
    with open(f'{fn}_mutation_table.txt', 'w') as fw:
        fw.write('Donor_del\tPos_diff\tins_len\tjunction_seq\tinversion\tinsert_position\tcount\t(%)\n')
        for x, y in sorted(filtered_mut.items(), key=lambda x: x[1], reverse=True):
            s = '\t'.join(map(str,x))
            percent = round(y*100/all_cnt,2) if all_cnt > 0 else 0
            fw.write(f'{s}\t{y}\t{percent}\n')

    # 3. Shorter Mutation Table
    with open(f'{fn}_mutation_shorter_1000bp_table_wo_min.txt', 'w') as fw:
        fw.write('Donor_del\tPos_diff\tins_len\tjunction_seq\tinversion\tinsert_position\tcount\t(%)\n')
        for x, y in sorted(small_mut.items(), key=lambda x: x[1], reverse=True):
            s = '\t'.join(map(str,x))
            percent = round(y*100/small_all_cnt,2) if small_all_cnt > 0 else 0
            fw.write(f'{s}\t{y}\t{percent}\n')


if __name__ == '__main__':
    main()