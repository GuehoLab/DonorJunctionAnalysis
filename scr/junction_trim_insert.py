import itertools
import sys
import os
import gzip
import argparse
from multiprocessing import Pool, cpu_count
import edlib

def fastq_parser_single_thread(fastq_file_path, max_reads):
    """
    Parses a gzipped FASTQ file, yielding reads up to max_reads limit.
    Accepts max_reads as an argument.
    """
    read_count = 0

    try:
        # Use 'rt' for reading text from a gzip file
        with gzip.open(fastq_file_path, 'rt') as f: 
            # Read 4 lines at a time
            for header, sequence, quality_header, quality_string in itertools.zip_longest(*([f]*4), fillvalue=None):
                if header is None:
                    break
                
                read_count += 1
                if max_reads > 0 and read_count > max_reads:
                    print(f"⚠️ Maximum number of reads ({max_reads:,}) reached. Additional reads will be skipped.")
                    return
                
                yield (
                    header.strip(), 
                    sequence.strip(), 
                    quality_header.strip(), 
                    quality_string.strip()
                )
                
    except FileNotFoundError:
        raise FileNotFoundError(f"Can not find file: {fastq_file_path}")
    except Exception as e:
        raise Exception(f"Error in parsing file: {e}")


def find_best_match_length_edlib(read_seq, adapter_seq, adapter_len, min_overlap_start, strict_end_length, max_error_rate):
    """
    Finds the best match length for trimming, accepting all necessary parameters.
    """
    read_len = len(read_seq)
    max_len = min(read_len, adapter_len)
    
    for current_len in range(max_len, min_overlap_start - 1, -1):
        
        target_read = read_seq[:current_len]
        target_adapter = adapter_seq[:current_len]
        
        # Strict End Check
        if current_len >= strict_end_length:
            insert_end = target_adapter[-strict_end_length:]
            read_end = target_read[-strict_end_length:]
            
            if read_end != insert_end:
                continue

        # Edlib Alignment Check
        result = edlib.align(target_adapter.encode('ascii'), target_read.encode('ascii'), 
                             mode="HW", task="distance")
        
        edit_distance = result["editDistance"]
        error_rate = edit_distance / current_len
        
        if error_rate <= max_error_rate:
            return current_len

    return 0 


def process_single_read(read_record, config):
    """
    Process a single read, taking configuration via a dictionary.
    """
    header, sequence, quality_header, quality_string = read_record

    # Call the matching function with parameters from config
    trim_len = find_best_match_length_edlib(
        sequence, 
        config['ADAPTER_SEQ'],
        config['ADAPTER_LEN'],
        config['MIN_OVERLAP_START'],
        config['STRICT_END_LENGTH'],
        config['MAX_ERROR_RATE']
    )
    
    # Filtering by overlap
    if trim_len < config['MIN_OVERLAP_START']:
        return None
    
    difference_in_bp = trim_len - config['ADAPTER_LEN']
    
    new_sequence = sequence[trim_len:]
    new_quality = quality_string[trim_len:]

    # Filtering by trimmed length
    if len(new_sequence) < config['MIN_TRIMMED_LENGTH']:
        return None
        
    # Header modification (Extract UMI and add difference_in_bp)
    header_parts = header.split()
    try:
        umi = header_parts[1].split('+')[1][:9]
    except IndexError:
        # Handle case where UMI format is unexpected
        umi = "NOUMI9999" 
        
    new_header = f"{header_parts[0]}_{umi}_{difference_in_bp}"
    
    return (
        f"{new_header}\n"
        f"{new_sequence}\n"
        f"{quality_header}\n"
        f"{new_quality}\n"
    )


def process_fastq_read_level_multiprocessing(input_fastq, output_fastq, num_processes, max_reads, config):
    """
    Orchestrates the multiprocessing pipeline, passing all configuration.
    """
    if num_processes is None or num_processes <= 0:
        # Fallback to a default or cpu_count
        num_processes = cpu_count()
    
    print(f"⚙️ Starting multiprocessing. Number of processes: {num_processes}")
    
    # 1. Read Generation (Single-threaded Gzip parsing)
    print(f"   Parsing FASTQ file ({input_fastq}) into individual reads...")
    try:
        read_generator = fastq_parser_single_thread(input_fastq, max_reads)
    except Exception as e:
        print(f"❌ Error occurred while reading the file: {e}", file=sys.stderr)
        return

    # 2. Multiprocessing Pool Setup
    with Pool(num_processes) as pool:
        
        from functools import partial
        process_func = partial(process_single_read, config=config)
        results_iterator = pool.imap(process_func, read_generator)

        print("   Processing reads and merging results...")
        
        # 3. Output Writing
        try:
            with gzip.open(output_fastq, 'wt') as out_f:
                for record in results_iterator:
                    if record is not None:
                        out_f.write(record)
                        
        except Exception as e:
            print(f"❌ Error occurred while writing the output file: {e}", file=sys.stderr)
            pool.terminate()
            return

    print(f"✅ Completed processing of the gzip FASTQ file (read-level multiprocessing).")
    print(f"   Output file: {output_fastq}")


def parse_args():
    parser = argparse.ArgumentParser(description="FASTQ trimming and insert-sequence filtering pipeline")

    io_group = parser.add_argument_group("Input/Output options")
    io_group.add_argument("-i", "--input_fastq", required=True, help="Input FASTQ (gzip) file")
    io_group.add_argument("-o", "--output_fastq", required=True, help="Output FASTQ file")

    insert_group = parser.add_argument_group("Insert trimming options")
    insert_group.add_argument("-s", "--insert_sequence", required=True, help="Insert sequence")
    insert_group.add_argument("--min_overlap_start", type=int, default=10, help="Minimum overlap")
    insert_group.add_argument("--min_trimmed_length", type=int, default=50, help="Minimum trimmed length")
    insert_group.add_argument("--strict_end_length", type=int, default=5, help="Strict matching length at read end")
    insert_group.add_argument("--max_error_rate", type=float, default=0.1, help="Maximum allowed error rate")

    perf_group = parser.add_argument_group("Performance options")
    perf_group.add_argument("-t", "--threads", type=int, default=8, help="Number of threads")

    limit_group = parser.add_argument_group("Processing limits")
    limit_group.add_argument("--max_reads", type=int, default=1_000_000, help="Maximum number of reads (0 for unlimited)")

    return parser.parse_args()


def main():
    args = parse_args()

    # --- Configuration Dictionary ---
    ADAPTER_SEQ = args.insert_sequence.upper()
    config = {
        'ADAPTER_SEQ': ADAPTER_SEQ,
        'ADAPTER_LEN': len(ADAPTER_SEQ),
        'MIN_OVERLAP_START': args.min_overlap_start,
        'MIN_TRIMMED_LENGTH': args.min_trimmed_length,
        'STRICT_END_LENGTH': args.strict_end_length,
        'MAX_ERROR_RATE': args.max_error_rate,
    }

    # --- Execution ---
    print("\n--- Starting Trimming Pipeline ---")
    process_fastq_read_level_multiprocessing(
        input_fastq=args.input_fastq,
        output_fastq=args.output_fastq,
        num_processes=args.threads,
        max_reads=args.max_reads,
        config=config 
    )

if __name__ == '__main__':
    main()