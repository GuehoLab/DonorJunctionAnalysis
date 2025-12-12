import argparse, sys
import os
import subprocess


def main():
	
	parser = argparse.ArgumentParser(
		description='Analysis of mutation at junction site in integration',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter 
	)
	
	required = parser.add_argument_group('Required arguments')
	required.add_argument('fastq_file', type=str, help='FASTQ file with donor integration sequence. fastq file should have UMI sequence at header.')
	required.add_argument('output_dir', type=str, help='Output directory name')
	required.add_argument('output_name', type=str, help='Output file name')
	required.add_argument('insert_sequence', type=str, help="Insert sequence")
	required.add_argument('genome_ref', type=str, help="Bowtie2 indexed file")
	required.add_argument('cleavage_pos', type=str, help="The position of cleavage site. e.g. chr14:22547458. the chromosome must be same with in the genome reference")

	cutadapt_illumina = parser.add_argument_group('Options for cutadapt to remove Illumina adapt')
	cutadapt_illumina.add_argument('--illumina_adaptor', type=str, default='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT,CTGTCTCTTATACACATCT',   help='Illumina adaptor sequence')
	cutadapt_illumina.add_argument('--illumina_cutadapt_option', type=str, default='--overlap 10 --error-rate 0.10 -q 20 -m 20', help='cutadapt option')

	insert_trim_group = parser.add_argument_group('Insert trimming options')
	insert_trim_group.add_argument("--min_overlap_start", type=int, default=10, help="Minimum overlap")
	insert_trim_group.add_argument("--min_trimmed_length", type=int, default=50, help="Minimum trimmed length")
	insert_trim_group.add_argument("--strict_end_length", type=int, default=5, help="Strict matching length at read end")
	insert_trim_group.add_argument("--max_error_rate", type=float, default=0.1, help="Maximum allowed error rate")
	insert_trim_group.add_argument("--max_reads", type=int, default=5_000_000, help="Maximum number of reads")

	bowtie_group = parser.add_argument_group('bowtie2 options')
	bowtie_group.add_argument("--bowtie_option", type=str, default='--very-sensitive-local', help="option for bowtie2")

	mut_analyze_group = parser.add_argument_group('mutation analysis options')
	mut_analyze_group.add_argument("--umi_size_min", type=int, default=1, help="option for bowtie2")
	mut_analyze_group.add_argument("--read_count_min", type=int, default=1, help="option for bowtie2")

	system_group = parser.add_argument_group('System options')
	system_group.add_argument('-t', '--threads', type=int, default=8, help='Number of threads')
	

	def rc(s): return s.translate(s.maketrans('ATGC','TACG'))[::-1]

	args = parser.parse_args()

	OUTPUT_DIR = args.output_dir
	OUTPUT_NAME = args.output_name
	THREADS = str(args.threads)
	INPUT_FILE = args.fastq_file

	try:
		CLEAVAGE_CHR, cleavage_pos_raw = args.cleavage_pos.split(':')
		CLEAVAGE_POS = int(cleavage_pos_raw)
	except (ValueError, AttributeError):
		print("FATAL ERROR: Invalid cleavage_pos format. Expected 'chr:pos'.")
		sys.exit(1)
	
	os.makedirs(OUTPUT_DIR, exist_ok=True)

	# --- Configuration ---
	# Define the actual cutadapt command parameters here.
	# Replace the placeholder values with your specific data/adapter sequences.
	ADAPTER_SEQUENCE = args.illumina_adaptor.upper().split(',')
	CUTADAPT_OUTPUT_FILE = os.path.join(OUTPUT_DIR, f"{OUTPUT_NAME}_trimmed_reads.fastq.gz")
	CUTADAPT_OPT = args.illumina_cutadapt_option.split(' ')

	# --- Cutadapt Command Setup ---
	# Arguments must be defined as a list of strings.
	cutadapt_command = ["cutadapt"]
	for i in CUTADAPT_OPT:
		cutadapt_command.append(i)
	for i in ADAPTER_SEQUENCE:
		cutadapt_command += ['-a', i]
	cutadapt_command.append(f'--cores={THREADS}')	
	cutadapt_command += ['-o', CUTADAPT_OUTPUT_FILE, INPUT_FILE]

	# --- Execution Block ---
	print(f"Starting Cutadapt execution: {' '.join(cutadapt_command)}")

	try:
		# Execute the command. This call blocks until cutadapt finishes.
		# check=True: Raises CalledProcessError if cutadapt exits with a non-zero status.
		# stdout/stderr=sys.stdout/sys.stderr: Directs cutadapt's real-time output/progress 
		#                                     messages to the terminal.
		subprocess.run(
			cutadapt_command,
			check=True,
			stdout=sys.stdout, 
			stderr=sys.stderr 
		)
		
		# Execution reaches here only upon successful completion.
		print(f"\nCutadapt finished successfully. Output saved to: {CUTADAPT_OUTPUT_FILE}")
		
	except subprocess.CalledProcessError as e:
		# Error during command execution (e.g., cutadapt reported an internal error)
		print(f"\nFATAL ERROR: Cutadapt execution failed. Exit code: {e.returncode}")
		sys.exit(1)
		
	except FileNotFoundError:
		# Error if the 'cutadapt' executable itself cannot be found in the system PATH
		print("\nFATAL ERROR: 'cutadapt' command not found. Check installation and PATH.")
		sys.exit(1)


	# --- insert trimming ---
	TRIM_OUTPUT_FILE = os.path.join(OUTPUT_DIR, f"{OUTPUT_NAME}_insert_trimmed_reads.fastq.gz")
	INSERT_SEQ = args.insert_sequence
	insert_command = ['junction_trim_insert', 
		'-i', CUTADAPT_OUTPUT_FILE, 
		'-o', TRIM_OUTPUT_FILE,
		'-s', INSERT_SEQ,
		'--min_overlap_start', str(args.min_overlap_start),
		'--min_trimmed_length', str(args.min_trimmed_length),
		'--strict_end_length', str(args.strict_end_length),
		'--max_error_rate', str(args.max_error_rate),
		'--max_reads', str(args.max_reads),
		'-t', THREADS]

	print(f"Starting trim_insert.py execution: {' '.join(insert_command)}")

	try:
		# Execute the command. This call blocks until the script finishes.
		# check=True: Raises CalledProcessError if the script exits with a non-zero status.
		# The script's output/progress is directed to the terminal in real-time.
		subprocess.run(
			insert_command,
			check=True,
			stdout=sys.stdout, 
			stderr=sys.stderr 
		)
		
		# Execution reaches here only upon successful completion.
		print(f"\ntrim_insert.py finished successfully. Output saved to: {TRIM_OUTPUT_FILE}")
		
	except subprocess.CalledProcessError as e:
		# Error during script execution (e.g., internal error in trim_insert.py)
		print(f"\nFATAL ERROR: trim_insert.py failed. Exit code: {e.returncode}")
		sys.exit(1)
		
	except FileNotFoundError:
		# Error if the 'python' executable or 'trim_insert.py' script file cannot be found
		print("\nFATAL ERROR: 'python' executable or 'trim_insert.py' not found. Check PATH and script location.")
		sys.exit(1)

	# --- bowtie2 alignment ---

	BOWTIE_OUTPUT_FILE = os.path.join(OUTPUT_DIR, f"{OUTPUT_NAME}.sam")
	SORTED_BAM_FILE = os.path.join(OUTPUT_DIR, f"{OUTPUT_NAME}.sorted.bam")

	print("--- 1. Bowtie2 Alignment Start ---")


	bowtie2_command = ['bowtie2',
					'-x', args.genome_ref,
					'-U', TRIM_OUTPUT_FILE,
					'-p', THREADS,
					'-S', BOWTIE_OUTPUT_FILE]
	for i in args.bowtie_option.split(' '):
		bowtie2_command.append(i)

	# Execute Bowtie2
	try:
		print(f"Executing: {' '.join(bowtie2_command)}")
		subprocess.run(
			bowtie2_command,
			check=True,
			stdout=sys.stdout, 
			stderr=sys.stderr 
		)
		print(f"Bowtie2 finished successfully. Output: {BOWTIE_OUTPUT_FILE}")
		
	except subprocess.CalledProcessError as e:
		print(f"\nFATAL ERROR: Bowtie2 alignment failed. Exit code: {e.returncode}")
		sys.exit(1)
	except FileNotFoundError:
		print("\nFATAL ERROR: 'bowtie2' command not found. Check installation and PATH.")
		sys.exit(1)

	# --- 2. SAMtools Sort (SAM -> BAM Conversion and Sorting) ---
	print("\n--- 2. SAMtools Sort Start (SAM to BAM) ---")

	# samtools sort input.sam -@ THREADS -o output.sorted.bam 
	samtools_sort_command = [
		'samtools',
		'sort',
		BOWTIE_OUTPUT_FILE,
		'-@', THREADS,                # Number of threads
		'-o', SORTED_BAM_FILE         # Sorted BAM output file
	]

	# Execute SAMtools Sort
	try:
		print(f"Executing: {' '.join(samtools_sort_command)}")
		# SAMtools uses stderr for progress output, so we direct both to the terminal
		subprocess.run(
			samtools_sort_command,
			check=True,
			stdout=sys.stdout, 
			stderr=sys.stderr 
		)
		print(f"SAMtools sorting finished successfully. Output: {SORTED_BAM_FILE}")
		
		# Optional: Clean up the intermediate SAM file
		# os.remove(BOWTIE_OUTPUT_FILE)
		# print(f"Cleaned up intermediate SAM file: {BOWTIE_OUTPUT_FILE}")
		
	except subprocess.CalledProcessError as e:
		print(f"\nFATAL ERROR: SAMtools sort failed. Exit code: {e.returncode}")
		sys.exit(1)
	except FileNotFoundError:
		print("\nFATAL ERROR: 'samtools' command not found. Check installation and PATH.")
		sys.exit(1)

	# --- 3. SAMtools Index ---
	print("\n--- 3. SAMtools Indexing Start ---")

	# samtools index input.sorted.bam
	samtools_index_command = [
		'samtools',
		'index',
		SORTED_BAM_FILE               # Sorted BAM file to index
	]

	# Execute SAMtools Index
	try:
		print(f"Executing: {' '.join(samtools_index_command)}")
		subprocess.run(
			samtools_index_command,
			check=True,
			stdout=sys.stdout, 
			stderr=sys.stderr 
		)
		print(f"SAMtools indexing finished successfully. Index file: {SORTED_BAM_FILE}.bai")

	except subprocess.CalledProcessError as e:
		print(f"\nFATAL ERROR: SAMtools index failed. Exit code: {e.returncode}")
		sys.exit(1)

	# --- Pipeline Completion ---
	print("\n--- Pipeline Completed Successfully ---")
	print(f"Final aligned and indexed file: {SORTED_BAM_FILE}")

	#--- Mutation analysis ---

	mut_analyze_command = ['junction_analyze_mut',
		SORTED_BAM_FILE.replace('.sorted.bam',''), 
		INSERT_SEQ, 
		CLEAVAGE_CHR,
		str(CLEAVAGE_POS),
		'--umi_size_min', str(args.umi_size_min),
		'--read_count_min', str(args.read_count_min)]
	
	# --- Execution Block ---
	print(f"--- Mutation Analysis Start ---")
	print(f"Executing: {' '.join(mut_analyze_command)}")

	try:
		# Executes the mutation analysis script and waits for completion
		subprocess.run(
			mut_analyze_command,
			check=True,
			stdout=sys.stdout,
			stderr=sys.stderr
		)
		
		print("\nMutation analysis finished successfully.")

	except subprocess.CalledProcessError as e:
		print(f"\nFATAL ERROR: analyze_mut.py failed. Exit code: {e.returncode}")
		sys.exit(1)
	except FileNotFoundError:
		print("\nFATAL ERROR: Python or 'analyze_mut.py' not found.")
		sys.exit(1)

	print("\n--- All Pipeline Steps Completed ---")
	return


if __name__=='__main__':
	main()
