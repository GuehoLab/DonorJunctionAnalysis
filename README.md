# DonorJunctionAnalysis


## 🧬 Donor Junction Analysis

**DonorJunctionAnalysis** is an analysis tool designed to determine the precision of integration/insertion events when performing Illumina sequencing (Donor-seq) following genomic modification techniques like Prime Assembly or general DNA integration. It focuses on analyzing the junction sequences where the donor sequence was inserted.

---

### 📥 Installation

The following core requirements are necessary to run the pipeline:

#### Requirements

| Requirement | Version |
| :--- | :--- |
| **Python** | 3.10 |
| **bowtie2** | 2.5.4 |
| **samtools** | 1.21 |
| **cutadapt** | 5.1 |
| **pysam** | 0.23.3 |
| **matplotlib** | 3.10.7 |
| **edlib** | 1.3.9 |

#### Setup Commands

Follow these commands to clone the repository, create the necessary Conda environment (using `mamba` for speed), and install the package in editable mode:

```bash
# 1. Clone the repository
git clone [https://github.com/GuehoLab/DonorJunctionAnalysis](https://github.com/GuehoLab/DonorJunctionAnalysis)

# 2. Navigate to the directory
cd DonorJunctionAnalysis

# 3. Create the Conda environment using the provided environment file
mamba env create -f environment.yml -n mutation_analysis_env

# 4. Activate the newly created environment
mamba activate mutation_analysis_env

# 5. Install the package in editable mode (-e)
pip install -e .
```

~5–10 minutes on a standard desktop computer

### 📝 Usage

The tool is executed using a single command that requires several positional arguments.

#### Command Syntax

```bash
donor_junction_analysis {input_fastq} {output_directory} {output_name} {insert_sequence} {bowtie2_indexed_reference_genome} {cleavage_site}
```

#### Positional Arguments

| Argument | Description | Important Note |
| :--- | :--- | :--- |
| **`input_fastq`** | FASTQ file containing the donor integration sequence. Supports `.fastq` and `.fastq.gz` formats. | The read header **MUST** contain the UMI sequence, typically formatted like `+ATTGGCAG` (9 bases) at the end of the read ID field. |
| **`output_directory`** | Name of the directory where all output files will be saved. | The directory will be created if it does not exist. |
| **`output_name`** | Base name for all output files generated within the `output_directory`. | |
| **`insert_sequence`** | The full donor sequence expected to be inserted into the genome. | |
| **`bowtie2_indexed_reference_genome`** | The genomic reference indexed using `bowtie2`. | |
| **`cleavage_site`** | The target cleavage site position. | Format: `{chromosome}:{position}` (e.g., `chr14:22547458`). The chromosome name must match the name in the Bowtie2 index. |

#### Example

You can run the provided example script within the cloned directory (ensure your reference files are correctly configured, especially if using a specific genome like GRCh38).

```bash
chmod +x ./example_run.sh
./example_run.sh
```

Expected runtime is ~2–5 minutes on a standard desktop computer

### ⚙️ Parameters (Full Options List)

| Group | Argument | Type | Default | Description |
| :--- | :--- | :--- | :--- | :--- |
| **Illumina Adaptor** | `--illumina_adaptor` | `str` | `AGATCGG...CT` | Illumina adaptor sequence(s) for `cutadapt` removal. |
| | `--illumina_cutadapt_option` | `str` | `--overlap 10 --error-rate 0.10 -q 20 -m 20` | Options passed directly to `cutadapt` for Illumina adaptor trimming. |
| **Insert Trimming** | `--min_overlap_start` | `int` | `10` | Minimum required overlap length for successful donor sequence trimming. |
| | `--min_trimmed_length` | `int` | `50` | Minimum required length of the read after trimming. |
| | `--strict_end_length` | `int` | `5` | Length of the read end required to have a strict (perfect) match with the donor sequence. |
| | `--max_error_rate` | `float` | `0.1` | Maximum allowed error rate (Edit Distance / Overlap Length) during alignment/trimming. |
| | `--max_reads` | `int` | `5_000_000` | Maximum number of reads to process (for fast testing/debugging). |
| **Bowtie2** | `--bowtie_option` | `str` | `--very-sensitive-local` | Options passed directly to the `bowtie2` aligner. |
| **Mutation Analysis** | `--umi_size_min` | `int` | `1` | Minimum number of reads in a UMI cluster to be included in the analysis. |
| | `--read_count_min` | `int` | `1` | Minimum read count for the *most frequent mutation* within a UMI cluster to be considered the consensus. |
| **System** | `-t, --threads` | `int` | `8` | Number of threads/processes to use for parallel processing. |

---

### 📊 Result File (`{output_name}_mutation_shorter_1000bp_table_wo_min.txt`)

The primary output file, showing UMI-deduplicated mutation patterns where the deletion/insertion size is less than 1000 bp.

| Column | Description | Example/Context |
| :--- | :--- | :--- |
| **Donor\_del** | **Deletion length from the donor sequence.** | `-20` means the last 20 bp of the donor sequence were lost during insertion. |
| **Pos\_diff** | **Difference in alignment position** relative to the expected `cleavage_site`. | `0` indicates a precise insertion at the expected cleavage site. |
| **ins\_len** | **Length of sequence inserted** between the donor sequence and the genome alignment. | This represents non-templated DNA (N-nucleotides) or short microhomology insertions. |
| **junction\_seq** | **Sequence surrounding the junction.** | Format: `[Donor Sequence End]/[Inserted/Genome Sequence Start]`. |
| **inversion** | Indicates whether the aligned read suggests an **inversion** event. | (`True`/`False`) |
| **insert\_position** | The **genomic position** where the sequence was aligned/inserted. | |
| **count** | **UMI-deduplicated count** for this specific mutation pattern. | The number of unique DNA molecules detected with this outcome. |
| **(%)** | **Percentage** of this mutation pattern among all filtered UMI clusters. | |


