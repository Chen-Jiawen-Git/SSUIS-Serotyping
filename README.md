# Streptococcus suis CPS Serotyping Tool

This is a Python-based serotyping tool for *Streptococcus suis*. It utilizes BLAST+ to align Capsular Polysaccharide (CPS) synthesis gene clusters in Whole Genome Sequencing (WGS) data. It incorporates a high-precision differentiation strategy based on **cpsK gene point mutations (SNPs)** for highly similar serotypes (such as Serotype 1 vs 14 and Serotype 2 vs 1/2).

## Key Features

*   **Comprehensive Serotype Coverage**: Supports prediction of standard *Streptococcus suis* serotypes.
*   **High-Precision Differentiation**:
    *   **Serotype 1 vs 14**: Differentiates using key Single Nucleotide Polymorphism (SNP) sites in the `cpsK` gene.
    *   **Serotype 2 vs 1/2**: Differentiates based on minor sequence differences (only 2 base differences) in the `cpsK` gene.
*   **Detailed Results**: Outputs include predicted serotype, coverage, identity, and detailed typing notes.

## Requirements

Before using this tool, please ensure the following software is installed in your environment:

1.  **Python 3.x**
    *   Requires the `biopython` library:
        ```bash
        pip install biopython
        ```
2.  **NCBI BLAST+**
    *   `blastn` and `makeblastdb` must be installed and added to the system PATH environment variable.

## Directory Structure

```
CPS_Serotyping_Tool/
├── cps_serotyping.py      # Main script
├── references/            # Reference sequence database (contains CPS references and marker genes)
├── utils/                 # Utility scripts (for development and debugging)
└── README.md              # Documentation
```

## Quick Start

Place your genome assembly files (`.fasta`, `.fna`, `.fa`) in a folder (e.g., `input_genomes`), then run the following command:

```bash
python cps_serotyping.py -f input_genomes -r references -o results.csv
```

## Parameters

| Parameter | Abbr. | Description | Default |
| :--- | :--- | :--- | :--- |
| `--fasta_dir` | `-f` | **(Required)** Path to the folder containing genome FASTA files | None |
| `--ref_dir` | `-r` | **(Required)** Path to the folder containing reference CPS sequences (the `references` folder in this package) | None |
| `--output` | `-o` | Path for the output CSV file | `serotyping_results_v1.csv` |
| `--marker_db` | None | Marker gene database filename (usually does not need modification) | `marker_genes.fasta` |
| `--min_co` | `-min_co` | Minimum coverage percentage to be considered positive | 95.0 |
| `--min_id` | `-min_id` | Minimum identity percentage to be considered positive | 95.0 |

**Example**:
```bash
# Run with custom thresholds
python cps_serotyping.py -f ./genomes -r ./references -o my_results.csv -min_co 90.0 -min_id 90.0
```

## Typing Principles & Challenges

This tool employs a two-step method for serotyping:

### 1. Initial Screening (CPS Gene Cluster Alignment)
First, sample genomes are aligned against CPS gene cluster reference sequences of all known serotypes using BLASTN. The best match that satisfies both coverage and identity thresholds (default 95%) is selected.

### 2. Differentiation of Difficult Serotypes (Based on cpsK Gene)

Some serotypes have highly homologous CPS gene clusters that are difficult to distinguish based on full-length alignment alone. This tool introduces specific logic for the following two groups:

#### A. Serotype 1 vs Serotype 14
The CPS gene clusters of these two are almost identical. Studies indicate that point mutations in the `cpsK` gene (encoding sialyltransferase) are key to determining the serotype.
*   **Differentiation Strategy**: Detect the base at position **492** of the `cpsK` gene (relative to the start codon).
    *   **Serotype 1**: This position is **T** (encoding Cysteine, Cys).
    *   **Serotype 14**: This position is **G** (encoding Tryptophan, Trp).
*   **Algorithm**: The tool extracts sequences homologous to `cpsK` from the sample and checks the base type at position 492 to make a "voting" decision.

#### B. Serotype 2 vs Serotype 1/2 (Chz)
The CPS gene clusters of these two are also extremely similar, and the `cpsK` gene (972 bp in length) differs by only **2 bases** (99.79% identity).
*   **Differentiation Strategy**: Detect key SNP sites in the `cpsK` gene.
*   **Detection Sites**:
    *   **Pos 603**: Serotype 2 is **A**, Serotype 1/2 is **G**.
    *   **Pos 714**: Serotype 2 is **G**, Serotype 1/2 is **A**.
*   **Algorithm**: Extract the `cpsK` sequence from the sample, detect the bases at the above two sites, and determine the classification through a voting mechanism. If the voting results are tied, it is marked as Ambiguous.

## Output Interpretation (Output CSV)

The output file contains the following columns:

*   **Sample**: Sample filename.
*   **Predicted_Serotype**: Predicted serotype result.
    *   e.g., `2`, `14`, `Putative 1` (Suspected).
    *   If classification is not possible, it shows `NA`.
*   **Coverage**: Coverage (%) of the best matching reference sequence.
*   **Identity**: Average identity (%) of the best matching reference sequence.
*   **Notes**: Detailed notes.
    *   Will indicate if SNP correction was applied (e.g., `Re-assigned to 14 (SNP Vote: 14=1, 1=0)`).
    *   For suspected results (Putative), specific metrics will be explained (e.g., `Identity between 90% and 95%`).

---
*Created for Streptococcus suis Analysis Pipeline.*
