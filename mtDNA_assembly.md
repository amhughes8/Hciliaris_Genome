# Assembling the Mitochondrial Genome from Whole Genome Sequencing Data
These instructions are from Giacomo Bernardi

## 1. Filter anything between 5kb-20kb
```
awk 'BEGIN {OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) < 20000) {print header, seq, qheader, qseq}}' /work/gatins/hci_genome/processing/hci1.fastq > HCI_1_dorado0.9.1_lessthan20000bp.fastq
awk 'BEGIN {OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) < 20000) {print header, seq, qheader, qseq}}' /work/gatins/hci_genome/processing/hci2.fastq > HCI_2_dorado0.9.1_lessthan20000bp.fastq
```
```
awk 'BEGIN {OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) > 5000) {print header, seq, qheader, qseq}}' HCI_1_dorado0.9.1_lessthan20000bp.fastq > HCI_1_dorado0.9.1_lessthan20000bp_morethan5000bp.fastq
awk 'BEGIN {OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) > 5000) {print header, seq, qheader, qseq}}' HCI_2_dorado0.9.1_lessthan20000bp.fastq > HCI_2_dorado0.9.1_lessthan20000bp_morethan5000bp.fastq
```
```
cat HCI_1_dorado0.9.1_lessthan20000bp_morethan5000bp.fastq HCI_2_dorado0.9.1_lessthan20000bp_morethan5000bp.fastq > HCI_dorado0.9.1_lessthan20000bp_morethan5000bp.fastq
```

## 2. Assemble
```
#!/bin/bash
#SBATCH -J assembly_flye_mtdna                        # Job name
#SBATCH -p short                            # Partition
#SBATCH -N 1                                # Number of nodes
#SBATCH -n 8                                # Number of tasks/threads
#SBATCH --mem=100G                          # Memory
#SBATCH -o output_%j.txt                    # Standard output file
#SBATCH -e error_%j.txt                     # Standard error file
#SBATCH --mail-user=hughes.annab@northeastern.edu  # Email
#SBATCH --mail-type=BEGIN,END                     # Email notification at job start and completion
#SBATCH --time=48:00:00                     # Maximum run time


# Your program/command here
module load anaconda3/2022.05 discovery
source activate /work/gatins/hci_genome/env
flye --nano-hq HCI_dorado0.9.1_lessthan20000bp_morethan5000bp.fastq --out-dir /work/gatins/hci_genome/processing/mtdna --genome-size 16.8k --threads 8
```

## 3. Look at assembly stats for a circular contig with high coverage
size should be betwen 16K - 32K
