#!/bin/bash
#SBATCH -J hifiasm_assembly_nomt                       # Job name
#SBATCH -p short                            # Partition
#SBATCH -N 1                                # Number of nodes
#SBATCH -n 32                               # Number of tasks/threads
#SBATCH --mem=500G                           # Memory
#SBATCH -o output_%j.txt                    # Standard output file
#SBATCH -e error_%j.txt                     # Standard error file
#SBATCH --mail-user=hughes.annab@northeastern.edu  # Email
#SBATCH --mail-type=END                     # Email notification at job completion
#SBATCH --time=48:00:00                     # Maximum run time


# Your program/command here
module load anaconda3/2022.05
module load samtools/1.19.2
source activate /work/gatins/hci_genome/env
samtools view -b -f 4 -@ 20 mito_aln.sorted.bam > unmapped.bam
samtools fastq unmapped.bam > reads_no_mito.fastq
hifiasm -o assembly_hifiasm_no_mito.asm --ont -t32 /work/gatins/hci_genome/processing/mtdna/removal/reads_no_mito.fastq
