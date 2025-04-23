#!/bin/bash
#SBATCH -J medaka_align                        # Job name
#SBATCH -p short                            # Partition
#SBATCH -N 1                                # Number of nodes
#SBATCH -n 32                                # Number of tasks/threads
#SBATCH -o output_%j.txt                    # Standard output file
#SBATCH -e error_%j.txt                     # Standard error file
#SBATCH --mail-user=hughes.annab@northeastern.edu  # Email
#SBATCH --mail-type=END                     # Email notification at job completion
#SBATCH --time=48:00:00                     # Maximum run time

module load anaconda3/2022.05
source activate medaka
mini_align -i /work/gatins/hci_genome/processing/hci_concat_noadapters.fastq -r /work/gatins/hci_genome/processing/assembly_Flye/assembly.fasta -P -m -p medaka_align.bam -t 32
