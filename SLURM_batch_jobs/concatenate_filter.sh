#!/bin/bash
#SBATCH -J concatenate_filter                        # Job name
#SBATCH -p short                            # Partition
#SBATCH -N 1                                # Number of nodes
#SBATCH -n 10                                # Number of tasks/threads
#SBATCH -o output_%j.txt                    # Standard output file
#SBATCH -e error_%j.txt                     # Standard error file
#SBATCH --mail-user=hughes.annab@northeastern.edu  # Email
#SBATCH --mail-type=END                     # Email notification at job completion
#SBATCH --time=48:00:00                     # Maximum run time


# Your program/command here
module load anaconda3/2022.05 discovery
source activate /work/gatins/hci_genome/env
cat /work/gatins/hci_genome/processing/hci1_noadapters.fastq /work/gatins/hci_genome/processing/hci2_noadapters.fastq > /work/gatins/hci_genome/processing/hci_concat_noadapters.fastq
cat /work/gatins/hci_genome/processing/hci_concat_noadapters.fastq | seqkit seq -m 2000 -Q 3 -j 10 > /work/gatins/hci_genome/processing/hci_filtered_2kQ3.fastq
