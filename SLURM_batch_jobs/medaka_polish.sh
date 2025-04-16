#!/bin/bash
#SBATCH -J medaka_polish                        # Job name
#SBATCH -p short                            # Partition
#SBATCH -N 1                                # Number of nodes
#SBATCH -n 32                                # Number of tasks/threads
#SBATCH -o output_%j.txt                    # Standard output file
#SBATCH -e error_%j.txt                     # Standard error file
#SBATCH --mail-user=hughes.annab@northeastern.edu  # Email
#SBATCH --mail-type=END                     # Email notification at job completion
#SBATCH --time=48:00:00                     # Maximum run time


# Your program/command here
module load anaconda3/2022.05
source activate medaka
medaka_consensus -i /work/gatins/hci_genome/processing/hci_concat_noadapters.fastq -d /work/gatins/hci_genome/processing/assembly_Flye/assembly.fasta -o medaka_out -t 32 -m r1041_e82_400bps_sup_v5.0.0
