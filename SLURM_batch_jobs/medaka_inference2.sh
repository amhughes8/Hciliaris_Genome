#!/bin/bash
#SBATCH -J medaka_inference2                        # Job name
#SBATCH -p short                            # Partition
#SBATCH -N 1                                # Number of nodes
#SBATCH -n 2                                # Number of tasks/threads
#SBATCH -o output_%j.txt                    # Standard output file
#SBATCH -e error_%j.txt                     # Standard error file
#SBATCH --mail-user=hughes.annab@northeastern.edu  # Email
#SBATCH --mail-type=END                     # Email notification at job completion
#SBATCH --time=48:00:00                     # Maximum run time


# Your program/command here
module load anaconda3/2022.05
source activate medaka
medaka inference /work/gatins/hci_genome/processing/medaka/medaka_align.bam.bam /work/gatins/hci_genome/processing/medaka/results_inference/contigs85-218.hdf --model r1041_e82_400bps_sup_v5.0.0 --threads 2 --region contig_85 contig_107 contig_125 contig_131 contig_141 contig_145 contig_148 contig_153 contig_158 contig_169 contig_170 contig_171 contig_174 contig_181 contig_197 contig_209 contig_214 contig_216 contig_217 contig_218
