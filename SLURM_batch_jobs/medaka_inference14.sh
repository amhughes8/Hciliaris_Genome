#!/bin/bash
#SBATCH -J medaka_inference14                        # Job name
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
medaka inference /work/gatins/hci_genome/processing/medaka/medaka_align.bam.bam /work/gatins/hci_genome/processing/medaka/results_inference/contigs1110-1153.hdf --model r1041_e82_400bps_sup_v5.0.0 --threads 2 --region contig_1110 contig_1111 contig_1115 contig_1116 contig_1117 contig_1118 contig_1119 contig_1120 contig_1123 contig_1124 contig_1126 contig_1127 contig_1129 contig_1130 contig_1146 contig_1147 contig_1149 contig_1151 contig_1152 contig_1153
