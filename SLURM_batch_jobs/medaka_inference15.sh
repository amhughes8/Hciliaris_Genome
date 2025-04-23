#!/bin/bash
#SBATCH -J medaka_inference15                        # Job name
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
medaka inference /work/gatins/hci_genome/processing/medaka/medaka_align.bam.bam /work/gatins/hci_genome/processing/medaka/results_inference/contigs1157-1201.hdf --model r1041_e82_400bps_sup_v5.0.0 --threads 2 --region contig_1157 contig_1158 contig_1163 contig_1164 contig_1169 contig_1174 contig_1176 contig_1180 contig_1184 contig_1187 contig_1188 contig_1189 contig_1190 contig_1191 contig_1193 contig_1194 contig_1195 contig_1196 contig_1199 contig_1201
