#!/bin/bash
#SBATCH -J medaka_inference3                        # Job name
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
medaka inference /work/gatins/hci_genome/processing/medaka/medaka_align.bam.bam /work/gatins/hci_genome/processing/medaka/results_inference/contigs220-250.hdf --model r1041_e82_400bps_sup_v5.0.0 --threads 2 --region contig_220 contig_221 contig_222 contig_223 contig_225 contig_230 contig_231 contig_232 contig_233 contig_234 contig_235 contig_236 contig_238 contig_241 contig_242 contig_244 contig_246 contig_247 contig_248 contig_250
