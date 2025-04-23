#!/bin/bash
#SBATCH -J medaka_inference4                        # Job name
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
medaka inference /work/gatins/hci_genome/processing/medaka/medaka_align.bam.bam /work/gatins/hci_genome/processing/medaka/results_inference/contigs256-287.hdf --model r1041_e82_400bps_sup_v5.0.0 --threads 2 --region contig_256 contig_257 contig_258 contig_259 contig_261 contig_263 contig_264 contig_265 contig_269 contig_272 contig_273 contig_274 contig_275 contig_276 contig_279 contig_280 contig_281 contig_284 contig_285 contig_287
