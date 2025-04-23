#!/bin/bash
#SBATCH -J medaka_inference5                        # Job name
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
medaka inference /work/gatins/hci_genome/processing/medaka/medaka_align.bam.bam /work/gatins/hci_genome/processing/medaka/results_inference/contigs289-332.hdf --model r1041_e82_400bps_sup_v5.0.0 --threads 2 --region contig_289 contig_290 contig_291 contig_294 contig_299 contig_300 contig_303 contig_304 contig_305 contig_307 contig_309 contig_313 contig_315 contig_320 contig_322 contig_323 contig_329 contig_330 contig_331 contig_332
