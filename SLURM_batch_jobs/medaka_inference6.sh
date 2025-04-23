#!/bin/bash
#SBATCH -J medaka_inference6                        # Job name
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
medaka inference /work/gatins/hci_genome/processing/medaka/medaka_align.bam.bam /work/gatins/hci_genome/processing/medaka/results_inference/contigs333-369.hdf --model r1041_e82_400bps_sup_v5.0.0 --threads 2 --region contig_333 contig_336 contig_339 contig_341 contig_342 contig_344 contig_345 contig_346 contig_347 contig_348 contig_349 contig_350 contig_351 contig_352 contig_353 contig_354 contig_356 contig_361 contig_366 contig_369
