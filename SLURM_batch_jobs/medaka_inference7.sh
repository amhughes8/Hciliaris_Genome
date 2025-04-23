#!/bin/bash
#SBATCH -J medaka_inference7                        # Job name
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
medaka inference /work/gatins/hci_genome/processing/medaka/medaka_align.bam.bam /work/gatins/hci_genome/processing/medaka/results_inference/contigs370-420.hdf --model r1041_e82_400bps_sup_v5.0.0 --threads 2 --region contig_370 contig_375 contig_376 contig_379 contig_389 contig_392 contig_396 contig_397 contig_399 contig_400 contig_401 contig_407 contig_408 contig_409 contig_412 contig_413 contig_414 contig_415 contig_418 contig_420
