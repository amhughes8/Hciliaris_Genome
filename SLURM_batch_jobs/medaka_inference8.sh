#!/bin/bash
#SBATCH -J medaka_inference8                        # Job name
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
medaka inference /work/gatins/hci_genome/processing/medaka/medaka_align.bam.bam /work/gatins/hci_genome/processing/medaka/results_inference/contigs422-501.hdf --model r1041_e82_400bps_sup_v5.0.0 --threads 2 --region contig_422 contig_423 contig_425 contig_426 contig_429 contig_430 contig_431 contig_432 contig_433 contig_434 contig_435 contig_436 contig_437 contig_439 contig_440 contig_441 contig_446 contig_480 contig_492 contig_501
