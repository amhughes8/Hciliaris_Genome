#!/bin/bash
#SBATCH -J medaka_inference9                        # Job name
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
medaka inference /work/gatins/hci_genome/processing/medaka/medaka_align.bam.bam /work/gatins/hci_genome/processing/medaka/results_inference/contigs554-756.hdf --model r1041_e82_400bps_sup_v5.0.0 --threads 2 --region contig_554 contig_578 contig_579 contig_580 contig_581 contig_582 contig_583 contig_585 contig_587 contig_602 contig_683 contig_693 contig_737 contig_748 contig_749 contig_750 contig_751 contig_752 contig_755 contig_756
