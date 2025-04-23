#!/bin/bash
#SBATCH -J medaka_inference11                        # Job name
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
medaka inference /work/gatins/hci_genome/processing/medaka/medaka_align.bam.bam /work/gatins/hci_genome/processing/medaka/results_inference/contigs809-919.hdf --model r1041_e82_400bps_sup_v5.0.0 --threads 2 --region contig_809 contig_810 contig_820 contig_822 contig_831 contig_852 contig_856 contig_860 contig_863 contig_880 contig_884 contig_895 contig_911 contig_912 contig_914 contig_915 contig_916 contig_917 contig_918 contig_919
