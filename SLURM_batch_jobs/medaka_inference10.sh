#!/bin/bash
#SBATCH -J medaka_inference10                        # Job name
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
medaka inference /work/gatins/hci_genome/processing/medaka/medaka_align.bam.bam /work/gatins/hci_genome/processing/medaka/results_inference/contigs757-807.hdf --model r1041_e82_400bps_sup_v5.0.0 --threads 2 --region contig_757 contig_763 contig_764 contig_765 contig_766 contig_769 contig_773 contig_774 contig_775 contig_779 contig_780 contig_784 contig_787 contig_789 contig_790 contig_792 contig_800 contig_804 contig_805 contig_807
