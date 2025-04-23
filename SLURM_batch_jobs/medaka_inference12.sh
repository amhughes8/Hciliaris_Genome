#!/bin/bash
#SBATCH -J medaka_inference12                        # Job name
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
medaka inference /work/gatins/hci_genome/processing/medaka/medaka_align.bam.bam /work/gatins/hci_genome/processing/medaka/results_inference/contigs920-1007.hdf --model r1041_e82_400bps_sup_v5.0.0 --threads 2 --region contig_920 contig_926 contig_927 contig_928 contig_937 contig_940 contig_941 contig_943 contig_948 contig_954 contig_968 contig_977 contig_980 contig_987 contig_991 contig_993 contig_994 contig_995 contig_1005 contig_1007
