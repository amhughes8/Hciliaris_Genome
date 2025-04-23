#!/bin/bash
#SBATCH -J medaka_inference21                        # Job name
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
medaka inference /work/gatins/hci_genome/processing/medaka/medaka_align.bam.bam /work/gatins/hci_genome/processing/medaka/results_inference/contigs1921-2088.hdf --model r1041_e82_400bps_sup_v5.0.0 --threads 2 --region contig_1921 contig_1924 contig_1932 contig_1933 contig_1942 contig_1956 contig_1958 contig_1967 contig_1968 contig_1969 contig_1972 contig_1974 contig_1983 contig_1984 contig_1987 contig_1995 contig_1996 contig_1999 contig_2000 contig_2004  contig_2008 contig_2009 contig_2011 contig_2014 contig_2019 contig_2024 contig_2029 contig_2032 contig_2033 contig_2046 contig_2049 contig_2057 contig_2069 contig_2070 contig_2071 contig_2073 contig_2074 contig_2082 contig_2083 contig_2088
