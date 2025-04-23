#!/bin/bash
#SBATCH -J medaka_inference22                        # Job name
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
medaka inference /work/gatins/hci_genome/processing/medaka/medaka_align.bam.bam /work/gatins/hci_genome/processing/medaka/results_inference/contigs2090-2235.hdf --model r1041_e82_400bps_sup_v5.0.0 --threads 2 --region contig_2090 contig_2091 contig_2092 contig_2095 contig_2101 contig_2102 contig_2106 contig_2109 contig_2113 contig_2123 contig_2124 contig_2126 contig_2128 contig_2139 contig_2147 contig_2148 contig_2151 contig_2154 contig_2156 contig_2174 contig_2175 contig_2176 contig_2178 contig_2179 contig_2180 contig_2181 contig_2184 contig_2185 contig_2186 contig_2187 contig_2198 contig_2200 contig_2201 contig_2205 contig_2206 contig_2209 contig_2212 contig_2213 contig_2218 contig_2235
