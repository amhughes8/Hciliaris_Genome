#!/bin/bash
#SBATCH -J stats_2.5kQ5                        # Job name
#SBATCH -p short                            # Partition
#SBATCH -N 1                                # Number of nodes
#SBATCH -n 10                                # Number of tasks/threads
#SBATCH -o output_%j.txt                    # Standard output file
#SBATCH -e error_%j.txt                     # Standard error file
#SBATCH --mail-user=hughes.annab@northeastern.edu  # Email
#SBATCH --mail-type=END                     # Email notification at job completion
#SBATCH --time=48:00:00                     # Maximum run time

module load anaconda3/2022.05
source activate /work/gatins/hci_genome/env
seqkit stats hci_filtered_2.5kQ5.fastq -a -j 10
