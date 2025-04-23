#!/bin/bash
#SBATCH -J busco_hifiasm                        # Job name
#SBATCH -p short                            # Partition
#SBATCH -N 1                                # Number of nodes
#SBATCH -n 25                                # Number of tasks/threads
#SBATCH -o output_%j.txt                    # Standard output file
#SBATCH -e error_%j.txt                     # Standard error file
#SBATCH --mail-user=hughes.annab@northeastern.edu  # Email
#SBATCH --mail-type=END                     # Email notification at job completion
#SBATCH --time=48:00:00                     # Maximum run time


# Your program/command here
module load anaconda3/2022.05
source activate /work/gatins/hci_genome/busco
busco -i /work/gatins/hci_genome/processing/test_hifiasm.fa --mode genome --lineage_dataset actinopterygii_odb12 --cpu 25 --out hifiasm_busco
