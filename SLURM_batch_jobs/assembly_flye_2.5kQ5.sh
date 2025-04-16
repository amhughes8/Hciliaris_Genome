#!/bin/bash
#SBATCH -J assembly_flye_2.5kQ5                        # Job name
#SBATCH -p short                            # Partition
#SBATCH -N 1                                # Number of nodes
#SBATCH -n 32                                # Number of tasks/threads
#SBATCH --mem=500G                          # Memory
#SBATCH -o output_%j.txt                    # Standard output file
#SBATCH -e error_%j.txt                     # Standard error file
#SBATCH --mail-user=hughes.annab@northeastern.edu  # Email
#SBATCH --mail-type=BEGIN,END                     # Email notification at job start and completion
#SBATCH --time=48:00:00                     # Maximum run time


# Your program/command here
module load anaconda3/2022.05 discovery
source activate /work/gatins/hci_genome/env
flye --nano-raw /work/gatins/hci_genome/processing/hci_filtered_2.5kQ5.fastq --threads 32 --out-dir /work/gatins/hci_genome/processing/assembly_Flye_2.5kQ5
