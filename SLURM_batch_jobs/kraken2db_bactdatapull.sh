#!/bin/bash
#SBATCH -J bacterial_data_kraken2                        # Job name
#SBATCH -p short                            # Partition
#SBATCH -N 1                                # Number of nodes
#SBATCH -n 32                                # Number of tasks/threads
#SBATCH --mem=200G                          # Memory
#SBATCH -o output_%j.txt                    # Standard output file
#SBATCH -e error_%j.txt                     # Standard error file
#SBATCH --mail-user=hughes.annab@northeastern.edu  # Email
#SBATCH --mail-type=END                     # Email notification at job completion
#SBATCH --time=48:00:00                     # Maximum run time


# Your program/command here
module load python/3.8.1
module load anaconda3/2022.05
source activate /work/gatins/hci_genome/env
cd /work/gatins/hci_genome/processing/kraken2_builtpython
python /work/gatins/hci_genome/kraken2/download_domain.py --domain bacteria --complete True --ext dna
