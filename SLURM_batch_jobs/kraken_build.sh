#!/bin/bash
#SBATCH -J kraken_build                        # Job name
#SBATCH -p short                            # Partition
#SBATCH -N 1                                # Number of nodes
#SBATCH -n 6                                # Number of tasks/threads
#SBATCH --mem=500G			    # Memory
#SBATCH -o output_%j.txt                    # Standard output file
#SBATCH -e error_%j.txt                     # Standard error file
#SBATCH --mail-user=hughes.annab@northeastern.edu  # Email
#SBATCH --mail-type=END                     # Email notification at job completion
#SBATCH --time=48:00:00                     # Maximum run time

module load gcc/9.2.0
/work/gatins/hci_genome/kraken2/kraken2-build --db /work/gatins/hci_genome/processing/kraken2_builtpython --build --threads 6
