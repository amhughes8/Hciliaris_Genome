#!/bin/bash
#SBATCH -J kraken_addtolib                        # Job name
#SBATCH -p short                            # Partition
#SBATCH -N 1                                # Number of nodes
#SBATCH -n 10                                # Number of tasks/threads
#SBATCH -o output_%j.txt                    # Standard output file
#SBATCH -e error_%j.txt                     # Standard error file
#SBATCH --mail-user=hughes.annab@northeastern.edu  # Email
#SBATCH --mail-type=END                     # Email notification at job completion
#SBATCH --time=48:00:00                     # Maximum run time


# Your program/command here
module load gcc/9.2.0
find /work/gatins/hci_genome/processing/kraken2_builtpython/genomes/ -name '*.fna' -print0 | xargs -0 -I{} -n1 -P10 /work/gatins/hci_genome/kraken2/kraken2-build --add-to-library {} --db /work/gatins/hci_genome/processing/kraken2_builtpython
