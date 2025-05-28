#!/bin/bash
#SBATCH -J add_fish_tolib_build                       # Job name
#SBATCH -p short                            # Partition
#SBATCH -N 1                                # Number of nodes
#SBATCH -n 40                               # Number of tasks/threads
#SBATCH --mem=500G                            # Memory
#SBATCH -o output_%j.txt                    # Standard output file
#SBATCH -e error_%j.txt                     # Standard error file
#SBATCH --mail-user=hughes.annab@northeastern.edu  # Email
#SBATCH --mail-type=BEGIN,END                     # Email notification at job completion
#SBATCH --time=48:00:00                     # Maximum run time

module load gcc/9.2.0
# Pinfish
/work/gatins/hci_genome/kraken2/kraken2-build --threads 20 --add-to-library /work/gatins/hci_genome/processing/krakendb_fish/GCA_039737535.1_Lrho_1.0_genomic.fna --db /work/gatins/hci_genome/processing/krakendb_fish
# Atlantic cod
/work/gatins/hci_genome/kraken2/kraken2-build --threads 20 --add-to-library /work/gatins/hci_genome/processing/krakendb_fish/GCF_902167405.1_gadMor3.0_genomic.fna --db /work/gatins/hci_genome/processing/krakendb_fish
# Bicolor damselfish
/work/gatins/hci_genome/kraken2/kraken2-build --threads 20 --add-to-library /work/gatins/hci_genome/processing/krakendb_fish/GCF_000690725.1_Stegastes_partitus-1.0.2_genomic.fna --db /work/gatins/hci_genome/processing/krakendb_fish
# build db
/work/gatins/hci_genome/kraken2/kraken2-build --db /work/gatins/hci_genome/processing/krakendb_fish --build --threads 40
