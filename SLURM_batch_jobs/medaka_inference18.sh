#!/bin/bash
#SBATCH -J medaka_inference18                        # Job name
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
medaka inference /work/gatins/hci_genome/processing/medaka/medaka_align.bam.bam /work/gatins/hci_genome/processing/medaka/results_inference/contigs1484-1624.hdf --model r1041_e82_400bps_sup_v5.0.0 --threads 2 --region contig_1484 contig_1489 contig_1493 contig_1495 contig_1500 contig_1501 contig_1502 contig_1510 contig_1511 contig_1516 contig_1523 contig_1524 contig_1530 contig_1532 contig_1540 contig_1541 contig_1542 contig_1544 contig_1552 contig_1554 contig_1558 contig_1560 contig_1572 contig_1576 contig_1578 contig_1579 contig_1580 contig_1581 contig_1586 contig_1588 contig_1596 contig_1599 contig_1602 contig_1603 contig_1604 contig_1612 contig_1613 contig_1620 contig_1621 contig_1624
