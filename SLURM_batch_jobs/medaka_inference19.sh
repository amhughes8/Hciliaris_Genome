#!/bin/bash
#SBATCH -J medaka_inference19                        # Job name
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
medaka inference /work/gatins/hci_genome/processing/medaka/medaka_align.bam.bam /work/gatins/hci_genome/processing/medaka/results_inference/contigs1625-1765.hdf --model r1041_e82_400bps_sup_v5.0.0 --threads 2 --region contig_1625 contig_1626 contig_1630 contig_1638 contig_1641 contig_1642 contig_1651 contig_1655 contig_1656 contig_1671 contig_1673 contig_1679 contig_1680 contig_1684 contig_1687 contig_1692 contig_1693 contig_1699 contig_1701 contig_1703 contig_1704 contig_1706 contig_1710 contig_1711 contig_1713 contig_1714 contig_1715 contig_1716 contig_1736 contig_1744 contig_1745 contig_1747 contig_1748 contig_1751 contig_1753 contig_1758 contig_1759 contig_1761 contig_1763 contig_1765
