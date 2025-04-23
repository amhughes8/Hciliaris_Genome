#!/bin/bash
#SBATCH -J medaka_inference20                        # Job name
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
medaka inference /work/gatins/hci_genome/processing/medaka/medaka_align.bam.bam /work/gatins/hci_genome/processing/medaka/results_inference/contigs1766-1916.hdf --model r1041_e82_400bps_sup_v5.0.0 --threads 2 --region contig_1766 contig_1768 contig_1769 contig_1772 contig_1780 contig_1781 contig_1784 contig_1786 contig_1788 contig_1789 contig_1791 contig_1792 contig_1794 contig_1801 contig_1803 contig_1808 contig_1810 contig_1811 contig_1825 contig_1827 contig_1835 contig_1837 contig_1841 contig_1843 contig_1859 contig_1860 contig_1862 contig_1864 contig_1865 contig_1866 contig_1868 contig_1870 contig_1872 contig_1895 contig_1896 contig_1899 contig_1900 contig_1903 contig_1910 contig_1916
