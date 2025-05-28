#!/bin/bash
#SBATCH -J convert_bam_sort_index                        # Job name
#SBATCH -p short                            # Partition
#SBATCH -N 1                                # Number of nodes
#SBATCH -n 30                                # Number of tasks/threads
#SBATCH -o output_%j.txt                    # Standard output file
#SBATCH -e error_%j.txt                     # Standard error file
#SBATCH --mail-user=hughes.annab@northeastern.edu  # Email
#SBATCH --mail-type=END                     # Email notification at job completion
#SBATCH --time=48:00:00                     # Maximum run time

module load samtools/1.19.2
samtools view -Sb -@ 30 aln_minimap2.sam > mtdna_aligned.bam
samtools sort -@ 20 mtdna_aligned.bam -o mito_aln.sorted.bam
samtools index mito_aln.sorted.bam
