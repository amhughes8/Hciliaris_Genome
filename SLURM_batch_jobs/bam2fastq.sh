#!/bin/bash
#SBATCH -J bam2fastq                        # Job name
#SBATCH -p short                            # Partition
#SBATCH -N 2                                # Number of nodes
#SBATCH -n 2                                # Number of tasks
#SBATCH -e error_%j.txt                     # Standard error file
#SBATCH --mail-user=hughes.annab@northeastern.edu  # Email
#SBATCH --mail-type=END                     # Email notification at job completion
#SBATCH --time=48:00:00                     # Maximum run time


# Your program/command here

module load samtools/1.9
samtools bam2fq /work/gatins/hci_genome/processing/fc1.bam > /work/gatins/hci_genome/processing/hci1.fastq
samtools bam2fq /work/gatins/hci_genome/processing/fc2_SS.bam > /work/gatins/hci_genome/processing/hci2.fastq
