# Annotation of the Queen Angelfish (*H. ciliaris*) Genome

Steps:
1. RepeatModeler to model repetitive elements in the genome

Build database with final genome assembly
```
module load singularity
singularity run dfam-tetools-latest.sif
BuildDatabase -name hci_genome_repeats /work/gatins/hci_genome/processing/assembly_FINAL.fasta
```

Run RepeatModeler
- job name: repeatmodeler
- job id: 48737549
- run time:
```
#!/bin/bash
#SBATCH -J repeatmodeler                    # Job name
#SBATCH -p short                            # Partition
#SBATCH -N 1                                # Number of nodes
#SBATCH -n 20                               # Number of tasks/threads
#SBATCH --mem=50G                           # Memory
#SBATCH -o output_%j.txt                    # Standard output file
#SBATCH -e error_%j.txt                     # Standard error file
#SBATCH --mail-user=hughes.annab@northeastern.edu  # Email
#SBATCH --mail-type=END                     # Email notification at job completion
#SBATCH --time=48:00:00                     # Maximum run time

module load singularity/3.5.3
singularity exec dfam-tetools-latest.sif RepeatModeler -LTRStruct -database hci_genome_repeats -threads 20
```

3. RepeatMasker to mask repetitive elements before annotating
4. Maker
