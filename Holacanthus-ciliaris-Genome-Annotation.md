# Annotation of the Queen Angelfish (*H. ciliaris*) Genome

Steps:
## 1. RepeatModeler to model repetitive elements in the genome

Build database with final genome assembly
```
module load singularity
singularity run dfam-tetools-latest.sif
BuildDatabase -name hci_genome_repeats /work/gatins/hci_genome/processing/assembly_FINAL.fasta
```

Run RepeatModeler

I actually requested and received access to the long partition to make this finally run because my previous runs were taking longer than the allotted 48 hours on the short partition. In the end, my final job ended up running for only 32 hours (!??), likely because I bumped up the threads and got a bit lucky with this specific run (each round can grab different portions of the genome and some rounds may take longer in separate runs). I also ran this on Discovery even though Northeastern is trying to move everything over to Explorer because I had Singularity already installed on Discovery...

- job name: repeatmodeler
- job id: 48960035
- run time: 1-08:29:28
```
#!/bin/bash
#SBATCH -J repeatmodeler                    # Job name
#SBATCH -p long                            # Partition
#SBATCH -N 1                                # Number of nodes
#SBATCH -n 50                               # Number of tasks/threads
#SBATCH --mem=50G                           # Memory
#SBATCH -o output_%j.txt                    # Standard output file
#SBATCH -e error_%j.txt                     # Standard error file
#SBATCH --mail-user=hughes.annab@northeastern.edu  # Email
#SBATCH --mail-type=END                     # Email notification at job completion
#SBATCH --time=48:00:00                     # Maximum run time
module load singularity/3.5.3
singularity exec dfam-tetools-latest.sif RepeatModeler -LTRStruct -database hci_genome_repeats -threads 50
```

## 2. RepeatMasker to mask repetitive elements before annotating

## 3. Gene prediction with [BRAKER3](https://github.com/Gaius-Augustus/BRAKER)

Download singularity container
```
module load singularity
singularity build braker3.sif docker://teambraker/braker3:latest
singularity exec braker3.sif braker.pl
```

test:
```
singularity exec -B $PWD:$PWD braker3.sif cp /opt/BRAKER/example/singularity-tests/test1.sh .
singularity exec -B $PWD:$PWD braker3.sif cp /opt/BRAKER/example/singularity-tests/test2.sh .
singularity exec -B $PWD:$PWD braker3.sif cp /opt/BRAKER/example/singularity-tests/test3.sh .
export BRAKER_SIF=/work/gatins/hci_genome/annotation/braker3.sif # may need to modify
bash test1.sh
bash test2.sh
bash test3.sh
```

There are 4 ways to run BRAKER3: 
- [ ] 1. genome only
- [ ] 2. genome + RNA-seq data
- [ ] 3. genome + protein database (when no RNA-seq data is available)
- [ ] 4. genome + RNA-seq data + protein database

Before I get RNA-seq data, I'm going to try to annotate the genome using method 3: genome + protein database.

First, I downloaded the Vertebrata OrthoDB (I may need to adjust this input database to be more marine fish-specific, but just to see what happens I'm using this for now)
```
wget https://bioinf.uni-greifswald.de/bioinf/partitioned_odb12/Vertebrata.fa.gz
gunzip Vertebrata.fa.gz
```

Next, I had to copy over my final genome assembly to the directory I'm working in (I had to specify my working directory in the command so I guess this is why? I tried providing a path but it didn't work).

From /work/gatins/hci_genome/processing:
```
cp assembly_FINAL.fasta ../annotation/braker
```

Now, let's run BRAKER
```
singularity exec -B ${PWD}:${PWD} braker3.sif braker.pl --genome=assembly_FINAL.fasta --prot_seq=Vertebrata.fa --workingdir=/work/gatins/hci_genome/annotation/braker/ --threads 8 --skipOptimize --busco_lineage actinopterygii_odb12 &> initial_run_braker3.log
```

**I need to rerun this once I mask the genome. This is just a trial.**
