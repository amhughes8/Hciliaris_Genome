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

Working directory:  /work/gatins/hci_genome/annotation/RM_44655.FriJun131036432025
RepeatModeler output:
```
hci_genome_repeats-families.fa  #Consensus sequences for each family identified.
hci_genome_repeats-families.stk  #Seed alignments for each family identified.
hci_genome_repeats-rmod.log  #Execution log.  Useful for reproducing results.
```

Following [this tutorial](https://darencard.net/blog/2022-07-09-genome-repeat-annotation/), we will now use the output FASTA file to get a better understanding of how many sequences were known vs unknown:
```
cat hci_genome_repeats-families.fa | seqkit fx2tab | grep -v "Unknown" | seqkit tab2fx > hci_genome_repeats-families.fa.known
cat hci_genome_repeats-families.fa | seqkit fx2tab | grep "Unknown" | seqkit tab2fx > hci_genome_repeats-families.fa.unknown
# quantify number of classified elements
grep -c ">" hci_genome_repeats-families.fa.known #811
# quantify number of unknown elements
grep -c ">" hci_genome_repeats-families.fa.unknown #1204
```

## 2. RepeatMasker to mask repetitive elements before annotating

I created a RepeatMasker conda environment: /projects/gatins/programs_explorer/RepeatMasker.

This environment really only contains h5py and all other requirements as well as the program itself are located in /projects/gatins/hci_genome/annotation.

To run RepeatMasker:
```
/projects/gatins/hci_genome/annotation/RepeatMasker/RepeatMasker -pa 10 -lib hci_genome_repeats-families.fa -xsmall -gff /projects/gatins/hci_genome/processing/assembly_FINAL.fasta
```
-pa: indicates the number of parallel sequence batch jobs running against the database at a time. **Since I designated 10 parallel runs and am using the RMBlast database, I am running this job with 40 threads (10*40).**
-xsmall: soft masks
-gff: creates an additional Gene Feature Finding format output

RepeatMasker run time: ~3 days
Output:
```
[hughes.annab@explorer-01 processing]$ cat assembly_FINAL.fasta.tbl
==================================================
file name: assembly_FINAL.fasta
sequences:           102
total length:  604914273 bp  (604914273 bp excl N/X-runs)
GC level:         41.41 %
bases masked:  114112123 bp ( 18.86 %)
==================================================
               number of      length   percentage
               elements*    occupied  of sequence
--------------------------------------------------
Retroelements        95872     18091213 bp    2.99 %
   SINEs:            13898      1507537 bp    0.25 %
   Penelope:          5482       453553 bp    0.07 %
   LINEs:            57154     11048343 bp    1.83 %
    CRE/SLACS            0            0 bp    0.00 %
     L2/CR1/Rex      43215      8039610 bp    1.33 %
     R1/LOA/Jockey    1189       330372 bp    0.05 %
     R2/R4/NeSL       1033       250114 bp    0.04 %
     RTE/Bov-B        5826       811435 bp    0.13 %
     L1/CIN4          3051       924195 bp    0.15 %
   LTR elements:     19338      5081780 bp    0.84 %
     BEL/Pao          1658       898697 bp    0.15 %
     Ty1/Copia           0            0 bp    0.00 %
     Gypsy/DIRS1      2848      1715396 bp    0.28 %
       Retroviral     8539      1312047 bp    0.22 %

DNA transposons     181164     26362430 bp    4.36 %
   hobo-Activator    60136      8744061 bp    1.45 %
   Tc1-IS630-Pogo    41633      7224934 bp    1.19 %
   En-Spm                0            0 bp    0.00 %
   MULE-MuDR           461        84157 bp    0.01 %
   PiggyBac            596       102241 bp    0.02 %
   Tourist/Harbinger  8982      1625562 bp    0.27 %
   Other (Mirage,     2957       670965 bp    0.11 %
    P-element, Transib)

Rolling-circles        977       191287 bp    0.03 %

Unclassified:       412076     52186314 bp    8.63 %

Total interspersed repeats:    96639957 bp   15.98 %


Small RNA:            8549      1799805 bp    0.30 %

Satellites:           2686      1660952 bp    0.27 %
Simple repeats:     286580     12544883 bp    2.07 %
Low complexity:      35635      1990489 bp    0.33 %
==================================================

* most repeats fragmented by insertions or deletions
  have been counted as one element


RepeatMasker version 4.1.9 , default mode
run with rmblastn version 2.14.1+
The query was compared to classified sequences in "hci_genome_repeats-families.fa"
FamDB:
```

## 3. Gene prediction with [BRAKER3](https://github.com/Gaius-Augustus/BRAKER)

Download singularity container
```
module load singularity
singularity build braker3.sif docker://teambraker/braker3:latest
singularity exec braker3.sif braker.pl
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

From /projects/gatins/hci_genome/processing:
```
cp assembly_FINAL.fasta.masked ../annotation/braker
```

Now, let's run BRAKER
- I think BRAKER ran for about 24 hrs?
```
apptainer exec -B /projects/gatins/hci_genome/annotation/braker /projects/gatins/hci_genome/annotation/braker/braker3.sif braker.pl \
--genome=/projects/gatins/hci_genome/annotation/braker/assembly_FINAL.fasta.masked \
--prot_seq=/projects/gatins/hci_genome/annotation/braker/Vertebrata.fa \
--threads=10 --species=Hciliaris --softmasking --busco_lineage=actinopterygii \
--AUGUSTUS_CONFIG_PATH=/projects/gatins/hci_genome/annotation/braker/config &> hci_braker.log
```
I ended up deleting the output from this run because it seems like it may be better to not run BRAKER with a BUSCO lineage. I can't figure out how to rerun BRAKER without deleting the initial run... like it wants to override the previous run but can't so it just doesn't work? So odd. Luckily I'm just troubleshooting so if this doesn't work I can just rerun what I initially ran!

started at 13:35 on July 29 -- finished on July 30 at 10:59 
```
apptainer exec -B /projects/gatins/hci_genome/annotation/braker /projects/gatins/hci_genome/annotation/braker/braker3.sif braker.pl \
--genome=/projects/gatins/hci_genome/annotation/braker/assembly_FINAL.fasta.masked \
--prot_seq=/projects/gatins/hci_genome/annotation/braker/Vertebrata.fa \
--threads=10 --species=Hciliaris --softmasking \
--AUGUSTUS_CONFIG_PATH=/projects/gatins/hci_genome/annotation/braker/config &> hci_nobusco_braker.log
```

NOW, let's try BRAKER with RNAseq! Bumping up the threads

Started running at 15:38 on August 6 -- finished 22:17 on August 7
```
apptainer exec -B /projects/gatins/hci_genome/annotation/braker /projects/gatins/hci_genome/annotation/braker/braker3.sif braker.pl \
--genome=/projects/gatins/hci_genome/annotation/braker/assembly_FINAL.fasta.masked \
--prot_seq=/projects/gatins/hci_genome/annotation/braker/Vertebrata.fa \
--bam=/projects/gatins/hci_genome/annotation/braker/hci_all-rnaseq.bam \
--threads=30 --species=Hciliaris --softmasking \
--AUGUSTUS_CONFIG_PATH=/projects/gatins/hci_genome/annotation/braker/config &> hci_nobusco_rnaseq_braker.log
```
Total number of protein-coding genes predicted:
```
grep -c "^>" braker.aa
```
**27032**

Initial protein BUSCO run after BRAKER:
C:95.0%[S:79.3%,D:15.7%],F:0.8%,M:4.2%,n:7207
	6846	Complete BUSCOs (C)
	5712	Complete and single-copy BUSCOs (S)
	1134	Complete and duplicated BUSCOs (D)
	61	Fragmented BUSCOs (F)
	300	Missing BUSCOs (M)
	7207	Total BUSCO groups searched

Lots of duplicated BUSCOs. Let's see how filtering with TSEBRA alters this.

## 4. Transcript filtering with [TSEBRA](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-04482-0)

```
apptainer exec braker3.sif tsebra.py --help
```

TSEBRA takes a list of gene prediciton files, a list of hintfiles and a configuration file as mandatory input.

### Step 1: Filter single-exon genes out
*Downloaded default.cfg from TSEBRA github*
```
apptainer exec braker3.sif tsebra.py \
-g /projects/gatins/hci_genome/annotation/braker/braker/GeneMark-ETP/genemark.gtf \
-k /projects/gatins/hci_genome/annotation/braker/braker/Augustus/augustus.hints.gtf \
-e /projects/gatins/hci_genome/annotation/braker/braker/hintsfile.gff \
--filter_single_exon_genes \
-c /projects/gatins/hci_genome/annotation/braker/default.cfg \
-o hci_braker_filtered.gtf
```

### Step 2: Getting the longest isoform of each gene loci from different gene sets
Combines multiple gene sets and reports the transcript with the longest coding region for each cluster of overlapping transcripts (one transcript per gene loci), e.g.
```
apptainer exec braker3.sif get_longest_isoform.py --gtf hci_braker_filtered.gtf --out hci_longest_insoforms.gtf
```

### Step 3: Extract protein sequences from filtered dataset using [gffread](https://ccb.jhu.edu/software/stringtie/gff.shtml#gffread)
```
/projects/gatins/programs_explorer/gffread/bin/gffread -w hci_transcripts_li_nseg.fa -y hci_transcripts_li_nseg.aa -g /projects/gatins/hci_genome/processing/assembly_FINAL.fasta.masked hci_longest_insoforms.gtf
```
> li = longest isoform
> nseg = no single exon genes

## 5. Protein BUSCO
```
module load anaconda3/2024.06
source activate /projects/gatins/programs_explorer/busco
busco -i hci_transcripts_li_nseg.aa --mode proteins --lineage_dataset actinopterygii_odb12 --cpu 10 --out hci_filtered_busco
```
C:99.1%[S:98.7%,D:0.4%],F:0.6%,M:0.2%,n:7207
	7143	Complete BUSCOs (C)
	7114	Complete and single-copy BUSCOs (S)
	29	Complete and duplicated BUSCOs (D)
	46	Fragmented BUSCOs (F)
	18	Missing BUSCOs (M)
	7207	Total BUSCO groups searched

I'm going to try running this all over again from the braker output but exclude the step where I filter out single-exon genes. I just don't really understand this step and I want to see what impact it has on these BUSCO results.

```
# generate consensus gtf
apptainer exec braker3.sif tsebra.py \
-g /projects/gatins/hci_genome/annotation/braker/braker/GeneMark-ETP/genemark.gtf \
-k /projects/gatins/hci_genome/annotation/braker/braker/Augustus/augustus.hints.gtf \
-e /projects/gatins/hci_genome/annotation/braker/braker/hintsfile.gff \
-c /projects/gatins/hci_genome/annotation/braker/default.cfg \
-o hci_braker_filtered_2.gtf

# longest isoform
apptainer exec braker3.sif get_longest_isoform.py --gtf hci_braker_filtered_2.gtf --out hci_longest_insoforms_2.gtf

# extract protein sequences
/projects/gatins/programs_explorer/gffread/bin/gffread -w hci_transcripts_li.fa -y hci_transcripts_li.aa -g /projects/gatins/hci_genome/processing/assembly_FINAL.fasta.masked hci_longest_insoforms_2.gtf

# busco
module load anaconda3/2024.06
source activate /projects/gatins/programs_explorer/busco
busco -i hci_transcripts_li.aa --mode proteins --lineage_dataset actinopterygii_odb12 --cpu 10 --out hci_filtered_busco_2
```

C:99.1%[S:98.7%,D:0.4%],F:0.6%,M:0.2%,n:7207
	7143	Complete BUSCOs (C)
	7114	Complete and single-copy BUSCOs (S)
	29	Complete and duplicated BUSCOs (D)
	46	Fragmented BUSCOs (F)
	18	Missing BUSCOs (M)
	7207	Total BUSCO groups searched

Hmm... so it doesn't change the completeness. I'm going to move on to the InterProScan step and see what happens. Jen said she doesn't always just filter out single exon genes because they may have function so we'll see which ones map. 

Got the heart RNAseq data! Re-running BRAKER with the full RNAseq dataset:
```
apptainer exec -B /projects/gatins/hci_genome/annotation/braker /projects/gatins/hci_genome/annotation/braker/braker3.sif braker.pl \
--genome=/projects/gatins/hci_genome/annotation/braker/assembly_FINAL.fasta.masked \
--prot_seq=/projects/gatins/hci_genome/annotation/braker/Vertebrata.fa \
--bam=/projects/gatins/hci_genome/annotation/braker/hci_all_trimmed_rnaseq_withheart.bam \
--workingdir=/projects/gatins/hci_genome/annotation/braker/hci_braker_final \
--threads=30 --species=Hciliaris_full --softmasking \
--AUGUSTUS_CONFIG_PATH=/projects/gatins/hci_genome/annotation/braker/config &> hci_nobusco_rnaseq_braker.log
```

## 6. Functional annotation with [InterProScan](https://www.ebi.ac.uk/interpro/) and [EnTAP](https://entap.readthedocs.io/en/latest/Getting_Started/introduction.html)


following formatting from red sea urchin genome annotation
```
/path/to/EnTAP --runP -i /path/to/input/transcriptome \
-d /path/to/diamond/database1.dmnd,/path/to/diamond/database2.dmnd \
--entap-ini /path/to/entap_config.ini -t 8
```

