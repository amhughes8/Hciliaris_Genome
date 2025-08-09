# Processing RNAseq data from Novogene

I think Novogene already performed filtering... here are the steps they list in their filtering process:

(1) Remove reads containing adapters.

(2) Remove reads containing N > 10% (N represents the base cannot be determined).

(3) Remove reads containing low quality (Qscore<= 5) base which is over 50% of the total base.

I checked a couple of samples with fastqc and it seems like the adapters are gone but they still have PolyA contamination?
![plot](photos/polyA_liver.png)

I'm going to use cutadapt to try and remove this. First, I need to create a conda environment for cutadapt:
```
# Create and activate conda env for cutadapt
conda create --prefix=/projects/gatins/programs_explorer/cutadapt python=3.13.5 anaconda
source activate /projects/gatins/programs_explorer/cutadapt

# Install cutadapt
pip install cutadapt
```

Now, let's run it. This took ~35 mins with 1 core:
```
# move to wd
cd /projects/gatins/hci_genome/rnaseq/fastqs

# make 'files_all'
ls *.fq.gz > files_all
sed -i -e 's/.fq.gz//g' files_all

# activate conda env
module load anaconda3/2024.06
source activate /projects/gatins/programs_explorer/cutadapt

# run cutadapt
for i in `cat files_all`;
  do cutadapt --poly-a -o ${i}_polyAremoved.fq.gz ${i}.fq.gz;
  done
```

Next, I'm going to trim the ends with trimgalore:
```
cd /projects/gatins/hci_genome/rnaseq/fastqs

# activate conda env
module load anaconda3/2024.06
source activate /projects/gatins/programs/trimgalore_ex

# make a list of all file names without extensions
ls *_1.fq.gz > files
sed -i -e 's/_1.fq.gz//g' files

# hard trim
for i in `cat files`;
  do trim_galore --fastqc --hardtrim3 142 -o /projects/gatins/hci_genome/rnaseq/fastqs/trimmed --paired --cores 2 ${i}_1_polyAremoved.fq.gz ${i}_2_polyAremoved.fq.gz;
  done
```

FastQC all files:
```
pwd
/projects/gatins/hci_genome/rnaseq/fastqs/trimmed

# name all files
ls *.fq.gz > all_trimmed

# load fastqc
module load fastqc/0.12.1

# run fastqc
for i in `cat all_trimmed`;
  do fastqc ${i} -o /projects/gatins/hci_genome/rnaseq/fastqs/fastqc/trimmed;
  done
```

## Mapping RNAseq data to reference genome
Index the genome -- total time = 00:04:50
```
cd /projects/gatins/hci_genome/rnaseq
hisat2-build -p 20 /projects/gatins/hci_genome/processing/assembly_FINAL.fasta.masked HCI_masked
```

Map
```
# making SAM files in the scratch directory
cd /hughes.annab/scratch

# -x indicates the reference genome index. hisat2 looks for the specified index first in the current directory, then in the directory specified in the HISAT2_INDEXES environment variable.
export HISAT2_INDEXES=/projects/gatins/hci_genome/rnaseq

# map to genome and create SAM file
for i in `cat files`; do hisat2 -x HCI_masked -1 ${i}_1.fq.gz -2 ${i}_2.fq.gz -S $i.sam; done

# load modules
module load samtools/1.21

# convert SAM to BAM
for i in `cat files`; do samtools view -u $i.sam | samtools sort -o $i.bam; done

# merge all sample BAM files
samtools merge -@ 32 hci_all_trimmed_rnaseq.bam ./*bam
```
