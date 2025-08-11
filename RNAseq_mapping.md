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
pwd
/projects/gatins/hci_genome/rnaseq

# -x indicates the reference genome index. hisat2 looks for the specified index first in the current directory, then in the directory specified in the HISAT2_INDEXES environment variable.
export HISAT2_INDEXES=/projects/gatins/hci_genome/rnaseq

# map to genome and create SAM file
cd /projects/gatins/hci_genome/rnaseq/fastqs/trimmed
for i in `cat files`; do hisat2 -x HCI_masked -1 ${i}_1_polyAremoved.142bp_3prime.fq.gz -2 ${i}_2_polyAremoved.142bp_3prime.fq.gz -S $i.sam; done

# load modules
module load samtools/1.21

# convert SAM to BAM
for i in `cat files`; do samtools view -u $i.sam | samtools sort -o $i.bam; done

# merge all sample BAM files
samtools merge -@ 32 hci_all_trimmed_rnaseq.bam ./*bam
```
Mapping has not improved following trimming steps. I think this is because there are so many short sequences (for example, here are the errors I've gotten only when mapping the trimmed reads):
```
Warning: skipping mate #1 of read 'LH00328:636:22V3VWLT4:8:2498:39476:29417 1:N:0:TAATGCCGAC+AGCGTTCTTG' because it was < 2 characters long
Warning: skipping mate #2 of read 'LH00328:636:22V3VWLT4:8:2498:7626:29431 2:N:0:TAATGCCGAC+AGCGTTCTTG' because length (0) <= # seed mismatches (0)
Warning: skipping mate #2 of read 'LH00328:636:22V3VWLT4:8:2498:7626:29431 2:N:0:TAATGCCGAC+AGCGTTCTTG' because it was < 2 characters long
Warning: skipping mate #1 of read 'LH00328:636:22V3VWLT4:8:2498:8564:29431 1:N:0:TAATGCCGAC+AGCGTTCTTG' because length (0) <= # seed mismatches (0)
Warning: skipping mate #1 of read 'LH00328:636:22V3VWLT4:8:2498:8564:29431 1:N:0:TAATGCCGAC+AGCGTTCTTG' because it was < 2 characters long
Warning: skipping mate #2 of read 'LH00328:636:22V3VWLT4:8:2498:20602:29431 2:N:0:TAATGCCGAC+AGCGTTCTTG' because length (0) <= # seed mismatches (0)
Warning: skipping mate #2 of read 'LH00328:636:22V3VWLT4:8:2498:20602:29431 2:N:0:TAATGCCGAC+AGCGTTCTTG' because it was < 2 characters long
Warning: skipping mate #2 of read 'LH00328:636:22V3VWLT4:8:2498:22722:29431 2:N:0:TAATGCCGAC+AGCGTTCTTG' because length (0) <= # seed mismatches (0)
Warning: skipping mate #2 of read 'LH00328:636:22V3VWLT4:8:2498:22722:29431 2:N:0:TAATGCCGAC+AGCGTTCTTG' because it was < 2 characters long
Warning: skipping mate #2 of read 'LH00328:636:22V3VWLT4:8:2498:8621:29417 2:N:0:TAATGCCGAC+AGCGTTCTTG' because it was < 2 characters long
Warning: skipping mate #1 of read 'LH00328:636:22V3VWLT4:8:2498:41846:29431 1:N:0:TAATGCCGAC+AGCGTTCTTG' because length (0) <= # seed mismatches (0)
```
Let's try using TrimGalore to remove short reads. Here is a blurb from there manual:
>## Paired-End Data
>Note that it is not recommended to remove too-short sequences if the analysed FastQ file is one of a pair of paired-end files, since this confuses the sequence-by-sequence order of paired-end reads which is again required by many aligners. For paired-end files, Trim Galore! has an option --paired which runs a paired-end validation on both trimmed _1 and _2 FastQ files once the trimming has completed. This step removes entire read pairs if at least one of the two sequences became shorter than a certain threshold. If only one of the two reads is longer than the set threshold, e.g. when one read has very poor qualities throughout, this singleton read can be written out to unpaired files (see option retain_unpaired) which may be aligned in a single-end manner.

```
# activate conda env
module load anaconda3/2024.06
source activate /projects/gatins/programs/trimgalore_ex

# remove sequences shorter than 20 bp (default)
trim_galore --length 20 --paired --cores 5 \
-o /projects/gatins/hci_genome/rnaseq/fastqs/trimmed/trimmed_length \
BRAIN_RNA_1_polyAremoved.142bp_3prime.fq.gz BRAIN_RNA_2_polyAremoved.142bp_3prime.fq.gz
