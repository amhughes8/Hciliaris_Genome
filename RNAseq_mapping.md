# Processing RNAseq data from Novogene

I think Novogene already performed filtering... here are the steps they list in their filtering process:

(1) Remove reads containing adapters.

(2) Remove reads containing N > 10% (N represents the base cannot be determined).

(3) Remove reads containing low quality (Qscore<= 5) base which is over 50% of the total base.

I checked a couple of samples with fastqc and it seems like the adapters are gone but they still have PolyA contamination?
![plot](photos/PolyA_liver.png)

I'm going to use cutadapt to try and remove this.
```
# Create and activate conda env for cutadapt
conda create --prefix=/projects/gatins/programs_explorer/cutadapt python=3.13.5 anaconda
source activate /projects/gatins/programs_explorer/cutadapt

# Install cutadapt
conda install bioconda::cutadapt

# Trim Poly-A
cutadapt --poly-a 
```

## Mapping RNAseq data to reference genome
Index the genome -- total time = 00:04:50
```
cd /projects/gatins/hci_genome/rnaseq
hisat2-build -p 20 /projects/gatins/hci_genome/processing/assembly_FINAL.fasta.masked HCI_masked
```

Map
```
cd /projects/gatins/hci_genome/rnaseq/HCI_CUR_092401_RNAseq/01.RawData
hisat2 -x HCI_masked -q 20 \
-1 ./BRAIN_RNA/BRAIN_RNA_1.fq.gz,./FIN_RNA/FIN_RNA_1.fq.gz,./GILL_RNA/GILL_RNA_1.fq.gz,./LIVER_RNA/LIVER_RNA_1.fq.gz,./MUSCLE_RNA/MUSCLE_RNA_1.fq.gz \
-2 ./BRAIN_RNA/BRAIN_RNA_2.fq.gz,./FIN_RNA/FIN_RNA_2.fq.gz,./GILL_RNA/GILL_RNA_2.fq.gz,./LIVER_RNA/LIVER_RNA_2.fq.gz,./MUSCLE_RNA/MUSCLE_RNA_2.fq.gz \
-S HCI_mapped_RNA.sam
```
okay this did not work because i need to run hisat2 for each tissue type separately i think... running the following instead:
```
for i in /projects/gatins/hci_genome/rnaseq/fastqs; do
hisat2 -x HCI_masked -1 ${i}_1.fastq -2 ${i}_2.fastq -S HCI_mapped_RNA.sam
done
```
this didn't work because i guess my first line is indicating that the command should be looking for fastqs_1.fastq? idk, trying again

based on Jen's code for red sea urchin:
```
cd /projects/gatins/hci_genome/rnaseq/fastqs
ls *_1.fq.gz > files
sed -i -e 's/_1.fq.gz//g' files
export HISAT2_INDEXES=/projects/gatins/hci_genome/rnaseq
for i in `cat files`; do hisat2 -x HCI_masked -1 ${i}_1.fq.gz -2 ${i}_2.fq.gz -S $i.sam; done
```

OK this is working! Now from within /projects/gatins/hci_genome/rnaseq/fastqs
```
# Convert SAM to BAM
module load samtools/1.21
for i in `cat files`; do samtools view -u $i.sam | samtools sort -o $i.bam; done
```
```
# Merge all sample BAM files
cd ../
samtools merge -@ 32 hci_all-rnaseq.bam ./fastqs/*bam
```
