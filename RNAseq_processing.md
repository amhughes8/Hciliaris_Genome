# Processing RNAseq data from Novogene

I think Novogene already performed filtering... here are the steps they list in their filtering process:

(1) Remove reads containing adapters.

(2) Remove reads containing N > 10% (N represents the base cannot be determined).

(3) Remove reads containing low quality (Qscore<= 5) base which is over 50% of the total base.

I checked a couple of samples with fastqc and they are looking ok, so let's proceed.

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
