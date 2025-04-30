# Removal of Mitochondrial DNA Before Genome Assembly

## Download mitochondrial sequence from NCBI
>NC_027595.1_HCI_mito.fasta

## Use minimap2 to map filtered fastq reads pre-assembly (hci_filtered_2.5kQ5.fastq) to reference mitochondrial sequence (NC_027595.1_HCI_mito.fasta)

```
module load minimap2/2.26
minimap2 -t 40 -ax map-ont NC_027595.1_HCI_mito.fasta hci_filtered_2.5kQ5.fastq > aln_minimap2.sam
```
## Use samtools to extract sequences that mapped to mitochondrial genome
Convert to BAM
```
module load samtools/1.19.2
samtools view -Sb -@ 30 aln_minimap2.sam > mito_aln.bam
```
Sort BAM
```
samtools sort -@ 20 mito_aln.bam -o mito_aln.sorted.bam
```
Index BAM
```
samtools index mito_aln.sorted.bam
```
Extract sequences that were unmapped and save them to a new BAM file
```
samtools view -b -f 4 -@ 20 mito_aln.sorted.bam > unmapped.bam
```
Convert BAM to FASTQ
```
samtools fastq unmapped.bam > reads_no_mito.fastq
```

## Re-assemble with hifiasm (32 threads, 500G memory)
```
module load anaconda3/2022.05
source activate /work/gatins/hci_genome/env
hifiasm -o assembly_hifiasm_no_mito.asm --ont -t32 /work/gatins/hci_genome/processing/mtdna/removal/reads_no_mito.fastq
```
