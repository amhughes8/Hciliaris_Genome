# PSMC Analysis: Inferring Population Size History using the Pairwise Sequentially Markovian Coalescent (PSMC) Model

## Download PSMC
```
git clone https://github.com/lh3/psmc.git
make; (cd utils; make)
```

## Prepare genome assembly data
I'm going to generate a highly-accurate dataset from my original fastq output by filtering to a **minimum length of 3000bp** and a **minimum quality score of Q10**. PSMC relies on accurate heterozygosity data, so we want the reads mapped to our reference to be as informative as possible with very little error.
- job name: filter_3kQ10
- job id: 48459438
- run time: 00:17:50
```
module load anaconda3/2022.05 discovery
source activate /work/gatins/hci_genome/env
cat /work/gatins/hci_genome/processing/hci_concat_noadapters.fastq | seqkit seq -m 3000 -Q 10 -j 10 > /work/gatins/hci_genome/processing/hci_filtered_3kQ10.fastq
```
Next, we can align these reads to our reference
```
module load minimap2/2.26
minimap2 -t 30 -ax map-ont assembly_no_contaminants.fasta hci_filtered_3kQ10.fastq > HCI_aligned.sam
```
Use samtools to convert sam->bam
```
samtools view -Sb -@ 30 -o HCI_aligned.bam HCI_aligned.sam
```
Sort bam file
```
samtools sort -o HCI_aligned_sorted.bam -O bam -@ 20 HCI_aligned.bam
```
Finally, index sorted bam file
```
samtools index -b -@ 20 HCI_aligned_sorted.bam
```

## Generate the whole-genome diploid consensus sequence for input
```
module load samtools/1.9
bcftools/1.21
samtools mpileup -C50 -f reference.fa aln.bam | bcftools view -c --threads 10 | vcfutils.pl vcf2fq -d 50 -D 300 | gzip > diploid_HCI_50_300.fq.gz
```
- -d sets minimum read depth and is recommended to be set to 1/3 of average read depth (in this case 50)
- -D sets the maximum read depth and is recommended to be set to 2x the average read depth (in this case 300)

## Run PSMC
```
/work/gatins/hci_genome/PSMC/psmc/utils/fq2psmcfa -q20 diploid_HCI_50_300.fq.gz > diploid_HCI.psmcfa
/work/gatins/hci_genome/PSMC/psmc/psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o diploid_HCI.psmc diploid_HCI.psmcfa
/work/gatins/hci_genome/PSMC/psmc/utils/psmc2history.pl diploid_HCI.psmc | /work/gatins/hci_genome/PSMC/psmc/utils/history2ms.pl > ms-cmd.sh
/work/gatins/hci_genome/PSMC/psmc/utils/psmc_plot.pl diploid diploid_HCI.psmc
```

