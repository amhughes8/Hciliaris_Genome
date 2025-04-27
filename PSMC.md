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
module load samtools/1.9
samtools sort -o HCI_aligned_sorted.bam -O bam -@ 20 HCI_aligned.bam
```
Finally, index sorted bam file
```
samtools index -b -@ 20 HCI_aligned_sorted.bam
```

## Generate the whole-genome diploid consensus sequence for input
```
module load samtools/1.9
module load bcftools/1.21
bcftools mpileup -C50 -f assembly_no_contaminants.fasta HCI_aligned_sorted.bam | bcftools view -c --threads 10 | vcfutils.pl vcf2fq -d 50 -D 300 | gzip > diploid_HCI_50_300.fq.gz
```
- -d sets minimum read depth and is recommended to be set to 1/3 of average read depth (in this case 50)
- -D sets the maximum read depth and is recommended to be set to 2x the average read depth (in this case 300)

This script was designed to take a BAM and a reference assembly to create a VCF with genotype likelihoods for a diploid individual --> consensus call --> translate VCF to FASTQ --> gzip.

It seems that this pipeline (specifically vcfutils.pl) has been deprecated and people have recently switched over to using *samtools consensus* to get a .fq.gz consensus sequence for the PSMC analysis. Let's try it:
```
module load samtools/1.19.2
samtools consensus --ambig -f fastq -d 50 /work/gatins/hci_genome/processing/HCI_aligned_sorted.bam -o consensus.fq
```
now gzip
```
gzip consensus.fq
```
output:
```
consensus.fq.gz
```

## Run PSMC
First, convert diploid FASTQ into a psmcfa file:
```
/work/gatins/hci_genome/PSMC/psmc/utils/fq2psmcfa -q20 consensus.fq.gz > diploid_HCI.psmcfa
```
Now, run PSMC:
```
/work/gatins/hci_genome/PSMC/psmc/psmc -N30 -t30 -r5 -p "4+30*2+4+6+10" -o diploid_HCI.psmc diploid_HCI.psmcfa
```
PSMC parameters:
- p STR pattern of parameters [4+5*3+4]
- t FLOAT maximum 2N0 coalescent time [15]
- N INT maximum number of iterations [30]
- r FLOAT initial theta/rho ratio [4]
- o FILE output file [stdout]

Just doing a quick test run, so using Remy's input (mutation rate=u, generation time in years=g)
```
/work/gatins/hci_genome/PSMC/psmc/utils/psmc_plot.pl -u 1e-08 -g 5 HCI_t30r5_plot_u1-8g5 diploid_HCI.psmc
```
Output:
```
HCI_t30r5_plot_u1-8g5.eps  HCI_t30r5_plot_u1-8g5.par
```
Download files to computer and visualize:

# Bootstrapping
Follow all same steps as above and use your .psmc file for bootstrapping.

Use the splitfa command to split long chromosome sequences found in diploid_HCI.psmcfa file to shorter segments for bootstrapping.
```
/work/gatins/hci_genome/PSMC/psmc/utils/splitfa diploid_HCI.psmcfa > diploid_HCI_split.psmcfa
```
Once you have your diploid_split.psmcfa file you will need to copy this file into 100 independent files. I personally like to do this in a separate directory, so mkdir bootstrapping. Now, copy (or move) diploid_HCI_split.psmcfa and your original psmc run outfile into your new bootstrap directory. The psmc file will be used after you run the bootstrap to concatenate with the other output files.
```
mkdir bootstrapping
cp diploid_HCI_split.psmcfa bootstrapping
cp diploid_HCI.psmc bootstrapping
```
