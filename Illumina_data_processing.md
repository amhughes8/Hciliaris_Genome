# Processing Illumina WGS data

Remy demultiplexed the files from her large BSB WGS run, but I am now going to quality trim them.

First, let's check them with fastqc
```
module load fastqc/0.12.1
```

```
source activate /projects/gatins/programs/trimgalore_ex

#hardtrim3 will keep sequence from the 3' end (thus removing bp's from the 5' end)
#all my sequences are 150bp's
trim_galore --fastqc --hardtrim3 142 -o /projects/gatins/hci_genome/illumina/trimmed --paired --cores 2 HCI_CUR_092401_merged.1.fq.gz HCI_CUR_092401_merged.2.fq.gz

#now remove remnant adapters and quality trim
trim_galore --fastqc -o /projects/gatins/hci_genome/illumina/clean/trimgalore_hard --paired --cores 2 /projects/gatins/hci_genome/illumina/trimmed/HCI_CUR_092401_merged.1.142bp_3prime.fq.gz /projects/gatins/hci_genome/illumina/trimmed/HCI_CUR_092401_merged.2.142bp_3prime.fq.gz
```

run multiqc from /projects/gatins/hci_genome/illumina/clean/trimgalore_hard

```
source activate /projects/gatins/programs_explorer
multiqc .
```
