#1. Get basecalled data from Miten

#2. run summary statistics with NanoStat on sequencing_summary.txt files from both flow cells
```
NanoStat --summary flowcell1_summary.txt flowcell2_summary.txt --outdir /work/gatins/hci_genome/statreports --name Nanostat_Stat_Report
```

#3. samtools bam --> fastq
```
module load samtools/1.9
samtools bam2fq Output_From_Miten1.bam > hci1.fastq
samtools bam2fq Output_From_Miten2.bam > hci2.fastq
```
##3a. check stats from each fastq using NanoStat OR seqkit (save output in excel file)
```
NanoStat --fastq hci1.fastq --outdir /work/gatins/hci_genome/statreports --name hci1_Nanostat_fastqreport
NanoStat --fastq hci2.fastq --outdir /work/gatins/hci_genome/statreports --name hci2_Nanostat_fastqreport
```
```
seqkit stat *.fastq
```
##3a. concatenate and gzip (not sure if it'll be better to run each flow cell through whole pipeline separately or together... revisit this)
```
cat hci1.fastq hci2.fastq > hci_concat.fastq.gz
```

#4. Porechop - trim adapters
```
porechop -i hci_concat.fastq.gz -o hci_concat_noadapters.fastq.gz
```

#5. seqkit - filtering
```
cat hci_concat_noadapters.fastq.gz | seqkit seq -m 2000 > hci_min2000.fast.gz
seqkit stats hci_min2000.fast.gz
```
https://bioinf.shenwei.me/seqkit/usage/

#6. Flye - assemble (maybe Hifiasm or shasta... try Flye first based on ONT recommendation)

#7. polish (Medaka - ONT recommendation)

#8. Blobtools - decontaminate and inspect
##8a. BUSCO

