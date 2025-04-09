#1. Get basecalled data from Miten
```
wget https://s3.amazonaws.com/gtl-public-data/miten/remy/03_24_25_R10_HCI_CUR_092401_GENOME.5mC_5hmC.sup.dorado.0.9.1.bam hughes.annab@xfer.discovery.neu.edu:/work/gatins/hci_genome/bams
wget https://s3.amazonaws.com/gtl-public-data/miten/remy/03_24_25_R10_HCI_CUR_092401_GENOME_SS.5mC_5hmC.sup.dorado.0.9.1.bam hughes.annab@xfer.discovery.neu.edu:/work/gatins/hci_genome/bams
```
let's store this data in an untouched 'bams' directory and also make a copy to work with in a 'processing' directory
```
cp /work/gatins/hci_genome/bams/03_24_25_R10_HCI_CUR_092401_GENOME.5mC_5hmC.sup.dorado.0.9.1.bam /work/gatins/hci_genome/processing
cp /work/gatins/hci_genome/bams/03_24_25_R10_HCI_CUR_092401_GENOME_SS.5mC_5hmC.sup.dorado.0.9.1.bam /work/gatins/hci_genome/processing
```
renaming the files so they are easier to work with
```
mv 03_24_25_R10_HCI_CUR_092401_GENOME.5mC_5hmC.sup.dorado.0.9.1.bam fc1.bam
mv 03_24_25_R10_HCI_CUR_092401_GENOME_SS.5mC_5hmC.sup.dorado.0.9.1.bam fc2_SS.bam
```

#2. run summary statistics with NanoStat on sequencing_summary.txt files from both flow cells
```
NanoStat --summary flowcell1_summary.txt flowcell2_summary.txt --outdir /work/gatins/hci_genome/statreports --name Nanostat_Stat_Report
```

#3. samtools bam --> fastq
```
module load samtools/1.9
samtools bam2fq fc1.bam > hci1.fastq
samtools bam2fq fc2_SS.bam > hci2.fastq
```
** started at 4:51pm

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

#9. NCBI adapter check
https://events.zoomgov.com/ej/Akmyb_uwsX0jDQtdW4EkddmRP2U7zDbJG3GwqFa2375b7pPpHMRS~A3Bfz_cM3xfiOOmt-OzkU9PoO_juK_MAF-VBZ_S2R9a6OUJSYJ-KCPNVdK2KaIQO3RZ-rva04f5PW1oUJEFC-AwyRCX9sVgyNkNgp-MLxkRBkdodc

