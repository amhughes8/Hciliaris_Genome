# 1. Get basecalled data from Miten
```
wget https://s3.amazonaws.com/gtl-public-data/miten/remy/03_24_25_R10_HCI_CUR_092401_GENOME.5mC_5hmC.sup.dorado.0.9.1.bam hughes.annab@xfer.discovery.neu.edu:/work/gatins/hci_genome/bams
wget https://s3.amazonaws.com/gtl-public-data/miten/remy/03_24_25_R10_HCI_CUR_092401_GENOME_SS.5mC_5hmC.sup.dorado.0.9.1.bam hughes.annab@xfer.discovery.neu.edu:/work/gatins/hci_genome/bams
```
let's store this data in an untouched 'bams' directory and also make a copy to work with in a 'processing' directory
```
cp /work/gatins/hci_genome/bams/03_24_25_R10_HCI_CUR_092401_GENOME.5mC_5hmC.sup.dorado.0.9.1.bam /work/gatins/hci_genome/processing
cp /work/gatins/hci_genome/bams/03_24_25_R10_HCI_CUR_092401_GENOME_SS.5mC_5hmC.sup.dorado.0.9.1.bam /work/gatins/hci_genome/processing
```
compressing original files to save storage space
- job name: tar_gzip_bams
- job id: 48012219
- run time:
```
tar -zcvf /work/gatins/hci_genome/HCI_CUR_092401_ONT_bams.gz /work/gatins/hci_genome/bams
```
renaming the files so they are easier to work with
```
mv 03_24_25_R10_HCI_CUR_092401_GENOME.5mC_5hmC.sup.dorado.0.9.1.bam fc1.bam
mv 03_24_25_R10_HCI_CUR_092401_GENOME_SS.5mC_5hmC.sup.dorado.0.9.1.bam fc2_SS.bam
```

# 2. samtools bam --> fastq
- job name: bam2fastq
- job id:48007076
- run time: 01:15:41
```
module load samtools/1.9
samtools bam2fq fc1.bam > hci1.fastq
samtools bam2fq fc2_SS.bam > hci2.fastq
```

## 2a. check stats from each fastq using NanoStat OR seqkit (save output in excel file)
- job name: flowcell_stats
- job id: 48011288
- run time:
```
module load anaconda3/2022.05
source activate /work/gatins/hci_genome/env
NanoStat --fastq /work/gatins/hci_genome/processing/hci1.fastq --outdir /work/gatins/hci_genome/processing/statreports --name hci1_Nanostat_fastqreport
NanoStat --fastq /work/gatins/hci_genome/processing/hci2.fastq --outdir /work/gatins/hci_genome/processing/statreports --name hci2_Nanostat_fastqreport
```
```
seqkit stat *.fastq
```
## 2b. concatenate and gzip (not sure if it'll be better to run each flow cell through whole pipeline separately or together... revisit this)
```
cat hci1.fastq hci2.fastq > hci_concat.fastq.gz
```

# 3. Porechop - trim adapters
```
porechop -i hci_concat.fastq.gz -o hci_concat_noadapters.fastq.gz
```

# 4. seqkit - filtering
```
cat hci_concat_noadapters.fastq.gz | seqkit seq -m 2000 > hci_min2000.fast.gz
seqkit stats hci_min2000.fast.gz
```
https://bioinf.shenwei.me/seqkit/usage/

# 5. Flye - assemble (maybe Hifiasm or shasta... try Flye first based on ONT recommendation)

# 6. polish (Medaka - ONT recommendation)

# 7. Blobtools - decontaminate and inspect
## 7a. BUSCO

# 8. NCBI adapter check
https://events.zoomgov.com/ej/Akmyb_uwsX0jDQtdW4EkddmRP2U7zDbJG3GwqFa2375b7pPpHMRS~A3Bfz_cM3xfiOOmt-OzkU9PoO_juK_MAF-VBZ_S2R9a6OUJSYJ-KCPNVdK2KaIQO3RZ-rva04f5PW1oUJEFC-AwyRCX9sVgyNkNgp-MLxkRBkdodc

