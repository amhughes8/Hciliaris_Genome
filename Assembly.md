# Get basecalled data from Miten
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
- run time: 02:03:15
```
tar -zcvf /work/gatins/hci_genome/HCI_CUR_092401_ONT_bams.gz /work/gatins/hci_genome/bams
```
renaming the bams I will work with for simplicity
```
mv 03_24_25_R10_HCI_CUR_092401_GENOME.5mC_5hmC.sup.dorado.0.9.1.bam fc1.bam
mv 03_24_25_R10_HCI_CUR_092401_GENOME_SS.5mC_5hmC.sup.dorado.0.9.1.bam fc2_SS.bam
```

# samtools bam --> fastq
- job name: bam2fastq
- job id:48007076
- run time: 01:15:41
```
module load samtools/1.9
samtools bam2fq fc1.bam > hci1.fastq
samtools bam2fq fc2_SS.bam > hci2.fastq
```

## check stats from each fastq using NanoStat OR seqkit (save output in excel file)
- job name: flowcell_stats
- job id: 48011288
- run time: 02:54:42
```
module load anaconda3/2022.05
source activate /work/gatins/hci_genome/env
NanoStat --fastq /work/gatins/hci_genome/processing/hci1.fastq --outdir /work/gatins/hci_genome/processing/statreports --name hci1_Nanostat_fastqreport
NanoStat --fastq /work/gatins/hci_genome/processing/hci2.fastq --outdir /work/gatins/hci_genome/processing/statreports --name hci2_Nanostat_fastqreport
```
I'm curious about the difference in outputs from NanoStat and seqkit so let's run seqkit too
```
module load anaconda3/2022.05
source activate /work/gatins/hci_genome/env
seqkit stat /work/gatins/hci_genome/processing/*.fastq
```

# Porechop - trim adapters
Porechop uses a lot of memory, so it is going to be really challenging to run on a concatenated file. I submitted a batch job on April 10 for hci1.fastq asking for 800G of memory (it is still in the queue) and I'm now interactively running Porechop on hci2.fastq with 300G of memory. **this step takes longer than 4 hours (short partition default), so make sure to indicate longer runtimes!**
```
srun --partition=short --nodes=1 --cpus-per-task=1 --mem=300G --time=48:00:00 --pty /bin/bash
porechop -i hci2.fastq -o hci2_noadapters.fastq
```
- job name: porechop1
- job id: 48049716
- run time: 20:39:38
```
module load anaconda3/2022.05
source activate /work/gatins/hci_genome/env
porechop -i /work/gatins/hci_genome/processing/hci1.fastq -o hci1_noadapters.fastq
```

# fastqc
running fastqc on each fastq before and after adapter removal ** I'm struggling to make this work. It keeps running out of memory but not alerting me... I need to figure out how to run this more efficiently since these files are all so large **
- job name: fastqc
- job id:
- run time: 
```
module load fastqc/0.11.9
fastqc /work/gatins/hci_genome/processing/*.fastq -o ./fastqc/
```

# concatenate to one big file
now that things look good after Porechop, let's concatenate before we filter and assemble
```
cat hci1_noadapters.fastq hci2_noadapters.fastq > hci_concat_noadapters.fastq
```

# [seqkit](https://bioinf.shenwei.me/seqkit/usage/) - filtering
we're filtering to a minimum sequence length of 2000 and minimum Q-score of 3 (I also allocated 10 threads to try and speed up the process). I also ran the concatenate function above with the seqkit filtering step as a batch job together:
- job name: concatenate_filter
- job id: 48175086
- run time: 00:33:29
```
module load anaconda3/2022.05 discovery
source activate /work/gatins/hci_genome/env
cat /work/gatins/hci_genome/processing/hci1_noadapters.fastq /work/gatins/hci_genome/processing/hci2_noadapters.fastq > /work/gatins/hci_genome/processing/hci_concat_noadapters.fastq
cat /work/gatins/hci_genome/processing/hci_concat_noadapters.fastq | seqkit seq -m 2000 -Q 3 -j 10 > /work/gatins/hci_genome/processing/hci_filtered_2kQ3.fastq
```
running stats on filtered output (saved output as .txt file)
```
seqkit stats hci_filtered_2kQ3.fastq -a
```

# Estimating genome size with Jellyfish (k=21)
- job name: jellyfish
- job id: 48198203
- run time:
```
module load jellyfish/2.2.10
jellyfish count -m 21 -s 100M -t 10 -C -o hci_21mer_output /work/gatins/hci_genome/processing/hci_filtered_2kQ3.fastq
```
```
jellyfish histo mer_counts.jf
```

# ASSEMBLIES
## Flye + polish with Medaka (ONT recommendation)
Flye is the method recommended by ONT for animal genome assemblies, so let's give it a shot. Requested 500G RAM and 32 threads:
- job name: assembly_flye
- job id: 48197184
- run time:
```
flye --nano-raw /work/gatins/hci_genome/processing/hci_filtered_2kQ3.fastq --threads 32 --out-dir /work/gatins/hci_genome/processing/assembly_Flye
```

## Hifiasm (newer and no polish required)

now let's check and make sure the adapters came off with FCS from NCBI
```
module load singularity/3.10.3
curl https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/releases/latest/fcs-adaptor.sif -Lo fcs-adaptor.sif
```
```
mkdir fcs_output
./run_fcsadaptor.sh --fasta-input hci_concat_noadapters.fastq.gz --output-dir /work/gatins/hci_genome/processing/fcs_output --euk --container-engine singularity --image fcs-adaptor.sif
```

# polish with Medaka
```
medaka_consensus -i reads from basecaller -d assembly_flye.fasta -o medaka_out -t 32 -m basecaller model
```
to find model list:
```
medaka tools list\_models
```

# Blobtools - decontaminate and inspect
## BUSCO

# 8. NCBI adapter check
https://events.zoomgov.com/ej/Akmyb_uwsX0jDQtdW4EkddmRP2U7zDbJG3GwqFa2375b7pPpHMRS~A3Bfz_cM3xfiOOmt-OzkU9PoO_juK_MAF-VBZ_S2R9a6OUJSYJ-KCPNVdK2KaIQO3RZ-rva04f5PW1oUJEFC-AwyRCX9sVgyNkNgp-MLxkRBkdodc

