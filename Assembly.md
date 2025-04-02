1. Get data from Miten
2. samtools bam --> fastq
  2a. concatenate if needed (not sure if it'll be better to run each flow cell through whole pipeline separately or together...?)
3. Porechop - trim adapters
from Remy's pipeline:
``` porechop -i HPA_fastq/HPA_big_01.fastq.gz -o HPA_fastq/HPA_big_01_pc.fastq.gz --threads 24```
4. seqkit - filtering (Remy used Nanofilt but i think seqkit looks ok)
`cat filename | seqkit seq -m 2000 | seqkit stats`
https://bioinf.shenwei.me/seqkit/usage/
5. Hifiasm - assemble (maybe Flye or shasta... try Hifiasm first)
  5a. estimate genome size with jellyfish (use Q5 reads)
6. Blobtools - decontaminate and inspect
  6a. BUSCO
7. polish
