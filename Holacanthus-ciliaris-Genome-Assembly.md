# Reference Genome Assembly for the Queen Angelfish (*Holacanthus ciliaris*)

## 1. Get basecalled data from Miten
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
renaming the bams I will work with in the /processing directory for simplicity
```
mv 03_24_25_R10_HCI_CUR_092401_GENOME.5mC_5hmC.sup.dorado.0.9.1.bam fc1.bam
mv 03_24_25_R10_HCI_CUR_092401_GENOME_SS.5mC_5hmC.sup.dorado.0.9.1.bam fc2_SS.bam
```

## 2. Use samtools to convert bam --> fastq
- job name: bam2fastq
- job id:48007076
- run time: 01:15:41
```
module load samtools/1.9
samtools bam2fq fc1.bam > hci1.fastq
samtools bam2fq fc2_SS.bam > hci2.fastq
```

## 3. Check stats from each fastq 
We can use [NanoStat](https://github.com/wdecoster/nanostat) OR [SeqKit](https://bioinf.shenwei.me/seqkit/usage/). Let's test NanoStat first:
- job name: flowcell_stats
- job id: 48011288
- run time: 02:54:42
```
module load anaconda3/2022.05
source activate /work/gatins/hci_genome/env
NanoStat --fastq /work/gatins/hci_genome/processing/hci1.fastq --outdir /work/gatins/hci_genome/processing/statreports --name hci1_Nanostat_fastqreport
NanoStat --fastq /work/gatins/hci_genome/processing/hci2.fastq --outdir /work/gatins/hci_genome/processing/statreports --name hci2_Nanostat_fastqreport
```
I'm curious about the difference in outputs from NanoStat and SeqKit so let's run SeqKit too
- job name: stats_seqkit
- job id: 48021445
- run time: 00:07:34
```
module load anaconda3/2022.05
source activate /work/gatins/hci_genome/env
seqkit stat /work/gatins/hci_genome/processing/*.fastq
```

## 4. Trim Adapters with [Porechop](https://github.com/rrwick/Porechop)
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

## 5. Concatenate to one big file
now that things look good after Porechop, let's concatenate before we filter and assemble
```
cat hci1_noadapters.fastq hci2_noadapters.fastq > hci_concat_noadapters.fastq
```

## 6. Filter with SeqKit
We're filtering to a minimum sequence length of 2500 bp and a minimum quality score of Q5.
- job name: filter_2.5kQ5
- job id: 48242201
- run time: 00:17:17
```
module load anaconda3/2022.05 discovery
source activate /work/gatins/hci_genome/env
cat /work/gatins/hci_genome/processing/hci_concat_noadapters.fastq | seqkit seq -m 2500 -Q 5 -j 10 > /work/gatins/hci_genome/processing/hci_filtered_2.5kQ5.fastq
```

## 7. Estimate genome size with [Jellyfish](https://github.com/gmarcais/Jellyfish) (k=21)
- job name: jellyfishQ5
- job id: 48316419
- run time: 02:12:55
```
module load jellyfish/2.2.10
jellyfish count -m 21 -s 500M -t 10 -C -o hciQ5_21mer_output /work/gatins/hci_genome/processing/hci_filtered_2.5kQ5.fastq
```
```
jellyfish histo hciQ5_21mer_output
```
Upload this histogram file to [GenomeScope](http://genomescope.org/genomescope2.0/) to look at sequencing coverage and get a genome size estimate.
- Description: HCI_CUR_092401
- K-mer length: 21
- Ploidy: 2
- Max k-mer coverage: -1 (left at default)
- Average k-mer coverage for polyploid genome: -1 (left at default)
![plot](photos/GenomeScope_Q5.png)
![plot](photos/GenomeScope_Q5_model.png)

## 8. Removal of mitochondrial DNA before genome assembly
My first attempt at assembling left me with >60 contigs assembled that mapped to the *H. ciliaris* mitogenome. In order to decrease mismappings and false assemblies, I'm removing the mtDNA before assembling.
### Download mitochondrial sequence from NCBI
>NC_027595.1_HCI_mito.fasta

### Use minimap2 to map filtered fastq reads pre-assembly (hci_filtered_2.5kQ5.fastq) to reference mitochondrial sequence (NC_027595.1_HCI_mito.fasta)
```
module load minimap2/2.26
minimap2 -t 40 -ax map-ont NC_027595.1_HCI_mito.fasta hci_filtered_2.5kQ5.fastq > aln_minimap2.sam
```
### Use samtools to extract sequences that mapped to mitochondrial genome
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

## 9. Assemble with [hifiasm](https://github.com/chhylp123/hifiasm)
- job name: hifiasm_assembly_nomt
- job id: 48558478
- run time: 14:14:41
```
module load anaconda3/2022.05
source activate /work/gatins/hci_genome/env
hifiasm -o assembly_hifiasm_no_mito.asm --ont -t32 /work/gatins/hci_genome/processing/mtdna/removal/reads_no_mito.fastq
```
### Running [BUSCO](https://busco.ezlab.org/busco_userguide.html#tips-for-running-busco)
first, find which dataset to run Busco against
```
busco --list-datasets
```
running against the **actinopterygii_odb12** dataset
- job name: busco_hifiasm_nomito
- job id: 48580072
- run time: 00:35:52
```
module load anaconda3/2022.05
source activate /work/gatins/hci_genome/busco
busco -i /work/gatins/hci_genome/processing/assembly_Flye/assembly.fasta --mode genome --lineage_dataset actinopterygii_odb12 --cpu 25 --out initial_assembly_flye_busco
```
### SeqKit stats on assembly + BUSCO:

|  file    | format | type | num_seqs  |  sum_len | min_len   |   avg_len  |   max_len   |    Q1   |     Q2    |     Q3 | sum_gap   |  N50 | N50_num | Q20(%) | Q30(%) | AvgQual | GC(%) | sum_n | BUSCO |
|----------|---------------|------|-----------|------------|-----------|------------|-------------|---------|-----------|--------|-----------|---------|---------|--------|--------|---------|-------|-------|------|
assembly_hifiasm_no_mito.fa | FASTA  | DNA    |    164 | 606,222,924  |  3,221 | 3,696,481.2 | 31,961,345 | 7,125 | 10,976 | 82,991.5   |     0 | 25,061,566    |   11   |    0  |     0    |    0  | 41.41   |   0 | 98.8% |

## 10. Adapter removal check with [FCS](https://github.com/ncbi/fcs) from NCBI
Download:
```
module load singularity/3.10.3
curl -LO https://github.com/ncbi/fcs/raw/main/dist/run_fcsadaptor.sh
chmod 755 run_fcsadaptor.sh
curl https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/releases/latest/fcs-adaptor.sif -Lo fcs-adaptor.sif
```
```
module load singularity/3.10.3
mkdir fcs_output_hifiasm_nomito
./run_fcsadaptor.sh --fasta-input /work/gatins/hci_genome/processing/mtdna/removal/hifiasm_nomito/assembly_hifiasm_no_mito.fa --output-dir /work/gatins/hci_genome/processing/fcs/fcs_output_hifiasm_nomito --euk --container-engine singularity --image fcs-adaptor.sif
```
Output is clean! No adapters detected.

## 11. Contamination Identification with [Kraken2](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown)
As I was setting up the Kraken2 standard database following the Git instructions, I kept receiving timeout errors `Step 1/2: Performing ftp file transfer of requested files
rsync_from_ncbi.pl: FTP connection error: Net::FTP: connect: timeout`... I found [this GitHub issues page](https://github.com/DerrickWood/kraken2/issues/272) where people were running into the same problem. I'm going to take the approach recommended in this thread and use a [custom-built python script](https://github.com/R-Wright-1/peptides/blob/master/download_domain.py) to pull all of the data from NCBI.

I started pulling the bacteria domain interactively and it's taking forever, so I'm submitting a batch job for the rest.
- job name: kraken2db_datapull (accidentally have it named the same as medaka_align because i never changed it in batch job oops)
- job id: 48321063
- run time: 00:43:01
```
module load python/3.8.1
conda install -c conda-forge biopython
conda install pandas
python /work/gatins/hci_genome/kraken2/download_domain.py --domain archaea --complete True --ext dna
python /work/gatins/hci_genome/kraken2/download_domain.py --domain viral --complete True --ext dna
python /work/gatins/hci_genome/kraken2/download_domain.py --domain plasmid --complete True --ext dna
python /work/gatins/hci_genome/kraken2/download_domain.py --domain vertebrate_mammalian --complete True --ext dna --human True
```

I ended up needing to pull the bacteria with a batch job too:
- job name: bacterial_data_kraken2
- job id: 48322385
- run time: 06:42:49
```
module load python/3.8.1
module load anaconda3/2022.05
source activate /work/gatins/hci_genome/env
cd /work/gatins/hci_genome/processing/kraken2_builtpython
python /work/gatins/hci_genome/kraken2/download_domain.py --domain bacteria --complete True --ext dna
```

Now that everything is downloaded into a /genomes directory in kraken2_builtpython, I need to --add-to-library:
- job name: kraken_addtolib
- job id: 48364282
- run time: 01:50:18
```
find /work/gatins/hci_genome/processing/kraken2_builtpython/genomes/ -name '*.fna' -print0 | xargs -0 -I{} -n1 -P10 /work/gatins/hci_genome/kraken2/kraken2-build --add-to-library {} --db /work/gatins/hci_genome/processing/kraken2_builtpython
```
Next, we build the Kraken2 database:
- job name: kraken_build
- job id: 48392591
- run time: 1-00:51:47
```
module load gcc/9.2.0
/work/gatins/hci_genome/kraken2/kraken2-build --db /work/gatins/hci_genome/processing/kraken2_builtpython --build --threads 6
```

You'll know the Kraken database has been built properly when you have these three files: 
- hash.k2d: Contains the minimizer to taxon mappings
- opts.k2d: Contains information about the options used to build the database
- taxo.k2d: Contains taxonomy information used to build the database

We have them! So let's go ahead and run the Hifiasm assembly against this database:
```
module load gcc/9.2.0
/work/gatins/hci_genome/kraken2/kraken2 --threads 20 --db /work/gatins/hci_genome/processing/kraken2_builtpython --use-names --report kraken2_hifiasm_nomito_report /work/gatins/hci_genome/processing/mtdna/removal/hifiasm_nomito/assembly_hifiasm_no_mito.fa
```

So, I initially built this **kraken2_builtpython** database but I am also interested in testing this analysis with a "positive control" approach, meaning adding in other fish sequences to the database and tossing all contaminants that are not of fish origin.

To do this, I copied over the library I created my first Kraken database with and am now adding three fish genomes to have positive controls: Atlantic cod (*Gadus morhua*), bicolor damselfish (*Stegastes partitus*), and pinfish (*Lagodon rhomboides*). 
```
cd /work/gatins/hci_genome/processing/krakendb_fish
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/902/167/405/GCF_902167405.1_gadMor3.0/GCF_902167405.1_gadMor3.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/690/725/GCF_000690725.1_Stegastes_partitus-1.0.2/GCF_000690725.1_Stegastes_partitus-1.0.2_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/039/737/535/GCA_039737535.1_Lrho_1.0/GCA_039737535.1_Lrho_1.0_genomic.fna.gz
```
Now, let's add them to our new library and build the database
- job name:
- job id: 48573973
- run time:
```
# Pinfish
/work/gatins/hci_genome/kraken2/kraken2-build --add-to-library GCA_039737535.1_Lrho_1.0_genomic.fna.gz --db /work/gatins/hci_genome/processing/krakendb_fish
# Atlantic cod
/work/gatins/hci_genome/kraken2/kraken2-build --add-to-library GCF_902167405.1_gadMor3.0_genomic.fna.gz --db /work/gatins/hci_genome/processing/krakendb_fish
# Bicolor damselfish
/work/gatins/hci_genome/kraken2/kraken2-build --add-to-library GCF_000690725.1_Stegastes_partitus-1.0.2_genomic.fna.gz --db /work/gatins/hci_genome/processing/krakendb_fish
/work/gatins/hci_genome/kraken2/kraken2-build --db /work/gatins/hci_genome/processing/krakendb_fish --build --threads 10
```
