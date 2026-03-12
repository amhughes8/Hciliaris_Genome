# Annotation of the Queen Angelfish (*Holacanthus ciliaris*) genome assembly
## Previously, I annotated the genome with 102 contigs. Since then, I have filtered the assembly to have 29 contigs and have done a test annotation of this filtered assembly to see if any gene content was removed. It doesn't seem like this is the case, so I'm now going to document all of the steps I took to annotate this FINAL assembly!

### 1. RepeatModeler to model repetitive elements in the genome
```
# RepeatModeler build database
cd /projects/gatins/2025_HCI_Genome/processing/final_filtering_blobtools
apptainer exec dfam-tetools-latest.sif BuildDatabase -name hci29_genome_repeats final_assembly_filtered.fasta

mv hci29_genome_repeats* /projects/gatins/2025_HCI_Genome/annotation/hci_29

# RepeatModeler run
cd /projects/gatins/2025_HCI_Genome/annotation/hci_29
apptainer exec dfam-tetools-latest.sif RepeatModeler -LTRStruct -database hci29_genome_repeats -threads 50
```

### 2. RepeatMasker to mask repetitive elements in the genome
```
module load anaconda3/2024.06
module load perl/5.40.0
source activate /projects/gatins/programs_explorer/RepeatMasker
/projects/gatins/2025_HCI_Genome/annotation/RepeatMasker/RepeatMasker -pa 10 -lib hci29_genome_repeats-families.fa -xsmall -gff /projects/gatins/2025_HCI_Genome/processing/final_filtering_blobtools/final_assembly_filtered.fasta
```
```
[hughes.annab@explorer-01 final_filtering_blobtools]$ cat final_assembly_filtered.fasta.tbl
==================================================
file name: final_assembly_filtered.fasta
sequences:            29
total length:  602025945 bp  (602025945 bp excl N/X-runs)
GC level:         41.37 %
bases masked:  112424191 bp ( 18.67 %)
==================================================
               number of      length   percentage
               elements*    occupied  of sequence
--------------------------------------------------
Retroelements        94545     19408870 bp    3.22 %
   SINEs:            22051      3332184 bp    0.55 %
   Penelope:          3310       523186 bp    0.09 %
   LINEs:            54979     10885284 bp    1.81 %
    CRE/SLACS            0            0 bp    0.00 %
     L2/CR1/Rex      43711      8101081 bp    1.35 %
     R1/LOA/Jockey    1665       384993 bp    0.06 %
     R2/R4/NeSL        903       214058 bp    0.04 %
     RTE/Bov-B        4364       770613 bp    0.13 %
     L1/CIN4          1638       924971 bp    0.15 %
   LTR elements:     14205      4668216 bp    0.78 %
     BEL/Pao          1221      1117108 bp    0.19 %
     Ty1/Copia           0            0 bp    0.00 %
     Gypsy/DIRS1      3237      1663161 bp    0.28 %
       Retroviral     5210       969380 bp    0.16 %

DNA transposons     165413     28506155 bp    4.74 %
   hobo-Activator    56266     12179658 bp    2.02 %
   Tc1-IS630-Pogo    40196      7187500 bp    1.19 %
   En-Spm                0            0 bp    0.00 %
   MULE-MuDR           353        65196 bp    0.01 %
   PiggyBac            928       145946 bp    0.02 %
   Tourist/Harbinger  8784      1589107 bp    0.26 %
   Other (Mirage,     2437       428659 bp    0.07 %
    P-element, Transib)

Rolling-circles        801       189462 bp    0.03 %

Unclassified:       405405     49316507 bp    8.19 %

Total interspersed repeats:    97231532 bp   16.15 %


Small RNA:           16293      2699075 bp    0.45 %

Satellites:            306       169943 bp    0.03 %
Simple repeats:     286760     12598240 bp    2.09 %
Low complexity:      35590      1996338 bp    0.33 %
==================================================

* most repeats fragmented by insertions or deletions
  have been counted as one element


RepeatMasker version 4.1.9 , default mode
run with rmblastn version 2.14.1+
The query was compared to classified sequences in "hci29_genome_repeats-families.fa"
FamDB:
```

### 3. Gene prediction with BRAKER3
```
#!/bin/bash
#SBATCH -J hci29_braker_nobusco_rnaseq                    # Job name
#SBATCH -p short                             # Partition
#SBATCH -N 1                                # Number of nodes
#SBATCH -n 30                               # Number of tasks/threads
#SBATCH -o output_%j.txt                    # Standard output file
#SBATCH -e error_%j.txt                     # Standard error file
#SBATCH --mail-user=hughes.annab@northeastern.edu  # Email
#SBATCH --mail-type=END                     # Email notification at job completion
#SBATCH --time=48:00:00                     # Maximum run time

cd /projects/gatins/2025_HCI_Genome/processing/final_filtering_blobtools
cp final_assembly_filtered.fasta.masked /projects/gatins/2025_HCI_Genome/annotation/braker/

cd /projects/gatins/2025_HCI_Genome/annotation/braker
mv final_assembly_filtered.fasta.masked final_assembly_filtered_hci29.fasta.masked
apptainer exec -B /projects/gatins/2025_HCI_Genome/annotation/braker /projects/gatins/2025_HCI_Genome/annotation/braker/braker3.sif braker.pl \
--genome=/projects/gatins/2025_HCI_Genome/annotation/braker/final_assembly_filtered_hci29.fasta.masked \
--prot_seq=/projects/gatins/2025_HCI_Genome/annotation/braker/Vertebrata.fa \
--bam=/projects/gatins/2025_HCI_Genome/annotation/braker/hci29_mapped_rnaseq.bam \
--threads=30 --species=Hciliaris29 --softmasking \
--AUGUSTUS_CONFIG_PATH=/projects/gatins/2025_HCI_Genome/annotation/braker/config &> hci29_nobusco_rnaseq_braker.log
```

[hughes.annab@explorer-02 hci29_braker]$ pwd
/projects/gatins/2025_HCI_Genome/annotation/braker/hci29_braker

```
grep -c "^>" braker.aa
```
27142

### 4. Transcript filtering with TSEBRA

Filter out single-exon genes:
```
apptainer exec braker3.sif tsebra.py \
-g /projects/gatins/2025_HCI_Genome/annotation/braker/hci29_braker/GeneMark-ETP/genemark.gtf \
-k /projects/gatins/2025_HCI_Genome/annotation/braker/hci29_braker/Augustus/augustus.hints.gtf \
-e /projects/gatins/2025_HCI_Genome/annotation/braker/hci29_braker/hintsfile.gff \
--filter_single_exon_genes \
-c /projects/gatins/2025_HCI_Genome/annotation/braker/default.cfg \
-o hci29_braker_nseg.gtf
```

Filter for longest isoform:
```
apptainer exec braker3.sif get_longest_isoform.py --gtf hci29_braker_nseg.gtf --out hci29_braker_nseg_li.gtf
```

### 5. BUSCO (protein mode)
Extract protein sequences:
```
/projects/gatins/programs_explorer/gffread/bin/gffread -w hci29_braker_nseg_li.fa -y hci29_braker_nseg_li.aa -g /projects/gatins/2025_HCI_Genome/annotation/braker/final_assembly_filtered_hci29.fasta.masked hci29_braker_nseg_li.gtf
```

```
module load anaconda3/2024.06
source activate /projects/gatins/programs_explorer/busco
busco -i hci29_braker_nseg_li.aa --mode proteins --lineage_dataset actinopterygii_odb12 --cpu 10 --out hci29_braker_busco
```

    ---------------------------------------------------
    |Results from dataset actinopterygii_odb12         |
    ---------------------------------------------------
    |C:99.1%[S:98.8%,D:0.4%],F:0.6%,M:0.3%,n:7207      |
    |7145    Complete BUSCOs (C)                       |
    |7118    Complete and single-copy BUSCOs (S)       |
    |27    Complete and duplicated BUSCOs (D)          |
    |43    Fragmented BUSCOs (F)                       |
    |19    Missing BUSCOs (M)                          |
    |7207    Total BUSCO groups searched               |
    ---------------------------------------------------

```
grep -c ">" hci29_braker_nseg_li.aa
```
29049
