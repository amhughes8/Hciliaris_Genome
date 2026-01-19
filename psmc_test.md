Align high quality dataset to ref genome
```
/projects/gatins/programs_explorer/minimap2/minimap2 -t 30 -ax map-ont /projects/gatins/2025_HCI_Genome/processing/final_filtering_blobtools/final_assembly_filtered.fasta /projects/gatins/2025_HCI_Genome/processing/hci_filtered_3kQ10.fastq > HCI_new_aligned.sam
```

```
module load samtools/1.21
samtools view -Sb -@ 30 -o HCI_new_aligned.bam HCI_new_aligned.sam
samtools sort -o HCI_new_aligned_sorted.bam -O bam -@ 20 HCI_new_aligned.bam
samtools index -b -@ 20 HCI_new_aligned_sorted.bam
```

Generate whole-genome consensus file for PSMC and gzip
```
samtools consensus --ambig -f fastq -d 50 HCI_new_aligned_sorted.bam -o consensus_new.fq
gzip consensus_new.fq
```

Run PSMC
```
/projects/gatins/2025_HCI_Genome/PSMC/psmc/utils/fq2psmcfa -q20 consensus_new.fq.gz > diploid_HCI_new.psmcfa
/projects/gatins/2025_HCI_Genome/PSMC/psmc/psmc -N30 -t30 -r5 -p "4+30*2+4+6+10" -o diploid_HCI_new_final.psmc diploid_HCI_new.psmcfa
```

Plot
```
source activate /projects/gatins/programs_explorer/gnuplot
/projects/gatins/2025_HCI_Genome/PSMC/psmc/utils/psmc_plot.pl -u 5.97e-09 -g 5 HCI_new_t30r5_plot_u597-9g5 diploid_HCI_new_final.psmc
```

Bootstrapping
```
/projects/gatins/2025_HCI_Genome/PSMC/psmc/utils/splitfa diploid_HBE.psmcfa > diploid_HBE_split.psmcfa
```
```
mkdir bootstrap
cp diploid_HBE_split.psmcfa bootstrap
cp diploid_HBE_final.psmc bootstrap
cd bootstrap
echo split_HBE_{001..100}.psmcfa| xargs -n 1 cp diploid_HBE_split.psmcfa
```

```
#!/bin/bash
#SBATCH -J psmc_array			    # Job name
#SBATCH -p short                            # Partition
#SBATCH -N 1                                # Number of nodes
#SBATCH -n 2                                # Number of tasks/threads
#SBATCH -o array_%A_%a.out    		    # Name of stdout output file
#SBATCH -e array_%A_%a.err    		    # Name of stdout output file
#SBATCH --array=1-100			    # Array index
#SBATCH --mem=6000MB 			    # Memory to be allocated PER NODE
#SBATCH --mail-user=hughes.annab@northeastern.edu  # Email
#SBATCH --mail-type=END                     # Email notification at job completion
#SBATCH --time=48:00:00                     # Maximum run time

echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
#
# ----------------Your Commands------------------- #
#
echo "This job in the array has:"
echo "- SLURM_JOB_ID=${SLURM_JOB_ID}"
echo "- SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"

# select our filename
N=${SLURM_ARRAY_TASK_ID}
# Comment one of the following two lines, depending on if the file names have leading zeros
#FILENAME=run-${N}.inp # without leading zeros
 FILENAME=split_HCI_$(printf "%03d" ${N}).psmcfa # with leading zeros
# adjust "%03d" to as many digits as are in the numeric part of the file name
echo "My input file is ${FILENAME}"

#
echo $P
#
/projects/gatins/2025_HCI_Genome/PSMC/psmc/psmc -N30 -t30 -r5 -b -p "4+30*2+4+6+10" -o /projects/gatins/2025_HBE_Genome/PSMC/bootstrap/${FILENAME}.psmc /projects/gatins/2025_HBE_Genome/PSMC/bootstrap/${FILENAME}
#

echo "Job finished" `date`
echo "My input file is ${FILENAME}"
```

```
cat *.psmc > HBE_combined.psmc
```

```
source activate /projects/gatins/programs_explorer/gnuplot
/projects/gatins/2025_HCI_Genome/PSMC/psmc/utils/psmc_plot.pl -u 5.97e-09 -g 5 HBE_t30r5_plot_u597-9g5 HBE_combined.psmc
```
