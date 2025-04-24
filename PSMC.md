# PSMC Analysis: Inferring Population Size History using the Pairwise Sequentially Markovian Coalescent (PSMC) Model

## Download PSMC
```
git clone https://github.com/lh3/psmc.git
make; (cd utils; make)
```
## Generate the whole-genome diploid consensus sequence for input
```
samtools mpileup -C50 -uf ref.fa aln.bam | bcftools view -c - | vcfutils.pl vcf2fq -d 10 -D 100 | gzip > diploid.fq.gz
```
## Run PSMC
```
utils/fq2psmcfa -q20 diploid.fq.gz > diploid.psmcfa
psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o diploid.psmc diploid.psmcfa
utils/psmc2history.pl diploid.psmc | utils/history2ms.pl > ms-cmd.sh
utils/psmc_plot.pl diploid diploid.psmc
```
