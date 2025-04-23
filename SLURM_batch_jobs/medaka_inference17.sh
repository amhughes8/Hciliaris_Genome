#!/bin/bash
#SBATCH -J medaka_inference17                        # Job name
#SBATCH -p short                            # Partition
#SBATCH -N 1                                # Number of nodes
#SBATCH -n 2                                # Number of tasks/threads
#SBATCH -o output_%j.txt                    # Standard output file
#SBATCH -e error_%j.txt                     # Standard error file
#SBATCH --mail-user=hughes.annab@northeastern.edu  # Email
#SBATCH --mail-type=END                     # Email notification at job completion
#SBATCH --time=48:00:00                     # Maximum run time


# Your program/command here
module load anaconda3/2022.05
source activate medaka
medaka inference /work/gatins/hci_genome/processing/medaka/medaka_align.bam.bam /work/gatins/hci_genome/processing/medaka/results_inference/contigs1373-1483.hdf --model r1041_e82_400bps_sup_v5.0.0 --threads 2 --region contig_1373 contig_1375 contig_1376 contig_1388 contig_1390 contig_1391 contig_1392 contig_1393 contig_1397 contig_1398 contig_1400 contig_1402 contig_1405 contig_1406 contig_1410 contig_1411 contig_1414 contig_1415 contig_1416 contig_1422 contig_1423 contig_1424 contig_1425 contig_1428 contig_1435 contig_1437 contig_1439 contig_1446 contig_1450 contig_1458 contig_1460 contig_1462 contig_1466 contig_1472 contig_1473 contig_1477 contig_1479 contig_1480 contig_1482 contig_1483
