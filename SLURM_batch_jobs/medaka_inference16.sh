#!/bin/bash
#SBATCH -J medaka_inference16                        # Job name
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
medaka inference /work/gatins/hci_genome/processing/medaka/medaka_align.bam.bam /work/gatins/hci_genome/processing/medaka/results_inference/contigs1203-1344.hdf --model r1041_e82_400bps_sup_v5.0.0 --threads 2 --region contig_1203 contig_1204 contig_1206 contig_1208 contig_1236 contig_1238 contig_1239 contig_1240 contig_1244 contig_1248 contig_1249 contig_1250 contig_1258 contig_1261 contig_1264 contig_1281 contig_1282 contig_1283 contig_1285 contig_1286 contig_1288 contig_1289 contig_1290 contig_1292 contig_1300 contig_1304 contig_1305 contig_1306 contig_1311 contig_1312 contig_1313 contig_1314 contig_1321 contig_1323 contig_1325 contig_1327 contig_1335 contig_1339 contig_1343 contig_1344
