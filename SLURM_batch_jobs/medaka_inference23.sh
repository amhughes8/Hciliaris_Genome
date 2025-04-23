#!/bin/bash
#SBATCH -J medaka_inference23                        # Job name
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
medaka inference /work/gatins/hci_genome/processing/medaka/medaka_align.bam.bam /work/gatins/hci_genome/processing/medaka/results_inference/contigs2239-2370.hdf --model r1041_e82_400bps_sup_v5.0.0 --threads 2 --region contig_2239 contig_2242 contig_2246 contig_2252 contig_2258 contig_2263 contig_2264 contig_2267 contig_2269 contig_2271 contig_2272 contig_2273 contig_2274 contig_2277 contig_2282 contig_2283 contig_2284 contig_2285 contig_2288 contig_2289 contig_2291 contig_2292 contig_2299 contig_2300 contig_2301 contig_2302 contig_2303 contig_2304 contig_2305 contig_2308 contig_2311 contig_2314 contig_2315 contig_2318 contig_2320 contig_2321 contig_2322 contig_2324 contig_2327 contig_2329 contig_2330 contig_2331 contig_2334 contig_2338 contig_2339 contig_2341 contig_2345 contig_2347 contig_2348 contig_2349 contig_2355 contig_2361 contig_2364 contig_2370
