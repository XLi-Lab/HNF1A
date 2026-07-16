#!/bin/bash -l
#SBATCH --job-name=hnf1a_closest
#SBATCH --output=closest_%j.log
#SBATCH --error=closest_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00

module load bedtools2/2.31.0-gcc-12.3.0-python-3.11.6
cd /scratch/users/k21049150/HNF1A

# sort the gene BED (closest requires both inputs coordinate-sorted)
sort -k1,1 -k2,2n gencode_protein_coding.bed > gencode_protein_coding_sorted.bed

# sort the variant BEDs (closest requires -a sorted too)
sort -k1,1 -k2,2n top_hits_pp05_dedup.bed > top_hits_pp05_dedup_sorted.bed
sort -k1,1 -k2,2n hep_top_hits_pp05_dedup.bed > hep_top_hits_pp05_dedup_sorted.bed

# islet
bedtools closest -a top_hits_pp05_dedup_sorted.bed \
                 -b gencode_protein_coding_sorted.bed \
                 -d > nearest_gene_Islets.bed

# hepatocyte
bedtools closest -a hep_top_hits_pp05_dedup_sorted.bed \
                 -b gencode_protein_coding_sorted.bed \
                 -d > nearest_gene_HepG2.bed

echo "Done"
