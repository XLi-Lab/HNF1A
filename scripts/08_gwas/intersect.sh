#!/bin/bash
#SBATCH --job-name=hnf1a_intersect
#SBATCH --output=intersect_%j.log
#SBATCH --error=intersect_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00

# Load bedtools
module load bedtools2/2.31.0-gcc-12.3.0-python-3.11.6

cd /scratch/users/k21049150/HNF1A

# Sort GWAS file using multiple cores
echo "Sorting GWAS file..."
sort -k1,1 -k2,2n --parallel=4 gwas_credible_set.bed > gwas_credible_set_sorted.bed

# Sort ChIP-seq files
echo "Sorting ChIP-seq files..."
sort -k1,1 -k2,2n GSM6248576_Islets_HNF1A_ab96777_peaks.bed > Islets_sorted.bed
sort -k1,1 -k2,2n GSM6248575_D20_EP_HNF1A_ab96777_peaks.bed > D20_EP_sorted.bed
sort -k1,1 -k2,2n GSM6248568_EndoC-betaH1_HNF1A_sc-6547_peaks.bed > EndoC_sorted.bed
sort -k1,1 -k2,2n GSM6248577_HepG2_HNF1A_ab96777_peaks.bed > HepG2_sorted.bed

# intersect each tissue
echo "Intersecting..."
bedtools intersect -a gwas_credible_set_sorted.bed -b Islets_sorted.bed -wa -wb > hits_Islets.bed
bedtools intersect -a gwas_credible_set_sorted.bed -b D20_EP_sorted.bed -wa -wb > hits_D20_EP.bed
bedtools intersect -a gwas_credible_set_sorted.bed -b EndoC_sorted.bed -wa -wb > hits_EndoC.bed
bedtools intersect -a gwas_credible_set_sorted.bed -b HepG2_sorted.bed -wa -wb > hits_HepG2.bed

# count hits
echo "Islets:" && wc -l hits_Islets.bed
echo "D20_EP:" && wc -l hits_D20_EP.bed
echo "EndoC:" && wc -l hits_EndoC.bed
echo "HepG2:" && wc -l hits_HepG2.bed

echo "Done!"
