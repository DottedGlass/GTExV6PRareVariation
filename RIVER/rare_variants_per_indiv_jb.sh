#!/bin/bash
#SBATCH --job-name=rare_var_per_indiv
#SBATCH --time=8:0:0
#SBATCH --partition=lrgmem
#SBATCH --nodes=1
# number of tasks (processes) per node
#SBATCH --ntasks-per-node=1
#SBATCH --mail-type=end
#SBATCH --mail-user=xwang145@jhu.edu


module load vcftools
module load bcftools
conda activate atg

gtex_code=/home-4/xwang145@jhu.edu/code/GTExV6PRareVariation/RIVER
ls $gtex_code
gtex_vcf=/scratch/groups/abattle4/jessica/ATG/Final/pseudofun/GTEx_v8_WGS_ARA.recode_SNP.vcf.gz
afr_1kg=/scratch/groups/abattle4/jessica/ATG/Final/1KG_AFR_biallelic
RAREVARDIR=/scratch/groups/abattle4/victor/GTExV6PRareVariationData
# # make 1KG frequency file
# for j in {1..22};do
# vcftools --gzvcf ${afr_1kg}/AFR_1KG_chr${j}.vcf.gz --freq --out $RAREVARDIR/download/1KG/AFR_1KG_chr${j}_af
# done
#
oneKG_vcf_frq=$RAREVARDIR/reference/AFR_1KG_chrALL_af.frq
# head -n1 $RAREVARDIR/download/1KG/AFR_1KG_chr1_af.frq > $oneKG_vcf_frq
# for j in {1..22}; do
# tail -n +2 $RAREVARDIR/download/1KG/AFR_1KG_chr${j}_af.frq >> $oneKG_vcf_frq
# done

# make gtex frequency file
gtex_vcf_frq=$RAREVARDIR/reference/GTEx_AFA_AF.frq
# vcftools --gzvcf ${gtex_vcf} --freq --out $RAREVARDIR/reference/GTEx_AFA_AF

#~/.conda/envs/genomics/bin/python $gtex_code/filter_for_rare_variants.py \
head -n1 $RAREVARDIR/download/1KG/AFR_1KG_chr${j}_af.frq > ~/work/jessica/ATG/Final/reference/GTEx_AFA_AF_rare.frq
for j in {1..22}; do
head -n1 $gtex_vcf_frq > ~/work/jessica/ATG/Final/reference/GTEx_tmp.frq
awk -v chrom=$j '$1=="chr"chrom {print $0}' $gtex_vcf_frq >> ~/work/jessica/ATG/Final/reference/GTEx_tmp.frq
python filter_for_rare_variants.py \
~/work/jessica/ATG/Final/reference/GTEx_tmp.frq \
$RAREVARDIR/download/1KG/AFR_1KG_chr${j}_af.frq \
~/work/jessica/ATG/Final/reference/GTEx_AFA_AF_rare_chr${j}.frq

cat ~/work/jessica/ATG/Final/reference/GTEx_AFA_AF_rare_chr${j}.frq >> ~/work/jessica/ATG/Final/reference/GTEx_AFA_AF_rare.frq
done


infile=~/work/jessica/ATG/Final/reference/GTEx_AFA_AF_rare.frq
tail -n +2 ${infile} | awk 'BEGIN {OFS="\t"} {split($NF,last,":");  if((last[2] < 0.5)) {print $1,$2-1,$2}}' > ~/work/jessica/ATG/Final/reference/out_file1.bed
tail -n +2 ${infile} | awk 'BEGIN {OFS="\t"} {split($NF,last,":");  if((last[2] > 0.5)) {print $1,$2-1,$2}}' > ~/work/jessica/ATG/Final/reference/out_file2.bed

bcftools --output-file ~/work/jessica/ATG/Final/reference/gtex_vcf_1.vcf.gz --output-type z --regions-file $RAREVARDIR/reference/out_file1.bed --no-header $gtex_vcf
bcftools --output-file ~/work/jessica/ATG/Final/reference/gtex_vcf_2.vcf.gz --output-type z --regions-file $RAREVARDIR/reference/out_file2.bed --no-header $gtex_vcf
