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
module load htslib
conda activate atg

gtex_code=/home-4/xwang145@jhu.edu/code/GTExV6PRareVariation/RIVER
ls $gtex_code
gtex_vcf=/scratch/groups/abattle4/jessica/ATG/Final/pseudofun/GTEx_v8_WGS_ARA.recode_SNP.vcf.gz
tabix $gtex_vcf

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

jb_out=~/work/jessica/ATG/Final/reference

#~/.conda/envs/genomics/bin/python $gtex_code/filter_for_rare_variants.py \
head -n1 $RAREVARDIR/download/1KG/AFR_1KG_chr${j}_af.frq > ${jb_out}/GTEx_AFA_AF_rare.frq
for j in {1..22}; do
head -n1 $gtex_vcf_frq > ${jb_out}/GTEx_tmp.frq
awk -v chrom=$j '$1=="chr"chrom {print $0}' $gtex_vcf_frq >> ${jb_out}/GTEx_tmp.frq
python filter_for_rare_variants.py \
${jb_out}/GTEx_tmp.frq \
$RAREVARDIR/download/1KG/AFR_1KG_chr${j}_af.frq \
${jb_out}/GTEx_AFA_AF_rare_chr${j}.frq

cat ~/work/jessica/ATG/Final/reference/GTEx_AFA_AF_rare_chr${j}.frq >> ~/work/jessica/ATG/Final/reference/GTEx_AFA_AF_rare.frq
done


infile=~/work/jessica/ATG/Final/reference/GTEx_AFA_AF_rare.frq
tail -n +2 ${infile} | awk 'BEGIN {OFS="\t"} {split($NF,last,":");  if((last[2] < 0.5)) {print $1,$2-1,$2}}' > ~/work/jessica/ATG/Final/reference/out_file1.bed
tail -n +2 ${infile} | awk 'BEGIN {OFS="\t"} {split($NF,last,":");  if((last[2] > 0.5)) {print $1,$2-1,$2}}' > ~/work/jessica/ATG/Final/reference/out_file2.bed

bcftools view --output-file ~/work/jessica/ATG/Final/reference/gtex_vcf_1.vcf.gz --output-type z --targets-file ~/work/jessica/ATG/Final/reference/out_file1.bed  $gtex_vcf
tabix ~/work/jessica/ATG/Final/reference/gtex_vcf_1.vcf.gz
bcftools view --output-file ~/work/jessica/ATG/Final/reference/gtex_vcf_2.vcf.gz --output-type z --targets-file ~/work/jessica/ATG/Final/reference/out_file2.bed  $gtex_vcf
tabix ~/work/jessica/ATG/Final/reference/gtex_vcf_2.vcf.gz

bcftools query -l ${jb_out}/gtex_vcf_1.vcf.gz | tr '\n' ' '   | awk '{print "CHRPOS",$0}' > ${jb_out}/header.txt


bcftools query -f'%CHROM:%POS [%GT ]\n' ${jb_out}/gtex_vcf_1.vcf.gz | awk '{for(i=2;i<=NF;i++){split($i,gt,"|");out=0;tmp=gt[1]+gt[2]; if (tmp > 0){out=1};$i=out} print $0}' > ${jb_out}/gtex_vcf_1_count.txt

bcftools query -f'%CHROM:%POS [%GT ]\n' ${jb_out}/gtex_vcf_2.vcf.gz | awk '{for(i=2;i<=NF;i++){split($i,gt,"|");out=0;tmp=gt[1]+gt[2]; if (tmp < 2){out=1};$i=out} print $0}'  > ${jb_out}/gtex_vcf_2_count.txt

cat ${jb_out}/header.txt ${jb_out}/gtex_vcf_1_count.txt ${jb_out}/gtex_vcf_2_count.txt > ${jb_out}/gtex_rare_alleles.txt


cut -d " " -f 2- ${jb_out}/gtex_rare_alleles.txt | awk 'BEGIN{OFS=" "}
     NR==1{print $0}
     NR>1{for (i=1;i<=NF;i++) a[i]+=$i}
     END{for (i=1;i<=NF;i++) printf a[i] OFS; printf "\n"}'  > ${jb_out}/gtex_rare_counts.txt




