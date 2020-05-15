#!/bin/bash

infolder=$RAREVARDIR/preprocessing/PEER_gtex8
outfolder=$RAREVARDIR/preprocessing/PEER_gtex8_wide
# outfolder=/scratch/groups/abattle4/jessica/ATG/Final/PEER_gtex8_wide
for file in $(ls ${infolder}/*globalOutliers.v8ciseqtls.removed.outlierGenes.zscores.txt); do
base=$(basename $file)
echo ${file}
Rscript --vanilla munge.R ${file} ${outfolder}/${base}
done