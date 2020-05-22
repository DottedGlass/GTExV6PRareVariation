#!/bin/bash

dir=$RAREVARDIR

cat ${dir}/archive/data/v8/outliers_medz_picked.txt | awk '$3==5' > possible.txt

# subset the normalized expression file to these example genes
expr=${dir}/preprocessing/gtex_v8_normalized_expression.txt
head -n1 $expr > expr.subset.by.genes.txt
cat possible.txt | cut -f1 | grep -f - $expr >> expr.subset.by.genes.txt


