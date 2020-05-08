# Dependencies
python2.7
* numpy
* pybedtools
* pysam
* scipy

vcftools

# Pipeline

## Preparing reference files used for RIVER

### Get list of subjects from v8
```
cat ${RAREVARDIR}/download/gtex8/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt | awk -F "\t" 'NR>1{print $1}' \
> ${RAREVARDIR}/preprocessing/gtex_v8_individuals_all_normalized_samples.txt
```

### Get list of African American subjects
```
cat ${RAREVARDIR}/download/gtex8/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt | awk -F '\t' '{if ($5 == 2) print $1;}' \
> ${RAREVARDIR}/preprocessing/gtex_v8_individuals_AFA.txt
```

### Get list of African American subjects with available WGS data
```
vcf-query -l ${RAREVARDIR}/download/gtex8/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz \
> ${RAREVARDIR}/preprocessing/gtex_v8_wgs_individuals.txt

awk 'NR==FNR { lines[$0]=1; next } $0 in lines' ${RAREVARDIR}/preprocessing/gtex_v8_wgs_individuals.txt ${RAREVARDIR}/preprocessing/gtex_v8_individuals_AFA.txt > ${RAREVARDIR}/preprocessing/gtex_v8_wgs_individuals_AFA.txt
```

### Get list of tissue names
```
ls ${RAREVARDIR}/preprocessing/PEER_gtex8 | cut -d. -f1 \
> ${RAREVARDIR}/preprocessing/gtex_v8_tissues_all_normalized_samples.txt
```

### Make flat files from normalized data
Flat file is `preprocessing/gtex_v8_normalized_expression.txt`
```
cd preprocessing
python gather_filter_normalized_expression_v8.py
```

### gene locations from gencode
`gencode.v26.genes.v8.patched_contigs_genetypes_autosomal.txt` is empty for some reason
```
cd preprocessing

bash process.reference.files_v8.sh ${RAREVARDIR}/download/gtex8/gencode.v26.GRCh38.genes.gtf.gz ${RAREVARDIR}/download/gencode/gencode.v26.annotation.gtf
```

## Outlier calling from normalized expression data
### Call multi-tissue outliers
```
Rscript call_outliers/call_outliers_medz_v8.R
```
Note that this version outputs `preprocessing/gtex_v8_wgs_ids_outlier_filtered.txt`, which contains individuals of multiple ancestries, not just European ancestry.

### Call single-tissue outliers
```
python call_outliers/call_outliers_single_tissue_v8.py
```

### Filter out African American outliers


## Making rare variant calls on African American individuals


### [debugging] Get TSS of relevant genes
Get list of relevant genes (lincRNA and protein coding only)
```
cut -f1 ${RAREVARDIR}/reference/gencode.v26.GRCh38.genes_genetypes_autosomal_PCandlinc_only.txt > ${RAREVARDIR}/reference/v8.genes.lincRNA.protein.txt
```
Use that list to subset `gencode.v26.genes.v8.patched_contigs_TSS.bed` and get TSS. Output file is called `reference/v8.genes.TSS.bed`
```
python preprocessing/get_TSS.py
```

### [debugging] Check for rare variants
#### Compute rare variants only from 103 African American individuals with WGS
Get subset of wgs with only African American individuals
```
bcftools view -S ${RAREVARDIR}/preprocessing/gtex_v8_wgs_individuals_AFA.txt GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz -o ${RAREVARDIR}/download/gtex8/GTEx_v8_WGS_ARA.vcf.gz
```
Compute allele frequency
```
vcftools --gzvcf ${RAREVARDIR}/download/gtex8/GTEx_v8_WGS_ARA.vcf.gz --freq --remove-indels --remove-filtered-all --out ${RAREVARDIR}/download/gtex8/GTEx_v8_WGS_ARA_2
```

### Get 10kb positions for TSS of genes
```
cat ${RAREVARDIR}/reference/gencode.v26.genes.v8.patched_contigs_TSS.bed | awk '{print $1, $2-1000, $3-1, $4}' OFS='\t' \
| grep -v '^chrM' > ${RAREVARDIR}/reference/gencode.v26.genes.v8.10kb_TSS.bed

grep -v '^chrM' ${RAREVARDIR}/reference/gencode.v26.genes.v8.10kb_TSS.bed > ${RAREVARDIR}/reference/gencode.v26.genes.v8.10kb_TSS.bed
```

Use bedtools to get intersection on TSS file with GTEx SNVs
```
bedtools intersect -header -a ${RAREVARDIR}/data/wgs/GTEx_AFA.vcf.gz -b ${RAREVARDIR}/reference/gencode.v26.genes.v8.10kb_TSS.bed > ${RAREVARDIR}/data/wgs/GTEx_AFA_10kb_TSS.vcf
```

## RIVER on African American individuals

#### [Step 0] Prepare  GTEx data files

##### Get list of African American individuals that have WGS in GTExv8
```
grep -f ${RAREVARDIR}/preprocessing/gtex_v8_individuals_AFA.txt ${RAREVARDIR}/preprocessing/gtex_v8_wgs_ids_outlier_filtered.txt \
> ${RAREVARDIR}/preprocessing/gtex_v8_wgs_ids_outlier_filtered_AFA.txt
```

#### Expression Ensembl_ID file
There exists a `gene_ensembl_ids.txt` in the repo, which is a list of genes we are restricting to and was used in v6p analysis. May need to reproduce a v8 specific version.

#### [Step 1] Generate matlab annotation files with both indivs and genes considered in GTEx v6p and gene expression matrix used for calling outliers later
```
cd RIVER

matlab -nodisplay -nodesktop -nosplash -nojvm -singleCompThread -r "run('generate_annotations_matlab_v8.m')"
matlab -nodisplay -nodesktop -nosplash -nojvm -singleCompThread -r "run('${RAREVARDIR}/RIVER/code/generate_expmat_44tissues.m')"
```
#### [Step 4] Extract a list of targeted regions for genes of interest
```
matlab -nodisplay -nodesktop -nosplash -nojvm -singleCompThread -r "run('${RAREVARDIR}/RIVER/code/generate_regions.m')"
```


# Figures and Analysis

## Get genes with African American individuals as multi-tissue outliers
```
grep -f ${RAREVARDIR}/preprocessing/gtex_v8_individuals_AFA.txt ${RAREVARDIR}/data/outliers_medz_picked.txt
```

## Compare number of genes for which a GTEx individual is a multi-tissue outlier

### subset v8 individuals by v6p individuals (441 subjects shared between v6p and v8)
```
grep -f ${RAREVARDIR}/preprocessing/gtex_v8_individuals_all_normalized_samples.txt ${RAREVARDIR}/data/outliers_medz_picked_counts_per_ind.txt \
> ${RAREVARDIR}/figures/outliers_medz_picked_counts_per_ind_v6p.txt

grep -f ${RAREVARDIR}/preprocessing/gtex_2015-01-12_individuals_all_normalized_samples.txt ${RAREVARDIR}/data/v8/outliers_medz_picked_counts_per_ind.txt \
> ${RAREVARDIR}/figures/outliers_medz_picked_counts_per_ind_v8.txt
```


# Pipeline
## Expression data correction and normalization
Generates processed data that can be downloaded from https://s3-us-west-2.amazonaws.com/gtex-v6p-rare-variation-data/GTExV6PRareVariationData.tar.gz. <br>

#### Generate rpkm and read count matrices from GTEx V6P combined file
Using v8 SampleAttributes instead of v6p
```
cat /work-zfs/abattle4/lab_data/GTEx_v8/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.txt | \
	cut -f1,14 | sed 's/ - /_/' | sed 's/ /_/g' | sed 's/(//' | sed 's/)//' | sed 's/c-1/c1/' | grep -v NA12878 > \
	${RAREVARDIR}/preprocessing/gtex_2017-06-05_samples_tissues.txt

```
Note that the following samples are missing tissue names from from samples_tissues.txt
GTEX-145LU-0004-SM-DLIOV
GTEX-14JFF-0002-SM-DLIOQ
GTEX-1JKYN-0004-SM-AYROE
GTEX-1K2DA-0002-SM-AYROF
GTEX-1LGRB-0003-SM-AYROG
GTEX-1LVAN-0004-SM-AYROD

Produced rpkm
```
python3.7 preprocessing/split_by_tissues.py \
    --GTEX ${RAREVARDIR}/download/GTEx_v6p/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct \
    --SAMPLE ${RAREVARDIR}/preprocessing/gtex_2017-06-05_samples_tissues.txt \
    --OUT ${RAREVARDIR}/preprocessing/PEER \
    --END .rpkm.txt
```

```
python preprocessing/split_by_tissues.py \
    --GTEX <path to downloaded file>/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct \
    --SAMPLE ${RAREVARDIR}/preprocessing/gtex_2017-06-05_samples_tissues.txt \
    --OUT ${RAREVARDIR}/preprocessing/PEER \
    --END .reads.txt
```

#### Run bash scripts to generate PEER corrected data (includes non-EAs) with covariates removed
(Uses multiple cores. Currently set to 10 cores. Can be altered by changing the number of cores <br>
specified by `parallel --jobs 10` in the scripts `preprocessing/PEER/calc.PEER.factors.all.tissues.sh` and <br>
`preprocessing/PEER/calc.residuals.sh`) <br>
```
bash preprocessing/PEER/PEER.pipeline.sh
```

#### Make list of individuals and tissues
```
bash preprocessing/get_tissue_by_individual.sh
```

#### Make flat files from raw RPKMs and PEER-corrected data for all tissues and individuals
```
python preprocessing/gather_filter_normalized_expression.py
python preprocessing/gather_filter_rpkm.py
```

#### Make list of expressed genes
```
cat preprocessing/PEER/*peer.ztrans.txt | cut -f1 | sort | uniq | grep -v Id > preprocessing/gtex.expressed.genes.txt
```

#### Make summary statistics for expressed genes
(Uses multiple cores. Can set the number at the top of the script. Currently set to 10 cores.)
```
Rscript preprocessing/rpkm.expression.analysis.R
```

## Preparing reference files used later
Generates processed data that can be downloaded from https://s3-us-west-2.amazonaws.com/gtex-v6p-rare-variation-data/GTExV6PRareVariationData.tar.gz. <br>

Copy or move the GTEx annotation file (`gencode.v19.genes.v6p_model.patched_contigs.gtf.gz`) to `${RAREVARDIR}/reference`.
```
bash preprocessing/process.reference.files.sh <path to>/gencode.v19.genes.v6p_model.patched_contigs.gtf.gz <path to>/gencode.v19.annotation.gtf
```
(relies on `pad.gtf.exons.py`, `gtf2TSS.sh`, and `gtf2genebed.sh`)


## Outlier calling
Generates processed data that can be downloaded from https://s3-us-west-2.amazonaws.com/gtex-v6p-rare-variation-data/GTExV6PRareVariationData.tar.gz. <br>

#### Call multi-tissue outliers
(Uses 12 cores. Can set the number at the top of the script.)
```
Rscript call_outliers/call_outliers_medz.R
```

#### Call single-tissue outliers
```
python call_outliers/call_outliers_single_tissue.py
```

#### Compare single-tissue and multi-tissue outliers as well as get stats on each
```
Rscript call_outliers/compare_single_multi_outliers.R
```

#### Run replication for single-tissue and multi-tissue outliers
(The multi-tissue replication uses 12 cores. Can set the number at the top of the script.)
```
Rscript call_outliers/multi_tissue_replication.R
Rscript call_outliers/single_tissue_replication.R
```

## Feature generation

#### Processing VCFs into bed files for each individual
```
bash feature_construction/vcf2bedfiles.sh \
	 <path to>/GTEx_Analysis_2015-01-12_WholeGenomeSeq_148Indiv_GATK_HaplotypeCaller.vcf.gz \
	 <path to>/gtex.lumpy.gs.svscore.low_conf.vcf.gz
```
(this script and some of its dependencies use multiple cores [set number at top of relevant scripts]; relies on :
* `vcf2bedfiles_helper_processVCF.sh`
* `vcf2bedfiles_helper_processVCF_SV.sh`
* `vcf2bedfiles_helper_processVCFtoolsOutput.sh`
* `vcf2bedfiles_helper_processVCFtoolsOutput_CNV.sh`
* `compileCADDscores.sh`
* `extractCADDscores_ekt.py`)

#### Extract features to be combined with the individual bed files
```
bash feature_construction/extract.1kg.AF.sh
```
(uses 15 cores, the number of which is set at the top of the script; relies on `process.1kg.AF.py`)
```
bash feature_construction/subset.CADD.features.sh
```
(uses 8 cores, set in the sort command)
```
bash feature_construction/TFBS_pipeline.sh
```
(relies on `pouya.raw.summary.py`)
```
bash feature_construction/ER_pipeline.sh
```

#### Add extracted features to individual bed files
```
bash run_add_features_variant_beds.sh
```
**Important:** Make sure to use bedtools version 2.26.0 or later.
Memory leak in previous versions causes the memory for this script to blow up. <br>
(uses 15+ cores. Set the number of processes at the top of the script. Relies on `add_features_variant_beds.sh`.)

#### Collapse site-level features created above into gene-level features
```
bash feature_construction/run_build_feature_count_summaries_all_genes.sh
```
(uses multiple cores; relies on :
* `build_count_summaries_all_genes.sh` set number of processes at top of script
* `build_feature_summaries_all_genes.sh` set number of processes at top of script
* `build_feature_set.py`)

#### Compile features for outliers and controls
```
bash feature_construction/run_compile_features_outliers.sh
```
(uses up to 10 cores; relies on:
* `compile_features_outliers.sh` set number of processes at top of script
* `compile_features_outliers_nothresh.sh`
* `compile_features_outliers_singletissue.sh` set number of processes at top of script
* `pick_outliers_controls_imbalanced.py`)

#### Compute minor allele frequencies of rare variants segregating in GTEx using the UK10K dataset
```
bash feature_construction/calc.uk10k.freqs.sh
```

## Disease gene annotations
We are providing the processed gene lists for the eight disease gene sets we analyzed for overlap with genes with multi-tissue outliers. <br>
We are also providing, where applicable, the commands and raw files needed to generate these processed lists.

#### ACMG
Source: http://www.ncbi.nlm.nih.gov/clinvar/docs/acmg/ <br>
Raw file: `acmg.csv` <br>
Download from https://s3-us-west-2.amazonaws.com/gtex-v6p-rare-variation-data/GTExV6PRareVariationData.tar.gz. <br>
Move the downloaded file to `${RAREVARDIR}/features/annotations/ACMG/`.

#### ClinVar
Source: ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/gene_condition_source_id <br>
Raw file: `gene_condition_source_id` <br>
Download from https://s3-us-west-2.amazonaws.com/gtex-v6p-rare-variation-data/GTExV6PRareVariationData.tar.gz. <br>
Move the downloaded file to `${RAREVARDIR}/features/annotations/ClinVar/`.

#### GWAS
Source: http://www.ebi.ac.uk/gwas/ <br>
Raw file: `gwas_catalog_v1.0-downloaded_2015-11-30.tsv` <br>
Download from https://s3-us-west-2.amazonaws.com/gtex-v6p-rare-variation-data/GTExV6PRareVariationData.tar.gz. <br>
Move the downloaded file to `${RAREVARDIR}/features/annotations/GWAS/`.

#### OMIM
Source: http://www.omim.org/ <br>
Raw files: `morbidmap.txt` and `mim2gene.txt` <br>
Processed file: `omim.genes.txt` <br>
Download the raw and processed files from https://s3-us-west-2.amazonaws.com/gtex-v6p-rare-variation-data/GTExV6PRareVariationData.tar.gz. <br>
Move these files to `${RAREVARDIR}/features/annotations/OMIM/`.

To produce the processed file from the raw files: <br>
```
grep '(3)' ${RAREVARDIR}/features/annotations/OMIM/morbidmap.txt | cut -f2 | sed 's/, /\n/g' | sort | uniq > ${RAREVARDIR}/features/annotations/OMIM/omim.genes.temp.txt

grep -wf ${RAREVARDIR}/features/annotations/OMIM/omim.genes.temp.txt ${RAREVARDIR}/features/annotations/OMIM/mim2gene.txt  > ${RAREVARDIR}/features/annotations/OMIM/temp.mim2gene.intersection.txt

cut -f4,5 ${RAREVARDIR}/features/annotations/OMIM/temp.mim2gene.intersection.txt | sort | uniq > ${RAREVARDIR}/features/annotations/OMIM/omim.genes.txt

rm ${RAREVARDIR}/features/annotations/OMIM/omim.genes.temp.txt ${RAREVARDIR}/features/annotations/OMIM/temp.mim2gene.intersection.txt
```

#### Orphanet
Source: http://www.orphadata.org/data/xml/en_product6.xml <br>
Raw file: `en_product6.xml` <br>
Processed file: `orphanet.genes.txt` <br>
Download the raw and processed files from https://s3-us-west-2.amazonaws.com/gtex-v6p-rare-variation-data/GTExV6PRareVariationData.tar.gz. <br>
Move these files to `${RAREVARDIR}/features/annotations/Orphanet/`.

To produce the processed file from the raw file: <br>
```
grep ENSG ${RAREVARDIR}/features/annotations/Orphanet/en_product6.xml | sort | uniq | grep -o 'ENSG[0-9]*' > ${RAREVARDIR}/features/annotations/Orphanet/orphanet.genes.txt
```

#### DDG2P
Source: http://www.ebi.ac.uk/gene2phenotype/downloads <br>
Raw file: `DDG2P_2_8_2017.csv.gz` <br>
Download the raw file from https://s3-us-west-2.amazonaws.com/gtex-v6p-rare-variation-data/GTExV6PRareVariationData.tar.gz. <br>
Move this file to `${RAREVARDIR}/features/annotations/DDG2P/`.

#### Other: Cardiovascular and Cancer disease genes
We assessed overlap of genes with multi-tissue outliers with two expert curated disease gene lists: one for heritable cancer predisposition and one for heritable cardiovascular disease. See the methods section of our manuscript for more information. <br>
Raw files: `cancer.genes.gold.standard.csv` (Cancer), `cardio.genes.gold.standard.csv` (Cardio)
Download the raw files from https://s3-us-west-2.amazonaws.com/gtex-v6p-rare-variation-data/GTExV6PRareVariationData.tar.gz. <br>
Move these files to `${RAREVARDIR}/features/annotations/Other/`.

## Shared eQTLs defined by METASOFT
Process the METASOFT results choosing the single best variant tested per gene as determined by P-value from the RE2 model. <br>
Also provide summary statistics regarding the number of tissues the eQTL is active in and the expression level for the gene across tissues.
```
python shared.eqtls/bf.metasoft.py --META ${RAREVARDIR}/data/metasoft/Metasoft_Output_v6p.txt \
    --TISS ${RAREVARDIR}/data/metasoft/Metasoft_tissue_order.txt \
    --OUT ${RAREVARDIR}/data/metasoft/gtex.metasoft.v6p.selected.txt
Rscript shared.eqtls/metasoft.summary.R
```

## Validation of large-effect rare variants using CRISPR-Cas9 genome editing
Prioritize variants for validation with CRISPR
```
Rscript crispr/prioritize.for.crispr.R
```

Prune the list of prioritized CRISPR variants down to a manageable size (N ~ 12)
```
Rscript crispr/prune.crispr.variants.R
```

Add major/minor allele info
```
python crispr/add.major.minor.alleles.py --IN ${RAREVARDIR}/data/CRISPR/crispr.overexpression.candidates.pruned.vcf \
    --OUT ${RAREVARDIR}/data/CRISPR/crispr.overexpression.candidates.pruned.maf.alleles.vcf \
    --FRQ ${RAREVARDIR}/features/variantBeds/GTEx_Analysis_2015-01-12_WholeGenomeSeq_148Indiv_GATK_HaplotypeCaller_123EAonly_SNPs.frq
```

Extract 100 bp sequences from the hg19 reference centered on each variant position
```
out=${RAREVARDIR}/data/CRISPR/crispr.candidates.donor.seq.raw.fa
if [ -e $out ]; then
    rm $out
fi
ref=${HG19}/hg19.fa
positions=`tail -n +2 ${RAREVARDIR}/data/CRISPR/crispr.candidates.pruned.vep.loftee.parsed.vcf | awk '{print "chr"$1":"$2-49"-"$2+49}'`
for line in $positions; do
    samtools faidx $ref $line >> $out
done
```

Generate donor sequences. One sequence for wild-type and one for rare allele.
```
python crispr/process.crispr.donor.seq.py \
    --VCF ${RAREVARDIR}/data/CRISPR/crispr.candidates.pruned.vep.loftee.parsed.vcf \
    --FASTA ${RAREVARDIR}/data/CRISPR/crispr.candidates.donor.seq.raw.fa \
    --OUT ${RAREVARDIR}/data/CRISPR/crispr.candidates.donor.seq.fa
```

Index the reference FASTA files for each amplicon region in the cDNA and gDNA
```
bwa index -p crispr/crispr.outlier.cdna.ref.fa crispr/crispr.outlier.cdna.ref.fa
bwa index -p crispr/crispr.control.cdna.ref.fa crispr/crispr.control.cdna.ref.fa
bwa index -p crispr/crispr.outlier.gdna.ref.fa crispr/crispr.outlier.gdna.ref.fa
bwa index -p crispr/crispr.control.gdna.ref.fa crispr/crispr.control.gdna.ref.fa
```

Perform mapping with BWA, sort and output to BAM format using samtools. <br>
Uses 10 cores. This is set at the top of the script.
```
bash crispr/crispr.bwa.aln.sort.sh
```

Summarize the CRISPR results
```
R CMD BATCH --no-save crispr/summarize.crispr.results.R
```
## RIVER
Here is an entire procedure of generating genomic features and expression values for running RIVER. To run these codes, you might need to move `gencode.v19.genes.v6p.patched_contigs.autosome.coding_linkRNA.gtf`,`gene_ensembl_ids.txt`,`tissue_names.txt` into `{RAREVARDIR}/reference and `an original WGS vcf file` into `{RAREVARDIR}/data/wgs`.
Note that you might need to change specified directory names into your own directories and revise some codes depending on which and how many features you would like to consider in your study. In order to run RIVER in bioconductor, you might also need to generate data following data structure and format in RIVER R package based upon final two files and a list of individual pairs having same rare variants within 10K from TSS of each gene.

#### [Step 1] Generate matlab annotation files with both indivs and genes considered in GTEx v6p and gene expression matrix used for calling outliers later
```
matlab -nodisplay -nodesktop -nosplash -nojvm -singleCompThread -r "run('${RAREVARDIR}/RIVER/code/generate_annotation_matlab.m')"
matlab -nodisplay -nodesktop -nosplash -nojvm -singleCompThread -r "run('${RAREVARDIR}/RIVER/code/generate_expmat_44tissues.m')"
```
#### [Step 2] In a WGS vcf file, all indels were removed, high quality variant calls (VQSLOD = PASS) were considered, only sites having <= 10 individuals in terms of missing genotypes were considered, and only autosomes were considered, and only European subjects were considered.
```
vcftools --gzvcf ${RAREVARDIR}/data/wgs/${an_original_GTEx_WGS_file} --keep ${RAREVARDIR}/preprocessing/gtex_2015-01-12_wgs_ids.txt --remove-indels --remove-filtered-all --max-missing-count 10 --not-chr X --recode --recode-INFO-all --stdout | bgzip -c > ${RAREVARDIR}/data/wgs/a_filtered_and_compressed_GTEx_WGS_vcf_file
tabix -p vcf ${RAREVARDIR}/data/wgs/filtered_and_compressed_GTEx_WGS_vcf_file
```
#### [Step 3] Compute GTEx allele frequencies with only subjects of interest
```
vcftools --gzvcf ${RAREVARDIR}/data/wgs/filtered_and_compressed_GTEx_WGS_vcf_file} --freq --out ${RAREVARDIR}/data/wgs/GTEx_af.vcf
bgzip -c ${RAREVARDIR}/data/wgs/GTEx_af.vcf > ${RAREVARDIR}/data/wgs/GTEx_af.vcf.gz
tabix -p vcf ${RAREVARDIR}/data/wgs/GTEx_af.vcf.gz
```
#### [Step 4] Extract a list of targeted regions for genes of interest
```
matlab -nodisplay -nodesktop -nosplash -nojvm -singleCompThread -r "run('${RAREVARDIR}/RIVER/code/generate_regions.m')"
```
#### [Step 5] In each subject of interest, extract a list of individual-specific rare variant sites based on AFs of both GTEx and EUR 1K population
```
count_ind=0
for ID in $(cat ${RAREVARDIR}/preprocessing/gtex_2015-01-12_wgs_ids_outlier_filtered.txt)
do
count_ind=$(( $count_ind + 1 ))
cat ${RAREVARDIR}/RIVER/data/rvsite/region.tss10k.txt | ${RAREVARDIR}/RIVER/code/extract_rvsites_ByInd.py -n $count_ind --id $ID --WGSvcf_in ${RAREVARDIR}/data/wgs/a_filtered_and_compressed_GTEx_WGS_vcf_file --GTExvcf_in ${RAREVARDIR}/data/wgs/a_compressed_GTEx_allele_frequency_file --EURvcf_in ${RAREVARDIR}/data/wgs/1KG/a_compressed_EUR_allele_freq_vcf_file --site_out ${RAREVARDIR}/RIVER/data/score/indiv/${ID}.${count_ind}.rvsites.txt
done
```
#### [Step 6] Extract all the features simulataneously (CADD, chromHMM, phylop, DANN).
```
count_ind=0
for ID in $(cat ${RAREVARDIR}/preprocessing/gtex_2015-01-12_wgs_ids_outlier_filtered.txt)
do
count_ind=$(( $count_ind + 1))
cat ${RAREVARDIR}/RIVER/data/score/indiv/${ID}.txt | ${RAREVARDIR}/RIVER/code/extract_scores_combined.py -n $count_ind --id $ID --af_in ${RAREVARDIR}/data/wgs/GTEx_af.vcf.gz --wgs_in filtered_and_compressed_GTEx_WGS_vcf_file --anno_in ${RAREVARDIR}/data/wgs/GTEx_vep.vcf.gz --cadd_in ${RAREVARDIR}/RIVER/data/whole_genome_SNVs_inclAnno.tsv.gz --dann_in ${RAREVARDIR}/RIVER/data/DANN_whole_genome_SNVs.tsv.bgz --chromHMM_in ${RAREVARDIR}/RIVER/data/wgEncodeBroadHmmGm12878HMM.sorted.hg19.bed.txt.gz --phylop_in ${RAREVARDIR}/RIVER/data/phyloP100way.txt.gz --score_out ${RAREVARDIR}/RIVER/data/score/indiv/${ID}.${count_ind}.score.nuc.txt
done
```
#### [Step 7] Extract quantitative values of genomic features from individual files generated by the previous script at gene level
```
count_ind=0
for ID in $(cat ${RAREVARDIR}/preprocessing/gtex_2015-01-12_wgs_ids_outlier_filtered.txt)
do
count_ind=$(( $count_ind + 1 ))
sed "s/order = 1;/order = $count_ind/" ${RAREVARDIR}/RIVER/code/extract_features_byInd.m > ${RAREVARDIR}/RIVER/code/mfiles/extract_features.$ID.$count_ind.m
matlab -nodisplay -nodesktop -nosplash -nojvm -singleCompThread -r "run('${RAREVARDIR}/RIVER/code/mfiles/extract_features.$ID.$count_ind.m')"
done
```
#### [Step 8] Compute normalized scores of genomic features
```
for order in {1..11}
do
sed "s/order = 1;/order = $order/" ${RAREVARDIR}/RIVER/code/compute_scores.m > ${RAREVARDIR}/RIVER/code/mfiles/compute_scores.$order.m
matlab -nodisplay -nodesktop -nosplash -nojvm -singleCompThread -r "run('${RAREVARDIR}/RIVER/code/mfiles/compute_scores.$order.m')"
done
```
#### [Step 9] Generate both "genomic_features.txt" and "zscores.txt" for running RIVER
```
matlab -nodisplay -nodesktop -nosplash -nojvm -singleCompThread -r "run('${RAREVARDIR}/RIVER/code/generate_data_RIVER.m')"
```


## Main figures
#### Figure 1
```
bash paper_figures/pick.cartoon.example.sh
Rscript paper_figures/figure1a.plot.cartoon.example.R
Rscript paper_figures/figure1b.outlier.sharing.R
Rscript paper_figures/figure1c.replication.rate.consistent.R
```

#### Figure 2
```
bash paper_figures/count_rarevars.sh
Rscript paper_figures/figure2a.count.enrichments.R
Rscript paper_figures/figure2b.threshold.enrichments.R
Rscript paper_figures/figure2c.ASE.enrichments.R
Rscript paper_figures/figure2.R
```

#### Figure 3
```
Rscript paper_figures/figure3a.rare.variant.class.enrichments.R
Rscript paper_figures/figure3b.feature.enrichments.R
Rscript paper_figures/figure3de.outlier.effect.size.R
Rscript paper_figures/figure3.R
```

#### Figure 4
```
Rscript paper_figures/figure4a.uk10k.R
Rscript paper_figures/figure4b.exac.enrichments.R
Rscript paper_figures/figure4c.gene.list.enrichments.R
Rscript paper_figures/figure4.R
```

#### Figure 5
```
Rscript getPosteriorsEval.R data/genomic_features.txt data/outliers.txt
Rscript paper_figures/figure5b.R
Rscript getPosteriorsApp.R data/genomic_features.txt data/outliers.txt data/postprobs_all.txt
Rscript paper_figures/figure5c.R
```

## Supplemental figures
#### EDF 1 and Supplementary Tables 1 and 2
```
## Adjusted R-squared values between top 15 PEER factors and top 20 sample and subject covariates in skeletal muscle
Rscript paper_figures/muscle_covariates_peerfactors.R
Rscript paper_figures/superheat_peer_muscle.R /data/muscle_samples_covariates_peerFactors.RData muscle.sample.peer.pdf 0.7
Rscript paper_figures/superheat_peer_muscle.R /data/muscle_subject_covariates_peerFactors.RData muscle.subject.peer.pdf 0.7

## Adjusted R-squared values between the total expression component removed by PEER in each of 44 tissues and top 20 sample and subject covariates
Rscript paper_figures/pve_samples_pertiss.R
Rscript paper_figures/pve_subject_pertiss.R
Rscript paper_figures/process_results.R
Rscript paper_figures/superheat_expression.R /data/superheat.subject.RData rv.subject.expression.pdf 0.25
Rscript paper_figures/superheat_expression.R /data/superheat.sample.RData rv.sample.expression.pdf 0.2

## Improvement of rare variant enrichments at varying levels of PEER correction
#Generate expression residuals with the top 0 or 5 PEER factors removed in addition to sex and genotype PCs
bash preprocessing/PEER/calc_residuals_covs_peerless.sh ${RAREVARDIR}/preprocessing/PEER 0
bash preprocessing/PEER/calc_residuals_covs_peerless.sh ${RAREVARDIR}/preprocessing/PEER 5

# Generate flat files for these new PEER corrected datasets to use for outlier calling
python preprocessing/gather_filter_normalized_expression_peerless.py \
       ${RAREVARDIR}/preprocessing/gtex_2015-01-12_tissues_all_normalized_samples.txt \
       ${RAREVARDIR}/preprocessing/gtex_2015-01-12_individuals_all_normalized_samples.txt \
       .peer.top0.ztrans.txt \
       ${RAREVARDIR}/preprocessing/gtex_2015-01-12_normalized_expression.peer.top0.txt
python preprocessing/gather_filter_normalized_expression_peerless.py \
       ${RAREVARDIR}/preprocessing/gtex_2015-01-12_tissues_all_normalized_samples.txt \
       ${RAREVARDIR}/preprocessing/gtex_2015-01-12_individuals_all_normalized_samples.txt \
       .peer.top5.ztrans.txt \
       ${RAREVARDIR}/preprocessing/gtex_2015-01-12_normalized_expression.peer.top5.txt

# Call Median Z-score outliers on data with top 0/5 PEER factors removed and compile features
R -f call_outliers/call_outliers_medz_peerless.R --slave --vanilla --args .peer.top0.txt
R -f call_outliers/call_outliers_medz_peerless.R --slave --vanilla --args .peer.top5.txt
bash feature_construction/run_compile_features_outliers_peerless.sh
```
Panel c is produced in `paper_figures/figure2a.count.enrichments.R`

#### EDF 2
You need to set the path to the subject annotations in the script below. <br>
```
Rscript paper_figures/suppfig.number.outliers.per.individual.R
```

#### EDF 3
```
Rscript paper_figures/suppfig.single.replication.compare.multi.R
```

#### EDF 4
You need to set the path to the downloaded expression covariates and subject annotations in the script below. <br>
```
Rscript paper_figures/suppfig.number.rare.vars.pca.R
```

#### EDF 5
```
Rscript paper_figures/suppfig.count.enrichments.R
```

#### EDF 6
Before running, set path to RNA-seq RPKMs in `preprocessing/get_median_rpkms.R ` correctly. <br>
```
Rscript preprocessing/get_median_rpkms.R
Rscript paper_figures/suppfig.over.under.expression.R
```

#### EDF 7
```
Rscript paper_figures/EDF7.R
```

#### EDF 8
Relies on the eGene and singificant variant-gene associations as downloaded for the GTEx v6p release from the portal.
```
bash paper_figures/annotate.variants.by.gene.sh
Rscript paper_figures/suppfig.rare.var.counts.disease.genes.gtex.cohort.R
Rscript paper_figures/suppfig.egene.enrichment.R
```

#### EDF 9
```
Rscript paper_figures/main_RIVER_VaryingThrds.R
Rscript paper_figures/Generate_figures_RIVER_VaryingThrds.R
Rscript paper_figures/main_RIVER_10CV.R
Rscript paper_figures/EDF9c.R
Rscript paper_figures/EDF9d.R
Rscript paper_figures/EDF9e.R
```

#### EDF 10
```
Rscript paper_figures/EDF10bd.R
Rscript paper_figures/EDF10ef.rpkm.R
Rscript paper_figures/EDF10ef.Zscores.R
```

### EDF 11
EDF11c is generated in `crispr/summarize.crispr.results.R`, run above.
