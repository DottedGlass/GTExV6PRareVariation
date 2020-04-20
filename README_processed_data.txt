#------ GTEx v6p Rare Variation Data for Public Release
Contacts: Alexis Battle (ajbattle@cs.jhu.edu), Stephen Montgomery (smontgom@stanford.edu)
Link to the Github repo: https://github.com/joed3/GTExV6PRareVariation
Download link for this processed dataset: https://s3-us-west-2.amazonaws.com/gtex-v6p-rare-variation-data/GTExV6PRareVariationData.tar.gz
Github repo states that users should download this processed data into a directory named <processed_data>

#------ File list
This processed_data directory contains 4 subdirectories:
	1. data
	2. features
	3. preprocessing
	4. reference

The files in each of these directories are listed as follows. Order follows alphabetically according to file name

#-- data subdirectory files

#-- data/genes.rpkm.mean.mat.txt
Description: Matrix of mean expression (RPKM) for each gene across tissues
File format: Tab delimited file where each row is a gene and each column is a tissue. Values represent mean RPKM values across individuals

#-- data/genes.rpkm.summary.stats.txt
Description: Summary statistics on mean, median, and sd of gene expression across tissues for each gene
File format: Tab delimited format where each row is a gene and each column is the mean, median, and standard deviation of mean RPKM values across the 44 tissues

#-- data/outliers_medz_counts.txt
Description: Matrix of genes by individuals representing the number of tissues used in computing the median Z-score for each (individual, gene) pair
File format: Tab delimited file where rows are genes (ENSG IDs) and columns are individuals (GTEx individual IDs). Values represent the number of tissues for the given (individual, gene) pair used in computation of the median Z-score

#-- data/outliers_medz_nothreshold_picked.txt
Description: Most extreme multi-tissue outliers for each gene with criterion |median Z-score| > 0
File format: Tab delimited file where first column is the gene ENSG ID, second column GTEx individual ID, third column the number of tissues included in the median Z-score computation, and the fourth column the median Z-score

#-- data/outliers_medz_picked.txt
Description: Most extreme multi-tissue outliers for each gene with criterion |median Z-score| > 2
File format: Tab delimited file where first column is the gene ENSG ID, second column GTEx individual ID, third column the number of tissues included in the median Z-score computation, and the fourth column the median Z-score

#-- data/outliers_medz_picked_counts_per_ind.txt
Description: Number of genes where the GTEx individual is a multi-tissue outlier
File format: Tab delimited file where the first column is the GTEx individual ID and the second column is the number of genes where that individual is a multi-tissue outlier

#-- data/outliers_medz_zscores.txt
Description: Matrix of genes by individuals holding the median Z-score for each (individual, gene) pair
File format: Tab delimited file where rows are genes (ENSG IDs) and columns are individuals (GTEx individual IDs). Values represent the median Z-score for the given (individual, gene) pair. For (individual, gene) pairs where the number of available tissues is less than 5, NA is used to indicate the outlier test was not performed

#-- data/outliers_singlez_picked.txt
Description: All single tissue outliers for each expressed autosomal lincRNA or protein coding gene in each tissue
File format: Tab delimited file where first column is gene ENSG ID, second column is GTEx individual ID, third column is tissue ID (spaces delimited by underscores), and fourth column is the expression Z-score

#-- data/singlez/outliers_singlez_nothreshold_{tissue name}_counts.txt
Description: Binary matrix of genes by individuals representing whether the individual has expression data for that gene in the given tissue
File format: Tab delimited file where rows are genes (ENSG IDs) and columns are GTEx individuals. Values are either 0 or 1 depending on whether expression was measured for that gene in that individuals in the given tissue. ${tissue name} represents 1 of the 44 GTEx tissue names

#-- data/singlez/outliers_singlez_nothreshold_${tissue name}_picked.txt
Description: Most extreme outlier for each expressed autosomal lincRNA or protein coding gene in each tissue with the criterion that |Z-score| > 0
File format: Tab delimited file where first column is gene ENSG ID, second column is GTEx individual ID, third column is the tissue ID, and fourth column is the Z-score. ${tissue name} represents 1 of the 44 GTEx tissue names

#------ features subdirectory files

#-- features/annotations/ACMG/completed_acmg.csv
Description: List of genes prioritized by the ACMG for having actionable, incidental findings
File format: Comma delimited file with format described in the header. Raw file downloaded from http://www.ncbi.nlm.nih.gov/clinvar/docs/acmg/

#-- features/annotations/DDG2P/DDG2P_2_8_2017.csv.gz
Description: List of disease associated genes as compiled by DDG2P
File format: Tab delimited file with format described in the header. Raw file downloaded from http://www.ebi.ac.uk/gene2phenotype/

#-- features/annotations/GWAS/gwas_catalog_v1.0-downloaded_2015-11-30.tsv
Description: List of disease and trait associated genes as compiled by the NHGRI GWAS catalog
File format: Tab delimited file with format described in the header. Downloaded from http://www.ebi.ac.uk/gwas

#-- features/annotations/GeneListsOther/cancer.genes.gold.standard.csv
Description: List of genes associated with heritable cancers
File format: Comma delimited file with format described in the header. Raw file obtained as descrived in the methods section of the manuscript

#-- features/annotations/GeneListsOther/cardio.genes.gold.standard.csv
Description: List of genes associated with heritable cardiovascular disease
File format: Comma delimited file with format described in the header. Raw file obtained as described in the methods section of the manuscript

#-- features/annotations/GeneListsOther/gene_condition_source_id
Description: List of disease associated genes as compiled by ClinVar
File format: Tab delimited file with format described in the header. Downloaded from http://www.ncbi.nlm.nih.gov/clinvar/

#-- features/annotations/OMIM/omim.genes.txt
Description: List of disease associated genes as compiled by OMIM
File format: Tab delimited file where first column is gene symbol and second column is the ENSG ID. Raw file was downloaded from http://omim.org/

#-- features/annotations/OrphaNet/orphanet.genes.txt
Description: List of disease associated genes (ENSG IDs) as compiled by OrphaNet
File format: List of ENSG IDs. Raw file was downloaded from  http://www.orpha.net/

#------ preprocessing subdirectory files

#-- preprocessing/PEER/${tissue name}.peer.ztrans.txt
Description: Expression residuals for each tissue with PEER factors and known covariates (sex, ancestry) removed
File format: Tab delimited file where each row represents a gene and each column represents an individual with expression data for the tissue. Values represent the expression Z-scores for each gene and individual with PEER factors, sex, and ancestry covariates removed.

#-- preprocessing/PEER/${tissue name}.reads.txt
GTEx v6p read count matrices for each tissue: preprocessing/PEER/${tissue name}.reads.txt
Description: GTEx v6p read count matrices for each tissue
File format: Tab delimited file. First column is the gene name (ENSG ID). Remaining columns represent the individuals with data for that tissue. Values for each gene and individual pair are read counts. ${tissue_name} represents one of the 44 tissues in this analysis set

#-- preprocessing/PEER/${tissue name}.rpkm.log2.ztrans.txt
Description: log2(RPKM + 2) transformed and centered and scaled expression matrices for each tissue
File format: Tab delimited file where each row represents a GTEx individual and column represents a gene. Values are log2(RPKM + 2) transformed followed by centering and scaling. ${tissue_name} represents one of the 44 tissues in this analysis set

#-- preprocessing/PEER/${tissue_name}.rpkm.txt
Description: GTEx v6p RPKM matrices for each tissue
File format: Tab delimited file. First column is the gene name (ENSG ID). Remaining columns represent the individuals with data for that tissue. Values for each gene and individual pair are RPKM values. ${tissue_name} represents one of the 44 tissues in this analysis set

#-- preprocessing/PEER/${tissue name}_Factors${number of factors}
Description: subdirectories for each tissue containing PEER factors, residuals, and QC plots
File format: ${tissue name} represents one of the 44 tissues in this analysis set. ${number of factors} depends on the sample size for the tissue with values in {15, 20, 25, 30, 35}

#-- preprocessing/gtex.expressed.genes.txt
Description: List of expressed genes (ENSG IDs)
File format: List of expressed genes (ENSG IDs) in any tissue in GTEx

#-- preprocessing/gtex_2015-01-12_individuals_all_normalized_samples.txt
Description: List of all individuals used in outlier calling
File format: List of individual IDs for individuals included in outlier calling with each ID on a separate row

#-- preprocessing/gtex_2015-01-12_normalized_expression.txt
Description: Flat file containing the expression Z-scores for each gene in each tissue across individuals
File format: Tab delimited file where first two columns represent tissue names and gene IDs (ENSG). Remaining columns hold the expression Z-scores for the 449 individuals with expression data in the 44 tissues

#-- preprocessing/gtex_2015-01-12_rpkm.txt
Description: Flat file containing the RPKM for each gene in each tissue across individuals
File format: Tab delimited file where first two columns represent tissue names and gene IDs (ENSG). Remaining columns hold the RPKM values for the 449 individuals with expression data in the 44 tissues

#-- preprocessing/gtex_2015-01-12_samples_tissues.txt
Description: Map of GTEx sample IDs to tissue type
File format: Tab delimited file with two columns. First column is GTEx sample ID, second column is tissue name (spaces delimited by underscores)

#-- preprocessing/gtex_2015-01-12_tissue_by_ind.txt
Description: Mapping of tissue names and individual IDs used in outlier calling
File format: Tab delimited file with two columns. First column is the tissue name (spaces delimited by underscores). Second colun is the GTEx individual ID

#-- preprocessing/gtex_2015-01-12_tissues_all_normalized_samples.txt
Description: List of all tissues used in outlier calling
File format: List of tissue names for tissues included in outlier calling with each name on a separate row

#-- preprocessing/gtex_2015-01-12_wgs_ids_HallLabSV_outlier_filtered.txt
Description: List of GTEx individuals with WGS data filtered for individuals who were multi-tissue outliers for 50 or more genes. Additional filter for individuals with large chromosomal aberrations. This is a filtered version of preprocessing/gtex_2015-01-12_wgs_ids_HallLabSV.txt
File format: List of GTEx individual IDs each on a separate line

#-- preprocessing/gtex_2015-01-12_wgs_ids_outlier_filtered.txt
Description: List of GTEx individuals with WGS data filtered for individuals who were multi-tissue outliers for 50 or more genes. This is a filtered version of preprocessing/gtex_2015-01-12_wgs_ids.txt
File format: List of GTEx individual IDs each on a separate line

#------ reference subdirectory files

#-- reference/gencode.v19.annotation_coding.lincRNA_padded.bed
Description: BED file with lincRNA and protein coding exons padded by 5 bp
File format: BED file representing the exons of autosomal lincRNA and protein coding genes with exons padded by 5bp

#-- reference/gencode.v19.genes.v6p.patched_contigs.bed
Description: BED file with gene TSS and TES
File format: BED file with 4th column giving the gene ENSG ID

#-- reference/gencode.v19.genes.v6p.patched_contigs_TSS.bed
Description: BED file with TSS positions for each gene with annotation in GENCODE V19
File format: BED file with 4th column giving the gene ENSG ID

#-- reference/gencode.v19.genes.v6p.patched_contigs_genetypes_autosomal.txt
Description: Mapping of gene names (ENSG IDs; autosomes only) to types (e.g. lincRNA, pseudogene, etc.)
File format: Tab delimited file with first column gene ENSG ID (autosomes only) and second column gene type (e.g. lincRNA, pseudogene, protein_coding, etc.)
