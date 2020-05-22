#!/usr/bin/env Rscript

# created the carcass of figure 1. the rest is done in inkscape.

library(ggplot2)
library(reshape2)
require(dplyr)

dir = Sys.getenv('RAREVARDIR')
GENE_CHOICE_v6='ENSG00000198610.6'
GENE_CHOICE='ENSG00000004139.13'

expr = read.table('expr.subset.by.genes.txt', header = T, stringsAsFactors = F)
examples = read.table('possible.txt', header = F, stringsAsFactors = F)
colnames(examples) = c('gene','ind','ntissue','medz')

# limit to tissue for which we have a pretty cartoon
tissues = c("Artery_Aorta","Artery_Coronary","Esophagus_Gastroesophageal_Junction","Heart_Atrial_Appendage","Heart_Left_Ventricle","Liver","Lung","Skin_Not_Sun_Exposed_Suprapubic","Skin_Sun_Exposed_Lower_leg","Stomach","Whole_Blood")

# for each gene, individual pair print the tissues if they are all contained in this set
get.tis = function(gene, ind) {
    indname = sub('-','.', ind, fixed = T)
    subset = expr[expr$Gene == gene, c('Tissue','Gene',indname)]
    subset = subset[!is.na(subset[,3]), ]
    stopifnot(nrow(subset) == 5)
    if (sum(subset$Tissue %in% tissues) == 4) {
        print(subset)
    }
}

dummy = mapply(get.tis, examples$gene, examples$ind)

expr.melted = melt(expr, id.vars = c('Tissue','Gene'))

expr.melted = expr.melted[!is.na(expr.melted$value), ]

# expr.melted.counts = expr.melted %>% group_by(Gene) %>% summarize(Tcount=n_distinct(Tissue))

#For Now we will look at a gene that has a unique cartoon for the exact same tissues as the v6 paper
indtissues = c("Adipose_Subcutaneous","Liver","Lung","Stomach","Whole_Blood")

# expr.melted.genes = expr.melted %>% filter(Tissue %in% indtissues) %>%
#   group_by(Gene) %>%
#   summarize(Tcount=n_distinct(Tissue)) %>% filter(Tcount == 5)


# picked one that would have a unique cartoon for each tissue

plot.data = expr.melted[expr.melted$Gene == GENE_CHOICE & expr.melted$Tissue %in% indtissues, ]

# for that gene, get the medz for all individual (with at least 5 tissues)
gene.expr = expr[expr$Gene == GENE_CHOICE, -c(1,2)]
gene.medz = apply(gene.expr, 2, function(x) ifelse(sum(!is.na(x)) >= 5, yes = median(x, na.rm = TRUE), no = NA))
gene.medz = gene.medz[!is.na(gene.medz)]
gene.medz.df = data.frame(Tissue = "Median", Gene = GENE_CHOICE, variable = names(gene.medz), value = gene.medz)

plot.data = rbind(plot.data, gene.medz.df)
plot.data$Tissue = factor(plot.data$Tissue, levels = c("Adipose_Subcutaneous","Liver","Stomach","Whole_Blood","Lung","Median"))

cols = c("#FF6600","#AABB66","#99FF00","#FFDD99","#FF00BB")
names(cols) = indtissues

tissue.labs <- gsub("_","\n",c("Adipose_Subcutaneous","Liver","Stomach","Whole_Blood","Lung","Median"))
names(tissue.labs) <- c("Adipose_Subcutaneous","Liver","Stomach","Whole_Blood","Lung","Median")

#pdf('./figure1a.cartoon.draft_v8.pdf', height = 7, width = 4.5)

pdf(paste0(dir, '/paper_figures/figure1a.cartoon.draft_v8.pdf'), height = 7, width = 4.5)
plot.data %>% 
ggplot(., aes(x = value)) +
    geom_histogram(binwidth = 0.15, colour = "white", fill = "darkgrey") + xlab('Z-score') + ylab('') +
    facet_grid(Tissue~.,scales = "free", labeller=labeller(Tissue=tissue.labs)) + theme_classic() + guides(fill = FALSE) + 
  xlim(c(-5,5)) +
    # geom_vline(xintercept = c(-0.9624364,-0.3787555,4.6037353,3.2303840,3.8570197), size = 1.1) +
    scale_fill_manual(values = cols) +
    theme(axis.ticks.y = element_blank(),
          #axis.text.y = element_blank(),
          strip.text = element_text(size = 8),
          strip.background = element_blank(),
          axis.text.x = element_text(size = 12),
          axis.title.x = element_text(size = 13),
          axis.title.y = element_text(size = 13))

dev.off()
