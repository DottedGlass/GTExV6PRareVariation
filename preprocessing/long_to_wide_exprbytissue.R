#!/usr/bin/env Rscript

require(dplyr)
require(tidyr)

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)<2) {
  stop("At least two argument must be supplied (input file, output file).n", call.=FALSE)
}

infile <- args[1]
outfile <- args[2]
orig <- read.table(infile,header = TRUE)
fix <- orig %>% spread(key = Ind,value = Zscore) %>% rename(Id=Gene)
write.table(fix,file=outfile,quote = FALSE,row.names = FALSE,col.names = TRUE,sep="\t")


