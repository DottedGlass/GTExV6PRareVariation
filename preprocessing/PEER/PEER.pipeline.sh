#!/bin/bash

# PEER correction pipeline for outlier detection and rare variant analysis
set -o nounset -o errexit -o pipefail

export PEER_DIR=${RAREVARDIR}/preprocessing/PEER

# Make RPKM and read count matrices for each tissue
Rscript preprocessing/PEER/preprocess.expr.data.R

# Estimate PEER factors
bash preprocessing/PEER/calc.PEER.factors.all.tissues.sh

# Remove residuals
bash preprocessing/PEER/calc.residuals.sh
