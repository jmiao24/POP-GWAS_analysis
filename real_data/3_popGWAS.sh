#!/bin/bash

# Prepare the Regenie input into the POP-TOOLS format
cd POP-TOOLS

trait=Head_BMD

python3 ./POP-GWAS.py \
--gwas-yhat-unlab ./data/${trait}_yhat_unlab.txt.gz \
--gwas-y-lab ./data/${trait}_y_lab.txt.gz \
--gwas-yhat-lab ./data/${trait}_yhat_lab.txt.gz \
--out ./result/${trait}