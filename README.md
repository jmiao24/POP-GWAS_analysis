# POP-GWAS analysis

This repository contains source code to reproduce analyses in "Valid inference for machine learning-assisted genome-wide association studies". 

The official software package is in [POP-TOOLS](https://github.com/qlu-lab/POP-TOOLS) GitHub repo.

## File structures
### Simulations
* `./simulation/Fun.R`: Functions for simulation and simple implementation of POP-GWAS in R
* `./simulation/qt.R`: simulation for quantitative phenotype
* `./simulation/bt.R`: simulation for binary phenotype
* `./simulation/imputation_r_FPR.R`: type-I error simulation for varying imputation correlation 
* `./simulation/imputation_r_power.R`: power simulation for varying imputation correlation
* `./simulation/vary_ratio.R`: simulation for varying sample size of unlabeled data

### UK Biobank data analysis
* `./real_data/1_softimpute.R`: phenotype imputation using softimpute
* `./real_data/2_regenie.sh`: run GWAS using regenie
* `./real_data/3_popGWAS.sh`: apply POP-GWAS to summary statistics
* `./real_data/4.1_post-GWAS.sh`: estimating heritability and genetic correlation using LD score regression
* `./real_data/4.2_coloc.R`: colocalization analysis

Example data for reproducing the results is available for download on [Box](https://uwmadison.box.com/s/iede9965ic49yol1r8w7bbx1ngkcrht0). After downloading, place the files in the 'data' folder and execute the provided codes.
  
## Contact
Feel free to reach out to Jiacheng at jmiao24@wisc.edu for questions.
