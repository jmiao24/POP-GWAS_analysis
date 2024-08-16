#!/bin/bash

## see more from: https://rgcgithub.github.io/regenie/recommendations/
## regenie documentations: https://rgcgithub.github.io/regenie/options/

## Replace the following variables with your variable of interest
# trait=Arms
# pheno=lab.y_hat

## Step1, stage1
mkdir -p ./data/gwas_regenie/DXA_BMD/regenie/step0

./plink \
--bfile ./UKB/genotype/ukb_EUR_pedigree_finetuned \
--keep ./data/pheno/DXA_BMD/${trait}.${pheno}_eur.txt \
--autosome \
--snps-only just-acgt \
--geno 0.01 \
--hwe 0.000001 \
--maf 0.01 \
--allow-no-sex \
--write-snplist \
--out ./data/gwas_regenie/DXA_BMD/regenie/step0/${trait}_${pheno}.qc.pass

mkdir -p ./data/gwas_regenie/DXA_BMD/regenie/step1

./regenie/v3.2.1/regenie \
--step 1 \
--gz \
--threads 8 \
--bed ./UKB/genotype/ukb_EUR_pedigree_finetuned \
--extract ./data/gwas_regenie/DXA_BMD/regenie/step0/${trait}_${pheno}.qc.pass.snplist \
--keep ./data/pheno/DXA_BMD/${trait}.${pheno}_eur.txt \
--phenoFile ./data/pheno/DXA_BMD/${trait}.${pheno}_eur.txt \
--phenoColList ${pheno} \
--covarFile ./data/covar/DXA_BMD/${trait}.${pheno}_eur.txt \
--covarColList Year,Sex,SexYear,Chip,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20 \
--catCovarList Sex,Chip \
--maxCatLevels 100 \
--bsize 1000 \
--split-l0 ./data/gwas_regenie/DXA_BMD/regenie/step1/${trait}_${pheno}_fit_bin_parallel,100 \
--out ./data/gwas_regenie/DXA_BMD/regenie/step1/${trait}_${pheno}_fit_bin_l0

## Step1, stage 2
./regenie/v3.2.1/regenie \
--step 1 \
--gz \
--threads 8 \
--bed ./UKB/genotype/ukb_EUR_pedigree_finetuned \
--extract ./data/gwas_regenie/DXA_BMD/regenie/step0/${trait}_${pheno}.qc.pass.snplist \
--keep ./data/pheno/DXA_BMD/${trait}.${pheno}_eur.txt \
--phenoFile ./data/pheno/DXA_BMD/${trait}.${pheno}_eur.txt \
--phenoColList ${pheno} \
--covarFile ./data/covar/DXA_BMD/${trait}.${pheno}_eur.txt \
--covarColList Year,Sex,SexYear,Chip,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20 \
--catCovarList Sex,Chip \
--maxCatLevels 100 \
--bsize 1000 \
--run-l0 ./data/gwas_regenie/DXA_BMD/regenie/step1/${trait}_${pheno}_fit_bin_parallel.master,${job} \
--out ./data/gwas_regenie/DXA_BMD/regenie/step1/${trait}_${pheno}_fit_bin_l0_${job}

## Step1, stage 3
./regenie/v3.2.1/regenie \
--step 1 \
--gz \
--threads 32 \
--bed ./UKB/genotype/ukb_EUR_pedigree_finetuned \
--extract ./data/gwas_regenie/DXA_BMD/regenie/step0/${trait}_${pheno}.qc.pass.snplist \
--keep ./data/pheno/DXA_BMD/${trait}.${pheno}_eur.txt \
--phenoFile ./data/pheno/DXA_BMD/${trait}.${pheno}_eur.txt \
--phenoColList ${pheno} \
--covarFile ./data/covar/DXA_BMD/${trait}.${pheno}_eur.txt \
--covarColList Year,Sex,SexYear,Chip,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20 \
--catCovarList Sex,Chip \
--maxCatLevels 100 \
--bsize 1000 \
--run-l1 ./data/gwas_regenie/DXA_BMD/regenie/step1/${trait}_${pheno}_fit_bin_parallel.master \
--out ./data/gwas_regenie/DXA_BMD/regenie/step1/${trait}_${pheno}

# Stage 2:
./plink \
--bfile ./genotype/chunks/chr${chr}.EUR.idupdated_snpid_${chunk}_50 \
--keep ./data/pheno/DXA_BMD/${trait}.${pheno}_eur.txt \
--geno 0.01 \
--hwe 0.000001 \
--maf 0.01 \
--write-snplist \
--out ./data/gwas_regenie/DXA_BMD/regenie/step1/${trait}_${pheno}_qc_pass_chr${chr}_chunk${chunk}

mkdir -p ./data/gwas_regenie/DXA_BMD/regenie/step2

./regenie/v3.2.1/regenie \
--step 2 \
--gz \
--threads 5 \
--bsize 400 \
--pred ./data/gwas_regenie/DXA_BMD/regenie/step1/${trait}_${pheno}_pred.list \
--bed ./genotype/chunks/chr${chr}.EUR.idupdated_snpid_${chunk}_50 \
--extract ./data/gwas_regenie/DXA_BMD/regenie/step1/${trait}_${pheno}_qc_pass_chr${chr}_chunk${chunk}.snplist \
--keep ./data/pheno/DXA_BMD/${trait}.${pheno}_eur.txt \
--phenoFile ./data/pheno/DXA_BMD/${trait}.${pheno}_eur.txt \
--phenoColList ${pheno} \
--covarFile ./data/covar/DXA_BMD/${trait}.${pheno}_eur.txt \
--covarColList Year,Sex,SexYear,Chip,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20 \
--catCovarList Sex,Chip \
--maxCatLevels 100 \
--out ./data/gwas_regenie/DXA_BMD/regenie/step2/gwas_linear_${trait}_${pheno}_chr${chr}_snpid_${chunk}_50