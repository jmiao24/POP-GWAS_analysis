#!/bin/bash

# Clumping
./plink/plink_1.9_linux_x86_64/plink \
--bfile ./UKB_LD/genotype/ukb_10k \
--clump ./result/${trait}_POP-GWAS.txt \
--clump-field P \
--clump-p1 1.4e-8 \
--clump-p2 1.4e-8 \
--clump-r2 0.01 \
--clump-kb 5000 \
--out ./results/clumping/${trait}_popgwas

# Heritability & Genetic correlation
## LDSC munge
traits=(
    Arms
    Femur_neck
    Femur_shaft
    Femur_total
    Femur_troch
    Femur_wards
    Head
    L1-L4
    Legs
    Pelvis
    Ribs
    Spine
    Total
    Trunk
)

OutDir="./data/munged_gwas"

for trait in "${traits[@]}"; do
    filepath="./DXA/cleaned/${trait}_pop.txt.gz"
    output="${OutDir}/munged_${trait}_pop"

    ./bin/python ./Software/ldsc/ldsc/munge_sumstats.py\
        --merge-alleles ./Software/ldsc/Inputs/w_hm3.snplist\
        --sumstats $filepath\
        --ignore "Z"\
        --out $output &
done


# Cross-site Genetic correlation
# Put the following code in R:
library(data.table)
library(dplyr)
library(glue)

pop_gwas_path <- list.files("../data/munged_gwas", pattern = "pop.sumstats.gz", full.names = TRUE)

gwas_list_collapsed <- paste(pop_gwas_path, collapse = ",")

for(i in seq_along(pop_gwas_path)){
    index_gwas <- pop_gwas_path[i]
    index_gwas_trait <- sub(".*/munged_(.*?)_pop.\\.sumstats\\.gz", "\\1", index_gwas)
    cmd <- paste0(
        "./bin/python ./Software/ldsc/ldsc/ldsc.py",
        " --rg ", index_gwas, ",", gwas_list_collapsed,
        " --ref-ld-chr ./Software/ldsc/Inputs/EUR_1KGphase1/eur_w_ld_chr/",
        " --w-ld-chr ./Software/ldsc/Inputs/EUR_1KGphase1/eur_w_ld_chr/",
        " --out ../data/gen_cor/", index_gwas_trait, "_vs_all.pop &")

    system(cmd)
}

traits=(
    Arms
    Femur_neck
    Femur_shaft
    Femur_total
    Femur_troch
    Femur_wards
    Head
    L1-L4
    Legs
    Pelvis
    Ribs
    Spine
    Total
    Trunk
)
# make directory for each trait

for trait in "${traits[@]}"; do

    mkdir $PWD/data/gen_cor_vs_other_traits/${trait}

done

# create .R file for each trait to run ldsc gencor
for trait in "${traits[@]}"; do
    # Create a unique .R file for each trait
    echo "source(\"./Tools/GeneticCorrelation/GenCorrLDSC.R\")" > "$PWD/data/gen_cor_vs_other_traits/${trait}/ldsc_gencor.R"
    echo "" >> "$PWD/data/gen_cor_vs_other_traits/${trait}/ldsc_gencor.R"
    echo "ss <- c(" >> "$PWD/data/gen_cor_vs_other_traits/${trait}/ldsc_gencor.R"
    echo "    \"./data/munged_gwas/munged_${trait}_pop.sumstats.gz\"" >> "$PWD/data/gen_cor_vs_other_traits/${trait}/ldsc_gencor.R"
    echo ")" >> "$PWD/data/gen_cor_vs_other_traits/${trait}/ldsc_gencor.R"
    echo "ss_name <- c(\"${trait}\")" >> "$PWD/data/gen_cor_vs_other_traits/${trait}/ldsc_gencor.R"
    echo "working_folder <- \"temp_ldsc\"" >> "$PWD/data/gen_cor_vs_other_traits/${trait}/ldsc_gencor.R"
    echo "multiple_ldsc_gencorr_func(working_folder, ss, ss_name)" >> "$PWD/data/gen_cor_vs_other_traits/${trait}/ldsc_gencor.R"
done

# SMR