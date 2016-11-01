
# Based on the methods from:
# /lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/pcs/README

# https://groups.google.com/forum/#!topic/plink2-users/gVhn1iu_X2k
#
# smartpca normally alternates between computing PCs and removing outliers from
# the dataset, where outlier is defined as "more than 6 sigma from the mean along
# one of the top 10 just-computed PCs", until either no more outliers exist or
# outliers have already been removed 5 times.  This is frequently a good idea,
# and a reason to prefer smartpca to plink's basic --pca.  However, if you want
# directly comparable results, you can disable outlier removal by adding "-m 0"
# to your smartpca.perl invocation.

# PC calculation for 3 GWAS cohorts
#WTCCC1
gwas1="/lustre/scratch114/teams/barrett/coreex_gaibdc/imputation/imputation_panel_A_WTCCC1.bcf"
#WTCCC2
gwas2="/lustre/scratch113/projects/crohns/GWAS3/imputation/imputation_panel_B_WTCCC2_reallyclean.bcf"
#GWAS3
gwas3="/lustre/scratch113/projects/crohns/GWAS3/imputation/imputation_panel_C_GWAS3.bcf"
#new GWAS3
gwas3="/lustre/scratch113/projects/crohns/GWAS3/imputation/imputation_panel_C_GWAS3_NEW.bcf"
#Convert to PLINK FORMAT
proj="gwas3"

IN_BCF=""

gwas3="/lustre/scratch113/projects/crohns/GWAS3/imputation/"

# Convert into plink format
plink \
    --bcf $'$proj'
    --const-fid 
    --biallelic-only
    --make-bed 
    --out 

# Prune SNPS
plink 
    --bcf "$prok"
    --bfile '$proj' 
    --biallelic-only 
    --exclude range "/lustre/scratch113/projects/tb/RUSSIA/data/combined/high-LD-regions.txt"
    --indep-pairwise 50 5 "0.2"
    --out '$proj|bsub -J "$proj-prune"

# Prune SNPs so that no pair of SNPs (within a given number of base pairs) has an R2 value greater than a given threshold
plink \
    --bfile "$OUT_DIR/$IN_DATA_PREFIX.qc3" \
    --exclude range "/lustre/scratch113/teams/barrett/coreex_gaibdc/refs/high-LD-regions.txt" \
    --indep-pairwise 50 5 "0.2" \
    --out "$OUT_DIR/$IN_DATA_PREFIX.qc3.prune_markers_by_ld"


plink 
--bfile '$proj' 
--maf 0.01 
--extract '$proj'.prune.in 
--pca 
--out '$proj|bsub -w "done($proj-prune)" 

