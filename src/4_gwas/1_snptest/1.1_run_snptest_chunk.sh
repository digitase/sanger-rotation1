#
# Run snptest on a specified chunk of the genome
#

OUT_DIR="$1"
dataset="$2"
assoc="$3"
chr="$4"
chunk_size="$5"
chunk_i="$LSB_JOBINDEX"

# Testing params
# OUT_DIR="/nfs/users/nfs_b/bb9/workspace/rotation1/crohns_workspace/4_gwas/1_snptest/"
# dataset="gwas3"
# assoc="ibd"
# chr=1
# chunk_size=200000
# chunk_i=20

# Imputed chromosome .gen.gz files
IN_GEN_FILES_DIR="/nfs/users/nfs_b/bb9/workspace/rotation1/crohns_workspace/4_gwas/1_snptest/data/"
# snptest binary
SNPTEST_BIN="/nfs/users/nfs_y/yl2/team143/software/snptest_v2.5-beta4_Linux_x86_64_dynamic/snptest_v2.5-beta4"
# .sample file with PCs for use as covariates
IN_PCS_SAMPLE_FILE="/lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc_NEW/GWAS3/refs/GWAS3.ibd.sample"
# samples to exclude when only testing for cd or uc specifically.
IN_EXCLUDE_SAMPLES_DIR="/lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc_NEW/GWAS3/refs/"
# snps to exclude per chromosome
IN_EXCLUDE_SNPS_DIR="/lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/info/GWAS3_new/"

# Create output dir for chromosome
CHUNK_OUT_DIR="$OUT_DIR/results/$dataset/$assoc/$chr/" 
mkdir -p "$CHUNK_OUT_DIR"

# Calculate start and end of chunk
start=$(echo "($chunk_i - 1) * $chunk_size" | bc)
end=$(echo "($chunk_i * $chunk_size) - 1" | bc)

#
# Run snptest
#
#            -data <a> <b>...: specify data files for analysis in .gen and .sample pairs. Automatic detection of .gz files.
#                      -o <a>: name of output file.
#     -frequentist <a> <b>...: specify which Frequentist tests to fit.
#                 -method <a>: method used to fit model, this can be one of "threshold", "expected", "score", "ml", "newml", or "em". 
#                  -pheno <a>: specify name of phenotype to use.
#         -cov_all_continuous: condition on all continuous covariates (C) in the sample files.
#           -range <a> <b>...: Operate only on SNPs within one of the given ranges.
# -exclude_samples <a> <b>...: list of samples to be excluded.
#           -exclude_snps <a>: list of snps to be excluded.
#
# The five different models are coded as 1=Additive, 2=Dominant, 3=Recessive, 4=General and 5=Heterozygote. 
# The additive model is the Cochran-Armitage test for an additive genetic effect. 
# The dominant and recessive models are specified using the AA genotype as the baseline genotype. 
# The general model is a the standard 2-df test of association.
#
# Some observations on runtime
    # -range "1-100" \
    # -frequentist 1 2 3 4 \

    # 14m to run just 1
    # 13m to run 1 2 3 4 

    # 65m to run a -range 200Kb chunk with 1 2 3 4 

date && echo Starting...
echo "LSB_JOBINDEX is $LSB_JOBINDEX"

chunk_output_file="$OUT_DIR/results/$dataset/$assoc/$chr/chunk_$chunk_i.range_${start}_$end.snptest.out" 

# Check if chunk already has output
# If chunk output file has non zero size, skip it.
# Note that regions with no association will have outputs with zero size.
# if [ -s "$chunk_output_file" ]; then
    # echo "Non-zero size output file detected: $chunk_output_file"
    # echo "Skipping chunk $LSB_JOBINDEX"
    # exit
# fi

# Run snptest
"$SNPTEST_BIN" \
    -data "$IN_GEN_FILES_DIR/$chr.gen.gz" "$IN_PCS_SAMPLE_FILE" \
    -exclude_samples "$IN_EXCLUDE_SAMPLES_DIR/gwas3-"$assoc"-assoc-sample-exclusion.txt" \
    -exclude_snps "$IN_EXCLUDE_SNPS_DIR/$chr-fail.list" \
    -range "$start-$end" \
    -frequentist 1 2 3 4 5 \
    -method score \
    -pheno bin1 \
    -cov_all_continuous \
    -o "$chunk_output_file" 

date && echo Done.

