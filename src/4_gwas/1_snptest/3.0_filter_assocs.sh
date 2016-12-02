#!/usr/local/bin/bash
#
# Concatenate together associations into one .assoc file per chromosome, filtering on:
#
# - all_maf >= 0.001 (overall, i.e. cases + controls)
# - info >= 0.4 (the imputation quality score)
# 
# INFO score
#
# The info measure takes the value 1 if all genotypes are completely certain, and
# the value 0 if the genotype probabilities for each sample are completely
# uncertain in Hardy-Weinberg proportions (i.e. they equal (1-θ)2, 2θ(1-θ), θ2).
# 
# MAF filtering overall vs. controls only:
#
# Yep, so just to clarify after crawling through old Slack discussions - we
# decided to filter on MAF_total > 0.1%, because we were attempting to capture
# all sites we thought we could reasonably impute. This would require the
# population-wide MAF to be 0.1% or greater (imputation is done on cases and
# controls together), and if we had limited it to controls only we are going to
# see some sites at < 0.1% (if it is at 0.1% in controls and, say, 0.05% in
# cases).
#
# The confusion came because for our sequencing-only data, where we have
# directly captured all variants and are not relying on imputation, we set a
# threshold for rare variant tests of  MAF <=0.5% in controls only. This is
# because we don’t want to exclude variants that might be rare in the
# population, but very elevated in IBD cases, because those are some of our
# most interesting findings.
#

# Which dataset was used
dataset="gwas3"
# Which samples to include
# Set to "cd", "uc", or "ibd"
# assoc="cd"
# assoc="uc"
assoc="ibd"

OUT_DIR="/nfs/users/nfs_b/bb9/workspace/rotation1/crohns_workspace/4_gwas/1_snptest/"
mkdir -p "$OUT_DIR/results.filtered/$dataset/$assoc/"

# Filtering thresholds
all_maf_thresh="0.001"
info_thresh="0.4"

# Concat chunks by chromosome and filter
assoc_header_line=$(grep '^alternate' $(find "$OUT_DIR/results/$dataset/$assoc/1/" -name "chunk_*" | head -n 1))
for ((chr = 1; chr <= 22; chr++)) {
    date && echo Merge and filter chr $chr associations...
    chr_assoc_file="$OUT_DIR/results.filtered/$dataset/$assoc/chr_$chr.snptest.filtered.out"
    # Add header
    echo "$assoc_header_line" > "$chr_assoc_file"
    # Filter and add chunks in order
    cat $(ls --sort=version "$OUT_DIR/results/$dataset/$assoc/$chr/chunk_"*) \
        | grep -v -e '^#' -e '^alternate_ids' \
        | awk -v all_maf_thresh="$all_maf_thresh" -v info_thresh="$info_thresh" '($9 >= info_thresh && $29 >= all_maf_thresh) {print $0}' \
            >> "$chr_assoc_file"
}

# Quick check of wc -l compared to katie's results
# for ((chr = 1; chr <= 22; chr++)) {
    # katie_assoc_wc=$(wc -l < "/lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc_NEW/GWAS3/ibd/$chr.assoc")
    # chr_assoc_file_wc=$(wc -l < "$OUT_DIR/results.filtered/$dataset/$assoc/chr_$chr.snptest.filtered.out")
    # echo katie $katie_assoc_wc, this $chr_assoc_file_wc, this-katie $(($chr_assoc_file_wc - $katie_assoc_wc))
# }

