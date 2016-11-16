#!/usr/local/bin/bash
#
# Output the positions of snps within known loci.
# Filter .gen files to leave only the genotypes for these loci.

out_dir="/nfs/users/nfs_b/bb9/workspace/rotation1/crohns_workspace/5_model_comparison/"
mkdir -p "$out_dir"

# The list of snps contained within known loci per chromosome
contained_snps="$out_dir/contained.snps.txt"

# Output positions of snps within each chr
mkdir -p "$out_dir/snps_pos"
for ((chr = 1; chr <= 22; chr++)) {
    echo Getting positions of contained snps from chr $chr...
    awk -F'\t' -v chr="$chr" '($2 == chr) {print $27}' "$contained_snps" > "$out_dir/snps_pos/$chr.txt"
}

# Get corresponding rows from gen files
in_gen_files_dir="/nfs/users/nfs_b/bb9/workspace/rotation1/crohns_workspace/4_gwas/1_snptest/data/" 

# Set stripe on the output dir
mkdir -p "$out_dir/gen"
lfs setstripe "$out_dir/gen" -c -1 

for ((chr = 1; chr <= 22; chr++)) {
    date && echo Getting genotypes from chr $chr gen file...
    zcat "$in_gen_files_dir/$chr.gen.gz" | awk '{if (NR == FNR) {pos[$1] = 1} else {if (pos[$3]) {print}}}' "$out_dir/snps_pos/$chr.txt" - > "$out_dir/gen/$chr.gen"
}

