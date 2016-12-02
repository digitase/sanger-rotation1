#!/usr/local/bin/bash
#
# Output the ids of significant snps within known loci.
# Filter .gen files to leave only the genotypes for these loci.
# Write filtered .gen files in chunks.

out_dir="/nfs/users/nfs_b/bb9/workspace/rotation1/crohns_workspace/4_gwas/2_R_glm/"
mkdir -p "$out_dir"

# The list of snps contained within known loci per chromosome
contained_snps="$out_dir/contained.snps.all.txt"
# contained_snps="$out_dir/contained.snps.signif.txt"

# Output unique snp ids within each chr
mkdir -p "$out_dir/snp_ids"
for ((chr = 1; chr <= 22; chr++)) {
    echo Getting ids of contained snps from chr $chr...
    awk -F'\t' -v chr="$chr" '($1 == chr) {print $24}' "$contained_snps" | sed 's/"//g' | sort -u > "$out_dir/snp_ids/$chr.txt"
}

# Get corresponding rows from gen files
in_gen_files_dir="/nfs/users/nfs_b/bb9/workspace/rotation1/crohns_workspace/4_gwas/1_snptest/data/" 

# Set stripe on the output dir
mkdir -p "$out_dir/gen/logs"
lfs setstripe "$out_dir/gen" -c -1 

# Get contained loci genotypes.
# Split output files.
chunk_size=100

for ((chr = 1; chr <= 22; chr++)) {
    bsub \
        -n 2 \
        -R "span[hosts=1] select[mem>1000] rusage[mem=1000]" -M 1000 \
        -o "$out_dir/gen/logs/$chr.gen.bsub_o.log" \
        -e "$out_dir/gen/logs/$chr.gen.bsub_e.log" \
        "zcat $in_gen_files_dir/$chr.gen.gz | awk '{if (NR == FNR) {id[\$1] = 1} else {if (id[\$2]) {print}}}' $out_dir/snp_ids/$chr.txt - | \
            split -a 10 -d -l $chunk_size - $out_dir/gen/$chr.gen."
}

