#!/usr/local/bin/bash
#
# Run 2.1_run_glm_models.R on a single chunk of snps

out_dir="/nfs/users/nfs_b/bb9/workspace/rotation1/crohns_workspace/5_model_comparison/"
mkdir -p "$out_dir/chunks"

# Which dataset was used
dataset="gwas3"
# Which samples to include
# Set to "cd", "uc", or "ibd"
assoc="ibd"
# Chunk size in number of SNPs
chunk_size=100

#
# Which chromosomes to run
#
for ((chr = 1; chr <= 22; chr++)) {

    gen_file="/nfs/users/nfs_b/bb9/workspace/rotation1/crohns_workspace/5_model_comparison/gen/$chr.gen"
    n_snps=$(cat $gen_file.[0-9]* | wc -l)
    n_chunks=$(echo $(printf %.0f $(echo $n_snps "/ $chunk_size" | bc)) "+1" | bc)
    date && echo "Chrom $chr ($n_snps snps) will run in $n_chunks chunks of $chunk_size snps."

    # Make logging dir
    mkdir -p "$out_dir/logs/$dataset/$assoc/$chr/"

    # bsub a job array over the chunks
    # Mem usage for chunk_size=100 ~ 300-400 MB
    # Restart on termination by memory limit
    bsub \
        -G team152 -q normal \
        -R "select[mem>500] rusage[mem=500]" -M 500 \
        -J "$dataset.$assoc.$chr[1-$n_chunks]" \
        -o "$out_dir/logs/$dataset/$assoc/$chr/jobid_%J.i_%I.bsub_o.log" \
        -e "$out_dir/logs/$dataset/$assoc/$chr/jobid_%J.i_%I.bsub_e.log" \
        "/software/R-3.3.0/bin/Rscript 2.1_run_glm_models.R $dataset $assoc $chr $chunk_size $gen_file $n_snps $out_dir"

}

