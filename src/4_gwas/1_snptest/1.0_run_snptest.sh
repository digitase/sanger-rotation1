#!/usr/local/bin/bash
#
# Run GWAS with snptest
# 
# Runs the analysis on chunks of each chromosome
#

OUT_DIR="/nfs/users/nfs_b/bb9/workspace/rotation1/crohns_workspace/4_gwas/1_snptest/"
mkdir -p "$OUT_DIR"

# Gauge resource requirements from the longest running of previously run jobs
# Running -frequentist 1, CHUNK_SIZE=2000000, max time~22h, max mem~1.6G
# zgrep -a "CPU time" /lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc_NEW/GWAS3/ibd/calling/logs.tar.gz | awk '{print $4}' | sort -n | tail
# zgrep -a "Max Memory" /lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc_NEW/GWAS3/ibd/calling/logs.tar.gz | awk '{print $4}' | sort -n | tail

# Which dataset was used
dataset="gwas3"

# Which samples to include
# Set to "cd", "uc", or "ibd"
assoc="ibd"

# Chunk size in bp
# Chunks will be run in parallel
# Orig, 2M = 22h. run 4 times the dfs, so div blocksize by 4, then 2 again to fit in 12h normal queue
# Max chunk runtime ~8h
chunk_size=200000

# Dir with info on marker positions, used to calculate chunk coords.
QC_SUMMARY_DIR="/lustre/scratch113/projects/crohns/RELEASE/QCsummaries/"

# Which chromosomes to run
#
for ((chr = 1; chr <= 22; chr++)) {

    # Calculate number of chunks based on the maximum position of a snp on the chromosome
    max_genome_pos="$(zcat $QC_SUMMARY_DIR/$chr.txt.gz | cut -d':' -f 2 | cut -d'_' -f 1 | sort -n | tail -n 1)"
    n_chunks=$(echo $(printf %.0f $(echo $max_genome_pos "/ $chunk_size" | bc)) "+1" | bc)

    date && echo "Chrom $chr will run in $n_chunks chunks of $chunk_size bp."

    # Make logging dir
    mkdir -p "$OUT_DIR/logs/$dataset/$assoc/$chr/"

    # bsub a job array over the chunks
    #
    # -Q to resubmit jobs failing due to std::bad_alloc
    # Seems to occur haphazardy, with different chunks failing each time. Possible node failures.
    # TODO Currently unsure if different -R requests change the failing chunks.
    #
    bsub \
        -G team152 -q normal \
        -Q "255" \
        -R "select[mem>2000] rusage[mem=2000]" -M 2000 \
        -J "$dataset.$assoc.$chr[1-$n_chunks]" \
        -o "$OUT_DIR/logs/$dataset/$assoc/$chr/jobid_%J.i_%I.bsub_o.log" \
        -e "$OUT_DIR/logs/$dataset/$assoc/$chr/jobid_%J.i_%I.bsub_e.log" \
        "bash 1.1_run_snptest_chunk.sh $OUT_DIR $dataset $assoc $chr $chunk_size"

}

