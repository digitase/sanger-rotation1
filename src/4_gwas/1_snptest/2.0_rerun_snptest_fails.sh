#!/usr/local/bin/bash
#
# Re-run GWAS with snptest
# 
# Runs the analysis on chunks of each chromosome that failed (non zero exit code)
# e.g. timed out, killed etc.
#

OUT_DIR="/nfs/users/nfs_b/bb9/workspace/rotation1/crohns_workspace/4_gwas/1_snptest/"
mkdir -p "$OUT_DIR"

# Which dataset was used
dataset="gwas3"
# Which samples to include
# Set to "cd", "uc", or "ibd"
# assoc="cd"
# assoc="uc"
assoc="ibd"
# Chunk size in bp
chunk_size=200000
# Dir with info on marker positions, used to calculate chunk coords.
QC_SUMMARY_DIR="/lustre/scratch113/projects/crohns/RELEASE/QCsummaries/"

# Get all logs with failures
echo Scanning log files in $OUT_DIR/logs/$dataset/$assoc for failed chunks...
failed_logs=($(grep -lr 'Exited with exit code' --include=*.bsub_o.log $OUT_DIR/logs/$dataset/$assoc))

# Rerun failed chunks for each chr
for ((chr = 1; chr <= 22; chr++)) {
    echo Processing re-runs for chr $chr...
    # Accumulate id of failed chunks for the chr
    failed_chunks=()
    for log in ${failed_logs[@]}; do
        if [[ "$log" =~ "/$dataset/$assoc/$chr/jobid_" ]]; then
            chunk_i=$(echo $(basename $log) | cut -d'.' -f 2)
            chunk_i=${chunk_i#i_}
            failed_chunks+=("$chunk_i")
            # Change file extension of logs to indicate rerun has been scheduled
            # This ensures a re-run of this re-run script will not pick up failed chunks that have already been re-run.
            mv "$log" "${log%.log}.reran.log"
            mv "${log%.bsub_o.log}.bsub_e.log" "${log%.bsub_o.log}.bsub_e.reran.log"
        fi
    done
    # Sort chunk numbers
    IFS=$'\n' failed_chunks=($(sort -n <<<"${failed_chunks[*]}"))
    unset IFS
    # echo ${failed_chunks[@]}
    # echo ${#failed_chunks[@]}
    if [[ ${#failed_chunks[@]} > 0 ]]; then
        # Change delimiter to comma
        failed_chunks=$(echo ${failed_chunks[@]} | tr ' ' ',' | sed 's/,$//')
        date && echo Re-running chunks: $failed_chunks 
        # bsub a job array over the failed chunks
        bsub \
            -G team152 -q long \
            -Q "255" \
            -R "select[mem>2000] rusage[mem=2000]" -M 2000 \
            -J "$dataset.$assoc.$chr[$failed_chunks]" \
            -o "$OUT_DIR/logs/$dataset/$assoc/$chr/jobid_%J.i_%I.bsub_o.log" \
            -e "$OUT_DIR/logs/$dataset/$assoc/$chr/jobid_%J.i_%I.bsub_e.log" \
            "bash 1.1_run_snptest_chunk.sh $OUT_DIR $dataset $assoc $chr $chunk_size"
    else
        echo No chunks to re-run.
    fi
}

