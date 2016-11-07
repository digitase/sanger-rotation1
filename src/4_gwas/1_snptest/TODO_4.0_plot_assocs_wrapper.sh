
# A wrapper script to bsub plot assocs.

# Template from snptest run:
# bsub \
    # -G team152 -q normal \
    # -R "select[mem>10000] rusage[mem=10000]" -M 2000 \
    # -J "$dataset.$assoc.$chr[1-$n_chunks]" \
    # -o "$OUT_DIR/logs/$dataset/$assoc/$chr/jobid_%J.i_%I.bsub_o.log" \
    # -e "$OUT_DIR/logs/$dataset/$assoc/$chr/jobid_%J.i_%I.bsub_e.log" \
    # "bash 1.1_run_snptest_chunk.sh $OUT_DIR $dataset $assoc $chr $chunk_size"

