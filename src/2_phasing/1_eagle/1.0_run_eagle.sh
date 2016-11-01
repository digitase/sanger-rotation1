#!/usr/local/bin/bash
#
# BSUB -J 1.0_run_eagle
# BSUB -n 16
# BSUB -R "select[mem>4000] rusage[mem=4000]" -M 4000
# BSUB -cwd .output
# BSUB -o 1.0_run_eagle.bsub_%J_%I_o.log
# BSUB -e 1.0_run_eagle.bsub_%J_%I_e.log
#
# Phase data with eagle2
#

# Only process chr20.
CHR=20

# Match this in the BSUB header
N_THREADS=16

IN_DATA_DIR="/nfs/users/nfs_b/bb9/workspace/rotation1/crohns_workspace/2_phasing/1_sanger_imputation_service/"
IN_DATA_PREFIX="coreex_gaibdc_usgwas_raw.qc6.maf_0.001.alleles_ordered.no_indel.no_ref_mismatch"

OUT_DIR="/nfs/users/nfs_b/bb9/workspace/rotation1/crohns_workspace/2_phasing/1_eagle/"
mkdir -p "$OUT_DIR"

#
# Use eagle2 to phase with a reference panel
#

# Convert 1000genomes reference panel data (IMPUTE format) to vcf
# REF_HAPLO_FILE="/lustre/scratch113/projects/crohns/2013Aug07/imputation/reference/ALL.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nosing/ALL.chr$CHR.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nosing.haplotypes.gz"
# REF_LEGEND_FILE="/lustre/scratch113/projects/crohns/2013Aug07/imputation/reference/ALL.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nosing/ALL.chr$CHR.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nosing.legend.gz"
# REF_PREFIX="$(basename $REF_HAPLO_FILE .haplotypes.gz)"

# gunzip -c "$REF_HAPLO_FILE" > "$OUT_DIR/$REF_PREFIX.haplotypes"
# gunzip -c "$REF_LEGEND_FILE" > "$OUT_DIR/$REF_PREFIX.legend"

# qctool \
    # -g "$OUT_DIR/$REF_PREFIX.haplotypes" -filetype impute_haplotypes \
    # -assume-chromosome "$CHR" \
    # -og "$OUT_DIR/$REF_PREFIX.vcf"

# bcftools convert -Ob - > "$OUT_DIR/$REF_PREFIX.bcf.gz"

# ABORT this pipeline implementation:
#
# NOTE: If the data set you wish to phase contains more than twice as many
# samples as the largest reference panel available to you, then using a reference
# panel is unlikely to give much of a boost in phasing accuracy. 
#    The 1000 genomes panel: ~ 2K samples
#    Our dataset: ~ 18K samples
#
# Therefore...

#
# Use eagle2 to perform phasing without a reference panel
#

# Compress vcf -> bcf.gz
ln -s "$IN_DATA_DIR/$IN_DATA_PREFIX.vcf" "$OUT_DIR/$IN_DATA_PREFIX.vcf"
bcftools convert -Ob <"$OUT_DIR/$IN_DATA_PREFIX.vcf" >"$OUT_DIR/$IN_DATA_PREFIX.bcf.gz"

# TODO can change output format using --vcfOutFormat
eagle \
    --vcf "$OUT_DIR/$IN_DATA_PREFIX.bcf.gz" \
    --geneticMapFile="/nfs/users/nfs_b/bb9/packages/Eagle_v2.3/tables/genetic_map_hg19_withX.txt.gz" \
    --chrom "$CHR" \
    --numThreads "$N_THREADS" \
    --outPrefix="$OUT_DIR/$IN_DATA_PREFIX.phased" \
        2>&1 | tee "$OUT_DIR/$IN_DATA_PREFIX.phased.log"

