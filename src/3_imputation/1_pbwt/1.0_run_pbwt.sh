#!/usr/local/bin/bash
#
# BSUB -J 1.0_run_pbwt
# BSUB -G team152
# BSUB -R "select[mem>10000] rusage[mem=10000]" -M 10000
# BSUB -cwd .output
# BSUB -o 1.0_run_pbwt.bsub_%J_%I_o.log
# BSUB -e 1.0_run_pbwt.bsub_%J_%I_e.log
#
# Impute haplotypes
#

# Which chromosome to impute
CHR=20

PHASED_BCF="/nfs/users/nfs_b/bb9/workspace/rotation1/crohns_workspace/2_phasing/1_eagle/coreex_gaibdc_usgwas_raw.qc6.maf_0.001.alleles_ordered.no_indel.no_ref_mismatch.phased.bcf"

OUT_DIR="/nfs/users/nfs_b/bb9/workspace/rotation1/crohns_workspace/3_imputation/1_pbwt/"
mkdir -p "$OUT_DIR"

#
# Correct the strand of the data
# TODO Currently does exactly the opposite (i.e. horribly broken).
#
# wget -O "$OUT_DIR/HumanCoreExome-12v1-1_B-b37.strand.RefAlt.zip" "http://www.well.ox.ac.uk/~wrayner/strand/RefAlt/HumanCoreExome-12v1-1_B-b37.strand.RefAlt.zip"
# unzip "$OUT_DIR/HumanCoreExome-12v1-1_B-b37.strand.RefAlt.zip" -d "$OUT_DIR"
# IN_PLINK_DIR="/nfs/users/nfs_b/bb9/workspace/rotation1/crohns_workspace/1_qc/"
# IN_PLINK_PREFIX="coreex_gaibdc_usgwas_raw.qc6.maf_0.001"
# plink \
    # --bfile "$IN_PLINK_DIR/$IN_PLINK_PREFIX" \
    # --chr "$CHR" \
    # --reference-allele "$OUT_DIR/HumanCoreExome-12v1-1_B-b37.strand.RefAlt" \
    # --recode-vcf \
    # --out "$OUT_DIR/$IN_PLINK_PREFIX.ref_allele_forced"
# export BCFTOOLS_PLUGINS="/nfs/users/nfs_b/bb9/packages/bcftools/plugins/"
# bcftools +fixref "$OUT_DIR/$IN_PLINK_PREFIX.ref_allele_forced.vcf" -- -f "../1_sanger_imputation_service/.output/human_g1k_v37.fasta"

# Create pbwt imputeRef files for the panel
PANEL_PREFIX="ALL.chr20.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nosing"

gunzip -c "/lustre/scratch113/projects/crohns/2013Aug07/imputation/reference/ALL.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nosing/$PANEL_PREFIX.haplotypes.gz" > "$OUT_DIR/$PANEL_PREFIX.haplotypes"
gunzip -c "/lustre/scratch113/projects/crohns/2013Aug07/imputation/reference/ALL.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nosing/$PANEL_PREFIX.legend.gz" > "$OUT_DIR/$PANEL_PREFIX.legend"

date && echo Building panel...
pbwt \
    -log "$OUT_DIR/$PANEL_PREFIX.imputeRef.log" \
    -readHapLegend \
        "$OUT_DIR/$PANEL_PREFIX.haplotypes" \
        "$OUT_DIR/$PANEL_PREFIX.legend" \
        "$CHR" \
    -writeAll "$OUT_DIR/$PANEL_PREFIX.imputeRef" \
    -writeImputeRef "$OUT_DIR/$PANEL_PREFIX.imputeRef"

# Convert vcf.gz to bcf
# gunzip -c "/nfs/users/nfs_b/bb9/workspace/rotation1/crohns_workspace/2_phasing/1_eagle/coreex_gaibdc_usgwas_raw.qc6.maf_0.001.alleles_ordered.no_indel.no_ref_mismatch.phased.vcf.gz" | \
    # bcftools convert -Ou - > "$OUT_DIR/coreex_gaibdc_usgwas_raw.qc6.maf_0.001.alleles_ordered.no_indel.no_ref_mismatch.phased.bcf"

# Perform imputation
# Runtime on chr20: 2h, 6Gb memory
date && echo Imputing chr $CHR...
pbwt \
    -log "$OUT_DIR/$CHR.imputation.log" \
    -checkpoint 10000 \
    -readVcfGT "$PHASED_BCF" \
    -referenceImpute "$OUT_DIR/$PANEL_PREFIX.imputeRef" \
    -referenceFasta "/nfs/users/nfs_b/bb9/workspace/rotation1/crohns_workspace/2_phasing/0_prep_data/human_g1k_v37.fasta" \
    -writeBcfGz "$OUT_DIR/$CHR.imputed.bcf.gz"

date && echo Finished.

