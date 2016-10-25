
CHR=20
IN_PLINK_DIR="/nfs/users/nfs_b/bb9/workspace/rotation1/crohns_workspace/1_qc/"
IN_PLINK_PREFIX="coreex_gaibdc_usgwas_raw.qc6.maf_0.001"

OUT_DIR="/lustre/scratch113/projects/crohns/bb9/2_phasing/2_pbwt/"
mkdir -p "$OUT_DIR"

# Correct the strand of the data
# TODO Does exactly the opposite (horribly broken).
# wget -O "$OUT_DIR/HumanCoreExome-12v1-1_B-b37.strand.RefAlt.zip" "http://www.well.ox.ac.uk/~wrayner/strand/RefAlt/HumanCoreExome-12v1-1_B-b37.strand.RefAlt.zip"
# unzip "$OUT_DIR/HumanCoreExome-12v1-1_B-b37.strand.RefAlt.zip" -d "$OUT_DIR"

# plink \
    # --bfile "$IN_PLINK_DIR/$IN_PLINK_PREFIX" \
    # --chr "$CHR" \
    # --reference-allele "$OUT_DIR/HumanCoreExome-12v1-1_B-b37.strand.RefAlt" \
    # --recode-vcf \
    # --out "$OUT_DIR/$IN_PLINK_PREFIX.ref_allele_forced"
# export BCFTOOLS_PLUGINS="/nfs/users/nfs_b/bb9/packages/bcftools/plugins/"
# bcftools +fixref "$OUT_DIR/$IN_PLINK_PREFIX.ref_allele_forced.vcf" -- -f "../1_sanger_imputation_service/.output/human_g1k_v37.fasta"

CHROM:POS_REF_ALT


# Convert from output gzipped Oxford HAPS/SAMPLE format (used by SHAPEIT2) to vcf
bcftools convert \
    -Oz \
    --hapsample2vcf "$OUT_DIR/$IN_PLINK_PREFIX.phased.haps.gz","$OUT_DIR/$IN_PLINK_PREFIX.phased.sample" \
        > "$OUT_DIR/$IN_PLINK_PREFIX.phased.vcf.gz"

pbwt -checkpoint 10000 -readVcfGT data.vcf -writeAll data
pbwt 
-log logfile
-readVcfGT - 
-referenceFasta resources/refs/human_g1k_v37.fasta 
-referenceImpute resources/refs/imputation/hrc.r1.1/pbwt/HRC.r1-1.GRCh37.chr20.shapeit3.mac5.aa.genotypes
-writeBcfGz adfsdfasdf.part

