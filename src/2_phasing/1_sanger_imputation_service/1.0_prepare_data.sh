
# For VCF, we require:
    # Valid VCF
    # All alleles on the forward strand
    # Coordinates are on GRCh37
    # REF allele matches GRCh37. See the resources for help checking and fixing the REF allele.
    # A single VCF file, not one file per-chromosome
    # Records are sorted by genomic position (chromosomal order is not important)
    # Chromosome names should be 1, 2, 3, etc… not chr1, chr2, chr3, etc… They
        # should match the names in this reference index file. Some programs will
        # represent X as 23, Y as 24, etc…. Please remove or replace these names. See
        # the resources for help renaming chromosomes in a VCF.  If not requesting
        # pre-phasing, then all sites and samples should be phased with no missing
        # data.

# Set paths and vars
PLINK_DIR="/software/hgi/pkglocal/plink-1.90b3w/bin/"
TABIX_DIR="/software/hgi/pkglocal/tabix-git-1ae158a/bin/"
BEDTOOLS_DIR="/software/hgi/pkglocal/bedtools-2.22.0/bin/"

export BCFTOOLS_PLUGINS="/nfs/users/nfs_b/bb9/packages/bcftools/plugins/"

PATH="$PLINK_DIR:$TABIX_DIR:$BEDTOOLS_DIR:$PATH"

IN_DATA_DIR="/nfs/users/nfs_b/bb9/workspace/rotation1/crohns_workspace/1_qc/"
IN_DATA_PREFIX="coreex_gaibdc_usgwas_raw.qc6.maf_0.001"

OUT_DIR="/nfs/users/nfs_b/bb9/workspace/rotation1/crohns_workspace/2_phasing/1_sanger_imputation_service/"
mkdir -p "$OUT_DIR"

# Only chr 20 for now
CHR=20

# Convert binary PLINK to vcf
plink \
    --bfile "$IN_DATA_DIR/$IN_DATA_PREFIX" \
    --chr "$CHR" \
    --recode-vcf \
    --out "$OUT_DIR/$IN_DATA_PREFIX"

#
# Check REF allele matches GRCh37 base and coordinate
#

# Alleles are already on the fwd strand.

# Fix allele ordering.
# plink automatically sets the major (common) allele as the reference allele for each population when generating the bim (map) files.
# Hence the true reference allele has been lost.
# Force the A2 allele to match the ref allele in GRCh37.

# Download reference
if [ ! -f "$OUT_DIR/human_g1k_v37.fasta" ]; then
    wget -O "$OUT_DIR/human_g1k_v37.fasta.gz" ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
    gunzip "$OUT_DIR/human_g1k_v37.fasta.gz"
fi

# Extract ref bases at the coords of the markers from the reference fasta.
# Extract marker ids from the vcf.
# Merge.
paste \
    <(bedtools getfasta -fi "$OUT_DIR/human_g1k_v37.fasta" -bed "$OUT_DIR/$IN_DATA_PREFIX.vcf" -fo - | grep -v '>') \
    <(awk '/^#/ {next} {print $3}' "$OUT_DIR/$IN_DATA_PREFIX.vcf") \
        > "$OUT_DIR/$IN_DATA_PREFIX.ref_alleles.txt"

# Regenerate vcf with correct allele ordering
plink \
    --bfile "$IN_DATA_DIR/$IN_DATA_PREFIX" \
    --chr "$CHR" \
    --a2-allele "$OUT_DIR/$IN_DATA_PREFIX.ref_alleles.txt" 1 2 \
    --recode-vcf \
    --out "$OUT_DIR/$IN_DATA_PREFIX.alleles_ordered"

# Remove indels
awk '/^#/ {print; next} {if ($4 !~ /[I|D]/) {print}}' "$OUT_DIR/$IN_DATA_PREFIX.alleles_ordered.vcf" > "$OUT_DIR/$IN_DATA_PREFIX.alleles_ordered.no_indel.vcf"

# Check for remaining mismatches
# which occur when the reference base does not match either of the vcf alleles.
bcftools +fixref "$OUT_DIR/$IN_DATA_PREFIX.alleles_ordered.no_indel.vcf" -- -f "$OUT_DIR/human_g1k_v37.fasta"

# Get pos of mismatches
bcftools norm --check-ref w -f "$OUT_DIR/human_g1k_v37.fasta" "$OUT_DIR/$IN_DATA_PREFIX.alleles_ordered.no_indel.vcf" -o /dev/null \
    2> >(awk '/^REF_MISMATCH/ {print "POS=="$3}' - > "$OUT_DIR/$IN_DATA_PREFIX.alleles_ordered.no_indel.ref_mismatches_pos.txt")

# Exclude the snps at those positions
bcftools view \
    --exclude \
        "$(python -c "import sys; print(' | '.join(l.strip() for l in sys.stdin))" < $OUT_DIR/$IN_DATA_PREFIX.alleles_ordered.no_indel.ref_mismatches_pos.txt)" \
    "$OUT_DIR/$IN_DATA_PREFIX.alleles_ordered.no_indel.vcf" \
        > "$OUT_DIR/$IN_DATA_PREFIX.alleles_ordered.no_indel.no_ref_mismatch.vcf"

# The resulting VCF file is uploaded to the online impuation service.
# Check job status here:
# https://imputation.sanger.ac.uk/?status=1&rid=a279d523d8b71059c38494fb118c8bfd

# Results were downloaded to /lustre/scratch113/projects/crohns/bb9/2_phasing/1_sanger_imputation_service/results_2016-10-25/

