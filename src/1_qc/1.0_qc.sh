#!/usr/local/bin/bash

#
# Following the QC protocol from Anderson et al 2010.
# Modifications suggested by Katie de Lange.
#

# Set paths and vars
PLINK_DIR="/software/hgi/pkglocal/plink-1.90b3w/bin/"
R_DIR="/software/R-3.3.0/bin/"

PATH="$PLINK_DIR:$R_DIR:$PATH"

IN_DATA_DIR="$HOME/workspace/rotation1/data/gwas3/"
IN_DATA_PREFIX="coreex_gaibdc_usgwas_raw"

OUT_DIR="$HOME/workspace/rotation1/crohns_workspace/1_qc/"
mkdir -p "$OUT_DIR"

# Symlink in raw opticalled data to output dir
ln -s "$IN_DATA_DIR/$IN_DATA_PREFIX.bed" "$OUT_DIR/$IN_DATA_PREFIX.bed"
ln -s "$IN_DATA_DIR/$IN_DATA_PREFIX.bim" "$OUT_DIR/$IN_DATA_PREFIX.bim"
ln -s "$IN_DATA_DIR/$IN_DATA_PREFIX.fam" "$OUT_DIR/$IN_DATA_PREFIX.fam"

#
# Part 1.
# Remove markers based on missingness.
#
# Turner, S., Armstrong, L. L., Bradford, Y., Carlson, C. S., Crawford, D. C., Crenshaw, A. T., … Ritchie, M. D. (2011). Quality control procedures for genome-wide association studies. Current Protocols in Human Genetics / Editorial Board, Jonathan L. Haines ... [et Al.], Chapter 1, Unit1.19. https://doi.org/10.1002/0471142905.hg0119s68
# We recommend removing poor quality SNPs before running the sample genotyping
# efficiency check discussed above, so that fewer samples will be dropped from
# the analysis simply because they were genotyped with SNP assays that had poor
# performance. 
#
# Markers can also be recovered via imputation, samples cannot.
#
MAX_MARKER_MISSINGNESS="0.05"

plink \
    --bfile "$OUT_DIR/$IN_DATA_PREFIX" \
    --geno "$MAX_MARKER_MISSINGNESS" \
    --make-bed \
    --out "$OUT_DIR/$IN_DATA_PREFIX.qc1" 

#
# Part 2.
# Remove samples based on missingness.
# also Identify samples with excess heterozygosity.
#
# The genotype failure rate and heterozygosity rate per individual are both measures of DNA sample quality.
# These two filters are independent, as they only consider info relevant to each a sample.
#

# Get sample missingness and heterozygosity rate
plink \
    --bfile "$OUT_DIR/$IN_DATA_PREFIX.qc1" \
    --missing \
    --het \
    --out "$OUT_DIR/$IN_DATA_PREFIX.qc2" 

# Plot heterozygosity, missingness and optimise thresholds
MAX_SAMPLE_MISSINGNESS="0.01"
MAX_HET_RATE_SD_THRESH=3

Rscript "1.1_plot_imiss_vs_het.R" \
    "$OUT_DIR/$IN_DATA_PREFIX.qc2.imiss" \
    "$OUT_DIR/$IN_DATA_PREFIX.qc2.het" \
    "$MAX_SAMPLE_MISSINGNESS" \
    "$MAX_HET_RATE_SD_THRESH" \
    "$OUT_DIR/$IN_DATA_PREFIX.qc2.pdf" \
    "$OUT_DIR/$IN_DATA_PREFIX.qc2.sample_fail_imiss.txt" \
    "$OUT_DIR/$IN_DATA_PREFIX.qc2.sample_fail_het.txt"

# Remove samples failing missingness check
plink \
    --bfile "$OUT_DIR/$IN_DATA_PREFIX.qc1" \
    --remove "$OUT_DIR/$IN_DATA_PREFIX.qc2.sample_fail_imiss.txt" \
    --make-bed \
    --out "$OUT_DIR/$IN_DATA_PREFIX.qc3" 

#
# Part 3.
#
# Other sample QC filters:
# Check sex concordancy, sample ancestry, relatedness.
# Apply filters as to minimise the number of samples lost.
#

#
# --check-sex
# Identification of individuals with discordant sex information
# Calculate the mean homozygosity rate across X chromosome markers for each individual in the study.
# When the homozygosity rate is more than 0.2 but less than 0.8 the genotype data is inconclusive regarding the sex of an individual and these are marked in column 4 with a 0.
#
plink \
    --bfile "$OUT_DIR/$IN_DATA_PREFIX.qc3" \
    --check-sex \
    --out "$OUT_DIR/$IN_DATA_PREFIX.qc3.check_sex" 

# Identification of duplicated or related individuals

# First prune to sites with >5% freq to save time
# Rare variation is generally not required to determine relatedness.
plink \
    --bfile "$OUT_DIR/$IN_DATA_PREFIX.qc3" \
    --maf "0.05" \
    --make-bed \
    --out "$OUT_DIR/$IN_DATA_PREFIX.qc3.maf_0.05"

# Use KING for relationship inference 
# up to and including third degree relatives
# Runtime: approx 2-5h
# Extra flag for farm use:
    # -G team152 \

bsub \
    -K \
    -J "$IN_DATA_PREFIX.qc3.maf_0.05.king" \
    -R "select[mem>4000] rusage[mem=4000]" \
    -M 4000 \
    -o "$OUT_DIR/$IN_DATA_PREFIX.qc3.maf_0.05.king.log" \
        "king \
            -b $OUT_DIR/$IN_DATA_PREFIX.qc3.maf_0.05.bed \
            --kinship --ibs \
            --related --degree 3 \
            --prefix $OUT_DIR/$IN_DATA_PREFIX.qc3.maf_0.05.king"

# Identification of individuals of divergent ancestry using PCA

# Merge genotypes to HapMap Phase III (HapMap3) data from four ethnic populations.
HAPMAP3R2_MARKERS="/lustre/scratch113/teams/barrett/coreex_gaibdc/refs/hapmap3r2_CEU.CHB.JPT.YRI.no-at-cg-snps.txt"
HAPMAP3R2_PREFIX="/lustre/scratch113/teams/barrett/coreex_gaibdc/refs/hapmap3r2_CEU.CHB.JPT.YRI.founders.no-at-cg-snps"

# Extract hapmap markers from our dataset
plink \
    --bfile "$OUT_DIR/$IN_DATA_PREFIX.qc3" \
    --extract "$HAPMAP3R2_MARKERS" \
    --make-bed \
    --out "$OUT_DIR/$IN_DATA_PREFIX.qc3.hapmap_markers"

# Prune SNPs so that no pair of SNPs (within a given number of base pairs) has an R2 value greater than a given threshold
plink \
    --bfile "$OUT_DIR/$IN_DATA_PREFIX.qc3" \
    --exclude "/lustre/scratch113/teams/barrett/coreex_gaibdc/refs/high-LD-regions.txt" \
    --range --indep-pairwise 50 5 "0.2" \
    --out "$OUT_DIR/$IN_DATA_PREFIX.qc3.prune_markers_by_ld"

# Merge our samples and hapmap samples, keeping only the LD pruned markers.
plink \
    --bfile "$OUT_DIR/$IN_DATA_PREFIX.qc3.hapmap_markers" \
    --bmerge "$HAPMAP3R2_PREFIX.bed" "$HAPMAP3R2_PREFIX.bim" "$HAPMAP3R2_PREFIX.fam" \
    --allow-no-sex \
    --extract "$OUT_DIR/$IN_DATA_PREFIX.qc3.prune_markers_by_ld.prune.in" \
    --make-bed \
    --out "$OUT_DIR/$IN_DATA_PREFIX.qc3.pruned.hapmap_merged"

# From the list of snps with strand problems caused by >2 alleles (listed in .missnp output of the merge), 
# change their alleles by flipping their strands.
plink \
    --bfile "$OUT_DIR/$IN_DATA_PREFIX.qc3" \
    --extract "$HAPMAP3R2_MARKERS" \
    --flip "$OUT_DIR/$IN_DATA_PREFIX.qc3.pruned.hapmap_merged-merge.missnp" \
    --make-bed \
    --out "$OUT_DIR/$IN_DATA_PREFIX.qc3.hapmap_markers.flipped"

# Repeat the merge with problematic markers flipped 
plink \
    --bfile "$OUT_DIR/$IN_DATA_PREFIX.qc3.hapmap_markers.flipped" \
    --bmerge "$HAPMAP3R2_PREFIX.bed" "$HAPMAP3R2_PREFIX.bim" "$HAPMAP3R2_PREFIX.fam" \
    --allow-no-sex \
    --extract "$OUT_DIR/$IN_DATA_PREFIX.qc3.prune_markers_by_ld.prune.in" \
    --make-bed \
    --out "$OUT_DIR/$IN_DATA_PREFIX.qc3.pruned.hapmap_merged.flipped"

# Shorten family/sample ids in fam file, as smartpca has a max length limit
    # -v map="$OUT_DIR/$IN_DATA_PREFIX.qc3.pruned.hapmap_merged.flipped.id_to_short_id.mappings" \
    # '{print "ID" NR FS $1 >> map; print "ID" NR FS "ID" NR FS $3 FS $4 FS $5 FS $6 >> fam}' \
awk -v fam="$OUT_DIR/$IN_DATA_PREFIX.qc3.pruned.hapmap_merged.flipped.short_id.fam" \
    '{print "ID" NR FS "ID" NR FS $3 FS $4 FS $5 FS $6 >> fam}' \
    "$OUT_DIR/$IN_DATA_PREFIX.qc3.pruned.hapmap_merged.flipped.fam" 

# Create symlinks to data with different extensions
# smartpca requires specific extensions to recognise input file formats
ln -s "$OUT_DIR/$IN_DATA_PREFIX.qc3.pruned.hapmap_merged.flipped.bim" "$OUT_DIR/$IN_DATA_PREFIX.qc3.pruned.hapmap_merged.flipped.pedsnp"
ln -s "$OUT_DIR/$IN_DATA_PREFIX.qc3.pruned.hapmap_merged.flipped.short_id.fam" "$OUT_DIR/$IN_DATA_PREFIX.qc3.pruned.hapmap_merged.flipped.pedind"

#
# Run smartpca perl script version. smartpca must be in $PATH
#
# -i, -a, -b specify input genotype, snp and indiv files
# -o, -e, -p, -l specify output file names
# -k    number of eigenvectors to output
# -t    number of principal components along which to remove outliers during each outlier removal iteration
# -w    If wishing to infer eigenvectors using only individuals from a 
#       subset of populations, and then project individuals from all populations 
#       onto those eigenvectors, this input file contains a list of population names,
#       one population name per line, which will be used to infer eigenvectors.  
#       It is assumed that the population of each individual is specified in the 
#       indiv file.
#       In this case, pops 3, 4, 5, 6 are the populations in the last column of the hapmap indiv file.
perl "$HOME/packages/eigensoft/bin/smartpca.perl" \
    -i "$OUT_DIR/$IN_DATA_PREFIX.qc3.pruned.hapmap_merged.flipped.bed" \
    -a "$OUT_DIR/$IN_DATA_PREFIX.qc3.pruned.hapmap_merged.flipped.pedsnp" \
    -b "$OUT_DIR/$IN_DATA_PREFIX.qc3.pruned.hapmap_merged.flipped.pedind" \
    -o "$OUT_DIR/$IN_DATA_PREFIX.qc3.pruned.hapmap_merged.flipped.pca" \
    -e "$OUT_DIR/$IN_DATA_PREFIX.qc3.pruned.hapmap_merged.flipped.eval" \
    -p "$OUT_DIR/$IN_DATA_PREFIX.qc3.pruned.hapmap_merged.flipped.plot" \
    -l "$OUT_DIR/$IN_DATA_PREFIX.qc3.pruned.hapmap_merged.flipped.log" \
    -k 2 -t 2 \
    -w "/lustre/scratch113/teams/barrett/coreex_gaibdc/refs/pca-populations.txt"

# Choose a PC2 threshold to exclude samples of non-European ancestry
MIN_PC2_THRESH="0.066"
Rscript "1.2_plot-pca-results_alpha.R" "$OUT_DIR/$IN_DATA_PREFIX.qc3.pruned.hapmap_merged.flipped.pca" "$MIN_PC2_THRESH"

# Remove samples failing any QC checks:
#   Samples with excess missingness have already been removed.
#   Remove samples with excess heterozygosity, 
#   discordant sex (accounting for known mishaps i.e. plate swap),
#   non-European ancestry.
#   For each pair of samples that are first degree relatives (kinship coefficient > 0.177) or closer,
#   where both samples pass other QC checks, 
#   remove the individual with the lower call rate.
#   See http://people.virginia.edu/~wc9c/KING/manual.html for info on thresholds on kinship.
#
MAX_KINSHIP_THRESH="0.177"

# Mark for removal
Rscript "1.3_identify_failed_samples.R" \
    "$OUT_DIR/$IN_DATA_PREFIX.qc1.fam" \
    "$OUT_DIR/$IN_DATA_PREFIX.qc2.imiss" \
    "$OUT_DIR/$IN_DATA_PREFIX.qc2.sample_fail_imiss.txt" \
    "$OUT_DIR/$IN_DATA_PREFIX.qc2.sample_fail_het.txt" \
    "$OUT_DIR/$IN_DATA_PREFIX.qc3.check_sex.sexcheck" \
    "$OUT_DIR/$IN_DATA_PREFIX.qc3.pruned.hapmap_merged.flipped.pca.evec" \
    "$OUT_DIR/$IN_DATA_PREFIX.qc3.pruned.hapmap_merged.flipped.fam" \
    "$OUT_DIR/$IN_DATA_PREFIX.qc3.maf_0.05.king.ibs0" \
    "$MIN_PC2_THRESH" \
    "$MAX_KINSHIP_THRESH" \
    "$OUT_DIR/$IN_DATA_PREFIX.qc3.summary.txt" \
    "$OUT_DIR/$IN_DATA_PREFIX.qc3.summary_table.txt" \
    "$OUT_DIR/$IN_DATA_PREFIX.qc3.sample_fail_any.txt"

# Remove
plink \
    --bfile "$OUT_DIR/$IN_DATA_PREFIX.qc3" \
    --remove "$OUT_DIR/$IN_DATA_PREFIX.qc3.sample_fail_any.txt" \
    --make-bed \
    --out "$OUT_DIR/$IN_DATA_PREFIX.qc4" 

#
# Part 4.
# Marker QC
#
# Recheck missingness, difference in case vs. control missingness, departure from HWE, batch effects.
#

# Recheck marker missingness after sample QC
plink \
    --bfile "$OUT_DIR/$IN_DATA_PREFIX.qc4" \
    --missing \
    --out "$OUT_DIR/$IN_DATA_PREFIX.qc4.missing" 

# Verify that the original missingness threshold is still sensible
MAX_MARKER_MISSINGNESS_1="$MAX_MARKER_MISSINGNESS"
Rscript "1.4_lmiss-hist.R" "$OUT_DIR/$IN_DATA_PREFIX.qc4.missing" "$MAX_MARKER_MISSINGNESS_1"

# Find markers with a significantly different in missingness rate (genotype call rate) between cases and controls
plink \
    --bfile "$OUT_DIR/$IN_DATA_PREFIX.qc4" \
    --test-missing \
    --out "$OUT_DIR/$IN_DATA_PREFIX.qc4.diffmiss"

DIFFMISS_PVAL_THRESH="0.00001"

# Write out failures
awk -v thresh="$DIFFMISS_PVAL_THRESH" '($5 < thresh) {print $2}' \
    "$OUT_DIR/$IN_DATA_PREFIX.qc4.diffmiss.missing" \
        > "$OUT_DIR/$IN_DATA_PREFIX.qc4.diffmiss.marker_fail_diffmiss.txt"

#
# Filter markers
#
# 1. Remove markers that depart significantly from Hardy-Weinburg equilibrium in controls.
#    Most GWA studies exclude markers that show extensive devia tion from HWE because this can be indicative of a genotyping or genotype-calling error.
#    The "midp" option applies the mid-p adjustment described in Graffelman J,
#        Moreno V (2013) The mid p-value in exact tests for Hardy-Weinberg
#        equilibrium. The mid-p adjustment tends to bring the null rejection rate
#        in line with the nominal p-value, and also reduces the filter's tendency
#        to favor retention of variants with missing data. 
#
# 2. Refilter for marker missingness.
# 
# 3. --exclude markers significantly different in missingness rate (genotype call rate) between cases and controls
#
HWE_PVAL_THRESH="0.00001"

plink \
    --bfile "$OUT_DIR/$IN_DATA_PREFIX.qc4" \
    --geno "$MAX_MARKER_MISSINGNESS_1" \
    --hwe "$HWE_PVAL_THRESH" midp \
    --exclude "$OUT_DIR/$IN_DATA_PREFIX.qc4.diffmiss.marker_fail_diffmiss.txt" \
    --make-bed \
    --out "$OUT_DIR/$IN_DATA_PREFIX.qc5"

#
# 4. --exclude markers with the batch effects identified by Yang's PC analysis. This involved:
#     - Computing within-sample PCs using common variants (MAF > 1%)
#     - A clear case outlier set was seen. These were all part of batch3
#     - Use PC1 to split cases into outliers and non-outliers. Perform an association
#     	test between these two groups, and identify 'significant' sites (P<1e-5)
#     - Because this problematic batch is all cases, these will generate false hits if
#    	left in. Thus all these sites (429) need to be removed.
#

plink \
    --bfile "$OUT_DIR/$IN_DATA_PREFIX.qc5" \
    --exclude "/lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/pcs/batch/batch-fail-1e-5.txt" \
    --make-bed \
    --out "$OUT_DIR/$IN_DATA_PREFIX.qc6"

#
# Part 5. 
# Before imputation, threshold on MAF
# 
# Rare variants are difficult to genotype, and power to detect associations is low.
#

IMPUTE_MAF_THRESH="0.001"
plink \
    --bfile "$OUT_DIR/$IN_DATA_PREFIX.qc6" \
    --maf "$IMPUTE_MAF_THRESH" \
    --make-bed \
    --out "$OUT_DIR/$IN_DATA_PREFIX.qc6.maf_$IMPUTE_MAF_THRESH"

