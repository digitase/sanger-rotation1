#!/usr/local/bin/bash

#
# Following the QC protocol from Anderson et al 2010.
# Modifications suggested by Katie de Lange.
#

# Set paths
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
# 1.
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

plink --noweb \
    --bfile "$OUT_DIR/$IN_DATA_PREFIX" \
    --geno "$MAX_MARKER_MISSINGNESS" \
    --make-bed \
    --out "$OUT_DIR/$IN_DATA_PREFIX.qc1" 

#
# 2.
# Remove samples based on missingness.
# Identify samples with excess heterozygosity.
#
# The genotype failure rate and heterozygosity rate per individual are both measures of DNA sample quality.
# These two filters are independent, as they only consider info relevant to each a sample.
#

# Get sample missingness and heterozygosity rate
plink --noweb \
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
    "$OUT_DIR/$IN_DATA_PREFIX.qc2.sample_fail_missingness_het.txt"

# Remove samples failing missingness check
plink --noweb \
    --bfile "$OUT_DIR/$IN_DATA_PREFIX.qc1" \
    --mind "$MAX_SAMPLE_MISSINGNESS" \
    --make-bed \
    --out "$OUT_DIR/$IN_DATA_PREFIX.qc3" 

#
# 3.
#
# Other sample filters:
# Check sex concordancy, sample ancestry, relatedness.
# Apply filters as to minimise the number of samples lost.
#

#
# --check-sex
# Identification of individuals with discordant sex information
# Calculate the mean homozygosity rate across X chromosome markers for each individual in the study.
# When the homozygosity rate is more than 0.2 but less than 0.8 the genotype data is inconclusive regarding the sex of an individual and these are marked in column 4 with a 0.
#
plink --noweb \
    --bfile "$OUT_DIR/$IN_DATA_PREFIX.qc3" \
    --check-sex \
    --out "$OUT_DIR/$IN_DATA_PREFIX.qc3.check_sex" 

# Identification of duplicated or related individuals

# First prune to sites with >5% freq to save time
plink --noweb \
    --bfile "$OUT_DIR/$IN_DATA_PREFIX.qc3" \
    --maf "0.05" \
    --make-bed \
    --out "$OUT_DIR/$IN_DATA_PREFIX.qc3.maf_0.05"

# Use KING for relationship inference 
# using Identity By State up to and including third degree relatives
# Runtime: approx 2-5h

bsub \
    -J "$IN_DATA_PREFIX.qc3.maf_0.05.king" \
    -R "select[mem>4000] rusage[mem=4000]" \
    -M 4000 \
    -o "$OUT_DIR/$IN_DATA_PREFIX.qc3.maf_0.05.king.log" \
    -G team152 \
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
plink --noweb \
    --bfile "$OUT_DIR/$IN_DATA_PREFIX.qc3" \
    --extract "$HAPMAP3R2_MARKERS" \
    --make-bed \
    --out "$OUT_DIR/$IN_DATA_PREFIX.qc3.hapmap_markers"

# Prune SNPs so that no pair of SNPs (within a given number of base pairs) has an R2 value greater than a given threshold
plink --noweb \
    --bfile "$OUT_DIR/$IN_DATA_PREFIX.qc3" \
    --exclude "/lustre/scratch113/teams/barrett/coreex_gaibdc/refs/high-LD-regions.txt" \
    --range --indep-pairwise 50 5 "0.2" \
    --out "$OUT_DIR/$IN_DATA_PREFIX.qc3.prune_markers_by_ld"

# Merge our samples and hapmap samples, keeping only the LD pruned markers.
plink --noweb \
    --bfile "$OUT_DIR/$IN_DATA_PREFIX.qc3.hapmap_markers" \
    --bmerge "$HAPMAP3R2_PREFIX.bed" "$HAPMAP3R2_PREFIX.bim" "$HAPMAP3R2_PREFIX.fam" \
    --allow-no-sex \
    --extract "$OUT_DIR/$IN_DATA_PREFIX.qc3.prune_markers_by_ld.prune.in" \
    --make-bed \
    --out "$OUT_DIR/$IN_DATA_PREFIX.qc3.pruned.hapmap_merged"

# From the list of snps with strand problems caused by >2 alleles (listed in .missnp output of the merge), 
# change their alleles by flipping their strands.
plink --noweb \
    --bfile "$OUT_DIR/$IN_DATA_PREFIX.qc3" \
    --extract "$HAPMAP3R2_MARKERS" \
    --flip "$OUT_DIR/$IN_DATA_PREFIX.qc3.pruned.hapmap_merged-merge.missnp" \
    --make-bed \
    --out "$OUT_DIR/$IN_DATA_PREFIX.qc3.hapmap_markers.flipped"

# Repeat the merge with problematic markers flipped 
plink --noweb \
    --bfile "$OUT_DIR/$IN_DATA_PREFIX.qc3.hapmap_markers.flipped" \
    --bmerge "$HAPMAP3R2_PREFIX.bed" "$HAPMAP3R2_PREFIX.bim" "$HAPMAP3R2_PREFIX.fam" \
    --allow-no-sex \
    --extract "$OUT_DIR/$IN_DATA_PREFIX.qc3.prune_markers_by_ld.prune.in" \
    --make-bed \
    --out "$OUT_DIR/$IN_DATA_PREFIX.qc3.pruned.hapmap_merged.flipped"

# Shorten family/sample ids in fam file, as smartpca has a max length limit
awk \
    -v fam="$OUT_DIR/$IN_DATA_PREFIX.qc3.pruned.hapmap_merged.flipped.short_id.fam" \
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
# -o, -p, -e, -l specify output file names
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
    -p "$OUT_DIR/$IN_DATA_PREFIX.qc3.pruned.hapmap_merged.flipped.plot" \
    -e "$OUT_DIR/$IN_DATA_PREFIX.qc3.pruned.hapmap_merged.flipped.eval" \
    -l "$OUT_DIR/$IN_DATA_PREFIX.qc3.pruned.hapmap_merged.flipped.log" \
    -k 2 -t 2 \
    -w "/lustre/scratch113/teams/barrett/coreex_gaibdc/refs/pca-populations.txt"

# Choose a PC2 threshold to exclude samples of non-European ancestry
MIN_PC2="0.066"
Rscript "1.2_plot-pca-results_alpha.R" "$OUT_DIR/$IN_DATA_PREFIX.qc3.pruned.hapmap_merged.flipped.pca" "$MIN_PC2"

# Remove samples failing QC

# Mark individuals with mismatched gender.
# Mark one of each pair of individuals that are first degree relatives (kinship coefficient > 0.177) or closer: the individual with the lower call rate.
# See http://people.virginia.edu/~wc9c/KING/manual.html for thresholds.

# TODO nothing below here works
exit

# 1) All those with non-European ancestry (PC1)
cat /lustre/scratch113/teams/barrett/coreex_gaibdc/QC/COMBINED/3_sample_QC/coreex_gaibdc_usgwas_qc2.fail-ancestry.txt > /lustre/scratch113/teams/barrett/coreex_gaibdc/QC/COMBINED/3_sample_QC/sample_QC_fail.samples
# 2) All with heterozygosity rate more than 3s from the mean
cat /lustre/scratch113/teams/barrett/coreex_gaibdc/QC/COMBINED/3_sample_QC/coreex_gaibdc_usgwas_qc2.fail-het-QC.txt >> /lustre/scratch113/teams/barrett/coreex_gaibdc/QC/COMBINED/3_sample_QC/sample_QC_fail.samples
# 3) All duplicates or 1st-degree relatives
cat /lustre/scratch113/teams/barrett/coreex_gaibdc/QC/COMBINED/3_sample_QC/coreex_gaibdc_usgwas_qc2.duplicates_to_remove | awk '{print $1 "\t" $2}' >> /lustre/scratch113/teams/barrett/coreex_gaibdc/QC/COMBINED/3_sample_QC/sample_QC_fail.samples
cat /lustre/scratch113/teams/barrett/coreex_gaibdc/QC/COMBINED/3_sample_QC/coreex_gaibdc_usgwas_qc2.first_degree_to_remove | awk '{print $1 "\t" $2}' >> /lustre/scratch113/teams/barrett/coreex_gaibdc/QC/COMBINED/3_sample_QC/sample_QC_fail.samples
# 4) All gender mismatches (except those from the known plate-swap)
cat /lustre/scratch113/teams/barrett/coreex_gaibdc/QC/COMBINED/3_sample_QC/fail-gender.txt >> /lustre/scratch113/teams/barrett/coreex_gaibdc/QC/COMBINED/3_sample_QC/sample_QC_fail.samples
# 5) Those with elevated missing rates have already been removed.



#----------------------------------------------- PART 3 -----------------------------------------------#

# Check the missing data rate across all the markers in this batch.
bsub -J "missingness_cases" -R"select[mem>2000] rusage[mem=2000]" -M2000 -o logs/QC_missing.log -G crohns "plink --noweb --bfile /lustre/scratch113/teams/barrett/coreex_gaibdc/QC/COMBINED/3_sample_QC/coreex_gaibdc_usgwas_qc3 --missing --out /lustre/scratch113/teams/barrett/coreex_gaibdc/QC/COMBINED/4_marker_QC/coreex_gaibdc_usgwas_qc3"

Rscript /lustre/scratch113/teams/barrett/coreex_gaibdc/QC/scripts/lmiss-hist.Rscript /lustre/scratch113/teams/barrett/coreex_gaibdc/QC/COMBINED/4_marker_QC/coreex_gaibdc_usgwas_qc3

# Check for markers with different genotype call rates between cases and controls.
bsub -J "diffmiss" -R"select[mem>2000] rusage[mem=2000]" -M2000 -o /lustre/scratch113/teams/barrett/coreex_gaibdc/QC/COMBINED/logs/QC_diffmiss.log -G crohns "plink --noweb --bfile /lustre/scratch113/teams/barrett/coreex_gaibdc/QC/COMBINED/3_sample_QC/coreex_gaibdc_usgwas_qc3 --test-missing --out /lustre/scratch113/teams/barrett/coreex_gaibdc/QC/COMBINED/4_marker_QC/coreex_gaibdc_usgwas_qc3"


#bsub -J "diffmiss" -R"select[mem>2000] rusage[mem=2000]" -M2000 -o /lustre/scratch113/teams/barrett/coreex_gaibdc/QC/COMBINED/logs/QC_diffmiss.log -G crohns "perl /lustre/scratch113/teams/barrett/coreex_gaibdc/QC/scripts/run-diffmiss-qc.pl /lustre/scratch113/teams/barrett/coreex_gaibdc/QC/COMBINED/4_marker_QC/coreex_gaibdc_usgwas_qc3"

# run-diffmiss-qc.pl didn't work.. just do it with awk.
cat /lustre/scratch113/teams/barrett/coreex_gaibdc/QC/COMBINED/4_marker_QC/coreex_gaibdc_usgwas_qc3.missing | awk '($5 < 0.00001){print $2}' > /lustre/scratch113/teams/barrett/coreex_gaibdc/QC/COMBINED/4_marker_QC/fail-diffmiss-qc.txt

# Prune out those with significant diffmiss, and do another double check on the missingness
# now that lots of dodgy samples have been removed.
bsub -J "prune" -R"select[mem>2000] rusage[mem=2000]" -M2000 -o /lustre/scratch113/teams/barrett/coreex_gaibdc/QC/COMBINED/logs/QC_marker_prune_a.log -G crohns "plink --noweb --bfile /lustre/scratch113/teams/barrett/coreex_gaibdc/QC/COMBINED/3_sample_QC/coreex_gaibdc_usgwas_qc3 --exclude /lustre/scratch113/teams/barrett/coreex_gaibdc/QC/COMBINED/4_marker_QC/fail-diffmiss-qc.txt --geno 0.05 --make-bed --out /lustre/scratch113/teams/barrett/coreex_gaibdc/QC/COMBINED/4_marker_QC/coreex_gaibdc_usgwas_qc3a"


# Remove markers with the batch effects identified by Yang's PC analysis. This involved:
# - Computing within-sample PCs using common variants (MAF > 1%)
# - A clear case outlier set was seen. These were all part of batch3
# - Use PC1 to split cases into outliers and non-outliers. Perform an association
# 	test between these two groups, and identify 'significant' sites (P<1e-5)
# - Because this problematic batch is all cases, these will generate false hits if
#	left in. Thus all these sites (429) need to be removed.
bsub -J "prune" -R"select[mem>2000] rusage[mem=2000]" -M2000 -o /lustre/scratch113/teams/barrett/coreex_gaibdc/QC/COMBINED/logs/QC_marker_prune_3b.log -G crohns "plink --noweb --bfile /lustre/scratch113/teams/barrett/coreex_gaibdc/QC/COMBINED/4_marker_QC/coreex_gaibdc_usgwas_qc3a --exclude /lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/pcs/batch/batch-fail-1e-5.txt --make-bed --out /lustre/scratch113/teams/barrett/coreex_gaibdc/QC/COMBINED/4_marker_QC/coreex_gaibdc_usgwas_qc3b"

# Find markers with a HWE P-value < 0.00001 (in controls)
bsub -J "prune" -R"select[mem>2000] rusage[mem=2000]" -M2000 -o /lustre/scratch113/teams/barrett/coreex_gaibdc/QC/COMBINED/logs/QC_hwe_markers.log -G crohns "plink --noweb --bfile /lustre/scratch113/teams/barrett/coreex_gaibdc/QC/COMBINED/4_marker_QC/coreex_gaibdc_usgwas_qc3b --hardy --out /lustre/scratch113/teams/barrett/coreex_gaibdc/QC/COMBINED/4_marker_QC/coreex_gaibdc_usgwas_qc3b"

# Pull out those with a hwe < 0.00001. 
#bsub -J "prune" -R"select[mem>2000] rusage[mem=2000]" -M2000 -o /lustre/scratch113/teams/barrett/coreex_gaibdc/QC/COMBINED/logs/QC_marker_prune_b.log -G crohns -q yesterday "plink --noweb --bfile /lustre/scratch113/teams/barrett/coreex_gaibdc/QC/COMBINED/4_marker_QC/coreex_gaibdc_usgwas_qc3b --maf 0.05 --make-bed --out /lustre/scratch113/teams/barrett/coreex_gaibdc/QC/COMBINED/4_marker_QC/coreex_gaibdc_usgwas_qc3c"

bsub -J "prune" -R"select[mem>2000] rusage[mem=2000]" -M2000 -o /lustre/scratch113/teams/barrett/coreex_gaibdc/QC/COMBINED/logs/QC_marker_prune_b.log -G crohns -q yesterday "plink --noweb --bfile /lustre/scratch113/teams/barrett/coreex_gaibdc/QC/COMBINED/4_marker_QC/coreex_gaibdc_usgwas_qc3b --hwe 0.00001 midp --make-bed --out /lustre/scratch113/teams/barrett/coreex_gaibdc/QC/COMBINED/4_marker_QC/coreex_gaibdc_usgwas_qc4"


# Before imputation, prune to SNPs with a MAF > 0.1%
bsub -J "prune" -R"select[mem>2000] rusage[mem=2000]" -M2000 -o /lustre/scratch113/teams/barrett/coreex_gaibdc/QC/COMBINED/logs/QC_marker_prune.log -G crohns "plink --noweb --bfile /lustre/scratch113/teams/barrett/coreex_gaibdc/QC/COMBINED/4_marker_QC/coreex_gaibdc_usgwas_qc4 --maf 0.001 --make-bed --out /lustre/scratch113/teams/barrett/coreex_gaibdc/QC/COMBINED/4_marker_QC/coreex_gaibdc_usgwas_qcIMP"



# Before imputation, may want to prune based on MAF??
bsub -J "prune" -R"select[mem>2000] rusage[mem=2000]" -M2000 -o /lustre/scratch113/teams/barrett/coreex_gaibdc/QC/COMBINED/logs/QC_marker_prune.log -G crohns "plink --noweb --bfile /lustre/scratch113/teams/barrett/coreex_gaibdc/QC/COMBINED/4_marker_QC/coreex_gaibdc_usgwas_qc4 --maf 0.01 --make-bed --out /lustre/scratch113/teams/barrett/coreex_gaibdc/QC/COMBINED/4_marker_QC/coreex_gaibdc_usgwas_qc5"

bsub -J "prune" -R"select[mem>2000] rusage[mem=2000]" -M2000 -o /lustre/scratch113/teams/barrett/coreex_gaibdc/QC/COMBINED/logs/QC_marker_prune.log -G crohns "plink --noweb --bfile /lustre/scratch113/teams/barrett/coreex_gaibdc/QC/COMBINED/4_marker_QC/coreex_gaibdc_usgwas_qc4 --maf 0.005 --make-bed --out /lustre/scratch113/teams/barrett/coreex_gaibdc/QC/COMBINED/4_marker_QC/coreex_gaibdc_usgwas_qc6"


bsub -J "prune" -R"select[mem>2000] rusage[mem=2000]" -M2000 -o /lustre/scratch113/teams/barrett/coreex_gaibdc/QC/COMBINED/logs/QC_marker_prune.log -G crohns "plink --noweb --bfile /lustre/scratch113/teams/barrett/coreex_gaibdc/QC/COMBINED/4_marker_QC/coreex_gaibdc_usgwas_qc4 --maf 0.001 --make-bed --out /lustre/scratch113/teams/barrett/coreex_gaibdc/QC/COMBINED/4_marker_QC/coreex_gaibdc_usgwas_qc7"
