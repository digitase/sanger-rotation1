#!/software/R-3.3.0/bin/Rscript
#
# Read in snptest model output
# Write out significant snps within known loci
#
library(data.table)
library(plyr)
library(doMC)

registerDoMC(cores=2)
options(stringsAsFactors=F)

out.dir <- "/nfs/users/nfs_b/bb9/workspace/rotation1/crohns_workspace/5_model_comparison/"
dir.create(out.dir)

# Load table of known loci
knownLoci.dt <- data.table(read.csv("/nfs/users/nfs_b/bb9/workspace/rotation1/data/known_ibd_loci/Table S2 - association statistics at all known loci.csv", header=T, comment.char="#"))
# Generate unique id for each locus
knownLoci.dt$locus.id <- paste(knownLoci.dt$Chr, knownLoci.dt$"LD_left", knownLoci.dt$"LD_right", sep="_")
setkey(knownLoci.dt, Chr, LD_left, LD_right)

# Read in snptest rows within known loci
knownLoci.snptest.results <- data.table(ldply(list.files("/nfs/users/nfs_b/bb9/workspace/rotation1/crohns_workspace/4_gwas/1_snptest/results.filtered/gwas3/ibd/", pattern="*.snptest.filtered.out", full.names=T), function(snptest.result.file) {
# knownLoci.snptest.results <- data.table(ldply(list.files("/nfs/users/nfs_b/bb9/workspace/rotation1/crohns_workspace/4_gwas/1_snptest/results.filtered/gwas3/ibd/", pattern="*.snptest.filtered.out", full.names=T)[14:15], function(snptest.result.file) {
    print(paste("Processing", snptest.result.file))
    snptest.result.dt <- fread(snptest.result.file, header=T)
    # Duplicate snp positions to create 0 length ranges
    snptest.result.dt$position_end <- snptest.result.dt$position
    # Find overlaps with known loci
    contained.snps <- foverlaps(snptest.result.dt, knownLoci.dt, by.x=c("alternate_ids", "position", "position_end"), which=F, nomatch=0)
    print(paste("Found", nrow(contained.snps), "snps within known loci."))
    return(contained.snps)
}, .parallel=F))

# Note that known loci bound can overlap, so a snp may fall into multiple loci

# Sanity checks
stopifnot(all(knownLoci.snptest.results$alternate_ids == knownLoci.snptest.results$Chr))
stopifnot(all(knownLoci.snptest.results$position == knownLoci.snptest.results$position_end))
stopifnot(all(knownLoci.snptest.results$LD_left <= knownLoci.snptest.results$position))
stopifnot(all(knownLoci.snptest.results$LD_right >= knownLoci.snptest.results$position))
print(paste("Total number of contained loci:", nrow(knownLoci.snptest.results)))

# Get maximum IBD topsnp additive snptest pvalues
# Ensures at least the topsnp for each IBD locus will be output, and any snps that fit at least as well.
# gwas.thresh <- max(knownLoci.snptest.results[grepl("IBD", Trait) & (topSNP.Position..bp. == position), frequentist_add_pvalue])

# ... or specify generic thresh if the above thresh is too high...
gwas.thresh <- 5e-8

# For each locus, determine the local p value threshold to use:
# max((if the topSNP exists in the snptest results, topSNP additive p value, else 0), generic gwas thresh)
#
# Determine if each locus has a topsnp
knownLoci.snptest.results[, is.topsnp := topSNP.Position..bp. == position]
knownLoci.snptest.results[, has.topsnp := any(is.topsnp), by=locus.id]
# Calculate locus thresholds
knownLoci.snptest.results[
    , 
    locus.gwas.thresh := {
        if (has.topsnp) {
            max(frequentist_add_pvalue[is.topsnp], gwas.thresh)
        } else {
            gwas.thresh
        }
    }, 
    by=.(locus.id, has.topsnp)
]

# Check that all thresholds are at least as large as the generic threshold
stopifnot(all(knownLoci.snptest.results$locus.gwas.thresh >= gwas.thresh))
# Check the correct number of unique locus thresholds has been returned
stopifnot(
    # Number of unique locus thresholds ==
    length(unique(knownLoci.snptest.results$locus.gwas.thresh)) == 
    # Number of loci -
    length(unique(knownLoci.snptest.results$locus.id)) - 
    # Number of loci without a topsnp -
    sum(!knownLoci.snptest.results[, NA, by=.(locus.id, has.topsnp)]$has.topsnp) -
    # Number of loci with topsnp with an ADD p value of < gwas.thresh +
    sum(knownLoci.snptest.results[, frequentist_add_pvalue[is.topsnp] < gwas.thresh, by=.(locus.id, has.topsnp)]$V1) +
    # The generic threshold
    1 
)

# Enforce the threshold at each locus
models.signif <- (
    knownLoci.snptest.results$frequentist_add_pvalue <= knownLoci.snptest.results$locus.gwas.thresh |
    knownLoci.snptest.results$frequentist_dom_pvalue <= knownLoci.snptest.results$locus.gwas.thresh |
    knownLoci.snptest.results$frequentist_rec_pvalue <= knownLoci.snptest.results$locus.gwas.thresh |
    knownLoci.snptest.results$frequentist_gen_pvalue <= knownLoci.snptest.results$locus.gwas.thresh |
    knownLoci.snptest.results$frequentist_het_pvalue <= knownLoci.snptest.results$locus.gwas.thresh 
)
print(paste("SNPs that reached the locus pvalue in at least one model:", sum(models.signif, na.rm=T)))

# ... or use the generic threshold only.
# models.signif <- (
    # knownLoci.snptest.results$frequentist_add_pvalue <= gwas.thresh |
    # knownLoci.snptest.results$frequentist_dom_pvalue <= gwas.thresh |
    # knownLoci.snptest.results$frequentist_rec_pvalue <= gwas.thresh |
    # knownLoci.snptest.results$frequentist_gen_pvalue <= gwas.thresh |
    # knownLoci.snptest.results$frequentist_het_pvalue <= gwas.thresh 
# )
# print(paste("SNPs that reached the generic pvalue threshold in at least one model:", sum(models.signif, na.rm=T), gwas.thresh))

# Filter and write out contained loci
print(paste("topSNPs in known loci table:", nrow(knownLoci.dt)))
print(paste("topSNPs found in snptest results:", sum(knownLoci.snptest.results$is.topsnp)))
knownLoci.snptest.results.filtered <- knownLoci.snptest.results[models.signif | is.topsnp]
print(paste("SNPs that reached the locus pvalue threshold in at least one model OR is topSNP:", nrow(knownLoci.snptest.results.filtered)))
write.table(knownLoci.snptest.results.filtered, file.path(out.dir, "contained.snps.txt"), sep="\t", row.names=F)

