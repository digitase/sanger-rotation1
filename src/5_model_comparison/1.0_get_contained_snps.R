#!/software/R-3.3.0/bin/Rscript
#
# Read in snptest model output
# 
#

library(data.table)

options(stringsAsFactors=F)

out.dir <- "/nfs/users/nfs_b/bb9/workspace/rotation1/crohns_workspace/5_model_comparison/"
dir.create(out.dir)

# Load table of known loci
knownLoci.df <- read.csv("/nfs/users/nfs_b/bb9/workspace/rotation1/data/known_ibd_loci/Table S2 - association statistics at all known loci.csv", header=T, comment.char="#")
knownLoci.dt <- data.table(knownLoci.df)
# Generate unique id for each locus
knownLoci.dt$locus.id <- paste(knownLoci.dt$Chr, knownLoci.dt$"LD_left", knownLoci.dt$"LD_right", sep="_")
setkey(knownLoci.dt, Chr, LD_left, LD_right)

knownLoci.snptest.results <- data.table()

# Read in snptest rows within known loci
# for (snptest.result.file in list.files("/nfs/users/nfs_b/bb9/workspace/rotation1/crohns_workspace/4_gwas/1_snptest/results.filtered/gwas3/ibd/", pattern="*.snptest.filtered.out", full.names=T)[14:15]) {
for (snptest.result.file in list.files("/nfs/users/nfs_b/bb9/workspace/rotation1/crohns_workspace/4_gwas/1_snptest/results.filtered/gwas3/ibd/", pattern="*.snptest.filtered.out", full.names=T)) {

    print(paste("Processing", snptest.result.file))
    # snptest.result.dt <- fread("/nfs/users/nfs_b/bb9/workspace/rotation1/crohns_workspace/4_gwas/1_snptest/results.filtered/gwas3/ibd/chr_22.snptest.filtered.out", header=T)
    snptest.result.dt <- fread(snptest.result.file, header=T)

    # Duplicate snp positions to create 0 length ranges
    snptest.result.dt$position_end <- snptest.result.dt$position

    # Find overlaps with known loci
    contained.snps <- foverlaps(snptest.result.dt, knownLoci.dt, by.x=c("alternate_ids", "position", "position_end"), which=F, nomatch=0)
    print(paste("Found", nrow(contained.snps), "snps within known loci."))
    knownLoci.snptest.results <- rbind(knownLoci.snptest.results, contained.snps)

}

# snptest.pvalue.names <- names(knownLoci.snptest.results)[grep('pvalue', names(knownLoci.snptest.results))]

# Sanity checks
stopifnot(all(knownLoci.snptest.results$alternate_ids == knownLoci.snptest.results$Chr))
stopifnot(all(knownLoci.snptest.results$position == knownLoci.snptest.results$position_end))
stopifnot(all(knownLoci.snptest.results$LD_left <= knownLoci.snptest.results$position))
stopifnot(all(knownLoci.snptest.results$LD_right >= knownLoci.snptest.results$position))

# Filter out models that failed to fit (singular values, convergence etc.)
models.failed <- which(
    is.na(knownLoci.snptest.results$frequentist_add_pvalue) &
    is.na(knownLoci.snptest.results$frequentist_dom_pvalue) &
    is.na(knownLoci.snptest.results$frequentist_rec_pvalue) &
    is.na(knownLoci.snptest.results$frequentist_gen_pvalue) &
    is.na(knownLoci.snptest.results$frequentist_het_pvalue)
)
# More restrictive version
# models.failed <- which(is.na(knownLoci.snptest.results$frequentist_add_pvalue))

knownLoci.snptest.results <- knownLoci.snptest.results[-models.failed]

write.table(knownLoci.snptest.results, file.path(out.dir, "contained.snps.txt"), sep="\t", row.names=F)

# TODO below here
stopifnot(F)

# Get rows of known topSNPs per loci
# Make sure they are all there
stopifnot(nrow(knownLoci.snptest.results[topSNP.Position..bp. == position, ]) == nrow(knownLoci.dt))

# Check if topSNPs in snptest output are consistent across models
knownLoci.snptest.results[
    , 
    .(
        add=which.min(frequentist_add_pvalue),
        dom=which.min(frequentist_dom_pvalue),
        rec=which.min(frequentist_rec_pvalue),
        gen=which.min(frequentist_gen_pvalue),
        het=which.min(frequentist_het_pvalue)
    ), 
    by=locus.id
]

