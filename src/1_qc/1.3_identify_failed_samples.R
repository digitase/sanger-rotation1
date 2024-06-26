#
# Identify and output samples that fail any QC checks.
#

library(plyr)

options(stringsAsFactors = F)

#
# Input files
#

print("Args:")
print(commandArgs(T))

in.fam.orig.file <- commandArgs(T)[1]
in.imiss.file <- commandArgs(T)[2]

in.imiss.fails.file <- commandArgs(T)[3]
in.het.fails.file <- commandArgs(T)[4]
in.sexcheck.file <- commandArgs(T)[5]
in.pca.file <- commandArgs(T)[6]
in.hapmap.fam.file <- commandArgs(T)[7]
in.king.ibs0.file <- commandArgs(T)[8]

in.pc2.thresh <- as.numeric(commandArgs(T)[9])
in.kinship.thresh <- as.numeric(commandArgs(T)[10])

out.summary.file <- commandArgs(T)[11]
out.summaryTable.file <- commandArgs(T)[12]
out.failedIds.file <- commandArgs(T)[13]

# in.fam.orig.file <- "../../crohns_workspace/1_qc/coreex_gaibdc_usgwas_raw.qc1.fam"
# in.imiss.file <- "../../crohns_workspace/1_qc/coreex_gaibdc_usgwas_raw.qc2.imiss"
# 
# in.imiss.fails.file <- "../../crohns_workspace/1_qc/coreex_gaibdc_usgwas_raw.qc2.sample_fail_imiss.txt"
# in.het.fails.file <- "../../crohns_workspace/1_qc/coreex_gaibdc_usgwas_raw.qc3.missing_het.sample_fail_het.txt"
# in.sexcheck.file <- "../../crohns_workspace/1_qc/coreex_gaibdc_usgwas_raw.qc3.check_sex.sexcheck"
# in.pca.file <- "../../crohns_workspace/1_qc/coreex_gaibdc_usgwas_raw.qc3.pruned.hapmap_merged.flipped.pca.evec"
# in.hapmap.fam.file <- "../../crohns_workspace/1_qc/coreex_gaibdc_usgwas_raw.qc3.pruned.hapmap_merged.flipped.fam"
# in.king.ibs0.file <- "../../crohns_workspace/1_qc/coreex_gaibdc_usgwas_raw.qc3.maf_0.05.king.ibs0"
# 
# in.pc2.thresh <- as.numeric(0.067)
# in.kinship.thresh <- as.numeric(0.177)
# 
# out.summary.file <- "../../crohns_workspace/1_qc/coreex_gaibdc_usgwas_raw.qc3.summary.txt"
# out.summaryTable.file <- "../../crohns_workspace/1_qc/coreex_gaibdc_usgwas_raw.qc3.summary_table.txt"
# out.failedIds.file <- "../../crohns_workspace/1_qc/coreex_gaibdc_usgwas_raw.qc3.sample_fail_any.txt"

#
# Identify samples failing missingness filter
# Identify samples failing heterozygosity filter
#

# read.table that handles empty input files
read.table.try <- function(filename, ...) {
    result <- tryCatch({
        read.table(filename, ...)
    }, error = function(e) {
        data.frame()
    })
}

# Read in previously identified individuals
in.imiss.fails.df <- read.table.try(in.imiss.fails.file)
in.het.fails.df <- read.table.try(in.het.fails.file)

#
# Identify samples with discordant sex information
#
in.sexcheck.df <- read.table.try(in.sexcheck.file, header = T)

# Account for known plate-swap:
#
# There was a plate swap, between samples on plates 333833 and 333835, 
# so you can check for those numbers in the sample names (PLATE_WELL_SAMPLE). 
# If you want to be super careful though, they fixed the genders in the 4th release of the data, 
# so you can check back against the fam file in here: 
# /lustre/scratch114/teams/barrett/coreex_gaibdc/release/coreex_gaibdc_20150304

# Swap the plate numbers
# Doesn't quite work, as the two plates have different numbers of samples
# in.sexcheck.df$PLATE_WELL <- sapply(strsplit(in.sexcheck.df$IID, split = "_"), function(x) paste(x[1], x[2], sep = "_"))
# in.sexcheck.df$PLATE_WELL.CORRECTED <- sub("^333833", "TMP333835", in.sexcheck.df$PLATE_WELL)
# in.sexcheck.df$PLATE_WELL.CORRECTED <- sub("^333835", "333833", in.sexcheck.df$PLATE_WELL.CORRECTED)
# in.sexcheck.df$PLATE_WELL.CORRECTED <- sub("^TMP333835", "333835", in.sexcheck.df$PLATE_WELL.CORRECTED)
# in.sexcheck.df$PEDSEX.CORRECTED <- in.sexcheck.df$PEDSEX[match(in.sexcheck.df$PLATE_WELL.CORRECTED, in.sexcheck.df$PLATE_WELL)]

# Ignore problems from these plates
in.sexcheck.df <- within(in.sexcheck.df, {
    CHECKSEX.FAIL <- (PEDSEX != SNPSEX) & (PEDSEX != 0) & (!grepl("^333833|^333835", IID))
})

# Also exclude "333835_G08_UC755045" due to mismatched sex compared to a later release
in.sexcheck.df$CHECKSEX.FAIL[in.sexcheck.df$IID == "333835_G08_UC755045"] <- T

# 
# Identify individuals with non-European ancestry
#
in.pca.df <- read.table.try(in.pca.file, skip = 1)
colnames(in.pca.df) <- c("SAMPLE", paste0("PC", 1:(ncol(in.pca.df)-2)), "POP")

# Unshorten sample IDs
in.hapmap.fam.df <- read.table.try(in.hapmap.fam.file)

# Assign original FID/IIDs
#
# Actually, not every hapmap sample will be included in the pca due to smartpca outlier removal
# stopifnot(nrow(in.pca.df) == nrow(in.hapmap.fam.df))
# in.pca.df$FID <- in.hapmap.fam.df$V1
# in.pca.df$IID <- in.hapmap.fam.df$V2
#
# So subset just the samples that are left
in.pca.df.indices <- unlist(Map(
    function (x) as.numeric(sub("ID", "", strsplit(x, split=':', fixed=T)[[1]][1])),
    in.pca.df$SAMPLE
))
in.pca.df$FID <- in.hapmap.fam.df$V1[in.pca.df.indices]
in.pca.df$IID <- in.hapmap.fam.df$V2[in.pca.df.indices]

# Identify failures based on PC2 coord
in.pca.df$ANCESTRY.FAIL <- in.pca.df$PC2 < in.pc2.thresh

#
# Summarise QC failures thus far
#

# Read in all original samples
# This is the fam file after qc1, which removed markers with excess missingness
in.fam.orig.df <- read.table.try(in.fam.orig.file)

# Add QC failures
failures.summary <- data.frame(FID = in.fam.orig.df$V1, IID = in.fam.orig.df$V2)

failures.summary$IMISS.FAIL <- F
failures.summary$IMISS.FAIL[failures.summary$IID %in% in.imiss.fails.df$V2] <- T

failures.summary$HET.FAIL <- F
failures.summary$HET.FAIL[failures.summary$IID %in% in.het.fails.df$V2] <- T

failures.summary$CHECKSEX.FAIL <- F
failures.summary$CHECKSEX.FAIL[failures.summary$IID %in% in.sexcheck.df$IID[in.sexcheck.df$CHECKSEX.FAIL]] <- T

failures.summary$ANCESTRY.FAIL <- F
failures.summary$ANCESTRY.FAIL[failures.summary$IID %in% in.pca.df$IID[in.pca.df$ANCESTRY.FAIL]] <- T

# 
# Identify duplicated and highly related individuals
#
in.king.ibs0.df <- read.table.try(in.king.ibs0.file, header = T)

# Read in sample call rates
in.imiss.df <- read.table.try(in.imiss.file, header=T)

# For each related pair, if both members passed other QC, mark the member with the lower call rate.
# Return whether each member of the pair fails.
process.related.pair <- function(x, kinship.thresh, failures.summary, in.imiss.df) {
    kinship <- x$Kinship
    # Threshold based on kinship
    if (kinship > kinship.thresh) {
        id1 <- x$ID1
        id2 <- x$ID2
        # Check for failures of other QC
        id1.prev_failed <- any(unlist(failures.summary[failures.summary$IID == id1, -c(1, 2)]))
        id2.prev_failed <- any(unlist(failures.summary[failures.summary$IID == id2, -c(1, 2)]))
        # If exactly one member failed previous QC, mark that member
        if (xor(id1.prev_failed, id2.prev_failed)) {
            return(c(id1.prev_failed, id2.prev_failed))
        # else mark the member with lower call rate.    
        } else {
            id1.fmiss <- in.imiss.df$F_MISS[in.imiss.df$IID == id1]
            id2.fmiss <- in.imiss.df$F_MISS[in.imiss.df$IID == id2]
            return(c(id1.fmiss >= id2.fmiss, id1.fmiss < id2.fmiss))
        }
    # Samples insufficiently closely related
    } else {
        return(c(F, F))
    }
}

pair.failures <- aaply(in.king.ibs0.df, 1, process.related.pair, in.kinship.thresh, failures.summary, in.imiss.df, .expand=F)

failures.summary$RELATEDNESS.FAIL <- F
failures.summary$RELATEDNESS.FAIL[failures.summary$IID %in% in.king.ibs0.df$ID1[pair.failures[, 1]]] <- T
failures.summary$RELATEDNESS.FAIL[failures.summary$IID %in% in.king.ibs0.df$ID2[pair.failures[, 2]]] <- T

# Write summary of failures
failures.summary <- within(failures.summary, {
    ANY.FAIL <- IMISS.FAIL | HET.FAIL | CHECKSEX.FAIL | ANCESTRY.FAIL | RELATEDNESS.FAIL
})

sink(out.summary.file)
print(summary(failures.summary))
sink()

write.table(failures.summary, out.summaryTable.file)

# Write list of failed ids to remove
write.table(
    failures.summary[failures.summary$ANY.FAIL, c("FID", "IID")],
    file = out.failedIds.file, 
    sep = "\t", 
    quote = F, 
    row.names = F, col.names = F
)
