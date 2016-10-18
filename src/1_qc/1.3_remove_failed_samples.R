
library(plyr)

options(stringsAsFactors = F)

in.imiss.file <- commandArgs(T)[1]
in.het.file <- commandArgs(T)[2]
in.imiss.thresh <- as.numeric(commandArgs(T)[3])
in.het.thresh <- as.numeric(commandArgs(T)[4])
out.pdf.file <- commandArgs(T)[5]
out.failedIds.file <- commandArgs(T)[6]

#
# Identify samples failing missingness filter
# Identify samples failing heterozygosity filter
#

# Read in previously identified individuals
in.imiss.fails.file <- "../../crohns_workspace/1_qc/coreex_gaibdc_usgwas_raw.qc2.sample_fail_imiss.txt"
in.het.fails.file <- "../../crohns_workspace/1_qc/coreex_gaibdc_usgwas_raw.qc2.sample_fail_het.txt"

in.imiss.fails.df <- read.table(in.imiss.fails.file)
in.het.fails.df <- read.table(in.het.fails.file)

#
# Identify samples with discordant sex information
#
in.sexcheck.file <- "../../crohns_workspace/1_qc/coreex_gaibdc_usgwas_raw.qc3.check_sex.sexcheck"
in.sexcheck.df <- read.table(in.sexcheck.file, header = T)

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

# 
# Identify individuals with non-European ancestry
#

in.pca.file <- "../../crohns_workspace/1_qc/coreex_gaibdc_usgwas_raw.qc3.pruned.hapmap_merged.flipped.pca.evec"
in.pca.df <- read.table(in.pca.file, skip = 1)
colnames(in.pca.df) <- c("SAMPLE", "PC1", "PC2", "POP")

# Restore original sample IDs
in.hapmap.fam.file <- "../../crohns_workspace/1_qc/coreex_gaibdc_usgwas_raw.qc3.pruned.hapmap_merged.flipped.fam"
in.hapmap.fam.df <- read.table(in.hapmap.fam.file)

stopifnot(nrow(in.pca.df) == nrow(in.hapmap.fam.df))
in.pca.df$FID <- in.hapmap.fam.df$V1
in.pca.df$IID <- in.hapmap.fam.df$V2

# Identify failures based on PC2 coord
in.pc2.thresh <- as.numeric(0.066)
in.pca.df$ANCESTRY.FAIL <- in.pca.df$PC2 < in.pc2.thresh

#
# Summarise QC failures thus far
#

# Read in all original samples
# This is the fam file after qc1, which removed markers with excess missingness
in.fam.orig.file <- "../../crohns_workspace/1_qc/coreex_gaibdc_usgwas_raw.qc1.fam"
in.fam.orig.df <- read.table(in.fam.orig.file)

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

in.king.ibs0.file <- "../../crohns_workspace/1_qc/coreex_gaibdc_usgwas_raw.qc3.maf_0.05.king.ibs0"
in.king.ibs0.df <- read.table(in.king.ibs0.file, header = T)

# Read in sample call rates
in.imiss.file <- "../../crohns_workspace/1_qc/coreex_gaibdc_usgwas_raw.qc1.imiss"
in.imiss.df <- read.table(in.imiss.file, header=T)

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
            return(c(id1.fmiss <= id2.fmiss, id1.fmiss > id2.fmiss))
        }
    # Samples insufficiently closely related
    } else {
        return(c(F, F))
    }
}

in.kinship.thresh <- as.numeric(0.177)
pair.failures <- aaply(in.king.ibs0.df, 1, process.related.pair, in.kinship.thresh, failures.summary, in.imiss.df, .expand=F)

failures.summary$RELATEDNESS.FAIL <- F
failures.summary$RELATEDNESS.FAIL[failures.summary$IID %in% in.king.ibs0.df$ID1[pair.failures[, 1]]] <- T
failures.summary$RELATEDNESS.FAIL[failures.summary$IID %in% in.king.ibs0.df$ID2[pair.failures[, 2]]] <- T

# Write summary of failures
failures.summary <- within(failures.summary, {
    ANY.FAIL <- IMISS.FAIL | HET.FAIL | CHECKSEX.FAIL | ANCESTRY.FAIL | RELATEDNESS.FAIL
})

out.summary.file <- "../../crohns_workspace/1_qc/coreex_gaibdc_usgwas_raw.qc3.summary.txt"
sink(out.summary.file)
print(summary(failures.summary))
sink()

out.summaryTable.file <- "../../crohns_workspace/1_qc/coreex_gaibdc_usgwas_raw.qc3.summary_table.txt"
write.table(failures.summary, out.summaryTable.file)

# Write list of failed ids to remove
out.failedIds.file <- "../../crohns_workspace/1_qc/coreex_gaibdc_usgwas_raw.qc3.sample_fail_any.txt"
write.table(
    failures.summary[failures.summary$ANY.FAIL, c("FID", "IID")],
    file = out.failedIds.file, 
    sep = "\t", 
    quote = F, 
    row.names = F, col.names = F
)
