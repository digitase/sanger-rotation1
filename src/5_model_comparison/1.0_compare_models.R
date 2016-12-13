#!/software/R-3.3.0/bin/Rscript
#
library(data.table)
library(gtools)
library(plyr)
library(ggplot2)
library(gplots)
library(gdata)
library(gridExtra)
library(irr)
library(doMC)

dataset <- "gwas3"
# assoc <- "cd"
# assoc <- "uc"
assoc <- "ibd"

outdir.base <- "/nfs/users/nfs_b/bb9/workspace/rotation1/crohns_workspace/5_model_comparison/1.0_compare_models/"
outdir <- file.path(outdir.base, dataset, assoc)
dir.create(outdir, recursive=T)
setwd(outdir)

# Start log file
sink("misc_output.txt")
sink()

# Read in snptest results for loci that were re run in glm
contained.snps.dt <- fread(sprintf("/nfs/users/nfs_b/bb9/workspace/rotation1/crohns_workspace/4_gwas/2_R_glm/contained.snps.%s.%s.all.txt", dataset, assoc))

knownLoci.dt <- data.table(read.csv("/nfs/users/nfs_b/bb9/workspace/rotation1/data/known_ibd_loci/Table_S2_-_association_statistics_at_all_known_loci.csv", header=T, comment.char="#"))
knownLoci.dt$locus.id <- paste(knownLoci.dt$Chr, knownLoci.dt$"LD_left", knownLoci.dt$"LD_right", sep="_")
# There is one locus (width 1) that did not have any contained snps.
# This dropped out at the post snptest MAF filter stage.
setdiff(knownLoci.dt$locus.id, unique(contained.snps.dt$locus.id))

# Verify snptest output filtering worked
stopifnot(all(contained.snps.dt$info > 0.4))
stopifnot(all(contained.snps.dt$all_maf > 0.001))

# Read in results from glm
all.glm.out.files <- mixedsort(list.files(file.path("/nfs/users/nfs_b/bb9/workspace/rotation1/crohns_workspace/4_gwas/2_R_glm/chunks/", dataset, assoc), include.dirs=F, recursive=T, full.names=T))
length(all.glm.out.files)
glm.results.dt <- data.table(ldply(all.glm.out.files, fread, .progress="text"))

# Merge in snptest p values
merged.results.dt <- merge(contained.snps.dt, glm.results.dt, by.x=c("rsid"), by.y=c("snp.id"), all.y=T)
dim(merged.results.dt)
# Find topsnps
merged.results.dt[, is.topsnp.meta := topSNP.Position..bp. == position.x, by=locus.id]
merged.results.dt[, is.topsnp.observed := frequentist_add_pvalue == min(frequentist_add_pvalue, na.rm=T), by=locus.id]
merged.results.dt[, is.locus.implicated := grepl(toupper(assoc), Trait)]

# Extra round of MAF filtering
merged.results.dt.orig <- merged.results.dt
# Allow through particular NOD2 snps that will be specifically checked downstream
merged.results.dt <- merged.results.dt[all_maf >= 0.01 | (alternate_ids == 16 & position.x %in% c(50745926, 50756540, 50763778))]

#
# Start comparisons
#

#
# Question 0.0
#
# Can we trust our association summary stats vs. snptest's?
#

if (F) {
# Coefficient p value is almost the same p value as the lrt
ggplot(data=merged.results.dt, aes(x=ADD_pValue.LRT, y=ADD_pValue.coef)) + 
    geom_abline(slope=1, intercept=0, linetype="dashed") +
    geom_point(aes(colour=frequentist_add_info), size=0.5) +
    scale_colour_gradient(low="red") +
    scale_x_log10() + scale_y_log10()

# Check that robust LRT is concordant with traditional LRT
ggplot(data=merged.results.dt[NULL_ADD_pValue.vuongtest.LRT > 0], aes(x=NULL_ADD_pValue.vuongtest.LRT, y=ADD_pValue.LRT)) +
    geom_abline(slope=1, intercept=0, linetype="dashed") +
    geom_point(aes(colour=frequentist_add_info), size=0.5) +
    scale_colour_gradient(low="red") +
    scale_x_log10() + scale_y_log10() 
}

# LRT value corresponds well to snptest values with high INFO
svg("0.0_frequentist_add-gen_pvalue_vs_ADD-GEN_pValue.LRT.svg", width=14, height=7)
    grid.arrange(
        # For ADD
        ggplot(data=merged.results.dt, aes(x=frequentist_add_pvalue, y=ADD_pValue.LRT)) +
            geom_abline(slope=1, intercept=0, linetype="dashed", col="grey40") +
            geom_hline(yintercept=5e-8, linetype="dotted", col="grey40") +
            geom_vline(xintercept=5e-8, linetype="dotted", col="grey40") +
            geom_point(aes(colour=frequentist_add_info), size=0.3) +
            scale_colour_gradient(low="red") +
            scale_x_log10() + scale_y_log10(),
        # For GEN
        ggplot(data=merged.results.dt, aes(x=frequentist_gen_pvalue, y=GEN_pValue.LRT)) +
            geom_abline(slope=1, intercept=0, linetype="dashed", col="grey40") +
            # geom_abline(slope=1.2, intercept=0, linetype="dashed", col="blue") +
            geom_hline(yintercept=5e-8, linetype="dotted", col="grey40") +
            geom_vline(xintercept=5e-8, linetype="dotted", col="grey40") +
            geom_point(aes(colour=frequentist_gen_info), size=0.3) +
            scale_colour_gradient(low="red") +
            # geom_text(data=merged.results.dt[log10(GEN_pValue.LRT) < log10(frequentist_gen_pvalue)*1.2 & frequentist_gen_pvalue < 5e-8 & GEN_pValue.LRT < 5e-8], aes(label="O"), size=3) +
            scale_x_log10() + scale_y_log10(),
        nrow=1, ncol=2 
    )
dev.off()

# Investigate points below 1.2x line in GEN plot
# Almost all MHC for IBD
merged.results.dt[log10(GEN_pValue.LRT) < log10(frequentist_gen_pvalue)*1.2 & frequentist_gen_pvalue < 5e-8 & GEN_pValue.LRT < 5e-8]

#
# Detecting evidence for non additive models being better
#
better <- data.table(rsid=merged.results.dt$rsid)

# Determine the significance threshold to use
# Bonferroni on the number of 'independent' loci
# TODO, use this threshold for icci
signif.thresh <- 0.05 / length(unique(merged.results.dt$locus.id))

#
# A decision based on raw p value
#
better[, DOM.raw.pValue := merged.results.dt[, frequentist_dom_pvalue < frequentist_add_pvalue]]
better[, REC.raw.pValue := merged.results.dt[, frequentist_rec_pvalue < frequentist_add_pvalue]]
better[, GEN.raw.pValue := merged.results.dt[, frequentist_gen_pvalue < frequentist_add_pvalue]]
better[, HET.raw.pValue := merged.results.dt[, frequentist_het_pvalue < frequentist_add_pvalue]]

#
# A decision based on raw BIC
#
raw.BIC.thresh <- 6
better[, DOM.raw.BIC := merged.results.dt[, NULL_ADD_BIC2 - ADD_DOM_BIC2 > raw.BIC.thresh]]
better[, REC.raw.BIC := merged.results.dt[, NULL_ADD_BIC2 - ADD_REC_BIC2 > raw.BIC.thresh]]
better[, GEN.raw.BIC := merged.results.dt[, NULL_ADD_BIC2 - ADD_GEN_BIC2 > raw.BIC.thresh]]
better[, HET.raw.BIC := merged.results.dt[, NULL_ADD_BIC2 - ADD_HET_BIC2 > raw.BIC.thresh]]

#
# Decision methods for ADD vs GEN (nested)
#
# Use the traditional LRT
better[, GEN.LRT := merged.results.dt[, ADD_GEN_pValue.LRT < signif.thresh]]
#
# Use vuong 2 step test (nested version)
#
better[
    , 
    GEN.vuong := merged.results.dt[
        , 
        ADD_GEN_pValue.vuongtest.varTest > 0 & 
        ADD_GEN_pValue.vuongtest.LRT > 0 &
        ADD_GEN_pValue.vuongtest.varTest < signif.thresh &
        ADD_GEN_pValue.vuongtest.LRT < 0.05
    ]
]

#
# Decision methods for non nested models
#
# Compare concordance of vuong LRT and BIC confidence intervals for non nested models
# More difficult to apply FDR to ci, so prefer LRT if possible
#
# Dodgy negative p values from vartest have been observed in conjunction with this error:
# 5: In imhof(n * omega.hat.2, lamstar^2) :
#    Note that Qq + abserr is positive.
#
# Use vuong 2 step test (non nested version)
better[
    , 
    DOM.vuong := merged.results.dt[
        , 
        ADD_DOM_pValue.vuongtest.varTest > 0 & 
        ADD_DOM_pValue.vuongtest.LRT > 0 &
        ADD_DOM_pValue.vuongtest.varTest < signif.thresh &
        ADD_DOM_pValue.vuongtest.LRT < 0.05
    ]
]
better[
    , 
    REC.vuong := merged.results.dt[
        , 
        ADD_REC_pValue.vuongtest.varTest > 0 & 
        ADD_REC_pValue.vuongtest.LRT > 0 &
        ADD_REC_pValue.vuongtest.varTest < signif.thresh &
        ADD_REC_pValue.vuongtest.LRT < 0.05
    ]
]
better[
    , 
    HET.vuong := merged.results.dt[
        , 
        ADD_HET_pValue.vuongtest.varTest > 0 & 
        ADD_HET_pValue.vuongtest.LRT > 0 &
        ADD_HET_pValue.vuongtest.varTest < signif.thresh &
        ADD_HET_pValue.vuongtest.LRT < 0.05
    ]
]

# A decision based on BICci (equiv. to AIC)
#
# Note: if models are nested or if the "variance test" from
     # ‘vuongtest()’ indicates models are indistinguishable, then the
     # intervals returned from ‘icci()’ will be incorrect.
# Requires models to be distinguishable
# Requires 0 to not be in the ci, and the decrease to be positive
#
better[, DOM.BIC := merged.results.dt[, ADD_DOM_pValue.vuongtest.varTest > 0 & ADD_DOM_pValue.vuongtest.varTest < signif.thresh & ADD_DOM_BICci.lower * ADD_DOM_BICci.upper > 0 & ADD_DOM_BICci.lower > 0 ]]
better[, REC.BIC := merged.results.dt[, ADD_REC_pValue.vuongtest.varTest > 0 & ADD_REC_pValue.vuongtest.varTest < signif.thresh & ADD_REC_BICci.lower * ADD_REC_BICci.upper > 0 & ADD_REC_BICci.lower > 0 ]]
better[, HET.BIC := merged.results.dt[, ADD_HET_pValue.vuongtest.varTest > 0 & ADD_HET_pValue.vuongtest.varTest < signif.thresh & ADD_HET_BICci.lower * ADD_HET_BICci.upper > 0 & ADD_HET_BICci.lower > 0 ]]

# A decision based on coxtest
better[, DOM.coxtest := merged.results.dt[, ADD_DOM_pValue.coxtest > 0 & ADD_DOM_pValue.coxtest < signif.thresh & ADD_DOM_pValue.coxtest.converse > signif.thresh]]
better[, REC.coxtest := merged.results.dt[, ADD_REC_pValue.coxtest > 0 & ADD_REC_pValue.coxtest < signif.thresh & ADD_REC_pValue.coxtest.converse > signif.thresh]]
better[, HET.coxtest := merged.results.dt[, ADD_HET_pValue.coxtest > 0 & ADD_HET_pValue.coxtest < signif.thresh & ADD_HET_pValue.coxtest.converse > signif.thresh]]
# A decision based on jtest
better[, DOM.jtest := merged.results.dt[, ADD_DOM_pValue.jtest > 0 & ADD_DOM_pValue.jtest < signif.thresh & ADD_DOM_pValue.jtest.converse > signif.thresh]]
better[, REC.jtest := merged.results.dt[, ADD_REC_pValue.jtest > 0 & ADD_REC_pValue.jtest < signif.thresh & ADD_REC_pValue.jtest.converse > signif.thresh]]
better[, HET.jtest := merged.results.dt[, ADD_HET_pValue.jtest > 0 & ADD_HET_pValue.jtest < signif.thresh & ADD_HET_pValue.jtest.converse > signif.thresh]]

# Filter data to leave only genome wide significant snps 
gwas.thresh <- 5e-8

merged.results.dt[
    , 
    signif.gwas.thresh :=
        frequentist_add_pvalue < gwas.thresh |
        frequentist_dom_pvalue < gwas.thresh |
        frequentist_rec_pvalue < gwas.thresh |
        frequentist_gen_pvalue < gwas.thresh |
        frequentist_het_pvalue < gwas.thresh
]
merged.results.dt$signif.gwas.thresh[is.na(merged.results.dt$signif.gwas.thresh)] <- F

#
# Question 1.0
# Basic numbers on the original and filtered analyses.
#
# How concordant is this analysis with the meta analysis?
#
sink("misc_output.txt", append=T, split=T)
#
print("How many snps total?")
nrow(merged.results.dt)
print("How many snps left?")
nrow(merged.results.dt[signif.gwas.thresh == T])
print("How many known loci remain?")
length(unique(merged.results.dt$locus.id))
print("How many known topsnps remain?")
sum(merged.results.dt$is.topsnp.meta)
#
print("How many known assoc specific loci remain?")
length(unique(merged.results.dt[grepl(toupper(assoc), Trait), locus.id]))
print("How many known assoc specific topsnps remain?")
sum(merged.results.dt[grepl(toupper(assoc), Trait), is.topsnp.meta])
#
print("How many known topsnps reached signif (ADD model only)?")
sum(merged.results.dt[signif.gwas.thresh == T, is.topsnp.meta])
print("How many known loci reached signif in at least one model?")
length(unique(merged.results.dt[signif.gwas.thresh == T, locus.id]))
print("How many gwas3 topsnps reached signif (ADD model only)?")
# There are loci where the gwas3 topsnp (based on frequentist_add_pvalue) did not reach significance,
# but another snp in the locus did (based on a non additive p value).
sum(merged.results.dt[signif.gwas.thresh == T, is.topsnp.observed])
print("Distribution of signif snps per locus:")
table(merged.results.dt[signif.gwas.thresh == T, locus.id])
#
sink()

plotVennWrapper <- function(dt, names) {
    subset.dt <- dt[, names, with=F]
    subset.dt[is.na(subset.dt)] <- F
    v <- venn(subset.dt, show.plot=T)
}

svg("1.0_distr_of_topsnps_in_signif_loci.svg", width=7, height=7)
    # plotVennWrapper(merged.results.dt[signif.gwas.thresh == T & (is.topsnp.meta|is.topsnp.observed)], c("is.topsnp.meta", "is.topsnp.observed", "is.locus.implicated"))
    plotVennWrapper(
        merged.results.dt[
            # Consider the topsnps in loci that contain snps significant in ANY gwas3 model...
            signif.gwas.thresh == T, 
            .(
                # Is the topsnp in a locus implicated in the IBD meta analysis?
                in.implicated.locus=all(is.locus.implicated), 
                # Is the topsnp significant in gwas3 using ADD model?
                is.signif.in.gwas3=any(is.topsnp.observed),
                #
                is.signif.in.gwas3.and.meta=any(is.topsnp.meta & is.topsnp.observed)
            ), 
            by=locus.id
        ],
        c("in.implicated.locus", "is.signif.in.gwas3", "is.signif.in.gwas3.and.meta")
    )
    # For those (3 for assoc=IBD) loci that are implicated but not significant according to the additive model, 
    # a non additive model caused the locus to pass the threshold.is.topsnp.meta
dev.off()

#
# Question 2.0
# How should we judge if a snp displays sufficient evidence for a non additive model being better?
#
better.onlySignif <- better[merged.results.dt$signif.gwas.thresh]

# Compare models
if (F) {
plotVennWrapper(better.onlySignif, c("DOM.raw.pValue", "REC.raw.pValue", "GEN.raw.pValue", "HET.raw.pValue"))
plotVennWrapper(better.onlySignif, c("DOM.raw.BIC", "REC.raw.BIC", "GEN.raw.BIC", "HET.raw.BIC"))
plotVennWrapper(better.onlySignif, c("DOM.BIC", "REC.BIC", "HET.BIC"))
plotVennWrapper(better.onlySignif, c("DOM.vuong", "REC.vuong", "HET.vuong"))
plotVennWrapper(better.onlySignif, c("DOM.coxtest", "REC.coxtest", "HET.coxtest"))
plotVennWrapper(better.onlySignif, c("DOM.jtest", "REC.jtest", "HET.jtest"))
}

# Same model, compare methods
svg("2.0_compare_model_selection_methods.DOM.svg", width=7, height=7)
    plotVennWrapper(better.onlySignif, c("DOM.raw.pValue", "DOM.vuong", "DOM.BIC", "DOM.raw.BIC"))
dev.off()
svg("2.0_compare_model_selection_methods.REC.svg", width=7, height=7)
    plotVennWrapper(better.onlySignif, c("REC.raw.pValue", "REC.vuong", "REC.BIC", "REC.raw.BIC"))
dev.off()
svg("2.0_compare_model_selection_methods.GEN.svg", width=7, height=7)
    plotVennWrapper(better.onlySignif, c("GEN.raw.pValue", "GEN.raw.BIC", "GEN.vuong", "GEN.LRT"))
dev.off()
svg("2.0_compare_model_selection_methods.HET.svg", width=7, height=7)
    plotVennWrapper(better.onlySignif, c("HET.raw.pValue", "HET.vuong", "HET.BIC", "HET.raw.BIC"))
dev.off()

# Heatmap of method concordance
#
# better.molten <- melt(better[order(rowSums(better[, -1, with=F], na.rm=T)), sort(names(better)), with=F], id.vars=c("rsid"))
# ggplot(better.molten, aes(rsid, variable)) + 
    # geom_tile(aes(fill = value)) + theme(axis.text.x = element_blank())

#
# Investigate properties of coxtest and j test
# TODO suspected breakdown for HET
#
# Plot coxtest and j test p values against decrease in BIC
# par(mfrow=c(1, 3))
# plot(merged.results.dt$NULL_ADD_BIC2 - merged.results.dt$ADD_DOM_BIC2, merged.results.dt$ADD_DOM_pValue.coxtest)
# plot(merged.results.dt$NULL_ADD_BIC2 - merged.results.dt$ADD_REC_BIC2, merged.results.dt$ADD_REC_pValue.coxtest)
# plot(merged.results.dt$NULL_ADD_BIC2 - merged.results.dt$ADD_HET_BIC2, merged.results.dt$ADD_HET_pValue.coxtest)

# par(mfrow=c(1, 3))
# plot(merged.results.onlySignif.dt$NULL_ADD_BIC2 - merged.results.onlySignif.dt$ADD_DOM_BIC2, merged.results.onlySignif.dt$ADD_DOM_pValue.jtest.FDR)
# plot(merged.results.onlySignif.dt$NULL_ADD_BIC2 - merged.results.onlySignif.dt$ADD_REC_BIC2, merged.results.onlySignif.dt$ADD_REC_pValue.jtest.FDR)
# plot(merged.results.onlySignif.dt$NULL_ADD_BIC2 - merged.results.onlySignif.dt$ADD_HET_BIC2, merged.results.onlySignif.dt$ADD_HET_pValue.jtest.FDR)

# Plot % significant as decided by the test for different bic decreases
# Poor performance at high BICs: few with that high of a BIC all get rejected
# merged.results.onlySignif.dt[NULL_ADD_BIC2 - ADD_HET_BIC2 > 45]
#
# BIC.breaks <- seq(min(merged.results.onlySignif.dt[, NULL_ADD_BIC2 - ADD_DOM_BIC2]), max(merged.results.onlySignif.dt[, NULL_ADD_BIC2 - ADD_DOM_BIC2]))
# grid.arrange(
    # qplot(
        # BIC.breaks, 
        # laply(BIC.breaks, function(x) {
            # signif <- better$DOM.coxtest[merged.results.onlySignif.dt[, NULL_ADD_BIC2 - ADD_DOM_BIC2 > x]]
            # sum(signif)/length(signif)
        # })
    # ),
    # qplot(
        # BIC.breaks, 
        # laply(BIC.breaks, function(x) {
            # signif <- better$DOM.jtest[merged.results.onlySignif.dt[, NULL_ADD_BIC2 - ADD_DOM_BIC2 > x]]
            # sum(signif)/length(signif)
        # })
    # ),
    # nrow=1, ncol=2
# )

# BIC.breaks <- seq(min(merged.results.onlySignif.dt[, NULL_ADD_BIC2 - ADD_REC_BIC2]), max(merged.results.onlySignif.dt[, NULL_ADD_BIC2 - ADD_REC_BIC2]))
# grid.arrange(
    # qplot(
        # BIC.breaks, 
        # laply(BIC.breaks, function(x) {
            # signif <- better$REC.coxtest[merged.results.onlySignif.dt[, NULL_ADD_BIC2 - ADD_REC_BIC2 > x]]
            # sum(signif)/length(signif)
        # })
    # ),
    # qplot(
        # BIC.breaks, 
        # laply(BIC.breaks, function(x) {
            # signif <- better$REC.coxtest[merged.results.onlySignif.dt[, NULL_ADD_BIC2 - ADD_REC_BIC2 > x]]
            # sum(signif)/length(signif)
        # })
    # ),
    # nrow=1, ncol=2
# )

# BIC.breaks <- seq(min(merged.results.onlySignif.dt[, NULL_ADD_BIC2 - ADD_GEN_BIC2]), max(merged.results.onlySignif.dt[, NULL_ADD_BIC2 - ADD_GEN_BIC2]))
# qplot(
    # BIC.breaks, 
    # laply(BIC.breaks, function(x) {
        # signif <- better$GEN.vuong[merged.results.onlySignif.dt[, NULL_ADD_BIC2 - ADD_GEN_BIC2 > x]]
        # sum(signif)/length(signif)
    # })
# )

# BIC.breaks <- seq(min(merged.results.onlySignif.dt[, NULL_ADD_BIC2 - ADD_HET_BIC2]), max(merged.results.onlySignif.dt[, NULL_ADD_BIC2 - ADD_HET_BIC2]))
# grid.arrange(
    # qplot(
        # BIC.breaks, 
        # laply(BIC.breaks, function(x) {
            # signif <- better$HET.coxtest[merged.results.onlySignif.dt[, NULL_ADD_BIC2 - ADD_HET_BIC2 > x]]
            # sum(signif)/length(signif)
        # })
    # ),
    # qplot(
        # BIC.breaks, 
        # laply(BIC.breaks, function(x) {
            # signif <- better$HET.jtest[merged.results.onlySignif.dt[, NULL_ADD_BIC2 - ADD_HET_BIC2 > x]]
            # sum(signif)/length(signif)
        # })
    # ),
    # nrow=1, ncol=2
# )

#
# Question 3.0
# Which snps show sufficient evidence for non additive effects?
# 
# The tests we trust should throw few positives specific to only that test.
# Using raw p value/BIC is too lenient.
# BIC and vuong tend to agree for non nested comparisons, 
# but vuong finds nothing for HET
# From the venn diagrams,
#   non-nested: BIC and vuong
#   nested: LRT
#
merged.results.dt[, better.DOM := better[, DOM.vuong & DOM.BIC]]
merged.results.dt[, better.REC := better[, REC.vuong & REC.BIC]]
merged.results.dt[, better.GEN := better[, GEN.LRT]]
merged.results.dt[, better.HET := better[, HET.vuong & HET.BIC]]

svg("3.0_venn_of_nonADD_betters.svg", width=7, height=7)

    plotVennWrapper(merged.results.dt[signif.gwas.thresh == T], c("better.DOM", "better.REC", "better.HET", "better.GEN"))

dev.off()

#
# Determine the "best model" for snps that show evidence for non additive effects
# This is the lowest BIC of significant non additive models.
#
non.add.BICs <- as.matrix(merged.results.dt[, .(ADD_DOM_BIC2, ADD_REC_BIC2, ADD_GEN_BIC2, ADD_HET_BIC2)])
non.add.pValues <- as.matrix(merged.results.dt[, .(frequentist_dom_pvalue, frequentist_rec_pvalue, frequentist_gen_pvalue, frequentist_het_pvalue)])
non.add.better <- as.matrix(merged.results.dt[, .(better.DOM, better.REC, better.GEN, better.HET)])
non.add.BICs.masked <- ifelse(non.add.better, non.add.BICs, NA)
best.models <- apply(non.add.BICs.masked, 1, function(x) {
    if(all(is.na(x))) {
        # If there is insufficient evidence for any nonADD model being better.
        return("ADD")
    } else {
        return(sub("^better.", "", names(x)[which.min(x)]))
    }
})
best.models.numeric <- as.numeric(mapvalues(best.models, from=c("ADD", "DOM", "REC", "GEN", "HET"), to=1:5))
merged.results.dt$best.model <- best.models
# Report the associated best model BIC and p value
merged.results.dt$best.BIC <- cbind(merged.results.dt$NULL_ADD_BIC2, non.add.BICs)[cbind(1:length(best.models), best.models.numeric)]
merged.results.dt$best.pValue <- cbind(merged.results.dt$frequentist_add_pvalue, non.add.pValues)[cbind(1:length(best.models), best.models.numeric)]

# Summary of which model is best
sink("misc_output.txt", append=T, split=T)
print("Tally of better models")
table(merged.results.dt$best.model)
print("Tally of better models (signif.gwas.thresh only)")
table(merged.results.dt[signif.gwas.thresh == T]$best.model)
sink()
# Distribution of how much of a BIC decrease is achieved with the best model
merged.results.dt$best.BIC.decrease <- merged.results.dt$NULL_ADD_BIC2 - merged.results.dt$best.BIC
svg("3.1_bic_decreases.png", width=7, height=7)
    hist(merged.results.dt[best.BIC.decrease > 0, best.BIC.decrease], breaks=100)
dev.off()

#
# Question 4.0
# Consider the known topsnps and the gwas3 topsnps.
# Do any of them show evidence for a non additive model being better?
#
write.fwf(
    merged.results.dt[(is.topsnp.observed | is.topsnp.meta) & best.model != "ADD" & signif.gwas.thresh],
    file = "4.0_better_nonADD_at_observed_or_meta_topsnp.txt",
    sep="\t"
)

# 
# Question 4.1
# Let each snp take the snptest p value of its best model (which may be non additive).
# Do the topsnps for the loci change?
#

# Check concordance of rankings of snps by BIC and p value
if (F) {
grid.arrange(
    qplot(rank(log10(merged.results.dt[signif.gwas.thresh == T, frequentist_add_pvalue])), rank(merged.results.dt[signif.gwas.thresh == T, NULL_ADD_BIC2])),
    qplot(rank(log10(merged.results.dt[signif.gwas.thresh == T, frequentist_dom_pvalue])), rank(merged.results.dt[signif.gwas.thresh == T, ADD_DOM_BIC2])),
    qplot(rank(log10(merged.results.dt[signif.gwas.thresh == T, frequentist_rec_pvalue])), rank(merged.results.dt[signif.gwas.thresh == T, ADD_REC_BIC2])),
    qplot(rank(log10(merged.results.dt[signif.gwas.thresh == T, frequentist_gen_pvalue])), rank(merged.results.dt[signif.gwas.thresh == T, ADD_GEN_BIC2])),
    qplot(rank(log10(merged.results.dt[signif.gwas.thresh == T, frequentist_het_pvalue])), rank(merged.results.dt[signif.gwas.thresh == T, ADD_HET_BIC2])),
    nrow=2, ncol=3
)
}

# Rank snps within loci by ADD and best model
merged.results.dt[, ADD.BIC.rank := rank(NULL_ADD_BIC2), by=locus.id]
merged.results.dt[, frequentist_add_pvalue.rank := rank(frequentist_add_pvalue), by=locus.id]
merged.results.dt[, best.BIC.rank := rank(best.BIC), by=locus.id]

# Find topsnps whose rank changes between the rankings
to.followup <- merged.results.dt[
    (frequentist_add_pvalue.rank <= 3 | best.BIC.rank <= 3) & 
    signif.gwas.thresh &
    frequentist_add_pvalue.rank != best.BIC.rank,
    .(
        alternate_ids, rsid, locus.id, position.x,
        is.topsnp.meta, is.topsnp.observed, 
        best.model, 
        frequentist_add_pvalue, best.pValue, 
        best.BIC.decrease, 
        frequentist_add_pvalue.rank, ADD.BIC.rank, best.BIC.rank, 
        rank.up=frequentist_add_pvalue.rank-best.BIC.rank,
        signif.gwas.thresh
    )
]

# Merge in rsids from bim file.
qcIMP.bim.dt <- fread('/lustre/scratch114/teams/barrett/coreex_gaibdc/QC/COMBINED/4_marker_QC/coreex_gaibdc_usgwas_qcIMP.bim')

write.fwf(
    to.followup,
    file = "4.1_better_nonADD_topsnp_appears_for_locus.txt",
    sep="\t"
)

write.fwf(
    qcIMP.bim.dt[paste(qcIMP.bim.dt[, V1], qcIMP.bim.dt[, V4]) %in% paste(to.followup[, alternate_ids], to.followup[, position.x])],
    file = "4.2_rsids_for_new_topsnps.txt",
    sep="\t"
)

# 
# Question 5.0
# Consider all snps for which there is a better non additive model.
#
# Where are they?
#

sink("misc_output.txt", append=T, split=T)
print("Tally of which loci that the significant non additive effects lie:")
table(merged.results.dt[(better.DOM | better.REC | better.HET | better.GEN) & signif.gwas.thresh, locus.id])
sink()

merged.results.dt[
    (better.DOM | better.REC | better.HET | better.GEN) & signif.gwas.thresh, 
    .(.N, mean.best.BIC.rank=mean(best.BIC.rank), mean.frequentist_add_pvalue.rank=mean(frequentist_add_pvalue.rank), mean.rank.up=mean(frequentist_add_pvalue.rank - best.BIC.rank)), 
    by=locus.id
]

#
# Output loci with a large BIC decrease
#
best.BIC.decrease.95.pc <- quantile(merged.results.dt[best.BIC.decrease > 0, best.BIC.decrease], 0.95)
write.fwf(
    merged.results.dt[
        best.BIC.decrease > best.BIC.decrease.95.pc, 
        .(
            rsid, locus.id, position.x, all_maf, info,
            is.topsnp.meta, is.topsnp.observed, 
            best.model, 
            frequentist_add_pvalue, best.pValue, 
            best.BIC.decrease, 
            ADD.BIC.rank, best.BIC.rank, 
            signif.gwas.thresh, best.BIC.decrease.95.pc.thresh=best.BIC.decrease.95.pc
        )
    ],
    file = "5.0_snps_with_best_BIC_decrease_95_pc.txt",
    sep="\t"
)

# What is the distribution of their info scores and MAFs?
svg("5.1_distr_of_info_and_maf_of_betters.svg", width=14, height=7)
    grid.arrange(
    ggplot(data=merged.results.dt[signif.gwas.thresh == T, .(better.DOM, better.REC, better.HET, better.GEN, info, all_maf)]) +
        geom_violin(aes(x=(better.DOM | better.REC | better.HET | better.GEN), y=info)),
    ggplot(data=merged.results.dt[signif.gwas.thresh == T, .(better.DOM, better.REC, better.HET, better.GEN, info, all_maf)]) +
        geom_violin(aes(x=(better.DOM | better.REC | better.HET | better.GEN), y=all_maf)),
        nrow=1, ncol=2 
    )
dev.off()

#
# Question 6.0
# Consider loci where there is previous evidence of nonadditive effects.
# e.g. Known possible recessive models @ NOD2?
#

# NOD2: containing locus 16:50511169_51009351 , topsnp 16:50745926_C_T
# Gene coords: chr 16, 50693581 to 50733081
# 3 key snps: rs2066844, rs2066845, rs2066847 @ pos c(50745926, 50756540, 50763778) 
write.fwf(
    merged.results.dt[
        (alternate_ids == 16 & position.x >= 50693581 & position.x <= 50733081 & signif.gwas.thresh) | 
        rsid == "16:50745926_C_T" |
        (alternate_ids == 16 & position.x %in% c(50745926, 50756540, 50763778)),
    ],
    file = "6.0_NOD2_snps.txt",
    sep="\t"
)
# TYK2: containing locus 19:10408439_10600418 , topsnp 19:10512911_G_A
# Gene coords: chr 19, 10350528 to 10380676
write.fwf(
    merged.results.dt[
        (alternate_ids == 19 & position.x >= 10350528 & position.x <= 10380676 & signif.gwas.thresh) | 
        rsid == "19:10512911_G_A",
    ],
    file = "6.1_TYK2_snps.txt",
    sep="\t"
)

#
# TODO unfinished things below here
#
stopifnot(F)

# Are the positives enriched for snps with high/low INFO?
# Plot non add p value against info
par(mfrow=c(2, 2))
plot(rank(log10(merged.results.onlySignif.dt[best.model == "DOM", frequentist_dom_pvalue])), rank(merged.results.onlySignif.dt[best.model == "DOM", frequentist_dom_info]))
plot(rank(log10(merged.results.onlySignif.dt[best.model == "REC", frequentist_rec_pvalue])), rank(merged.results.onlySignif.dt[best.model == "REC", frequentist_rec_info]))
plot(rank(log10(merged.results.onlySignif.dt[best.model == "GEN", frequentist_gen_pvalue])), rank(merged.results.onlySignif.dt[best.model == "GEN", frequentist_gen_info]))
plot(rank(log10(merged.results.onlySignif.dt[best.model == "HET", frequentist_het_pvalue])), rank(merged.results.onlySignif.dt[best.model == "HET", frequentist_het_info]))

# Kendall rank correlation coefficient, commonly referred to as Kendall's tau coefficient (after the Greek letter τ), is a statistic used to measure the ordinal association between two measured quantities. 
cor.test(log10(merged.results.onlySignif.dt[best.model == "DOM", frequentist_dom_pvalue]), merged.results.onlySignif.dt[best.model == "DOM", frequentist_dom_info], method="kendall")
cor.test(log10(merged.results.onlySignif.dt[best.model == "REC", frequentist_rec_pvalue]), merged.results.onlySignif.dt[best.model == "REC", frequentist_rec_info], method="kendall")
cor.test(log10(merged.results.onlySignif.dt[best.model == "GEN", frequentist_gen_pvalue]), merged.results.onlySignif.dt[best.model == "GEN", frequentist_gen_info], method="kendall")
cor.test(log10(merged.results.onlySignif.dt[best.model == "HET", frequentist_het_pvalue]), merged.results.onlySignif.dt[best.model == "HET", frequentist_het_info], method="kendall")

#
# Is the known meta-analysis topsnp per locus changed by allowing for non additive models?
#

(merged.results.onlySignif.dt[better.DOM | better.REC | better.HET | better.GEN, .(rsid, locus.id)])

merged.results.onlySignif.dt[is.topsnp.meta == T, ]
merged.results.onlySignif.dt[is.topsnp.observed == T, ]

x <- "8_49047317_49206289"
merged.results.onlySignif.dt[locus.id == x & ADD.BIC.rank < 20, ]

# Within loci, do the rankings change significantly comparing additive snptest vs additive LRT p values/BIC
# If not, do the rankings change significantly comparing additive BIC vs best BIC

merged.results.dt[
    ,
    # kendall(cbind(frequentist_add_pvalue, ADD_pValue.LRT))$value,
    kendall(cbind(frequentist_add_pvalue, ADD_pValue.LRT), correct=T)$value,
    by=locus.id
]

with(merged.results.dt,{
    kendall(cbind(frequentist_add_pvalue, ADD_pValue.LRT), correct=T)
})

x <- "8_49047317_49206289"
kendall(cbind(merged.results.onlySignif.dt[locus.id == x, frequentist_add_pvalue], merged.results.onlySignif.dt[locus.id == x, ADD_pValue.LRT]), correct=T)
kendall(cbind(merged.results.onlySignif.dt[locus.id == x, ADD_pValue.LRT], merged.results.onlySignif.dt[locus.id == x, NULL_ADD_BIC2]), correct=T)
kendall(cbind(merged.results.onlySignif.dt[locus.id == x, NULL_ADD_BIC2], merged.results.onlySignif.dt[locus.id == x, best.BIC]), correct=T)

merged.results.onlySignif.dt[locus.id == x & ADD.BIC.rank < 20]

#
# Possiblities
#
    # Interaction effects?
    # Why does interaction appear as dominant
    #
    # Implication of loci in UC/CD that were not implicated before due to wrong model choice?
    #
    # Bayesian model comparison

