#!/software/R-3.3.0/bin/Rscript
#
library(data.table)
library(gtools)
library(plyr)
library(ggplot2)
library(gplots)
library(gridExtra)
library(irr)

dataset <- "gwas3"
assoc <- "ibd"

# Read in snptest results for loci that were re run in glm
contained.snps.dt <- fread(sprintf("/nfs/users/nfs_b/bb9/workspace/rotation1/crohns_workspace/4_gwas/2_R_glm/contained.snps.%s.%s.all.txt", dataset, assoc))

# QC checks
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
merged.results.dt[, is.topsnp.known := topSNP.Position..bp. == position.x, by=locus.id]
merged.results.dt[, is.topsnp.observed := frequentist_add_pvalue == min(frequentist_add_pvalue, na.rm=T), by=locus.id]

#
# Start comparisons
#

# pdf("1.0_compare_models.pdf")

# Can we trust our analysis techniques?

# Coefficient p value is almost the same p value as the lrt
ggplot(data=merged.results.dt, aes(x=ADD_pValue.LRT, y=ADD_pValue.coef)) + 
    geom_point(aes(colour=frequentist_add_info), size=0.5) +
    scale_x_log10() + scale_y_log10()

# Check that robust LRT is concordant with traditional LRT
ggplot(data=merged.results.dt[NULL_ADD_pValue.vuongtest.LRT > 0], aes(x=NULL_ADD_pValue.vuongtest.LRT, y=ADD_pValue.LRT)) +
    geom_point(aes(colour=frequentist_add_info), size=0.1) +
    scale_x_log10() + scale_y_log10() 

# @@@@@
# LRT value corresponds well to snptest values with high INFO
grid.arrange(
    # For ADD
    ggplot(data=merged.results.dt, aes(x=frequentist_add_pvalue, y=ADD_pValue.LRT)) +
        geom_point(aes(colour=frequentist_add_info), size=0.5) +
        geom_abline(slope=1, intercept=0) +
        scale_colour_gradient(low="red") +
        scale_x_log10() + scale_y_log10(),
        # geom_text(data=merged.results.dt[abs(1 - frequentist_add_pvalue/ADD_pValue.LRT) > 0.99], aes(label=round(frequentist_add_info, digits=2)), size=3),
    # For GEN
    ggplot(data=merged.results.dt, aes(x=frequentist_gen_pvalue, y=GEN_pValue.LRT)) +
        geom_point(aes(colour=frequentist_gen_info), size=0.5) +
        geom_abline(slope=1, intercept=0) +
        scale_colour_gradient(low="red") +
        scale_x_log10() + scale_y_log10(), 
    nrow=1, ncol=2 
)

#
# Detecting evidence for non additive models being better
#
better <- data.table(rsid=merged.results.dt$rsid)

# Determine the significance threshold to use
# Bonferroni on the number of 'independent' loci
# TODO, consider subsetting, then holm for more power
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
better[, DOM.raw.BIC := merged.results.dt[, ADD_DOM_BIC2 < NULL_ADD_BIC2]]
better[, REC.raw.BIC := merged.results.dt[, ADD_REC_BIC2 < NULL_ADD_BIC2]]
better[, GEN.raw.BIC := merged.results.dt[, ADD_GEN_BIC2 < NULL_ADD_BIC2]]
better[, HET.raw.BIC := merged.results.dt[, ADD_HET_BIC2 < NULL_ADD_BIC2]]

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

# A decision based on BIC (equiv. to AIC)
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
dim(merged.results.dt[signif.gwas.thresh == T])

#
# Question 0
# Basic numbers on the original and filtered analysis.
# Consider the known topsnps and the gwas3 topsnps.
# How concordant are the two sets?
#
# How many known loci remain?
length(unique(merged.results.dt$locus.id))
# How many known topsnps remain?
sum(merged.results.dt$is.topsnp.known)
#
# How many known assoc specific loci remain?
length(unique(merged.results.dt[grepl(toupper(assoc), Trait), locus.id]))
# How many known assoc specific topsnps remain?
sum(merged.results.dt[grepl(toupper(assoc), Trait), is.topsnp.known])
#
# How many known loci reached signif?
length(unique(merged.results.dt[signif.gwas.thresh == T, locus.id]))
# How many known topsnps reached signif?
sum(merged.results.dt[signif.gwas.thresh == T, is.topsnp.known])
# How many gwas3 topsnps reached signif?
sum(merged.results.dt[signif.gwas.thresh == T, is.topsnp.observed])
# Distribution of signif snps per locus
# table(merged.results.onlySignif.dt$locus.id)

plotVennWrapper <- function(dt, names) {
    subset.dt <- dt[, names, with=F]
    subset.dt[is.na(subset.dt)] <- F
    v <- venn(subset.dt, show.plot=T)
}
# TODO

#
# Question 0.5
# How should we judge if a snp displays sufficient evidence for a non additive model?
#

better.onlySignif <- better[merged.results.dt$signif.gwas.thresh]

# Compare models
plotVennWrapper(better.onlySignif, c("DOM.raw.pValue", "REC.raw.pValue", "GEN.raw.pValue", "HET.raw.pValue"))
plotVennWrapper(better.onlySignif, c("DOM.raw.BIC", "REC.raw.BIC", "GEN.raw.BIC", "HET.raw.BIC"))
plotVennWrapper(better.onlySignif, c("DOM.BIC", "REC.BIC", "HET.BIC"))
plotVennWrapper(better.onlySignif, c("DOM.vuong", "REC.vuong", "HET.vuong"))
plotVennWrapper(better.onlySignif, c("DOM.coxtest", "REC.coxtest", "HET.coxtest"))
plotVennWrapper(better.onlySignif, c("DOM.jtest", "REC.jtest", "HET.jtest"))

# Same model, compare methods
plotVennWrapper(better.onlySignif, c("DOM.raw.pValue", "DOM.vuong", "DOM.BIC", "DOM.coxtest", "DOM.jtest"))
plotVennWrapper(better.onlySignif, c("REC.raw.pValue", "REC.vuong", "REC.BIC", "REC.coxtest", "REC.jtest"))
plotVennWrapper(better.onlySignif, c("GEN.raw.pValue", "GEN.raw.BIC", "GEN.vuong", "GEN.LRT"))
plotVennWrapper(better.onlySignif, c("HET.raw.pValue", "HET.vuong", "HET.BIC", "HET.coxtest", "HET.jtest"))

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

# @@@@@
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
merged.results.dt[, better.HET := better[, HET.BIC]]

# Determine the "best model" for snps that show evidence for non additive effects
# This is the lowest BIC of significant non additive models.
# Report this lowest BIC.
#
# best.BIC.models <- alply(merged.results.dt, 1, function(x) {
    # modelNames <- c("DOM", "REC", "GEN", "HET")
    # BICs <- with(x, c(ADD_DOM_BIC2, ADD_REC_BIC2, ADD_GEN_BIC2, ADD_HET_BIC2))
    # betters <- with(x, c(better.DOM, better.REC, better.GEN, better.HET))
    # if(length(bestBIC)) {
        # result <- min(bestBIC)
        # names(result) <- modelNames[which(BICs == min(bestBIC))][1]
    # } else {
        # result <- x$NULL_ADD_BIC2
        # names(result) <- "ADD"
    # }
    # return(result)
# }, .progress="text")

# Get best non additive model, BIC, and evidence for non additivity with that model
non.add.BICs <- merged.results.dt[, .(ADD_DOM_BIC2, ADD_REC_BIC2, ADD_GEN_BIC2, ADD_HET_BIC2)]
min.BIC.col <- max.col(-non.add.BICs, ties.method="first")
min.BIC.model <- c("DOM", "REC", "GEN", "HET")[min.BIC.col]
min.BIC.BICs <- data.frame(merged.results.dt[, .(ADD_DOM_BIC2, ADD_REC_BIC2, ADD_GEN_BIC2, ADD_HET_BIC2)])[cbind(1:nrow(merged.results.dt), min.BIC.col)]
min.BIC.pValues <- data.frame(merged.results.dt[, .(frequentist_add_pvalue, frequentist_rec_pvalue, frequentist_gen_pvalue, frequentist_het_pvalue)])[cbind(1:nrow(merged.results.dt), min.BIC.col)]
min.BIC.better <- data.frame(merged.results.dt[, .(better.DOM, better.REC, better.GEN, better.HET)])[cbind(1:nrow(merged.results.dt), min.BIC.col)]

# Set best models for each snp
merged.results.dt$best.model <- "ADD"
merged.results.dt$best.model[min.BIC.better] <- min.BIC.model[min.BIC.better]
merged.results.dt$best.BIC <- merged.results.dt$NULL_ADD_BIC2
merged.results.dt$best.BIC[min.BIC.better] <- min.BIC.BICs[min.BIC.better]
merged.results.dt$best.pValue <- merged.results.dt$frequentist_add_pvalue
merged.results.dt$best.pValue[min.BIC.better] <- min.BIC.pValues[min.BIC.better]
merged.results.dt$best.BIC.decrease <- merged.results.dt$NULL_ADD_BIC2 - merged.results.dt$best.BIC

# @@@@@
# Summary of which model is best
table(merged.results.dt$best.model)
# Distribution of how much of a BIC decrease is achieved with the best model
hist(merged.results.dt[best.BIC.decrease > 0, best.BIC.decrease], breaks=50)

#
# Question 1
# Consider the known topsnps and the gwas3 topsnps.
# Do any of them show evidence for a non additive model being better?
#
merged.results.dt[is.topsnp.known == T & best.model != "ADD"]
merged.results.dt[is.topsnp.observed == T & best.model != "ADD"]

# What if we require genome wide significance?
merged.results.dt[is.topsnp.known == T & best.model != "ADD" & signif.gwas.thresh]
merged.results.dt[is.topsnp.observed == T & best.model != "ADD" & signif.gwas.thresh]
# 161479745_A_G rs1801274

# 
# Question 2
# Let each snp take the snptest p value of its best model (which may be non additive).
# Do the topsnps for the loci change?
#

# Rankings of snps by BIC and p value are concordant
grid.arrange(
    qplot(log10(merged.results.dt[signif.gwas.thresh == T, frequentist_add_pvalue]), merged.results.dt[signif.gwas.thresh == T, NULL_ADD_BIC2]),
    qplot(log10(merged.results.dt[signif.gwas.thresh == T, frequentist_dom_pvalue]), merged.results.dt[signif.gwas.thresh == T, ADD_DOM_BIC2]),
    qplot(log10(merged.results.dt[signif.gwas.thresh == T, frequentist_rec_pvalue]), merged.results.dt[signif.gwas.thresh == T, ADD_REC_BIC2]),
    qplot(log10(merged.results.dt[signif.gwas.thresh == T, frequentist_gen_pvalue]), merged.results.dt[signif.gwas.thresh == T, ADD_GEN_BIC2]),
    qplot(log10(merged.results.dt[signif.gwas.thresh == T, frequentist_het_pvalue]), merged.results.dt[signif.gwas.thresh == T, ADD_HET_BIC2]),
    nrow=2, ncol=3
)

# Rank snps within loci by ADD and best model
merged.results.dt[, ADD.BIC.rank := rank(NULL_ADD_BIC2), by=locus.id]
merged.results.dt[, best.BIC.rank := rank(best.BIC), by=locus.id]

# Find topsnps whose rank changes between the rankings
merged.results.dt[
    (ADD.BIC.rank == 1 | best.BIC.rank == 1) & 
    signif.gwas.thresh &
    ADD.BIC.rank != best.BIC.rank,
    .(
        rsid, locus.id, 
        is.topsnp.known, is.topsnp.observed, 
        best.model, 
        frequentist_add_pvalue, best.pValue, 
        best.BIC.decrease, 
        ADD.BIC.rank, best.BIC.rank, 
        signif.gwas.thresh
    )
]

# There are 2 snps in the output.
# One is a known topsnp.
# In both cases, in gwas3:
# - there is enough evidence to suggest the best model is not ADD
# - there is a large BIC decrease using the best model
# - their frequentist_add_pvalue is not genome wide significant, but their best p value is
# - they are not the lead snp with ranked using frequentist_add_pvalue, but they are when ranked using BIC
# 

# 
# Question 3
# Consider all snps for which there is a better non additive model.
# What is the distribution of their info scores and MAFs?
#

#
# TODO below here
#
dev.off()
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

# @@@@@
(merged.results.onlySignif.dt[better.DOM | better.REC | better.HET | better.GEN, .(rsid, locus.id)])

# @@@@@
# Loci that the non additive effects lie in
unique(merged.results.onlySignif.dt[better.DOM | better.REC | better.HET | better.GEN, locus.id])
# [1] "1_161460211_161638410" "6_31011373_32778656"   "8_49047317_49206289"

merged.results.onlySignif.dt[is.topsnp.known == T, ]
merged.results.onlySignif.dt[is.topsnp.observed == T, ]

# @@@@@
x <- "8_49047317_49206289"
merged.results.onlySignif.dt[locus.id == x & ADD.BIC.rank < 20, ]

# Within loci, do the rankings change significantly comparing additive snptest vs additive LRT p values/BIC
# If not, do the rankings change significantly comparing additive BIC vs best BIC

#
# Known possible recessive models?
#

# @@@@@
# NOD2: containing locus 16:50511169_51009351 , topsnp 16:50745926_C_T
# Gene coords: chr 16, 50693581 to 50733081
better[rsid ==  "16:50745926_C_T", ]
NOD2.rsid <- merged.results.dt[alternate_ids == 16 & position.x >= 50693581 & position.x <= 50733081, rsid]
better[rsid %in% NOD2.rsid, ]

# TYK2: containing locus 19:10408439_10600418 , topsnp 19:10512911_G_A
# Gene coords: chr 19, 10350528 to 10380676
better[rsid ==  "19:10512911_G_A", ]
TYK2.rsid <- merged.results.dt[alternate_ids == 19 & position.x >= 10350528 & position.x <= 10380676, rsid]
better[rsid %in% TYK2.rsid, ]

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

# Interaction effects?
# Why does interaction appear as dominant

# Implication of loci in UC/CD that were not implicated before due to wrong model choice?

# Bayesian model comparison

