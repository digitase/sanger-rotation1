#!/software/R-3.3.0/bin/Rscript
#
library(data.table)
library(gtools)
library(plyr)
library(ggplot2)
library(gplots)
library(gridExtra)
library(irr)

# 
# Read in snptest results for loci that were re run in glm
# contained.snps.dt <- fread("/nfs/users/nfs_b/bb9/workspace/rotation1/crohns_workspace/4_gwas/2_R_glm/contained.snps.signif.txt")
contained.snps.dt <- fread("/nfs/users/nfs_b/bb9/workspace/rotation1/crohns_workspace/4_gwas/2_R_glm/contained.snps.all.txt")

# Read in results from glm
all.glm.out.files <- mixedsort(list.files("/nfs/users/nfs_b/bb9/workspace/rotation1/crohns_workspace/4_gwas/2_R_glm/chunks/gwas3/ibd/", include.dirs=F, recursive=T, full.names=T))
length(all.glm.out.files)
glm.results.dt <- data.table(ldply(all.glm.out.files, fread, .progress="text"))

# Merge in snptest p values
merged.results.dt <- merge(contained.snps.dt, glm.results.dt, by.x=c("rsid"), by.y=c("snp.id"), all.y=T)
dim(merged.results.dt)
# Find topsnps
merged.results.dt[, is.topsnp.known := topSNP.Position..bp. == position.x, by=locus.id]
merged.results.dt[, is.topsnp.observed := frequentist_add_pvalue == min(frequentist_add_pvalue, na.rm=T), by=locus.id]
# Get model number for lowest p value for each snp
merged.results.dt[, lowest.pValue.modelNumber := which.min(c(frequentist_add_pvalue, frequentist_dom_pvalue, frequentist_rec_pvalue, frequentist_gen_pvalue, frequentist_het_pvalue)), by=1:nrow(merged.results.dt)]

#
# Start comparisons
#

pdf("plots.pdf")

# Can we trust our analysis techniques?

# Coefficient p value is NOT the same p value as the lrt
ggplot(data=merged.results.dt, aes(x=ADD_pValue.LRT, y=ADD_pValue.coef)) + 
    geom_point(aes(colour=frequentist_add_info), size=0.1)

# @@@@@
# LRT value corresponds well to snptest values with high INFO
# For ADD
ggplot(data=merged.results.dt, aes(x=frequentist_add_pvalue, y=ADD_pValue.LRT)) +
    geom_point(aes(colour=frequentist_add_info), size=0.1) +
    scale_colour_gradient(low="red") 
    # scale_x_log10() + scale_y_log10() +
    # geom_text(data=merged.results.dt[abs(1 - frequentist_add_pvalue/ADD_pValue.LRT) > 0.99], aes(label=round(frequentist_add_info, digits=2)), size=3) 

# For GEN
summary(merged.results.dt$frequentist_gen_pvalue)
summary(merged.results.dt$GEN_pValue.LRT)
ggplot(data=merged.results.dt, aes(x=frequentist_gen_pvalue, y=GEN_pValue.LRT)) +
    geom_point(aes(colour=frequentist_gen_info), size=0.1) +
    scale_colour_gradient(low="red") 
    # scale_x_log10() + scale_y_log10() +
    # geom_text(data=merged.results.dt[abs(1 - frequentist_gen_pvalue/GEN_pValue.LRT) > 0.99], aes(label=round(frequentist_gen_info, digits=2)), size=3) 

# @@@@@
# Check that robust LRT is concordant with traditional LRT
ggplot(data=merged.results.dt, aes(x=NULL_ADD_pValue.vuongtest.LRT, y=ADD_pValue.LRT)) +
    geom_point(aes(colour=frequentist_add_info), size=0.1)

# Prefilter data to leave only genome wide significant snps 
gwas.thresh <- 5e-8
signif.thresh <- 0.05

# Prefilter on snptest p values from score test
merged.results.onlySignif.dt <- merged.results.dt[
    frequentist_add_pvalue < gwas.thresh |
    frequentist_dom_pvalue < gwas.thresh |
    frequentist_rec_pvalue < gwas.thresh |
    frequentist_gen_pvalue < gwas.thresh |
    frequentist_het_pvalue < gwas.thresh
]
dim(merged.results.onlySignif.dt)

# Check the power of the cohort analysis
# How many known IBD loci remain?
length(unique(contained.snps.dt$locus.id))
length(unique(merged.results.dt$locus.id))
length(unique(merged.results.onlySignif.dt$locus.id))
# How many of the corresponding topsnps remain?
nrow(contained.snps.dt[topSNP.Position..bp. == position])
sum(merged.results.dt$is.topsnp.known)
sum(merged.results.dt$is.topsnp.observed, na.rm=T) # TODO why are there multiple best snps
sum(merged.results.onlySignif.dt$is.topsnp.known)
sum(merged.results.onlySignif.dt$is.topsnp.observed)

# Distribution of snps per chromosome...
table(merged.results.onlySignif.dt$alternate_ids)
# ... and per locus
table(merged.results.onlySignif.dt$locus.id)

# Get hist of proportion of lowest p values are from each type of model
hist(merged.results.onlySignif.dt$lowest.pValue.modelNumber)

#
# Detecting snps where non additive models are better
#

# Apply FDR corrections to relevant cols
merged.results.onlySignif.dt[, ADD_DOM_pValue.coxtest.FDR := p.adjust(ADD_DOM_pValue.coxtest, method="holm")]
merged.results.onlySignif.dt[, ADD_REC_pValue.coxtest.FDR := p.adjust(ADD_REC_pValue.coxtest, method="holm")]
merged.results.onlySignif.dt[, ADD_HET_pValue.coxtest.FDR := p.adjust(ADD_HET_pValue.coxtest, method="holm")]
merged.results.onlySignif.dt[, ADD_DOM_pValue.coxtest.converse.FDR := p.adjust(ADD_DOM_pValue.coxtest.converse, method="holm")]
merged.results.onlySignif.dt[, ADD_REC_pValue.coxtest.converse.FDR := p.adjust(ADD_REC_pValue.coxtest.converse, method="holm")]
merged.results.onlySignif.dt[, ADD_HET_pValue.coxtest.converse.FDR := p.adjust(ADD_HET_pValue.coxtest.converse, method="holm")]
#
merged.results.onlySignif.dt[, ADD_DOM_pValue.jtest.FDR := p.adjust(ADD_DOM_pValue.jtest, method="holm")]
merged.results.onlySignif.dt[, ADD_REC_pValue.jtest.FDR := p.adjust(ADD_REC_pValue.jtest, method="holm")]
merged.results.onlySignif.dt[, ADD_HET_pValue.jtest.FDR := p.adjust(ADD_HET_pValue.jtest, method="holm")]
merged.results.onlySignif.dt[, ADD_DOM_pValue.jtest.converse.FDR := p.adjust(ADD_DOM_pValue.jtest.converse, method="holm")]
merged.results.onlySignif.dt[, ADD_REC_pValue.jtest.converse.FDR := p.adjust(ADD_REC_pValue.jtest.converse, method="holm")]
merged.results.onlySignif.dt[, ADD_HET_pValue.jtest.converse.FDR := p.adjust(ADD_HET_pValue.jtest.converse, method="holm")]
#
merged.results.onlySignif.dt[, ADD_DOM_pValue.vuongtest.varTest.FDR := p.adjust(ADD_DOM_pValue.vuongtest.varTest, method="holm")]
merged.results.onlySignif.dt[, ADD_REC_pValue.vuongtest.varTest.FDR := p.adjust(ADD_REC_pValue.vuongtest.varTest, method="holm")]
merged.results.onlySignif.dt[, ADD_GEN_pValue.vuongtest.varTest.FDR := p.adjust(ADD_GEN_pValue.vuongtest.varTest, method="holm")]
merged.results.onlySignif.dt[, ADD_HET_pValue.vuongtest.varTest.FDR := p.adjust(ADD_HET_pValue.vuongtest.varTest, method="holm")]
#
merged.results.onlySignif.dt[, ADD_DOM_pValue.vuongtest.LRT.FDR := p.adjust(ADD_DOM_pValue.vuongtest.LRT, method="none")]
merged.results.onlySignif.dt[, ADD_REC_pValue.vuongtest.LRT.FDR := p.adjust(ADD_REC_pValue.vuongtest.LRT, method="none")]
merged.results.onlySignif.dt[, ADD_GEN_pValue.vuongtest.LRT.FDR := p.adjust(ADD_GEN_pValue.vuongtest.LRT, method="none")]
merged.results.onlySignif.dt[, ADD_HET_pValue.vuongtest.LRT.FDR := p.adjust(ADD_HET_pValue.vuongtest.LRT, method="none")]
#
merged.results.onlySignif.dt[, ADD_GEN_pValue.LRT.FDR := p.adjust(ADD_GEN_pValue.LRT, method="holm")]

# Ensure rsids are unique in this subset
stopifnot(length(merged.results.onlySignif.dt$rsid) == length(unique(merged.results.onlySignif.dt$rsid)))

#
# A decision based on raw p value
#
# Find snps with better non additive p values
better <- data.table(rsid=merged.results.onlySignif.dt$rsid)
better[, DOM.raw.pValue := merged.results.onlySignif.dt[, frequentist_dom_pvalue < frequentist_add_pvalue]]
better[, REC.raw.pValue := merged.results.onlySignif.dt[, frequentist_rec_pvalue < frequentist_add_pvalue]]
better[, GEN.raw.pValue := merged.results.onlySignif.dt[, frequentist_gen_pvalue < frequentist_add_pvalue]]
better[, HET.raw.pValue := merged.results.onlySignif.dt[, frequentist_het_pvalue < frequentist_add_pvalue]]

#
# Decision methods for non nested models
#
# Compare concordance of vuong LRT and BIC confidence intervals for non nested models
# More difficult to apply FDR to ci, so prefer LRT if possible

# Investigate dodgy negative p values from vartest
# Has been observed in conjunction with this error:
# 5: In imhof(n * omega.hat.2, lamstar^2) :
  # Note that Qq + abserr is positive.
#
# Is it the case that a negative p value always leads to failure to reject the null via the BIC criterion?
length(
    intersect(
        merged.results.onlySignif.dt[ADD_DOM_pValue.vuongtest.varTest < 0, rsid],
        merged.results.onlySignif.dt[ADD_DOM_BICci.lower * ADD_DOM_BICci.upper > 0 & ADD_DOM_BICci.lower > 0, rsid]
    )
) > 0
#
# This is not the case, thus better to just exclude negative p values.

# Use vuong 2 step test (non nested version)
better[
    , 
    DOM.vuong := merged.results.onlySignif.dt[
        , 
        ADD_DOM_pValue.vuongtest.varTest.FDR > 0 & 
        ADD_DOM_pValue.vuongtest.varTest.FDR < signif.thresh &
        ADD_DOM_pValue.vuongtest.LRT.FDR > 0 &
        ADD_DOM_pValue.vuongtest.LRT.FDR < signif.thresh
    ]
]
better[
    , 
    REC.vuong := merged.results.onlySignif.dt[
        , 
        ADD_REC_pValue.vuongtest.varTest.FDR > 0 &
        ADD_REC_pValue.vuongtest.varTest.FDR < signif.thresh &
        ADD_REC_pValue.vuongtest.LRT.FDR > 0 &
        ADD_REC_pValue.vuongtest.LRT.FDR < signif.thresh
    ]
]
better[
    , 
    HET.vuong := merged.results.onlySignif.dt[
        , 
        ADD_HET_pValue.vuongtest.varTest.FDR > 0 &
        ADD_HET_pValue.vuongtest.varTest.FDR < signif.thresh &
        ADD_HET_pValue.vuongtest.LRT.FDR > 0 &
        ADD_HET_pValue.vuongtest.LRT.FDR < signif.thresh
    ]
]

# A decision based on BIC (equiv. to AIC)
# Note: if models are nested or if the "variance test" from
     # ‘vuongtest()’ indicates models are indistinguishable, then the
     # intervals returned from ‘icci()’ will be incorrect.
# Requires models to be distinguishable
# Requires 0 to not be in the ci, and the decrease to be positive
better[, DOM.BIC := merged.results.onlySignif.dt[, ADD_DOM_pValue.vuongtest.varTest.FDR > 0 & ADD_DOM_pValue.vuongtest.varTest.FDR < signif.thresh & ADD_DOM_BICci.lower * ADD_DOM_BICci.upper > 0 & ADD_DOM_BICci.lower > 0 ]]
better[, REC.BIC := merged.results.onlySignif.dt[, ADD_REC_pValue.vuongtest.varTest.FDR > 0 & ADD_REC_pValue.vuongtest.varTest.FDR < signif.thresh & ADD_REC_BICci.lower * ADD_REC_BICci.upper > 0 & ADD_REC_BICci.lower > 0 ]]
better[, HET.BIC := merged.results.onlySignif.dt[, ADD_HET_pValue.vuongtest.varTest.FDR > 0 & ADD_HET_pValue.vuongtest.varTest.FDR < signif.thresh & ADD_HET_BICci.lower * ADD_HET_BICci.upper > 0 & ADD_HET_BICci.lower > 0 ]]
# A decision based on coxtest
better[, DOM.coxtest := merged.results.onlySignif.dt[, ADD_DOM_pValue.coxtest.FDR > 0 & ADD_DOM_pValue.coxtest.FDR < signif.thresh & ADD_DOM_pValue.coxtest.converse.FDR > signif.thresh]]
better[, REC.coxtest := merged.results.onlySignif.dt[, ADD_REC_pValue.coxtest.FDR > 0 & ADD_REC_pValue.coxtest.FDR < signif.thresh & ADD_REC_pValue.coxtest.converse.FDR > signif.thresh]]
better[, HET.coxtest := merged.results.onlySignif.dt[, ADD_HET_pValue.coxtest.FDR > 0 & ADD_HET_pValue.coxtest.FDR < signif.thresh & ADD_HET_pValue.coxtest.converse.FDR > signif.thresh]]
# A decision based on jtest
better[, DOM.jtest := merged.results.onlySignif.dt[, ADD_DOM_pValue.jtest.FDR > 0 & ADD_DOM_pValue.jtest.FDR < signif.thresh & ADD_DOM_pValue.jtest.converse.FDR > signif.thresh]]
better[, REC.jtest := merged.results.onlySignif.dt[, ADD_REC_pValue.jtest.FDR > 0 & ADD_REC_pValue.jtest.FDR < signif.thresh & ADD_REC_pValue.jtest.converse.FDR > signif.thresh]]
better[, HET.jtest := merged.results.onlySignif.dt[, ADD_HET_pValue.jtest.FDR > 0 & ADD_HET_pValue.jtest.FDR < signif.thresh & ADD_HET_pValue.jtest.converse.FDR > signif.thresh]]

#
# Decision methods for ADD vs GEN (nested)
#
# Use the traditional LRT
better[, GEN.LRT := merged.results.onlySignif.dt[, ADD_GEN_pValue.LRT.FDR < signif.thresh]]
#
# Use vuong 2 step test (nested version)
#
better[
    , 
    GEN.vuong := merged.results.onlySignif.dt[
        , 
        ADD_GEN_pValue.vuongtest.varTest.FDR > 0 &
        ADD_GEN_pValue.vuongtest.varTest.FDR < signif.thresh &
        ADD_GEN_pValue.vuongtest.LRT.FDR > 0 &
        ADD_GEN_pValue.vuongtest.LRT.FDR < signif.thresh
    ]
]

#
# Venn diagrams
#

plotVennWrapper <- function(dt, names) {
    subset.dt <- dt[, names, with=F]
    subset.dt[is.na(subset.dt)] <- F
    v <- venn(subset.dt, show.plot=T)
}

# @@@@@
# Compare models

plotVennWrapper(better, c("DOM.raw.pValue", "REC.raw.pValue", "GEN.raw.pValue", "HET.raw.pValue"))

plotVennWrapper(better, c("DOM.BIC", "REC.BIC", "HET.BIC"))

plotVennWrapper(better, c("DOM.vuong", "REC.vuong", "HET.vuong"))

plotVennWrapper(better, c("DOM.coxtest", "REC.coxtest", "HET.coxtest"))

plotVennWrapper(better, c("DOM.jtest", "REC.jtest", "HET.jtest"))

# Same model, compare methods

plotVennWrapper(better, c("DOM.raw.pValue", "DOM.vuong", "DOM.BIC", "DOM.coxtest", "DOM.jtest"))

plotVennWrapper(better, c("REC.raw.pValue", "REC.vuong", "REC.BIC", "REC.coxtest", "REC.jtest"))

plotVennWrapper(better, c("GEN.raw.pValue", "GEN.vuong", "GEN.LRT"))

plotVennWrapper(better, c("HET.raw.pValue", "HET.vuong", "HET.BIC", "HET.coxtest", "HET.jtest"))

# @@@@@
# Heatmap of method concordance
#
better.molten <- melt(better[order(rowSums(better[, -1, with=F], na.rm=T)), sort(names(better)), with=F], id.vars=c("rsid"))
ggplot(better.molten, aes(rsid, variable)) + 
    geom_tile(aes(fill = value)) + theme(axis.text.x = element_blank())

#
# Investigate properties of coxtest and j test
#
# Plot coxtest and j test p values against decrease in BIC
par(mfrow=c(1, 3))
plot(merged.results.onlySignif.dt$NULL_ADD_BIC2 - merged.results.onlySignif.dt$ADD_DOM_BIC2, merged.results.onlySignif.dt$ADD_DOM_pValue.coxtest.FDR)
plot(merged.results.onlySignif.dt$NULL_ADD_BIC2 - merged.results.onlySignif.dt$ADD_REC_BIC2, merged.results.onlySignif.dt$ADD_REC_pValue.coxtest.FDR)
plot(merged.results.onlySignif.dt$NULL_ADD_BIC2 - merged.results.onlySignif.dt$ADD_HET_BIC2, merged.results.onlySignif.dt$ADD_HET_pValue.coxtest.FDR)

par(mfrow=c(1, 3))
plot(merged.results.onlySignif.dt$NULL_ADD_BIC2 - merged.results.onlySignif.dt$ADD_DOM_BIC2, merged.results.onlySignif.dt$ADD_DOM_pValue.jtest.FDR)
plot(merged.results.onlySignif.dt$NULL_ADD_BIC2 - merged.results.onlySignif.dt$ADD_REC_BIC2, merged.results.onlySignif.dt$ADD_REC_pValue.jtest.FDR)
plot(merged.results.onlySignif.dt$NULL_ADD_BIC2 - merged.results.onlySignif.dt$ADD_HET_BIC2, merged.results.onlySignif.dt$ADD_HET_pValue.jtest.FDR)

# @@@@@
# Plot % significant as decided by the test for different bic decreases
# Poor performance at high BICs: few with that high of a BIC all get rejected
merged.results.onlySignif.dt[NULL_ADD_BIC2 - ADD_HET_BIC2 > 45]
#
BIC.breaks <- seq(min(merged.results.onlySignif.dt[, NULL_ADD_BIC2 - ADD_DOM_BIC2]), max(merged.results.onlySignif.dt[, NULL_ADD_BIC2 - ADD_DOM_BIC2]))
grid.arrange(
    qplot(
        BIC.breaks, 
        laply(BIC.breaks, function(x) {
            # signif <- merged.results.onlySignif.dt[NULL_ADD_BIC2 - ADD_DOM_BIC2 > x, ADD_DOM_pValue.coxtest.FDR] < signif.thresh
            signif <- better$DOM.coxtest[merged.results.onlySignif.dt[, NULL_ADD_BIC2 - ADD_DOM_BIC2 > x]]
            sum(signif)/length(signif)
        })
    ),
    qplot(
        BIC.breaks, 
        laply(BIC.breaks, function(x) {
            # signif <- merged.results.onlySignif.dt[NULL_ADD_BIC2 - ADD_DOM_BIC2 > x, ADD_DOM_pValue.jtest.FDR] < signif.thresh
            signif <- better$DOM.jtest[merged.results.onlySignif.dt[, NULL_ADD_BIC2 - ADD_DOM_BIC2 > x]]
            sum(signif)/length(signif)
        })
    ),
    nrow=1, ncol=2
)

BIC.breaks <- seq(min(merged.results.onlySignif.dt[, NULL_ADD_BIC2 - ADD_REC_BIC2]), max(merged.results.onlySignif.dt[, NULL_ADD_BIC2 - ADD_REC_BIC2]))
grid.arrange(
    qplot(
        BIC.breaks, 
        laply(BIC.breaks, function(x) {
            signif <- better$REC.coxtest[merged.results.onlySignif.dt[, NULL_ADD_BIC2 - ADD_REC_BIC2 > x]]
            sum(signif)/length(signif)
        })
    ),
    qplot(
        BIC.breaks, 
        laply(BIC.breaks, function(x) {
            signif <- better$REC.coxtest[merged.results.onlySignif.dt[, NULL_ADD_BIC2 - ADD_REC_BIC2 > x]]
            sum(signif)/length(signif)
        })
    ),
    nrow=1, ncol=2
)

BIC.breaks <- seq(min(merged.results.onlySignif.dt[, NULL_ADD_BIC2 - ADD_GEN_BIC2]), max(merged.results.onlySignif.dt[, NULL_ADD_BIC2 - ADD_GEN_BIC2]))
qplot(
    BIC.breaks, 
    laply(BIC.breaks, function(x) {
        signif <- better$GEN.vuong[merged.results.onlySignif.dt[, NULL_ADD_BIC2 - ADD_GEN_BIC2 > x]]
        sum(signif)/length(signif)
    })
)

BIC.breaks <- seq(min(merged.results.onlySignif.dt[, NULL_ADD_BIC2 - ADD_HET_BIC2]), max(merged.results.onlySignif.dt[, NULL_ADD_BIC2 - ADD_HET_BIC2]))
grid.arrange(
    qplot(
        BIC.breaks, 
        laply(BIC.breaks, function(x) {
            signif <- better$HET.coxtest[merged.results.onlySignif.dt[, NULL_ADD_BIC2 - ADD_HET_BIC2 > x]]
            sum(signif)/length(signif)
        })
    ),
    qplot(
        BIC.breaks, 
        laply(BIC.breaks, function(x) {
            signif <- better$HET.jtest[merged.results.onlySignif.dt[, NULL_ADD_BIC2 - ADD_HET_BIC2 > x]]
            sum(signif)/length(signif)
        })
    ),
    nrow=1, ncol=2
)

# The tests we trust should throw few positives specific to only that test
# From the venn diagrams, generally vuong tests, bic and coxtest.
#
better.DOM <- merged.results.onlySignif.dt[better[, DOM.vuong & DOM.BIC & DOM.coxtest]][, .(rsid, locus.id)]
better.REC <- merged.results.onlySignif.dt[better[, REC.vuong & REC.BIC & REC.coxtest]][, .(rsid, locus.id)]
better.GEN <- merged.results.onlySignif.dt[better[, GEN.vuong & GEN.LRT]][, .(rsid, locus.id)]
better.HET <- merged.results.onlySignif.dt[better[, HET.vuong & HET.BIC & HET.coxtest]][, .(rsid, locus.id)]

# Find snps that seem to have multiple modes of effect
merged.results.onlySignif.dt <- merged.results.onlySignif.dt
merged.results.onlySignif.dt$better.DOM <- merged.results.onlySignif.dt$rsid %in% better.DOM$rsid
merged.results.onlySignif.dt$better.REC <- merged.results.onlySignif.dt$rsid %in% better.REC$rsid
merged.results.onlySignif.dt$better.GEN <- merged.results.onlySignif.dt$rsid %in% better.GEN$rsid
merged.results.onlySignif.dt$better.HET <- merged.results.onlySignif.dt$rsid %in% better.HET$rsid

plotVennWrapper(merged.results.onlySignif.dt, c("better.DOM", "better.REC", "better.GEN", "better.HET"))

#
# Known possible recessive models?
#

# NOD2: containing locus 16:50511169_51009351 , topsnp 16:50745926_C_T
# Gene coords: chr 16, 50693581 to 50733081
better[rsid ==  "16:50745926_C_T", ]
NOD2.rsid <- merged.results.onlySignif.dt[alternate_ids == 16 & position.x >= 50693581 & position.x <= 50733081, rsid]
better[rsid %in% NOD2.rsid, ]

# TYK2: containing locus 19:10408439_10600418 , topsnp 19:10512911_G_A
# Gene coords: chr 19, 10350528 to 10380676
better[rsid ==  "19:10512911_G_A", ]
TYK2.rsid <- merged.results.onlySignif.dt[alternate_ids == 19 & position.x >= 10350528 & position.x <= 10380676, rsid]
better[rsid %in% TYK2.rsid, ]

# Call a "best model" for snps that show evidence for non additive effects
# Based on lowest BIC of significant non additive models.
# Report this lowest BIC

best.BIC.models <- alply(merged.results.onlySignif.dt, 1, function(x) {
    modelNames <- c("DOM", "REC", "GEN", "HET")
    BICs <- with(x, c(ADD_DOM_BIC2, ADD_REC_BIC2, ADD_GEN_BIC2, ADD_HET_BIC2))
    betters <- with(x, c(better.DOM, better.REC, better.GEN, better.HET))
    bestBIC <- BICs[betters]
    if(length(bestBIC)) {
        result <- min(bestBIC)
        names(result) <- modelNames[which(BICs == min(bestBIC))][1]
    } else {
        result <- x$NULL_ADD_BIC2
        names(result) <- "ADD"
    }
    return(result)
}, .progress="text")

merged.results.onlySignif.dt$best.BIC <- unlist(best.BIC.models)
merged.results.onlySignif.dt$best.model <- substr(names(unlist(best.BIC.models)), start=nchar(names(unlist(best.BIC.models))) - 2, stop=nchar(names(unlist(best.BIC.models))))
merged.results.onlySignif.dt$best.BIC.decrease <- merged.results.onlySignif.dt$NULL_ADD_BIC2 - merged.results.onlySignif.dt$best.BIC

table(merged.results.onlySignif.dt$best.model)
hist(merged.results.onlySignif.dt[best.BIC.decrease > 0, best.BIC.decrease], breaks=50)

# Compare rankings of snps by BIC and p value
grid.arrange(
    qplot(log10(merged.results.onlySignif.dt$ADD_pValue.LRT), merged.results.onlySignif.dt$NULL_ADD_BIC2),
    qplot(log10(merged.results.onlySignif.dt$DOM_pValue.LRT), merged.results.onlySignif.dt$ADD_DOM_BIC2),
    qplot(log10(merged.results.onlySignif.dt$REC_pValue.LRT), merged.results.onlySignif.dt$ADD_REC_BIC2),
    qplot(log10(merged.results.onlySignif.dt$GEN_pValue.LRT), merged.results.onlySignif.dt$ADD_GEN_BIC2),
    qplot(log10(merged.results.onlySignif.dt$HET_pValue.LRT), merged.results.onlySignif.dt$ADD_HET_BIC2), 
    nrow=2, ncol=3
)

# Close enough
merged.results.onlySignif.dt[, ADD.BIC.rank := rank(NULL_ADD_BIC2), by=locus.id]
merged.results.onlySignif.dt[, best.BIC.rank := rank(best.BIC), by=locus.id]

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

# Loci that the non additive effects lie in
unique(merged.results.onlySignif.dt[better.DOM | better.REC | better.HET | better.GEN, locus.id])
# [1] "1_161460211_161638410" "6_31011373_32778656"   "8_49047317_49206289"

merged.results.onlySignif.dt[is.topsnp.known == T, ]
merged.results.onlySignif.dt[is.topsnp.observed == T, ]

# Within loci, do the rankings change significantly comparing additive snptest vs additive LRT p values/BIC
# If not, do the rankings change significantly comparing additive BIC vs best BIC

# TODO below here
dev.off()
stopifnot(F)

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

