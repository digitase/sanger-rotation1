#!/software/R-3.3.0/bin/Rscript
#
library(data.table)
library(gtools)
library(plyr)
library(ggplot2)
library(gplots)

# 
# Read in snptest results for loci that were re run in glm
contained.snps.dt <- fread("/nfs/users/nfs_b/bb9/workspace/rotation1/crohns_workspace/4_gwas/2_R_glm/contained.snps.txt")

# Read in results from glm
all.glm.out.files <- mixedsort(list.files("/nfs/users/nfs_b/bb9/workspace/rotation1/crohns_workspace/4_gwas/2_R_glm/chunks/gwas3/ibd/", include.dirs=F, recursive=T, full.names=T))
glm.results.dt <- data.table(ldply(all.glm.out.files, fread))

# Merge in snptest p values
merged.results.dt <- merge(contained.snps.dt, glm.results.dt, by.x=c("rsid"), by.y=c("snp.id"), all.y=T)
# Find topsnps
merged.results.dt[, is.topsnp := topSNP.Position..bp. == position.x, by=locus.id]
# Get model number for lowest p value for each snp
merged.results.dt[, lowest.pValue.modelNumber := which.min(c(frequentist_add_pvalue, frequentist_dom_pvalue, frequentist_rec_pvalue, frequentist_gen_pvalue, frequentist_het_pvalue)), by=1:nrow(merged.results.dt)]

#
# Start comparisons
#

# Can we trust our analysis techniques?

# Coefficient p value is NOT the same p value as the lrt
ggplot(data=merged.results.dt, aes(x=ADD_pValue.LRT, y=ADD_pValue.coef)) + 
    geom_point(aes(colour=frequentist_add_info), size=0.1)

# LRT value corresponds well to snptest values with high INFO
# For ADD
ggplot(data=merged.results.dt, aes(x=frequentist_add_pvalue, y=ADD_pValue.LRT)) +
    geom_point(aes(colour=frequentist_add_info), size=0.1) +
    scale_colour_gradient(low="red") +
    # scale_x_log10() + scale_y_log10() +
    geom_text(data=merged.results.dt[abs(1 - frequentist_add_pvalue/ADD_pValue.LRT) > 0.95], aes(label=round(frequentist_add_info, digits=2)), size=3) 

# For GEN
summary(merged.results.dt$frequentist_gen_pvalue)
summary(merged.results.dt$GEN_pValue.LRT)
ggplot(data=merged.results.dt, aes(x=frequentist_gen_pvalue, y=GEN_pValue.LRT)) +
    geom_point(aes(colour=frequentist_gen_info), size=0.1) +
    scale_colour_gradient(low="red") +
    # scale_x_log10() + scale_y_log10() +
    geom_text(data=merged.results.dt[abs(1 - frequentist_gen_pvalue/GEN_pValue.LRT) > 0.95], aes(label=round(frequentist_gen_info, digits=2)), size=3) 

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

# Alternate filtering based on LRT vs null glm.
# merged.results.onlySignif.dt <- merged.results.dt[
    # ADD_pValue.LRT < gwas.thresh |
    # DOM_pValue.LRT < gwas.thresh |
    # REC_pValue.LRT < gwas.thresh |
    # GEN_pValue.LRT < gwas.thresh |
    # HET_pValue.LRT < gwas.thresh 
# ]
# dim(merged.results.onlySignif.dt)

# Get stats on what proportion of lowest p values are from each type of model
hist(merged.results.onlySignif.dt$lowest.pValue.modelNumber)

# Check the power of the cohort analysis
# How many known IBD loci remain?
length(unique(contained.snps.dt$locus.id))
length(unique(merged.results.onlySignif.dt$locus.id))
# How many of the corresponding topsnps remain?
nrow(contained.snps.dt[topSNP.Position..bp. == position])
sum(merged.results.onlySignif.dt$is.topsnp)

#
# Detecting snps where non additive models are better
#

# Apply FDR corrections to relevant cols
merged.results.onlySignif.dt[, ADD_DOM_pValue.coxtest.FDR := p.adjust(ADD_DOM_pValue.coxtest, method="holm")]
merged.results.onlySignif.dt[, ADD_REC_pValue.coxtest.FDR := p.adjust(ADD_REC_pValue.coxtest, method="holm")]
merged.results.onlySignif.dt[, ADD_HET_pValue.coxtest.FDR := p.adjust(ADD_HET_pValue.coxtest, method="holm")]

merged.results.onlySignif.dt[, ADD_DOM_pValue.jtest.FDR := p.adjust(ADD_DOM_pValue.jtest, method="holm")]
merged.results.onlySignif.dt[, ADD_REC_pValue.jtest.FDR := p.adjust(ADD_REC_pValue.jtest, method="holm")]
merged.results.onlySignif.dt[, ADD_HET_pValue.jtest.FDR := p.adjust(ADD_HET_pValue.jtest, method="holm")]

merged.results.onlySignif.dt[, ADD_DOM_pValue.vuongtest.varTest.FDR := p.adjust(ADD_DOM_pValue.vuongtest.varTest, method="holm")]
merged.results.onlySignif.dt[, ADD_REC_pValue.vuongtest.varTest.FDR := p.adjust(ADD_REC_pValue.vuongtest.varTest, method="holm")]
merged.results.onlySignif.dt[, ADD_GEN_pValue.vuongtest.varTest.FDR := p.adjust(ADD_GEN_pValue.vuongtest.varTest, method="holm")]
merged.results.onlySignif.dt[, ADD_HET_pValue.vuongtest.varTest.FDR := p.adjust(ADD_HET_pValue.vuongtest.varTest, method="holm")]

merged.results.onlySignif.dt[, ADD_DOM_pValue.vuongtest.LRT.FDR := p.adjust(ADD_DOM_pValue.vuongtest.LRT, method="holm")]
merged.results.onlySignif.dt[, ADD_REC_pValue.vuongtest.LRT.FDR := p.adjust(ADD_REC_pValue.vuongtest.LRT, method="holm")]
merged.results.onlySignif.dt[, ADD_GEN_pValue.vuongtest.LRT.FDR := p.adjust(ADD_GEN_pValue.vuongtest.LRT, method="holm")]
merged.results.onlySignif.dt[, ADD_HET_pValue.vuongtest.LRT.FDR := p.adjust(ADD_HET_pValue.vuongtest.LRT, method="holm")]

merged.results.onlySignif.dt[, ADD_GEN_pValue.LRT.FDR := p.adjust(ADD_GEN_pValue.LRT, method="holm")]

# Ensure rsids are unique in this subset
stopifnot(length(merged.results.onlySignif.dt$rsid) == length(unique(merged.results.onlySignif.dt$rsid)))
better <- data.table(rsid=merged.results.onlySignif.dt$rsid)

#
# A decision based on raw p value
#
# Find snps with better non additive p values
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
intersect(
    merged.results.onlySignif.dt[ADD_DOM_pValue.vuongtest.varTest < 0, rsid],
    merged.results.onlySignif.dt[ADD_DOM_BICci.lower * ADD_DOM_BICci.upper > 0 & ADD_DOM_BICci.lower > 0, rsid]
)
#
# This is not the case, thus better to just exclude negative p values.

# Use vuong 2 step test (non nested version)
better[
    , 
    DOM.vuong := merged.results.onlySignif.dt[
        , 
        ADD_DOM_pValue.vuongtest.varTest.FDR > 0 & ADD_DOM_pValue.vuongtest.varTest.FDR < signif.thresh &
        ADD_DOM_pValue.vuongtest.LRT > 0 &
        ADD_DOM_pValue.vuongtest.LRT < signif.thresh
    ]
]
better[
    , 
    REC.vuong := merged.results.onlySignif.dt[
        , 
        ADD_REC_pValue.vuongtest.varTest.FDR > 0 &
        ADD_REC_pValue.vuongtest.varTest.FDR < signif.thresh &
        ADD_REC_pValue.vuongtest.LRT > 0 &
        ADD_REC_pValue.vuongtest.LRT < signif.thresh
    ]
]
better[
    , 
    HET.vuong := merged.results.onlySignif.dt[
        , 
        ADD_HET_pValue.vuongtest.varTest.FDR > 0 &
        ADD_HET_pValue.vuongtest.varTest.FDR < signif.thresh &
        ADD_HET_pValue.vuongtest.LRT > 0 &
        ADD_HET_pValue.vuongtest.LRT < signif.thresh
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
better[, DOM.coxtest := merged.results.onlySignif.dt[, ADD_DOM_pValue.coxtest.FDR > 0 & ADD_DOM_pValue.coxtest.FDR < signif.thresh]]
better[, REC.coxtest := merged.results.onlySignif.dt[, ADD_REC_pValue.coxtest.FDR > 0 & ADD_REC_pValue.coxtest.FDR < signif.thresh]]
better[, HET.coxtest := merged.results.onlySignif.dt[, ADD_HET_pValue.coxtest.FDR > 0 & ADD_HET_pValue.coxtest.FDR < signif.thresh]]

#
# Decision methods for ADD vs GEN (nested)
#
# Use the traditional LRT
better[, GEN.LRT := merged.results.onlySignif.dt[, ADD_GEN_pValue.LRT.FDR < signif.thresh]]

# Use vuong 2 step test (nested version)
#
better[
    , 
    GEN.vuong := merged.results.onlySignif.dt[
        , 
        ADD_GEN_pValue.vuongtest.varTest.FDR > 0 &
        ADD_GEN_pValue.vuongtest.varTest.FDR < signif.thresh &
        ADD_GEN_pValue.vuongtest.LRT > 0 &
        ADD_GEN_pValue.vuongtest.LRT < signif.thresh
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

plotVennWrapper(better, c("DOM.raw.pValue", "REC.raw.pValue", "GEN.raw.pValue", "HET.raw.pValue"))

plotVennWrapper(better, c("DOM.BIC", "REC.BIC", "HET.BIC"))

plotVennWrapper(better, c("DOM.vuong", "REC.vuong", "HET.vuong"))

plotVennWrapper(better, c("DOM.coxtest", "REC.coxtest", "HET.coxtest"))


Same model, diff methods

foo=as.matrix(better[, sort(names(better)), with=F][, -ncol(better), with=F])
foo[foo] <- 1

image(foo, col=c("red", "blue"), xlab="genome-wide signif in at least 1 model snps")



bar <- melt(better[, sort(names(better)), with=F], id.vars=c("rsid"))

ggplot(bar, aes(rsid, reorder(variable, variable))) + 
    geom_tile(aes(fill = value)) + theme(axis.text.x = element_blank())


# Known possible recessive models?
#
# NOD2: 16:50511169_51009351 , 16:50745926_C_T
merged.results.onlySignif.dt[rsid ==  "16:50745926_C_T", ]
better[rsid ==  "16:50745926_C_T", ]
# TYK2: 19:10408439_10600418 , 19:10512911_G_A
merged.results.onlySignif.dt[rsid ==  "19:10512911_G_A", ]
better[rsid ==  "19:10512911_G_A", ]

# Is the end result enriched for snps with high INFO?

