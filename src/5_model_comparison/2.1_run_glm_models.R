#!/software/R-3.3.0/bin/Rscript
#
library(data.table)
library(plyr)
library(lmtest)
library(nonnest2)

# dataset <- "gwas3"
# assoc <- "ibd"
# chrom.i <- 1
# chunk.i <- 1
# chunk.size <- 100
# gen.file.dir <- "/nfs/users/nfs_b/bb9/workspace/rotation1/crohns_workspace/5_model_comparison/gen/"
# gen.file <- file.path(gen.file.dir, paste(chrom.i, "gen", sep="."))
# gen.file.chunk.file <- paste(gen.file, formatC(chunk.i-1, width=10, flag="0"), sep=".")
# n.snps <- as.numeric(system(paste("wc -l <", gen.file.chunk.file), intern=T))
# out.dir.base <- "/nfs/users/nfs_b/bb9/workspace/rotation1/crohns_workspace/5_model_comparison/"

print("----- Args: ------")
print(commandArgs(T))
print("------------------")
dataset <- commandArgs(T)[1]
assoc <- commandArgs(T)[2]
chrom.i <- commandArgs(T)[3]
chunk.i <- as.numeric(Sys.getenv("LSB_JOBINDEX"))
chunk.size <- as.numeric(commandArgs(T)[4])
gen.file <- commandArgs(T)[5]
n.snps <- as.numeric(commandArgs(T)[6])
out.dir.base <- commandArgs(T)[7]
print(paste("LSB_JOBINDEX (chunk.i + 1) is:", chunk.i))

# Convert gen genotype encoding (3 col)
# > gen.dt[1, 1:10, with=F]
   # V1              V2       V3 V4 V5 V6 V7 V8 V9  V10
# 1: 22 22:21910280_T_C 21910280  T  C  1  0  0  0 0.65
# to expected number of copies
genToDosage <- function(gen.dt, models=c("add", "dom", "rec", "het", "gen")) {
    refHomo.probs <- gen.dt[, seq(6, ncol(gen.dt), 3), with=F]
    het.probs <- gen.dt[, seq(7, ncol(gen.dt), 3), with=F]
    altHomo.probs <- gen.dt[, seq(8, ncol(gen.dt), 3), with=F]
    # Get locations of missing data, denoted by 0 0 0 code
    total.probs <- refHomo.probs + het.probs + altHomo.probs
    # Convert to expected number of copies of the alt 
    dosage <- list()
    if ("gen" %in% models) {
        # General model requires these to be built
        models <- c(models, "add", "rec")
    }
    if ("add" %in% models) {
        dosage$add <- het.probs + 2*altHomo.probs
        is.na(dosage$add) <- total.probs == 0
    }
    if ("dom" %in% models) {
        dosage$dom <- het.probs + altHomo.probs
        is.na(dosage$dom) <- total.probs == 0
    }
    if ("rec" %in% models) {
        dosage$rec <- altHomo.probs
        is.na(dosage$rec) <- total.probs == 0
    }
    if ("het" %in% models) {
        dosage$het <- het.probs
        is.na(dosage$het) <- total.probs == 0
    }
    # Return metadata also.
    dosage$metadata <- gen.dt[, 1:5, with=F]
    names(dosage$metadata) <- c("chr", "snp.id", "position", "allele1", "allele2")
    return(dosage)
}

# Likelihood ratio test for nested GLMs
lrTest <- function(m, null.model) {
    1 - pchisq(null.model$deviance - m$deviance, null.model$df.residual - m$df.residual)
}

getModelSummary <- function(model, null.model, terms=c("geno.add"), prefix="ADD") {
    result <- c()
    for (term in terms) {
        # Get p value, coefficient estimates with std err
        if (term %in% rownames(coef(summary(model)))) {
            result <- c(
                result, 
                estimate.coef=coef(summary(model))[term, "Estimate"],
                stdErr.coef=coef(summary(model))[term, "Std. Error"],
                # z.coef=coef(summary(model))[term, "z value"],
                pValue.coef=coef(summary(model))[term, "Pr(>|z|)"]
            )
        } else {
            result <- c(
                result, 
                estimate.coef=NA,
                stdErr.coef=NA,
                # z.coef=NA,
                pValue.coef=NA
            )
        }
    }
    nested <- all(labels(terms(null.model)) %in% labels(terms(model)))
    result <- c(
        result,
        # Report deviances
        deviance=model$deviance,
        # deviance.null=null.model$deviance,
        # Report LRT against the provided null model, only valid if models are nested
        pValue.LRT=ifelse(nested, lrTest(model, null.model), NA),
        # Report ICs
        # AIC.null=AIC(null.model),
        # AIC=AIC(model),
        # BIC.null=BIC(null.model),
        # BIC=BIC(model),
        # Report if the model converged, or if fitted values are on the boundary
        converged=model$converged,
        boundary=model$boundary
    )
    # Add model name to colnames
    names(result) <- paste(prefix, names(result), sep="_")
    return(data.frame(t(result)))
}

# Perform test of model1 vs model2.
# H0: Model fits are equal
# H1: Alternate model 2 fits better than model 1
compareModels <- function(model1, model2, prefix1="ADD", prefix2="DOM", nested=F) {
    #
    # Handler to suppress warnings due to glm objects having the lm class too.
    #
    handler <- function(w) {
        if (grepl("the condition has length > 1 and only the first element will be used", w)) invokeRestart("muffleWarning")
    }
    # Tests based on Vuong’s (1989) theory of non-nested model comparison 
    # The test is whether the first argument fits better than the second argument.
    # vuongtest.result <- suppressWarnings(vuongtest(model2, model1, nested=nested))
    vuongtest.result <- withCallingHandlers(vuongtest(model2, model1, nested=nested), warning=handler)
    # Variance test to determine if models are distinguishable...
    pValue.vuongtest.varTest <- vuongtest.result$p_omega
    # Perform LRT.
    # ... if models are nested, this is a robust LRT:
    # The traditional LRT uses a chi-square null distribution, whereas the Vuong LRT uses a weighted sum of chi-square
    # distributions. The latter distribution converges to the traditional chi-square distribution when the full model
    # is the true model.
    # ... else this is a non-nested likelihood ratio test to compare model fit.
    pValue.vuongtest.LRT <- vuongtest.result$p_LRT$A
    #
    # Also obtain CI on the difference in AIC/BIC
    # Note: if models are nested or if the "variance test" from
    # ‘vuongtest()’ indicates models are indistinguishable, then the
    # intervals returned from ‘icci()’ will be incorrect.
    #
    icci.result <- icci(model1, model2, conf.level=0.95)
    AIC2 <- icci.result$AIC$AIC2
    BIC2 <- icci.result$BIC$BIC2
    if (nested) {
        AICci.lower <- NA
        AICci.upper <- NA
        BICci.lower <- NA
        BICci.upper <- NA
    } else {
        AICci.lower <- icci.result$AICci[1]
        AICci.upper <- icci.result$AICci[2]
        BICci.lower <- icci.result$BICci[1]
        BICci.upper <- icci.result$BICci[2]
    }
    #
    # Perform coxtest and jtests to attempt to reject model1:
    # Under the assumption that the first argument contains the contains the correct set of regressors,
    # tests whether considering the second argument provides a significant improvement.
    #
    pValue.coxtest <- coxtest(model1, model2)["fitted(M1) ~ M2", "Pr(>|z|)"]
    pValue.jtest <- jtest(model1, model2)["M1 + fitted(M2)", "Pr(>|t|)"]
    #
    # Collate test results
    #
    result <- c(
        AIC2=AIC2,
        AICci.lower=AICci.lower,
        AICci.upper=AICci.upper,
        BIC2=BIC2,
        BICci.lower=BICci.lower,
        BICci.upper=BICci.upper,
        pValue.vuongtest.varTest=pValue.vuongtest.varTest,
        pValue.vuongtest.LRT=pValue.vuongtest.LRT,
        pValue.coxtest=pValue.coxtest,
        pValue.jtest=pValue.jtest
    )
    names(result) <- paste(prefix1, prefix2, names(result), sep="_")
    return(data.frame(t(result)))
}

# Read in phenotype and covariates .sample file
sample.dt <- fread("/lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc_NEW/GWAS3/refs/GWAS3.ibd.sample", header=F, skip=2)
n.samples <- nrow(sample.dt)
colnames(sample.dt) <- colnames(fread("/lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc_NEW/GWAS3/refs/GWAS3.ibd.sample", header=T, nrows=0))

# Standarise covariates
for (n in colnames(sample.dt)) {
    if (grepl("^cov", sample.dt[, n])) {
        sample.dt[[n]] <- scale(sample.dt[[n]])
    }
}

# Create outdir for chunk
print(paste("Processing chr", chrom.i, "in", ceiling(n.snps/chunk.size), "chunks..."))
chunk.out.dir <- file.path(out.dir.base, "chunks", dataset, assoc, chrom.i)
dir.create(chunk.out.dir, recursive=T)

# Read in genotypes .gen file for chunk
gen.file.chunk.file <- paste(gen.file, formatC(chunk.i-1, width=10, flag="0"), sep=".")
print(paste("Reading chunk file:", gen.file.chunk.file))
# chunk.starts <- seq(1, n.snps, chunk.size)
# chunk.start <- chunk.starts[chunk.i]
# gen.dt <- fread(gen.file.chunk.file, sep=" ", skip=chunk.start-1, nrow=chunk.size)
gen.dt <- fread(gen.file.chunk.file, sep=" ")

# Check that we have all genotypes
stopifnot((ncol(gen.dt) - 5)/3 == n.samples)

# Convert to dosages
dosage.dt <- genToDosage(gen.dt)

# Fit a logistic regression with covariates for each site
# Note the the default glm logit is the logistic function
model.null <- glm(bin1 ~ cov1+cov2+cov3+cov4+cov5+cov6+cov7+cov8+cov9+cov10, family=binomial, data=sample.dt)

chunk.results <- ldply(1:nrow(gen.dt), function(snp.i) {

    print(paste(date(), " | Testing snp: ", snp.i, " for: chr_", chrom.i, ".chunk_i_", chunk.i, ".glm_R.out", sep=""))

    geno.add <- unlist(dosage.dt$add[snp.i])
    geno.dom <- unlist(dosage.dt$dom[snp.i])
    geno.rec <- unlist(dosage.dt$rec[snp.i])
    geno.het <- unlist(dosage.dt$het[snp.i])

    model.add <- glm(bin1 ~ cov1+cov2+cov3+cov4+cov5+cov6+cov7+cov8+cov9+cov10 + geno.add, family=binomial, data=sample.dt)
    model.dom <- glm(bin1 ~ cov1+cov2+cov3+cov4+cov5+cov6+cov7+cov8+cov9+cov10 + geno.dom, family=binomial, data=sample.dt)
    model.rec <- glm(bin1 ~ cov1+cov2+cov3+cov4+cov5+cov6+cov7+cov8+cov9+cov10 + geno.rec, family=binomial, data=sample.dt)
    model.het <- glm(bin1 ~ cov1+cov2+cov3+cov4+cov5+cov6+cov7+cov8+cov9+cov10 + geno.het, family=binomial, data=sample.dt)
    model.gen <- glm(bin1 ~ cov1+cov2+cov3+cov4+cov5+cov6+cov7+cov8+cov9+cov10 + geno.add + geno.rec, family=binomial, data=sample.dt)

    # Report coef estimates with stderr for each model
    # Report LRT against appropriate null for each model
    # Report 95% conf. intervals of difference in AIC/BIC to additive model
    # Report model distinguishability and vuong robust LRT 
    cbind(
        NULL_deviance=model.null$deviance,
        getModelSummary(model.add, null.model=model.null, terms=c("geno.add"), prefix="ADD"),
        getModelSummary(model.dom, null.model=model.null, terms=c("geno.dom"), prefix="DOM"),
        getModelSummary(model.rec, null.model=model.null, terms=c("geno.rec"), prefix="REC"),
        getModelSummary(model.het, null.model=model.null, terms=c("geno.het"), prefix="HET"),
        getModelSummary(model.gen, null.model=model.null, terms=c("geno.add", "geno.rec"), prefix="GEN"),
        NULL_AIC1=AIC(model.null),
        NULL_BIC1=BIC(model.null),
        # Test whether model2 fits better
        compareModels(model.null, model.add, prefix1="NULL", prefix2="ADD", nested=T),
        compareModels(model.add, model.dom, prefix1="ADD", prefix2="DOM", nested=F),
        compareModels(model.add, model.rec, prefix1="ADD", prefix2="REC", nested=F),
        compareModels(model.add, model.het, prefix1="ADD", prefix2="HET", nested=F),
        compareModels(model.add, model.gen, prefix1="ADD", prefix2="GEN", nested=T),
        ADD_GEN_pValue.LRT=lrTest(model.gen, null.model=model.add)
    )

}, .parallel=F)

write.table(
    cbind(dosage.dt$metadata, chunk.results), row.names=F, sep="\t",
    file.path(chunk.out.dir, paste("chr_", chrom.i, ".chunk_i_", chunk.i, ".glm_R.out", sep=""))
)

