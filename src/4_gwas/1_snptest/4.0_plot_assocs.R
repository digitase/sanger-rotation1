#!/software/R-3.3.0/bin/Rscript
#
# Read in snptest association tables.
# Plot QQ plots and Manhattan plots for each chromosome, for each model.
#
# Based on /lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc_NEW/GWAS3/scripts/plot_manhattan.R
#
library(data.table)
library(plyr)
library(doMC)

source("/lustre/scratch113/projects/crohns/2013Aug07/assoc/qqman.r")

# bsub wrapper will be required to run chroms in parallel
# Approx. 8GB required per core to read in snptest data
# Approx. 6GB required per core to plot models in parallel
# bsub -n 5 -R "span[hosts=1] select[mem>10000] rusage[mem=10000]" -M 10000 -o out.log -e err.log "./4.0_plot_assocs.R"
doMC::registerDoMC(cores=5)

snptest.out.dir <- "/nfs/users/nfs_b/bb9/workspace/rotation1/src/4_gwas/1_snptest/.output/results.filtered/gwas3/ibd"
# snptest.out.dir <- commandArgs(T)[1]

date()

plots.out.dir <- file.path(snptest.out.dir, "plots")
dir.create(plots.out.dir)
setwd(plots.out.dir)

l_ply(list.files(snptest.out.dir, full.names=T), function(snptest.result.file) {
    
    # Get chr number from filename
    snptest.result.file.basename <- basename(snptest.result.file)
    chr <- as.numeric(sub("chr_", "", strsplit(snptest.result.file.basename, split=".", fixed=T)[[1]][1]))

    # Read in snptest assoc values
    print(paste("Reading in snptest assocs for chr", chr, ":", snptest.result.file))
    data <- fread(snptest.result.file, header=T, showProgress=F)

    # Plot QQ plot and manhattan plot for each type of model
    pvalue_colnames <- c("frequentist_add_pvalue", "frequentist_dom_pvalue", "frequentist_rec_pvalue", "frequentist_gen_pvalue", "frequentist_het_pvalue")

    l_ply(pvalue_colnames, function(pvalue_colname) {

        print(paste("Plotting QQ plot for model:", pvalue_colname, ", chr:", chr))
        png(paste(snptest.result.file.basename, pvalue_colname, "qq.png", sep="."), h=1000, w=1000)
            lambda <- median(qchisq(na.omit(data[[pvalue_colname]]), df=1, lower.tail=F), na.rm=T)/0.456 
            qq(na.omit(data[[pvalue_colname]]), main=lambda)
        dev.off()

        print(paste("Plotting Manhattan plot for model:", pvalue_colname, ", chr:", chr))
        png(paste(snptest.result.file.basename, pvalue_colname, "manhattan.png", sep="."), h=1000, w=1000)
            subset.nonMissing <- na.omit(data.frame(cbind(data$position, data[[pvalue_colname]])))
            colnames(subset.nonMissing) <- c("BP", "P")
            subset.nonMissing$CHR <- chr
            # Color scheme values from http://colorbrewer2.org/
            manhattan(subset.nonMissing, col=c("#08306b", "#4292c6"), 
                      suggestiveline=-log10(1e-5), genomewideline=-log10(5e-8),
                      pch=20, cex=2)
        dev.off()

    }, .parallel=T)

    rm(data)
    gc()

}, .parallel=F)

date()

