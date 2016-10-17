
library(ggplot2)

in.imiss.file <- commandArgs(T)[1]
in.het.file <- commandArgs(T)[2]
in.imiss.thresh <- as.numeric(commandArgs(T)[3])
in.het.thresh <- as.numeric(commandArgs(T)[4])
out.pdf.file <- commandArgs(T)[5]
out.failure.file <- commandArgs(T)[6]

# in.imiss.file <- "output/coreex_gaibdc_usgwas_raw.qc2.imiss"
# in.het.file <- "output/coreex_gaibdc_usgwas_raw.qc2.het"
# in.imiss.thresh <- 0.01
# in.het.thresh <- 3
# out.pdf.file <- "output/coreex_gaibdc_usgwas_raw.qc2.pdf"
# out.failure.file <- "output/coreex_gaibdc_usgwas_raw.qc2.txt"

imiss.df <- read.table(in.imiss.file, header=T)
het.df <- read.table(in.het.file, header=T)
stopifnot(all(imiss.df$IID == het.df$IID))

# Calculate the observed het rate
het.df$HET.RATE <- (het.df$N.NM. - het.df$O.HOM.) / het.df$N.NM.

het.rate.upperThresh <- mean(het.df$HET.RATE) + in.het.thresh*sd(het.df$HET.RATE)
het.rate.lowerThresh <- mean(het.df$HET.RATE) - in.het.thresh*sd(het.df$HET.RATE)

pdf(out.pdf.file)
    # Inform choice of marker missingness threshold
    hist(log10(imiss.df$F_MISS), main="log10(marker missingness rate)")
    abline(v=log10(in.imiss.thresh), lty=2)
    
    # Also examine choice of excess het rate threshold
    ggplot(mapping=aes(x=imiss.df$F_MISS, y=het.df$HET.RATE)) +
        geom_point() +
        geom_vline(xintercept = in.imiss.thresh) +
        geom_hline(yintercept = het.rate.lowerThresh) +
        geom_hline(yintercept = het.rate.upperThresh) +
        scale_x_log10()
dev.off()

# Write list of failures
merged.imiss.het.df <- merge(imiss.df, het.df)
merged.imiss.het.df <- within(merged.imiss.het.df, {
   IMISS.FAIL <- F_MISS > in.imiss.thresh 
   HET.FAIL <- (HET.RATE < het.rate.lowerThresh) | (HET.RATE > het.rate.upperThresh)
})

print(paste(sum(merged.imiss.het.df$IMISS.FAIL), "samples failed missingness filter."))
print(paste(sum(merged.imiss.het.df$HET.FAIL), "samples failed heterozygosity rate filter."))
print(paste(sum(merged.imiss.het.df$IMISS.FAIL | merged.imiss.het.df$HET.FAIL), "samples failed at least one filter."))

write.table(
    merged.imiss.het.df[merged.imiss.het.df$IMISS.FAIL | merged.imiss.het.df$HET.FAIL, c("FID", "IID", "HET.FAIL", "IMISS.FAIL")], 
    file=out.failure.file, 
    sep="\t", 
    quote=F, 
    row.names=F, col.names=F
)
