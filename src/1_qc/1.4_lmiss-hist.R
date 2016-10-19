#
# Plot histogram of marker missingness
#
# Modified from /lustre/scratch113/teams/barrett/coreex_gaibdc/QC/scripts/lmiss-hist.Rscript
#

prefix <- commandArgs(T)[1]
thresh <- as.numeric(commandArgs(T)[2])

x <- read.table(paste(prefix,".lmiss",sep=""),header=T)

ylabels=c("0","20K","40K","60K","80K","100K")
xlabels=c("0.0001","0.001","0.01","0.1","1")

pdf(paste(prefix,".lmiss.pdf",sep=""))
    hist(log10(x$F_MISS),axes=F,xlim=c(-4,0),col="RED",ylab="Number of SNPs",xlab="Fraction of missing data",main="All SNPs",ylim=c(0,100000))

    axis(side=2,labels=F)
    mtext(ylabels,side=2,las=2, at=c(0,20000,40000,60000,80000,100000),line=1)

    axis(side=1,labels=F)
    mtext(xlabels,side=1,at=c(-4,-3,-2,-1,0),line=1)

    abline(v=log10(thresh),lty=2)
dev.off()

