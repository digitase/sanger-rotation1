#
# Plot first two PCs
# 
# A modification of the script from:
# /lustre/scratch113/teams/barrett/coreex_gaibdc/QC/scripts/plot-pca-results_alpha.Rscript
#

.libPaths("/software/vertres/lib/R/")

# library(ggplot2)
library(plyr)

args <- commandArgs(trailingOnly = TRUE)

data <- read.table(paste(args[1],".evec",sep=""),h=F,skip=1)
colnames(data) <- c("Sample","PC1","PC2","Population")
data$Population <- revalue(data$Population, c("3"="CEU", "4"="CHB", "5"="JPT", "6"="YRI"))

# group.colors <- c(JPT="PURPLE",CHB="PURPLE",YRI="GREEN",CEU="RED",BATCH="ORANGE",Case="BLUE",Control="BLACK")

# ggplot(data,aes(x=PC1,y=PC2,color=Population)) + geom_point(alpha=0.75,size=1) + geom_hline(yintercept=as.numeric(args[2]),linetype="dashed") + geom_vline(xintercept=as.numeric(args[3]),linetype="dashed")  + scale_colour_manual(values=c("PURPLE","YELLOW","GREEN","RED","ORANGE","BLUE","BLACK")) + xlim(0,0.02) + ylim(0.06,0.08)
# ggsave(paste(args[1],".evec.pdf",sep=""), width=8,height=8) 

pdf(paste(args[1],".evec.pdf",sep=""),width=8,height=8)
    plot(data$PC1,data$PC2,col=data$Population)
    abline(h=as.numeric(args[2]),lty=2)
dev.off()

