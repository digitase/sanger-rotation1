#!/software/R-3.3.0/bin/Rscript
#
# Output the positions of snps in each chr
#

library(data.table)

out.dir <- "/nfs/users/nfs_b/bb9/workspace/rotation1/crohns_workspace/5_model_comparison/"
contained.snps.file <- "/nfs/users/nfs_b/bb9/workspace/rotation1/crohns_workspace/5_model_comparison/contained.snps.csv"

# Read in list of snps to test
# These snps that fall within known loci
snps.dt <- fread(contained.snps.file, drop=1)

dir.create(file.path(out.dir, "snps_pos"))
setwd(file.path(out.dir, "snps_pos"))

for (chr in 1:22) {
    write(snps.dt[alternate_ids == chr, position], paste(chr, "txt", sep=".")) 
}

