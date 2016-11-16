# Check that snptest output for frequentist additive model has been replicated against:
# /lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc_NEW/GWAS3
# NOTE this has essentially been confirmed.

katie.file <- "gwas3.ibd.0-1999999"
my.file <- "chunk_i.1800000_1999999.snptest.out"

options(stringsAsFactors = F)

katie.df <- read.table(katie.file, skip = 58000, nrows = 10000)
my.df <- read.table(my.file, header = T)
colnames(katie.df) <- colnames(my.df)[1:ncol(katie.df)]
colnames(katie.df)[ncol(katie.df)] <- "comment"

katie.df.filtered <- katie.df[katie.df$rsid %in% my.df$rsid, ]
stopifnot(all(my.df$rsid == katie.df.filtered$rsid))

for (n in colnames(katie.df.filtered)) {
    if (n == "index" || n == "comment") {
        next
    }
    if (any(na.omit(my.df[[n]]) != na.omit(katie.df.filtered[[n]]))) {
        print(n)
    }
}
