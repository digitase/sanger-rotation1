# Stripe .gen.gz files for faster access.

IN_GEN_FILES_DIR="/lustre/scratch113/projects/crohns/RELEASE/v1/GWAS3/"
DATA_DIR="/lustre/scratch113/projects/crohns/bb9/4_gwas/1_snptest/data/"

mkdir -p "$DATA_DIR"

lfs setstripe "$DATA_DIR" -c -1 

lfs getstripe "$DATA_DIR"

for ((chr = 1; chr <= 22; chr++)) {
    echo "Copying chr $chr gen file..."
    cp "$IN_GEN_FILES_DIR/$chr.gen.gz" "$DATA_DIR/$chr.gen.gz"
}

lfs getstripe "$DATA_DIR"

echo "Data source: $IN_GEN_FILES_DIR" > "$DATA_DIR/README.txt"

