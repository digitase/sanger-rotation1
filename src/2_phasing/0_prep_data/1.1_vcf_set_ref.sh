#
# Use Shane's bcftools plugin to switch reference of a VCF
#
# Modified from: /lustre/scratch113/teams/barrett/coreex_gaibdc/imputation/scripts/prep_bcfs.sh
# See plugin code at /nfs/users/nfs_s/sm15/dev/bcftools/misc/set-ref.c
#

# Requires "hgi" in .softwarerc
module add hgi/htslib/latest
module add hgi/bcftools/latest

# export BCFTOOLS_PLUGINS="~sm15/libexec/bcftools_plugins"

export BCFTOOLS_PLUGINS="/nfs/users/nfs_s/sm15/libexec/bcftools_plugins/"

IN_VCF="$1"
IN_REF_FA="$2"
OUT_VCF="$3"

bcftools +set-ref "$IN_VCF" -- -v -f "$IN_REF_FA" 2>"$OUT_VCF.log" >"$OUT_VCF"

