#!/usr/local/bin/bash
#
# Evaluate results from QC pipeline reimplementation
# 
# Katie de Lange's pipeline is described in "/lustre/scratch114/teams/barrett/coreex_gaibdc/QC/COMBINED/qc_steps.txt"
#

OUT_DIR="/nfs/users/nfs_b/bb9/workspace/rotation1/crohns_workspace/1_qc/"

# plink file prefix from QC pipeline 
IN_DATA_PREFIX_1="/nfs/users/nfs_b/bb9/workspace/rotation1/crohns_workspace/1_qc/coreex_gaibdc_usgwas_raw.qc6.maf_0.001"
IN_DATA_PREFIX_1="/nfs/users/nfs_b/bb9/workspace/rotation1/crohns_workspace/1_qc/coreex_gaibdc_usgwas_raw.qc6"
IN_DATA_PREFIX_1="/nfs/users/nfs_b/bb9/workspace/rotation1/crohns_workspace/1_qc/coreex_gaibdc_usgwas_raw.qc3"
# plink file prefix from Katie de Lange's QC pipeline 
IN_DATA_PREFIX_2="/lustre/scratch114/teams/barrett/coreex_gaibdc/QC/COMBINED/4_marker_QC/coreex_gaibdc_usgwas_qcIMP"
IN_DATA_PREFIX_2="/lustre/scratch114/teams/barrett/coreex_gaibdc/QC/COMBINED/4_marker_QC/coreex_gaibdc_usgwas_qc4"
IN_DATA_PREFIX_1="/nfs/users/nfs_b/bb9/workspace/rotation1/crohns_workspace/1_qc/coreex_gaibdc_usgwas_raw.qc3"
IN_DATA_PREFIX_2="/lustre/scratch114/teams/barrett/coreex_gaibdc/QC/COMBINED/1_initial_marker_QC/coreex_gaibdc_usgwas_qc1"

# Raw
IN_DATA_PREFIX_1="/nfs/users/nfs_b/bb9/workspace/rotation1/crohns_workspace/1_qc/coreex_gaibdc_usgwas_raw"
python 2.1_compare_plink_files.py "$IN_DATA_PREFIX_1" "$IN_DATA_PREFIX_1" "fam"
python 2.1_compare_plink_files.py "$IN_DATA_PREFIX_1" "$IN_DATA_PREFIX_1" "bim"

# Marker missingness
IN_DATA_PREFIX_1="/nfs/users/nfs_b/bb9/workspace/rotation1/crohns_workspace/1_qc/coreex_gaibdc_usgwas_raw.qc1"
IN_DATA_PREFIX_2="/lustre/scratch114/teams/barrett/coreex_gaibdc/QC/COMBINED/1_initial_marker_QC/coreex_gaibdc_usgwas_qc1"
python 2.1_compare_plink_files.py "$IN_DATA_PREFIX_1" "$IN_DATA_PREFIX_2" "fam"
python 2.1_compare_plink_files.py "$IN_DATA_PREFIX_1" "$IN_DATA_PREFIX_2" "bim"

# Sample missingness
IN_DATA_PREFIX_1="/nfs/users/nfs_b/bb9/workspace/rotation1/crohns_workspace/1_qc/coreex_gaibdc_usgwas_raw.qc3"
IN_DATA_PREFIX_2="/lustre/scratch114/teams/barrett/coreex_gaibdc/QC/COMBINED/2_initial_sample_QC/coreex_gaibdc_usgwas_qc2"
python 2.1_compare_plink_files.py "$IN_DATA_PREFIX_1" "$IN_DATA_PREFIX_2" "fam"
python 2.1_compare_plink_files.py "$IN_DATA_PREFIX_1" "$IN_DATA_PREFIX_2" "bim"

# Other sample QC
IN_DATA_PREFIX_1="/nfs/users/nfs_b/bb9/workspace/rotation1/crohns_workspace/1_qc/coreex_gaibdc_usgwas_raw.qc4"
IN_DATA_PREFIX_2="/lustre/scratch114/teams/barrett/coreex_gaibdc/QC/COMBINED/3_sample_QC/coreex_gaibdc_usgwas_qc3"
python 2.1_compare_plink_files.py "$IN_DATA_PREFIX_1" "$IN_DATA_PREFIX_2" "fam"
python 2.1_compare_plink_files.py "$IN_DATA_PREFIX_1" "$IN_DATA_PREFIX_2" "bim"

# Other marker QC
IN_DATA_PREFIX_1="/nfs/users/nfs_b/bb9/workspace/rotation1/crohns_workspace/1_qc/coreex_gaibdc_usgwas_raw.qc6"
IN_DATA_PREFIX_2="/lustre/scratch114/teams/barrett/coreex_gaibdc/QC/COMBINED/4_marker_QC/coreex_gaibdc_usgwas_qc4"
python 2.1_compare_plink_files.py "$IN_DATA_PREFIX_1" "$IN_DATA_PREFIX_2" "fam"
python 2.1_compare_plink_files.py "$IN_DATA_PREFIX_1" "$IN_DATA_PREFIX_2" "bim"

# Prune MAF
IN_DATA_PREFIX_1="/nfs/users/nfs_b/bb9/workspace/rotation1/crohns_workspace/1_qc/coreex_gaibdc_usgwas_raw.qc6.maf_0.001"
IN_DATA_PREFIX_2="/lustre/scratch114/teams/barrett/coreex_gaibdc/QC/COMBINED/4_marker_QC/coreex_gaibdc_usgwas_qcIMP"
python 2.1_compare_plink_files.py "$IN_DATA_PREFIX_1" "$IN_DATA_PREFIX_2" "fam"
python 2.1_compare_plink_files.py "$IN_DATA_PREFIX_1" "$IN_DATA_PREFIX_2" "bim"

# Alternate version. Easier to output ids.
# function compare_fam_files {
    # echo Present in both
    # comm -12 <(cut -d' ' -f1 "$1.fam" | sort -u) <( cut -d' ' -f1 "$2.fam" | sort -u) | wc -l
    # echo Present in only file 1
    # comm -23 <(cut -d' ' -f1 "$1.fam" | sort -u) <( cut -d' ' -f1 "$2.fam" | sort -u) | wc -l
    # echo Present in only file 2
    # comm -13 <(cut -d' ' -f1 "$1.fam" | sort -u) <( cut -d' ' -f1 "$2.fam" | sort -u) | wc -l
# }

