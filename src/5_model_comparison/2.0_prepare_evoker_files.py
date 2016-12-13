#!/usr/bin/env python
#
# Source: /nfs/users/nfs_d/dr9/prepare_evoker_files.py
# Author: Dan Rice
#

from plumbum import local
import os
import re
import shutil
import pandas as pd
from glob import glob
import tarfile


plink = local['/nfs/users/nfs_d/dr9/software/plink/plink']
int2bnt = local['/nfs/users/nfs_d/dr9/evoker/int2bnt.pl']


def extract_and_write_intensity(original_intensity_file, new_intensity_file, rsids):
    regex_or = '|'.join(rsids)
    rsids = set(rsids)
    p = re.compile(r'({})\s+'.format(regex_or))
    rsids_seen_so_far = set()
    with open(original_intensity_file) as i, open(new_intensity_file, 'w') as o:
        header = i.readline()
        o.write(header)
        for line in i:
            if rsids_seen_so_far == rsids:
                return
            m = p.match(line)
            if m:
                v = m.groups()[0]
                rsids_seen_so_far.add(v)
                o.write(line)

                
def sync_intensity_fam_and_write(original_intensity_file, new_intensity_file, fam_file):
    """
    Sometimes the fam rows and the intensity file columns are not in agreement.
    """
    samples = []
    with open(fam_file) as i:
        for line in i:
            tokens = line.split()
            assert tokens[0] == tokens[1]
            samples.extend([tokens[0] + 'A', tokens[0] + 'B'])
    df = pd.read_csv(original_intensity_file, sep='\t', dtype='str')
    fam_order = list(df.columns[:3]) + samples
    if fam_order == list(df.columns):
        print('\nOrder is the same, skipping: \n{}\n{}\n'.format(intensity, fam))
        return
    df = df[fam_order]
    df.to_csv(new_intensity_file, sep='\t', index=False)
    
    
def make_tarfile(output_filename, source_dir):
    with tarfile.open(output_filename, "w:gz") as tar:
        tar.add(source_dir, arcname=os.path.basename(source_dir))


def prepare_evoker_files(work_path, snps):
    
    #-----------------------------------------------------------------------------------------
    # These are correct for GWAS3
    CASE = 'coreex_gaibdc_raw'
    CONTROL = 'coreex_usgwas_raw'
    basenames = {
        CASE: 'coreex_gaibdc',
        CONTROL: 'coreex_usgwas'
    }
    intensities_path = {
        CASE: '/lustre/scratch115/teams/barrett/users/dr9/coreex_gaibdc_20150304/coreex_gaibdc_20150304_{}.int',
        CONTROL: '/lustre/scratch115/teams/barrett/users/dr9/coreex_usgwas_20131119.int/coreex_usgwas_20131119_{}.int'
    }
    plink_bfile_path = '/lustre/scratch113/teams/barrett/coreex_gaibdc/evoker/data/{chrom}/{collection}.{chrom}{suffix}'
    qc_fail_samples_path = '/lustre/scratch113/teams/barrett/coreex_gaibdc/evoker/qc_fail_samples.txt'
    #-----------------------------------------------------------------------------------------

    snps_path = os.path.join(work_path, 'snp_keep')

    if not os.path.exists(work_path):
        os.mkdir(work_path)

    snp_chroms = set([el['chrom'] for el in snps])
    for chrom in snp_chroms:
        rsids_with_chrom = [el['rsid'] for el in snps if el['chrom'] == chrom]
        rsids_with_chrom = '\n'.join(rsids_with_chrom)
        with open(snps_path, 'w') as o:
            o.write(rsids_with_chrom)
        for collection in [CASE, CONTROL]:
            out = os.path.join(work_path, basenames[collection] + '.' + chrom)
            bfile = plink_bfile_path.format(chrom=chrom, collection=collection, suffix='')
            print(plink([
                '--bfile', bfile,
                '--extract', snps_path,
                '--make-bed',
                '--out', out,
            ]))

            print(plink([
                '--bfile', out,
                '--remove', qc_fail_samples_path,
                '--make-bed',
                '--out', out,
            ]))

            extracted_intensity_path = out + '.int'
            variant_rsids = [el['rsid'] for el in snps if el['chrom'] == chrom]
            intensity_path = intensities_path[collection].format(chrom)
            extract_and_write_intensity(intensity_path, extracted_intensity_path, variant_rsids)
            synced_intensity_path = out + '.synced.int'
            sync_intensity_fam_and_write(original_intensity_file=extracted_intensity_path, \
                                         new_intensity_file=synced_intensity_path, \
                                         fam_file=out + '.fam' \
                                        )

            # Convert int to bnt
            print(int2bnt([
                '-f', 'illuminus',
                '-i', synced_intensity_path,
                '-o', out + '.bnt'
            ]))

            # Can only have one fam
            src = out + '.fam'
            dst = out[:-2] + '.fam'
            shutil.move(src, dst)

    to_remove = ('~', '.nosex', '.log', '.synced.int', '.int', 'snp_keep')
    for suffix in to_remove:
        for f in glob(work_path + '*' + suffix):
            os.remove(f)

    # Create marker list.
    rsids = '\n'.join([el['rsid'] for el in snps])
    with open(os.path.join(work_path, 'markers.txt'), 'w') as o:
        o.write(rsids)

    evoker_directory = os.path.join(work_path, 'evoker')
    if not os.path.exists(evoker_directory):
        os.makedirs(evoker_directory)
    for f in glob(os.path.join(work_path, '*')):
        shutil.move(f, evoker_directory)

    make_tarfile(os.path.join(work_path, 'evoker.tar.gz'), evoker_directory)
    
    
def main():
    
    # This directory should be empty because temp files are deleted.
    #  work_path = '/lustre/scratch115/teams/barrett/users/dr9/gwas3_ben/'
    work_path = '/lustre/scratch113/projects/crohns/bb9/5_model_comparison/2.0_prepare_evoker_files/'

    # The snps you want to view.
    snps = [
        {
            'chrom': '1',
            'coor': '161479745',
            'rsid': 'exm117933'
        },
        {
            'chrom': '1',
            'coor': '161479745',
            'rsid': 'rs1801274'
        },
        {
            'chrom': '1',
            'coor': '169511555',
            'rsid': 'exm121879'
        },
        {
            'chrom': '8',
            'coor': '49053165',
            'rsid': 'exm-rs11778329'
        },
    ]
    
    prepare_evoker_files(work_path, snps)
    
    
if __name__ == '__main__':
    main()
