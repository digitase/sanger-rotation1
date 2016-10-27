
Phase genotypes -> haplotypes.

Operate on a small chromosome e.g. chr20 to save time.

The bigger plan for phasing and imputation:
    Phasing with EAGLE2
    Phasing with pbwt

    Try imputation with the Sanger imputation service: https://imputation.sanger.ac.uk/
        which implements this pipeline, imputing against the 
            Haplotype Reference Consortium (release 1.1)    32,470  39M     autosomes only; SNPs only
                S. McCarthy et al. (2016) A reference panel of 64,976 haplotypes for genotype imputation, Nature Genetics. 48(10):1279-83 [27548312]

    Also try imputation with pbwt against the 1000 genomes reference panel

        The 1000 genomes reference panel (2504? individuals) is located at:
        /lustre/scratch113/projects/crohns/2013Aug07/imputation/reference/ALL.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nosing

    and then comparison of results between panels.

