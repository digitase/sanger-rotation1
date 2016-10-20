
QC pipeline for GWAS chip data.
Currently applied to GWAS3 data.

A merger of several different pipelines/recommendations:

1.
    213 GWAS quality control. We removed variants that did not overlap between the two versions of the chip, 
    214 had missingness > 5%, a significant difference in call rate between cases and controls (P < 1x10-?5), 
    215 deviated from Hardy-?Weinberg equilibrium (HWE) in controls (P < 1x10-?5), or that were affected by a 
    216 genotyping batch effect (significant association [P < 1x10-?5] between an outlier group of cases discovered 
    217 using principal component analysis [PC1 < -?0.005], and the remainder of the samples). We then removed 
    218 samples with missingness > 1%, heterozygosity ±3 standard deviations from the mean, mismatch 
    219 between reported and genotypic sex, first-?degree relatives or closer (kinship coefficient > 0.177), and 
    220 non-?European samples identified through principal component analysis with HapMap3 populations. After 
    221 quality control, data were available for 4,474 Crohn’s disease, 4,173 ulcerative colitis, 592 IBD-? 
    222 unclassified cases and 9,500 controls for 296,203 variants.

    Ref:
        de Lange, K. M. K. M. et al., Moutsianas, L., Lee, J. C., Lamb, C. A.,
        Luo, Y., Kennedy, N. A., … Barrett, J. (2016). Genome-wide association study
        implicates immune activation of multiple integrin genes in inflammatory bowel
        disease. BIORXIV, 1–19. https://doi.org/10.1101/058255

2.
   Description from Katie de Lange: 

   1. Remove markers with extreme missingness (will mess up individual statistics if kept)
   2. Remove individuals with extreme missingness (will mess up downstream QC metrics if kept)
   3. Check all other individual filters. For related pairs, only remove those with lowest call rate or mismatched gender as applicable. Remove everyone who fails one or more filters.
   4. Check all other marker filters.

   Order details: /lustre/scratch114/teams/barrett/coreex_gaibdc/QC/COMBINED/qc_steps.txt

3.
    Recommendations from:
    
    Anderson, C. A., Pettersson, F. H., Clarke, G. M., Cardon, L. R., Morris, A.
    P., & Zondervan, K. T. (2010). Data quality control in genetic case-control
    association studies. Nature Protocols, 5(9), 1564–73.
    https://doi.org/10.1038/nprot.2010.116

Notes:
    On the known plate swap 

        There was a plate swap, between samples on plates
        333833 and 333835, so you can check for those numbers in the sample names
        (PLATE_WELL_SAMPLE). If you want to be super careful though, they fixed the
        genders in the 4th release of the data, so you can check back against the
        fam file in here:
        /lustre/scratch114/teams/barrett/coreex_gaibdc/release/coreex_gaibdc_20150304

        de Lange, 2016-10-13

    On the order of filters

        I applied most of the filters all at once - mostly because some of them
        were inter-dependent. For example, when you are removing one individual
        from a related pair, you want to make sure that if one individual fails (or
        is worse, even if it doesn’t reach the failure threshold) on other
        filter(s), then that is the one you remove. Otherwise you may end up
        throwing away both individuals, even though one is fine. (edited) For most
        of the other filters, it doesn’t matter if you apply them together or
        separately, because failure is failure, regardless of their status
        elsewhere (e.g. if it is a missingness fail, then you are going to remove
        it no matter what) The only one I would really suggest making sure you do
        separately first is extreme missingness, as these can affect your
        calculations of the mean and sd for heterozygosity rate.

        de Lange, 2016-10-14

Cleaned dataset stats after filtering for MAF from the final plink log file.

    296117 variants and 18846 people pass filters and QC.  
    Among remaining phenotypes, 9324 are cases and 9522 are controls.
    
    Comparison with existing results in de Lange 2016 BIORXIV:
        9239 cases, 9500 controls over 296,203 variants (filtering for MAF > 0.01% assumed).

