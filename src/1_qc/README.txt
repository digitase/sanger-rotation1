
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

A comparison of sample-related QC between this pipeline and de Lange (/lustre/scratch114/teams/barrett/coreex_gaibdc/QC/COMBINED/qc_steps.txt):

Number of fails under each criteria from this pipeline:
    
     FID                IID            IMISS.FAIL       HET.FAIL       CHECKSEX.FAIL   ANCESTRY.FAIL   RELATEDNESS.FAIL  ANY.FAIL      
 Length:22252       Length:22252       Mode :logical   Mode :logical   Mode :logical   Mode :logical   Mode :logical    Mode :logical  
 Class :character   Class :character   FALSE:21932     FALSE:21495     FALSE:22001     FALSE:21069     FALSE:20418      FALSE:18755    
 Mode  :character   Mode  :character   TRUE :320       TRUE :757       TRUE :251       TRUE :1183      TRUE :1834       TRUE :3497     
                                       NA's :0         NA's :0         NA's :0         NA's :0         NA's :0          NA's :0 

Number of fails from de Lange indicated with <-:
    ANCESTRY (clear)
    /lustre/scratch113/teams/barrett/coreex_gaibdc/QC/COMBINED/3_sample_QC/coreex_gaibdc_usgwas_qc2.fail-ancestry.txt
    1446 using a PC2 threshold of 0.066
    1466 using a PC2 threshold of 0.067 <-
        We observe 1183, as the 1466 figure includes hapmap samples.

    HET.FAIL (clear)
    /lustre/scratch113/teams/barrett/coreex_gaibdc/QC/COMBINED/3_sample_QC/coreex_gaibdc_usgwas_qc2.fail-het-QC.txt
    758 <-
    but the script /lustre/scratch113/teams/barrett/coreex_gaibdc/QC/scripts/imiss-vs-het.BATCH.Rscript uses 2SDs instead of 3.

    RELATEDNESS.FAIL 
    /lustre/scratch113/teams/barrett/coreex_gaibdc/QC/COMBINED/3_sample_QC/coreex_gaibdc_usgwas_qc2.duplicates_to_remove
    1008 duplicates (kinship > 0.354) (but only 967 unique samples) <-
    /lustre/scratch113/teams/barrett/coreex_gaibdc/QC/COMBINED/3_sample_QC/coreex_gaibdc_usgwas_qc2.first_degree_to_remove
    1049 first degree (kinship > 0.177, < 0.354) (but only 868 unique samples) <-
        1824 unique samples in the union of these two files.

    SEX.FAIL (clear)
    /lustre/scratch113/teams/barrett/coreex_gaibdc/QC/COMBINED/3_sample_QC/fail-gender.txt
    251 <-
        250 + 333835_G08_UC755045 was also removed due to mismatched sex compared to a later release

    # 5) Those with elevated missing rates have already been removed.
    threshold of 0.01 was used

    Total unique samples to remove:
        3185 in /lustre/scratch113/teams/barrett/coreex_gaibdc/QC/COMBINED/3_sample_QC/sample_QC_fail.samples

Explanation of differences
    Samples differing between the pipelines are either in a related pair, or part of a second, extra round of sample missingness filtering that was done in this pipeline.
    e.g. > all(setdiff(failures.summary$FID[failures.summary$ANY.FAIL], katie$V1) %in% union(failures.summary$IID[failures.summary$RELATEDNESS.FAIL], failures.summary$IID[failures.summary$IMISS.FAIL]))
        is TRUE

    Consider discrepencies in related pairs.
        The member of the pair retained by the pipelines differs.
        e.g. 386938_A10_gaibdc6035222 and 386936_E05_gaibdc6034982 are both members of a related pair where neither member fails other filters.
            In these cases, the member with the lower missingness rate should be retained.
                de Lange's pipeline retains the other member.

