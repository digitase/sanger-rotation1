
library(coin)


# alternate_ids rsid chromosome position alleleA alleleB index average_maximum_posterior_call info cohort_1_AA cohort_1_AB cohort_1_BB cohort_1_NULL all_AA all_AB all_BB all_NULL all_total 
# 1 1:161479745_A_G NA 161479745 A G 2203 0.983032 0.968462 3910.33 9537.19 4907.48 0 3910.33 9537.19 4907.48 0 18355 
# cases_AA cases_AB cases_BB cases_NULL cases_total 
# 1797.43 4857.38 2205.18 0 8860 
# controls_AA controls_AB controls_BB controls_NULL controls_total 
# 2112.89 4679.81 2702.29 0 9495 
# all_maf cases_maf controls_maf missing_data_proportion 
# 0.472837 0.476989 0.468963 1.09011e-15 
# het_OR het_OR_lower het_OR_upper 
# 1.22011 1.13239 1.31463 
# hom_OR hom_OR_lower hom_OR_upper 
# 0.959262 0.881645 1.04371 
# all_OR all_OR_lower all_OR_upper 
# 0.968312 0.92941 1.00884 
# frequentist_add_pvalue frequentist_add_info frequentist_add_beta_1 frequentist_add_se_1 
# 0.0706157 0.966752 -0.0396967 0.0219569 
# frequentist_dom_pvalue frequentist_dom_info frequentist_dom_beta_1 frequentist_dom_se_1 
# 0.00181831 0.954272 0.116256 0.0372805 
# frequentist_rec_pvalue frequentist_rec_info frequentist_rec_beta_1 frequentist_rec_se_1 
# 1.14749e-08 0.958788 -0.196243 0.0343842 
# frequentist_gen_pvalue frequentist_gen_info frequentist_gen_beta_1 frequentist_gen_se_1 frequentist_gen_beta_2 frequentist_gen_se_2 
# 7.94825e-14 0.905965 -0.0254919 0.0220373 0.233251 0.0308791 
# frequentist_het_pvalue frequentist_het_info frequentist_het_beta_1 frequentist_het_se_1 
# 1.58606e-14 0.93758 0.236299 0.0307665 
# comment
# NA

# http://stats.stackexchange.com/questions/8774/what-is-the-difference-between-independence-test-in-r-and-cochrane-and-armitage
# As a follow-up to my comment, if independence.test refers to
# coin::independence_test, then you can reproduce a Cochrane and Armitage trend
# test, as it is used in GWAS analysis, as follows:

foo = as.table(rbind(
    round(c(1797.43, 4857.38, 2205.18)),
    round(c(2112.89, 4679.81, 2702.29))
))
rownames(foo) <- c("case", "control")
names(dimnames(foo)) <- c("pheno", "geno")

foo
foo/rowSums(foo)

independence_test(foo, teststat="quad", scores=list(geno=c(0, 1, 2)))
independence_test(foo, teststat="quad", scores=list(geno=c(0, 1, 1)))
independence_test(foo, teststat="quad", scores=list(geno=c(0, 0, 1)))
independence_test(foo, teststat="quad", scores=list(geno=c(0, 1, 0)))

A.freq.case <- 0.5 * (2*foo[1, 1] + foo[1, 2])/sum(foo[1, ])
B.freq.case <- 1 - A.freq.case
case.expected <- round(c(A.freq.case*A.freq.case, 2*A.freq.case*B.freq.case, B.freq.case*B.freq.case)*sum(foo[1, ]))
A.freq.control <- 0.5 * (2*foo[2, 1] + foo[2, 2])/sum(foo[2, ])
B.freq.control <- 1 - A.freq.control
control.expected <- round(c(A.freq.control*A.freq.control, 2*A.freq.control*B.freq.control, B.freq.control*B.freq.control)*sum(foo[2, ]))

# A: 0.489 (492)
# G: 0.511 (514)
# A|A: 0.256 (129)
# A|G: 0.465 (234)
# G|G: 0.278 (140)
eur.1000.genomes.expected <- c(129, 234, 140)

chisq.test(rbind(foo[1, ], case.expected))
chisq.test(rbind(foo[2, ], control.expected))
chisq.test(rbind(foo[2, ], eur.1000.genomes.expected))

