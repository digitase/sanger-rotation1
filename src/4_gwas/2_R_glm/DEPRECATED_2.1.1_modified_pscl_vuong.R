
source("./nonnest2/R/vuongtest.R")

vuong.modified <- function (m1, m2, digits = getOption("digits")) {
    m1y <- m1$y
    m2y <- m2$y
    m1n <- length(m1y)
    m2n <- length(m2y)
    if (m1n == 0 | m2n == 0) 
        stop("Could not extract dependent variables from models.")
    if (m1n != m2n) 
        stop(paste("Models appear to have different numbers of observations.\n", 
            "Model 1 has ", m1n, " observations.\n", "Model 2 has ", 
            m2n, " observations.\n", sep = ""))
    if (any(m1y != m2y)) {
        stop(paste("Models appear to have different values on dependent variables.\n"))
    }
    p1 <- predprob(m1)
    p2 <- predprob(m2)
    if (!all(colnames(p1) == colnames(p2))) {
        stop("Models appear to have different values on dependent variables.\n")
    }
    whichCol <- match(m1y, colnames(p1))
    whichCol2 <- match(m2y, colnames(p2))
    if (!all(whichCol == whichCol2)) {
        stop("Models appear to have different values on dependent variables.\n")
    }
    m1p <- rep(NA, m1n)
    m2p <- rep(NA, m2n)
    for (i in 1:m1n) {
        m1p[i] <- p1[i, whichCol[i]]
        m2p[i] <- p2[i, whichCol[i]]
    }
    k1 <- length(coef(m1))
    k2 <- length(coef(m2))
    lm1p <- log(m1p)
    lm2p <- log(m2p)
    m <- lm1p - lm2p
    bad1 <- is.na(lm1p) | is.nan(lm1p) | is.infinite(lm1p)
    bad2 <- is.na(lm2p) | is.nan(lm2p) | is.infinite(lm2p)
    bad3 <- is.na(m) | is.nan(m) | is.infinite(m)
    bad <- bad1 | bad2 | bad3
    neff <- sum(!bad)
    if (any(bad)) {
        cat("NA or numerical zeros or ones encountered in fitted probabilities\n")
        cat(paste("dropping these", sum(bad), "cases, but proceed with caution\n"))
    }
    aic.factor <- (k1 - k2)/neff
    bic.factor <- (k1 - k2)/(2 * neff) * log(neff)
    v <- rep(NA, 3)
    arg1 <- matrix(m[!bad], nrow = neff, ncol = 3, byrow = FALSE)
    arg2 <- matrix(c(0, aic.factor, bic.factor), nrow = neff, 
        ncol = 3, byrow = TRUE)
    num <- arg1 - arg2
    s <- apply(num, 2, sd)
    numsum <- apply(num, 2, sum)
    v <- numsum/(s * sqrt(neff))
    names(v) <- c("Raw", "AIC-corrected", "BIC-corrected")
    pval <- rep(NA, 3)
    msg <- rep("", 3)
    for (j in 1:3) {
        # if (v[j] > 0) {
            # pval[j] <- 1 - pnorm(v[j])
            # msg[j] <- "model1 > model2"
        # }
        # else {
            pval[j] <- pnorm(v[j])
            msg[j] <- "model2 > model1"
        # }
    }
    # out <- data.frame(v, msg, format.pval(pval))
    # names(out) <- c("Vuong z-statistic", "H_A", "p-value")
    out <- data.frame(v, msg, pval, format.pval(pval))
    names(out) <- c("Vuong z-statistic", "H_A", "p-value", "p-value.formatted")
    # cat(paste("Vuong Non-Nested Hypothesis Test-Statistic:", 
        # "\n"))
    # cat("(test-statistic is asymptotically distributed N(0,1) under the\n")
    # cat(" null that the models are indistinguishible)\n")
    # cat("-------------------------------------------------------------\n")
    # print(out)
    # return(invisible(NULL))
    return(out)
}

distinguish.test <- function (object1, object2, nested = FALSE, adj = "none") {
    classA <- class(object1)[1L]
    classB <- class(object2)[1L]
    callA <- if (isS4(object1)) 
        object1@call
    else object1$call
    callB <- if (isS4(object2)) 
        object2@call
    else object2$call
    llA <- llcont(object1)
    llB <- llcont(object2)
    if (nested) {
        if (sum(llB) > sum(llA)) {
            tmp <- object1
            object1 <- object2
            object2 <- tmp
            tmp <- llA
            llA <- llB
            llB <- tmp
        }
    }
    n <- NROW(llA)
    omega.hat.2 <- (n - 1)/n * var(llA - llB)
    lamstar <- calcLambda(object1, object2, n)
    pOmega <- imhof(n * omega.hat.2, lamstar^2)$Qq
    return(pOmega)
}

# function (object1, object2, nested = FALSE, adj = "none") 
# {
    # classA <- class(object1)[1L]
    # classB <- class(object2)[1L]
    # callA <- if (isS4(object1)) 
        # object1@call
    # else object1$call
    # callB <- if (isS4(object2)) 
        # object2@call
    # else object2$call
    # llA <- llcont(object1)
    # llB <- llcont(object2)
    # if (nested) {
        # if (sum(llB) > sum(llA)) {
            # tmp <- object1
            # object1 <- object2
            # object2 <- tmp
            # tmp <- llA
            # llA <- llB
            # llB <- tmp
        # }
    # }
    # n <- NROW(llA)
    # omega.hat.2 <- (n - 1)/n * var(llA - llB)
    # lamstar <- calcLambda(object1, object2, n)
    # pOmega <- imhof(n * omega.hat.2, lamstar^2)$Qq
    # lr <- sum(llA - llB)
    # teststat <- (1/sqrt(n)) * lr/sqrt(omega.hat.2)
    # if (adj == "aic") {
        # teststat <- teststat - (length(coef(object1)) - length(coef(object2)))
    # }
    # if (adj == "bic") {
        # teststat <- teststat - (length(coef(object1)) - length(coef(object2))) * 
            # log(n)/2
    # }
    # if (nested) {
        # teststat <- 2 * lr
        # pLRTA <- imhof(teststat, -lamstar)[[1]]
        # pLRTB <- NA
    # }
    # else {
        # pLRTA <- pnorm(teststat, lower.tail = FALSE)
        # pLRTB <- pnorm(teststat)
    # }
    # rval <- list(omega = omega.hat.2, p_omega = pOmega, LRTstat = teststat, 
        # p_LRT = list(A = pLRTA, B = pLRTB), class = list(class1 = classA, 
            # class2 = classB), call = list(call1 = callA, call2 = callB), 
        # nested = nested)
    # class(rval) <- "vuongtest"
    # return(rval)
# }
