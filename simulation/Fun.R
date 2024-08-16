### INT
int = function(x){return(qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x))))}

# Simulate distribution
grt_err = function(n, dis){
  ### dis=1: normal distribution
  ### dis=2-4: t distirbution with df=10,5,3
  ### dis=5-7; chisq distribution with df=15,5,1
  if(dis==1){e = rnorm(n);escaled=(e-0)/1
  }else if(dis==2){e = rt(n,10);escaled=(e-0)/sqrt(10/8)
  }else if(dis==3){e = rt(n,5);escaled=(e-0)/sqrt(5/3)
  }else if(dis==4){e = rt(n,3);escaled=(e-0)/sqrt(3/1)
  }else if(dis==5){e = rchisq(n,15);escaled=(e-15)/sqrt(2*15)
  }else if(dis==6){e = rchisq(n,5);escaled=(e-5)/sqrt(2*5)
  }else if(dis==7){e = rchisq(n,1);escaled=(e-1)/sqrt(2*1)
  }else if(dis==8){e =  rexp(n, rate=1);escaled=(e-1)/sqrt(2*1)
  }else{stop("ERROR: dis option can only be 1-7.")}
  return(escaled)
}

require(data.table)
require(plyr)
require(dplyr)
require(optparse)

read_sumstats <- function(ss, snp_col="SNP", a1_col="A1", a2_col="A2", z_col="Z") {

    # ss could be either path or data.frame
    if (is.character(ss)) {
        if (file.exists(ss)) {
            ss <- fread(ss)
        } else {
            cat("ERROR: \n")
            q()
        }
    }

    ss <- as.data.frame(ss)
    header <- c(snp_col, a1_col, a2_col, z_col)
    if (sum(header %in% colnames(ss)) != length(header)) {cat("ERROR: line91\n"); q()}
    ss <- ss[, header]
    colnames(ss) <- c("SNP", "A1", "A2", "Z")
    ss <- ss[!is.na(ss$Z),]

    return(as.data.frame(ss))
}

read_sumstats_noallele <- function(ss, snp_col="SNP", z_col="Z") {

    # ss could be either path or data.frame
    if (is.character(ss)) {
        if (file.exists(ss)) {
            ss <- fread(ss)
        } else {
            cat("ERROR: \n")
            q()
        }
    }

    ss <- as.data.frame(ss)
    header <- c(snp_col, z_col)
    if (sum(header %in% colnames(ss)) != length(header)) {cat("ERROR: line91\n"); q()}
    ss <- ss[, header]
    colnames(ss) <- c("SNP", "Z")
    ss <- ss[!is.na(ss$Z),]

    return(as.data.frame(ss))
}


## no direct use, if you need, just input sumstats data frame with colnames including SNP and Z at least
merge_match_a1a2 <- function(ss1, ss2) {

    ss1 <- as.data.frame(ss1)
    ss2 <- as.data.frame(ss2)

    snp_ovp <- intersect(ss1$SNP, ss2$SNP)
    ss1 <- ss1[match(snp_ovp, ss1$SNP),]
    ss2 <- ss2[match(snp_ovp, ss2$SNP),]

    not.match.not.flip <- ((ss1$A1 %in% c("A","T")) & (ss1$A2 %in% c("A","T")) & (ss2$A1 %in% c("G", "C") | ss2$A2 %in% c("G", "C"))) | ((ss1$A1 %in% c("G","C")) & (ss1$A2 %in% c("G","C")) & (ss2$A1 %in% c("A", "T") | ss2$A2 %in% c("A", "T"))) | ((ss1$A1 %in% c("A","T") | ss1$A2 %in% c("A","T")) & (ss2$A1 %in% c("G", "C")) & (ss2$A2 %in% c("G", "C"))) | ((ss1$A1 %in% c("G","C") | ss1$A2 %in% c("G","C")) & (ss2$A1 %in% c("A", "T")) & (ss2$A2 %in% c("A", "T")))
    cat("  Num of SNPs whose A1, A2 neither matched nor flipped:", paste0(sum(not.match.not.flip), ". "))

    if (sum(not.match.not.flip) > 0) {
        cat("Remove these SNPs.")
        ss1 <- ss1[!not.match.not.flip,]
        ss2 <- ss2[!not.match.not.flip,]
    }
    cat("\n")

    con <- ss1$A1==ss2$A1 | ((ss1$A1 %in% c("A","T")) & (ss2$A1 %in% c("A","T"))) | ((ss1$A1 %in% c("G","C")) & (ss2$A1 %in% c("G","C")))
    ss2$Z <- ifelse(con, ss2$Z, -ss2$Z)
    ss <- cbind(ss1, ss2$Z)

    return(as.data.frame(ss))
}

merge_nomatch_a1a2 <- function(ss1, ss2) {

    ss1 <- as.data.frame(ss1)
    ss2 <- as.data.frame(ss2)

    snp_ovp <- intersect(ss1$SNP, ss2$SNP)
    ss1 <- ss1[match(snp_ovp, ss1$SNP),]
    ss2 <- ss2[match(snp_ovp, ss2$SNP),]

    ss <- cbind(ss1, ss2$Z)

    return(as.data.frame(ss))
}


## sumstats path/data.frame -> Z table
# input: a vector/list of sumstats path/data.frame in order of (1) predicted-outcome on labeled dataset (2) predicted-outcome on unlabeled dataset (3) labeled outcome in labeled dataset
# optional input: vector/list of column names if different from our usual format
read_z <- function(sumstats, snp_cols=rep("SNP",3), z_cols=rep("Z",3), no.allele=F, a1_cols=rep("A1", 3), a2_cols=rep("A2", 3)) {

    if (no.allele) {
        sumstats <- mapply(read_sumstats_noallele, sumstats, snp_cols, z_cols, SIMPLIFY=F)
    } else {
        sumstats <- mapply(read_sumstats, sumstats, snp_cols, a1_cols, a2_cols, z_cols, SIMPLIFY=F)
    }

    sumstats <- lapply(sumstats, as.data.frame, SIMPLIFY=F)

    if (no.allele) {
        z_df <- as.data.frame(Reduce(merge_nomatch_a1a2, sumstats))
    } else {
        z_df <- as.data.frame(Reduce(merge_match_a1a2, sumstats))
    }

    colnames(z_df)[(ncol(z_df)-2):ncol(z_df)] <- paste0("Z", 1:3)
    return(as.data.frame(z_df))
}

z_to_beta <- function(z, n){
    beta <- z / sqrt(n)
    return(beta)
}

## Z table -> BETA, SE, P
# input: read_z(), n, N, c
# c is the weights, r is the correaltion
cal_qt_wtd <- function(z_df, n, N, c, r, eaf) {

    ## order of naive z: predicted-outcome on labeled dataset, predicted-outcome on unlabeled dataset, labeled outcome in labeled dataset
    ## Z1, Z2, Z3 -> beta1, beta2, beta3

    cat("Start correcting:\n")

    # 1 Transform the Z-score into effect size
    beta.lab.yhat <- z_to_beta(z_df$Z1, n)
    beta.unlab.yhat <- z_to_beta(z_df$Z2, N)
    beta.lab.y <- z_to_beta(z_df$Z3, n)

    # 2 Calculate the corrected effect size and its correspond standard error
    beta.pop <- beta.lab.y + c * (beta.unlab.yhat  - beta.lab.yhat)
    se.pop <- sqrt((1 + c^2 - 2 * c * r)/n + c^2/N)
    z.pop <- beta.pop/se.pop
    p.pop <- 2 * pnorm(abs(z.pop), lower.tail = F)
    n.eff <- 1/se.pop^2

    # Transform the standardized scale into the allele scale
    factor <- sqrt(2 * eaf * (1 - eaf))
    df.out <- data.frame(
    BETA = beta.pop / factor,
    SE = se.pop / factor,
    Z = z.pop,
    P = p.pop,
    N.eff = n.eff
  )

    cat("Done.\n\n")

    return(as.data.frame(df.out))
}

cal_qt_wtd_ipw <- function(z_df, n_eff, n, N, c, r, eaf) {

    ## order of naive z: predicted-outcome on labeled dataset, predicted-outcome on unlabeled dataset, labeled outcome in labeled dataset
    ## Z1, Z2, Z3 -> beta1, beta2, beta3

    cat("Start correcting:\n")

    # 1 Transform the Z-score into effect size
    beta.lab.yhat <- z_to_beta(z_df$Z1, n)
    beta.unlab.yhat <- z_to_beta(z_df$Z2, N)
    beta.lab.y <- z_to_beta(z_df$Z3, n_eff)

    # 2 Calculate the corrected effect size and its correspond standard error
    beta.pop <- beta.lab.y + c * (beta.unlab.yhat  - beta.lab.yhat)
    se.pop <- sqrt((1 + c^2 - 2 * c * r)/n + c^2/N)
    z.pop <- beta.pop/se.pop
    p.pop <- 2 * pnorm(abs(z.pop), lower.tail = F)
    n.eff <- 1/se.pop^2

    # Transform the standardized scale into the allele scale
    factor <- sqrt(2 * eaf * (1 - eaf))
    df.out <- data.frame(
    BETA = beta.pop / factor,
    SE = se.pop / factor,
    Z = z.pop,
    P = p.pop,
    N.eff = n.eff
  )

    cat("Done.\n\n")

    return(as.data.frame(df.out))
}

cal_qt_wtd_ovp <- function(z_df, n, N, n_ovp, c, r, eaf) {

    ## order of naive z: predicted-outcome on labeled dataset, predicted-outcome on unlabeled dataset, labeled outcome in labeled dataset
    ## Z1, Z2, Z3 -> beta1, beta2, beta3

    cat("Start correcting:\n")

    # 1 Transform the Z-score into effect size
    beta.lab.yhat <- z_to_beta(z_df$Z1, n)
    beta.unlab.yhat <- z_to_beta(z_df$Z2, N)
    beta.lab.y <- z_to_beta(z_df$Z3, n)

    # 2 Calculate the corrected effect size and its correspond standard error
    beta.pop <- beta.lab.y + c * (beta.unlab.yhat  - beta.lab.yhat)
    se.pop <- sqrt((1 + c^2 - 2 * c * r)/n + c^2/N)
    z.pop <- beta.pop/se.pop
    p.pop <- 2 * pnorm(abs(z.pop), lower.tail = F)
    n.eff <- 1/se.pop^2

    # Transform the standardized scale into the allele scale
    factor <- sqrt(2 * eaf * (1 - eaf))
    df.out <- data.frame(
    BETA = beta.pop / factor,
    SE = se.pop / factor,
    Z = z.pop,
    P = p.pop,
    N.eff = n.eff
  )

    cat("Done.\n\n")

    return(as.data.frame(df.out))
}

cal_qt_1 <- function(z_df, n, N, r, eaf) {
    cal_qt_wtd(z_df, n, N, c = 1, r = r, eaf = eaf)
}

cal_qt_opt <- function(z_df, n, N, r, eaf) {
    cal_qt_wtd(z_df, n, N, c = r * N / (N+n), r = r, eaf = eaf)
}

cal_qt_opt_ipw <- function(z_df, n_eff, n, N, r, eaf) {
    cal_qt_wtd_ipw(z_df, n_eff, n, N, c = r * N / (N+n), r = r, eaf = eaf)
}

cal_qt_opt_ovp <- function(z_df, n, N, n_ovp, r, eaf) {
    cal_qt_wtd_ovp(z_df, n, N, n_ovp, c = r * (N -n_ovp)  / (N+n - 2*n_ovp), r = r, eaf = eaf)
}

effect_std_to_allele <- function(effect_std, eaf){
    effect_allele <- effect_std /(sqrt(2*eaf*(1-eaf)))
}

linear_to_log_or <- function(beta, case_prop, eaf) {
    # log_or: log odds ratio
    # beta: estimated linear regression coefficients with centered genotype and outcome
    # case_prop: propotion of cases in the sample
    # eaf: effect allele frequency
    log_or <- beta * (case_prop * (1 - case_prop) + 
                    0.5 * (1 - 2 * case_prop) * (1 - 2 * eaf) * beta - 
                    (0.084 + 0.9 * case_prop * (1 - 2 * case_prop) * eaf * (1 - eaf))/(case_prop * (1 - case_prop)) * beta^2)^(-1)
    ## Account for variance in outcome: beta = beta * sqrt(var(Y)) = beta * sqrt(case_prop * (1 - case_prop))
    log_or <- log_or * sqrt(case_prop * (1 - case_prop))
    log_or
}

cal_bt_1 <- function(z_df, n, N, r, eaf, case_prop) {

    cat("Start correcting:\n")
    df.out <- cal_qt_1(z_df, n, N, r = r, eaf = eaf)
    df.out$BETA <- linear_to_log_or(df.out$BETA, case_prop, eaf)
    df.out$SE <- linear_to_log_or(df.out$SE, case_prop, eaf)

    cat("Done.\n\n")

    return(as.data.frame(df.out))
}

cal_bt_opt <- function(z_df, n, N, r, eaf, case_prop) {

    cat("Start correcting:\n")
    df.out <- cal_qt_opt(z_df, n, N, r = r, eaf = eaf)
    df.out$BETA <- linear_to_log_or(df.out$BETA, case_prop, eaf)
    df.out$SE <- linear_to_log_or(df.out$SE, case_prop, eaf)

    cat("Done.\n\n")
    return(as.data.frame(df.out))
}