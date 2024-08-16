rm(list = ls())
library(data.table)
library(openxlsx)
library(janitor)
library(rpart)
library(softImpute)
library(caret)
library(optparse)

option_list = list(
    make_option("--t", action="store", type='character')
)

opt = parse_args(OptionParser(option_list=option_list))

t <- opt$t

# Y
pheno <- fread("./data/train/DXA/dxa_bone_size_mineral_density.txt.gz")
covar <- fread("./Resource/Phenotype/EUR_covar.txt.gz")
pheno.indep <- pheno[pheno$IID %in% covar$IID, ]
covar <- covar[match(pheno.indep$IID, covar$IID),]

# Z
pred <- fread(paste0("./data/train/DXA_BMD/",t,"_eur.Z.txt"))
colnames(pred)[1] <- "IID"
pred <- pred[match(pheno.indep$IID, pred$IID), ]

pred_tmp <- cbind(pred, as.data.frame(pheno.indep)[, t, drop = F])
colnames(pred_tmp)[ncol(pred_tmp)] <- "y"
pred_tmp.1 <- pred_tmp[,"IID",drop=F]
for (i in 2:ncol(pred_tmp)) {
    c <- covar
    c$est <- as.data.frame(pred_tmp)[, i]
    model <- lm(est ~ Sex+Year+SexYear+Chip+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20, data = c)
    res <- rep(NA, nrow(pred_tmp))
    res[as.integer(names(residuals(model)))] <- residuals(model)
    pred_tmp.1 <- cbind(pred_tmp.1, res)
}
colnames(pred_tmp.1)[2:ncol(pred_tmp.1)] <- paste0(colnames(pred_tmp)[2:ncol(pred_tmp)], ".res")
pred_tmp.2 <- pred_tmp.1[as.vector(!is.na(pred_tmp.1[, ncol(pred_tmp.1), with = F])), ]
pred_tmp.2 <- clean_names(pred_tmp.2)
dat <- pred_tmp.2


# Cross-fitting
cv_folds <- createFolds(dat$y_res, k = 10)  # Assuming you're doing 10-fold CV
test_all <- c()
r_vec <- c()

set.seed(1234)
for (i in 1:10) {
    print(i)
    # i <- 1
    # Split the data into training and test based on the folds
    train <- as.matrix(dat[-cv_folds[[i]], -1])
    test_tmp <- dat[cv_folds[[i]],]
    test <- as.matrix(test_tmp[, -1])
    test[, "y_res"] <- NA
    incomplete <- as.matrix(rbind(train, test))

    ## softImpute
    lambda.len <- 100
    maxrank <- T
    fixed.maxrank <- ncol(incomplete)-1
    verbose <- F
    lam0 <- lambda0(incomplete)
    lamseq <- exp(seq(from=log(lam0+.2),to=log(.001),length=lambda.len))
    ranks <- as.integer( lamseq )
    rank.max <- ifelse(maxrank, 2, fixed.maxrank )
    warm <- NULL
    for(j in seq(along=lamseq)){
        if( verbose ) cat( j, ' ' )
        out <- softImpute(x=incomplete, lambda=lamseq[j], rank=rank.max, warm=warm, maxit=1000)
        complete <- complete(incomplete, out)
        warm <- out
        if( maxrank ){
            ranks[j]  <- sum(round(out$d,4)>0)
            rank.max  <- min( ranks[j]+2, fixed.maxrank ) ### grows by at most 2, bounded by P/2
        }
        if( verbose ) cat( '\n' )
    }
    test <- cbind(test_tmp, complete[(nrow(train)+1):nrow(complete), "y_res"])
    colnames(test)[ncol(test)] <- "y_hat"
    r_vec <- c(r_vec, cor(test$y_res, test$y_hat, use = "pairwise.complete.obs"))
    test_all <- rbind(test_all, test)
}

r_vec

# Prediction on all data:
pred_tmp.3 <- clean_names(as.data.frame(pred_tmp.1)[, -ncol(pred_tmp.1)])
pred <- as.matrix(clean_names(as.data.frame(pred_tmp.1)[, -1]))
lambda.len <- 100
maxrank <- T
fixed.maxrank <- ncol(pred)-1
verbose <- F
lam0 <- lambda0(pred)
lamseq <- exp(seq(from=log(lam0+.2),to=log(.001),length=lambda.len))
ranks <- as.integer( lamseq )
rank.max <- ifelse(maxrank, 2, fixed.maxrank )
warm <- NULL
for(j in seq(along=lamseq)){
    if( verbose ) cat( j, ' ' )
    out <- softImpute(x=pred, lambda=lamseq[j], rank=rank.max, warm=warm, maxit=1000)
    pred_complete <- complete(pred, out)
    warm <- out
    if( maxrank ){
        ranks[j]  <- sum(round(out$d,4)>0)
        rank.max  <- min( ranks[j]+2, fixed.maxrank ) ### grows by at most 2, bounded by P/2
    }
    if( verbose ) cat( '\n' )
}
pred_all <- cbind(pred_tmp.3[, 1], as.data.frame(pred_complete))
colnames(pred_all)[ncol(pred_all)] <- "y_hat"
colnames(pred_all)[1] <- "iid"
pred_all_out <- pred_all[!(pred_all$iid %in% test_all$iid), ]

fwrite(pred_all_out, paste0("./data/train/DXA_BMD/",t,".unlab_eur.pop.txt"), sep = "\t", quote = F, row.names = F, col.names = T)
fwrite(test_all, paste0("./data/train/DXA_BMD/",t,".lab_eur.pop.txt"), sep = "\t", quote = F, row.names = F, col.names = T)
