\name{gskat_seq_opt}
\alias{gskat_seq_opt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Perform GEE_SKAT Optimal Test
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gskat_seq_opt(y, XC, Z, ID, NB=10000, impute.method = "fixed", SNP.weights = NULL, w_a = 1, w_b = 25, resampling = TRUE, pw = "Rade", Uc = TRUE, sW = FALSE, np = 10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
binary phenotype coded as 0, 1
}
  \item{XC}{
covaraite matrix, not including the intercept column
}
  \item{Z}{
SNP genotypes coded 0/1/2 (minor allele count).
}
  \item{ID}{
Pedigree ID matrix, including Family ID (FID) and Individual ID (IID)
}
  \item{NB}{
number of comparisons to get p-value, default=10000
}
  \item{impute.method}{
default is fixed method i.e. fill with means
}
  \item{SNP.weights}{
If NULL, the default beta (1,25) desensity will be used, or a custimoized weight vector may be used
}
  \item{w_a}{
The first parameter of the beta density in the weight function
}
  \item{w_b}{
The second parameter of the beta density in the weight function
}
  \item{resampling}{
If TRUE, resampling will be applied
}
  \item{pw}{
r.v. used in the perturbation, "Norm"=Normal , "Rade"=Rademacher
}
  \item{Uc}{
Score centered or not
}
  \item{sW}{
standardize weights
}
  \item{np}{
number of perturbed samples. default=10000
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{pfinal}{Returns p-value}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Xuefeng Wang
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (y, XC, Z, ID, NB, impute.method = "fixed", SNP.weights = NULL, 
    w_a = 1, w_b = 25, resampling = TRUE, pw = "Rade", Uc = TRUE, 
    sW = FALSE, np = 10000) 
{
    maf <- apply(Z, 2, function(x) {
        sum(x, na.rm = T)/((length(x) - sum(is.na(x))) * 2)
    })
    if (impute.method == "fixed") 
        Z <- apply(Z, 2, function(x) {
            x[is.na(x)] = mean(x, na.rm = TRUE)
            return(x)
        })
    Z <- scale(Z, scale = F)
    if (sum(XC[, 1] != 1) != 0) {
        XC = cbind(1, data.matrix(XC))
    }
    Z <- data.matrix(Z)
    y <- as.vector(y)
    FID <- as.factor(ID$FID)
    FID.old <- FID
    FID <- FID[order(FID.old)]
    y <- y[order(FID.old)]
    XC <- XC[order(FID.old), ]
    Z <- data.matrix(Z[order(FID.old), ])
    ID <- ID[order(FID.old), ]
    p = ncol(Z)
    q = ncol(XC) - 1
    if (is.null(SNP.weights)) {
        if (length(maf) == 1) {
            W = 1
        }
        else {
            W <- diag(dbeta(maf, w_a, w_b))
        }
        if (sW == TRUE) {
            W <- W/max(W)
        }
    }
    else {
        W <- diag(SNP.weights)
    }
    null <- fit_FSKAT_IC(y, XC, FID)
    alpha <- null[[1]]
    mu <- as.vector(plogis(XC \%*\% alpha))
    U = t(Z) \%*\% t(t((y - mu)))
    V = diag((mu) * (1 - mu))
    Covy <- tcrossprod(y - mu)
    Covy = Covy * blockMatrixDiagonal(lapply(split(FID, FID), 
        function(alist) {
            matrix(1, length(alist), length(alist))
        }))
    XZ <- cbind(XC, Z)
    B <- t(XZ) \%*\% Covy \%*\% XZ
    A <- t(XZ) \%*\% V \%*\% XZ
    Azx <- A[(q + 1 + 1):ncol(A), 1:(q + 1)]
    Axx <- A[1:(q + 1), 1:(q + 1)]
    C <- cbind(-Azx \%*\% solve(Axx), diag(p))
    rho = seq(0, 1, 0.1)
    TS = rho
    pval = rho
    pvalb = rho
    onev = matrix(1, nrow(W), 1)
    for (i in 1:11) {
        TS[i] = t(U) \%*\% ((1 - rho[i]) * W + rho[i] * tcrossprod(sqrt(W) \%*\% 
            onev)) \%*\% U
    }
    Bsqrt = mysqrt(B)
    BC <- Bsqrt \%*\% t(C) \%*\% ((1 - rho[i]) * W + rho[i] * tcrossprod(sqrt(W) \%*\% 
        onev)) \%*\% C \%*\% Bsqrt
    Lamda <- eigen(BC, only.values = T)$values
    results <- davies(TS, Lamda, rep(1, length(Lamda)))
    results
    if (resampling == TRUE | results$Qq > 1 | results$ifault == 
        1) {
        a = diag((mu) * (1 - mu))
        R <- diag(dim(a)[1])
        iV <- solve(sqrt(a) \%*\% (R) \%*\% sqrt(a))
        y_mu_t = t(t(y - mu))
        ZaV <- t(Z) \%*\% a \%*\% iV
        ni <- do.call(c, lapply(split(ID$FID, ID$FID), function(x) {
            length(x)
        }))
        Ub.boot <- matrix(ncol = np, nrow = nrow(U))
        n = length(unique(ID$FID))
        if (pw == "Norm") {
            Ub.boot <- apply(Ub.boot, 2, function(x) {
                res.pert <- y_mu_t * rep(rnorm(n), ni)
                return(ZaV \%*\% res.pert)
            })
        }
        else if (pw == "Rade") {
            Ub.boot <- apply(Ub.boot, 2, function(x) {
                res.pert <- y_mu_t * rep(sample(c(1, -1), n, 
                  replace = T), ni)
                return(ZaV \%*\% res.pert)
            })
        }
        if (Uc == TRUE) {
            Ub.boot <- apply(Ub.boot, 1, scale, scale = FALSE)
        }
        else if (length(maf) == 1) {
            Ub.boot <- t(t(Ub.boot))
        }
        else {
            Ub.boot <- t(Ub.boot)
        }
        MEANTS = rep(1, 11)
        VARTS = rep(1, 11)
        DF = rep(1, 11)
        for (i in 1:11) {
            Ts_boot <- apply(Ub.boot, 1, function(x) {
                t(x) \%*\% ((1 - rho[i]) * W + rho[i] * tcrossprod(sqrt(W) \%*\% 
                  onev)) \%*\% x
            })
            MEANTS[i] = mean(Ts_boot)
            VARTS[i] = var(Ts_boot)
            DF[i] = 12/kurtosis(Ts_boot)
            df = 12/kurtosis(Ts_boot)
            if (df < 0) {
                df = 100
            }
            mu_Ts <- sum(diag(B \%*\% t(C) \%*\% W \%*\% C))
            if (df < 0) {
                pval = NA
            }
            else {
                pval[i] = 1 - pchisq((TS[i] - mean(Ts_boot)) * 
                  sqrt(2 * df)/sqrt(var(Ts_boot)) + df, df = df)
            }
        }
    }
    else {
        Ts_boot = mu_Ts = var_theory = df = pval = NA
    }
    pmin = min(pval)
    Ub.h <- matrix(ncol = NB, nrow = nrow(U))
    Ub.h <- apply(Ub.h, 2, function(x) {
        res.pert <- y_mu_t * rep(rnorm(n), ni)
        return(ZaV \%*\% res.pert)
    })
    pminb = rep(1, NB)
    for (j in 1:NB) {
        for (i in 1:11) {
            TSh = t(Ub.h[, j]) \%*\% ((1 - rho[i]) * W + rho[i] * 
                tcrossprod(sqrt(W) \%*\% onev)) \%*\% Ub.h[, j]
            pvalb[i] = 1 - pchisq((TSh - MEANTS[i]) * sqrt(2 * 
                DF[i])/sqrt(VARTS[i]) + DF[i], df = DF[i])
        }
        pminb[j] = min(pvalb)
    }
    pfinal = sum(pminb < pmin)/NB
    return(pfinal)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
