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
       attach(gdata)
       gskat_seq_opt(y,X[,-1],Z,ID=data.frame(FID=ID$FID,IID=ID$IID))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
