\name{score_FSKAT_IC_pertu}
\alias{score_FSKAT_IC_pertu}

\title{
GEE_SKAT for rare variants
}
\description{
This function performs GEE_SKAT score test based on the davies and perturbation method
}
\usage{
score_FSKAT_IC_pertu(y,XC,Z,ID,w_a=1,w_b=5,pw="Rade",Uc=FALSE,sW=TRUE,np=10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{y}{binary phenotype coded as 0, 1}
	\item{XC}{covaraite matrix, WITHOUT the intercept column.}
	\item{Z}{SNP genotypes coded 0/1/2 (minor allele count).}
	\item{ID}{Pedigree ID matrix. The four columns are the Family ID (FID), Individual ID (IID),  
Paternal ID (FAT),  Maternal ID (MAT), respectively.}
	\item{w_a}{The first parameter of the beta density in the weight function}
	\item{w_b}{The second parameter of the beta density in the weight function}
	\item{pw}{perturbatio method: "Rade"=Rademacher; "Norm"=Normal}
	\item{Uc}{Score centered or not}
	\item{sW}{standardize weights or not}
	\item{np}{No. of perturbed samples}	
}

\value{
	\item{pval_davies}{Returns p-value from davies method}
  	\item{ifault}{fault indicator from the davies method: 0: no error, 1: requested accuracy
    could not be obtained, 2: round-off error possibly significant, 3:
    invalid parameters, 4: unable to locate integration parameters}
   	\item{PM}{Empirical mean of Ts based on perturbation}
   	\item{mu_Ts}{Theoretical mean of Ts}
   	\item{PV}{Empirical variance of Ts based on perturbation}
   	\item{var_Ts}{Theoretical variance of Ts}
   	\item{pval_pert}{Returns p-value from perturbation method}
}
%% \references{
%% ~put references to the literature/web site here ~
%% }
\author{
Xuefeng Wang
}

\seealso{
score.RFAM_IC_burden
}
\examples{
	score_FSKAT_IC_pertu(y=gdata$y,XC=gdata$X[,-1],Z=gdata$Z,ID=gdata$ID)
}

