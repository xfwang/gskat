\name{gskat_score_pert}
\alias{gskat_score_pert}

\title{
Perform GEE_SKAT Score Test (perturbation)
}
\description{
This function performs family based marker set based association (score test) based on the GEE_SKAT.
}
\usage{
gskat_score_pert(pedDat,F1=FALSE,pw="Rade",np=10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pedDat}{
Input Family Data:  
A list including four data matrix: ID, y, X, Z. 
ID: Pedigree ID matrix. The four columns are the Family ID (FID), Individual ID (IID),  
Paternal ID (FAT),  Maternal ID (MAT), respectively.
y: binary phenotype coded as 0, 1
X: covaraite matrix, including the intercept column
Z: SNP genotypes coded 0/1/2 (minor allele count). 	  
}
  \item{F1}{
If TRUE, use identity working corr. matrix	
}
	\item{pw}{
Perturbation method."Rade":Rademacher; "Norm": Normal distribution
}
\item{np}{
number of perturbed samples. default=10000
}
}

\value{
	\item{pval}{Returns p-value}
  	\item{Ts}{Score statistic}
  	\item{mu_Ts}{Theoretical mean of Ts}
  	\item{var_Ts}{Theoretical variance of Ts}
  	\item{PM}{Empirical mean of Ts based on the perturbation}
  	\item{PV}{Empirical variance of Ts based on the perturbation}
}
%% \references{
%% ~put references to the literature/web site here ~
%% }
\author{
Xuefeng Wang
}


\examples{
		gskat_score_pert(gdata)
}

