\name{gskat_seq_cont}
\alias{gskat_seq_cont}

\title{
Perform GEE_KM Score Test with sequencing data for continuous traits
}
\description{
This function performs family based marker set based association (score test) based on the GEE_KM.
}
\usage{
gskat_seq_cont(y,XC,Z,ID,impute.method="fixed",SNP.weights=NULL,
	w_a=1,w_b=25,resampling=TRUE,pw="Rade",Uc=TRUE,sW=FALSE,np=10000)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{y}{binary phenotype coded as 0, 1}
	\item{XC}{covaraite matrix, including the intercept column}
	\item{Z}{SNP genotypes coded 0/1/2 (minor allele count)}
	\item{ID}{ID: Pedigree ID matrix, having four columns are the Family ID (FID), Individual ID (IID),  
	Paternal ID (FAT),  Maternal ID (MAT), respectively.}
	\item{impute.method}{default is fixed method i.e. fill with means}
	\item{SNP.weights}{If NULL, the default beta (1,25) desensity will be used, or a custimoized weight vector may be used}
	\item{w_a}{The first parameter of the beta density in the weight function}
	\item{w_b}{The second parameter of the beta density in the weight function}
	\item{resampling}{If TRUE, resampling will be applied}
	\item{pw}{r.v. used in the perturbation, "Norm"=Normal , "Rade"=Rademacher}
	\item{Uc}{Score centered or not}
	\item{sW}{standardize weights or not}
	\item{np}{No. of perturbed samples}
}

\value{
	\item{p1}{Returns asymptotic p-value}
	\item{p2}{Returns resampling p-value}
}
%% \references{
%% ~put references to the literature/web site here ~
%% }
\author{
Xuefeng Wang
}

\examples{
		# gskat_seq_cont(y,XC=X,Z,ID)
}

