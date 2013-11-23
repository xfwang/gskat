\name{gskat_score}
\alias{gskat_score}

\title{
Perform GEE_SKAT Score Test
}
\description{
This function performs family based marker set based association (score test) based on the GEE_SKAT.
}
\usage{
gskat_score(pedDat,F1=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pedDat}{
Input Family Data:  
A list including four data matrix: ID, y, X, Z. 
ID: Pedigree ID matrix, having four columns are the Family ID (FID), Inidvidual ID (IID),  
Paternal ID (FAT),  Maternal ID (MAT), respectively.
y: binary phenotype coded as 0, 1
X: covaraite matrix, including the intercept column
Z: SNP genotypes coded 0/1/2 (minor allele count). 	  
}
  \item{F1}{
If TRUE, use identity working corr. matrix	
}

}

\value{
	\item{pval}{Returns p-value}
  	\item{ifault}{fault indicator from the davies method: 0: no error, 1: requested 
  	accuracy could not be obtained, 2: round-off error possibly significant, 3:
    invalid parameters, 4: unable to locate integration parameters}
}
%% \references{
%% ~put references to the literature/web site here ~
%% }
\author{
Xuefeng Wang
}


\examples{
		gskat_score(gdata)
}

