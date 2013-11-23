\name{score.RFAM_IC_burden}
\alias{score.RFAM_IC_burden}

\title{
Burden test based on GEE_SKAT
}
\description{
This function performs burden test based on the GEE_SKAT.
}
\usage{
score.RFAM_IC_burden(y,XC,Z,FID)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{y}{binary phenotype coded as 0, 1}
	\item{XC}{covaraite matrix, WITHOUT the intercept column.}
	\item{Z}{SNP genotypes coded 0/1/2 (minor allele count).}
	\item{FID}{Family ID}
}

\value{
	Returns p-value
}
%% \references{
%% ~put references to the literature/web site here ~
%% }
\author{
Xuefeng Wang
}


\seealso{
geeglm_wrap; gee_wrap
}
\examples{
	Z1=genoC(gdata$Z)
	score.RFAM_IC_burden(y=gdata$y,XC=gdata$X[,-1],Z=Z1,FID=gdata$ID[,1])
}

