\name{matchEffectAllele}
\alias{matchEffectAllele}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
check and update signs of estimated effect sizes in study 2 by matching alleles between study 1 and study 2
}
\description{
Compare effect alleles and other alleles and change the sign of the effect in study2 in cases were the effect alleles are different between studies. 
}
\usage{
matchEffectAllele(snpID, study1.direction, study2.direction, study1.alleleA, study2.alleleA, study1.alleleB = NULL, study2.alleleB = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{snpID}{
  a vector of snp identifiers (e.g. rsIDs)
}
  \item{study1.direction}{
  Numeric vector of the directions of effects of the SNPs from study1. These could also be the effect estimates themselves, 
  i.e. any numeric vector with signes matching the signes of the effects in study1. 
}
  \item{study2.direction}{
  Numeric vector of the directions of effects of the SNPs from study2. These could also be the effect estimates themselves, 
  i.e. any numeric vector with signes matching the signes of the effects in study2.
}
  \item{study1.alleleA}{
  Character vector providing the effect alleles corresponding to snpID in study1. 
}
  \item{study2.alleleA}{
  Character vector providing the effect alleles corresponding to snpID in study2. 
}
  \item{study1.alleleB}{
  Character vector providing the non-effect (other) alleles corresponding to snpID in study1
}
  \item{study2.alleleB}{
  Character vector providing the non-effect (other) alleles corresponding to snpID in study2.
}
}
\details{
  Returns a data frame matching each snpID with its estimated effect size is study1 and study2, after making sure that 
  both studies use the same effect alleles for all given SNPs. If there is enough information, strand-abgiuouos SNPs flag is also given.
}

\value{
Returns a data frame with SNP identifiers, direction of effects of study1 and study2 after possibly flipping 
the directions of effects of study2 so the same effect allele is used in both studyes. If effect estimates were given in study1.direction and study2.direction, then effect estimates are returned accordingly.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Tamar Sofer
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(dat)
head(dat)
## match effects on alleles:

dat.matched <- matchEffectAllele(dat$rsID, study1.direction = dat$study1.beta, study2.direction = dat$study2.beta, study1.alleleA = dat$study1.alleleA, study2.alleleA = dat$study2.alleleA, study1.alleleB = dat$study1.alleleB, study2.alleleB = dat$study2.alleleB)
head(dat.matched)

}
