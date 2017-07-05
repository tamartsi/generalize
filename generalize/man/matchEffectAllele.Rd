\name{matchEffectAllele}
\alias{matchEffectAllele}
update directions of effect of study2 to be according to the effect allele in study1
\title{
check and update signs of estimated effect sizes in study 2 by matching alleles between study 1 and study 2
}
\description{
Compare effect alleles and other alleles and change the sign of the effect in study2 in cases were the effect alleles are different between studies. 
}
\usage{
matchEffectAllele(snpID, study2.effect, study1.alleleA, study2.alleleA, 
	study1.alleleB = NULL, study2.alleleB = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{snpID}{
  a vector of snp identifiers (e.g. rsIDs)
}
  \item{study2.effect}{
  Numeric vector of the SNP effects from study2, or their directions given as 1, -1. 
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
  Returns a data frame matching each snpID with its estimated effect size in study2, after making sure that 
  both studies use the same effect alleles for all given SNPs. If there is enough information, strand-abgiuouos SNPs flag is also given.
}

\value{
Returns a data frame with SNP identifiers, direction of effects of study2 after possibly flipping 
the directions of effects of study2 so the same effect allele is used in both studyes. 
The returned data frame also has an indicator of whether the direction of association was flipped
based on allele matched. If it was, and you would like to provide the correct effect allele 
frequencies associated with the follow-up study effect sizes, make sure to change them 
(by eaf[which(flip)] <- 1-eaf[which(flip)]) for the flipped effect sizes.
If effect estimates were given in study1.direction and study2.direction, then effect estimates are returned accordingly.
}

\author{
Tamar Sofer
}
other sections like Warning with \section{Warning }{....} ~
\examples{
data(dat)
head(dat)
## match effects on alleles:

dat.matched <- matchEffectAllele(dat$rsID, study2.effect = dat$study2.beta, 
	study1.alleleA = dat$study1.alleleA, study2.alleleA = dat$study2.alleleA, 
	study1.alleleB = dat$study1.alleleB, study2.alleleB = dat$study2.alleleB)
	
head(dat.matched)

}
