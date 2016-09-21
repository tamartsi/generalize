\name{testGeneralization}
\alias{testGeneralization}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Test for generalization of associations from study1 to study2 
}
\description{
  The function calls another function that computes r-values given two vectors of p-values from primary 
  (study1) and follow-up (study2) studies and determine generalization (TRUE/FALSE) for each of the associations. The
 computes either FDR or FWER r-values.  Directional control is applied by default, so that one-sided p-values are computed
  and used.only and associations with the same sign in the two studies are generalized. 
}
\usage{
testGeneralization(snpID, study1.pval, study2.pval, study1.n.test, study1.direction = NULL,  study2.direction = NULL, directional.control = TRUE, control.measure = "FDR", q = 0.05, l00 = 0.8, c2 = 0.5, variation = c("none","use.m.star","use.t"), tt = NULL,  verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{snpID}{
  a vector of snp identifiers (e.g. rsIDs)
}
  \item{study1.pval}{
  Numeric vector of the two sided p-values from study1.
}

  \item{study2.pval}{
  Numeric vector of the two sided p-values from study2.
}
  \item{study1.n.test}{
 Number of SNPs tested in study 1.
}
  \item{study1.direction}{
   An optional numeric vector of the directions of effects of the SNPs from study1. These could also be the effect estimates themselves, 
  i.e. any numeric vector with signes matching the signes of the effects in study1. 
	Required for directional control using 1-sided p-values. 
}
  \item{study2.direction}{
   An optional numeric vector of the directions of effects of the SNPs from study2. These could also be the effect estimates themselves, 
  i.e. any numeric vector with signes matching the signes of the effects in study2.
	Required for directional control using 1-sided p-values. 
}
  \item{directional.control}{
  should directional controlled by applied? i.e., a SNP could be generalized only if it has 
  the same direction of association in both study1 and study2. 
}
  \item{control.measure}{
 Either "FDR" or "FWER", according to the desired control of false generalizations.
}
  \item{q}{
  The level or false generalization (e.g. control FDR at the q level or FWER at the q level).
}
  \item{l00}{
 Lower bound of the fraction of features (out of m) with true null hypotheses in both studies.
    For example, for GWAS on the whole genome, the choice of 0.8 is conservative
    in typical applications.
}
  \item{c2}{
 Parameter for relative boost to the p-values from the primary study.
    0.5 (default) is recommended, since was observed in simulations to yield
    similar power to procedure with the optimal value (which is unknown for real data).

}
  \item{variation}{
  When 'use.m.star' is selected m* is used. m* is defined as follows:
   \eqn{m^*=m\sum_{i=1}^{m}\frac{1}{i}}{m*=m*sum(1/i)}.
   When 'use.t' is selected c1 is computed given the threshold tt.
   Both variations guarantee that the procedure that decleares all r-values below  q as replicability claims,
   controls the FDR at level q, for any type of dependency of the p-values in the primary study.
   default is 'none'.
}
  \item{tt}{
   The selection rule threshold for p-values from the primary study. must be supplied when
   variation 'use.t' is selected.
}
  \item{verbose}{
}
}
\details{
Returns a data frame matching each snpID with its computed r-value and generalization status.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
## first match effects on alleles:

dat.matched <- matchEffectAllele(dat$rsID, study1.direction = dat$study1.beta, study2.direction = dat$study2.beta, study1.alleleA = dat$study1.alleleA, study2.alleleA = dat$study2.alleleA, study1.alleleB = dat$study1.alleleB, study2.alleleB = dat$study2.alleleB)


head(dat.matched)

dat$study2.beta <- dat.matched$study2.direction

## test for generalization:

gen.res <- testGeneralization(dat$rsID, dat$study1.pval, dat$study2.pval, dat$study1.n.test[1], study1.direction = dat$study1.beta, study2.direction = dat$study2.beta, directional.control = TRUE, control.measure = "FDR" )
}
