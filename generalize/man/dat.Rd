\name{dat}
\alias{dat}
\docType{data}
\title{
A data set of previouly published SNP-platelet count associations in Gieger et al., 2011, and association results for the same SNPs in the HCHS/SOL. 
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("dat")}
\format{
  A data frame with 42 observations on the following 15 variables.
  \describe{
    \item{\code{rsID}}{SNP identifier}
    \item{\code{chromosome}}{chromosome}
    \item{\code{position}}{position in human build 37}
    \item{\code{study1.alleleA}}{effect Allele in study 1}
    \item{\code{study1.alleleB}}{other Allele in study 1}
    \item{\code{study1.beta}}{estimated effect size is study 1}
    \item{\code{study1.se}}{SE of beta in study 1}
    \item{\code{study1.pval}}{two-sided p-value in study 1}
    \item{\code{study1.n.test}}{number of tested SNPs in study 1}
    \item{\code{study2.alleleA}}{effect Allele in study 2}
    \item{\code{study2.alleleB}}{other Allele in study 2}
    \item{\code{study2.beta}}{estimated effect size is study 2}
    \item{\code{study2.se}}{SE of beta in study 2}
    \item{\code{study2.pval}}{two-sided p-value in study 2}
    \item{\code{Ref}}{Reference for the platelet GWAS paper from which the study1 results are taken - Gieger,2011.}
  }
}
\details{
This data set provides the top results from a platelet GWAS reported in Gieger et al., 2011,
and the association testing results for the same SNP in the HCHS/SOL stuy. 
These results will soon be published in the American Journal of Human Genetics (as of Dec 2015), by Schick et al.
}
\source{
Gieger, Christian, et al. "New gene functions in megakaryopoiesis and platelet formation." Nature 480.7376 (2011): 201-208,
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(dat)
head(dat)
}
\keyword{datasets}
