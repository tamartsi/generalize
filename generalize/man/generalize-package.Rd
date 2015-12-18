\name{generalize-package}
\alias{generalize-package}
\alias{generalize}
\docType{package}
\title{
\packageTitle{generalize}
}
\description{
\packageDescription{generalize}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{generalize}
\packageIndices{generalize}
Code used to perform generalization testing for a set of association from study 1 (the discovery study) to study 2 (the follow-up, generalization study).
}
\author{
\packageAuthor{generalize}

Maintainer: Tamar Sofer \packageMaintainer{generalize}
}
\references{
~~ Literature or other references for background information ~~
}
~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation directory ~~
\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
data(dat)
head(dat)
## first match effects on alleles:
dat.matched <- matchEffectAllele(dat$rsID, study1.direction = dat$study1.beta, study2.direction = dat$study2.beta, study1.alleleA = dat$study1.alleleA, study2.alleleA = dat$study2.alleleA, study1.alleleB = dat$study1.alleleB, study2.alleleB = dat$study2.alleleB)


dat$study2.beta <- dat.matched$study2.direction

## test for generalization:

gen.res <- testGeneralization(dat$rsID, dat$study1.pval, dat$study2.pval, dat$study1.n.test[1], study1.direction = dat$study1.beta, study2.direction = dat$study2.beta, directional.control = TRUE, control.measure = "FDR" )
head(gen.res)

require(ggplot2)
require(gridExtra)
require(RColorBrewer)


## prepare figure
prepareGenResFigure(dat$rsID, dat$study1.beta, dat$study1.se, dat$study2.beta, dat$study2.se, gen.res$generalized, gen.res$gen.rvals, dat$study1.n.test[1], output.file = paste0(getwd(), "/Generalization_example.pdf"), study1.name = "Study1", study2.name = "Study2")

## marking genome-wide significant hits in study 2 with red stars (here nothing is genome-wide significant so there will be no difference!)
prepareGenResFigure(dat$rsID, dat$study1.beta, dat$study1.se, dat$study2.beta, dat$study2.se, gen.res$generalized, gen.res$gen.rvals, dat$study1.n.test[1], output.file = paste0(getwd(), "/Generalization_example.pdf"), study1.name = "Study1", study2.name = "Study2", add.stars = 5e-8) 
}
