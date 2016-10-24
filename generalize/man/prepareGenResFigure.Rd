\name{prepareGenResFigure}
\alias{prepareGenResFigure}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prepare a figure illustrating generalization testing results
}
\description{
The function takes vectors of effect sizes and standard errors from study1 and study2, and generalization 
  testing results in the form of r-values  and generalization status, corresponding to a set of SNPs 
  selected for follow-up after association testing in study1. 
   The function creates a figure to visualize these results, and the observed effect sizes and confidence intervals. 
  The confidence intervals are scaled to control the false coverage rate (FCR) at the 95% level in both studies. 
}
\usage{
prepareGenResFigure(snpID, study1.beta, study1.se, study2.beta, study2.se, 
	generalized, gen.rvals, study1.n.test, output.file, study1.name = "Study1", 
	study2.name = "Study2", make.title = FALSE, legend.position = c(0.9, 0.9), 
	CI.line.width = 1.7, x.ticks.size = 15, plot.width = 8, plot.height = 8, 
	trait = NULL, plot.title = NULL, add.stars = NULL, alpha.level.follow = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{snpID}{
  a vector of snp identifiers (e.g. rsIDs)
}
  \item{study1.beta}{
  Numeric vector of the effect estimates from study1.
}
  \item{study1.se}{
 Numeric vector of the SE of the effect estimates from study1.
}
  \item{study2.beta}{
    Numeric vector of the effect estimates from study2.
}
  \item{study2.se}{
  Numeric vector of the SE of the effect estimates from study2.
}
  \item{generalized}{
  Generalization status, return by the function testGeneralization
}
  \item{gen.rvals}{
  r-values from generalization testing. Returned by the function testGeneralization.
}
  \item{study1.n.test}{
 Number of SNPs tested in study 1.
}
  \item{output.file}{
 Name of file for the figure.
}
  \item{study1.name}{
  A string, name of study1.
}
  \item{study2.name}{
  A string, name of study2.
}
  \item{make.title}{
  Should the plot have a title?
}
  \item{legend.position}{
  Legend position.
}
  \item{CI.line.width}{
  the width of the line of the confidence intervals.
}
  \item{x.ticks.size}{
  The size of the snp identifiers strings given at the bottom of the figure (at the X-axis).
}
  \item{plot.width}{
  Figure width.
}
  \item{plot.height}{
 Figure height.
}
  \item{trait}{
  An optional string with the name of the trait associated (or not) with the SNPs. Used in the figure title if given. 
}
  \item{plot.title}{
  An optional string of a spefic title
}
  \item{add.stars}{
  An optional numerical value, corrsponding to a p-value threshold. If a SNP from study2 has a two-sided 
   p-value that is smaller than this threshold, than a star is plotted next to the estimated effect of the appropriate SNP.
}
  \item{alpha.level.folow}{
  Significance level for plotting confidence interval. Defult level is 0.05 (plots 95\% confidence intervals)
}

}
\details{
A figure is created in output.file.
}
\examples{
data(dat)
head(dat)

## first match effects on alleles:
dat.matched <- matchEffectAllele(dat$rsID, study2.effect = dat$study2.beta, 
	study1.alleleA = dat$study1.alleleA, study2.alleleA = dat$study2.alleleA, 
	study1.alleleB = dat$study1.alleleB, study2.alleleB = dat$study2.alleleB)

head(dat.matched)

dat$study2.beta <- dat.matched$study2.effect

## test for generalization:

gen.res <- testGeneralization(dat$rsID, dat$study1.pval, dat$study2.pval, 
	dat$study1.n.test[1], study1.effect = dat$study1.beta, 
	study2.effect = dat$study2.beta, directional.control = TRUE, 
	control.measure = "FDR" )

## prepare figure
require(ggplot2)
require(gridExtra)
require(RColorBrewer)

## prepare figure
prepareGenResFigure(dat$rsID, dat$study1.beta, dat$study1.se, dat$study2.beta, 
	dat$study2.se, gen.res$generalized, gen.res$gen.rvals, dat$study1.n.test[1], 
	output.file = paste0(getwd(), "/Generalization_example.pdf"), 
	study1.name = "Study1", study2.name = "Study2")

## marking genome-wide significant hits in study 2 with red stars 
## (here nothing is genome-wide significant so there will be no difference!)
prepareGenResFigure(dat$rsID, dat$study1.beta, dat$study1.se, dat$study2.beta, 
	dat$study2.se, gen.res$generalized, gen.res$gen.rvals, dat$study1.n.test[1], 
	output.file = paste0(getwd(), "/Generalization_example.pdf"), 
	study1.name = "Study1", study2.name = "Study2", add.stars = 5e-8) 

}
