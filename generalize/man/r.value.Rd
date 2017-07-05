\name{r.value}
\alias{r.value}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 r-value computation 
}
\description{
The function computes r-values given two vectors of p-values from primary and
    follow-up studies. r-values assess the False Discovery Rate (FDR) of repilcability
    claims across the primary and follow-up studies. 
  This is a function from Ruth Heller, adapted to compute FWER r-values in addition to FDR r-values. 
}
\usage{
r.value(p1, p2, m, c2 = 0.5, control.measure = "FDR", l00 = 0.8 , variation = c("none","use.m.star","use.t"), tt = NULL, Q= 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p1}{
  Numeric vector of the p-values from study 1
}
  \item{p2}{
  Numeric vector  of the p-values from study 2
}
  \item{m}{
  Number of features examined in the primary study.
}
  \item{c2}{
  Parameter for relative boost to the p-values from the primary study.
    0.5 (default) is recommended, since was observed in simulations to yield
    similar power to procedure with the optimal value (which is unknown for real data).

}
  \item{control.measure}{
  A sting, either FDR or FWER, depending on the desired measure of control on false generalizations. 
  }
  \item{l00}{
 Lower bound of the fraction of features (out of m) with true null hypotheses in both studies.
    For example, for GWAS on the whole genome, the choice of 0.8 is conservative
    in typical applications.
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
  \item{Q}{
  The level or false generalization (e.g. control FDR at the q level or FWER at the q level).
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
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
Ruth Heller, Shay Yaacoby (shay66@gmail.com), small adaptation by Tamar Sofer. 
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# General example from Ruth Heller's website: 
  pv <- read.csv("http://www.math.tau.ac.il/~ruheller/Software/CrohnExample.csv")
  rv <- r.value(p1=pv$p1,p2=pv$p2,m=635547,c2=0.5,l00=0.8)
  rv2 <- r.value(p1=pv$p1,p2=pv$p2,m=635547,c2=0.5,l00=0.8,variation="use.t",tt=1e-5)

#### in this package, the function is called by testGenerelization.
}

