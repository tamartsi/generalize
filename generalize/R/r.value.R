r.value <-
function (p1, p2, m, c2 = 0.5, control.measure = "FDR", l00= 0.8, variation = c("none","use.m.star","use.t"), tt = NULL, Q = 0.05)
{
  variation <- match.arg(variation)
  
  if (variation != "use.t" & !is.null(tt)) 
    warning("threshold tt is ignored")
  if (variation == "use.t" & !is.null(tt))
  {
    if (tt <= (1-c2)/(1-l00*(1-c2*Q)) * Q/m )
    {
      warning("since t < c(q)q/m, no modification to the original r-value computation was necessary (see section Derivation & Properties in the article)")
      variation <- "none"
    }
    if (tt >= (1-c2)/(1-l00*(1-c2*Q))*Q/(1+sum(1/(1:(m-1)))))
      warning("for the selected threshold t, the 'use.t' variation won't lead\nto more discoveries than the 'use.m.star' variation.")
    
  }
  if (variation == "use.m.star") m <- m*sum(1L/(1L:m)) 
  k <- length(p1) 
  #------- input checks: -------
  if (k!=length(p2))
    stop("p-value vectors must be of equallengths")
  if ((1<=c2)|(c2<=0))
    stop("c2 value should be in the interval (0,1)")
  if ((1<=l00)|(l00<0))
    stop("l00 value should be in the interval [0,1)")
  if (m<k)
    stop("Number of features in the primary stage (m) must be at least as large as the number of features followed-up.")
  if (any(is.na(c(p1,p2))))
    stop("NA's are not allowed as p-values")
  if (any(c(p1,p2)>1) | any(c(p1,p2)<=0))
    stop("p-values must be in the interval (0,1]")
  if (variation == "use.t" & is.null(tt)) 
    stop("specify threshold tt")
  if (!is.element(control.measure, c("FDR", "FWER")))
  	stop("specificy control measure FDR or FWER")  
    
  #------- function definition: compute r-value of given value of x: -------
  
  if (control.measure == "FWER"){
	  r.value.x <- function (x) {
	    
	    c1 <- switch(variation,
	                 none       = (1-c2)/(1-l00*(1-c2*x)),
	                 use.m.star = (1-c2)/(1-l00*(1-c2*x)),
	                 use.t = {
	                   if (tt <= (1-c2)/(1-l00*(1-c2*x)) * x/m)
	                     (1-c2)/(1-l00*(1-c2*x))
	                   else
	                   {
	                     lower <- 1e-6 ; upper <- (1-c2)/(1-l00*(1-c2*x))
	                     f <- function(a)
	                     {
	                       if (tt*m/(a*x) < 10)
	                       {
	                         a*(1+sum(1/(1:(ceiling(tt*m/(a*x)-1))))) - (1-c2)/(1-l00*(1-c2*x))
	                       }
	                       else
	                       {
	                         a*(1-digamma(1)+1/(2*ceiling(tt*m/(a*x)-1))+log(ceiling(tt*m/(a*x)-1))) - (1-c2)/(1-l00*(1-c2*x))
	                       }
	                     }
	                     
	                     c1_sol1 <- uniroot(f,c(lower,upper))
	                     next_step <- tt*m/(ceiling(tt*m/(c1_sol1$root*x))*x)
	                     if (f(next_step)<=0)
	                     {
	                       c1_sol2 <- uniroot(f,c(next_step,upper))
	                       c1_sol2$root
	                     }
	                     else
	                     {
	                       c1_sol1$root 
	                     }
	                   }
	                 })
	    
	    E   <- pmax(m*p1/c1, k*p2/c2)
	    return(E)
	  }
  } else if (control.measure == "FDR"){
  	  r.value.x <- function (x) {
	    c1 <- switch(variation,
	                 none       = (1-c2)/(1-l00*(1-c2*x)),
	                 use.m.star = (1-c2)/(1-l00*(1-c2*x)),
	                 use.t = {
	                   if (tt <= (1-c2)/(1-l00*(1-c2*x)) * x/m)
	                     (1-c2)/(1-l00*(1-c2*x))
	                   else
	                   {
	                     lower <- 1e-6 ; upper <- (1-c2)/(1-l00*(1-c2*x))
	                     f <- function(a)
	                     {
	                       if (tt*m/(a*x) < 10)
	                       {
	                         a*(1+sum(1/(1:(ceiling(tt*m/(a*x)-1))))) - (1-c2)/(1-l00*(1-c2*x))
	                       }
	                       else
	                       {
	                         a*(1-digamma(1)+1/(2*ceiling(tt*m/(a*x)-1))+log(ceiling(tt*m/(a*x)-1))) - (1-c2)/(1-l00*(1-c2*x))
	                       }
	                     }
	                     
	                     c1_sol1 <- uniroot(f,c(lower,upper))
	                     next_step <- tt*m/(ceiling(tt*m/(c1_sol1$root*x))*x)
	                     if (f(next_step)<=0)
	                     {
	                       c1_sol2 <- uniroot(f,c(next_step,upper))
	                       c1_sol2$root
	                     }
	                     else
	                     {
	                       c1_sol1$root 
	                     }
	                   }
	                 })
	    
	    E   <- pmax(m*p1/c1, k*p2/c2)
	    oe  <- order(E, decreasing =TRUE)
	    oer <- order(oe)
	    r   <- cummin((E/rank(E, ties.method= "max"))[oe])[oer]
	    r   <- pmin(r,1)
	    return(r)
	  }	
  }
  
  rv <- rep(NA,length(p2))
  tol = min(c(p1,p2)[c(p1,p2)!=0], 0.0001)
  for (i in 1:length(p2))
  {
    aux <- function(x) return(r.value.x(x)[i]-x)
    if (aux(1)>=0) rv[i] <- 1
    else 
    {
      sol <- uniroot(aux,c(tol,1),tol = tol)
      rv[i] <- sol$root
    }
  }  
  return(rv)
}
