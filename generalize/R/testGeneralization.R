testGeneralization <-
function(snpID,  study1.pval, study2.pval, study1.n.test, study1.direction = NULL,  study2.direction = NULL, directional.control = TRUE, control.measure = "FDR", q = 0.05, l00 = 0.8, c2 = 0.5, variation = c("none","use.m.star","use.t"), tt = NULL,  verbose = FALSE ){
	
	stopifnot(length(snpID) == length(study1.pval),  length(snpID) == length(study2.pval), length(snpID) < study1.n.test, l00>=0, l00 < 1, is.element(control.measure, c("FDR", "FWER")))
	
	if (!is.null(study1.direction)) if (length(snpID) != length(study1.direction)) stop("Number of SNP identifier different than number of SNP directions of effects for study 1")
	if (!is.null(study2.direction)) if (length(snpID) != length(study2.direction)) stop("Number of SNP identifier different than number of SNP directions of effects for study 2")
	if (directional.control & (is.null(study1.direction) | is.null(study2.direction) )) stop("effect directions must be given for directional control")
	
	if (directional.control){
		message("Generating one-sided p-values guided by study1's directions of effects...")
		study1.pval.1sided <- study1.pval/2
		study2.pval.1sided <- study2.pval/2
		inds.diff.sign <- which(sign(study1.direction) != sign(study2.direction))
		if (length(inds.diff.sign) > 0 ) study2.pval.1sided[inds.diff.sign] <- 1- study2.pval.1sided[inds.diff.sign]
		
		study1.pval <- study1.pval.1sided
		study2.pval <- study2.pval.1sided
	}
	
	
	message(paste("Calcluating", control.measure, "r-values..."))
	gen.rvals <- r.value(p1 = study1.pval, p2 = study2.pval, m = study1.n.test, control.measure = control.measure, c2 = c2, l00 = l00, variation = variation, tt = tt, Q=q )
	
	generalized <- (gen.rvals <= q)
	
	return(data.frame(snpID, gen.rvals, generalized , stringsAsFactors = F))
	
}
