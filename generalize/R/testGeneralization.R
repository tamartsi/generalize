testGeneralization <-
function(snpID,  study1.pval, study2.pval, study1.n.test, study1.effect = NULL,  study2.effect = NULL, directional.control = TRUE, control.measure = "FDR", q = 0.05, l00 = 0.8, c2 = 0.5, variation = c("none","use.m.star","use.t"), tt = NULL,  verbose = FALSE ){
	
	stopifnot(length(snpID) == length(study1.pval),  length(snpID) == length(study2.pval), length(snpID) < study1.n.test, l00>=0, l00 < 1, is.element(control.measure, c("FDR", "FWER")))
	
	inds.na <- unique(c(which(is.na(snpID)), which(is.na(study1.pval)), which(is.na(study2.pval))))
	if (length(inds.na) > 0){
		message(paste0(length(inds.na), " SNPs are missing either SNP ID, study1.pval, or study2.pval. Removing these SNPs."))
		snpID <- snpID[-inds.na]
		study1.pval <- study1.pval[-inds.na]
		study2.pval <- study2.pval[-inds.na]
		if (!is.null(study1.effect)){
			study1.effect <- study1.effect[-inds.na]
			study2.effect <- study2.effect[-inds.na]
		}
	}
	
	## checks if directional control is performed:
	if (directional.control){
		if (is.null(study1.effect) | is.null(study2.effect)) stop("effects (or directions as 1 or -1) must be given for directional control")
		if (sum(is.na(study1.effect)) > 0) stop("Some study1 effects are missing")
		if (sum(is.na(study2.effect)) > 0) stop("Some study2 effects are missing")	
		
	}
		
	
	message(paste0("Controlling ", control.measure, "at the ", q, " level"))
	
	if (directional.control){
		message("Generating one-sided p-values guided by study1's directions of effects...")
		study1.pval.1sided <- study1.pval/2
		study2.pval.1sided <- study2.pval/2
		inds.diff.sign <- which(sign(study1.effect) != sign(study2.effect))
		if (length(inds.diff.sign) > 0 ) study2.pval.1sided[inds.diff.sign] <- 1- study2.pval.1sided[inds.diff.sign]
		
		study1.pval <- study1.pval.1sided
		study2.pval <- study2.pval.1sided
	}
	
	
	message(paste("Calcluating", control.measure, "r-values..."))
	gen.rvals <- r.value(p1 = study1.pval, p2 = study2.pval, m = study1.n.test, control.measure = control.measure, c2 = c2, l00 = l00, variation = variation, tt = tt, Q=q )
	
	generalized <- (gen.rvals <= q)
	
	return(data.frame(snpID, gen.rvals, generalized , stringsAsFactors = F))
	
}
