matchEffectAllele <-
function(snpID, study1.direction, study2.direction, study1.alleleA, study2.alleleA, study1.alleleB = NULL, study2.alleleB = NULL){
	stopifnot(length(snpID) == length(study1.direction),  length(snpID) == length(study2.direction), length(snpID) == length(study1.alleleA), length(snpID) == length(study2.alleleA))
	if (!is.null(study1.alleleB)) stopifnot(length(snpID) == length(study1.alleleB))
	if (!is.null(study2.alleleB)) stopifnot(length(snpID) == length(study2.alleleB))	
	
	complement <- function(x, y){
 		if ((x == "A" & y == "T") | (x == "T" & y == "A") | (x == "G" & y == "C") | (x == "C" & y == "G")) return(TRUE) else return(FALSE) 
 	}

	
	
	flip <- rep(NA, length = length(snpID))
	strand.ambiguous <- rep(NA, length = length(snpID))
	
	if (is.null(study1.alleleB) & is.null(study2.alleleB)){
		strand.ambiguous <- NULL
		message("Allele matching based on partial information, only effect alleles are provided")
		for (i in 1:length(snpID)){
			if (study1.alleleA[i] == study2.alleleA[i]) {flip[i] <- FALSE} else{ flip[i] <- TRUE}
		}
	}

	if (is.null(study1.alleleB) & !is.null(study2.alleleB)){
		message("Allele matching based on partial information, other alleles are not provided for study1")
		for (i in 1:length(snpID)){
			if (study1.alleleA[i] == study2.alleleA[i]) {flip[i] <- FALSE} else{
				if (study1.alleleA[i] == study2.alleleB[i]) {flip[i] <- TRUE} else{
					if (complement(study1.alleleA[i], study2.alleleA[i])) {flip[i] <- FALSE} else{
						if (complement(study1.alleleA[i], study2.alleleB[i])) {flip[i] <- TRUE} else flip[i] <- NA
					}	
				}
			}
			if (complement(study2.alleleA[i], study2.alleleB[i])) strand.ambiguous[i] <- TRUE else strand.ambiguous[i] <- FALSE
		}	
	}
	
	if (is.null(study2.alleleB) & !is.null(study1.alleleB)){
		message("Allele matching based on partial information, other alleles are not provided for study2")
		for (i in 1:length(snpID)){
			if (study2.alleleA[i] == study1.alleleA[i]) {flip[i] <- FALSE} else{
				if (study2.alleleA[i] == study1.alleleB[i]) {flip[i] <- TRUE} else{
					if (complement(study2.alleleA[i], study1.alleleA[i])) {flip[i] <- FALSE} else{
						if (complement(study2.alleleA[i], study1.alleleB[i])) {flip[i] <- TRUE} else flip[i] <- NA
					}	
				}
			}
		if (complement(study1.alleleA[i], study1.alleleB[i])) strand.ambiguous[i] <- TRUE else strand.ambiguous[i] <- FALSE	
		}	
	}
	
	if (!is.null(study1.alleleB) & !is.null(study2.alleleB)){
		for (i in 1:length(snpID)){
			if (study2.alleleA[i] == study1.alleleA[i] & study2.alleleB[i] == study1.alleleB[i]) {flip[i] <- FALSE} else{
				if (study2.alleleA[i] == study1.alleleB[i] & study2.alleleB[i] == study1.alleleA[i]) {flip[i] <- TRUE} else{
					if (complement(study2.alleleA[i], study1.alleleA[i]) & complement(study2.alleleB[i], study1.alleleB[i])) {flip[i] <- FALSE} else{
						if (complement(study2.alleleA[i], study1.alleleB[i]) & complement(study2.alleleB[i], study1.alleleA[i])) {flip[i] <- TRUE} else flip[i] <- NA
					}	
				}
			}
		if (complement(study1.alleleA[i], study1.alleleB[i])) strand.ambiguous[i] <- TRUE else strand.ambiguous[i] <- FALSE	
		}	
		
	}
	
	inds.flip <- which(flip)
	if (length(inds.flip) > 0){
		study2.direction[inds.flip] <- -study2.direction[inds.flip]
	}
	
	
	if (is.null(strand.ambiguous)) return( data.frame(snpID, study1.direction, study2.direction, study1.alleleA)) else{
		return( data.frame(snpID, study1.direction, study2.direction, study1.alleleA, strand.ambiguous))
	}
}
