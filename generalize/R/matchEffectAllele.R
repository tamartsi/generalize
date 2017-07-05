matchEffectAllele <-
function(snpID, study2.effect, study1.alleleA, study2.alleleA, study1.alleleB = NULL, study2.alleleB = NULL){
	stopifnot( length(snpID) == length(study2.effect), length(snpID) == length(study1.alleleA), length(snpID) == length(study2.alleleA))
	
	inds.na <- unique(c(which(is.na(snpID)), which(is.na(study2.effect)), which(is.na(study1.alleleA)), which(is.na(study2.alleleA))))
	if (length(inds.na) > 0){
		message("Some of the SNPs had missing data, removing them")
		snpID <- snpID[-inds.na]
		study2.effect <- study2.effect[-inds.na]
		study1.alleleA <- study1.alleleA[-inds.na]
		study2.alleleA <- study2.alleleA[-inds.na]
		if (!is.null(study1.alleleB)) study1.alleleB <- study1.alleleB[-inds.na]
		if (!is.null(study2.alleleB)) study2.alleleB <- study2.alleleB[-inds.na]
	}
	
	study1.alleleA <- toupper(study1.alleleA)
	study2.alleleA <- toupper(study2.alleleA)
	
	if (!is.null(study1.alleleB)) {
		stopifnot(length(snpID) == length(study1.alleleB))
		study1.alleleB <- toupper(study1.alleleB)
		if (sum(is.na(study1.alleleB)) > 0) {
			message(paste0(sum(is.na(study1.alleleB)), " alleleB from study1 are NAs, removing information for study1.alleleB"))
			study1.alleleB <- NULL
			
		}
		if (sum(study1.alleleA == study1.alleleB) > 0){
			message("alleleA and alleleB of study1 are the same for some SNP, discarding alleleB info for study1")
			study1.alleleB <- NULL
		}
		}
	if (!is.null(study2.alleleB)) {
		stopifnot(length(snpID) == length(study2.alleleB))	
		study2.alleleB <- toupper(study2.alleleB)
		if (sum(is.na(study2.alleleB)) > 0) {
			message(paste0(sum(is.na(study2.alleleB)), " alleleB from study1 are NAs, removing information for study1.alleleB"))
			study2.alleleB <- NULL	
		}
		if (sum(study2.alleleA == study2.alleleB) > 0){
			message("alleleA and alleleB of study2 are the same for some SNP, discarding alleleB info for study2")
			study2.alleleB <- NULL
		}
		}
		
		
	message("passed data entry checks, orienting the effects of study2 to the correct effect allele. Assuming same strand.")
	
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
		study2.effect[inds.flip] <- -study2.effect[inds.flip]
	}
	
	
	if (is.null(strand.ambiguous)) return( data.frame(snpID, study2.effect, study1.alleleA, flip)) else{
		return( data.frame(snpID, study2.effect, study1.alleleA, flip, strand.ambiguous))
	}
}
