prepareGenResFigure <-
function(snpID, study1.beta, study1.se, study2.beta, study2.se,generalized, gen.rvals, study1.n.test, output.file, study1.name = "Study1", study2.name = "Study2", make.title = FALSE, legend.position = c(.9, .9), CI.line.width = 1.7, x.ticks.size = 15, plot.width = 8, plot.height = 8, trait = NULL, plot.title = NULL, add.stars = NULL, alpha.level.follow = NULL){
	
	stopifnot(length(snpID) == length(study1.beta), length(study1.beta) == length(study1.se), length(study1.beta) == length(study2.beta), length(snpID) == length(study2.se))
	
	if (is.null(alpha.level.follow)) alpha.level.follow <- 0.05
	
	alpha.level.dscvr <- 0.05*length(snpID)/study1.n.test
	
	ind.rev <- which(sign(study1.beta) < 0)
	if (length(ind.rev) > 0){
		study1.beta[ind.rev] <- -study1.beta[ind.rev]
		study2.beta[ind.rev] <- -study2.beta[ind.rev]
	}
	

	colors <- c("#7570b3", "#bcbddc","#e6550d", "#fdae6b")
	names(colors) <- c(paste(study1.name, "(generalized)") , paste(study2.name), paste(study1.name, "(not generalized)"), paste(study2.name, " "))
	
		
	colScale <- scale_color_manual(values=colors, breaks=names(colors))
	dodge <- position_dodge(width=0.2)

	Z.follow <- qnorm(alpha.level.follow/2, lower.tail = F)
	Z.dscvr <- qnorm(alpha.level.dscvr/2, lower.tail = F)

	for.bar.plot <- data.frame(snpID = c(snpID, snpID), beta = c(study1.beta, study2.beta), se = c(study1.se, study2.se))
	for.bar.plot$study <- c(rep(study1.name, length(snpID)), rep(study2.name, length(snpID)))
	inds.generalized <- c(which(generalized == TRUE), which(generalized == TRUE) + length(snpID)) 
	inds.not.generalized <- setdiff(1:nrow(for.bar.plot), inds.generalized)
	
	if (length(inds.generalized) > 0){
		inds.gen.dscvr <- intersect(inds.generalized, grep(study1.name, for.bar.plot$study))
		for.bar.plot$study[inds.gen.dscvr] <- paste(for.bar.plot$study[inds.gen.dscvr], "(generalized)")	
	} 
	if (length(inds.not.generalized) > 0) {
		inds.not.gen.dscvr <- intersect(inds.not.generalized, grep(study1.name, for.bar.plot$study))
		inds.not.gen.follow <- intersect(inds.not.generalized, grep(study2.name, for.bar.plot$study))
		for.bar.plot$study[inds.not.gen.dscvr] <- paste(for.bar.plot$study[inds.not.gen.dscvr],"(not generalized)")
		for.bar.plot$study[inds.not.gen.follow] <- paste(for.bar.plot$study[inds.not.gen.follow]," ")
		}
	
		
	for.bar.plot$snpID <- factor(for.bar.plot$snpID, levels = for.bar.plot[order(gen.rvals, decreasing = F), "snpID"])
	
	for.bar.plot$Z <- c(rep(Z.dscvr, length(snpID)), rep(Z.follow, length(snpID))) 	
	

	
	if (is.null(plot.title)){
		if (make.title & study1.name != "Study1") plot.title <- paste0(trait, " SNPs from ", study1.name)
		if (make.title & study1.name == "Study1") plot.title <- paste0(trait, " SNPs")
		if (!make.title) plot.title <- paste0("") 
	}
	
	p1 <- ggplot(for.bar.plot, aes(x = snpID, y= beta, colour = study)) + ggtitle(plot.title) +
		geom_errorbar(aes(ymin=beta - se*Z, ymax = beta +se*Z), width = 0, position = dodge, lwd = CI.line.width) +
		geom_point(position = dodge, size = 2.5) +  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size= x.ticks.size))  + colScale  +  theme(legend.position = legend.position) +   theme(legend.title=element_blank()) + labs(x = "SNP ID", y = "Effect size")

 

	
	if (!is.null(add.stars)){
		for.bar.plot$beta.star <- for.bar.plot$beta
		for.bar.plot$beta.star[grep(study1.name, for.bar.plot$study)] <- NA
		for.bar.plot$beta.star[which(pchisq((for.bar.plot$beta/for.bar.plot$se)^2, df = 1, lower.tail = F) > add.stars)] <- NA
		
		p1 <- p1 + geom_point(data = for.bar.plot[which(!is.na(for.bar.plot$beta.star)),], position = position_dodge(width=0.3), size = 2.5, colour = "red", shape = 8) 
		}
	
	p1
	
	ggsave(filename = output.file, width = plot.width, height = plot.height)


	
}
