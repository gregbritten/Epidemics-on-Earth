

NE <- read.csv('d:/dropbox/working/covid19/urop/github/greg1/northeastern.csv')

NEposrate <- NE$postest/(NE$postest + NE$negtest)

par(mfrow=c(2,3))
	y1 <- rev(NE$postest)
	y2 <- rev(NE$ntest)
	x <- 1:114
	lwy1 <- loess(y1 ~ x)$fitted
	g1   <- log(lwy1[2:114]/lwy1[1:113])

	lwy2 <- loess(y2 ~ x)$fitted
	g2   <- log(lwy2[2:114]/lwy2[1:113])

	plot(rev(NE$postest),type='l')
	
	lines(lwy1)
	
	plot(rev(NE$ntest),type='l')
	lines(lwy2)
	
		plot(g2,g1)
	abline(h=0,lty=2)
	abline(v=0,lty=2)
	abline(0,1)
	
	plot(g1,type='l')
	plot(g2,type='l')


	