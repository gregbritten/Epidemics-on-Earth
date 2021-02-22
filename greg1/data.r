library(colorRamps)
library(mgcv)
library(viridis)
setwd('d:/dropbox/working/covid19/urop/github/greg1/' )

##############################################################
## READ DATA #################################################
##############################################################
dat    <- read.csv('daily.csv',stringsAsFactor=FALSE)
abbrev <- read.csv('state_abbrev.csv',stringsAsFactor=FALSE)
pops   <- read.csv('co-est2019-alldata.csv',stringsAsFactor=FALSE)

states <- unique(dat$state)

#############################
## PLOT #####################
#############################
DATA <- list()
ex <- c(4,13,28,43,51)


for(i in (1:length(states))[-ex]){
	states[i]
	name <- abbrev$name[abbrev$abbrev==states[i]]	
	POP <- pops$POPESTIMATE2019[pops$STNAME==name & pops$COUNTY==0]
	
	st <- dat[dat$state==states[i],]
	st$doy <- as.numeric(strftime(paste(substr(as.character(st$date),1,4), 
										 substr(as.character(st$date),5,6), 
										 substr(as.character(st$date),7,8), sep='-'), format = "%j")) 
	dPOS <- rev(st$positiveIncrease);     dPOS[1] <- 0
		dPOS[dPOS<0] <- 0
	dNEG <- rev(st$negativeIncrease);     dNEG[1] <- 0
		dNEG[dNEG<0] <- 0
	doy  <- rev(st$doy)
	dHOS <- rev(st$hospitalizedIncrease); dHOS[1] <- 0
	dDEAD <- rev(st$deathIncrease)     

	Ntest <- dPOS + dNEG
	N_obs <- length(dPOS)
	
	d <- dPOS
	d[d<0] <- 0

	if(length(POP)>0){
	DATA[[states[i]]] <- list(POP=POP,
					  y=d,
					  dPOS=dPOS,
					  dNEG=dNEG,
					  N_obs=N_obs, 
					  doy=doy,
					  state=states[i],
					  Ntest=Ntest,
					  dDEAD=dDEAD,
					  dHOS=dHOS,
					  states=states[i])				  
	}


	lw                           <- loess(dPOS  ~ doy,span=0.3)
	DATA[[states[i]]]$loess_case <- lw$fitted
	
	lw                           <- loess(Ntest  ~ doy,span=0.3)
	DATA[[states[i]]]$loess_test <- lw$fitted

	lw                           <- loess(dNEG ~ doy,span=0.3)
	DATA[[states[i]]]$loess_neg  <- lw$fitted
	
	lw                           <- loess(dDEAD ~ doy,span=0.3)
	DATA[[states[i]]]$loess_dead <- lw$fitted
	
	lw                           <- loess(dHOS ~ doy,span=0.3)
	DATA[[states[i]]]$loess_hos  <- lw$fitted		
}

pdf('d:/dropbox/working/covid19/urop/daily_cases_12.17.2020.pdf',height=9,width=14)
par(mfrow=c(8,7),mar=c(0,2,0,0),oma=c(4,2,1,1) ,cex.axis=0.7)
	for(i in 1:length(DATA)){
		dPOS <- DATA[[i]]$dPOS
		doy  <- DATA[[i]]$doy
		lw   <- DATA[[i]]$loess_case
		plot(doy,dPOS,type='l',xlim=c(65,365),xaxt='n',bty='n')
		lines(doy,lw,col='red')
		mtext(states[i],line=-1.5,adj=0.05)
		if(i>48) axis(side=1)
	mtext(side=2,outer=TRUE,'Daily Cases')
	mtext(side=1,outer=TRUE,'Day of Year',line=2.5)
}

###########################################################
## GROWTH RATES OVER TIME #################################
###########################################################
pdf('d:/dropbox/working/covid19/urop/daily_growthrates_12.17.2020.pdf',height=9,width=14)
par(mfrow=c(8,7),mar=c(0,2,0,0),oma=c(4,2,1,1) ,cex.axis=0.7)
	for(i in 1:length(DATA)){
		tmp <- DATA[[i]]$loess_case
		n <- length(tmp)
		plot(DATA[[i]]$doy[-1],log(tmp[2:n]/tmp[1:(n-1)]),xaxt='n',type='l',bty='n',xlim=c(63,360),ylim=c(-0.1,0.1))
		mtext(DATA[[i]]$state,line=-4.5,adj=0.05)
		if(i>44) axis(side=1)
	}
	mtext(side=1,outer=TRUE,'Day of Year',line=1)
	mtext(side=2,outer=TRUE,'Growth Rate in Cases',line=0.5)
dev.off()


pdf('d:/dropbox/working/covid19/urop/daily_testgrowthrates_12.17.2020.pdf',height=9,width=14)
par(mfrow=c(8,7),mar=c(0,2,0,0),oma=c(4,2,1,1) ,cex.axis=0.7)
	for(i in 1:length(DATA)){
		tmp <- DATA[[i]]$loess_test
		n <- length(tmp)
		plot(DATA[[i]]$doy[-1],log(tmp[2:n]/tmp[1:(n-1)]),xaxt='n',type='l',bty='n',xlim=c(63,360),ylim=c(-0.1,0.1))
		mtext(DATA[[i]]$state,line=-4.5,adj=0.05)
		if(i>44) axis(side=1)
	}
	mtext(side=1,outer=TRUE,'Day of Year',line=1)
	mtext(side=2,outer=TRUE,'Growth Rate in Tests',line=0.5)
dev.off()


###########################################################
## GROWTH RATES AGAINST GROWTH RATES ######################
###########################################################
pdf('d:/dropbox/working/covid19/urop/growthrates_VS_12.17.2020.pdf',height=9,width=14)
par(mfrow=c(8,7),mar=c(0,2,0,0),oma=c(4,2,1,1) ,cex.axis=0.7)
	for(i in 1:length(DATA)){
		tmp <- DATA[[i]]$loess_test
		tmp2 <- DATA[[i]]$loess_case		
		n <- length(tmp)
		cols <- viridis(n)
		gtmp <- log(tmp[2:n]/tmp[1:(n-1)])
		gtmp2 <- log(tmp2[2:n]/tmp2[1:(n-1)])
		plot(gtmp,gtmp2,xaxt='n',bty='n',ylim=c(-0.1,0.1),xlim=c(-0.1,0.1),col=cols,pch=19)
		mtext(DATA[[i]]$state,line=-1.5,adj=0.05)
		if(i>44) axis(side=1)
			abline(h=0,lty=2)
	abline(v=0,lty=2)
	abline(0,1)
	}
	mtext(side=2,outer=TRUE,'Growth Rate in Cases',line=0.5)
	mtext(side=1,outer=TRUE,'Growth Rate in Tests',line=0.5)
dev.off()


pdf('d:/dropbox/working/covid19/urop/growthrates_VS_hospitalizations_01.13.2020.pdf',height=9,width=14)
par(mfrow=c(8,7),mar=c(0,2,0,0),oma=c(4,2,1,1) ,cex.axis=0.7)
	for(i in 1:length(DATA)){
		tmp <- DATA[[i]]$loess_test
		tmp2 <- DATA[[i]]$loess_hos		
		n <- length(tmp)
		cols <- viridis(n)
		gtmp <- log(tmp[2:n]/tmp[1:(n-1)])
		gtmp2 <- log(tmp2[2:n]/tmp2[1:(n-1)])
		plot(gtmp,gtmp2,xaxt='n',bty='n',ylim=c(-0.1,0.1),xlim=c(-0.1,0.1),col=cols,pch=19)
		mtext(DATA[[i]]$state,line=-1.5,adj=0.05)
		if(i>44) axis(side=1)
			abline(h=0,lty=2)
	abline(v=0,lty=2)
	abline(0,1)
	}
	mtext(side=1,outer=TRUE,'Growth Rate in Tests',line=0.5)
	mtext(side=2,outer=TRUE,'Growth Rate in Hospitalizations',line=0.5)
dev.off()



pdf('d:/dropbox/working/covid19/urop/growthrates_VS_deaths_01.13.2020.pdf',height=9,width=14)
par(mfrow=c(8,7),mar=c(0,2,0,0),oma=c(4,2,1,1) ,cex.axis=0.7)
	for(i in 1:length(DATA)){
		tmp <- DATA[[i]]$loess_test
		tmp2 <- DATA[[i]]$loess_dead		
		n <- length(tmp)
		cols <- viridis(n)
		gtmp <- log(tmp[2:n]/tmp[1:(n-1)])
		gtmp2 <- log(tmp2[2:n]/tmp2[1:(n-1)])
		plot(gtmp,gtmp2,xaxt='n',bty='n',ylim=c(-0.1,0.1),xlim=c(-0.1,0.1),col=cols,pch=19)
		mtext(DATA[[i]]$state,line=-1.5,adj=0.05)
		if(i>44) axis(side=1)
			abline(h=0,lty=2)
	abline(v=0,lty=2)
	abline(0,1)
	}
	mtext(side=1,outer=TRUE,'Growth Rate in Tests',line=0.5)
	mtext(side=2,outer=TRUE,'Growth Rate in Deaths',line=0.5)
dev.off()


pdf('d:/dropbox/working/covid19/urop/cases_VS_deaths_01.13.2020.pdf',height=9,width=14)
par(mfrow=c(8,7),mar=c(0,2,0,0),oma=c(4,2,1,1) ,cex.axis=0.7)
	for(i in 1:length(DATA)){
		tmp <- DATA[[i]]$loess_case
		tmp2 <- DATA[[i]]$loess_dead		
		n <- length(tmp)
		cols <- viridis(n)
		gtmp <- log(tmp[2:n]/tmp[1:(n-1)])
		gtmp2 <- log(tmp2[2:n]/tmp2[1:(n-1)])
		plot(gtmp,gtmp2,xaxt='n',bty='n',ylim=c(-0.1,0.1),xlim=c(-0.1,0.1),col=cols,pch=19)
		mtext(DATA[[i]]$state,line=-1.5,adj=0.05)
		if(i>44) axis(side=1)
			abline(h=0,lty=2)
	abline(v=0,lty=2)
	abline(0,1)
	}
	mtext(side=1,outer=TRUE,'Growth Rate in Cases',line=0.5)
	mtext(side=2,outer=TRUE,'Growth Rate in Deaths',line=0.5)
dev.off()



####################################################

##--PLOT TESTS--###########
pdf('d:/dropbox/working/covid19/urop/daily_tests_12.17.2020.pdf',height=9,width=14)
par(mfrow=c(8,7),mar=c(0,2,0,0),oma=c(4,2,1,1) ,cex.axis=0.7)
	for(i in 1:length(DATA)){
		ntest <- DATA[[i]]$Ntest
		plot(DATA[[i]]$doy,ntest,xaxt='n',type='l',bty='n',xlim=c(63,360))
		mtext(DATA[[i]]$state,line=-1.5,adj=0.05)
		if(i>44) axis(side=1)
	}
	mtext(side=1,outer=TRUE,'Day of Year',line=1)
	mtext(side=2,outer=TRUE,'Number of Tests',line=0.5)
dev.off()




par(mfrow=c(8,7),mar=c(0,2,0,0),oma=c(4,2,1,1) ,cex.axis=0.7)
	for(i in 1:length(DATA)){
		ntest <- DATA[[i]]$dPOS + DATA[[i]]$dNEG
		ntest[ntest<0] <- 0
#		n <- length(ntest)
		#plot(DATA[[i]]$doy,ntest,xaxt='n',type='l',bty='n',xlim=c(63,360))
		#lw <- loess(DATA[[i]]$dPOS + DATA[[i]]$dNEG  ~ DATA[[i]]$doy,span=0.6)
		ntestr <- rollapply(ntest,width=7,FUN=mean,fill=0,align='r')
		n <- length(ntestr)
		tg <- log(ntest[2:n]/ntest[1:(n-1)])

#		tg <- log(ntestr[2:n]/ntestr[1:(n-1)])
		troll <- rollapply(tg, width=7, FUN=sd, fill=0, align="r")
		#lines(DATA[[i]]$doy,lw$fitted,lwd=1.5,col='red')
		#DATA[[i]]$loess_test <- lw$fitted
		plot(ntestr[-1]/troll)
		
		mtext(DATA[[i]]$state,line=,-1.5,adj=0.05)
		if(i>44) axis(side=1)
	}
	mtext(side=1,outer=TRUE,'Day of Year',line=1)
	mtext(side=2,outer=TRUE,'Number of Tests',line=0.5)


library(mgcv)

par(mfrow=c(8,7),mar=c(0,2,0,0),oma=c(4,2,1,1) ,cex.axis=0.7)
	for(i in 1:length(DATA)){
		#ntest <- DATA[[i]]$dPOS + DATA[[i]]$dNEG
		#ntest[ntest<0] <- 0
#		n <- length(ntest)
		#plot(DATA[[i]]$doy,ntest,xaxt='n',type='l',bty='n',xlim=c(63,360))
		#lw <- loess(DATA[[i]]$dPOS + DATA[[i]]$dNEG  ~ DATA[[i]]$doy,span=0.6)
		
		Iobs <- rollapply(DATA[[i]]$dPOS,width=7,FUN=mean,fill=0,align='r',trim=0.1)
		#Iobs <- DATA[[i]]$dPOS
		n <- length(Iobs)
		#ts <- 1:n
		#fit <- gam(Iobs ~ s(ts))
		#Ism <- as.numeric(predict(fit,newdata=list(ts=ts)))

		Ig   <- log(Iobs[2:n]/Iobs[1:(n-1)])
#		Ig   <- log(Ism[2:n]/Ism[1:(n-1)])
		
		Igsd <- rollapply(Ig,width=7,FUN=sd,fill=0,align='r')
		
		plot(Ig,type='l')
		
		#ntestr <- rollapply(ntest,width=7,FUN=mean,fill=0,align='r')
		#n <- length(ntestr)
		#tg <- log(ntest[2:n]/ntest[1:(n-1)])

#		tg <- log(ntestr[2:n]/ntestr[1:(n-1)])
		#troll <- rollapply(tg, width=7, FUN=sd, fill=0, align="r")
		#lines(DATA[[i]]$doy,lw$fitted,lwd=1.5,col='red')
		#DATA[[i]]$loess_test <- lw$fitted
		#plot(ntestr[-1]/troll)
		
		mtext(DATA[[i]]$state,line=-1.5,adj=0.05)
		if(i>44) axis(side=1)
	}
	mtext(side=1,outer=TRUE,'Day of Year',line=1)
	mtext(side=2,outer=TRUE,'Number of Tests',line=0.5)




##--PLOT TESTS--###########
pdf('d:/dropbox/working/covid19/urop/posrate_ntest_12.17.2020.pdf',height=9,width=11)
par(mfrow=c(8,7),mar=c(0,2,1,0),oma=c(4,2,1,1) ,cex.axis=0.7)
	for(i in 1:length(DATA)){
		ntest <- DATA[[i]]$dPOS + DATA[[i]]$dNEG
		ntest[ntest<0] <- 0
		posrate <- DATA[[i]]$dPOS/(DATA[[i]]$dPOS+DATA[[i]]$dNEG)
		posrate[posrate<0] <- 0
		posrate[!is.finite(posrate)] <- 0
		plot(ntest,posrate,xaxt='n',bty='n',pch=19,cex=0.4,ylim=c(0,1))

	#	lw <- loess(DATA[[i]]$dPOS + DATA[[i]]$dNEG  ~ DATA[[i]]$doy,span=0.6)
	#	lines(DATA[[i]]$doy,lw$fitted,lwd=1.5,col='red')
	#	DATA[[i]]$loess_test <- lw$fitted
		mtext(DATA[[i]]$state,line=-1.5,adj=0.1,cex=0.7)
		if(i>44) axis(side=1)
		#lm <- lm(posrate  ~ ntest)
		#abline(lm,col='red')
		fit <- gam(posrate ~ s(ntest))
		lines(seq(1,max(ntest)),predict(fit,newdata=list(ntest=seq(1,max(ntest)))),col='red')
	}
	mtext(side=1,outer=TRUE,'Number of Tests',line=1)
	mtext(side=2,outer=TRUE,'Positivity Rate',line=0.5)
dev.off()



##--PLOT TESTS--###########
pdf('d:/dropbox/working/covid19/urop/posrate_ntest_12.17.2020.pdf',height=9,width=11)
par(mfrow=c(8,7),mar=c(0,2,1,0),oma=c(4,2,1,1) ,cex.axis=0.7)
	for(i in 1:length(DATA)){
		ntest <- DATA[[i]]$loess_test
		ntest[ntest<0] <- 0
		n <- length(ntest)
		tmp <- log(ntest[2:n]/ntest[1:(n-1)])
		posrate <- DATA[[i]]$dPOS/(DATA[[i]]$dPOS+DATA[[i]]$dNEG)
		posrate[posrate<0] <- 0
		posrate[!is.finite(posrate)] <- 0
		plot(tmp,posrate[-1],xaxt='n',bty='n',pch=19,cex=0.4,ylim=c(0,1))

	#	lw <- loess(DATA[[i]]$dPOS + DATA[[i]]$dNEG  ~ DATA[[i]]$doy,span=0.6)
	#	lines(DATA[[i]]$doy,lw$fitted,lwd=1.5,col='red')
	#	DATA[[i]]$loess_test <- lw$fitted
		mtext(DATA[[i]]$state,line=-1.5,adj=0.1,cex=0.7)
		if(i>44) axis(side=1)
		#lm <- lm(posrate  ~ ntest)
		#abline(lm,col='red')
		#fit <- gam(posrate ~ s(ntest))
		#lines(seq(1,max(ntest)),predict(fit,newdata=list(ntest=seq(1,max(ntest)))),col='red')
	}
	mtext(side=1,outer=TRUE,'Testing growth rate',line=1)
	mtext(side=2,outer=TRUE,'Positivity Rate',line=0.5)
dev.off()



##--PLOT TESTS--###########
# par(mfrow=c(8,7),mar=c(0,2,0,0),oma=c(4,2,1,1) ,cex.axis=0.7)
# for(i in 1:length(DATA)){
    # ntest <- DATA[[i]]$dPOS + DATA[[i]]$dNEG
	# ntest[ntest<0] <- 0
	# plot(DATA[[i]]$doy[-1],log(ntest[2:length(ntest)]/ntest[1:(length(ntest)-1)]),xaxt='n',type='l',bty='n',xlim=c(63,360))
	# lw <- loess(DATA[[i]]$dPOS + DATA[[i]]$dNEG  ~ DATA[[i]]$doy,span=0.6)
	# lines(DATA[[i]]$doy,lw$fitted,lwd=1.5,col='red')
	# DATA[[i]]$loess_test <- lw$fitted
	# mtext(DATA[[i]]$state,line=-1.5,adj=0.05)
	# if(i>44) axis(side=1)
# }
# mtext(side=1,outer=TRUE,'Day of Year',line=1)
# mtext(side=2,outer=TRUE,'Number of Tests',line=0.5)


pdf('d:/dropbox/working/covid19/urop/test_case_smooth_12.17.2020.pdf',height=9,width=11)
par(mfrow=c(8,7),mar=c(1,2,1,1),oma=c(3,3,1,1),cex.axis=0.7)
for(i in 1:length(DATA)){
	x <- DATA[[i]]$loess_test
	y <- DATA[[i]]$loess_case
	rng <- range(c(x,y))
	plot(x,y,type='l',bty='n')#,xlim=rng,ylim=rng)
		mtext(DATA[[i]]$state,line=-1.5,adj=0.1,cex=0.7)
	#abline(0,1)
}
dev.off()




