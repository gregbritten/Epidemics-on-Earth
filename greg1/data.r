library(colorRamps)

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
par(mfrow=c(8,7),mar=c(0,2,0,0),oma=c(4,2,1,1) ,cex.axis=0.7)
DATA <- list()
k <- 1
ex <- c(4,13,28,43,51)

#for(i in c(1:3,5:56)){
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
	doy  <- rev(st$doy)
	dHOS <- rev(st$hospitalizedIncrease); dHOS[1] <- 0

	dDEAD <- rev(st$deathIncrease)

	N_obs <- length(dPOS)
	
	d <- dPOS
	#d <- round(approx(POP*(dPOS/(dPOS+dNEG)),n=length(dPOS))$y)
	d[d<0] <- 0
	if(length(POP)>0){
	DATA[[states[i]]] <- list(POP=POP,
				      y=d,
					  dPOS=dPOS,
					  dNEG=dNEG,
				      N_obs=N_obs, 
				      x0 = x0,
				      doy=doy,
					  state=states[i],
					  Ntest=dPOS+dNEG,
					  dDEAD=dDEAD,
					  E0=E0,states=states[i])
	k <- k+1				  
	}

	plot(doy,dPOS,type='l',xlim=c(65,365),xaxt='n',bty='n')
	mtext(states[i],line=-1.5,adj=0.05)
	if(i>49) axis(side=1)
	lw <- loess(dPOS  ~ doy,span=0.6)
	DATA[[states[i]]]$loess_case <- lw$fitted
	lines(doy,lw$fitted,lwd=1.5,col='red')
}
mtext(side=2,outer=TRUE,'Daily Cases')
mtext(side=1,outer=TRUE,'Day of Year',line=2.5)


##--PLOT TESTS--###########
par(mfrow=c(8,7),mar=c(0,2,0,0),oma=c(4,2,1,1) ,cex.axis=0.7)
for(i in 1:length(DATA)){
    ntest <- DATA[[i]]$dPOS + DATA[[i]]$dNEG
	ntest[ntest<0] <- 0
	plot(DATA[[i]]$doy,ntest,xaxt='n',type='l',bty='n',xlim=c(63,360))
	lw <- loess(DATA[[i]]$dPOS + DATA[[i]]$dNEG  ~ DATA[[i]]$doy,span=0.6)
	lines(DATA[[i]]$doy,lw$fitted,lwd=1.5,col='red')
	DATA[[i]]$loess_test <- lw$fitted
	mtext(DATA[[i]]$state,line=-1.5,adj=0.05)
	if(i>44) axis(side=1)
}
mtext(side=1,outer=TRUE,'Day of Year',line=1)
mtext(side=2,outer=TRUE,'Number of Tests',line=0.5)


##--PLOT TESTS--###########
par(mfrow=c(8,7),mar=c(0,2,0,0),oma=c(4,2,1,1) ,cex.axis=0.7)
for(i in 1:length(DATA)){
    ntest <- DATA[[i]]$dPOS + DATA[[i]]$dNEG
	ntest[ntest<0] <- 0
	posrate <- DATA[[i]]$dPOS/(DATA[[i]]$dPOS+DATA[[i]]$dNEG)
	posrate[posrate<0] <- 0
	plot(ntest,posrate,xaxt='n',bty='n',pch=19,cex=0.4,ylim=c(0,1))

#	lw <- loess(DATA[[i]]$dPOS + DATA[[i]]$dNEG  ~ DATA[[i]]$doy,span=0.6)
#	lines(DATA[[i]]$doy,lw$fitted,lwd=1.5,col='red')
#	DATA[[i]]$loess_test <- lw$fitted
	mtext(DATA[[i]]$state,line=-1.5,adj=0.1,cex=0.7)
	if(i>44) axis(side=1)
}
mtext(side=1,outer=TRUE,'Number of Tests',line=1)
mtext(side=2,outer=TRUE,'Positivity Rate',line=0.5)



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



par(mfrow=c(8,7),mar=c(1,2,1,1),oma=c(3,3,1,1),cex.axis=0.7)
for(i in 1:length(DATA)){
	x <- DATA[[i]]$loess_test
	y <- DATA[[i]]$loess_case
	rng <- range(c(x,y))
	plot(x,y,type='l',bty='n')#,xlim=rng,ylim=rng)
	#abline(0,1)
}

