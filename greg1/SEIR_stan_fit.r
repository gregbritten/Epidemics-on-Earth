library(rstan)
library(lubridate)
library(lmodel2)
options(mc.cores = parallel::detectCores())
Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')

setwd('d:/dropbox/working/covid19/synthetic_observations/' )

#####################################################################################
## ESTIMATE PARAMETERS ##############################################################
#####################################################################################
##--COMPILE STAN CODE--###############
mod_betaAR1 <- stan_model('stan_SEIR_discrete_beta_AR1.stan')

mod2 <- stan_model('stan_SEIR_discrete_fixed.stan')

modrw <- stan_model('stan_SEIR_discrete_rw.stan')

#######################################
## READ DATA ##########################
#######################################
dat <- read.csv('d:/dropbox/working/covid19/urop/github/greg1/northeastern_uni.csv')
dat$time <- decimal_date(mdy(dat$date))
dat <- dat[rev(1:nrow(dat)),]
mat <- read.csv('d:/google/working/covid/massachusetts-history.csv')
mat$positiveIncrease[mat$positiveIncrease<0] <- 0
mat$time <- decimal_date(mdy(mat$date))
mat <- mat,[rev(1:nrow(mat)),]
mat <- mat[mat$time>2020.2,]

#####################################
## CASES & POSITIVITY ###############
#####################################
pdf('d:/dropbox/working/covid19/urop/northeastern_mass_cases_pos.pdf',height=5,width=6)
xlims=c(2020.2,2021)
par(mfrow=c(2,1),mar=c(2,4,2,4))
plot(dat$time,dat$postest,type='l',xlab='',ylab='',xlim=xlims); 
	#mtext(side=1,'Days',line=2.5)
	mtext(side=2,'Positive Cases',line=2.5)
par(new=TRUE)
plot(dat$time,dat$postest/(dat$postest + dat$negtest),type='l',col='blue',xlab='',ylab='',xaxt='n',yaxt='n',xlim=xlims)
axis(side=4)
	mtext(side=4,'Positive Rate',line=2.5,col='blue')
	mtext(adj=0.1,'Northeastern')

plot(mat$time,mat$positiveIncrease,type='l',xlab='',ylab='',xlim=xlims); 
	mtext(side=2,'Positive Cases',line=2.5)
par(new=TRUE)
plot(mat$time,mat$positiveIncrease/(mat$positiveIncrease + mat$negativeIncrease),type='l',col='blue',xlab='',ylab='',xaxt='n',yaxt='n',xlim=xlims,ylim=c(0,0.4))
axis(side=4)
	mtext(side=4,'Positive Rate',line=2.5,col='blue')
	mtext(adj=0.1,'Massachusetts')
dev.off()

###############################
## Ntest ######################
###############################

pdf('d:/dropbox/working/covid19/urop/northeastern_mass_ntest.pdf',height=5,width=6)
xlims=c(2020.2,2021)
par(mfrow=c(2,1),mar=c(2,4,2,4))
plot(dat$time,dat$postest+dat$negtest,type='l',xlab='',ylab='',xlim=xlims); 
	mtext(side=2,'Ntest',line=2.5)
	mtext(adj=0.1,'Northeastern')

plot(mat$time,mat$positiveIncrease+mat$negativeIncrease,type='l',xlab='',ylab='',xlim=xlims); 
	mtext(side=2,'Ntest',line=2.5)
	mtext(adj=0.1,'Massachusetts')
dev.off()


###########################################
## EFFECTIVE POPULATION SIZE ##############
###########################################
pdf('d:/dropbox/working/covid19/urop/northeastern_mass_effectivepop.pdf',height=4,width=10)
par(mfrow=c(1,2),mar=c(2,2,2,2),oma=c(2,2,2,2))
NeffNEmean <- mean(dat$postest/(dat$postest/(dat$postest + dat$negtest)),na.rm=TRUE)
NeffNE     <- dat$postest/(dat$postest/(dat$postest + dat$negtest))
plot(dat$time,dat$postest/(dat$postest/(dat$postest + dat$negtest)),xlim=xlims,ylab='')
	mtext(adj=0.1,'Northeastern')
abline(h=NeffNEmean,lty=2)
mtext(side=2,line=2.5,'Npop')
#plot(dat$postest/(dat$postest + dat$negtest)*Neff,dat$postest)
#abline(0,1)

posrateMA  <- mat$positiveIncrease/(mat$positiveIncrease + mat$negativeIncrease)
#posrateMA[posrateMA>0.35] <- 0.1
NeffMAmean <- mean(mat$positiveIncrease/posrateMA,na.rm=TRUE)
NeffMA     <- mat$positiveIncrease/posrateMA
plot(mat$time,NeffMA,xlim=xlims,ylab='')
abline(h=NeffMAmean,lty=2)
	mtext(adj=0.1,'Massachusetts')
mtext(side=2,line=2.5,'Npop')
dev.off()


yNE  <- dat$postest 
popNE <- NeffNEmean
yNEp <- dat$postest/(dat$postest+dat$negtest)*popNE

yMA <- mat$positiveIncrease
yMA[yMA>3600] <- 1000 
popMA <- 1E7
yMAp <- posrateMA*popMA
yMAp[yMAp>3.6E6] <- 1E6


###############################################################
## GROWTH RATES ###############################################
###############################################################
winsiz <- 14
NEroll  <- rollapply(yNE,width=winsiz,FUN=mean,align='right')
NErollp <- rollapply(yNEp,width=winsiz,FUN=mean,align='right')
MAroll <- rollapply(yMA,width=winsiz,FUN=mean,align='right')
MArollp <- rollapply(yMAp,width=winsiz,FUN=mean,align='right')

pdf('d:/dropbox/working/covid19/urop/MA_time.pdf',height=5,width=7)
par(mfrow=c(2,2),mar=c(2,2,2,2),oma=c(2,2,2,2))
plot(MAroll,type='l')
plot(MArollp,type='l')
NMA <- length(MAroll)
plot(log(MAroll[2:NMA]/MAroll[1:(NMA-1)]),ylim=c(-0.2,0.2),type='l')
abline(h=0)
plot(log(MArollp[2:NMA]/MArollp[1:(NMA-1)]),ylim=c(-0.2,0.2),type='l')
abline(h=0)
dev.off()

pdf('d:/dropbox/working/covid19/urop/NE_time.pdf',height=5,width=7)
par(mfrow=c(2,2),mar=c(2,2,2,2),oma=c(2,2,2,2))
plot(NEroll,type='l')
plot(NErollp,type='l')
NMA <- length(NEroll)
plot(log(NEroll[2:NMA]/NEroll[1:(NMA-1)]),ylim=c(-0.4,0.4),type='l')
abline(h=0)
plot(log(NErollp[2:NMA]/NErollp[1:(NMA-1)]),ylim=c(-0.4,0.4),type='l')
abline(h=0)
dev.off()
#################################
## SCATTER ######################
#################################
pdf('d:/dropbox/working/covid19/urop/regressions.pdf',height=4,width=8)
par(mfrow=c(1,2))

x <- log(MAroll[2:NMA]/MAroll[1:(NMA-1)])
y <- log(MArollp[2:NMA]/MArollp[1:(NMA-1)])

plot(x,y,ylim=c(-0.4,0.4),xlim=c(-0.4,0.4),pch=19,cex=0.6,xlab='',ylab='')
mtext(side=1,line=2.5,'Observed Growth Rates')
mtext(side=2,line=2.5,'Growth Rate from Positivity')
abline(0,1)
#abline(lm(y~x),lty=2)
lines(lmodel2(y~x),col=c('blue',NA))

x <- log(NEroll[2:NMA]/NEroll[1:(NMA-1)])
y <- log(NErollp[2:NMA]/NErollp[1:(NMA-1)])
plot(x,y,ylim=c(-0.4,0.4),xlim=c(-0.4,0.4),pch=19,cex=0.6,xlab='',ylab='')
mtext(side=1,line=2.5,'Observed Growth Rates')
mtext(side=2,line=2.5,'Growth Rate from Positivity')
abline(0,1)
#abline(lm(y~x),lty=2)
abline(lmodel2(y~x),lty=2)
lines(lmodel2(y~x),col=c('blue',NA))

dev.off()

# N <- nrow(dat)
# par(mfrow=c(2,1),mar=c(2,4,2,4))
# plot(dat$time[-c(1:7)],diff(NEroll),type='l',xlab='',ylab='',xlim=xlims); 
# abline(h=0)
# plot(dat$time[-1],log(dat$postest[2:N]/dat$postest[1:(N-1)]),type='l',xlab='',ylab='',xlim=xlims); 
# abline(h=0)
	# mtext(side=1,'Days',line=2.5)
	# mtext(side=2,'Observed Positive Cases',line=2.5)
# par(new=TRUE)
# plot(dat$time,dat$postest/(dat$postest + dat$negtest),type='l',col='blue',xlab='',ylab='',xaxt='n',yaxt='n',xlim=xlims)
# axis(side=4)



#plot(posrateMA*NeffMAmean,mat$positiveIncrease)
#abline(0,1)

#################################
## FIT TO NY DATA ###############
#################################
E0 <- 0.05
S0 <- 1.0 - E0
R0 <- 0.0
I0 <- 0.0

x0 <- c(S0,E0,I0,R0) #initial conditions

yNE  <- dat$postest 
popNE <- NeffNEmean
yNEp <- dat$postest/(dat$postest+dat$negtest)*popNE
plot(dat$time,yNE,xlim=xlims)
plot(dat$time,yNEp,xlim=xlims)

yMA <- mat$positiveIncrease
yMA[yMA>3600] <- 1000 
popMA <- 1E7
yMAp <- posrateMA*popMA
yMAp[yMAp>3.6E6] <- 1E6
plot(mat$time,yMA,xlim=xlims)
plot(mat$time,yMAp,xlim=xlims)


dataNE <- list(POP=popNE,
			 y=yNE,
			 N_obs=length(yNE),
			 t_obs=1:length(yNE),
			 x0 = x0)
dataNEp <- list(POP=popNE,
			 y=round(yNEp),
			 N_obs=length(yNEp),
			 t_obs=1:length(yNEp),
			 x0 = x0)
dataMA <- list(POP=1E5,
			 y=yMA,
			 N_obs=length(yMA),
			 t_obs=1:length(yMA),
			 x0 = x0)
dataMAp <- list(POP=popMA,
			 y=yMAp,
			 N_obs=length(yMAp),
			 t_obs=1:length(yMAp),
			 x0 = x0)

mcmcNE <- sampling(mod_betaAR1, data=dataNE ,open_progress=TRUE)
mcmcNE <- sampling(mod2, data=dataNE ,open_progress=TRUE)
mcmcNErw <- sampling(modrw, data=dataNE ,open_progress=TRUE)

postNErw <- extract(mcmcNErw)
postNE <- extract(mcmcNE)

#mcmcNEp <- sampling(mod_betaAR1, data=dataNEp ,open_progress=TRUE)
mcmcNEp <- sampling(mod2, data=dataNEp ,open_progress=TRUE)
postNEp <- extract(mcmcNEp)

mcmcMA <- sampling(mod_betaAR1, data=dataMA ,open_progress=TRUE)
mcmcMArw <- sampling(modrw, data=dataMA ,open_progress=TRUE)
#mcmcMA <- sampling(mod2, data=dataMA ,open_progress=TRUE)
postMA <- extract(mcmcMArw)

mcmcMAp <- sampling(mod_betaAR1, data=dataMAp ,open_progress=TRUE)
#mcmcMAp <- sampling(mod2, data=dataMAp ,open_progress=TRUE)
postMAp <- extract(mcmcMAp)

par(mfrow=c(2,1),mar=c(2,4,2,4))
plot(colMeans(postNErw$beta),type='l')
lines(colMeans(postNEp$beta))
plot(colMeans(postNE$beta),type='l')
lines(colMeans(postNEp$beta))

plot(colMeans(postMA$beta),type='l')
lines(colMeans(postMAp$beta))

######################################################i#############
## OPTIMIZING #####################################################
###################################################################
optNE <- optimizing(mod_betaAR1,data=data,algorithm="Newton",hessian=TRUE)

optMA <- optimizing(mod2,data=dataMA,algorithm="Newton",hessian=TRUE,init='0')

optMA <- optimizing(mod_betaAR1,data=dataMA,algorithm="Newton",hessian=TRUE,init='0')


plot(optMA$par[1:113])
abline(h=1)
plot(opt$par[1:38])
lines(opt$par[1:38] + 2*sqrt(diag(solve(-opt$hessian)))[1:38])
lines(opt$par[1:38] - 2*sqrt(diag(solve(-opt$hessian)))[1:38])

opt <- optimizing(mod_betaAR1,data=data,algorithm="Newton",hessian=TRUE)

mcmc <- sampling(mod_betaAR1,data=data,open_progress=TRUE)
