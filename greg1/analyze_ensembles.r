
library(zoo)
setwd('d:/dropbox/working/covid19/urop/github/greg1/')


bs <- c(0.10,0.15,0.20,0.25,0.30)
ns <- c(1000,10000,100000,1000000)

b <- bs[4]
n <- ns[4]
N <- 365*5



for(k in 1:4){
	n <- ns[k]
	for(j in 1:length(bs)){
		b <- bs[j]
		nm <- paste('SEIR_ensemble_n=',n,'_b=',format(b,nsmall=2),'.csv',sep='')
		dat <- read.csv(nm)
			
		plot(-999,xlim=c(0,900),ylim=c(0,2.0))
		for(i in 1:10){
			d <- dat[dat$trial==i,]


pdf('d:/dropbox/working/covid19/urop/growth_rate_variance.pdf',height=7,width=9)
par(mfrow=c(4,5),mar=c(2,2,1,1))

for(k in 1:4){
	n <- ns[k]
	for(j in 1:length(bs)){
		b <- bs[j]
		nm <- paste('SEIR_ensemble_n=',n,'_b=',format(b,nsmall=2),'.csv',sep='')
		dat <- read.csv(nm)
			
		plot(-999,xlim=c(0,900),ylim=c(0,2.0))
		for(i in 1:10){
			d <- dat[dat$trial==i,]
			I <- d$X3
			R <- d$X4
			Inew <- diff(I+R)

			#plot(I)
			g <- log(I[2:N]/I[1:(N-1)])
			g <- log(Inew[2:N]/Inew[1:(N-1)])
			#lines(g)
			groll <- rollapply(g/Inew[-N], width=7, FUN=sd, fill=0, align="r")
			lines(groll,col=i)
		}
	}
}

dev.off()



cs <- c(1,10,50,100,200,500)

TS=MAX <- array(NA,dim=c(1000,length(bs),length(ns),length(cs)))

for(k in 1:4){
print(k)
	n <- ns[k]
	for(j in 1:length(bs)){
		b <- bs[j]
		nm <- paste('SEIR_ensemble_n=',n,'_b=',format(b,nsmall=2),'.csv',sep='')
		dat <- read.csv(nm)
			
		for(p in 1:6){
			ii <- cs[p]

			for(i in 1:1000){
				d <- dat[dat$trial==i,]
				I <- d$X3

				MAX[i,j,k,p] <- max(I)
				
				TS[i,j,k,p] <- which(I>ii)[1]
			}
		}
	}
}

pdf('d:/dropbox/working/covid19/urop/outbreak_times_simulations.pdf',height=6,width=8)
par(mfrow=c(4,5),mar=c(2,2,1,1))
for(k in 1:4){
	for(j in 1:5){
		xx <- TS[1:1000,j,k,2]
		xx <- xx[!is.na(xx)]
		if(sum(!is.na(xx))>0) hist(xx,breaks=20,xlim=c(0,365),main='')
	}
}
dev.off()

loc=scl=shp <- matrix(NA,4,5)
pdf('d:/dropbox/working/covid19/urop/outbreak_times_simulations.pdf',height=6,width=8)
par(mfrow=c(4,5),mar=c(2,2,1,1))
for(k in 1:4){
	for(j in 1:5){
		xx <- TS[1:1000,j,k,2]
		xx <- xx[!is.na(xx)]
		fit <- fevd(xx)
		plot(fit,type='hist',main='',ylim=c(0,0.05),col='grey',legend=c())
		shp[k,j] <- fit$results$par[3]
		#if(sum(!is.na(xx))>0) hist(xx,breaks=20,xlim=c(0,365),main='')
	}
}
dev.off()


par(mfrow=c(4,5),mar=c(2,2,1,1))
for(k in 1:4){
	for(j in 1:5){
		xx <- TS[1:1000,j,k,2]
		xx <- xx[!is.na(xx)]
		fit <- fevd(xx)
		plot(fit,type='hist',main='',ylim=c(0,0.05),col='grey',legend=c())
		shp[k,j] <- fit$results$par[3]
		#if(sum(!is.na(xx))>0) hist(xx,breaks=20,xlim=c(0,365),main='')
	}
}

skew <- matrix(NA,4,5)
for(k in 1:4){
	for(j in 1:5){
		xx <- TS[1:1000,j,k,2]
		xx <- xx[!is.na(xx)]
		skew[k,j] <- skewness(xx)
		#if(sum(!is.na(xx))>0) hist(xx,breaks=20,xlim=c(0,365),main='')
	}
}



MAX <- array(NA,dim=c(1000,length(bs),length(ns)))

par(mfrow=c(4,5),mar=c(2,2,1,1))
for(k in 1:4){
	n <- ns[k]
	for(j in 1:length(bs)){
		b <- bs[j]

		nm <- paste('SEIR_ensemble_n=',n,'_b=',format(b,nsmall=2),'.csv',sep='')
		dat <- read.csv(nm)

		for(i in 1:1000){
				d <- dat[dat$trial==i,]
				I <- d$X3

				MAX[i,j,k] <- max(I)
		}
	}
}




par(mfrow=c(4,5),mar=c(2,2,1,1))
for(k in 1:4){
	for(j in 1:length(bs)){
		xx <- MAX[,j,k]
		xx <- xx[!is.na(xx) & xx > 5]
		fit <- fevd(xx)
		plot(fit,type='hist',main='',ylim=c(0,0.05),col='grey',legend=c())
	}
}





EX <- array(NA,dim=c(1000,length(bs),length(ns)))

for(k in 1:4){
print(k)
	n <- ns[k]
	for(j in 1:length(bs)){
		b <- bs[j]
		nm <- paste('SEIR_ensemble_n=',n,'_b=',format(b,nsmall=2),'.csv',sep='')
		dat <- read.csv(nm)
			
		for(i in 1:1000){
			d <- dat[dat$trial==i,]
			I <- d$X3
			
			EX[i,j,k] <- which(I==0)[1]

		}
	}
}


shp <- matrix(NA,4,5)
par(mfrow=c(4,5),mar=c(2,2,1,1))
for(k in 1:4){
	for(j in 1:length(bs)){
		xx <- EX[,j,k]
		xx <- xx[!is.na(xx) & xx > 150]
		fit <- fevd(xx)

		shp[k,j] <- fit$results$par[3]
		plot(fit,type='hist',main='',ylim=c(0,0.05),col='grey',legend=c(),xlim=c(100,600))
	}
}


skew <- matrix(NA,4,5)
par(mfrow=c(4,5),mar=c(2,2,1,1))
for(k in 1:4){
	for(j in 1:length(bs)){
		xx <- EX[,j,k]
		xx <- xx[!is.na(xx) & xx > 150]
		
		skew[k,j] <- skewness(xx)
	}
}







# If we perform binomial testing of the population, we have a variance n*p*(1-p)
# We have observed growth rate variability; 
# What explains the growth rate variance?
# What is the mean-variance relationship for growth rate?
# Contributions: 1) determinitic change in the true underlying growth rate
#                2) binomial sampling variance (decreases with increasing ntest; evaluate on deterministic sims)
#				 3) demographic stochasticity

# We can sample the deterministic solution with binomial sampling and evaluate the contribution of ntest 
# We can evaluate the variance of the stochastic simulation growth rates

#########################################################################
sims <- read.csv('d:/dropbox/working/covid19/urop/github/greg1/SEIR_ensemble_n=10000_b=0.30.csv')
sims2 <- read.csv('d:/dropbox/working/covid19/urop/github/greg1/SEIR_ensemble_n=1000000_b=0.30.csv')

sim <- sims[sims$trial==1,]
sim2 <- sims2[sims2$trial==2,]

I1 <- sim$X3
I2 <- sim2$X3

par(mfrow=c(2,2),mar=c(2,2,2,2))
	plot(I1,type='l')
	plot(I2,type='l')

	plot(log(I1[2:300]/I1[1:299]),ylim=c(-0.3,0.3),type='l')
	abline(h=0)
	plot(log(I2[2:300]/I2[1:299]),ylim=c(-0.3,0.3),type='l')
	abline(h=0)

p1 <- I1/10000
p2 <- I2/100000


sizes = round(seq(10,100,length.out=301))

Iobs1_100 <- rollapply(rbinom(301,size=100,prob=p1), 7, mean)
Iobs2_100 <- rollapply(rbinom(301,size=100,prob=p2), 7, mean)
Iobs1_1000 <- rollapply(rbinom(301,size=1000,prob=p1), 7, mean)
Iobs2_1000 <- rollapply(rbinom(301,size=1000,prob=p2), 7, mean)
#Iobs1 <- rollapply(rbinom(301,size=sizes,prob=p1), 7, mean)
#Iobs2 <- rollapply(rbinom(301,size=1000,prob=p2), 7, mean)


par(mfrow=c(4,2))
plot(Iobs1_100,type='l')
plot(Iobs2_100,type='l')
plot(log(Iobs1_100[2:300]/Iobs1_100[1:299]),type='l'); abline(h=0)
plot(log(Iobs2_100[2:300]/Iobs2_100[1:299]),type='l'); abline(h=0)
plot(Iobs1_1000,type='l')
plot(Iobs2_1000,type='l')
plot(log(Iobs1_1000[2:300]/Iobs1_1000[1:299]),type='l'); abline(h=0)
plot(log(Iobs2_1000[2:300]/Iobs2_1000[1:299]),type='l'); abline(h=0)


plot(I1[3:297],Iobs1_100)


plot(Iobs1)


plot(smooth(log(Iobs[2:300]/Iobs[1:299])))
abline(h=0)

