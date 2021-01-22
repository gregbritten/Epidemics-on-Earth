library(lubridate)
library(zoo)
library(dplyr)
library(plyr)
library(mgcv)

dat <- read.csv('d:/dropbox/working/covid19/urop/github/greg1/wastewater.csv',skip=1)
dat$time <- decimal_date(mdy(dat$date))

suf <- read.csv('d:/dropbox/working/covid19/urop/github/greg1/covid_confirmed_usafacts.csv')


suf_cases <- as.numeric(suf[suf[,1]==25025,5:ncol(suf)])
suf_time <- decimal_date(mdy(substring(colnames(suf[,5:ncol(suf)]),2)))

suf_I <- diff(rollapply(suf_cases,width=7,FUN=mean,align='center',fill=NA))

mar_cases <- as.numeric(suf[suf[,1]==4013,5:ncol(suf)])


# par(mfrow=c(2,2),mar=c(2,2,2,2),oma=c(2,2,2,2))
	# plot(suf_time[-1],suf_I,type='l',xlim=c(2020,2021),ylim=c(0,700))
		# points(suf_time[-1],diff(suf_cases),pch=19,cex=0.2)
	# plot(dat$time-3/365,dat$south_7avg,type='l',col='blue',xlim=c(2020,2021))
		# points(dat$time,dat$south,pch=19,cex=0.2,col='blue')
	# plot(suf_time[-1],suf_I,type='l',xlim=c(2020,2021),ylim=c(0,700))
	# par(new=TRUE)
	# plot(dat$time-3/365,dat$south_7avg,type='l',col='blue',xlim=c(2020,2021),xaxt='n',yaxt='n')
		# axis(side=4,col='blue')

	# plot(suf_time[-1],suf_I,type='l',xlim=c(2020,2021),ylim=c(0,700))
		# points(suf_time[-1],diff(suf_cases),pch=19,cex=0.2)
	# par(new=TRUE)
	# plot(dat$time-3/365,dat$south_7avg,type='l',col='blue',xlim=c(2020,2021),xaxt='n',yaxt='n')
		# axis(side=4,col='blue')
	# points(dat$time,dat$south,pch=19,cex=0.2,col='blue')

pdf('d:/dropbox/working/covid19/urop/RNA_cases.pdf',width=5,height=6)
par(mfrow=c(3,1),mar=c(2,2,2,2),oma=c(2,2,2,2))
	plot(suf_time[-1],suf_I,type='l',xlim=c(2020,2021),ylim=c(0,700))
		points(suf_time[-1],diff(suf_cases),pch=19,cex=0.2)
	plot(dat$time-3/365,na.approx(dat$south_7avg,na.rm=FALSE),type='l',col='blue',xlim=c(2020,2021))
		points(dat$time,dat$south,pch=19,cex=0.2,col='blue')
	plot(suf_time[-1],suf_I,type='l',xlim=c(2020,2021),ylim=c(0,700))
	par(new=TRUE)
	plot(dat$time-3/365,na.approx(dat$south_7avg,na.rm=FALSE),type='l',col='blue',xlim=c(2020,2021),xaxt='n',yaxt='n')
		axis(side=4,col='blue')
dev.off()


X1 <- data.frame(time = substring(date_decimal(dat$time - 3/365),1,10), x1=na.approx(dat$south_7avg,na.rm=FALSE))
X2 <- data.frame(time = substring(date_decimal(suf_time[-1]),1,10), x2 = suf_I)

X <- merge(X1,X2,by='time')
ii <- 29:nrow(X)

plot(decimal_date(ymd(X$time[ii])), X$x2[ii]/X$x1[ii],type='l',xlab='',ylab='# Cases per RNA copy')
abline(lm(X$x2[ii]/X$x1[ii] ~ decimal_date(ymd(X$time[ii]))))


###################################################
fit <- gam(diff(suf_cases) ~ s(suf_time[-1],k=100,sp=0.005))

pred <- predict(fit)
pred[pred<0] <- 0


par(mfrow=c(3,1),mar=c(2,2,2,2),oma=c(2,2,2,2))
plot(suf_time[-1],pred,type='l',xlim=c(2020,2021),ylim=c(0,700))
points(suf_time[-1],diff(suf_cases),pch=19,cex=0.2)

fit2 <- gam(south ~ s(time,k=100,sp=0.5), data=dat)
pred2 <- predict(fit2, newdata=list(time=dat$time))

plot(dat$time,pred2,type='l',xlim=c(2020,2021),col='blue')
points(dat$time,dat$south,pch=19,cex=0.2,col='blue')

plot(dat$time,pred2,type='l',xlim=c(2020,2021),col='blue')
lines(suf_time[-1],pred,type='l',xlim=c(2020,2021))




n <- length(dat$north)

ng <- log(dat$north_7avg[2:n]/dat$north_7avg[1:(n-1)])

par(mfrow=c(2,2),mar=c(2,2,2,2))
plot(1:315,dat$north_7avg,type='l')
plot(1:315,dat$south_7avg,type='l')

plot(1:314,ng,ylim=c(-0.5,0.5))


################################################################
tempe <- read.csv('d:/dropbox/working/covid19/urop/github/greg1/covid_wastewater_tempe.csv')
tempe$date <- unlist(strsplit(tempe$Sample_Date, " "))[seq(1,564,2)]
tempe$time <- decimal_date(mdy(tempe$date))

sum_tempe <- ddply(tempe,.(Site,time),summarize,mean=mean(Gene_Copies_per_Liter,na.rm=TRUE))
sum_tempe <- sum_tempe[2:nrow(sum_tempe),]

sum_tempe$mean[sum_tempe$mean > 4177154] <- NA

plot(sum_tempe$time,sum_tempe$mean)
fit <- loess(mean ~ time, data=sum_tempe,sp=0.2)
isort <- order(fit$x)
lines(fit$x[isort],fit$fitted[isort])


fit <- gam(log10(mean) ~ s(time),data=sum_tempe)
plot(10^predict(fit,newdata=list(time=seq(2020,2021,length.out=200))))


pdf('d:/dropbox/working/covid19/urop/tempe.pdf',height=4,width=5)
par(mfrow=c(1,1),cex.axis=0.8)
regions <- unique(sum_tempe$Site)
for(i in 1:length(regions)){
	dd <- sum_tempe[sum_tempe$Site==regions[i],]
	if(i==1){
		plot(dd$time,dd$mean,type='l',xlim=c(2020,2021),col=i+1,xlab='',ylab='')
	}else{
		lines(dd$time,dd$mean,col=i+1)
	}
}
#lines(fit$x[isort],fit$fitted[isort],lwd=2)

#par(new=TRUE)
#plot(fit$x[isort],fit$fitted[isort],lwd=2,xlim=c(2020,2021),type='l',col='blue')
par(new=TRUE)
plot(suf_time[-1],diff(mar_cases),xlim=c(2020,2021),type='l',lwd=2,lty=1,yaxt='n',xlab='',ylab='',col=adjustcolor('black',alpha.f=0.7))
axis(side=4)

dev.off()


plot(fit$x[isort],fit$fitted[isort],lwd=2,type='l',xlim=c(2020,2021))
par(new=TRUE)
plot(suf_time[-1],diff(mar_cases),xlim=c(2020,2021),type='l',lwd=2,lty=2,col=)












