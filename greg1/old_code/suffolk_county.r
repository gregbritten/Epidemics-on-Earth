library(zoo)

dat <- read.csv('d:/dropbox/working/covid19/urop/github/greg1/wastewater.csv',skip=1)
dat$time <- decimal_date(mdy(dat$date))

suf <- read.csv('d:/dropbox/working/covid19/urop/github/greg1/covid_confirmed_usafacts.csv')


suf_cases <- as.numeric(suf[suf[,1]==25025,5:ncol(suf)])
suf_time <- decimal_date(mdy(substring(colnames(suf[,5:ncol(suf)]),2)))

suf_I <- diff(rollapply(suf_cases,width=7,FUN=mean,align='center',fill=NA))

mar_cases <- as.numeric(suf[suf[,1]==4013,5:ncol(suf)])


plot(suf_I)
SUF_I <- round(suf_I)
SUF_I <- rollapply(suf_I[43:330],width=7,FUN=mean,align='right',fill=NA)
SUF_I <- SUF_I[6:length(SUF_I)]
SUF_I <- round(SUF_I)


WW <- rollapply(na.approx(dat$south_7avg),width=7,FUN=mean,align='right',fill=NA)
plot(WW)

round(SUF_I)


par(mfrow=c(2,2),mar=c(2,2,2,2),oma=c(2,2,2,2))
plot(suf_time[-1],suf_I,type='l',xlim=c(2020,2021),ylim=c(0,700))
points(suf_time[-1],diff(suf_cases),pch=19,cex=0.2)

plot(dat$time,dat$south_7avg,type='l',col='blue',xlim=c(2020,2021))
points(dat$time,dat$south,pch=19,cex=0.2,col='blue')














