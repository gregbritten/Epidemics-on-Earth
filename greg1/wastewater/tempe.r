
d <- read.csv('d:/dropbox/working/covid19/urop/github/greg1/covid-19_cases_by_zip_code.csv')

zips <- c(85202,85280,85281,85282,85283,85284,85285,85287)	

d <- d[d$POSTCODE %in% zips,]
d$time <- decimal_date(ymd(substring(d$create_date,1,10)))

dates <- unique(d$create_date)
ddf <- data.frame(date=character(),sum=numeric(),time=numeric())

for(i in 1:length(dates)){
	dd <- d[d$create_date==dates[i],]
	ddf <- rbind(ddf,data.frame(date=dates[i],sum=sum(dd$CaseCount,na.rm=TRUE),time=dd$time[1]))
}	

par(new=TRUE)
plot(ddf$time[-1],diff(ddf$sum),xlim=c(2020,2021),ylim=c(0,1000),type='l')


