suf <- read.csv('d:/dropbox/working/covid19/urop/github/greg1/covid_confirmed_usafacts.csv')

d <- suf[,5:ncol(suf)]

states <- unique(suf$State)

pdf('d:/dropbox/working/covid19/urop/outbreak_times.pdf',height=6,width=8)
breaks <- seq(1,365,10)
par(mfrow=c(7,8),mar=c(1,1,1,1),cex.axis=0.7,oma=c(2,2,2,2))
for(i in 1:51){
	d <- suf[suf$State==states[i],5:ncol(suf)]
	hist(apply(d,1,function(x) which(x>10)[1]),xlim=c(0,365),main='',breaks=breaks,xaxt='n')
	mtext(adj=0.8,states[i],line=-1,cex=0.7)
	if(i>43) axis(side=1)
}
dev.off()

pdf('d:/dropbox/working/covid19/urop/outbreak_times_combined.pdf',height=4,width=5)
par(mfrow=c(1,1))
	hist(apply(suf[,5:ncol(suf)],1,function(x) which(x>10)[1]),xlim=c(0,365),main='',breaks=breaks,xlab='')
dev.off()


cs <- c(1,10,50,100,200,500)
par(mfrow=c(2,3))
for(i in 1:length(cs)){
	d <- suf[,5:ncol(suf)]
	hist(apply(d,1,function(x) which(x>cs[i])[1]),xlim=c(0,365),main='',breaks=breaks,xaxt='n')
}


