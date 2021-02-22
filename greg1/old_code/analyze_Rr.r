
ma <- read.csv('d:/dropbox/working/covid19/urop/github/greg1/MA_try.csv')
mapos <- read.csv('d:/dropbox/working/covid19/urop/github/greg1/MA_try_pos.csv')

ma <- ma[,colnames(ma) %in% colnames(mapos)]

Rs <- mapos[,1]

rmod = rmodpos <- numeric(ncol(mapos))
for(i in 3:ncol(mapos)){
	rmod[i] = Rs[ma[,i]==max(ma[,i])]
	rmodpos[i] = Rs[mapos[,i]==max(mapos[,i])]
}


par(mfrow=c(1,1))
plot(rmod,type='l',xaxt='n')
lines(rmodpos,type='l',col='red')
abline(h=1,lty=2)
legend('topleft',legend=c('Raw','From positive rate'),lty=1,col=c('black','red'),bty='n')

labi <- seq(1,length(rmod),14)
axis(side=1,at=labi,labels=substring(colnames(ma),2)[labi],cex.axis=0.5)




