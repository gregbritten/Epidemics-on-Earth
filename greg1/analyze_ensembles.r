library(zoo)
library(extRemes)
library(dplyr)
library(viridis)
library(RColorBrewer)
#setwd('~/dropbox/working/covid19/urop/simulations/')
setwd('~/dropbox/working/covid19/urop/simulations_SIR/')
files <- list.files()

#betas <- as.numeric(substr(files,26,29))
betas <- as.numeric(substr(files,25,28))
gamma <- 0.1
R0s <- betas/gamma

#EX <- data.frame(trial=1:10000)
EXsir <- data.frame(trial=1:10000)
for(i in 1:length(files)){
print(i)
  dat <- read.csv(files[i])
  #ex <- ddply(dat, .(trial), summarize, col=which(X3==0)[1]) #X3 is I in SEIR
  ex <- ddply(dat, .(trial), summarize, col=which(X2==0)[1]) #X2 is I in SIR
  colnames(ex)[2] <- R0s[i]

  #EX <- merge(EX,ex)
  EXsir <- merge(EXsir,ex)
}

# write.csv(file='~/dropbox/working/COVID19/urop/results/EX.csv',EX)
# save(file='~/dropbox/working/COVID19/urop/results/EX.rdata',EX)
write.csv(file='~/dropbox/working/COVID19/urop/results/EXsir.csv',EXsir)
save(file='~/dropbox/working/COVID19/urop/results/EXsir.rdata',EXsir)

#################################################
## SUMMARY STATISTICS ###########################
#################################################
EX <- EXsir

df <- data.frame(iep=numeric(),iemean=numeric(),iesd=numeric(),ieskew=numeric(),iekurt=numeric(),
                 eep=numeric(),eemean=numeric(),eesd=numeric(),eeskew=numeric(),eekurt=numeric(),
                 obp=numeric(),obmean=numeric(),obsd=numeric(),obskew=numeric(),obkurt=numeric())

Qie=Qee=Qob <- matrix(NA,nrow=ncol(EX)-1,ncol=5)

for(i in 2:ncol(EX)){
  ex <- EX[,i]
  ob <- ex[ex>70 & !is.na(ex)]
  ie <- ex[ex==2 & !is.na(ex)]
  ee <- ex[ex>2 & ex<=70 & !is.na(ex)]
  
  dftmp <- data.frame(iep=length(ie)/1E4, iemean=mean(ie), iesd=sd(ie), ieskew=skewness(ie), iekurt=kurtosis(ie),
                      eep=length(ee)/1E4, eemean=mean(ee), eesd=sd(ee), eeskew=skewness(ee), eekurt=kurtosis(ee),
                      obp=length(ob)/1E4, obmean=mean(ob), obsd=sd(ob), obskew=skewness(ob), obkurt=kurtosis(ob))
  
  Qie[i-1,] <- quantile(ie,probs=c(0.05,0.25,0.5,0.75,0.95))
  Qee[i-1,] <- quantile(ee,probs=c(0.05,0.25,0.5,0.75,0.95))
  Qob[i-1,] <- quantile(ob,probs=c(0.05,0.25,0.5,0.75,0.95))
  
  df <- rbind(df,dftmp)  
}

#cols <- viridis(3)
cols <- brewer.pal(3, 'Dark2')
pdf('~/dropbox/working/COVID19/UROP/plots/class_probs.pdf',height=5,width=6)
par(mfrow=c(1,1),cex.axis=0.8)
  plot(R0s,df$iep,type='l',ylim=c(0,1),xaxt='n',col=cols[1],lwd=1.5,xlab='',ylab='')
  axis(side=1,at=R0s,cex.axis=0.4,padj=0)
  lines(R0s,df$eep,type='l',ylim=c(0,1),col=cols[2],lwd=1.5)
  lines(R0s,df$obp,type='l',ylim=c(0,1),col=cols[3],lwd=1.5)
  abline(v=1,lty=2)
  legend(1.5,1.05,legend=c('Immediate Extinction','Early Extinction','Outbreak'),bty='n',lty=1,col=cols)
  mtext(side=1,expression(italic(R[0])),line=2.5)
  mtext(side=2,expression('Probability'),line=2.5)
dev.off()
  
  
pdf('~/dropbox/working/COVID19/UROP/plots/quantiles_early_exintinctions.pdf',height=4,width=6)
par(mfrow=c(1,1),cex.axis=0.8)
plot(R0s,Qee[,3],type='l',ylim=c(0,55),xaxt='n',col=cols[1],xlab='',ylab='',,lwd=1)
axis(side=1,at=R0s,cex.axis=0.4)
matplot(R0s,Qee[,c(1,5)],type='l',add=TRUE,lty=3,col=cols[2],lwd=1.5)
matplot(R0s,Qee[,c(2,4)],type='l',add=TRUE,lty=2,col=cols[3],lwd=1.5)
mtext(side=1,expression(italic(R[0])),line=2.5)
mtext(side=2,expression('Day of Extinction'),line=2.5)
abline(v=1,lty=2)
legend(11,55,legend=c('0.95','0.75','0.50','0.25','0.05'),bty='n',col=c(cols[2],cols[3],cols[1],cols[3],cols[2]),
        lty=c(3,2,1,2,3))
dev.off()

pdf('~/dropbox/working/COVID19/UROP/plots/quantiles_outbreaks.pdf',height=4,width=6)
par(mfrow=c(1,1),cex.axis=0.8)
plot(R0s,Qob[,3],type='l',ylim=c(70,1200),xaxt='n',col=cols[1],xlab='',ylab='')
axis(side=1,at=R0s,cex.axis=0.4)
matplot(R0s,Qob[,c(1,5)],type='l',add=TRUE,lty=3,col=cols[2],lwd=1)
matplot(R0s,Qob[,c(2,4)],type='l',add=TRUE,lty=2,col=cols[3],lwd=1)
legend(11,1200,legend=c('0.95','0.75','0.50','0.25','0.05'),bty='n',col=c(cols[2],cols[3],cols[1],cols[3],cols[2]),
       lty=c(3,2,1,2,3))
mtext(side=2,expression('Day of Extinction'),line=2.5)
mtext(side=1,expression(italic(R[0])),line=2.5)
abline(v=1,lty=2)
dev.off()



pdf('~/dropbox/working/COVID19/UROP/plots/skewness_kurtosis.pdf',height=5,width=8)
par(mfrow=c(2,3),mar=c(1,1,1,1),oma=c(3,3,3,3))
plot(R0s,df$eeskew,type='l',ylim=c(0,5),xaxt='n')
  axis(side=1,at=R0s,labels=NA)
  mtext(adj=0.1,'Early Extinctions')
  mtext(side=2,'Skewnwss',line=2.5)
  abline(v=1,lty=2)
plot(R0s,df$obskew,type='l',ylim=c(-4,3),xaxt='n')
  axis(side=1,at=R0s,labels=NA)
  mtext(adj=0.1,'Outbreaks')
  abline(v=1,lty=2)
plot(R0s,df$obskew,type='l',ylim=c(0.4,1.1),xaxt='n')
  axis(side=1,at=R0s,labels=NA)
  mtext(adj=0.1,'Outbreaks (zoomed in)')
  abline(v=1,lty=2)
plot(R0s,df$eekurt,type='l',ylim=c(0,35),xaxt='n')
  mtext(side=2,'Kurtosis',line=2.5)
  axis(side=1,at=R0s,cex.axis=0.5)
plot(R0s,df$obkurt,type='l',ylim=c(-2,15),xaxt='n')
  axis(side=1,at=R0s,cex.axis=0.5)
plot(R0s,df$obkurt,type='l',ylim=c(1,6),xaxt='n')
  axis(side=1,at=R0s,cex.axis=0.5)
  mtext(side=1,expression(italic(R[0])),line=2,outer=TRUE)
dev.off()

# 
# par(mfrow=c(1,1))
# plot(R0s,df$eemean,type='l',ylim=c(0,35))
# lines(R0s,df$eemean + df$eesd,type='l',ylim=c(0,1),lty=2)
# lines(R0s,df$eemean - df$eesd,type='l',ylim=c(0,1),lty=2)
# 
# plot(R0s,df$obmean,type='l')
# lines(R0s,df$obmean + df$obsd,type='l',lty=2)
# lines(R0s,df$obmean - df$obsd,type='l',lty=2)
# 
# 
# plot(R0s,df$eeskew,type='l',ylim=c(0,5))
# plot(R0s,df$obskew,type='l',ylim=c(0,1))


##################################################
## FULL HISTOGRAM ################################
##################################################

pdf('~/dropbox/working/COVID19/UROP/plots/full_histograms.pdf',height=7,width=9)
breaks=seq(-1,2000,length.out=200)
par(mfrow=c(7,8),mar=c(1,1,0,0),oma=c(4,4,2,2),cex.axis=0.7,cex.lab=0.7)
for(i in 2:ncol(EX)){
  hist(EX[,i],xlim=c(0,1000),ylim=c(0,3000),main='',breaks=breaks,xaxt='n',yaxt='n')
  mtext(bquote('R'[0]~'='~.(R0s[i-1])),cex=0.6,line=-2,adj=0.8)
  if(i>(7*7-1)){axis(side=1,cex.axis=0.7)}
  if(i %in% seq(2,7*8,8)){axis(side=2,cex.axis=0.7)}
}
mtext(side=1,outer=TRUE,'Extinction Day',line=1.5)
mtext(side=2,outer=TRUE,'Number of Stochastic Trials',line=1.5)
dev.off()


## OUTBREAKS #######################
pdf('~/dropbox/working/COVID19/UROP/plots/outbreak_histograms.pdf',height=7,width=9)
breaks=seq(-1,2000,length.out=200)
par(mfrow=c(7,8),mar=c(1,1,0,0),oma=c(4,4,2,2))
for(i in 2:ncol(EX)){
  hist(EX[,i][EX[,i]>100],xlim=c(100,1000),ylim=c(0,2000),main='',breaks=breaks,yaxt='n',xaxt='n')
  mtext(bquote('R'[0]~'='~.(R0s[i-1])),cex=0.6,line=-2,adj=0.8)
  if(i>(7*7-1)){axis(side=1,cex.axis=0.7)}
  if(i %in% seq(2,7*8,8)){axis(side=2,cex.axis=0.7)}
}
mtext(side=1,outer=TRUE,'Extinction Day',line=1.5)
mtext(side=2,outer=TRUE,'Number of Stochastic Trials',line=1.5)
dev.off()

## IMMEDIATE EXTINCTION #######################
pdf('~/dropbox/working/COVID19/UROP/plots/immediate_histograms.pdf',height=7,width=9)
breaks=seq(-1,200,length.out=50)
par(mfrow=c(7,8),mar=c(1,1,0,0),oma=c(4,4,2,2))
for(i in 2:ncol(EX)){
  hist(EX[,i][EX[,i]==2],ylim=c(0,10000),xlim=c(0,100),main='',breaks=breaks,yaxt='n',xaxt='n')
  if(i>(7*7-1)){axis(side=1,cex.axis=0.7)}
  if(i %in% seq(2,7*8,8)){axis(side=2,cex.axis=0.7)}
}
mtext(side=1,outer=TRUE,'Extinction Day',line=1.5)
mtext(side=2,outer=TRUE,'Number of Stochastic Trials',line=1.5)
dev.off()

## EARLY EXTINCTION #######################
pdf('~/dropbox/working/COVID19/UROP/plots/early_histograms.pdf',height=7,width=9)
breaks=seq(-1,200,length.out=70)
par(mfrow=c(7,8),mar=c(1,1,0,0),oma=c(4,4,2,2))
for(i in 2:ncol(EX)){
  hist(EX[,i][EX[,i]>2 & EX[,i]<100],ylim=c(0,1000),xlim=c(0,100),main='',breaks=breaks,yaxt='n',xaxt='n')
  if(i>(7*7 -1 )){axis(side=1,cex.axis=0.7)}
  if(i %in% seq(2,7*8,8)){axis(side=2,cex.axis=0.7)}
}
mtext(side=1,outer=TRUE,'Extinction Day',line=1.5)
mtext(side=2,outer=TRUE,'Number of Stochastic Trials',line=1.5)
dev.off()

## EARLY EXTINCTION #######################
pdf('~/dropbox/working/COVID19/UROP/plots/immediate_early_histograms.pdf',height=7,width=9)
breaks=seq(-1,200,length.out=70)
par(mfrow=c(7,8),mar=c(1,1,0,0),oma=c(4,4,2,2))
for(i in 2:ncol(EX)){
  hist(EX[,i][EX[,i]>1 & EX[,i]<100],ylim=c(0,2000),xlim=c(0,100),main='',breaks=breaks,yaxt='n',xaxt='n')
  if(i>(7*7 -1 )){axis(side=1,cex.axis=0.7)}
  if(i %in% seq(2,7*8,8)){axis(side=2,cex.axis=0.7)}
}
mtext(side=1,outer=TRUE,'Extinction Day',line=1.5)
mtext(side=2,outer=TRUE,'Number of Stochastic Trials',line=1.5)
dev.off() 