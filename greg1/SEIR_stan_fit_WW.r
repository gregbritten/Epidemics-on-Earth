library(zoo)
library(rstan)
options(mc.cores = parallel::detectCores())
Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')

range01 <- function(x){(x-min(x))/(max(x)-min(x))}
#########################################################
dat <- read.csv('d:/dropbox/working/covid19/urop/github/greg1/wastewater.csv',skip=1)
suf <- read.csv('d:/dropbox/working/covid19/urop/github/greg1/covid_confirmed_usafacts.csv')

suf_cases <- as.numeric(suf[suf[,1]==25025,5:ncol(suf)])
suf_I     <- diff(rollapply(suf_cases,width=7,FUN=mean,align='center',fill=NA))
SUF_I     <- round(suf_I)
SUF_I     <- rollapply(suf_I[43:330],width=7,FUN=mean,align='right',fill=NA)
SUF_I     <- SUF_I[6:length(SUF_I)]
SUF_I     <- round(SUF_I)
#plot(SUF_I)
SUF_I[1]  <- 1
ySUF      <- SUF_I


WW <- rollapply(na.approx(c(1,dat$south_7avg)),width=7,FUN=mean,align='right',fill=NA)
WW <- WW[!is.na(WW)]
WW <- round(range01(WW)*max(ySUF))
#plot(WW)
yWW <- WW

POP <- 8E5
x0  <- c((POP-2)/POP,1/POP,1/POP,0)

dataSUF <- list(POP=POP,
			    y=ySUF,
			    N_obs=length(ySUF),
			    t_obs=1:length(ySUF),
			    x0 = x0)

dataWW <- list(POP=POP,
			    y=yWW,
			    N_obs=length(ySUF),
			    t_obs=1:length(ySUF),
			    x0 = x0)

########################################################

mod <- stan_model('d:/dropbox/working/covid19/urop/github/greg1/stan_SEIR.stan')



