dat <- read.csv('d:/dropbox/working/covid19/urop/github/greg1/wastewater.csv',skip=1)
suf <- read.csv('d:/dropbox/working/covid19/urop/github/greg1/covid_confirmed_usafacts.csv')

suf_cases <- as.numeric(suf[suf[,1]==25025,5:ncol(suf)])
suf_date <- mdy(substring(colnames(suf[,5:ncol(suf)]),2))


suf_I     <- rollapply(suf_cases,width=7,FUN=mean,align='center')
suf_t     <- suf_date[3:337]
suf_I     <- round(suf_I)
SUF <- data.frame(cases=suf_I,time=suf_t)

ww_I <- round(rollapply(na.approx(c(1,dat$south)),width=7,FUN=mean))
ww_t <- mdy(dat$date[1:288])
WW <- data.frame(cases=ww_I,time=ww_t)

X <- merge(SUF,WW,by='time',all=TRUE)
plot(X$cases.y/X$cases.x,ylim=c(0,0.1))
plot(log(X$cases.y/X$cases.x))


df <- data.frame(date=c(suf_t,ww_t),state=c(rep("SUF",length(suf_I)),rep("WW",length(ww_I))),cases=c(suf_I,cumsum(ww_I)))
dftib <- as_tibble(df)




#estimates_all <- covid_cases %>%
estimates_all <- dftib %>%
  filter(date >= "2020-03-01") %>%
  group_by(state) %>%
  # Ignore states that have not reached 100 infections
  filter(max(cases) > 100 ) %>%
  group_split() %>%
  map_df(~ {
    .x %>%
      smooth_new_cases() %>%
      compute_likelihood() %>%
      compute_posterior() %>%
      estimate_rt()
  }) %>%
  ungroup()


estimates_all %>%
  plot_estimates() +
  facet_wrap(~ state, ncol = 4) +
  labs(subtitle = "")

