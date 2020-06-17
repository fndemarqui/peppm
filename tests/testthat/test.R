

library(peppm)

data(telecom)
head(telecom)

fit <- with(telecom, peppm(time, status, cohesion=1, nburnin = 0, nlag = 1, npost = 100))
names(fit)
