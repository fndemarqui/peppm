

library(peppm)

data(telecom)

fit1 <- with(telecom, peppm(time, status, cohesion=1, nburnin = 0, nlag = 1, npost = 100))
fit2 <- with(telecom, peppm(time, status, cohesion=2, nburnin = 0, nlag = 1, npost = 100))
fit3 <- with(telecom, peppm(time, status, cohesion=3, nburnin = 0, nlag = 1, npost = 100))
fit4 <- with(telecom, peppm(time, status, cohesion=4, nburnin = 0, nlag = 1, npost = 100))

