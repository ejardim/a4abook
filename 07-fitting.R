## ----echo=FALSE---------------------------------------------------------------
knitr::opts_chunk$set(warning = FALSE, message=FALSE)


## -----------------------------------------------------------------------------
library(FLa4a)
data(ple4)
data(ple4.indices)
data(ple4.index)
fmod <- ~ s(age, k=8) + s(year, k=30) + te(age, year, k = c(5, 15))
fit <- sca(ple4, ple4.indices, fmodel=fmod)
stk <- ple4 + fit


## ----plt01, fig.cap="Stock summary", echo=FALSE-------------------------------
plot(stk)


## ----echo=FALSE---------------------------------------------------------------
fit


## ----fit0---------------------------------------------------------------------
fmod0 <- ~s(age, k=6)+s(year, k=10)+te(age, year, k=c(3,8))
qmod0 <- list(~s(age, k = 4),
       ~s(age, k = 3),
       ~s(age, k = 3) + year,
       ~s(age, k = 3),
       ~s(age, k = 4),
       ~s(age, k = 6))
srmod0 <- ~ s(year, k=20)
vmod0 <- list(~s(age, k=4), ~1,  ~1, ~1, ~1, ~1, ~1, ~1)
n1mod0 <- ~ s(age, k=3)
fit0 <- sca(ple4, ple4.indices,
       fmodel=fmod0,
       qmodel=qmod0,
       srmodel=srmod0,
       n1model=n1mod0,
       vmodel=vmod0)
stk0 <- ple4 + fit0


## ----plt02, fig.cap="Stock summary - close to official assessment", echo=FALSE----
plot(stk0)


## ----echo=FALSE---------------------------------------------------------------
fit0


## -----------------------------------------------------------------------------
fmod1 <- ~ factor(age) + factor(year)
fit1 <- sca(ple4, ple4.indices, fmodel=fmod1, fit="MP")


## -----------------------------------------------------------------------------
fmod2 <- ~ s(age, k=4) + s(year, k=20)
fit2 <- sca(ple4, ple4.indices, fmodel=fmod2, fit="MP")


## ----warning=FALSE------------------------------------------------------------
fmod3 <- ~ s(age, k=4) + s(year, k=20) + s(as.numeric(year-age), k=10)
fit3 <- sca(ple4, ple4.indices, fmodel=fmod3, fit="MP")


## ----sep00, fig.cap="Selection pattern of separable models. Each line represents the selection pattern in a specific year. Independent age and year effects (factor), internally dependent age and year (smooth), double separable (double).", echo=FALSE----
flqs <- FLQuants(factor=harvest(fit1), smooth=harvest(fit2), double=harvest(fit3))
pset <- list(strip.background=list(col="gray90"))
xyplot(data~age|qname, groups=year, data=flqs, type="l", col=1, layout=c(3,1), ylab="fishing mortality", par.settings=pset)


## ----sep01, fig.cap="Fishing mortality of separable models. Independent age and year effects (factor), internally dependent age and year (smooth), double separable (double).", echo=FALSE----
wireframe(data~age+year|qname, data=as.data.frame(flqs), layout=c(3,1))


## -----------------------------------------------------------------------------
fmod <- ~ te(age, year, k = c(4,20))
fit <- sca(ple4, ple4.indices, fmodel=fmod)


## ----te1, fig.cap="Fishing mortality smoothed non-separable model", echo=FALSE----
wireframe(harvest(fit), zlab="F")


## -----------------------------------------------------------------------------
fmod <- ~ te(age, year, k = c(4,20)) + s(year, k = 5, by = as.numeric(age==1))
fit2 <- sca(ple4, ple4.indices, fmodel=fmod)


## ----age1, fig.cap="Fishing mortality age-year interaction model with extra age 1 smoother.", echo=FALSE----
wireframe(harvest(fit2), zlab="F")


## -----------------------------------------------------------------------------
age <- 1:10
# last age same as previous
replace(age, age>9, 9)
# all ages after age 6
replace(age, age>6, 6)
year <- 1950:2010
replace(year, year>2005, 2005)


## -----------------------------------------------------------------------------
fmod <- ~ s(replace(age, age>9, 9), k=4) + s(year, k=20)
fit <- sca(ple4, ple4.indices, fmod)


## ----ctsselage, fig.cap="F-at-age fixed above age 9", echo=FALSE--------------
wireframe(harvest(fit), zlab="F")


## -----------------------------------------------------------------------------
fmod <- ~ s(age, k=4) + s(replace(year, year>2013, 2013), k=20)
fit <- sca(ple4, ple4.indices, fmod)


## ----ctsselyear, fig.cap="F-at-age fixed for the most recent 5 years", echo=FALSE----
wireframe(data~age+year, data=harvest(fit), screen=c(z=-130, y=0, x=-60), zlab="F")


## -----------------------------------------------------------------------------
year <- 1950:2010
# two levels separated in 2000
breakpts(year, 2000)
# five periods with equal interval
breakpts(year, seq(1949, 2010, length=6))


## -----------------------------------------------------------------------------
fmod <- ~s(age, k = 3, by = breakpts(year, 1990))
fit <- sca(ple4, ple4.indices, fmod)


## ----brk, echo=FALSE, fig.cap="F-at-age in two periods using in both cases a thin plate spline with 3 knots", echo=FALSE----
wireframe(harvest(fit), zlab="F")


## -----------------------------------------------------------------------------
fmodel <- ~ s(year, k = 15, by = factor(age)) + s(age, k = 4)


## -----------------------------------------------------------------------------
fmodel <- ~ te(age, year, k = c(5, 15))


## -----------------------------------------------------------------------------
fmodel <- ~ s(age, k = 4) + s(year, k = 15) + te(age, year, k = c(3, 5))


## -----------------------------------------------------------------------------
fmod <- ~ I(1/(1+exp(-age)))
fit <- sca(ple4, ple4.indices, fmod)


## ----logistic, fig.cap="F-at-age logistic", echo=FALSE------------------------
wireframe(harvest(fit), zlab="F")


## -----------------------------------------------------------------------------
qmod <- list(~factor(age))
fit <- sca(ple4, ple4.index, qmodel=qmod)


## ----dummyage, fig.cap="Catchability age independent model", echo=FALSE-------
qhat <- predict(fit)$qmodel[[1]]
wireframe(qhat, zlab="q")


## -----------------------------------------------------------------------------
qmod <- list(~ s(age, k=4))
fit <- sca(ple4, ple4.indices[1], qmodel=qmod)


## ----smoothage, fig.cap="Catchability smoother age model", echo=FALSE---------
qhat <- predict(fit)$qmodel[[1]]
wireframe(qhat, zlab="q")


## -----------------------------------------------------------------------------
qmod <- list( ~ s(age, k=4) + year)
fit <- sca(ple4, ple4.indices[1], qmodel=qmod)


## ----qtrend, fig.cap="Catchability with a linear trend in year", echo=FALSE----
qhat <- predict(fit)$qmodel[[1]]
wireframe(qhat, zlab="q")


## -----------------------------------------------------------------------------
# simulating a biomass index (note the name of the first dimension element) using
# the ple4 biomass and an arbritary catchability of 0.001 plus a lognormal error.
dnms <- list(age="all", year=range(ple4)["minyear"]:range(ple4)["maxyear"])
bioidx <- FLIndexBiomass(FLQuant(NA, dimnames=dnms))
index(bioidx) <- stock(ple4)*0.001
index(bioidx) <- index(bioidx)*exp(rnorm(index(bioidx), sd=0.1))
range(bioidx)[c("startf","endf")] <- c(0,0)
# note the name of the first dimension element
index(bioidx)
# fitting the model
fit <- sca(ple4, FLIndices(bioidx), qmodel=list(~1))


## -----------------------------------------------------------------------------
predict(fit)$qmodel[[1]][1,drop=TRUE]


## -----------------------------------------------------------------------------
# creating the index
dnms <- list(age="all", year=range(ple4)["minyear"]:range(ple4)["maxyear"])
bioidx <- FLIndexBiomass(FLQuant(NA, dimnames=dnms))
# but now use only ages 2:4
index(bioidx) <- tsb(ple4[ac(2:4)])*0.001
index(bioidx) <- index(bioidx)*exp(rnorm(index(bioidx), sd=0.1))
range(bioidx)[c("startf","endf")] <- c(0,0)
# to pass this information to the model one needs to specify an age range
range(bioidx)[c("min","max")] <- c(2,4)
# fitting the model
fit <- sca(ple4, FLIndices(bioidx), qmodel=list(~1))


## -----------------------------------------------------------------------------
predict(fit)$qmodel[[1]][1,drop=TRUE]


## -----------------------------------------------------------------------------
idx <- ple4.index[1]
fit <- sca(ple4, FLIndices(recidx=idx), qmodel=list(~1))
# the estimated catchability is
predict(fit)$qmodel[[1]][1,drop=TRUE]


## -----------------------------------------------------------------------------
srmod <- ~ factor(year)
fit <- sca(ple4, ple4.indices, srmodel=srmod)
srmod <- ~ s(year, k=15)
fit1 <- sca(ple4, ple4.indices, srmodel=srmod)
srmod <- ~ ricker(CV=0.1)
fit2 <- sca(ple4, ple4.indices, srmodel=srmod)
srmod <- ~ bevholt(CV=0.1)
fit3 <- sca(ple4, ple4.indices, srmodel=srmod)
srmod <- ~ hockey(CV=0.1)
fit4 <- sca(ple4, ple4.indices, srmodel=srmod)
srmod <- ~ geomean(CV=0.1)
fit5 <- sca(ple4, ple4.indices, srmodel=srmod)


## ----srmod, fig.cap="Recruitment estimates since 1960 by each stock-recruitment model.", echo=FALSE----
flqs <- FLQuants(factor=stock.n(fit)[1], smother=stock.n(fit1)[1], ricker=stock.n(fit2)[1], bevholt=stock.n(fit3)[1], hockey=stock.n(fit4)[1], geomean=stock.n(fit5)[1])
flqs <- lapply(flqs, "[", j=ac(1960:2017))
xyplot(data~year, groups=qname, data=flqs, type="l", auto.key=list(points=FALSE, lines=TRUE, columns=3), ylab="No. recruits")


## -----------------------------------------------------------------------------
# reference model with constant variance for the survey index
vmod <- list(~s(age, k=3), ~1)
fit1 <- sca(ple4, ple4.index, vmodel=vmod)
# to compare - survey index variance modelled has a U-shape smoother
vmod <- list(~s(age, k=3), ~s(age, k=3))
fit2 <- sca(ple4, ple4.index, vmodel=vmod)


## ----vmod, fig.cap="Abundance index observation variance estimate", echo=FALSE----
wireframe(predict(fit1)$vmodel[[1]], zlab="variance")


## ----vmodimpact, fig.cap="Population estimates using two different variance models for the survey", echo=FALSE----
flqs <- FLQuants(constant=stock.n(window(fit1, start=1990)), smoother=stock.n(window(fit2, start=1990)))
xyplot(data~year|factor(age), groups=qname, data=flqs, type="l",
       scales=list(y=list(relation="free", draw=FALSE)),
       auto.key=list(points=FALSE, lines=TRUE, columns=2),
       par.settings=list(superpose.line=list(col=c("red", "blue"), lwd=1.5),
       strip.background=list(col="gray90")), ylab="", layout=c(5,2))


## -----------------------------------------------------------------------------
# model with smoother
n1mod <- ~s(age, k=4)
fit1 <- sca(ple4, ple4.indices, n1model=n1mod)
# model with factor
n1mod <- ~factor(age)
fit2 <- sca(ple4, ple4.indices, n1model=n1mod)


## ----ny1, fig.cap="Nay=1 models", echo=FALSE----------------------------------
flqs <- FLQuants(smother=stock.n(fit1)[,1], factor=stock.n(fit2)[,1])
pset <- list(superpose.line=list(col=c("red", "blue")))
lgnd <- list(points=FALSE, lines=TRUE, space='right')
xyplot(data~age, groups=qname, data=flqs, type="l", auto.key=lgnd, par.settings=pset, ylab="")


## ----n1modimpact, fig.cap="Population estimates using two different variance models", echo=FALSE----
flqs <- FLQuants(smother=stock.n(fit1), factor=stock.n(fit2))
xyplot(data~year|factor(age), groups=qname, data=flqs, type="l",
       scales=list(y=list(relation="free", draw=FALSE)),
       auto.key=list(points=FALSE, lines=TRUE, columns=2),
       par.settings=list(superpose.line=list(col=c("red", "blue"), lwd=1.5),
       strip.background=list(col="gray90")), ylab="", layout=c(5,2))


## -----------------------------------------------------------------------------
fmod <- ~ factor(age) + s(year, k=10, by = breakpts(age, c(0:8)))
fit <- sca(ple4, ple4.indices, fmod)


## ----ageind, fig.cap="F-at-age as thin plate spline with 3 knots for each age", echo=FALSE----
wireframe(harvest(fit), zlab="F")


## -----------------------------------------------------------------------------
fmodel <- ~ s(age, k = 4) + s(pmax(year - age, 1957), k = 10) + s(year, k = 10)
fit <- sca(ple4, ple4.indices, fmodel=fmodel)


## ----coh, echo=FALSE, fig.cap="F-at-age with a cohort effect.", echo=FALSE----
wireframe(harvest(fit), zlab="F")


## -----------------------------------------------------------------------------
vmod <- list(
       ~ s(age, k = 3) + year,
       ~1, ~1, ~1, ~1, ~1, ~1
       )
fit <- sca(ple4, ple4.indices, vmodel=vmod)


## ----vm, echo=FALSE, fig.cap="Catch at age variance model with a year effect.", echo=FALSE----
wireframe(predict(fit)$vmodel[[1]], zlab="variance")


## -----------------------------------------------------------------------------
fmod <- ~s(age, k = 3, by = breakpts(age, 5)) + s(year, k = 10)
fit <- sca(ple4, ple4.indices, fmodel = fmod)


## ----danai01, echo=FALSE, fig.cap="Smoothers fitted to two sets of ages, 1 to 4 and 5 to 10.", echo=FALSE----
wireframe(harvest(fit), zlab="F")


## -----------------------------------------------------------------------------
stk <- ple4
idx <- ple4.indices[1]
# cv of observed catches
varslt <- catch.n(stk)
varslt[] <- 0.4
catch.n(stk) <- FLQuantDistr(catch.n(stk), varslt)
# cv of observed indices
varslt <- index(idx[[1]])
varslt[] <- 0.2
index.var(idx[[1]]) <- varslt
# run
fit1 <- sca(stk, idx, fmodel=fmod0, qmodel=qmod0, srmodel=srmod0, vmodel=vmod0, n1model=n1mod0)
flqs <- FLQuants(nowgt=stock.n(fit0), extwgt=stock.n(fit1))


## ----likwgt, fig.cap="Stock summary of distinct likelihood weightings", echo=FALSE----
#xyplot(data~year|factor(age), groups=qname, data=flqs, type="l", scales=scl, auto.key=lgnd, par.settings=pset, ylab="")


## ----likwgtimpact, fig.cap="Population estimates using two different variance models", eho=FALSE----
flsts <- FLStocks(nowgt=ple4+fit0, wgt=ple4 + fit1)
plot(flsts)


## ----echo=FALSE---------------------------------------------------------------
fit0 <- sca(ple4, ple4.indices[1], fmodel=fmod0, qmodel=qmod0, srmodel=srmod0, vmodel=vmod0, n1model=n1mod0)
(fitSumm(fit1)/fitSumm(fit0))[c(2,8,9),]


## -----------------------------------------------------------------------------
nao <- read.table("https://www.cpc.ncep.noaa.gov/products/precip/CWlink/pna/norm.nao.monthly.b5001.current.ascii.table", skip=1, fill=TRUE, na.strings="-99.90")
dnms <- list(quant="nao", year=1950:2024, unit="unique", season=1:12, area="unique")
nao <- FLQuant(unlist(nao[,-1]), dimnames=dnms, units="nao")
nao <- seasonMeans(trim(nao, year=dimnames(stock.n(ple4))$year))


## -----------------------------------------------------------------------------
srmod <- ~ s(nao, k=10)
fit2 <- sca(ple4, ple4.indices[1], qmodel=list(~s(age, k=4)), srmodel=srmod, covar=FLQuants(nao=nao))


## ----naor, echo=FALSE, fig.cap="Recruitment model with covariates. Using the NAO index as a recruitment index.", echo=FALSE----
flqs <- FLQuants(simple=stock.n(fit)[1], covar=stock.n(fit2)[1])
xyplot(data~year, groups=qname, data=flqs, type="l",
       auto.key=list(points=FALSE, lines=TRUE, columns=2),
       par.settings=list(superpose.line=list(col=c("red", "blue"), lwd=1.5),
       strip.background=list(col="gray90")), ylab="")


## -----------------------------------------------------------------------------
srmod <- ~ ricker(a=~nao, CV=0.25)
fit3 <- sca(ple4, ple4.indices[1], qmodel=list(~s(age, k=4)), srmodel=srmod, covar=FLQuants(nao=nao))


## ----naor2, echo=FALSE, fig.cap="Recruitment model with covariates. Using the NAO index as a covariate for the stock-recruitment model parameters.", echo=FALSE----
flqs <- FLQuants(simple=stock.n(fit)[1], covar=stock.n(fit3)[1])
xyplot(data~year, groups=qname, data=flqs, type="l",
       auto.key=list(points=FALSE, lines=TRUE, columns=2),
       par.settings=list(superpose.line=list(col=c("red", "blue"), lwd=1.5),
       strip.background=list(col="gray90")), ylab="")


## ----eval=FALSE---------------------------------------------------------------
# fit1 <- sca(ple4, ple4.indices, wkdir="fit1run")


## ----missing obs--------------------------------------------------------------
fit <- sca(ple4, ple4.indices)
ple4_missing <- ple4
catch.n(ple4_missing)[ac(1:2), "2013"] <- NA
fit_missing <- sca(ple4_missing, ple4.indices)


## ----obsmissing, echo = FALSE, fig.cap="Stock estimates with missing observations."----
pred <- ple4 + simulate(fit, 1000)
pred_missing <- ple4 + simulate(fit_missing, 1000)
flqs <- FLQuants(base = catch.n(pred)[ac(1:2), ac(2011:2015)], missing = catch.n(pred_missing)[ac(1:2), ac(2011:2015)])
bwplot(data ~ qname | factor(year) + factor(age), data = as.data.frame(flqs), scales = "free", auto.key = T, ylab = "Catch at age", layout=c(5,2), par.settings=list(box.rectangle = list(col = "black"), box.umbrella = list(col = "black"), plot.symbol = list(col = "black")))


## -----------------------------------------------------------------------------
# bevholt s/r CV was tweaked to give best results for the example
fit2 <- sca(ple4, ple4.indices, srmodel=~bevholt(CV=0.16))
fit_missing2 <- sca(ple4_missing, ple4.indices, srmodel=~bevholt(CV=0.16))


## ----obsmissing2, echo = FALSE, fig.cap="Stock estimates with missing observations."----
pred <- ple4 + simulate(fit2, 1000)
pred_missing <- ple4 + simulate(fit_missing2, 1000)
flqs2 <- FLQuants(base = catch.n(pred)[ac(1:2), ac(2011:2015)], missing = catch.n(pred_missing)[ac(1:2), ac(2011:2015)])
bwplot(data ~ qname | factor(year) + factor(age), data = as.data.frame(flqs2), scales = "free", auto.key = T, ylab = "Catch at age", layout=c(5,2), par.settings=list(box.rectangle = list(col = "black"), box.umbrella = list(col = "black"), plot.symbol = list(col = "black")))

