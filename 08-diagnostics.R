## ----echo=FALSE, message=FALSE, warning=FALSE---------------------------------
library(FLa4a)
data(ple4)
data(ple4.indices)


## -----------------------------------------------------------------------------
# use single indices and set plus group at 9
idx <- ple4.indices[c("BTS-Combined (all)")]
idx[[1]] <- idx[[1]][1:9]
stk <- setPlusGroup(ple4, 9)
iy <- 2008
# fit
fmod <- ~s(age, k=6)+s(year, k=10)+te(age, year, k=c(3,8))
qmod <- list(~s(age, k=4))
n1mod <- ~s(age, k=5)
vmod <- list(~s(age, k=4), ~1)
fit <- sca(stk, idx, fmodel=fmod, qmodel=qmod, n1model=n1mod, vmodel=vmod)


## -----------------------------------------------------------------------------
# residuals
d_s <- residuals(fit, stk, idx)


## ----res, fig.cap="Standardized residuals for abundance indices and catch numbers (catch.n). Each panel is coded by age class, dots represent standardized residuals and lines a simple smoother."----
plot(d_s)


## ----resaux, fig.cap="Standardized residuals for catch numbers (catch.n). Each panel is coded by age class, dots represent standardized residuals and lines a regression fit."----
# shorten time series
plot(d_s[1], auxline=c("r","g"))


## ----resp, fig.cap="Pearson residuals for abundance indices and catch numbers (catch.n). Each panel is coded by age class, dots represent standardized residuals and lines a simple smoother."----
d_p <- residuals(fit, stk, idx, type='pearson')
plot(d_p)


## ----resr, fig.cap="Raw residuals for abundance indices and catch numbers (catch.n). Each panel is coded by age class, dots represent standardized residuals and lines a simple smoother."----
d_r <- residuals(fit, stk, idx, type='deviances')
plot(d_r)


## ----resy, fig.cap="Standardized residuals for abundance indices and catch numbers (catch.n). Each panel is coded by year, dots represent standardized residuals and lines a simple smoother."----
# shorten time series for demonstration purposes
d_ss <- window(d_s, start=iy)
plot(d_ss, by='age', auxline=c("h", "g"))


## ----bub, fig.cap="Bubbles plot of standardized residuals for abundance indices and for catch numbers (catch.n)."----
bubbles(d_s)


## ----qq, fig.cap="Quantile-quantile plot of standardized residuals for abundance indices and catch numbers (catch.n). Each panel is coded by age class, dots represent standardized residuals and lines the normal distribution quantiles."----
qqmath(d_s)


## ----echo=FALSE---------------------------------------------------------------
fit <- window(fit, start=iy)
stk <- window(stk, start=iy)
idx <- window(idx, start=iy)


## ----selplt, fig.cap="Predict and observed catch-at-age"----------------------
plot(fit, stk)


## ----idxplt, fig.cap="Predict and observed abundance-at-age"------------------
plot(fit, idx)


## ----catchdiag, fig.cap="Diagnostics for age aggregated catch in weight", fig.height=10, fig.asp=1, out.width = '100%', warning=FALSE----
fmod <- ~ factor(age) + factor(year) + te(age, year, k = c(5, 15))
fit <- sca(ple4, ple4.indices, fmodel=fmod)
c_d <- computeCatchDiagnostics(fit, ple4)
plot(c_d)


## ----cpred, fig.cap="Prediction of aggregated catch in weight"----------------
plot(c_d, type="prediction", probs=c(0.025, 0.975))


## -----------------------------------------------------------------------------
fit <- sca(ple4, ple4.indices, srmodel=~bevholt(CV=0.2))
fitSumm(fit)


## -----------------------------------------------------------------------------
AIC(fit)
BIC(fit)


## ----echo=FALSE, eval=FALSE---------------------------------------------------
# library(a4adiags)
# theme_set(theme_bw())
# fit <- sca(mut09, mut09.idx, fmod = ~factor(age) + s(year, k = 8))
# res <- residuals(fit, mut09, mut09.idx)


## ----idxrunstest, fig.cap="Runstest for the abundance index", echo=FALSE, eval=FALSE----
# plotRunstest(fit, mut09.idx, combine = F) + theme_bw() + facet_wrap(~age)


## ----catchrunstest, fig.cap="Runstest for the catch by age", echo=FALSE, eval=FALSE----
# plotRunstest(catch.n(mut09), catch.n(mut09 + fit), combine = F) + theme_bw() + facet_wrap(~age)


## -----------------------------------------------------------------------------
fit0 <- sca(ple4, ple4.indices)
n <- 5
nret <- as.list(1:n)
stks <- FLStocks(lapply(nret, function(x){window(ple4, end=(range(ple4)["maxyear"]-x))}))
idxs <- lapply(nret, function(x){window(ple4.indices, end=(range(ple4)["maxyear"]-x))})
fits <- scas(stks, idxs, fmodel=list(fmodel(fit0)))
stks <- stks + fits
stks[[6]] <- ple4 + simulate(fit0, 250)


## ----echo=FALSE---------------------------------------------------------------
lapply(fits, fmodel)


## ----retro, fig.cap="Retrospective analysis of the plaice in ICES area IV stock. Fixed F model."----
plot(window(stks, start=2005))


## -----------------------------------------------------------------------------
n <- 5
nret <- as.list(1:n)
stks <- FLStocks(lapply(nret, function(x){window(ple4, end=(range(ple4)["maxyear"]-x))}))
idxs <- lapply(nret, function(x){window(ple4.indices, end=(range(ple4)["maxyear"]-x))})
# each model will have smootheness scaled to length of time series
fmod <- lapply(stks, defaultFmod)
fits <- scas(stks, idxs, fmodel=fmod)
stks <- stks + fits
stks[[6]] <- ple4 + simulate(fit0, 250)


## ----echo=FALSE---------------------------------------------------------------
lapply(fits, fmodel)


## ----retro2, fig.cap="Retrospective analysis of the plaice in ICES area IV stock. Updating F model."----
plot(window(stks, start=2005))


## ----warning=FALSE, message=FALSE---------------------------------------------
library(a4adiags)
theme_set(theme_bw())
nyears <- 5
# set number of year for average biology and selectivity
nsq <- 3
hc <- a4ahcxval(ple4, ple4.indices, nyears = nyears, nsq = nsq)


## ----hc, echo=FALSE, fig.cap="Survey predictions of year ahead indices in hindcast process. The MASE is presented in the strip about the index and is related to the predictive skill of the index.", warning=FALSE, message=FALSE----
plotXval2(hc$indices)

