## ----echo=FALSE, message=FALSE, warning=FALSE---------------------------------
library(knitr)
opts_chunk$set(dev='png', dev.args=list(type="cairo"))


## -----------------------------------------------------------------------------
# load libraries and data
library(FLa4a)
library(ggplotFL)
data(hke1567)
data(hke1567.idx)
nsim <- 250
# MLE estimate
fmod <- ~s(age, k = 4) +
    s(year, k = 8) +
    s(year, k = 8, by = as.numeric(age == 0)) +
    s(year, k = 8, by = as.numeric(age == 4))
qmod <- list(~I(1/(1 + exp(-age))))
fit <- sca(hke1567, hke1567.idx, fmodel=fmod, qmodel=qmod)
fit <- simulate(fit, nsim)


## -----------------------------------------------------------------------------
# mcmc
mc <- SCAMCMC()


## -----------------------------------------------------------------------------
# fit the model
fitmc00 <- sca(hke1567, hke1567.idx, fmodel=fmod, qmodel=qmod, fit = "MCMC", mcmc=mc)
# check acceptance rate
fitSumm(fitmc00)


## ----fig.cap="Stock assessment summaries of maximum likelihood (mle) and monte carlo (mc) fits."----
plot(FLStocks(mle=hke1567 + fit, mc=hke1567 + fitmc00))


## -----------------------------------------------------------------------------
library(coda)


## -----------------------------------------------------------------------------
# update initial fit, control random seed
mc <- SCAMCMC(mcmc=100000, mcsave=100, mcseed=10)
fitmc01 <- sca(hke1567, hke1567.idx, fmodel=fmod, qmodel=qmod, fit = "MCMC", mcmc=mc)
fitmc01.mc <- FLa4a::as.mcmc(fitmc01)
# highly correlated fit, control random seed
mc <- SCAMCMC(mcmc=1000, mcsave=1, mcseed=10)
fitmc02 <- sca(hke1567, hke1567.idx, fmodel=fmod, qmodel=qmod, fit = "MCMC", mcmc=mc)
fitmc02.mc <- FLa4a::as.mcmc(fitmc02)


## ----chain01, fig.cap="MCMC chain's trace for the first parameter. High correlation chain in blue, low correlation chain in red."----
# use mcmc.list() to create a list with both runs so they're plot together
traceplot(mcmc.list(mc01=fitmc01.mc[,1], mc02=fitmc02.mc[,2]), lwd=1.5, col=c(2,4), lty=1)


## ----chain01b, fig.cap="MCMC chain with high autocorrelation after removing the initial 250 samples (burnin period)."----
traceplot(FLa4a::as.mcmc(burnin(fitmc02, 250))[,1], lwd=1.5, col=4, lty=1)


## ----acf01hc, fig.cap="Autocorrelation plot of the first parameter in the MCMC chain", fig.show="hold"----
acfplot(fitmc02.mc[,1], lwd=3, main="High correlation chain", ylim=c(-1, 1))


## ----acf01, fig.cap="Autocorrelation plot of the first parameter in the MCMC chain", fig.show="hold"----
acfplot(fitmc01.mc[,1], lwd=3, main="Low correlation chain", ylim=c(-1, 1))


## ----ccr01, fig.cap="Crosscorrelation plots", fig.show="hold", out.width="50%"----
crosscorr.plot(fitmc01.mc, main="Low correlation chain")
crosscorr.plot(fitmc02.mc, main="High correlation chain")


## ----gew01, fig.cap="Geweke plot of the first parameter in the MCMC chains", fig.show="hold", out.width="50%"----
geweke.plot(fitmc01.mc[,1], main="Low correlation chain")
geweke.plot(fitmc02.mc[,1], main="High correlation chain")


## ----cmean01, fig.cap="Cumulative mean plots of the first parameter in the MCMC chains", fig.show="hold", out.width="50%"----
cm01 <- fitmc01.mc[,1]
cm01 <- cumsum(cm01) / seq_along(cm01)
cm02 <- fitmc02.mc[,1]
cm02 <- cumsum(cm02) / seq_along(cm02)
plot(cm01, type="l", xlab="samples", ylab="mean", main="Low correlation chain", ylim=c(-0.52, -0.42))
plot(cm02, type="l", xlab="samples", ylab="mean", main="High correlation chain", ylim=c(-0.52, -0.42))


## ----dens01, fig.cap="Density plots of the first parameter in the MCMC chains", fig.show="hold", out.width="50%"----
densplot(fitmc01.mc[,1], main="Low correlation chain")
densplot(fitmc02.mc[,1], main="High correlation chain")


## -----------------------------------------------------------------------------
# low correlation
mc <- SCAMCMC(mcmc=100000, mcsave=100, mcseed=30)
fitmc01b <- sca(hke1567, hke1567.idx, fmodel=fmod, qmodel=qmod, fit = "MCMC", mcmc=mc)
fitmc01b.mc <- FLa4a::as.mcmc(fitmc01b)
# highly correlated fit
mc <- SCAMCMC(mcmc=1000, mcsave=1, mcseed=30)
fitmc02b <- sca(hke1567, hke1567.idx, fmodel=fmod, qmodel=qmod, fit = "MCMC", mcmc=mc)
fitmc02b.mc <- FLa4a::as.mcmc(fitmc02b)
# create lists for comparison
mclst01 <- mcmc.list(a=fitmc01.mc, b=fitmc01b.mc)
mclst02 <- mcmc.list(a=fitmc02.mc, b=fitmc02b.mc)


## ----echo=FALSE---------------------------------------------------------------
gelman.diag(mclst01)


## ----echo=FALSE---------------------------------------------------------------
gelman.diag(mclst02)


## ----gelm01, fig.cap="Gelman-Rubin's diagnostic plots for the first parameter.", fig.show="hold", out.width="50%"----
mclst01 <- mcmc.list(a=fitmc01.mc[,1], b=fitmc01b.mc[,1])
mclst02 <- mcmc.list(a=fitmc02.mc[,1], b=fitmc02b.mc[,1])
gelman.plot(mclst01, main="Low correlation chain", ylim=c(1, 5))
gelman.plot(mclst02, main="High correlation chain", ylim=c(1, 5))


## -----------------------------------------------------------------------------
cbind(fitSumm(fitmc01), fitSumm(fitmc02))


## -----------------------------------------------------------------------------
# no scale
mc <- SCAMCMC(mcmc=100000, mcsave=100, mcseed=10, mcnoscale=TRUE)
fitmc01ns <- sca(hke1567, hke1567.idx, fmodel=fmod, qmodel=qmod, fit = "MCMC", mcmc=mc)
fitmc01ns.mc <- FLa4a::as.mcmc(fitmc01ns)
data.frame(hessian_scale=c(fitSumm(fitmc01)), hessian_noscale=c(fitSumm(fitmc01ns)))


## -----------------------------------------------------------------------------
data.frame(hessian_scale=c(autocorr.diag(fitmc01.mc[,1])), hessian_noscale=c(autocorr.diag(fitmc01ns.mc[,1])))


## -----------------------------------------------------------------------------
# more Cauchy
mc <- SCAMCMC(mcmc=1000, mcsave=1, mcseed=10, mcprobe=0.45)
fitmc02p <- sca(hke1567, hke1567.idx, fmodel=fmod, qmodel=qmod, fit = "MCMC", mcmc=mc)
fitmc02p.mc <- FLa4a::as.mcmc(fitmc02p)
data.frame(probe_0.05=c(fitSumm(fitmc02)), probe_0.45=c(fitSumm(fitmc02p)))


## -----------------------------------------------------------------------------
data.frame(probe_0.05=c(autocorr.diag(fitmc02.mc[,1])), probe_0.45=c(autocorr.diag(fitmc02p.mc[,1])))


## -----------------------------------------------------------------------------
# reduce correlation
mc <- SCAMCMC(mcmc=10000, mcsave=10, mcseed=10, mcrb=1)
fitmc02r <- sca(hke1567, hke1567.idx, fmodel=fmod, qmodel=qmod, fit = "MCMC", mcmc=mc)
fitmc02r.mc <- FLa4a::as.mcmc(fitmc02r)
data.frame(rslbound_no=c(fitSumm(fitmc02)), rslbound_high=c(fitSumm(fitmc02r)))


## -----------------------------------------------------------------------------
data.frame(rslbound_no=c(autocorr.diag(fitmc02.mc[,1])), rslbound_high=c(autocorr.diag(fitmc02r.mc[,1])))

