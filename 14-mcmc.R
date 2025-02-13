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
# check the default pars
mc


## -----------------------------------------------------------------------------
# fit the model
fitmc00 <- sca(hke1567, hke1567.idx, fmodel=fmod, qmodel=qmod, fit = "MCMC", mcmc=mc)
# check acceptance rate
fitSumm(fitmc00)


## -----------------------------------------------------------------------------
plot(hke1567 + fitmc00)


## -----------------------------------------------------------------------------
library(coda)


## -----------------------------------------------------------------------------
# update initial fit
mc <- SCAMCMC(mcmc=100000, mcsave=100)
fitmc01 <- sca(hke1567, hke1567.idx, fmodel=fmod, qmodel=qmod, fit = "MCMC", mcmc=mc)
fitmc01.mc <- FLa4a::as.mcmc(fitmc01)
# highly correlated fit
mc <- SCAMCMC(mcmc=1000, mcsave=1)
fitmc02 <- sca(hke1567, hke1567.idx, fmodel=fmod, qmodel=qmod, fit = "MCMC", mcmc=mc)
fitmc02.mc <- FLa4a::as.mcmc(fitmc02)


## ----chain01, fig.cap="MCMC chains trace. Correlated chain in blue, uncorrelated chain in red."----
traceplot(mcmc.list(mc01=fitmc01.mc[,1], mc02=fitmc02.mc[,2]), lwd=2, col=c(2,4), lty=1)


## ----chain01b, fig.cap="MCMC chain with high autocorrelation after removing the initial 250 samples (burnin period)."----
traceplot(FLa4a::as.mcmc(burnin(fitmc02, 250))[,1], lwd=2, col=4, lty=1)


## ----acf01, fig.cap="Autocorrelation plot of the first parameter in the MCMC chain. Correlated chain in blue, uncorrelated chain in red."----
acfplot(mcmc.list(mc01=fitmc01.mc[,1], mc02=fitmc02.mc[,2]), type="p", pch=19, col=c(2,4))


## ----gew01, fig.cap="Geweke plot of the first parameter in the MCMC chains", fig.show="hold", out.width="50%"----
geweke.plot(fitmc01.mc[,1], main="Uncorrelated chain")
geweke.plot(fitmc02.mc[,1], main="Correlated chain")


## ----cmean01, fig.cap="Cumulative mean plots of the first parameter in the MCMC chains", fig.show="hold", out.width="50%"----
cm01 <- fitmc01.mc[,1]
cm01 <- cumsum(cm01) / seq_along(cm01)
cm02 <- fitmc02.mc[,1]
cm02 <- cumsum(cm02) / seq_along(cm02)
plot(cm01, type="l", xlab="samples", ylab="mean", main="Uncorrelated chain")
plot(cm02, type="l", xlab="samples", ylab="mean", main="Correlated chain")


## ----dens01, fig.cap="Density plots of the first parameter in the MCMC chains", fig.show="hold", out.width="50%"----
densplot(fitmc01.mc[,1], main="Uncorrelated chain")
densplot(fitmc02.mc[,1], main="Correlated chain")


## -----------------------------------------------------------------------------
mc <- SCAMCMC(mcmc=200000, mcsave=200)
fitmc03 <- sca(hke1567, hke1567.idx, fmodel=fmod, qmodel=qmod, fit = "MCMC", mcmc=mc)
fitmc03.mc <- FLa4a::as.mcmc(fitmc03)


## ----chain02, "MCMC chain trace with a thining of 200 samples. No autocorrelation and burnin period apears to exist (left panel). Parameter density shows a well behaved simetric distribution (right panel)."----
plot(fitmc02.mc[,1])


## ----acf02, "Autocorrelation plot of the first parameter in the MCMC chain, with a thining of 200."----
acf(fitmc03.mc[,1])


## -----------------------------------------------------------------------------
mc <- SCAMCMC(mcprobe=0.1)
fitmc01 <- sca(hke1567, hke1567.idx, fmodel=fmod, qmodel=qmod, fit = "MCMC", mcmc=SCAMCMC())
fitSumm(fitmc01)
fitmc01p <- sca(hke1567, hke1567.idx, fmodel=fmod, qmodel=qmod, fit = "MCMC", mcmc=mc)
fitSumm(fitmc01p)
mc <- SCAMCMC(mcmc=1000, mcsave=1, mcprobe=0.0005)
fitmc02p <- sca(hke1567, hke1567.idx, fmodel=fmod, qmodel=qmod, fit = "MCMC", mcmc=mc)
fitSumm(fitmc02p)
mc <- SCAMCMC(mcmc=200000, mcsave=200, mcprobe=0.0005)
fitmc03p <- sca(hke1567, hke1567.idx, fmodel=fmod, qmodel=qmod, fit = "MCMC", mcmc=mc)
fitSumm(fitmc03p)




## -----------------------------------------------------------------------------
plot(FLStocks(ll=hke1567 + fit, mc=hke1567 + fitmc00, mc_alt=hke1567 + fitmc01))


## -----------------------------------------------------------------------------
# update initial fit
mc <- SCAMCMC(mcmc=100000, mcsave=100, mcseed=10)
fitmc01 <- sca(hke1567, hke1567.idx, fmodel=fmod, qmodel=qmod, fit = "MCMC", mcmc=mc)
fitmc01.mc <- FLa4a::as.mcmc(fitmc01)

mc <- SCAMCMC(mcmc=100000, mcsave=100, mcseed=30)
fitmc01b <- sca(hke1567, hke1567.idx, fmodel=fmod, qmodel=qmod, fit = "MCMC", mcmc=mc)
fitmc01b.mc <- FLa4a::as.mcmc(fitmc01b)

# highly correlated fit
mc <- SCAMCMC(mcmc=1000, mcsave=1, mcseed=10)
fitmc02 <- sca(hke1567, hke1567.idx, fmodel=fmod, qmodel=qmod, fit = "MCMC", mcmc=mc)
fitmc02.mc <- FLa4a::as.mcmc(fitmc02)

mc <- SCAMCMC(mcmc=1000, mcsave=1, mcseed=30)
fitmc02b <- sca(hke1567, hke1567.idx, fmodel=fmod, qmodel=qmod, fit = "MCMC", mcmc=mc)
fitmc02b.mc <- FLa4a::as.mcmc(fitmc02b)


