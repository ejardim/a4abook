## ----warning=FALSE, message=FALSE---------------------------------------------
library(FLa4a)
library(FLasher)
library(FLBRP)
library(ggplotFL)
data(ple4)
data(ple4.indices)


## ----warning=FALSE, message=FALSE---------------------------------------------
# fit model
fit00 <- sca(ple4, ple4.indices, srmodel=~geomean()) 
stk00 <- ple4 + fit00
# create stock recruitment model object
sr00 <- as(fit00, "FLSR")
# set projection
# number of iterations
nsim <- 250
# most recent year in the data
maxy <- range(ple4)["maxyear"]
# number of years to project
projy <- 5
# last year for projections
endpy <- maxy + projy
# initial year for projections
inipy <- maxy + 1
# extend stock object to store projection's results
stk00 <- fwdWindow(stk00, end = endpy)
# set the control for the projections, in this case a fixed f of 0.3
trg00 <- fwdControl(year = inipy:endpy, quant = "f", value = 0.3)
# project
stk01 <- fwd(stk00, control=trg00, sr=sr00)


## ----echo=FALSE, fig.cap="Projection of stock for 5 years with fixed fishing mortality and recruitment"----
plot(stk01)


## ----warning=FALSE, message=FALSE---------------------------------------------
stk00 <- ple4 + simulate(fit00, nsim)
stk00 <- fwdWindow(stk00, end = endpy)
res00 <- residuals(sr00)
rec00 <- window(rec(stk00), 2018, 2022)
rec00 <- rlnorm(rec00, mean(res00), sqrt(var(res00)))
stk02 <- fwd(stk00, control=trg00, sr=sr00, residuals=rec00)


## ----echo=FALSE, fig.cap="Stochastic projection of stock for 5 years with fixed fishing mortality and recruitment"----
plot(stk02)


## ----warning=FALSE, message=FALSE---------------------------------------------
fit00 <- sca(ple4, ple4.indices) 
stk00 <- ple4 + simulate(fit00, nsim)
sr00 <- as.FLSR(stk00, model="geomean")
sr00 <- fmle(sr00, control = list(trace = 0))
stk00 <- fwdWindow(stk00, end = endpy)
res00 <- residuals(sr00)
rec00 <- window(rec(stk00), inipy, endpy)
rec00 <- rlnorm(rec00, c(yearMeans(res00)), sqrt(c(yearVars(res00))))
stk03 <- fwd(stk00, control=trg00, sr=sr00, residuals=rec00)


## ----echo=FALSE, fig.cap="Stochastic projection of stock for 5 years with fixed fishing mortality and recruitment"----
plot(stk03)


## ----echo=FALSE, fig.cap="Stochastic projection of stock for 5 years with fixed fishing mortality and recruitment. 01: projection without uncertainty, stock recruitment model fit within the stock assessment model; 02: projection with uncertainty, stock recruitment model fit within the stock assessment model; 03: projection with uncertainty, stock recruitment model fit after the stock assessment model"----
plot(window(FLStocks('01' = stk01, '02' = stk02, '03' = stk03), start = 2000))


## ----warning=FALSE, message=FALSE---------------------------------------------
fit00 <- sca(ple4, ple4.indices, srmodel=~geomean()) 
sr00 <- as(fit00, "FLSR")
stk00 <- ple4 + fit00
stk00 <- fwdWindow(stk00, end = endpy, years = list(wt = 20, catch.sel = 10, discards.ratio = 10), fun = list(wt = "sample"))
trg00 <- fwdControl(year = inipy:endpy, quant = "f", value = 0.3)
stk04 <- fwd(stk00, control=trg00, sr=sr00)


## ----echo=FALSE, fig.cap="Stochastic projections of stock for 5 years with fixed fishing mortality and recruitment. Two scenarios with different assumptions about initial conditions"----
plot(FLStocks(default=stk01, alternative=stk04))


## ----warning=FALSE, message=FALSE---------------------------------------------
fit00 <- sca(ple4, ple4.indices, srmodel=~geomean())
sr00 <- as(fit00, "FLSR")
stk00 <- ple4 + simulate(fit00, nsim)
# set projection
projy <- 25
endpy <- maxy + projy
inipy <- maxy + 1
stk00 <- fwdWindow(stk00, end = endpy)
trg00 <- fwdControl(year = inipy:endpy, quant = "f", value = 0)
# recruitment uncertainty
res00 <- residuals(sr00)
rec00 <- window(rec(stk00), inipy, endpy)
rec00 <- rlnorm(rec00, mean(res00), sqrt(var(res00)))
# project
stk05 <- fwd(stk00, control=trg00, sr=sr00, residuals=rec00)


## ----echo=FALSE, fig.cap="Stochastic projection of stock for 25 years in the absence of fishing"----
plot(stk05)


## ----warning=FALSE, message=FALSE---------------------------------------------
trg00 <- fwdControl(year = inipy:endpy, quant = c(rep("ssb_end", 15), rep("f", 10)), value = c(rep(2000000, 15), rep(0.3, 10)))
stk06 <- fwd(stk00, control=trg00, sr=sr00, residuals=rec00)


## ----echo=FALSE, fig.cap="Stochastic projection of stock for 25 years with fixed SSB for 15 years followed by fixed fishing mortality for 10 years and constant recruitment"----
plot(stk06)


## ----warning=FALSE, message=FALSE---------------------------------------------
fit00 <- sca(ple4, ple4.indices, srmodel=~geomean())
sr00 <- as(fit00, "FLSR")
stk00 <- ple4 + simulate(fit00, nsim)
# set projection
projy <- 5
endpy <- maxy + projy
inipy <- maxy + 1
stk00 <- fwdWindow(stk00, end = endpy)
trg00 <- fwdControl(year = inipy:endpy, quant = "ssb_end", value = 1.1, relYear = inipy:endpy-1)
# recruitment uncertainty
res00 <- residuals(sr00)
rec00 <- window(rec(stk00), inipy, endpy)
rec00 <- rlnorm(rec00, mean(res00), sqrt(var(res00)))
# project
stk07 <- fwd(stk00, control=trg00, sr=sr00, residuals=rec00)


## ----warning=FALSE, message=FALSE---------------------------------------------
trg00 <- fwdControl(year = inipy:endpy, quant = "ssb_end", value = 1.1, relYear = maxy)
stk08 <- fwd(stk00, control=trg00, sr=sr00, residuals=rec00)


## ----echo=FALSE, fig.cap="Stochastic projection of stock for 25 years with fixed SSB for 15 years followed by fixed fishing mortality for 10 years and constant recruitment"----
plot(window(FLStocks("10% SSB growth relative to previous year" = stk07, "10% higher SSB relative to most recent estimate" = stk08), start = 2000))


## ----warning=FALSE, message=FALSE---------------------------------------------
minc <- 0.2*mean(catch(stk00), na.rm=TRUE)
trg00 <- fwdControl(year = inipy:endpy, quant = rep(c("ssb_end", "catch"), projy), value = rep(c(1500000, NA), projy), min=rep(c(NA, minc), projy))
stk09 <- fwd(stk00, control=trg00, sr=sr00, residuals=rec00)


## ----echo=FALSE, fig.cap="Stochastic projection of stock for 25 years with SSB target of 1500000t and catch limit of 50% historical catches"----
plot(window(stk09, start = 2000))


## ----warning=FALSE, message=FALSE---------------------------------------------
# fit model
fit00 <- sca(ple4, ple4.indices, srmodel=~geomean())
stk00 <- ple4 + fit00
# create stock recruitment model object
sr00 <- as(fit00, "FLSR")
# estimate reference points
brp00 <- FLBRP(stk00, sr=sr00)
brp00 <- brp(brp00)
ftrg <- refpts(brp00)['msy','harvest']
blim <- 0.5*refpts(brp00)['msy','biomass']
# set projection
# most recent year in the data
maxy <- range(ple4)["maxyear"]
# number of years to project
projy <- 2
# last year for projections
endpy <- maxy + projy
# initial year for projections
inipy <- maxy + 1
# extend stock object to store projection's results
stk00 <- fwdWindow(stk00, end = endpy)
# set the controls for the projections
trg00 <- fwdControl(year = inipy:endpy, quant = "f", value = ftrg)
trg01 <- fwdControl(year = inipy:endpy, quant = "f", value = 0.8*ftrg)
# project
if(ssb(stk00)[,ac(maxy)] > blim){
        stk10 <- fwd(stk00, control=trg00, sr=sr00)
    } else {
        stk10 <- fwd(stk00, control=trg01, sr=sr00)
}


## ----echo=FALSE, fig.cap="Projection of stock for 2 years following a HCR with a target of FMSY and limit of 50% SSBMSY"----
plot(window(stk10, start = 2000))

