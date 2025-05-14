## ----echo=FALSE, message=FALSE, warning=FALSE---------------------------------
library(FLa4a)
library(ggplot2)
data(hke1567)
data(hke1567.idx)


## ----predsim_fit0-------------------------------------------------------------
nsim <- 250
fmod <- ~s(age, k = 4) +
    s(year, k = 8) +
    s(year, k = 8, by = as.numeric(age == 0)) +
    s(year, k = 8, by = as.numeric(age == 4))
qmod <- list(~I(1/(1 + exp(-age))))
fit0 <- sca(hke1567, hke1567.idx, fmodel=fmod, qmodel=qmod)
stk0 <- hke1567 + fit0


## ----predsim_fit_pred---------------------------------------------------------
fit.pred <- predict(fit0)
lapply(fit.pred, names)


## ----echo=FALSE---------------------------------------------------------------
fit.pred$stkmodel$ny1


## ----simny1, fig.cap='Simulations from the model prediction of initial age structure', fig.pos = 'H', fig.height = 4, echo=FALSE, message=FALSE, warning=FALSE----
fit.sim <- simulate(fit0, nsim = 250)
sim.pred <- predict(fit.sim)
fit_sim_ny1 <- sim.pred$stkmodel$ny1

# reduce to quantiles
fit_sim_ny1 <- quantile(fit_sim_ny1, prob = c(0.025, 0.50, 0.975))

# reshape
dat <- reshape(
  as.data.frame(fit_sim_ny1, drop = TRUE),
  timevar = "iter", idvar = c("age"), direction = "wide"
)

# plot
ggplot(data = dat, aes(x = age, y = `data.50%`)) +
  geom_ribbon(aes(ymin = `data.2.5%`, ymax = `data.97.5%`),
    fill = "red", alpha = .15
  ) +
  geom_point() +
  geom_line() +
  ylab("Estimated initial age structure (numbers)") +
  scale_x_continuous(breaks = pretty(dat$age))


## ----predsimhist, fig.cap="Histogram of 250 draws from the approximate distribution of the estimate of survey observation error.", echo=FALSE----
hist(
  exp(coef(fit.sim)$vmodel[[2]]),
  main = "250 draws of a model parameter",
  nclass = 10,
  xlab = "Survey index observation error"
)


## -----------------------------------------------------------------------------
set.seed(1234)
fit.sim1 <- simulate(fit0, nsim = 250)
set.seed(1234)
fit.sim2 <- simulate(fit0, nsim = 250)
all.equal(fit.sim1, fit.sim2)


## ----sim2, fig.cap="Stock summary of the simulated and fitted data"-----------
stk.pred <- hke1567 + fit.sim
plot(FLStocks(simulated=stk.pred, fitted=stk0))


## ----eval=FALSE, echo=FALSE---------------------------------------------------
# nsim <- 250
# fmod <- ~s(age, k = 4) +
#     s(year, k = 8) +
#     s(year, k = 8, by = as.numeric(age == 0)) +
#     s(year, k = 8, by = as.numeric(age == 4))
# qmod <- list(~I(1/(1 + exp(-age))))
# fit0 <- sca(hke1567, hke1567.idx, fmodel=fmod, qmodel=qmod)
# stk0 <- hke1567 + fit0
# 
# # compute deviances
# res0 <- residuals(fit0, hke1567, hke1567.idx, type="deviances")
# # build object with estimation uncertainty
# stk.eu <- hke1567 + simulate(fit0, nsim)
# # build object with residual uncertainty
# stk.ru <- propagate(stk0, nsim) + res0
# # build object with prediction uncertainty
# stk.pu <- hke1567 + simulate(fit0, nsim) + res0


## ----preduncert, fig.cap="Prediction (pu), residual (ru) and estimation (eu) uncertainty", echo=FALSE, eval=FALSE----
# plot(FLStocks(pu=stk.pu, eu=stk.eu, ru=stk.ru)) + facet_grid(qname~stock, scales="free") + theme(legend.position = "top")


## ----message=FALSE, warning=FALSE---------------------------------------------
stk00 <- readRDS("data/MUT1_stk.rds")
idx00 <- readRDS("data/MUT1_idx.rds")


## -----------------------------------------------------------------------------
nits <- 25

shape <- FLModelSim(model=~exp(-age-0.5))
level <- FLModelSim(model=~k^0.66*t^0.57, params = FLPar(k=0.4, t=10),
                     vcov=matrix(c(0.001, 0.01,0.01, 1), ncol=2))
#trend <- FLModelSim(model=~b, params=FLPar(b=0.5), vcov=matrix(0.02))

m4 <- a4aM(shape=shape, level=level)
m4 <- mvrnorm(nits, m4)
range(m4)[] <- range(stk00)[]
range(m4)[c("minmbar","maxmbar")]<-c(1,1)
flq <- m(m4)[]
quant(flq) <- "age"
stk0 <- propagate(stk00, nits)
m(stk0) <- flq


## ----m, fig.cap="Natural mortality generated from M model's parameter uncertainty", echo=FALSE, message=FALSE, warning=FALSE----
bwplot(data~factor(age), data=m(stk0))


## -----------------------------------------------------------------------------
# create objects to store the results
stk01 <- stk0
stk02 <- stk0
stk03 <- propagate(stk00, nits*nits)

# run without estimation uncertainty
stk04 <- stk00 + sca(stk0, idx00)
# update M, the "+" method doesn't do it automatically
m(stk04) <- flq

for(i in 1:nits){
    stk <- iter(stk0, i)
    fit <- sca(stk, idx00)
    # Method 1
    iter(stk01, i) <- stk + simulate(fit, 1)
    # Method 2
    iter(stk02, i) <- qapply(stk + simulate(fit, nits), iterMedians)
    # Method 3
    iter(stk03, (nits*(i-1)+1):(nits*i)) <- stk + simulate(fit, nits)
}



## ----mprop, fig.cap="Stock summary. Stock metrics computed over fits including uncertainty in M and estimation uncertainty", echo=FALSE----
plot(FLStocks("M"=stk04, "M + 1 estimation sample"=stk01, "M + estimation median"=stk02, "M + n estimation samples"=stk03)) + facet_grid(qname~stock, scales="free") + theme(legend.position = "top")

