
library(fda)
library(lattice)
library(INLA)

##################################
## get smooth design matrix
##################################

load("B:/Loch_Ard_Felling/Package/rivertemp/data/fullwater.rda")
# a quick fix
stream $ year[stream $ Sensor.Model == "Squirrel"] <- stream $ year[stream $ Sensor.Model == "Squirrel"] + 1900 
# set up a bspline basis of dimension 48 and order 6 
hourbasis <- create.bspline.basis(c(0, 24), 48, 6)
# and a harmonic accelerator penalty... penalises deviations from sinusiodal curves
penalty <- vec2Lfd(c(0, (2*pi/24)^2, 0), c(0, 24))

# select a few days data from control burn and add in decimal hour
dat <- with(subset(stream, site == 2 & year %in% 2003:2003), 
    data.frame(temp = Original.Value, 
               dhour = hour + min/60, 
               month = factor(month.abb[mon + 1], levels = month.abb), 
               yday = yday, 
               year = year, 
               site = site))
dat <- subset(dat, yday < 365) # remove last day of leap year

getcoef <- function(lambda, .yday, .year = 2009, .site = 2) {
  sdat <- subset(dat, yday == .yday & year == .year & site == .site)
  #if (nrow(sdat) < 23) return (NULL)
  if (sum(unique(round(sdat $ dhour)) %in% 0:24) < 24) return(NULL)  
  c(coef(smooth.basis(sdat $ dhour, sdat $ temp, fdPar(hourbasis, penalty, lambda))))
}

yday <- rep(0:365, length(2003:2009))
year <- rep(2003:2009, each = 366)

ck <- mapply(getcoef, .yday = yday, .year = year, MoreArgs = list(lambda = exp(-8)))
whichNULL <- sapply(ck, is.null)
ck <- simplify2array(ck[!whichNULL])

# is this sensitive to missing values?...
pcs <- pca.fd(fd(ck, hourbasis), nharm = 4)


#####################################################


setupStuff <- function(yr, control = 2, felled = 10) {
  # we need a different S matrix by day and for each site
  # and we can only work with paired data for now
  Xdf <- with(subset(stream, site == control & year == yr), {
             data.frame(x = Original.Value, 
                        hr = hour + min/60,
                        yday = yday)})

  Ydf <- with(subset(stream, site == felled & year == yr), {
             data.frame(y = Original.Value, 
                        hr = hour + min/60,
                        yday = yday,
                        year = yr)})

  # which days to keep in the analysis..
  # make sure there are sufficient observations.. one per hour at least?
  days <- intersect(Xdf $ yday, Ydf $ yday)

  Ydf <- subset(Ydf, yday %in% days)
  Xdf <- subset(Xdf, yday %in% days)

  # extract hours by day
  Yhrs <- lapply(days, function(i) Ydf $ hr[Ydf $ yday == i])
  Xhrs <- lapply(days, function(i) Xdf $ hr[Xdf $ yday == i])

  # extract temperature obs by day
  Yobs <- lapply(days, function(i) Ydf $ y[Ydf $ yday == i])
  Xobs <- lapply(days, function(i) Xdf $ x[Xdf $ yday == i])
 
  # get Y summaries by day
  Sy <- lapply(Yhrs, eval.fd, fdobj = pcs $ harmonics)
  muy <- lapply(Yhrs, eval.fd, fdobj = pcs $ meanfd)
  Ybeta <- sapply(seq_along(days), function(i) solve(crossprod(Sy[[i]]), t(Sy[[i]])) %*% (Yobs[[i]] - muy[[i]]))

  # hat matrix to convert X into curve space (using the same functional transformation for now)
  Sx <- lapply(Xhrs, eval.fd, fdobj = pcs $ harmonics)
  mux <- lapply(Xhrs, eval.fd, fdobj = pcs $ meanfd)
  Xbeta <- sapply(seq_along(days), function(i) solve(crossprod(Sx[[i]]), t(Sx[[i]])) %*% (Xobs[[i]] - mux[[i]]))

  list(muy = muy, Ybeta = Ybeta, Xbeta = Xbeta,
       yr = yr, days = days)
}


felleffect <- function(yrs, days, theta = 0) {
  days <- days + (yrs - min(yrs)) * 366
#  ifelse(yrs == min(yrs), 0, 
#    ifelse(yrs == min(yrs) + 1, 1,
#      pmax(1-1/(theta * 366) * (days - 366*2), 0)))
  ifelse(yrs == min(yrs), 0, 
    ifelse(yrs == min(yrs) + 1, 1,
      dnorm(days - 366*2, 0, theta*366) * sqrt(2*pi) * theta * 366
    ))
}
plot(1:(365*10), felleffect(yrs = rep(1:10, each = 365), days = rep(1:365, 10), theta = 2), type = "l")

###############################################################


dff <- function(theta = 4, pc = 1) {
  y <- unlist(sapply(allstuff, function(x) x $ Ybeta[pc,]))
  X <- do.call(rbind, lapply(allstuff, function(x) t(x $ Xbeta)))
  colnames(X) <- paste0("PC", 1:4)
  days <- unlist(sapply(allstuff, function(x) x $ days))
  yr <- unlist(sapply(allstuff, function(x) rep(x $ yr, length(x $ days))))

  df <- data.frame(y = y, days = days, yr = yr,
                   felled = yr > min(yr))
  df <- cbind(df, X)
 
  within(df, {fell = felleffect(yr, days, theta)})
}

# doone <- function(form) {
# # estimate recovery period
# theta <- seq(2, 15, length = 50)
# llik <- sapply(theta, 
#          function(.th) { 
#            logLik(
#              gam(form, data = dff(.th), method = "ML")
#             )
#          }) 
# plot(theta, llik)
# thest <- theta[which.max(aic)]
# abline(v = thest)

# m <- gam(form, data = dff(thest), method = "ML")
# attr(m, "theta") <- thest
# m
# }


form <- function(pc = 1, type = "inla", smoothpc = pc) {
  lpc <- setdiff(1:4, pc)
  if (type == "inla") {
    out <- 
      paste0("y ~ ", 
            "PC", lpc[1], ":prefell + PC", lpc[1], ":fell + ",
            "PC", lpc[2], ":prefell + PC", lpc[2], ":fell + ",
            "PC", lpc[3], ":prefell + PC", lpc[3], ":fell - 1 + ",
            "f(days.fell.PC   , PC.fell   , model = 'rw2', cyclic = TRUE) + ",
            "f(days.prefell.PC, PC.prefell, model = 'rw2', cyclic = TRUE) + ",
            "f(days.fell      , fell   , model = 'rw2', cyclic = TRUE) + ", 
            "f(days.prefell   , prefell, model = 'rw2', cyclic = TRUE) + ",
            "f(days.ar1       , model = 'ar1')")
  } else if (type == "mgcv") {
    out <- 
      paste0("y ~ fell + ", 
            "PC", lpc[1], " + PC", lpc[1], ":fell + ",
            "PC", lpc[2], " + PC", lpc[2], ":fell + ",
            "PC", lpc[3], " + PC", lpc[3], ":fell + ",
            "s(days, by = PC", smoothpc, ".fell, bs = 'cc') + ",
            "s(days, by = PC", smoothpc, ".prefell, bs = 'cc') + ",
            "s(days, by = fell, bs = 'cc') + ", 
            "s(days, by = prefell, bs = 'cc')")
  } 
  as.formula(out)
}
form(2)
form(3, type = "mgcv")


data <- function(pc = 1, theta = 8, smoothpc = pc) {
  # data
  out <- within(dff(theta, pc = pc), {
      prefell <- 1 - fell
      days.fell.PC <- days
      days.prefell.PC <- days
      days.fell <- days
      days.prefell <- days
      prefell <- 1 - fell
      days.ar1 <- days + (yr - min(yr)) * 365
     })
  out[paste0("PC", smoothpc, ".fell")] <- with(out, get(paste0("PC", smoothpc)) * fell)
  out[paste0("PC", smoothpc, ".prefell")] <- with(out, get(paste0("PC", smoothpc)) * (1-fell))
  out
}
head(data(4, 5))

inlafit <- function(form, pc = 1, theta = 8) {
  # data
  dat <- data(pc = pc, theta = theta)

  r <- inla(form, data = dat, control.compute=list(cpo=TRUE))
  r
}

magicfit <- function(pc = 1, theta = 8, rho = 0.9, smoothpc = pc) { 
  ## fit using magic
  dat <- data(pc = pc, theta = theta, smoothpc = smoothpc)

  ## produce scaled covariance matrix for AR1 errors...
  V <- corMatrix(Initialize(corAR1(rho, form = ~ days | yr), dat))
  V <- as.matrix(do.call(bdiag, V))
  Cv <- chol(V)  # t(Cv)%*%Cv=V

  ## Fit smooth, taking account of *known* correlation...
  w <- solve(t(Cv)) # V^{-1} = w'w

  ## Use `gam' to set up model for fitting...
  G <- gam(form(pc = pc, type = "mgcv"), data = dat, fit = FALSE)

  ## fit using magic, with weight *matrix*
  mgfit <- magic(G$y, G$X, G$sp, G$S, G$off, rank = G$rank, C = G$C, w = w)

  ## GAM ignoring correlation
  b <- gam(form(pc = pc, type = "mgcv"), data = dat)
  ## Modify previous gam object using new fit, for plotting...    
  mg.stuff <- magic.post.proc(G $ X, mgfit, w)
  b $ edf[] <- mg.stuff $ edf
  b $ Vp[] <- mg.stuff $ Vb
  b $ coefficients[] <- mgfit $ b 
  b
}

##################################################################################

 
# try some models

with(subset(stream, yday < 45), table(year, yday, site))

# do all stuff once
#yrs <- c(1988, 1990, 1992:1995)
yrs <- c(1988, 1990, 1992:1995)
# need to run this bit as it is used in the setup function
allstuff <- lapply(yrs, setupStuff, control = 11, felled = 10)


# the theta values to try
thetas <- seq(2, 10, length = 30) # should be an hour per model

# PC 1
# ######

pc <- 1
system.time(
i1s <-
  lapply(thetas, function(th.) {
    inlafit(form(pc), pc = pc, theta = th.)
  })
)

# save stuff
save(i1s, file = "i1s.rdata")


# PC 2
# ######

pc <- 2
system.time(
i2s <-
  lapply(thetas, function(th.) {
    inlafit(form(pc), pc = pc, theta = th.)
  })
)
# save stuff
save(i2s, file = "i2s.rdata")


# PC 3
# ######

pc <- 3
system.time(
i3s <-
  lapply(thetas, function(th.) {
    inlafit(form(pc), pc = pc, theta = th.)
  })
)
# save stuff
save(i3s, file = "i3s.rdata")


# PC 4
# ######

pc <- 4
system.time(
i4s <-
  lapply(thetas, function(th.) {
    inlafit(form(pc), pc = pc, theta = th.)
  })
)
# save stuff
save(i4s, file = "i4s.rdata")



par(mfrow = c(2,2))
LS1 <- sapply(i1s, function(x) ifelse(sum(x $ cpo $ fail)>20, NA, mean(-log(x $ cpo $ cpo[x $ cpo $ fail == 0]))))
plot(thetas, LS1)

LS2 <- sapply(i2s, function(x) ifelse(sum(x $ cpo $ fail)>20, NA, mean(-log(x $ cpo $ cpo[x $ cpo $ fail == 0]))))
plot(thetas, LS2)

LS3 <- sapply(i3s, function(x) ifelse(sum(x $ cpo $ fail)>20, NA, mean(-log(x $ cpo $ cpo[x $ cpo $ fail == 0]))))
plot(thetas, LS3)

LS4 <- sapply(i4s, function(x) ifelse(sum(x $ cpo $ fail)>20, NA, mean(-log(x $ cpo $ cpo[x $ cpo $ fail == 0]))))
plot(thetas, LS4)


# extract ar1 parameter estimates and look at how the change with presistence variable
rhos1 <- sapply(i1s, function(i) i $ summary.hyperpar[7,])
rhos2 <- sapply(i2s, function(i) i $ summary.hyperpar[7,])
rhos3 <- sapply(i3s, function(i) i $ summary.hyperpar[7,])
rhos4 <- sapply(i4s, function(i) i $ summary.hyperpar[7,])

# plot
par(mfrow = c(2,2))
plot(LS1, rhos1[1,])
plot(LS2, rhos2[1,])
plot(LS3, rhos3[1,])
plot(LS4, rhos4[1,])

rho1 <- rhos1[which.min(LS1)]
rho2 <- rhos2[which.min(LS2)]
rho3 <- rhos3[which.min(LS3)]
rho4 <- rhos4[which.min(LS4)]

# test 
b <- magicfit(pc = 1, theta = 8, rho = 0.9)
plot(b, pages = 1, scale = 0)
summary(b)

# the theta values to try
thetas <- seq(0.1, 10, length = 30) # should be an hour per model

# Fit final model to PC1
system.time(
b1s <- lapply(thetas, 
  function(th.) magicfit(pc = 1, theta = th., rho = rho1))
)

# Fit final model to PC2
system.time(
b2s <- lapply(thetas, 
  function(th.) magicfit(pc = 2, theta = th., rho = rho2))
)

# Fit final model to PC3
system.time(
b3s <- lapply(thetas, 
  function(th.) magicfit(pc = 3, theta = th., rho = rho3))
)

# Fit final model to PC4
thetas4 <- seq(1.5, 10, length = 30)
system.time(
b4s <- lapply(thetas4, 
  function(th.) magicfit(pc = 4, theta = th., rho = rho4))
)

# plot up theta estimates
par(mfrow = c(2,2))

ll1s <- sapply(b1s, logLik)
plot(thetas, ll1s, type = "o", main = "PC1")

ll2s <- sapply(b2s, logLik)
plot(thetas, ll2s, type = "o", main = "PC2")

ll3s <- sapply(b3s, logLik)
plot(thetas, ll3s, type = "o", main = "PC3")

ll4s <- sapply(b4s, logLik)
plot(thetas4, ll4s, type = "o", main = "PC4")

# now collect best models together...
theta1 <- thetas[which.max(ll1s)]
theta2 <- thetas[which.max(ll2s)]
theta3 <- thetas[which.max(ll3s)]
theta4 <- thetas4[which.max(ll4s)]

b1 <- magicfit(pc = 1, theta = theta1, rho = rho1)
b2 <- magicfit(pc = 2, theta = theta2, rho = rho2)
b3 <- magicfit(pc = 3, theta = theta3, rho = rho3)
b4 <- magicfit(pc = 4, theta = theta4, rho = rho4)

# check them out for sanity
par(mfrow = c(4,4), oma = c(0,0,0,0), mar = c(3, 2, 2, 1) + .1)
summary(b1)
plot(b1, scale = 0)

summary(b2)
plot(b2, scale = 0)

summary(b3)
plot(b3, scale = 0)

summary(b4)
plot(b4, scale = 0)

# use fitted values to build the predicted PC scores
Bpred <- cbind(fitted(b1), fitted(b2), fitted(b3), fitted(b4))

nhrs <- 24
hrs <- seq(0, 24, length = nhrs)
Sy <- eval.fd(hrs, fdobj = pcs $ harmonics)
muy <- eval.fd(hrs, fdobj = pcs $ meanfd)

Ypred <- Sy %*% t(Bpred) + c(muy)
#matplot(Ypred, type = "l")

df <- data(theta = 3)
preddf <- data.frame(hrs = hrs, pred = c(Ypred),
            days = rep(df $ days, each = nhrs),
            yr = rep(df $ yr, each = nhrs),
            fell = rep(df $ fell, each = nhrs))
preddf $ ddays <- with(preddf, days + hrs / 24)
preddf $ col <- colorRampPalette(c("blue", "red"))(101)[round(preddf $ fell, 2) * 100 + 1]

#xyplot(pred ~ I(ddays + yr*365), data = preddf, type = "p", cex = 0.5, col = preddf $ col, pch = 16)

# so now some plots of the effect of felling

## lets choose a mean day from each month
df $ imonth <- as.numeric(cut(df $ days, cumsum(c(-1, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 32))))
df $ month <- factor(month.abb[df $ imonth], month.abb)

meanPCs <- function(felled) {

  out <- with(df, cbind.data.frame( 
                    PC1 = c(tapply(PC1, list(month, yr), mean)),
                    PC2 = c(tapply(PC2, list(month, yr), mean)),
                    PC3 = c(tapply(PC3, list(month, yr), mean)),
                    PC4 = c(tapply(PC4, list(month, yr), mean))))

  mdays <- c(31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  out <-
    within(out, {
      days <- cumsum(mdays) - mdays/2
      fell <- felled
      prefell <- 1 - fell
      PC1.fell <- PC1 * fell
      PC1.prefell <- PC1 * prefell
      PC2.fell <- PC2 * fell
      PC2.prefell <- PC2 * prefell
      PC3.fell <- PC3 * fell
      PC3.prefell <- PC3 * prefell
      PC4.fell <- PC4 * fell
      PC4.prefell <- PC4 * prefell
    })
  out
}

felleffect <- with(df, c(tapply(fell, list(month, yr), mean)))

felledPCs <- sapply(list(b1, b2, b3, b4), predict, newdata = meanPCs(felled = felleffect))
felledCurv <- Sy %*% t(felledPCs) + c(muy)

prefellPCs <- sapply(list(b1, b2, b3, b4), predict, newdata = meanPCs(felled = 0))
prefellCurv <- Sy %*% t(prefellPCs) + c(muy)

# confidence intervals for the felling effect
Curvsim <- function(n) {
  PCs <-
  sapply(list(b1, b2, b3, b4), function(x) {
    bsim <- mvrnorm(n, coef(x), x $ Vp)
    Xf <- predict(x, newdata = meanPCs(felled = felleffect), type = "lpmatrix")
    Xb <- predict(x, newdata = meanPCs(felled = 0), type = "lpmatrix")

    (Xf - Xb) %*% t(bsim)
  })

  curves <- Sy %*% t(PCs)
  dim(curves) <- c(24, nrow(meanPCs(felled = 0)), n)
  curves
}

curvsim <- Curvsim(1000)

felledCurv <- prefellCurv + apply(curvsim, 1:2, median)
felledCurv.cil <- prefellCurv + apply(curvsim, 1:2, quantile, 0.025, na.rm = TRUE)
felledCurv.ciu <- prefellCurv + apply(curvsim, 1:2, quantile, 0.975, na.rm = TRUE)

preddf <- data.frame(hrs = hrs, felled = c(felledCurv), 
            fell.ciu =  c(felledCurv.ciu), fell.cil = c(felledCurv.cil),
            prefell = c(prefellCurv),
            month = factor(rep(month.abb, each = nhrs), month.abb),
            year = rep(unique(df$yr), each = nhrs*12))

# ta da!
xyplot(prefell + felled + fell.cil + fell.ciu ~ hrs | month + factor(year), 
       data = preddf, type = "l", as.table = TRUE, lty = c(1, 1, 2, 2), col = c(1, 2, 2, 2))

## now to demonstrate the effect on an extreme day rather than the mean day

Curvsim2 <- function(n) {
  PCs <-
  sapply(list(b1, b2, b3, b4), function(x) {
    bsim <- mvrnorm(n, coef(x), x $ Vp)
    Xf <- predict(x, type = "lpmatrix")
    Xb <- Xf
    Xb. <- predict(x, newdata = meanPCs(felled = 0), type = "lpmatrix")
    Xb[,colSums(Xb., na.rm = TRUE) == 0] <- 0 

    (Xf - Xb) %*% t(bsim)
  })

  curves <- Sy %*% t(PCs)
  dim(curves) <- c(24, nrow(model.frame(b1)), n)
  curves
}

curvsim <- Curvsim2(1000)

maxCurvs <- sapply(1:1000, function(i) t(apply(curvsim[,,i], 1, function(x) tapply(x, list(df $ month, df $ yr), max))))
dim(maxCurvs) <- c(24, 72, 1000)

felledCurv <- prefellCurv + apply(maxCurvs, 1:2, median)
felledCurv.cil <- prefellCurv + apply(maxCurvs, 1:2, quantile, 0.025, na.rm = TRUE)
felledCurv.ciu <- prefellCurv + apply(maxCurvs, 1:2, quantile, 0.975, na.rm = TRUE)

preddf <- data.frame(hrs = hrs, felled = c(felledCurv), 
            fell.ciu =  c(felledCurv.ciu), fell.cil = c(felledCurv.cil),
            prefell = c(prefellCurv),
            month = factor(rep(month.abb, each = nhrs), month.abb),
            year = rep(unique(df$yr), each = nhrs*12))

# ta da!
xyplot(prefell + felled + fell.cil + fell.ciu ~ hrs | month + factor(year), 
       data = preddf, type = "l", as.table = TRUE, lty = c(1, 1, 2, 2), col = c(1, 2, 2, 2))

## and minimum effects

minCurvs <- sapply(1:1000, function(i) t(apply(curvsim[,,i], 1, function(x) tapply(x, list(df $ month, df $ yr), min))))
dim(minCurvs) <- c(24, 72, 1000)

felledCurv <- prefellCurv + apply(minCurvs, 1:2, median)
felledCurv.cil <- prefellCurv + apply(minCurvs, 1:2, quantile, 0.025, na.rm = TRUE)
felledCurv.ciu <- prefellCurv + apply(minCurvs, 1:2, quantile, 0.975, na.rm = TRUE)

preddf <- data.frame(hrs = hrs, felled = c(felledCurv), 
            fell.ciu =  c(felledCurv.ciu), fell.cil = c(felledCurv.cil),
            prefell = c(prefellCurv),
            month = factor(rep(month.abb, each = nhrs), month.abb),
            year = rep(unique(df$yr), each = nhrs*12))

# ta da!
xyplot(prefell + felled + fell.cil + fell.ciu ~ hrs | month + factor(year), 
       data = preddf, type = "l", as.table = TRUE, lty = c(1, 1, 2, 2), col = c(1, 2, 2, 2))

# So that is one felling event finished!  At least in draft.
# now for the other felling events....
