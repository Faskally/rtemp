
if (0) {
setwd("B:/Loch_Ard_Felling/")

# attach functions as if it were a package
local({
  # load required libraries
  library(lattice)
  library(fda)
  library(reshape)

  # remove old instances
  sapply(rev(grep("rivertempFunctions", search())), function(i) detach(pos = i))

  # attach a named empty environment and source in functions to it
  env <- attach(NULL, name = "rivertempFunctions")
  files <- dir("Package/rivertemp/R/", full = TRUE)
  lapply(files, sys.source, envir = env)
  cat("This is rivertemp version 1.0\n")
})

search()
ls(2)

# get full data set
load("Package/rivertemp/data/fullwater.rda")
stream $ year[stream $ Sensor.Model == "Squirrel"] <- stream $ year[stream $ Sensor.Model == "Squirrel"] + 1900 
head(stream)

# select a few days data
dat <- subset(stream, site == 2 & year == 2009 & yday %in% (100 + 0:10))
rownames(dat) <- NULL
dat <- within(dat, {
                val <- Original.Value
                dhour <- hour + min/60
            })[c("dhour","val")]

hour <- unique(dat $ dhour)
tempmat <- matrix(dat $ val, nrow = length(hour))
matplot(tempmat, type = "l")

# nbasis must be odd for fourier bases
tempbasis <- create.bspline.basis(c(0, 24), 96, 6)

# derivative penalties
harmaccelLfd = vec2Lfd(c(0, (2*pi/24)^2, 0), c(0, 24))

# get an appropriate penalty:
loglambdas <- seq(-12, -4, length = 50)
gcvs <- 
  sapply(loglambdas, 
    function(i) 
      sum(
        smooth.basis(hour, 
                 tempmat, 
                 fdPar(tempbasis, harmaccelLfd, exp(i))) $ gcv
        ))
# try later: lambda2gcv?
plot(loglambdas, gcvs, type = "l")

loglambdas[which.min(gcvs)]
loglambda <- 0
tempbasispen <- fdPar(tempbasis, harmaccelLfd, exp(loglambda))

tempSmooth <- smooth.basis(hour, tempmat, tempbasispen)
tempfd <- tempSmooth $ fd

plot(tempfd)

hourfine <- seq(0, 24, length = 1000)
velffine <- eval.fd(hourfine, tempfd, 1)
accffine <- eval.fd(hourfine, tempfd, 2)
hourpt <- seq(0, 24, length = 96)
velpt <- eval.fd(hourpt, tempfd, 1)
accpt <- eval.fd(hourpt, tempfd, 2)
par(mfrow = c(3,4))
for (i in 1:11) {
  plot(velffine[,i], accffine[,i], type = "l", 
       xlab = "Velocity", ylab = "Acceleration", main = paste("day", i))
  points(velpt[,i], accpt[,i], pch = 21, cex = 0.8, col = 1, bg = rainbow(96))
  abline(h = 0, v = 0, lty = 3)
}
plot(tempfd)

# Landmarkregistration
accelfdUN <- tempfd
accelmeanfdUN <- mean(accelfdUN)

MINs <- sapply(1:11, function(i) seq(0, 24, length = 100)[which.max(predict(accelfdUN[i], seq(0, 24, length = 100)))])
# see page 123 or there abouts in FDA in R
wbasisLM <- create.bspline.basis(c(0,24), 4, 3, c(0, mean(MINs),24))
WfdLM <- fd(matrix(0, 4, 1), wbasisLM)
WfdParLM <- fdPar(WfdLM, 1, 1e-12)

regListLM <- landmarkreg(accelfdUN, MINs, mean(MINs), WfdParLM, TRUE)
accelfdLM <- regListLM $ regfd
accelmeanfdLM <- mean(accelfdLM)
warpfdLM <- regListLM $ warpfd
WfdLM <- regListLM $ Wfd

plot(accelfdLM)

wbasisCR <- create.bspline.basis(c(0,24), 15, 5)
Wfd0CR <- fd(matrix(0, 15, 11), wbasisCR)
WfdParCR <- fdPar(Wfd0CR, 1, 1)

regList <- register.fd(mean(accelfdLM), accelfdLM, WfdParCR)
accelfdCR <- regList $ regfd
warpfdCR <- regList $ warpfd
WfdCR <- regList $ Wfd

par(mfrow = c(2,2))
plot(accelfdUN)
plot(accelfdLM)
plot(accelfdCR)

plot(warpfdLM)
plot(warpfdCR)


# now try some regressions
# set up response
daymeans <- colSums(predict(tempfd)) + rnorm(11, 0, .1)

# set up predictors
templist <- vector("list",2)
templist[[1]] <- rep(1,11)
templist[[2]] <- tempfd

# set up parameter bases
betalist = list(
  create.constant.basis(c(0,24)),
  fdPar(create.bspline.basis(c(0,24),24, 6), harmaccelLfd, exp(10)))

# do regression
annPrecTemp = fRegress(daymeans, templist, betalist)

# confidence intervals on smooth
resid = daymeans - daymeanshat
SigmaE.= sum(resid^2)/(11-annPrecTemp$df)
SigmaE = SigmaE.*diag(rep(1,11))
stderrList = fRegress.stderr(annPrecTemp, tempSmooth$y2cMap, SigmaE)

betafd = annPrecTemp $ betaestlist[[2]] $ fd
betastderrfd = stderrList$betastderrlist[[2]]
plot(betafd, xlab="Day", ylab="Temperature Reg. Coeff.", ylim=c(3.5,4.2), lwd=2)
lines(betafd+2*betastderrfd, lty=2, lwd=1)
lines(betafd-2*betastderrfd, lty=2, lwd=1)

# assess fit and check that FDA is appropriate
# Ftests... and fit with a constant basis for 2nd element

}