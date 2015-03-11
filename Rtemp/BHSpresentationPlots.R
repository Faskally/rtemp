



# load required libraries
library(lattice)
library(fda)
#library(reshape)
library(mgcv)
library(MASS)
library(nlme)



# get full data set
load("Z:/Loch_Ard_Felling/Package/rivertemp/data/fullwater.rda")
# a quick fix
stream $ year[stream $ Sensor.Model == "Squirrel"] <- stream $ year[stream $ Sensor.Model == "Squirrel"] + 1900 



# set up a bspline basis of dimension 48 and order 6 
hourbasis <- create.bspline.basis(c(0, 24), 48, 6)
# and a harmonic accelerator penalty... penalises deviations from sinusiodal curves
penalty <- vec2Lfd(c(0, (2*pi/24)^2, 0), c(0, 24)) 



# select a few days data and add in decimal hour
dat <- with(subset(stream, site == 2 & year %in% 2003:2009), 
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

info2 <-data.frame(yday = yday[!whichNULL], year = year[!whichNULL])

tempfd2 <- fd(ck, hourbasis)
pc2 <- pca.fd(tempfd2, nharm = 4) # is this sensitive to missing values?...

info2[pc2 $ harmonics $ fdnames[[2]]] <- as.data.frame(pc2 $ scores)



# select a few days data and add in decimal hour
dat <- with(subset(stream, site == 10 & year %in% 2003:2009), 
    data.frame(temp = Original.Value, 
               dhour = hour + min/60, 
               month = factor(month.abb[mon + 1], levels = month.abb), 
               yday = yday, 
               year = year, 
               site = site))  
dat <- subset(dat, yday < 365) # remove last day of leap year

getcoef <- function(lambda, .yday, .year = 2009, .site = 10) {
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

info10 <-data.frame(yday = yday[!whichNULL], year = year[!whichNULL])

tempfd10 <- fd(ck, hourbasis)
pc10 <- pca.fd(tempfd10, nharm = 4) # is this sensitive to missing values?...

info10[pc10 $ harmonics $ fdnames[[2]]] <- as.data.frame(pc10 $ scores)



info <- merge(info2, info10, by = c("yday", "year"), suffixes = c(".2", ".10"))
info $ day <- with(info, yday + (year - 2003) * 365)
info $ month <- month.abb[strptime(paste(info $ yday + 1, info $ year), format = "%j %Y") $ mon + 1]
info $ season <- with(info, ifelse(month %in% c("Nov", "Dec", "Jan", "Feb"), "Winter", 
                              ifelse(month %in% c("Mar", "Apr"), "Spring",
                                ifelse(month %in% c("May", "Jun", "Jul", "Aug"), "Summer", "Autumn"))))
info $ season <- factor(info $ season, levels = c("Winter", "Spring", "Summer", "Autumn"))
info <- info[order(info $ day),]  



plot.pcs <- function(infodat, main = NULL) {
  par(mfrow = c(4, 4), mar = c(0,0,0,0), oma = c(2, 2, 5, 4))
  for (i in 1:4) {
    for (j in 1:4) {
      plot(infodat[[paste0("PC",i, ".2")]], infodat[[paste0("PC",j, ".10")]], ann = FALSE, axes = FALSE)
      box()
    } 
  }
  mtext(paste("PC", 1:4), side = 3, line = .5, at = 1:4 / 4 - .125, outer = TRUE)
  mtext(paste("PC", 1:4), side = 4, line = .5, at = 4:1 / 4 - .125, outer = TRUE, las = 1)
  mtext("Burn 2", side = 1, line = .5, at = .5, outer = TRUE)
  mtext("Burn 10", side = 2, line = .5, at = .5, outer = TRUE)
  if (!is.null(main)) mtext(main, side = 3, line = 3, at = 0.5, outer = TRUE, font = 2, cex = 1.5)
}
plot.pcs(info, main = "All data") 



Qfunc <- function(n, phi, subset = 1:n) {
  Q <- Diagonal(n)
  diag(Q[2:(n-1),2:(n-1)]) <- 1 + phi^2
  diag(Q[,-1]) <- diag(Q[-1,]) <- -phi
  # work out Q for a subset of 1:n
  # hack:
  V <- solve(Q)
  Vs <- V[subset,subset]
  Qs <- solve(Vs)
  out <- Diagonal(nrow(Qs))
  #Qs
  diag(out) <- diag(Qs)
  diag(out[,-1]) <- diag(out[-1,]) <- diag(Qs[-1,])  
  out
}

## produce scaled covariance matrix for AR1 errors...
#Q <- Qfunc(max(info $ day + 1), .6, subset = info $ day + 1)
#w <- chol(Q)

V <- corMatrix(Initialize(corAR1(.9), data.frame(x = info $ day)))
Cv <- chol(V)  # t(Cv)%*%Cv=V
w <- solve(t(Cv)) # V^{-1} = w'w
# check w and w2 are identical!




# felling effect
feffect <- function(day, decay = 1, b1 = 100, a1 = 0.25, type = 2) {
  if (type == 1) {
    ifelse(day < 365, 0, 
      exp(ifelse(day < 365 * 2, 0,
        -decay/365 * (day - 2 * 365))))
  } else if (type == 2) {
    a1 <- a1 * 365
    a2 <- decay * 365
    ifelse(day < (1.25*365),  
      exp(-(abs(day - 1.25*365)/a1)^b1),
       ifelse( day >= (1.25*365) & day < (1.75*365),
       1, exp(-(abs(day - 1.75*365)/a2)^2)))
  }
} 

 
fitDecay <- function(decay, gcv = FALSE, data = info) {
  tmp <- data
  tmp $ fell <- feffect(tmp $ day, decay = decay, type = 1)

  form <- y ~ (PC1.x + PC2.x + PC3.x + PC4.x) * fell +
            s(yday, bs = "cc", k = 6) + 
            s(yday, by = PC1.x, bs = "cc", k = 4) + 
            s(yday, by = PC2.x, bs = "cc", k = 4) + 
            s(yday, by = PC3.x, bs = "cc", k = 4) + 
            s(yday, by = PC4.x, bs = "cc", k = 4) +
            s(yday, bs = "cc", k = 6, by = fell) +
            s(yday, by = PC1.x * fell, bs = "cc", k = 4) + 
            s(yday, by = PC2.x * fell, bs = "cc", k = 4) + 
            s(yday, by = PC3.x * fell, bs = "cc", k = 4) + 
            s(yday, by = PC4.x * fell, bs = "cc", k = 4) +
            s(day, k = 10)

  full <- gam(form, data = tmp)

  # Use `gam' to set up model for fitting...
  G <- gam(form, data = tmp, fit = FALSE)
  # fit using magic, with weight *matrix*
  mgfit <- magic(G $ y, G $ X, G $ sp, G $ S, G $ off, rank = G $ rank, C = G $ C, w = w)
  # Modify previous gam object using new fit, for plotting...    
  mg.stuff <- magic.post.proc(G $ X, mgfit, w)
  full $ edf[] <- mg.stuff $ edf
  full $ Vp[] <- mg.stuff $ Vb
  full $ coefficients[] <- mgfit $ b 

  if (gcv) full $ gcv.ubre else full
}

#full <- fitDecay(1)

dat <- info
names(dat) <- gsub("[.]2", ".x", names(dat))

decays <- c(seq(0, 0.5, length = 4), 1:5)
gcv1 <- sapply(decays, fitDecay, gcv = TRUE, data = within(dat, {y = PC1.10}))
gcv2 <- sapply(decays, fitDecay, gcv = TRUE, data = within(dat, {y = PC2.10}))
gcv3 <- sapply(decays, fitDecay, gcv = TRUE, data = within(dat, {y = PC3.10}))
gcv4 <- sapply(decays, fitDecay, gcv = TRUE, data = within(dat, {y = PC4.10}))
gcvs <- gcv1 + gcv2 + gcv3 + gcv4
plot(decays, gcvs)

full1 <- fitDecay(decays[which.min(gcvs)], data = within(dat, {y = PC1.10}))
full2 <- fitDecay(decays[which.min(gcvs)], data = within(dat, {y = PC2.10}))
full3 <- fitDecay(decays[which.min(gcvs)], data = within(dat, {y = PC3.10}))
full4 <- fitDecay(decays[which.min(gcvs)], data = within(dat, {y = PC4.10}))

plot(full1, scale = 0, page = 1)
summary(full1)

plot(full2, scale = 0, page = 1)
summary(full2)

plot(full3, scale = 0, page = 1)
summary(full3)

plot(full4, scale = 0, page = 1)
summary(full4)



getFelledPC <- function(model, sim = FALSE, n = 1) {
  if (!sim) {
    b <- coef(model)
    b[-grep("fell", names(b))] <- 0
    c(predict(model, type = "lpmatrix") %*% b)
  } else {
    b <- coef(model)
    V <- model $ Vp
    bsim <- mvrnorm(n, b, V)
    bsim[,!grepl("fell", names(b))] <- 0
    predict(model, type = "lpmatrix") %*% t(bsim)
  }
} 

X <- eval.fd(1:24, pc10 $ harmonics)

fits <- X %*% rbind(getFelledPC(full1), getFelledPC(full2), getFelledPC(full3), getFelledPC(full4))

plot(c(fits[, 365 * 4 + 200 + 1]), type = "l")

# we can compute daily summaries of felling effects
dat $ dailyMax <- apply(fits, 2, max)
xyplot(dailyMax ~ yday | year, data = dat, pch = ".")

dat $ dailyMean <- apply(fits, 2, mean)
xyplot(dailyMax ~ yday | year, data = dat, pch = ".")

# what about confidence intervals...
sim1 <- getFelledPC(full1, sim = TRUE, n = 100)
sim2 <- getFelledPC(full2, sim = TRUE, n = 100)
sim3 <- getFelledPC(full3, sim = TRUE, n = 100)
sim4 <- getFelledPC(full4, sim = TRUE, n = 100)

simfits <- sapply(1:100, function(i) X %*% rbind(sim1[,i], sim2[,i], sim3[,i], sim4[,i]))
dim(simfits) <- c(24, 2138, 1000)

# we can compute daily summaries of felling effects
simdailyMax <- apply(simfits, 2:3, max)
dat $ dailyMax.cil <- apply(simdailyMax, 1, quantile, .025)
dat $ dailyMax.ciu <- apply(simdailyMax, 1, quantile, .975)
xyplot(dailyMax.cil + dailyMax.ciu ~ yday | factor(year), data = dat, pch = ".", as.table = TRUE, grid = TRUE)

simdailyMean <- apply(simfits, 2:3, mean)
dat $ dailyMean.cil <- apply(simdailyMean, 1, quantile, .025)
dat $ dailyMean.ciu <- apply(simdailyMean, 1, quantile, .975)
xyplot(dailyMean.cil + dailyMean.ciu ~ yday | factor(year), data = dat, pch = ".", as.table = TRUE, grid = TRUE)




