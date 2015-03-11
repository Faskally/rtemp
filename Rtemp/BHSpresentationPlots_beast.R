

plot.pcs <- function(infodat, main = NULL) {
  par(mfrow = c(4, 4), mar = c(0,0,0,0), oma = c(2, 2, 5, 4))
  for (i in 1:4) {
    for (j in 1:4) {
      plot(infodat[[paste0("PC",i, ".x")]], infodat[[paste0("PC",j, ".y")]], ann = FALSE, axes = FALSE)
      box()
    } 
  }
  mtext(paste("PC", 1:4), side = 3, line = .5, at = 1:4 / 4 - .125, outer = TRUE)
  mtext(paste("PC", 1:4), side = 4, line = .5, at = 4:1 / 4 - .125, outer = TRUE, las = 1)
  mtext("Control", side = 1, line = .5, at = .5, outer = TRUE)
  mtext("Treatment", side = 2, line = .5, at = .5, outer = TRUE)
  if (!is.null(main)) mtext(main, side = 3, line = 3, at = 0.5, outer = TRUE, font = 2, cex = 1.5)
}

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
    #s(yday, by = PC1.x, bs = "cc", k = 4) + 
    #s(yday, by = PC2.x, bs = "cc", k = 4) + 
    #s(yday, by = PC3.x, bs = "cc", k = 4) + 
    #s(yday, by = PC4.x, bs = "cc", k = 4) +
    s(yday, bs = "cc", k = 6, by = fell) +
    s(yday, by = PC1.x * fell, bs = "cc", k = 4) + 
    s(yday, by = PC2.x * fell, bs = "cc", k = 4) + 
    s(yday, by = PC3.x * fell, bs = "cc", k = 4) + 
    s(yday, by = PC4.x * fell, bs = "cc", k = 4)# +
  #s(day, k = 10)
  
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


getFelledPC <- function(model, sim = FALSE, n = 1, which = grepl("fell", names(coef(model)))) {
  if (!sim) {
    b <- coef(model)
    if (is.null(which)) which <- rep(TRUE, length(b))
    if (!is.logical(which)) which <- 1:length(b) %in% which
    b[!which] <- 0
    c(predict(model, type = "lpmatrix") %*% b)
  } else {
    b <- coef(model)
    V <- model $ Vp
    bsim <- mvrnorm(n, b, V)
    if (is.null(which)) which <- rep(FALSE, length(b))    
    if (!is.logical(which)) which <- 1:length(b) %in% which
    bsim[,!which] <- 0
    predict(model, type = "lpmatrix") %*% t(bsim)
  }
} 




# load required libraries
library(lattice)
library(fda)
#library(reshape)
library(mgcv)
library(MASS)
library(nlme)


# get the data for each site:


# get full data set

choices <- data.frame(burn = c(10, 11, 7, 10, 11), year = c(1989, 1997, 1994, 2004, 2004))


for (ichoices in 1:5) {
  
  fname <- paste0("Z:/Loch_Ard_Felling/BHS/data/B", choices $ burn[ichoices], "_", choices $ year[ichoices], ".rda")
  load(fname)
  
  #load("Z:/Loch_Ard_Felling/BHS/data/B10_2004.rda")
  #load("Z:/Loch_Ard_Felling/BHS/data/B11_2004.rda")
  #load("Z:/Loch_Ard_Felling/BHS/data/B7_1994.rda")
  #load("Z:/Loch_Ard_Felling/BHS/data/B10_1989.rda")
  #load("Z:/Loch_Ard_Felling/BHS/data/B11_1997.rda")
  
  # sort out day
  info $ day <- info $ day - min(info $ day) + info $ yday[1]
  
  
  #plot.pcs(info, main = "All data") 
  
  #xyplot(PC1.y ~ yday | factor(year), data = info, as.table = TRUE)
  
  decays <- seq(0, 1, length = 11)
  gcv1 <- sapply(decays, fitDecay, gcv = TRUE, data = within(info, {y = PC1.y}))
  gcv2 <- sapply(decays, fitDecay, gcv = TRUE, data = within(info, {y = PC2.y}))
  gcv3 <- sapply(decays, fitDecay, gcv = TRUE, data = within(info, {y = PC3.y}))
  gcv4 <- sapply(decays, fitDecay, gcv = TRUE, data = within(info, {y = PC4.y}))
  gcvs <- gcv1 + gcv2 + gcv3 + gcv4
  
  decays2 <- seq(decays[max(1, which.min(gcvs) - 1)], decays[max(1, which.min(gcvs) + 1)], length = 11)
  gcv1 <- sapply(decays2, fitDecay, gcv = TRUE, data = within(info, {y = PC1.y}))
  gcv2 <- sapply(decays2, fitDecay, gcv = TRUE, data = within(info, {y = PC2.y}))
  gcv3 <- sapply(decays2, fitDecay, gcv = TRUE, data = within(info, {y = PC3.y}))
  gcv4 <- sapply(decays2, fitDecay, gcv = TRUE, data = within(info, {y = PC4.y}))
  gcvs2 <- gcv1 + gcv2 + gcv3 + gcv4
  
  gcvs <- c(gcvs, gcvs2)
  decays <- c(decays, decays2)
  
  plot(decays, gcvs)
  
  decay <- decays[which.min(gcvs)]
  full1 <- fitDecay(decay, data = within(info, {y = PC1.y}))
  full2 <- fitDecay(decay, data = within(info, {y = PC2.y}))
  full3 <- fitDecay(decay, data = within(info, {y = PC3.y}))
  full4 <- fitDecay(decay, data = within(info, {y = PC4.y}))
  
  info $ fell <- feffect(info $ day, decay = decay, type = 1)
  xyplot(fell ~ day, data = info)
  
  X <- eval.fd(seq(0, 24, length = 96), pc10 $ harmonics)
  Beta <- rbind(0, getFelledPC(full1), getFelledPC(full2), getFelledPC(full3), getFelledPC(full4))
  
  Beta.all <- rbind(1, getFelledPC(full1, which = NULL), getFelledPC(full2, which = NULL), getFelledPC(full3, which = NULL), getFelledPC(full4, which = NULL))
  
  Xmu <- eval.fd(seq(0, 24, length = 96), pc10 $ mean)
  X <- cbind(Xmu, X)
  
  
  fits <- X %*% Beta
  fits.all <- X %*% Beta.all
  
  # we can compute daily summaries of felling effects
  info $ dailyMax <- apply(fits, 2, max)
  info $ dailyMax.all <- apply(fits.all, 2, max)
  
  
  # what about confidence intervals...
  n <- 1000
  sim1 <- getFelledPC(full1, sim = TRUE, n = n)
  sim2 <- getFelledPC(full2, sim = TRUE, n = n)
  sim3 <- getFelledPC(full3, sim = TRUE, n = n)
  sim4 <- getFelledPC(full4, sim = TRUE, n = n)
  
  simfits <- sapply(1:n, function(i) X %*% rbind(0, sim1[,i], sim2[,i], sim3[,i], sim4[,i]))
  dim(simfits) <- c(nrow(X), ncol(Beta), n)
  
  # we can compute daily summaries of felling effects
  simdailyMax <- apply(simfits, 2:3, max)
  info $ dailyMax.cil <- apply(simdailyMax, 1, quantile, .025)
  info $ dailyMax.ciu <- apply(simdailyMax, 1, quantile, .975)
  
  
  p <- xyplot(dailyMax + dailyMax.cil  + dailyMax.ciu ~ yday | factor(year), data = info, 
         pch = 19, cex = .7, as.table = TRUE, grid = TRUE, type = c("p"),
         main = paste0(choices $ burn[ichoices], "_", choices $ year[ichoices]))
  
  plot(p)

  info $ event <- do.call(paste, c(choices[ichoices,], list(sep = "_"))) 
  
  
  vname <- do.call(paste, c(choices[ichoices,], list(sep = "_"))) 
  assign(paste0("info", vname), info, envir = .GlobalEnv)
  #assign(paste0("simfits", vname), simfits, envir = .GlobalEnv)
  
  fname <- paste0("Z:/Loch_Ard_Felling/BHS/data/B", choices $ burn[ichoices], "_", choices $ year[ichoices], "fits.rda")
  save(list = paste0(c("info"), vname),  file = fname)
  
}


info.all <- rbind(info10_1989, info10_2004, info11_1997, info11_2004, info7_1994)


xyplot(dailyMax ~ yday | year * event, data = info.all)
fname <- paste0("Z:/Loch_Ard_Felling/BHS/data/B", choices $ burn[ichoices], "_", choices $ year[ichoices], "fits.rda")
save(info.all,  file = "Z:/Loch_Ard_Felling/BHS/data/infoall.rda")

