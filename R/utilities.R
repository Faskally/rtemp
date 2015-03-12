
#' @import fda
NULL


#' Run a full calibration
#' 
#' This function does everything from read in files, select calibtation period and
#' thens runs a calibration and saves the file.
#'
#' @param object The number of simulations to use to evaluate confidence intervals.
#' @param nbsplines A number.
#' @param lambda A number.
#' @param breaks A number.
#' @return NULL
#' @export
#' @examples
#' 1 + 1
getfdyr <- function(object, nbsplines = 24, lambda = 0, breaks = NULL) {

  # create some basis functions
  daybasis365 <- create.bspline.basis(c(0, 365), nbsplines, norder = 6, breaks = breaks)

  # create roughness penalty
  harmaccelLfd <- vec2Lfd(c(0, (2*pi/365)^2, 0), c(0, 365))

  # so now estimate the coefficients for each day, year and site combo
  ck <-
    by(object, interaction(object $ year), 
      function(x) {
        if (sum(0:365 %in% x $ yday) < 200) return(NULL)
        sm <- try(smooth.basis(x $ yday, x $ control, fdPar(daybasis365, harmaccelLfd, lambda)))
        if (inherits(sm, "try-error")) return(NULL)
        list(coef = c(coef(sm)), gcv = sm $ gcv)
      })
  # next line removes NULL elements
  ck <- ck[!sapply(ck, is.null)]
  if (length(ck) == 0) stop("no days had sufficient coverage to fit ", nbsplines, " bsplines")
  gcv <- sapply(ck, "[[", "gcv")
  names(gcv) <- names(ck)
  ck <- sapply(ck, "[[", "coef")
  
  # now create a functional data object and fill it with the relavent stuff
  list(fd = fd(ck, daybasis365, fdnames= list(time = 0:365, reps = colnames(ck), values = "value")),
       gcv = gcv)
}





#' Run a full calibration
#' 
#' This function does everything from read in files, select calibtation period and
#' thens runs a calibration and saves the file.
#'
#' @param object The number of simulations to use to evaluate confidence intervals.
#' @param nbsplines A number.
#' @param lambda A number.
#' @param breaks A number.
#' @return NULL
#' @export
#' @examples
#' 1 + 1
getfd <- function(object, nbsplines = 24, lambda = 0, breaks = NULL) {

  # add in a column for decimal hour
  object $ dhour <- with(object, hour + min/60)

  # create some basis functions
  hourbasis24 <- create.bspline.basis(c(0, 24), nbsplines, norder = 6, breaks = breaks)

  # create roughness penalty
  harmaccelLfd <- vec2Lfd(c(0, (2*pi/24)^2, 0), c(0, 24))
  Rmat <- eval.penalty(hourbasis24, harmaccelLfd)

  # so now estimate the coefficients for each day, year and site combo
  ck <-
    by(object, interaction(object $ yday, object $ site, object $ year), 
      function(x) {
        if (!all(0:23 %in% x $ hour)) return(NULL)
        #Phi <- eval.basis(x $ dhour, hourbasis24)
        #XX <- crossprod(Phi) + lambda * Rmat
        #if (any(eigen(XX, only.values = TRUE) $ values <= 0)) return(NULL)
        #drop(solve(XX, crossprod(Phi, x $ Original.Value)))
        sm <- try(smooth.basis(x $ dhour, x $ Original.Value, fdPar(hourbasis24, harmaccelLfd, lambda)))
        if (inherits(sm, "try-error")) return(NULL)
        list(coef = coef(sm), gcv = sm $ gcv)
      })
  # next line removes NULL elements
  ck <- ck[!sapply(ck, is.null)]
  if (length(ck) == 0) stop("no days had sufficient coverage to fit ", nbsplines, " bsplines")
  gcv <- sapply(ck, "[[", "gcv")
  names(gcv) <- names(ck)
  ck <- sapply(ck, "[[", "coef")
  
  # now create a functional data object and fill it with the relavent stuff
  list(fd = fd(ck, hourbasis24, fdnames= list(time = 0:24, reps = colnames(ck), values = "value")),
       gcv = gcv)
}


#' Run a full calibration
#' 
#' This function does everything from read in files, select calibtation period and
#' thens runs a calibration and saves the file.
#'
#' @param object The number of simulations to use to evaluate confidence intervals.
#' @param nbsplines A number.
#' @param lambda A number.
#' @param breaks A number.
#' @return NULL
#' @export
#' @examples
#' 1 + 1
getgcv <- function(object, nbsplines = 24, lambda = 0, breaks = NULL) {

  # add in a column for decimal hour
  object $ dhour <- with(object, hour + min/60)

  # create some basis functions
  hourbasis24 <- create.bspline.basis(c(0, 24), nbsplines, norder = 6, breaks = breaks)

  # create roughness penalty
  harmaccelLfd <- vec2Lfd(c(0, (2*pi/24)^2, 0), c(0, 24))

  # so now estimate the coefficients for each day, year and site combo
  ck <-
    by(object, interaction(object $ yday, object $ site, object $ year), 
      function(x) {
        if (!all(0:23 %in% x $ hour)) return(NULL)
        sm <- try(smooth.basis(x $ dhour, x $ Original.Value, fdPar(hourbasis24, harmaccelLfd, lambda)))
        if (inherits(sm, "try-error")) return(NULL)
        sm $ gcv
      })
  # next line removes NULL elements
  ck <- ck[!sapply(ck, is.null)]
  gcv <- unlist(ck)
  names(gcv) <- names(ck)
  gcv
}


#' Run a full calibration
#' 
#' This function does everything from read in files, select calibtation period and
#' thens runs a calibration and saves the file.
#'
#' @param object The number of simulations to use to evaluate confidence intervals.
#' @param pc A number.
#' @return NULL
#' @export
#' @examples
#' 1 + 1
reduceData <- function(object, pc) {

  # now convert newdata to basis coefficients using the new bases.
  object $ dhour <- with(object, hour + min/60)

  dk <-
    by(object, interaction(object $ yday, object $ site, object $ year)[drop = TRUE], 
      function(x) {
        if (is.null(x)) return(rep(NA, length(pc $ meanfd $ basis $ params)))
        Phi <- eval.fd(x $ dhour, pc $ harmonics)
        # substract off mean function
        y <- x $ Original.Value - eval.fd(x $ dhour, pc $ meanfd)
        coef(lm( y ~ Phi - 1))
      })
  dk <- simplify2array(unclass(dk))

  # so now we have our new data!!
  out <- unique(object[c("year", "yday", "mon", "mday")])
  out <- out[rep(1:nrow(out), each = length(pc $ harmonics $ fdnames[[2]])), ]
  out $ pc <- pc $ harmonics $ fdnames[[2]]
  out $ coef <- c(dk)
  rownames(out) <- NULL

  out
}


#' Run a full calibration
#' 
#' This function does everything from read in files, select calibtation period and
#' thens runs a calibration and saves the file.
#'
#' @param dat The number of simulations to use to evaluate confidence intervals.
#' @param PC A number.
#' @param ncuts A number.
#' @return NULL
#' @export
#' @examples
#' 1 + 1
plotPCbyblock <- function(dat, PC = "PC1", ncuts = 4) {
#PC <- "PC1"
dat1 <- dat[dat $ pc == PC, ]
if (ncuts == 1) {
  dat1 $ wk <- "1"
} else {
  dat1 $ wk <- cut(dat1 $ yday, ncuts)
}
m1 <- lm(felled ~ factor(year)*wk*control, data = dat1)

# predict intercept
preddat <- expand.grid(wk = levels(dat1 $ wk[drop = TRUE]),
                       year  = unique(dat1 $ year))

preddat $ intercept <- predict(m1, newdata = cbind(preddat, control = 0))
preddat $ slope <- predict(m1, newdata = cbind(preddat, control = 1)) - preddat $ intercept

preddat <- reshape::melt(preddat, id.vars = c("wk", "year"))


ylim <- switch(PC, PC1 = list(c(-15, 25), c(0, 1.5)),
                   PC2 = list(c(-2, 2), c(0, 1.5)),
                   PC3 = list(c(-1.5, 1), c(0, 1.5)),
                   PC4 = list(c(-1, 1), c(0, 1.5)))

lattice::xyplot(value ~ wk | factor(year) + variable, data = preddat, 
                     type = c("p", "l", "g"), ylim = ylim[rep(1:2, each = length(unique(preddat $ year)))], 
                     as.table = TRUE, col = 1, pch = 16, cex = .65,
                     scales = list(relation = "free"), 
                     ylab = "monthly variation", xlab = "", main = PC,
                     axis = function(side, ...) {
                       if (lattice::current.column() == 1 & side == "left") {
                        lattice::panel.axis(side, outside = TRUE, tck = .3, text.cex = .7)
                       }
                     },
                     par.settings = list(layout.heights = list(axis.panel = 0), 
                                         layout.widths  = list(axis.panel = 0, 
                                                               ylab.axis.padding = 4),
                                         axis.components = list(left = list(pad1 = .2, pad2 = 2-0.2))))
}





#' Run a full calibration
#' 
#' This function does everything from read in files, select calibtation period and
#' thens runs a calibration and saves the file.
#'
#' @param infodat The number of simulations to use to evaluate confidence intervals.
#' @param main A number.
#' @return NULL
#' @export
#' @examples
#' 1 + 1
plot_pcs <- function(infodat, main = NULL) {
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



 
#' felling effect
#' 
#' This function does everything from read in files, select calibtation period and
#' thens runs a calibration and saves the file.
#'
#' @param day The number of simulations to use to evaluate confidence intervals.
#' @param decay A number.
#' @param b1 A number.
#' @param a1 A number.
#' @param type A number.
#' @return NULL
#' @export
#' @examples
#' 1 + 1
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



#' felling effect
#' 
#' This function does everything from read in files, select calibtation period and
#' thens runs a calibration and saves the file.
#'
#' @param decay The number of simulations to use to evaluate confidence intervals.
#' @param gcv A number.
#' @param data A number.
#' @return NULL
#' @export
#' @examples
#' 1 + 1
fitDecay <- function(decay, gcv = FALSE, data) {
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
  
  # generate AR1 structure
  V <- nlme::corMatrix(nlme::Initialize(nlme::corAR1(.9), data.frame(x = tmp $ day)))
  Cv <- chol(V)  # t(Cv)%*%Cv=V
  w <- solve(t(Cv)) # V^{-1} = w'w

  full <- mgcv::gam(form, data = tmp)
  
  # Use `gam' to set up model for fitting...
  G <- mgcv::gam(form, data = tmp, fit = FALSE)
  # fit using magic, with weight *matrix*
  mgfit <- mgcv::magic(G $ y, G $ X, G $ sp, G $ S, G $ off, rank = G $ rank, C = G $ C, w = w)
  # Modify previous gam object using new fit, for plotting...    
  mg.stuff <- mgcv::magic.post.proc(G $ X, mgfit, w)
  full $ edf[] <- mg.stuff $ edf
  full $ Vp[] <- mg.stuff $ Vb
  full $ coefficients[] <- mgfit $ b 
  
  if (gcv) full $ gcv.ubre else full
}



#' felling effect
#' 
#' This function does everything from read in files, select calibtation period and
#' thens runs a calibration and saves the file.
#'
#' @param model The number of simulations to use to evaluate confidence intervals.
#' @param sim A number.
#' @param n A number.
#' @param which A number.
#' @return NULL
#' @export
#' @examples
#' 1 + 1
getFelledPC <- function(model, sim = FALSE, n = 1, 
                        which = grepl("fell", names(coef(model)))) {
  if (!sim) {
    b <- coef(model)
    if (is.null(which)) which <- rep(TRUE, length(b))
    if (!is.logical(which)) which <- 1:length(b) %in% which
    b[!which] <- 0
    c(predict(model, type = "lpmatrix") %*% b)
  } else {
    b <- coef(model)
    V <- model $ Vp
    bsim <- MASS::mvrnorm(n, b, V)
    if (is.null(which)) which <- rep(FALSE, length(b))    
    if (!is.logical(which)) which <- 1:length(b) %in% which
    bsim[,!which] <- 0
    predict(model, type = "lpmatrix") %*% t(bsim)
  }
} 


