
#' @import fda
NULL

#' @import reshape
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

preddat <- melt(preddat, id.vars = c("wk", "year"))


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
