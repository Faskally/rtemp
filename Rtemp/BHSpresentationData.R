
# load required libraries
library(lattice)
library(fda)
#library(reshape)
library(mgcv)
library(MASS)
library(nlme)


# get the data for each site:


# get full data set
load("B:/Loch_Ard_Felling/Package/rivertemp/data/fullwater.rda")
# a quick fix
stream $ year[stream $ Sensor.Model == "Squirrel"] <- stream $ year[stream $ Sensor.Model == "Squirrel"] + 1900 



# set up a fourier basis of dimension 48 and order 6 
hourbasis <- create.fourier.basis(c(0, 24), 47, 24)
# and a harmonic accelerator penalty... penalises deviations from sinusiodal curves
penalty <- vec2Lfd(c(0, (2*pi/24)^2, 0), c(0, 24)) 

# would be good to use the same bases for all experiments...

for (choice in 1:5) {

cat("choice:", choice, "\n"); flush.console()

if (choice == 1) {
  years <- 2003:2009
  site2 <- 2
  site10 <- 10
} else if (choice == 2) {
  years <- 2003:2009
  site2 <- 2
  site10 <- 11
} else if (choice == 3) {
  years <- 1988:1994
  site2 <- 11
  site10 <- 10  
} else if (choice == 4) {
  years <- 1996:2002
  site2 <- 10
  site10 <- 11  
} else if (choice == 5) {
  years <- 1993:1995
  site2 <- 11
  site10 <- 7  
}


# select a few days data and add in decimal hour
dat <- with(subset(stream, site == site2 & year %in% years), 
            data.frame(temp = Original.Value, 
                       dhour = hour + min/60, 
                       month = factor(month.abb[mon + 1], levels = month.abb), 
                       yday = yday, 
                       year = year,  
                       site = site)) 
dat <- subset(dat, yday < 365) # remove last day of leap year

getcoef <- function(lambda, .yday, .year = 2009, .site = site2) {
  sdat <- subset(dat, yday == .yday & year == .year & site == .site)
  #if (nrow(sdat) < 23) return (NULL)
  if (sum(unique(round(sdat $ dhour)) %in% 0:24) < 24) return(NULL)  
  c(coef(smooth.basis(sdat $ dhour, sdat $ temp, fdPar(hourbasis, penalty, lambda))))
}

yday <- rep(0:365, length(years))
year <- rep(years, each = 366)

ck <- mapply(getcoef, .yday = yday, .year = year, MoreArgs = list(lambda = exp(5)))
whichNULL <- sapply(ck, is.null)
ck <- simplify2array(ck[!whichNULL])

info2 <-data.frame(yday = yday[!whichNULL], year = year[!whichNULL])

tempfd2 <- fd(ck, hourbasis)
pc2 <- pca.fd(tempfd2, nharm = 4) # is this sensitive to missing values?...

par(mfrow = c(2,2))
plot(pc2)

info2[pc2 $ harmonics $ fdnames[[2]]] <- as.data.frame(pc2 $ scores)



# select a few days data and add in decimal hour
dat <- with(subset(stream, site == site10 & year %in% years), 
            data.frame(temp = Original.Value, 
                       dhour = hour + min/60, 
                       month = factor(month.abb[mon + 1], levels = month.abb), 
                       yday = yday, 
                       year = year, 
                       site = site))  
dat <- subset(dat, yday < 365) # remove last day of leap year

getcoef <- function(lambda, .yday, .year = 2009, .site = site10) {
  sdat <- subset(dat, yday == .yday & year == .year & site == .site)
  #if (nrow(sdat) < 23) return (NULL)
  if (sum(unique(round(sdat $ dhour)) %in% 0:24) < 24) return(NULL)  
  c(coef(smooth.basis(sdat $ dhour, sdat $ temp, fdPar(hourbasis, penalty, lambda))))
  #c(coef(smooth.basis(sdat $ dhour, sdat $ temp, pc2 $ harmonics)))
}

yday <- rep(0:365, length(years))
year <- rep(years, each = 366)

ck <- mapply(getcoef, .yday = yday, .year = year, MoreArgs = list(lambda = exp(5)))
whichNULL <- sapply(ck, is.null)
ck <- simplify2array(ck[!whichNULL])

info10 <-data.frame(yday = yday[!whichNULL], year = year[!whichNULL])

tempfd10 <- fd(ck, hourbasis)
pc10 <- pca.fd(tempfd10, nharm = 4) # is this sensitive to missing values?...

par(mfrow = c(2,2))
plot(pc10)

info10[pc10 $ harmonics $ fdnames[[2]]] <- as.data.frame(pc10 $ scores)

info <- merge(info2, info10, by = c("yday", "year"), suffixes = c(".2", ".10"))
info $ day <- with(info, yday + (year - 2003) * 365)
info $ month <- month.abb[strptime(paste(info $ yday + 1, info $ year), format = "%j %Y") $ mon + 1]
info $ season <- with(info, ifelse(month %in% c("Nov", "Dec", "Jan", "Feb"), "Winter", 
                                   ifelse(month %in% c("Mar", "Apr"), "Spring",
                                          ifelse(month %in% c("May", "Jun", "Jul", "Aug"), "Summer", "Autumn"))))
info $ season <- factor(info $ season, levels = c("Winter", "Spring", "Summer", "Autumn"))
info <- info[order(info $ day),]  

names(info) <- gsub("[.]2", ".x", names(info))
names(info) <- gsub("[.]10", ".y", names(info))



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



fname <- paste0("Z:/Loch_Ard_Felling/BHS/data/B", site10, "_", years[2], ".rda")
save(info, pc10, w, file = fname)

}

