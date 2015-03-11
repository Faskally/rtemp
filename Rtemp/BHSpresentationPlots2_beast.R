
# load required libraries
library(lattice)
library(latticeExtra)
library(fda)

load("B:/Loch_Ard_Felling/BHS/data/infoall.rda")



info.all $ cevent <- factor(info.all $ event, levels = c("10_1989", "7_1994", "11_1997", "11_2004", "10_2004"))

p <-
xyplot(dailyMax ~ day | cevent, data = subset(info.all, day > 365), layout = c(1,5), as.table = TRUE, pch = 19, groups = dailyMax.cil > 0, cex = .5, col = c(1, 2),
       xlim = c(365, 365*7),
       xlab = "", ylab = "Effect on Daily Maximum", strip = FALSE, strip.left = TRUE,
       scale = list(x = list(at = 365 * 0:(4*8)/4, labels = rep(c("Jan", "Apr", "Jul", "Oct"), 8)))) +
xyplot(dailyMax ~ day | cevent, data = subset(info.all, day > 365), 
       panel = function(x, y) {
         panel.abline(h = c(-2.5, 2.5, 7.5), col = grey(0.9))
         panel.abline(h = c(-5, 5, 10), v = 365 * 1:(4*8)/4, col = grey(0.8))
         panel.abline(h = 0, v = 365 * 1:8, col = grey(0.5), lwd = 2)
       })

ppi <- 150
width <- 60
height <- 27.5
png(file = "B:/Loch_Ard_Felling/BHS/FellingEffect.png", width / 2.54 * ppi, height / 2.54 * ppi, res = ppi)
print(p)
dev.off()

#

load("B:/Loch_Ard_Felling/BHS/data/B10_2004.rda")
library(fda)

ppi <- 150
width <- 40
height <- 12
png(file = "B:/Loch_Ard_Felling/BHS/FPCA.png", width / 2.54 * ppi, height / 2.54 * ppi, res = ppi)
par(mfrow = c(1,4))
plot(pc10)
dev.off()

#
# get full data set
load("B:/Loch_Ard_Felling/Package/rivertemp/data/fullwater.rda")
# a quick fix
stream $ year[stream $ Sensor.Model == "Squirrel"] <- stream $ year[stream $ Sensor.Model == "Squirrel"] + 1900 

stream $ felled <- FALSE
stream $ felled[with(stream, (site == 10 & year %in% c(1989, 2004)) | 
                             (site == 11 & year %in% c(1997, 2004)) |
                             (site == 7 & year == 1994))] <- TRUE
p1 <- xyplot(Original.Value ~ date | paste0("Burn ", site), data = stream, groups = felled, layout = c(1,4), 
             pch = ".", col = c(1, 2), strip = FALSE, strip.left = TRUE, 
             ylab = "Raw temperature", xlab = "", grid = TRUE) 

ppi <- 150
width <- 40
height <- 12
png(file = "B:/Loch_Ard_Felling/BHS/data.png", width / 2.54 * ppi, height / 2.54 * ppi, res = ppi)
plot(p1)
dev.off()



load("B:/Loch_Ard_Felling/BHS/data/infoall.rda")

head(info.all)



with(info.all, tapply(dailyMax, list(event, year), max))  

with(info.all, tapply(fell, list(event, year), min))  
