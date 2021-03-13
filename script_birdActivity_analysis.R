#### aves bosques de la chec

library(tidyverse)


chec <- read.delim("bird_activity_data.txt")

str(chec)

for (i in 1:ncol(chec)) {
  if (is.character(chec[[i]])) {
    chec[[i]] <- as.factor(chec[[i]])
  } else {
    next
  }
}

library(dplyr)
chec %>%
filter(species == "mille") %>%
  arrange(forest)# 28 observaciones, 2 en reforestacion y 26 en BA

chec %>%
  filter(species == "rufi") %>%
  arrange(forest)# 39 observaciones, 4 en reforestación y 35 en BA
  
chec %>%
  filter(species == "frena") %>%
  arrange(forest)

chec %>%
  filter(species == "goud") %>%
  arrange(forest)


head(chec)
summary(chec)
range(chec$time)

timeRad <- chec$time * 2 * pi
table(chec$species,chec$forest)

library(boot)
library(circular)
library(overlap)
library(cowplot)
library(gridExtra)
library(png)


#
#patrón de actividad por especie----
#

#Arremon torquatus
torquatus <- timeRad[chec$species == "torq"]
densityPlot(torquatus, rug=TRUE, ylim=c(0,0.22), ylab="Density of activity", cex.axis=1.1, cex.lab=1.2, extend = NULL)
abline(v=6, lty="dashed", col="black", lwd=2); abline(v=18, lty="dashed", col="black", lwd=2)

#Arremon brunneinucha
brunneinucha <- timeRad[chec$species == "brun"]
densityPlot(brunneinucha, rug=TRUE, ylim=c(0,0.22), ylab="Density of activity", cex.axis=1.1, cex.lab=1.2, extend = NULL)
abline(v=6, lty="dashed", col="black", lwd=2); abline(v=18, lty="dashed", col="black", lwd=2)

#Grallaria ruficapilla
ruficapilla <- timeRad[chec$species == "rufi"]
densityPlot(ruficapilla, rug=TRUE, ylim=c(0,0.22), ylab="Density of activity", cex.axis=1.1, cex.lab=1.2, extend = NULL)
abline(v=6, lty="dashed", col="black", lwd=2); abline(v=18, lty="dashed", col="black", lwd=2) 

#Grallaria milleri
milleri <- timeRad[chec$species == "mille"]
densityPlot(milleri, rug=TRUE, ylim=c(0,0.22), ylab="Density of activity", cex.axis=1.1, cex.lab=1.2, extend = NULL)
abline(v=6, lty="dashed", col="black", lwd=2); abline(v=18, lty="dashed", col="black", lwd=2) 

#Chamaepetes goudotii
goudotii <- timeRad[chec$species == "goud"]
densityPlot(goudotii, rug=TRUE, ylim=c(0,0.22), ylab="Density of activity", cex.axis=1.1, cex.lab=1.2, extend = NULL)
abline(v=6, lty="dashed", col="black", lwd=2); abline(v=18, lty="dashed", col="black", lwd=2) 

#Zentrigon frenata
frenata <- timeRad[chec$species == "frena"]
densityPlot(frenata, rug=TRUE, ylim=c(0,0.22), ylab="Density of activity", cex.axis=1.1, cex.lab=1.2, extend = NULL)
abline(v=6, lty="dashed", col="black", lwd=2); abline(v=18, lty="dashed", col="black", lwd=2) 

levels(chec$species)
table(chec$species, chec$forest)

#
#
#
#coeficiente de solapamiento----

#C. goudotti ----
#overlap calculation 
bosq <- timeRad[chec$species == "goud" & chec$forest == "andean"]
plant <- timeRad[chec$species == "goud" & chec$forest == "reforestation"]
min(length(bosq), length(plant))
BosqPlantMILLE <- overlapEst(bosq, plant)
BosqPlantMILLE

#plot activity overlap
#png("c_goudotti1.png", height = 10 , width = 15, units = 'cm', res = 300)
par(mar=c(2.8,4.1,1,1))
torq <- overlapPlot(bosq, plant, ylim=c(0,0.29), ylab="Density of activity", cex.axis=1.4, cex.lab=1.4, cex.main= 1.6, xlab = '',
                    rug = T, linecol = c("black", "black"), linewidth = c(1,1), main='', extend = NULL)
a <- readPNG(paste(getwd(), "/c_goudotti.png", sep = ""))
rasterImage(a, xleft=18.2, xright=24.9, ybottom=0.155, ytop=0.27)
abline(v=6, lty="dashed", col="black", lwd=1); abline(v=18, lty="dashed", col="black", lwd=1)
text(21.5,0.125,expression(n == 53), cex = 0.95);text(21.5,0.148,expression(paste(italic("C. goudotti"))), cex = 0.95)
text(12,0.28, expression(hat(Delta) == 0.71  (0.58 - 0.96)), cex = 1.2)
#dev.off()
#legend("topleft", c("Andean ", "G. ruficapilla"), lty=c(1,2), col=c(1,4), bty="n")
#bootstrap
bosqboot <- resample(bosq, 10000)
dim(bosqboot)
plantboot <- resample(plant, 10000)
dim(plantboot)

#overlap estimates by boostrap
BosqPlant <- bootEst(bosqboot, plantboot, adjust = c(1, NA, NA))
dim(BosqPlant)
BSmean <- colMeans(BosqPlant)
BSmean

#calculation of CI
tmp <- BosqPlant[, 1] 
bootCI(BosqPlantMILLE[1], tmp)

tmp <- BosqPlant[, 1] # Extract the required column of the matrix
bootCIlogit(BosqPlantMILLE[1], tmp)




#Z. frenata ----
#overlap calculation 
bosq <- timeRad[chec$species == "frena" & chec$forest == "andean"]
plant <- timeRad[chec$species == "frena" & chec$forest == "reforestation"]
min(length(bosq), length(plant))
BosqPlantMILLE <- overlapEst(bosq, plant)
BosqPlantMILLE

#plot activity overlap
#png("z_frenata1.png", height = 10 , width = 15, units = 'cm', res = 300)
par(mar=c(2.8,4.1,1,1))
torq <- overlapPlot(bosq, plant, ylim=c(0,0.29), ylab="", cex.axis=1.4, cex.lab=1.4, cex.main= 1.4, xlab = '',
                    rug = T, linecol = c("black", "black"), linewidth = c(1,1), main = "", extend = NULL)
a <- readPNG(paste(getwd(), "/z_frenata.png", sep = ""))
rasterImage(a, xleft=19, xright=24, ybottom=0.155, ytop=0.265)
abline(v=6, lty="dashed", col="black", lwd=1); abline(v=18, lty="dashed", col="black", lwd=1)
text(21.5,0.12,expression(n == 87), cex = 0.95);text(21.5,0.143,expression(paste(italic("Z. frenata"))), cex = 0.95)
text(12,0.28,expression(hat(Delta) == 0.75  (0.69 - 0.98)), cex = 1.2)
#dev.off()
#legend("topleft", c("Andean ", "G. ruficapilla"), lty=c(1,2), col=c(1,4), bty="n")
#bootstrap
bosqboot <- resample(bosq, 10000)
dim(bosqboot)
plantboot <- resample(plant, 10000)
dim(plantboot)

#overlap estimates by boostrap
BosqPlant <- bootEst(bosqboot, plantboot, adjust = c(1, NA, NA))
dim(BosqPlant)
BSmean <- colMeans(BosqPlant)
BSmean

#calculation of CI
tmp <- BosqPlant[, 1] 
bootCI(BosqPlantMILLE[1], tmp)

tmp <- BosqPlant[, 1] # Extract the required column of the matrix
bootCIlogit(BosqPlantMILLE[1], tmp)




#A. brunneinucha ----
#overlap calculation 
bosq <- timeRad[chec$species == "brun" & chec$forest == "andean"]
plant <- timeRad[chec$species == "brun" & chec$forest == "reforestation"]
min(length(bosq), length(plant))
BosqPlantMILLE <- overlapEst(bosq, plant)
BosqPlantMILLE


#plot activity overlap
#png('a_brunneinucha1.png', width = 15, height = 10, units = "cm", res = 300)
par(mar=c(2.8,4.1,1,1))
torq <- overlapPlot(bosq, plant, ylim=c(0,0.29), ylab="Density of activity", cex.axis=1.4, cex.lab=1.4, cex.main= 1.4, xlab = '',
                    rug = T, linecol = c("black", "black"), linewidth = c(1,1), main='', extend = NULL)
a <- readPNG(paste(getwd(), "/a_brunneinucha.png", sep = ""))
rasterImage(a, xleft=18.7, xright=24.4, ybottom=0.21, ytop=0.28)
abline(v=6, lty="dashed", col="black", lwd=1); abline(v=18, lty="dashed", col="black", lwd=1)
text(21.5,0.175,expression(n == 49), cex = 0.95);text(21.5,0.20,expression(paste(italic("A. brunneinucha"))), cex = 0.95)
text(12,0.28,expression(hat(Delta) == 0.59  (0.39 - 0.76)), cex = 1.2)
#dev.off()
#legend("topleft", c("Andean ", "G. ruficapilla"), lty=c(1,2), col=c(1,4), bty="n")
#bootstrap
bosqboot <- resample(bosq, 10000)
dim(bosqboot)
plantboot <- resample(plant, 10000)
dim(plantboot)

#overlap estimates by boostrap
BosqPlant <- bootEst(bosqboot, plantboot, adjust = c(1, NA, NA))
dim(BosqPlant)
BSmean <- colMeans(BosqPlant)
BSmean

#calculation of CI
tmp <- BosqPlant[, 1] 
bootCI(BosqPlantMILLE[1], tmp)

tmp <- BosqPlant[, 1] # Extract the required column of the matrix
bootCIlogit(BosqPlantMILLE[1], tmp)



#A. torquatus ----
#overlap calculation 
bosq <- timeRad[chec$species == "torq" & chec$forest == "andean"]
plant <- timeRad[chec$species == "torq" & chec$forest == "reforestation"]
min(length(bosq), length(plant))
BosqPlantMILLE <- overlapEst(bosq, plant)
BosqPlantMILLE


#plot activity overlap
#png('a_torquatus1.png', width = 15, height = 10, units = 'cm', res = 300)
par(mar=c(2.8,4.1,1,1))
torq <- overlapPlot(bosq, plant, ylim=c(0,0.29), ylab="", cex.axis=1.4, cex.lab=1.4, cex.main= 1.4, xlab = '',
                    rug = T, linecol = c("black", "black"), linewidth = c(1,1), main='', extend = NULL)
a <- readPNG(paste(getwd(), "/a_torquatus.png", sep = ""))
rasterImage(a, xleft=18.3, xright=24.8, ybottom=0.21, ytop=0.28)
abline(v=6, lty="dashed", col="black", lwd=1); abline(v=18, lty="dashed", col="black", lwd=1)
text(21.5,0.175,expression(n == 33), cex = 0.95);text(21.5,0.20,expression(paste(italic("A. torquatus"))), cex = 0.95)
text(12,0.28,expression(hat(Delta) == 0.63  (0.53 - 0.94)), cex = 1.2)
#dev.off()
#legend("topleft", c("Andean ", "G. ruficapilla"), lty=c(1,2), col=c(1,4), bty="n")

#bootstrap
bosqboot <- resample(bosq, 10000)
dim(bosqboot)
plantboot <- resample(plant, 10000)
dim(plantboot)

#overlap estimates by boostrap
BosqPlant <- bootEst(bosqboot, plantboot, adjust = c(1, NA, NA))
dim(BosqPlant)
BSmean <- colMeans(BosqPlant)
BSmean

#calculation of CI
tmp <- BosqPlant[, 1] 
bootCI(BosqPlantMILLE[1], tmp)

tmp <- BosqPlant[, 1] # Extract the required column of the matrix
bootCIlogit(BosqPlantMILLE[1], tmp)



#G. ruficapilla ----
#overlap calculation 
bosq <- timeRad[chec$species == "rufi"]
#plant <- timeRad[chec$species == "mille" & chec$forest == 2]
#min(length(bosq), length(plant))
#BosqPlantMILLE <- overlapEst(bosq, plant)
#BosqPlantMILLE

#plot activity overlap
#png('g_ruficapilla1.png', height = 10, width = 15, units = 'cm', res = 300)
par(mar=c(4,4.1,1,1))
rufi <- densityPlot(bosq, ylim=c(0,0.195), ylab="Density of activity", cex.axis=1.4, cex.lab=1.4, cex.main= 1.4,
                     rug = T, linecol = c("black", "black"), linewidth = c(1,1), main='', extend = NULL)
abline(v=6, lty="dashed", col="black", lwd=1); abline(v=18, lty="dashed", col="black", lwd=1)
a <- readPNG(paste(getwd(), "/g_ruficapilla.png", sep = ""))
rasterImage(a, xleft=18.8, xright=24, ybottom=0.10, ytop=0.18)
text(21.6,0.075,expression(n == 39), cex = 0.95);text(21.6,0.09,expression(paste(italic("G. ruficapilla"))), cex = 0.95)
#dev.off()



#G. milleri ----
#overlap calculation 
bosq <- timeRad[chec$species == "mille"]
#plant <- timeRad[chec$species == "mille" & chec$forest == 2]
#min(length(bosq), length(plant))
#BosqPlantMILLE <- overlapEst(bosq, plant)
#BosqPlantMILLE

#plot activity overlap
#png('g_milleri1.png', width = 15, height = 10, units = 'cm', res = 300)
par(mar=c(4,4,1,1))
mille <- densityPlot(bosq, ylim=c(0,0.195), ylab="", cex.axis=1.4, cex.lab=1.4, cex.main= 1.4,
            rug = T, linecol = c("black", "black"), linewidth = c(1,1), main='', extend = NULL)
abline(v=6, lty="dashed", col="black", lwd=1); abline(v=18, lty="dashed", col="black", lwd=1)
a <- readPNG(paste(getwd(), "/g_milleri.png", sep = ""))
rasterImage(a, xleft=18.4, xright=24.4, ybottom=0.10, ytop=0.18)
text(21.5,0.075,expression(n == 28), cex = 0.95);text(21.5,0.09,expression(paste(italic("G. milleri"))), cex = 0.95)
#dev.off()

