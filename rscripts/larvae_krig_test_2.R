############################  Fit All Krigging Models - SGSL Icthyoplankton Survey   ##############################################
##  Andrew Smith & Jean Martin Chamberland
##  Written on Jan 31, 2020
##  
## Input: 1) dataframe to be transformed to spatial points dataframe with at minimum, year, value of interest, and lat long. 
##        2) Grid on which to interpolated the data
## Output: 1) models.xval.r2 - Fitted Krig models & second version which removes outliers above 95% quantile if fit is bad
##         2) results - best fits of the above 
##  Functions used in this script:
##  - sourced "moy.var.Krigeage.R" - provided by Hugo Bourdages
## -  gstat.xv.r2 - calculates R2 from fitted krigged models
## To Do: Turn loops into functions
## Update... feb 13 2020 cleaned and wrangled laraval data from scratch
############################  Load and clean data   ##############################################

# Load packages
library(sf)
library(tidyverse)
library(sp)
library(gstat)
library(spacetime)
library(raster)
library(mapview)

# load data and interpolation grid
load("./data/larvae/df.Rdata")
load("./data/larvae/grid.Rdata")
source("./Rscripts/moy.var.Krigeage.R")


df$lat <- ifelse(df$lat <45,48.4583,df$lat)

# Sensitivity run cut off at longitude -63
df2<- df %>% dplyr::filter(lon < -63)
save(df2, file = "./data/larvae/df2.Rdata")
df<-df2
grid_sens <- grid %>% dplyr::filter(lon < -63)
grid <- grid_sens
# create spatial objects
# cap_larvae$lon <- -cap_larvae$lon
xy<- subset(df, select = c("lon","lat"))
names(xy)<- c("lon","lat")
cap_larvae_sp <- SpatialPointsDataFrame(coords = xy, data = df, proj4string = CRS("+init=epsg:4326")) # wgs84 xy in deg.dec
plot(cap_larvae_sp)
# transform to projection of choice
cap_larvae_sp <- spTransform(cap_larvae_sp, CRSobj = CRS("+init=epsg:32196")) # now in xy in metres quebec lambert nad83

# do same for grid
xy <- subset(grid, select=c("lon", "lat"))
names(xy) <- c("lon", "lat")
grid.sp <- SpatialPointsDataFrame(coords=xy, data = grid, proj4string = CRS("+init=epsg:4326") )
grid.sp <- spTransform(grid.sp, CRSobj = CRS("+init=epsg:32196"))
plot(grid.sp, pch=1, cex=0.1, col="blue")
plot(cap_larvae_sp, pch = 21, col="black", add=T)

# calculate dist between stations
maxDist <- max(dist(coordinates(df_year)))
maxDist <- max(dist(coordinates(cap_larvae_sp)))
# create list of all krigging model types
model.list <- list(nug = vgm(model = "Nug"), 
                   sph = vgm(model = "Sph"),
                   exp = vgm(model = "Exp"),
                   sph.nug = vgm(model = "Sph", psill = NA, range = NA, nugget = NA),
                   exp.nug = vgm(model = "Exp", psill = NA, range = NA, nugget = NA)
                   
)

# function to calculate R2
gstat.xv.r2 <- function(cv.o){
  SSres <- sum(cv.o@data$residual^2)
  SStot <- sum( (cv.o@data$observed - mean(cv.o@data$observed))^2)
  rsquared <- 1-(SSres/SStot)
  return(rsquared)
}

# df to store results
models.xval.r2  <- data.frame(year=NA, model=NA, psill.nug=NA, range.nug=NA, psill.S.E = NA, range.S.E=NA, xval.r2 =NA)[0,]
compteur=1

# loop to run through each year and compare different krigging models 
# Note: this takes ~ 5-10 minutes to run and could be improved
for(i in unique(cap_larvae_sp$year)){
# Create subset and change year 
df_year <- cap_larvae_sp[cap_larvae_sp$year==i,]


# compute variogram
varioCloud <- variogram(Nm2 ~ 1,
                        data = df_year,
                        cloud= F,
                        cutoff = maxDist)
# Plot the cloud variogram
#plot(varioCloud$dist, varioCloud$gamma,     xlab = "Distance (m)",     ylab = "Variance",     las = 1)

# fit variograms:
nug  <- fit.variogram(varioCloud, model = model.list$nug, fit.method = 7)
sph  <- fit.variogram(varioCloud, model = model.list$sph, fit.method = 7)
exp  <- fit.variogram(varioCloud, model = model.list$exp, fit.method = 7)
sph.nug  <- fit.variogram(varioCloud, model = model.list$sph.nug, fit.method = 7)
exp.nug  <- fit.variogram(varioCloud, model = model.list$exp.nug, fit.method = 7)

#xval for each models: 
# Nug
cv.o <- krige.cv(Nm2~1, df_year, nug, nmin = 3, nmax = 15, nfold = nrow(df_year))
# store results
models.xval.r2[compteur, ] <- c(i, "Nug", nug$psill[1], nug$range[1], nug$psill[2], nug$range[2], r2 = gstat.xv.r2(cv.o))
compteur = compteur+1
# sph
cv.o <- krige.cv(Nm2~1, df_year, sph, nmin = 3, nmax = 15, nfold = nrow(df_year))
# store results
models.xval.r2[compteur, ] <- c(i, "Sph", sph$psill[1], sph$range[1], sph$psill[2], sph$range[2], r2 = gstat.xv.r2(cv.o))
compteur = compteur+1
# exp
cv.o <- krige.cv(Nm2~1, df_year, exp, nmin = 3, nmax = 15, nfold = nrow(df_year))
# store results
models.xval.r2[compteur, ] <- c(i, "Exp", exp$psill[1], exp$range[1], exp$psill[2], exp$range[2], r2 = gstat.xv.r2(cv.o))
compteur = compteur+1
# sph.nug
cv.o <- krige.cv(Nm2~1, df_year, sph.nug, nmin = 3, nmax = 15, nfold = nrow(df_year))
# store results
models.xval.r2[compteur, ] <- c(i, "Sph.nug", sph.nug$psill[1], sph.nug$range[1], sph.nug$psill[2], sph.nug$range[2], r2 = gstat.xv.r2(cv.o))
compteur = compteur+1
# exp.nug
cv.o <- krige.cv(Nm2~1, df_year, exp.nug, nmin = 3, nmax = 15, nfold = nrow(df_year))
# store results
models.xval.r2[compteur, ] <- c(i, "Exp.nug", exp.nug$psill[1], exp.nug$range[1], exp.nug$psill[2], exp.nug$range[2], r2 = gstat.xv.r2(cv.o))
compteur = compteur+1


if(nrow(df_year[df_year$Nm2 < quantile(df_year$Nm2, .95), ])>1 & sum(df_year[df_year$Nm2 < quantile(df_year$Nm2, .95), ]$Nm2)>0)
{
# Remove outliers if variogram noisy
varioCloud.nooutlier <- variogram(Nm2 ~ 1,
                        data = df_year[df_year$Nm2 < quantile(df_year$Nm2, .95), ],
                        cloud= F,
                        cutoff = maxDist)

# Normalise the variogram by dividing the data used by the variance if it is determined to be noisy
varioCloud.nooutlier$gamma <- varioCloud.nooutlier$gamma / var(df_year[df_year$Nm2<quantile(df_year$Nm2,.95),]$Nm2)
#plot(varioCloud.nooutlier$dist, varioCloud.nooutlier$gamma,     xlab = "Distance (m)",     ylab = "Variance",     las = 1)

# fit variograms for datasets with no outliers
nug2  <- fit.variogram(varioCloud.nooutlier, model = model.list$nug, fit.method = 7)
sph2  <- fit.variogram(varioCloud.nooutlier, model = model.list$sph, fit.method = 7)
exp2  <- fit.variogram(varioCloud.nooutlier, model = model.list$exp, fit.method = 7)
sph.nug2  <- fit.variogram(varioCloud.nooutlier, model = model.list$sph.nug, fit.method = 7)
exp.nug2  <- fit.variogram(varioCloud.nooutlier, model = model.list$exp.nug, fit.method = 7)

# rescale models to true variance in order to get proper estimate of variance in kriging
sph2$psill  <- sph2$psill * var(df_year$Nm2)
exp2$psill  <- exp2$psill * var(df_year$Nm2)
sph.nug2$psill  <- sph.nug2$psill * var(df_year$Nm2)
exp.nug2$psill  <- exp.nug2$psill * var(df_year$Nm2)


  # same but for 2nd set of models (without extreme noisy values)
  # nug
  cv.o <- krige.cv(Nm2~1, df_year, nug2, nmin = 3, nmax = 15, nfold = nrow(df_year))
  # store results
  models.xval.r2[compteur, ] <- c(i, "Nug2", nug2$psill[1], nug2$range[1], nug2$psill[2], nug2$range[2], r2 = gstat.xv.r2(cv.o))
  compteur = compteur+1
  # sph
  cv.o <- krige.cv(Nm2~1, df_year, sph2, nmin = 3, nmax = 15, nfold = nrow(df_year))
  # store results
  models.xval.r2[compteur, ] <- c(i, "Sph2", sph2$psill[1], sph2$range[1], sph2$psill[2], sph2$range[2], r2 = gstat.xv.r2(cv.o))
  compteur = compteur+1
  # exp
  cv.o <- krige.cv(Nm2~1, df_year, exp2, nmin = 3, nmax = 15, nfold = nrow(df_year))
  # store results
  models.xval.r2[compteur, ] <- c(i, "Exp2", exp2$psill[1], exp2$range[1], exp2$psill[2], exp2$range[2], r2 = gstat.xv.r2(cv.o))
  compteur = compteur+1
  # sph.nug
  cv.o <- krige.cv(Nm2~1, df_year, sph.nug2, nmin = 3, nmax = 15, nfold = nrow(df_year))
  # store results
  models.xval.r2[compteur, ] <- c(i, "Sph.nug2", sph.nug2$psill[1], sph.nug2$range[1], sph.nug2$psill[2], sph.nug2$range[2], r2 = gstat.xv.r2(cv.o))
  compteur = compteur+1
  # exp.nug
  cv.o <- krige.cv(Nm2~1, df_year, exp.nug2, nmin = 3, nmax = 15, nfold = nrow(df_year))
  # store results
  models.xval.r2[compteur, ] <- c(i, "Exp.nug2", exp.nug2$psill[1], exp.nug2$range[1], exp.nug2$psill[2], exp.nug2$range[2], r2 = gstat.xv.r2(cv.o))
  compteur = compteur+1

  } 
}
# models.xval.r2_new <- models.xval.r2
models.xval.r2_sens <- models.xval.r2

# write.csv(models.xval.r2, file = "data/larvae/models.xval.r2.csv")
write.csv(models.xval.r2_new, file = "data/larvae/models.xval.r2_new.csv")
write.csv(models.xval.r2_sens, file = "data/larvae/models.xval.r2_sens.csv")

# models.xval.r2 <- read.csv(file.choose(), row.names = 1)
# models.xval.r2 <- models.xval.r2_new
models.xval.r2 <- models.xval.r2_sens


models.xval.r2$xval.r2 <-as.numeric(as.character(models.xval.r2$xval.r2))
models.xval.r2$psill.nug <-as.numeric(as.character(models.xval.r2$psill.nug))
models.xval.r2$range.nug <-as.numeric(as.character(models.xval.r2$range.nug))
models.xval.r2$psill.S.E <-as.numeric(as.character(models.xval.r2$psill.S.E))
models.xval.r2$range.S.E <-as.numeric(as.character(models.xval.r2$range.S.E))
str(models.xval.r2)
hist(models.xval.r2$xval.r2)
models.xval.r2 %>% ggplot(aes(year, xval.r2, colour = model ))+geom_point(size = 4) + geom_line(aes(group = model)) + facet_wrap(vars(model))

# plots
par(mfrow=c(2,1))
plot(x=models.xval.r2$year, y = models.xval.r2$moy.krig)# the index
lines(x=models.xval.r2$year, y = models.xval.r2$moy.krig) 
plot(x=models.xval.r2$year, y = models.xval.r2$xval.r2) # how well is my data spatially structered? if r2 high then more structured that given year. if r2 negative it is due to low number of points or absence of spatial structure eg 1988 where only station has presence of capelin
lines(x=models.xval.r2$year, y = models.xval.r2$xval.r2)
abline(h=0, lty=2)
par(mfrow=c(1,1))
plot(x=results$xval.r2, y = models.xval.r2$moy.krig)

# Now that we have generated multiple fitted models for a given year, choose only the ones with the highest R2 per year

#  df to store data: 
results <- models.xval.r2
results$moy.krig <- NA
results$var.krig <- NA
results <- results[0,]

compteur <- 1

# moy et variance krigee: 
for(i in unique(cap_larvae_sp$year)){
  df_year <- cap_larvae_sp[cap_larvae_sp$year==i,]
  
  # select best model for given year
  best.model <- models.xval.r2[models.xval.r2$year==i,]
  best.model <- best.model[which(best.model$xval.r2==max(best.model$xval.r2,na.rm=T))[1],]
  
  # store info into df: 
  results[compteur,c(1:7)] <- best.model
  
  #  create variogram object by replacing parameters in existing variogram model object.
  if(sum(is.na(best.model[,c("psill.S.E","range.S.E")]))==2)
  {
    vario.obj <- nug
    vario.obj$psill <- best.model$psill.nug
    vario.obj$range <- best.model$range.nug
  }else{
    vario.obj <- sph.nug
    vario.obj$model <- c("Nug", substr(best.model$model,1,3)) 
    vario.obj$psill <- c(best.model$psill.nug, best.model$psill.S.E)
    vario.obj$range <- c(best.model$range.nug, best.model$range.S.E)
  }

  # moy var krig
  moy_var <- moy.var.Krigeage(xCoords = coordinates(df_year), xValues = df_year$Nm2, kCoords = coordinates(grid.sp), vario = vario.obj)
  
  results[compteur,]$moy.krig <- moy_var$Z
  results[compteur,]$var.krig <- moy_var$VarKrigeage
  compteur <-compteur+1 
}
results_sens <- results
results_new <- results

write.csv(results, file = "./results.capelin.larvae.kriging.csv", row.names = F)
write.csv(results_new, file = "./results.capelin.larvae.kriging_new.csv", row.names = F)
write.csv(results_sens, file = "./results.capelin.larvae.kriging_sens.csv", row.names = F)

summary(results_new$xval.r2)
par(mfrow=c(2,1))
plot(x=results$year, y = results$moy.krig)# the index
lines(x=results$year, y = results$moy.krig) 
plot(x=results$year, y = results$xval.r2) # how well is my data spatially structered? if r2 high then more structured that given year. if r2 negative it is due to low number of points or absence of spatial structure eg 1988 where only station has presence of capelin
lines(x=results$year, y = results$xval.r2)
abline(h=0, lty=2)
par(mfrow=c(1,1))
plot(x=results$xval.r2, y = results$moy.krig)

# explain negative rsquared: 1988, 1993, 1994, 2006 with data exploration: 
df %>% ggplot(aes(lon, lat, colour = presence, size = Nm2))+geom_point()+ facet_wrap(vars(year)) + theme_classic(base_size = 18)

# 1988 has 1 point with presence of capelin
df %>% ggplot(aes(log(Nm2)))+geom_density()+ facet_wrap(vars(year))


hist(results$xval.r2)
h = hist(results$xval.r2) # or hist(x,plot=FALSE) to avoid the plot of the histogram
h$density = h$counts/sum(h$counts)*100
plot(h,freq=FALSE)  # en y = pourcentage des r2 xval du meillleru model de krig de l<annee
# so it looks like roughly 60% of all years have very weak spatial structure... this justifies the use of a rough spatial model "trend surface analysis" ie adding 
# go back to neg binom zeroinfl and reproject coordinates to proper space and then put back in data frame to run the model 

hist(results$xval.r2, freq=F, ylim=c(0,1))

# figures
results <- read_csv("data/larvae/results.capelin.larvae.kriging.csv")
par(mfrow(c(2,1)))
results_capelin_larvae_kriging %>% ggplot(aes(year, lat, colour = presence, size = Nm2))+geom_point()+ facet_wrap(vars(year)) + theme_classic(base_size = 18)

