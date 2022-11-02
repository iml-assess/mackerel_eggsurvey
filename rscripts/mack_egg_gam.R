library(gratia)
library(mgcv)
library(gamair)
library(tidyverse)
library(magrittr)
library(broom)
data("mack", "mackp")
data("med")
data("meh")
data("coast")


plot(med$lo,med$la,cex=0.2+med$count^.5/10,col="red",
     pch=21,xlab="lo",ylab="la",main="mackerel")
ind <- med$count==0
points(med$lo[ind],med$la[ind],cex=0.1,pch=19)
lines(coast)



mack$log.net.area <- log(mack$net.area)
gmtw <- gam(egg.count ~ s(lon,lat,k=100) + s(I(b.depth^.5))+
                s(c.dist) + s(salinity) + s(temp.surf) + s(temp.20m)+
                offset(log.net.area),data=mack,family=tw,method="REML",
            select=TRUE)

gm2 <- gam(egg.count ~ s(lon,lat,k=100) + s(I(b.depth^.5)) +
               s(c.dist) + s(temp.20m) + offset(log.net.area),
           data=mack,family=tw,method="REML")

mackp$log.net.area <- rep(0,nrow(mackp))
lon <- seq(-15,-1,1/4); lat <- seq(44,58,1/4)
zz<-array(NA,57*57); zz[mackp$area.index]<-predict(gm2,mackp)
image(lon,lat,matrix(zz,57,57),col=gray(0:32/32),
      cex.lab=1.5,cex.axis=1.4)
contour(lon,lat,matrix(zz,57,57),add=TRUE)
lines(coast$lon,coast$lat,col=1)

set.seed(4) ## make reproducible
br1 <- rmvn(n=1000,coef(gm2),vcov(gm2))
Xp <- predict(gm2,newdata=mackp,type="lpmatrix")
mean.eggs1 <- colMeans(exp(Xp%*%t(br1)))
hist(mean.eggs1)



data(med); head(med) ## look at data
data(coast)
## initial plots...
plot(med$lo,med$la,cex=0.2+med$count^.5/10,col="red",
     pch=21,xlab="lo",ylab="la",main="mackerel")
ind <- med$count==0
points(med$lo[ind],med$la[ind],cex=0.1,pch=19)
lines(coast)
## ... survey seems to cover spawning area this time!
require(mgcv)
m1 <- gam(count~s(lo,la,k=100)+s(T.surf)+s(T.20)+s(I(b.depth^.5))+s(Sal20)+
              s(ship,bs="re")+offset(log(vol)),data=med,select=TRUE,family=tw)
gam.check(m1) ## mean variance relationship not quite right?
m2 <- gam(count~s(lo,la,k=100)+s(T.surf)+s(T.20)+s(I(b.depth^.5))+s(Sal20)+
              s(ship,bs="re")+offset(log(vol)),data=med,select=TRUE,family=nb)
gam.check(m2)
par(mfrow=c(1,2)) ## re-check residuals vs fitted
plot(fitted(m1)^.5,residuals(m1));plot(fitted(m2)^.5,residuals(m2))
AIC(m1,m2) ## neg bin much better
plot(m2,pages=1) ## effects


draw(m1)
draw(m2)
appraise(m1)
appraise(m2)


Prediction=predict.gam(m1,mackp,type='response')

ggplot(aes(lon,lat), data = mackp) + geom_point() +
    geom_point(aes(lon,lat, colour = "red"),data = mack)


ggplot(aes(lon,lat, colour = ship), data = med) + geom_point() 
names(med)
names(mackp)

mackp$vol <- rep(0,nrow(mackp))
med %<>% mutate(lon = lo, lat = la, egg.count = count, temp.surf = T.surf, temp.20m = T.20, salinity = Sal20 )

m3 <- gam(egg.count ~ 
              s(lon, lat, k=100) +
              s(temp.20m) +
              s(I(b.depth^.5)) +
              s(salinity) + 
              offset(log(vol)),
          data = med,
          select = TRUE,
          family = tw)

med %>% ggplot(aes(lon,lat))+geom_point() 
mackp %>% ggplot(aes(lon,lat))+geom_point() 
fitm3 <- predict(m3)
ffitm3<-predict.gam(m3)

df<-augment(m3, type.predict = "response")
df$predicted <- df$.fitted
df %>% ggplot(aes(egg.count,predicted)) + geom_point()


names(med)
library(lubridate)
med %<>% mutate(time = hour(DT), day = yday(DT))
medsub <- med %>% dplyr::filter(count<3000)

# null tweedie model
m0 <- gam(egg.count ~ 
              s(lon, lat, k=100) +
              offset(log(vol)),
          data = med,
          select = TRUE,
          family = tw)

m1 <- gam(egg.count ~ 
              s(lon, lat, k=100) +
              s(temp.20m, k = 3) +
              s(I(b.depth^.5)) +
              s(salinity, k = 3) + 
              s(ship,bs="re") +
              s(time, bs = "cc", k = 12) +
              s(day, bs = "cc", k = 3) +
              offset(log(vol)),
          data = med,
          select = TRUE,
          family = tw)

# subsetted tweedie
m1s <- gam(egg.count ~ 
              s(lon, lat, k=100) +
              s(temp.20m) +
              s(I(b.depth^.5)) +
              s(salinity) + 
              s(ship,bs="re") +
              s(time, bs = "cc") +
              s(day, bs = "cc") +
              offset(log(vol)),
          data = medsub,
          select = TRUE,
          family = tw)

# full negbin
m2 <- gam(egg.count ~ 
              s(lon, lat, k=100) +
              s(temp.20m) +
              s(I(b.depth^.5)) +
              s(salinity) + 
              s(ship,bs="re") +
              s(time, bs = "cc") +
              s(day, bs = "cc") +
              offset(log(vol)),
          data = med,
          select = TRUE,
          family = nb)

# subset negbin
m2s <- gam(egg.count ~ 
               s(lon, lat, k=100) +
               s(temp.20m) +
               s(I(b.depth^.5)) +
               s(salinity) + 
               s(ship,bs="re") +
               s(time, bs = "cc") +
               s(day, bs = "cc") +
               offset(log(vol)),
           data = medsub,
           select = TRUE,
           family = nb)

m3 <- gam(egg.count ~ 
               s(lon, lat, k=100) +
               s(temp.20m) +
               s(ship,bs="re") +
               s(day, bs = "cc") +
               offset(log(vol)),
           data = medsub,
           select = TRUE,
           family = nb)

m4 <- gam(egg.count ~ 
              s(lon, lat, k=100) +
              s(temp.20m) +
              s(ship,bs="re") +
              s(day, bs = "cc") +
              offset(log(vol)),
          data = medsub,
          select = TRUE,
          family = tw)

gam.check(m0)
gam.check(m1)
gam.check(m1s)
gam.check(m2)
gam.check(m2s)
gam.check(m3)
gam.check(m4)
summary(m1)
summary(m1s)
summary(m2)
summary(m2s)
summary(m3)
summary(m4)
par(mfrow=c(2,2)) ## re-check residuals vs fitted
plot(fitted(m1)^.5,residuals(m1));plot(fitted(m1s)^.5,residuals(m1s));plot(fitted(m2)^.5,residuals(m2));plot(fitted(m2s)^.5,residuals(m2s))
AIC(m0,m1,m1s,m2,m2s,m3,m4) ## neg bin much better
plot(m2,pages=1) ## effects


draw(m1)
draw(m1s)
draw(m3)
appraise(m1)
appraise(m1s)
appraise(m3)
draw(m2)
draw(m2s)
appraise(m2)
appraise(m2s)

performance(m1,m1s,m2,m2s)



# add vars to prediction grid (smaller but try)
med_pred <- med %>% dplyr::select(count, lon, lat, ship, vol, time, day, salinity, b.depth, temp.20m)
# med_pred$vol <- rep(0,nrow(med_pred))
med_pred$ship <- "29CS"
med_pred$day <- 80
med_pred$time <- 12
med_pred$area.index <- 1
med_pred %<>% dplyr::filter(count<3000)

med_pred2 <- med_pred %>% mutate(lat = round(lat,0), lon = round(lon,0)) %>% 
    group_by(lat, lon, ship) %>% 
    dplyr::summarise(count = mean(count), 
                     vol = mean(vol),
                     time = mean(time),
                     day = mean(day),
                     salinity = mean(salinity, na.rm = T),
                     b.depth = mean(b.depth),
                     temp.20m = mean(temp.20m, na.rm = T)) %>% 
    mutate(b.depth = ifelse(is.na(b.depth), 500, b.depth), temp.20m = ifelse(is.na(temp.20m), 13, temp.20m))

df<-augment(m2s, type.predict = "response", newdata = med_pred)
df$pred_std <- df$.fitted

df2<-augment(m2s, type.predict = "response", newdata = med_pred2)
df2$pred_std <- df2$.fitted


df %>% ggplot(aes(lon,lat, size = count)) +
    geom_point(alpha = 0.5) +
    geom_point(aes(lon,lat, size = pred_std), colour = "red", alpha = 0.5)+
    scale_size_area()

df2 %>% ggplot(aes(lon,lat, size = count)) +
    geom_point(alpha = 0.5) +
    geom_point(aes(lon,lat, size = pred_std), colour = "red", alpha = 0.5)+
    scale_size_area()


df2$df <- "grid"
df$df <- "raw"
DF <- bind_rows(df,df2)
DF
DF %>% ggplot(aes(lon,lat)) +
    geom_point(size = 0.1) +
    geom_point(aes(lon,lat, colour = "predicted", size = pred_std), alpha = 0.4, colour = "blue") +
    scale_size_area() +
    facet_wrap(vars(df))

DF %>% ggplot(aes(lon,lat)) +
    geom_point() +
    scale_size_area() +
    facet_wrap(vars(df))


DF %>% group_by(df) %>% dplyr::summarise(m = mean(pred_std, na.rm = T),
                                         v = var(pred_std, na.rm = T))

med %>% dplyr::summarise(m = median(count))

dd<-augment(m3, newdata = medsub, type.predict = "response")
dd
mean(dd$count)
mean(dd$.fitted,na.rm = T)


####
load("~/Data/Ichthyoplankton/Rdata/sgsl_ichthyo_2021.Rdata")
station_key <- read_csv("~/Data/Ichthyoplankton/data/station_key.csv") # incomplete
sgsl_env <- read_csv("~/Data/Ichthyoplankton/data/sgsl_env.csv")
sgsl_ichthyo_temp <- read_csv("~/Data/Ichthyoplankton/data/sgsl_ichthyo_temp.csv")
grille_pred <- read_delim("~/Data/Ichthyoplankton/data/grille_pred.csv",delim = ";", escape_double = FALSE, trim_ws = TRUE)

# wrangle
sgsl_ichthyo_temp %<>% pivot_longer(cols = 5:40, names_to = "year", values_to = "temp1")
sgsl_env <- left_join(sgsl_env, station_key)
sgsl_env %<>% transmute(year = str_sub(date_utc,8,11),
                       station_raw = station_raw, 
                       station_name = station_name, # incomplete
                       time = str_sub(time, 1,2),
                       lon = longitude, 
                       lat = as.numeric(latitude), 
                       temp2 = temp,
                       sal = sal,
                       sigmaT = as.numeric(sigma))
# subset
target <- c("Scomber scombrus (larva)",
            "Scomber scombrus (oeuf stade 1)",
            "Scomber scombrus (oeuf stade 2)",
            "Scomber scombrus (oeuf stade 3)",
            "Scomber scombrus (oeuf stade 4)",
            "Scomber scombrus (oeuf stade 5)")

dat <- sgsl_ichthyo %>% 
  dplyr::filter(taxons %in% target,
                trajet %in% c("T1","RT"),
                year != 2009, 
                pass == 1,
                extra_1991a == 0,
                extra_1991b == 0,
                extra_1986 == 0,
                extra_northumberland == 0,
                extra_ns_1998 == 0,
                extra_st_georges_bay_1994 == 0,
                extra_2000 == 0,
                extra_scotian_shelf_2009 == 0) %>% 
  dplyr::select(year,
                month,
                day,
                doy,
                hour_start,
                stratum,
                lat,
                lon,
                station,
                station_name,
                pass,
                pl_headr_start_depth,
                pl_headr_volume,
                pl_gen_split_fraction,
                pl_gen_counts,
                taxons,
                n, 
                n_m3,
                n_m2) %>% 
  mutate(scalar = pl_headr_volume*pl_headr_start_depth,
         f_year = as.factor(year), 
         lon = as.numeric(lon),
         lat = as.numeric(lat),
         vol = pl_headr_volume,
         depth = pl_headr_start_depth)

dat <- dat[complete.cases(dat), ]

# free up space
rm(sgsl_ichthyo)

glimpse(dat)
dat %>% 
  ggplot(aes(lon,lat, colour = taxons)) +
  geom_point() +
  facet_wrap(vars(year))

ml <- dat %>% dplyr::filter(taxons == "Scomber scombrus (larva)") %>% as_tibble()
m1 <- dat %>% dplyr::filter(taxons == "Scomber scombrus (oeuf stade 1)") %>% as_tibble()
m2 <- dat %>% dplyr::filter(taxons == "Scomber scombrus (oeuf stade 2)") %>% as_tibble()
m3 <- dat %>% dplyr::filter(taxons == "Scomber scombrus (oeuf stade 3)") %>% as_tibble()
m4 <- dat %>% dplyr::filter(taxons == "Scomber scombrus (oeuf stade 4)") %>% as_tibble()
m5 <- dat %>% dplyr::filter(taxons == "Scomber scombrus (oeuf stade 5)") %>% as_tibble()

### models
# larvae
# tweedie vs negative binomial
larvae_null_tw <- gam(n ~ f_year +
                     s(lon, lat, k = 100) +
                     offset(log(vol)),
                   data = ml,
                   select = TRUE,
                   family = tw,
                   method = "REML")

larvae_null_nb <- gam(n ~ f_year +
                        s(lon, lat, k = 100) +
                        offset(log(vol)),
                      data = ml,
                      select = TRUE,
                      family = nb,
                      method = "REML")

# remove outlier otherwise similar performance between tw and nb with tw having lower AIC 9484 vs 9619 and greater deviance explained ~ 66% vs 59% 
ml %<>% dplyr::filter(n < 20000)
ml2019 <- ml %>% dplyr::filter(year == 2019)

# larvae

larvae <- gam(n ~ f_year +
                s(lon, lat, k = 10) +
                s(depth, k = 3) + 
                s(hour_start, bs = "cc", k = 12) + 
                s(doy, bs = "cc", k = 12) + 
                offset(log(vol)),
              data = ml,
              select = TRUE,
              family = tw,
              method = "REML")

library(performance)
library(ggeffects)


plot(DF)
# stage 1 eggs
egg1 <- gam(n ~ f_year +
              s(lon, lat, k = 10) +
              s(depth, k = 3) + 
              s(hour_start, bs = "cc", k = 12) + 
              s(doy, bs = "cc", k = 12) + 
              offset(log(vol)),
            data = m1,
            select = TRUE,
            family = tw,
            method = "REML")

# stage 2 eggs

egg2 <- gam(n ~ f_year +
              s(lon, lat, k = 10) +
              s(depth, k = 3) + 
              s(hour_start, bs = "cc", k = 12) + 
              s(doy, bs = "cc", k = 12) + 
              offset(log(vol)),
            data = m2,
            select = TRUE,
            family = tw,
            method = "REML")

# stage 3 eggs

egg3 <- gam(n ~ f_year +
              s(lon, lat, k = 10) +
              s(depth, k = 3) + 
              s(hour_start, bs = "cc", k = 12) + 
              s(doy, bs = "cc", k = 12) + 
              offset(log(vol)),
            data = m3,
            select = TRUE,
            family = tw,
            method = "REML")

# stage 4 eggs

egg4 <- gam(n ~ f_year +
              s(lon, lat, k = 10) +
              s(depth, k = 3) + 
              s(hour_start, bs = "cc", k = 12) + 
              s(doy, bs = "cc", k = 12) + 
              offset(log(vol)),
            data = m4,
            select = TRUE,
            family = tw,
            method = "REML")

# stage 5 eggs

egg5 <- gam(n ~ f_year +
              s(lon, lat, k = 10) +
              s(depth, k = 3) + 
              s(hour_start, bs = "cc", k = 12) + 
              s(doy, bs = "cc", k = 12) + 
              offset(log(vol)),
            data = m5,
            select = TRUE,
            family = tw,
            method = "REML")


# model diagnostics
AIC(larvae, egg1, egg2, egg3, egg4, egg5)
summary(egg5)
gam.check(larvae)
draw(egg5)
appraise(larvae)
performance(larvae)

larvae_pred_df <- ml %>% dplyr::select(n, f_year, depth, hour_start, doy, vol)
larvae_fitted <- augment(larvae, type.predict = "response") %>% mutate(stage = "larvae")
egg1_fitted <- augment(egg1, type.predict = "response") %>% mutate(stage = "egg1")
egg2_fitted <- augment(egg2, type.predict = "response") %>% mutate(stage = "egg2")
egg3_fitted <- augment(egg3, type.predict = "response") %>% mutate(stage = "egg3")
egg4_fitted <- augment(egg4, type.predict = "response") %>% mutate(stage = "egg4")
egg5_fitted <- augment(egg5, type.predict = "response") %>% mutate(stage = "egg5")

fitted_df <- bind_rows(larvae_fitted, egg1_fitted)
fitted_df <- bind_rows(fitted_df, egg2_fitted)
fitted_df <- bind_rows(fitted_df, egg3_fitted)
fitted_df <- bind_rows(fitted_df, egg4_fitted)
fitted_df <- bind_rows(fitted_df, egg5_fitted)

fitted_df$fitted_n = fitted_df$.fitted

fitted_df %>% ggplot(aes(n, fitted_n, colour = stage)) + geom_point()


fitted_df %>% ggplot(aes(f_year, n, colour = stage)) + geom_point()
fitted_df %>% ggplot(aes(f_year, fitted_n, colour = stage)) + geom_point()
fitted_df %>% ggplot(aes(f_year, fitted_n, colour = stage)) + stat_summary(fun.data = "mean_cl_boot")


index <- fitted_df %>% group_by(f_year, stage) %>% dplyr::summarise(n = mean(fitted_n, na.rm = T))
### Figures
# Maps
# interactive map by year 
canada <- map_data("world", "Canada") %>% mutate(lon = long)
canada_map <- ggplot() + geom_polygon(data = canada,aes(x = lon, y = lat, group = group),fill = "gray20",color = "black") + theme_bw() 
canada_atlantic <- canada_map + 
  coord_map(xlim = c(-69, -57.5),  ylim = c(43.5, 50), projection = "lambert", parameters = c(45,50))       
sGSL <- canada_map + 
  coord_map(xlim = c(-66.5, -59),  ylim = c(45.6, 49), projection = "lambert", parameters = c(46.5 ,48))  
rm(canada);rm(canada_map)

sGSL + 
  geom_point(data = ml %>%
               dplyr::filter(year == 2019),
             aes(lon,lat, size = n_m2, fill = n_m2),
             shape = 21) + 
  scale_fill_viridis_c(option = "H") +
  scale_size_area()

sGSL + 
  geom_point(data = ml,
             aes(lon,lat, size = n_m2, fill = n_m2),
             shape = 21) + 
  scale_fill_viridis_c(option = "H") +
  scale_size_area() +
  facet_wrap(vars(year))
