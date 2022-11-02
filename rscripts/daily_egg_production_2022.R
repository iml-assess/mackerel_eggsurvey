######### Mackerel Egg Index
##
##  Script to subset southern Gulf of St Lawrence spring ichthyoplankton survey data (1979-2020) 
##  and calculate mackerel daily and total egg production for input into mackerel stock assessment model CCAM
##  
##  bY ANDREW D SMITH
##  
##  Change Log:
##  v.0 = pre 2017: all data on annual, differently formatted spreadsheets (.xls or .xlsx) in S/Pélagiques/Plancton/Relevés
##  v.1 = 2017: Baye Mbaye and Elisabeth van Beveren wrote eggSurvey.R but data was messy and only covered years 2012 onwards. published values were used otherwise.
##  v.1.5 = August 2018: Andrew Smith, Jean Martin Chamberland, and Mélanie Boudreau all make attempts to clean and merge the data. The end result, largely rewritten by ADS is SSB_index_egg_survey.R. Still many problems with formatting and bugs. raw data used for quality control and summary statistics. published values used in model with updates from raw for years 2015:2018
##  v.2.0 = November 2020: Following the validation and transfer of the data onto the national BioChem database by the DAISS team at IML (Isabelle St. Pierre and Hélène Dionne), 
##        - ADS writes read_ichthyo.R to merge and format the BioChem data and names the data sgsl_ichthyo_full. Some formatting still needed for station names etc. 1982-1985 hard coded in for the instant (original mission sheets MIA). 
##        - January 2021: ADS writes this script (mackerel_egg_index.R) to read and parse new format
##                          - Added summarised mission metadata as well gsi, krigging, and other varia model params
##                          - Added temperature data from Caroline. Lafleur - Francois Grégoire calculated these a different way as values differ. Decision: years(1979-2014 = published set, 2015-2020 C. Lafleur DAISS)
##                          - Added figures to do quality control
## Last full DFO Research Documents detailing the ichthyoplankton survey were:
#   - Grégoire, F., Girard, L. et Boudreau, M. 2014. Résultats des relevés du programme de monitorage zonal atlantique (PMZA)-maquereau bleu (Scomber scombrus L.) réalisés dans le sud du golfe du Saint-Laurent en 2012 et 2013. Secr. can. de consult. sci. du MPO. Doc. de rech. 2014/075. v + 82 p.
##  and
##  - Grégoire, F., Girard, L., Beaulieu, J.-L., Lussier, J.-F., et Gendron, M.-H. 2014. Détection des tendances communes dans les abondances d’oeufs et de larves de poissons récoltés dans le sud du golfe du Saint-Laurent entre 1983 et 2012. Secr. can. de consult. sci. du MPO. Doc. de rech. 2014/074. v + 34 p.
##  and
##  - Grégoire, F., Girard, L., et Beaulieu, J.-L. 2014. Analyses de similarité appliquées sur les abondances de larves de poissons récoltées dans le sud du golfe du Saint-Laurent entre 1983 et 2012. Secr. can. de consult. sci. du MPO. Doc. de rech. 2014/080. v + 16 p.
## 
##  - Update: 23/09/2022: Preparing transfer of code to Caroline Lehoux. Checking if everything runs correctly up to 2021 data
    # - reprojecting to Quebec Lambert does little in the way of analysis. It changes the mean and var of the interpolated values very little (Decimal order) and by keeping the data in WGS84 you have the advantage that creating figures is far easier as lat lon stay in deg.dec
    # - recalculated variograms and fit models for all years instead of using published results. Chose to use only exponential model. except 1982 refit (1982 temperature data not found...). 
    # - comparison with automap::autokrige function yields similar results to the custom function moy.var.krige but the estimated variances are less similar
#########

############################  Load and clean data   ##############################################

#####################  LOAD PACKAGES  ####################
packages = c('sf', 'sp','rgdal','rgeos','readxl','spacetime','automap', 'raster','mapview', 'gstat', 'magrittr', 'lubridate', 'readr', 'tidyverse', 'mapproj')
invisible(lapply(packages, function(x) {
  if (!require(x, character.only = T)) {
    install.packages(x)
    require(x)
  }
}))

source("./functions/moy.var.Krigeage.R") # This function is quite old now and maybe gstat or automap can do it? Results have never been compared
source("./functions/mackerel_fun_incubation.R") # Mackerel incubation eqns.

#####################  LOAD AND WRANGLE DATA  ####################

# Stage 1 and 5 mackerel eggs. Response variable of interest is n_m2, i.e. the number of eggs per square metre after standardizing for station depth and volume
load("./rdata/egg_index_stage15_1979-2022.Rdata")
# Classic interpolation grid
load("./rdata/grid.Rdata")
# Denser grid made in qgis for plots (optional)
# grid_dense <- readOGR(dsn = "C:/Users/SmithAND/Documents/My_Stocks/Mackerel/Assessments/2019/2018_Mackerel_Evaluation/Data/sGSL_grid/sGSL_grid_1NM.shp", layer = "sGSL_grid_1NM")

# Subset and filter data
# Only trajet 1
egg %<>% dplyr::filter(trajet == "T1", 
                       !is.na(lat)) %>% 
  droplevels()

# Check data set dimensions. There should be only one obs per station per year. If this is not the case it is likely due to parsing station names from biochem dataset (e.g. 2.1 cor, 2.1a etc.) -> either delete one or take the mean
temp <- egg %>% count(year, station) 
temp <- egg %>% group_by(year, station) %>% mutate(n = n()) %>% 
  dplyr::filter(n>1)

# Calculate incubation time. A number of equations exist. Could be good sensitivity test to compare Mendiola 2016 eqn to our classic Lockwood 1977
source("./functions/mackerel_fun_incubation.R")

egg %<>% 
  mutate(I = I_lockwood(temperature), # unit = hours
         I2 = I_mendiola(temperature)) # unit = hours


# Calculate daily egg production (DEP) by station
egg$dep = egg$n_m2/egg$I * 24
egg$dep2 = egg$n_m2/egg$I2 * 24
summary(egg$dep)

# no temperature data easily available for 1979, 1982, and none yet for 2022 so drop for now
egg %<>% dplyr::filter(!is.na(dep))

# Specify official year set rules (1979-2020). The specified years are usually not included in the time series due to varied reasons in previous peer reviews. Could try and use them all the same and have it as a sensitivity run / extended model. 
egg %<>% mutate(mainset = ifelse(year %in% c(1982, 1999, 2001, 2006), 0, 1))

#####################  SUMMARY FIGURES  ####################
# Maps
canada <- map_data("world", "Canada") %>% mutate(lon = long)
canada_map <- ggplot() + geom_polygon(data = canada,aes(x = lon, y = lat, group = group),fill = "gray20",color = "black") + theme_bw() 
canada_atlantic <- canada_map + 
  coord_map(xlim = c(-69, -57.5),  ylim = c(43.5, 50), projection = "lambert", parameters = c(45,50))       
sGSL <- canada_map + 
  coord_map(xlim = c(-66.5, -60),  ylim = c(45.6, 49), projection = "lambert", parameters = c(46.5 ,48))  
rm(canada);rm(canada_map);rm(canada_atlantic)

# egg stage maps
library(ggforce)
a <- sGSL +  geom_point(data = egg, aes(lon,lat, size = dep, fill = dep), shape = 21, colour = "black", alpha = 0.8) +
  labs(y = "Latitude", x = "Longitude") +
  theme_minimal() +
  scale_fill_viridis_c(option = "magma", end = 0.85) +
  scale_size_area(max_size = 6) +
  theme(legend.position = "top")
a + facet_wrap_paginate(vars(year), ncol = 3, nrow = 4, page = 1)
a + facet_wrap_paginate(vars(year), ncol = 3, nrow = 4, page = 2)
a + facet_wrap_paginate(vars(year), ncol = 3, nrow = 4, page = 3)

# stratum map
sGSL + geom_point(aes(lon, lat, fill = stratum), size = 6, shape = 21, data = egg) + scale_size_area() +
  theme_minimal(base_size = 14) + labs(x = "Longitude", y = "Latitude") +
  scale_fill_viridis_d(end = 0.9, option = "B") + theme(legend.position = "top") 
  # ggsave("dep_stratum_map.png", device = "png", path = figures_path, width = 6.5, height = 7, units = "in", dpi = 300)   

# station map
sGSL + geom_point(aes(lon,lat),fill = "white", colour ="black", shape = 21, size = 10, data = egg) +
  geom_text(aes(lon, lat, label = station), size = 3.5, colour = "black", data = egg) +
  theme_minimal(base_size = 14) + labs(x = "Longitude", y = "Latitude") +
  scale_fill_viridis_d(end = 0.9, option = "B") + theme(legend.position = "top") 
  # ggsave("station_map.png", device = "png", path = figures_path, width = 6.5, height = 7, units = "in", dpi = 300)   

# annual proportion of dep by temperature and ecdf
temp <- egg %>%
  group_by(year, temperature) %>%
  summarise(n = sum(dep,na.rm=T)) %>%
  group_by(year) %>% 
  mutate(freq = n / sum(n), temp = temperature) %>%
  dplyr::filter(!is.na(temp),!is.na(freq)) %>%
  as_tibble()

temp <- egg %>%
  group_by(year, temperature) %>%
  mutate(n = sum(dep,na.rm=T)) %>%
  group_by(year) %>% 
  dplyr::summarise(freq = n / sum(n), temp = temperature) %>% 
  dplyr::filter(!is.na(temp),!is.na(freq)) %>%
  as_tibble()

temp <- egg %>%
  group_by(year, round(temperature)) %>%
  summarise(n = n(), dep = sum(dep)) %>%
  mutate(prop = n/sum(n), w_prop = prop *dep)

egg %>% group_by(round(temperature,2)) %>% 
  dplyr::summarise(n_m2 = gm_mean(n_m2)) %>% 
  ggplot(aes(`round(temperature, 2)`,n_m2))+geom_point()

egg %>% ggplot(aes(temperature,dep))+ 
  geom_point() + 
  facet_wrap(vars(year), scales = "free_y") +
  xlab("temperature") +
  ylab("daily egg production")+
  labs(title = "Relationship between temperature and daily egg production", subtitle = "note that y axis scale varies")

egg %>% ggplot(aes(round(temperature))) +
  stat_ecdf() + 
  geom_point(aes(`round(temperature)`, prop),data = temp) +
  facet_wrap(vars(round(year))) + xlab("temperature") + ylab("cumulative distribution (ECDF)") +
  labs(title = "Relationship between temperature and daily egg production", subtitle = "points = % annual dep, line = ecdf of temperature")


#####################  CREATE SPATIAL OBJECTS FOR KRIGGING  ####################
egg2 <- egg %>% dplyr::group_by(year, station) %>% 
  dplyr::summarise(lon = mean(round(lon,2)), lat = mean(round(lat,2)),
                   dep = mean(round(dep,2)), temperature = mean(round(temperature,2)))
# for data set
xy <- subset(egg2, select = c("lon","lat"))
names(xy)<- c("lon","lat")
egg_sp <- SpatialPointsDataFrame(coords = xy, data = egg2, proj4string = CRS("+init=epsg:4326")) # wgs84 xy in deg.dec
plot(egg_sp)
plot(SGSL,add=T)

# transform to projection of choice (optional. It is fine in WGS84)
# egg_sp <- spTransform(egg_sp, CRSobj = CRS("+init=epsg:32196")) # now in xy in metres quebec lambert nad83
# CRS_LCC_83 <- CRS("+proj=lcc +lat_0=46.5 +lat_1=48 +lat_2=50 +lon_0=-70 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs") #ESRI:102002
# egg_sp <- spTransform(egg_sp, CRSobj = CRS_LCC_83)

# do same for classic interpolation grid
xy <- subset(grid, select=c("lon", "lat"))
names(xy) <- c("lon", "lat")
grid_sp <- SpatialPointsDataFrame(coords=xy, data = grid, proj4string = CRS("+init=epsg:4326") )
# grid_sp <- spTransform(grid_sp, CRSobj = CRS("+init=epsg:32196"))
# grid_sp <- spTransform(grid_sp, CRSobj = CRS_LCC_83)

# make sure they overlap
plot(grid_sp, pch=1, cex=0.1, col="blue")
plot(egg_sp, pch = 21, col="black", add=T)

# do same for denser interpolation grid (for smoother looking plots but not for main calculation)
# grid_dense_df <- as.data.frame(grid_dense) # to remove the default projection of +longlat wgs83 turn to dataframe (also good to have for plotting)
# grid_dense_df$lon <- grid_dense_df$coords.x1 
# grid_dense_df$lat <- grid_dense_df$coords.x2
# grid_dense_df <- grid_dense_df %>% dplyr::select(id,lon,lat)
# xy <- subset(grid_dense_df, select=c("lon", "lat"))
# names(xy) <- c("lon", "lat")
# grid_dense_sp <- SpatialPointsDataFrame(coords=xy, data = grid_dense_df, proj4string = CRS("+init=epsg:4326") )
# # grid_dense_sp <- spTransform(grid_dense_sp, CRSobj = CRS("+init=epsg:32196"))
# # grid_dense_sp <- spTransform(grid_dense_sp, CRSobj = CRS_LCC_83)
# rm(grid_dense)
# make sure they overlap
# plot(grid_dense_sp, pch=1, cex=0.1, col="blue")
# plot(egg_sp, pch = 21, col="black", add=T)

# VARIOGRAM BY YEAR (This could and should be turned into a function or for loop. But due to the list format of kriging output it might need some finesse. Especially for parameter starting values)
# This next section is just copy paste for making variograms for each year. Starting values and outlier removals that work are indicated
######################################################
d79 <- egg_sp[egg_sp@data$year == '1979', ]
# Remove outliers
my.vgm <- variogram(dep ~ 1, data = d79[d79$dep < quantile(d79$dep, .975), ]) 
plot(my.vgm)

# Normalise the variogram by dividing the data used by the variance
my.vgm$gamma <- my.vgm$gamma / var(d79[d79$dep<quantile(d79$dep,.975),]$dep)
plot(my.vgm)

# Set initial variogram parameters by eye 
model.i = vgm(psill= 1.2, model="Exp", range=75, nugget=0.1)
plot(my.vgm, model=model.i, plot.numbers=T, pch=16)

# Fit variogram from initial values (If no convergence, check tolerence of outlier threshhold. Above it is set to .99 in the quantile function)
model.f = fit.variogram(my.vgm, model=model.i,fit.ranges = T, fit.method=7) 
plot(my.vgm, model=model.f, plot.numbers=T, pch=16)

# Once the variogram is fit, rescale the semivariogram to the true values by multiplying by the variance of the values (undoing earlier division)
model.f$psill <- model.f$psill * var(d79$dep)
m79<-model.f

######################################################
d83 <- egg_sp[egg_sp@data$year == '1983', ]
# Remove outliers
my.vgm <- variogram(dep ~ 1, data = d83[d83$dep < quantile(d83$dep, .975), ]) 
plot(my.vgm)

# Normalise the variogram by dividing the data used by the variance
my.vgm$gamma <- my.vgm$gamma / var(d83[d83$dep<quantile(d83$dep,.975),]$dep)
plot(my.vgm)

# Set initial variogram parameters by eye 
model.i = vgm(psill= 2, model="Exp", range=75, nugget=0.1)
plot(my.vgm, model=model.i, plot.numbers=T, pch=16)

# Fit variogram from initial values (If no convergence, check tolerence of outlier threshhold. Above it is set to .99 in the quantile function)
model.f = fit.variogram(my.vgm, model=model.i,fit.ranges = T, fit.method=7) 
plot(my.vgm, model=model.f, plot.numbers=T, pch=16)

# Once the variogram is fit, rescale the semivariogram to the true values by multiplying by the variance of the values (undoing earlier division)
model.f$psill <- model.f$psill * var(d83$dep)
m83<-model.f

######################################################
d84 <- egg_sp[egg_sp@data$year == '1984', ]
# Remove outliers
my.vgm <- variogram(dep ~ 1, data = d84[d84$dep < quantile(d84$dep, .975), ]) 
plot(my.vgm)

# Normalise the variogram by dividing the data used by the variance
my.vgm$gamma <- my.vgm$gamma / var(d84[d84$dep<quantile(d84$dep,.975),]$dep)
plot(my.vgm)

# Set initial variogram parameters by eye 
model.i = vgm(psill= 1.2, model="Exp", range=75, nugget=0.1)
plot(my.vgm, model=model.i, plot.numbers=T, pch=16)

# Fit variogram from initial values (If no convergence, check tolerence of outlier threshhold. Above it is set to .99 in the quantile function)
model.f = fit.variogram(my.vgm, model=model.i,fit.ranges = T, fit.method=7) 
plot(my.vgm, model=model.f, plot.numbers=T, pch=16)

# Once the variogram is fit, rescale the semivariogram to the true values by multiplying by the variance of the values (undoing earlier division)
model.f$psill <- model.f$psill * var(d84$dep)
m84<-model.f


######################################################
d85 <- egg_sp[egg_sp@data$year == '1985', ]
# Remove outliers
my.vgm <- variogram(dep ~ 1, data = d85[d85$dep < quantile(d85$dep, .975), ]) 
plot(my.vgm)

# Normalise the variogram by dividing the data used by the variance
my.vgm$gamma <- my.vgm$gamma / var(d85[d85$dep<quantile(d85$dep,.975),]$dep)
plot(my.vgm)

# Set initial variogram parameters by eye 
model.i = vgm(psill= 1, model="Exp", range=75, nugget=0.4)
plot(my.vgm, model=model.i, plot.numbers=T, pch=16)

# Fit variogram from initial values (If no convergence, check tolerence of outlier threshhold. Above it is set to .99 in the quantile function)
model.f = fit.variogram(my.vgm, model=model.i,fit.ranges = T, fit.method=7) 
plot(my.vgm, model=model.f, plot.numbers=T, pch=16)

# Once the variogram is fit, rescale the semivariogram to the true values by multiplying by the variance of the values (undoing earlier division)
model.f$psill <- model.f$psill * var(d85$dep)
m85<-model.f

######################################################
d86 <- egg_sp[egg_sp@data$year == '1986', ]
# Remove outliers
my.vgm <- variogram(dep ~ 1, data = d86[d86$dep < quantile(d86$dep, .975), ]) 
plot(my.vgm)

# Normalise the variogram by dividing the data used by the variance
my.vgm$gamma <- my.vgm$gamma / var(d86[d86$dep<quantile(d86$dep,.975),]$dep)
plot(my.vgm)

# Set initial variogram parameters by eye 
model.i = vgm(psill= 1, model="Exp", range=75, nugget=0.4)
plot(my.vgm, model=model.i, plot.numbers=T, pch=16)

# Fit variogram from initial values (If no convergence, check tolerence of outlier threshhold. Above it is set to .99 in the quantile function)
model.f = fit.variogram(my.vgm, model=model.i,fit.ranges = T, fit.method=7) 
plot(my.vgm, model=model.f, plot.numbers=T, pch=16)

# Once the variogram is fit, rescale the semivariogram to the true values by multiplying by the variance of the values (undoing earlier division)
model.f$psill <- model.f$psill * var(d86$dep)
m86<-model.f

######################################################
d87 <- egg_sp[egg_sp@data$year == '1987', ]
# Remove outliers
my.vgm <- variogram(dep ~ 1, data = d87[d87$dep < quantile(d87$dep, .975), ]) 
plot(my.vgm)

# Normalise the variogram by dividing the data used by the variance
my.vgm$gamma <- my.vgm$gamma / var(d87[d87$dep<quantile(d87$dep,.975),]$dep)
plot(my.vgm)

# Set initial variogram parameters by eye 
model.i = vgm(psill= 1, model="Exp", range=75, nugget=0.1)
plot(my.vgm, model=model.i, plot.numbers=T, pch=16)

# Fit variogram from initial values (If no convergence, check tolerence of outlier threshhold. Above it is set to .99 in the quantile function)
model.f = fit.variogram(my.vgm, model=model.i,fit.ranges = T, fit.method=7) 
plot(my.vgm, model=model.f, plot.numbers=T, pch=16)

# Once the variogram is fit, rescale the semivariogram to the true values by multiplying by the variance of the values (undoing earlier division)
model.f$psill <- model.f$psill * var(d87$dep)
m87<-model.f

######################################################
d88 <- egg_sp[egg_sp@data$year == '1988', ]
# Remove outliers
my.vgm <- variogram(dep ~ 1, data = d88[d88$dep < quantile(d88$dep, .975), ]) 
plot(my.vgm)

# Normalise the variogram by dividing the data used by the variance
my.vgm$gamma <- my.vgm$gamma / var(d88[d88$dep<quantile(d88$dep,.975),]$dep)
plot(my.vgm)

# Set initial variogram parameters by eye 
model.i = vgm(psill= 1, model="Exp", range=75, nugget=0.1)
plot(my.vgm, model=model.i, plot.numbers=T, pch=16)

# Fit variogram from initial values (If no convergence, check tolerence of outlier threshhold. Above it is set to .99 in the quantile function)
model.f = fit.variogram(my.vgm, model=model.i,fit.ranges = T, fit.method=7) 
plot(my.vgm, model=model.f, plot.numbers=T, pch=16)

# Once the variogram is fit, rescale the semivariogram to the true values by multiplying by the variance of the values (undoing earlier division)
model.f$psill <- model.f$psill * var(d88$dep)
m88<-model.f

######################################################
d89 <- egg_sp[egg_sp@data$year == '1989', ]
# Remove outliers
my.vgm <- variogram(dep ~ 1, data = d89[d89$dep < quantile(d89$dep, .99), ]) 
plot(my.vgm)

# Normalise the variogram by dividing the data used by the variance
my.vgm$gamma <- my.vgm$gamma / var(d89[d89$dep<quantile(d89$dep,.99),]$dep)
plot(my.vgm)

# Set initial variogram parameters by eye 
model.i = vgm(psill= 1, model="Exp", range=30, nugget=0)
plot(my.vgm, model=model.i, plot.numbers=T, pch=16)

# Fit variogram from initial values (If no convergence, check tolerence of outlier threshhold. Above it is set to .99 in the quantile function)
model.f = fit.variogram(my.vgm, model=model.i,fit.ranges = T, fit.method=7) 
plot(my.vgm, model=model.f, plot.numbers=T, pch=16)

# Once the variogram is fit, rescale the semivariogram to the true values by multiplying by the variance of the values (undoing earlier division)
model.f$psill <- model.f$psill * var(d89$dep)
m89<-model.f

######################################################
d90 <- egg_sp[egg_sp@data$year == '1990', ]
# Remove outliers
my.vgm <- variogram(dep ~ 1, data = d90[d90$dep < quantile(d90$dep, .975), ]) 
plot(my.vgm)

# Normalise the variogram by dividing the data used by the variance
my.vgm$gamma <- my.vgm$gamma / var(d90[d90$dep<quantile(d90$dep,.975),]$dep)
plot(my.vgm)

# Set initial variogram parameters by eye 
model.i = vgm(psill= 1, model="Exp", range=50, nugget=0)
plot(my.vgm, model=model.i, plot.numbers=T, pch=16)

# Fit variogram from initial values (If no convergence, check tolerence of outlier threshhold. Above it is set to .99 in the quantile function)
model.f = fit.variogram(my.vgm, model=model.i,fit.ranges = T, fit.method=7) 
plot(my.vgm, model=model.f, plot.numbers=T, pch=16)

# Once the variogram is fit, rescale the semivariogram to the true values by multiplying by the variance of the values (undoing earlier division)
model.f$psill <- model.f$psill * var(d90$dep)
m90<-model.f

######################################################
######################################################
d91 <- egg_sp[egg_sp@data$year == '1991', ]
# Remove outliers
my.vgm <- variogram(dep ~ 1, data = d91[d91$dep < quantile(d91$dep, .99), ]) 
plot(my.vgm)

# Normalise the variogram by dividing the data used by the variance
my.vgm$gamma <- my.vgm$gamma / var(d91[d91$dep<quantile(d91$dep,.99),]$dep)
plot(my.vgm)

# Set initial variogram parameters by eye 
model.i = vgm(psill= 1.2, model="Exp", range=50, nugget=0)
plot(my.vgm, model=model.i, plot.numbers=T, pch=16)

# Fit variogram from initial values (If no convergence, check tolerence of outlier threshhold. Above it is set to .99 in the quantile function)
model.f = fit.variogram(my.vgm, model=model.i,fit.ranges = T, fit.method=7) 
plot(my.vgm, model=model.f, plot.numbers=T, pch=16)

# Once the variogram is fit, rescale the semivariogram to the true values by multiplying by the variance of the values (undoing earlier division)
model.f$psill <- model.f$psill * var(d91$dep)
m91<-model.f

######################################################
######################################################
d92 <- egg_sp[egg_sp@data$year == '1992', ]
# Remove outliers
my.vgm <- variogram(dep ~ 1, data = d92[d92$dep < quantile(d92$dep, .99), ]) 
plot(my.vgm)

# Normalise the variogram by dividing the data used by the variance
my.vgm$gamma <- my.vgm$gamma / var(d92[d92$dep<quantile(d92$dep,.99),]$dep)
plot(my.vgm)

# Set initial variogram parameters by eye 
model.i = vgm(psill= 1, model="Exp", range=50, nugget=0)
plot(my.vgm, model=model.i, plot.numbers=T, pch=16)

# Fit variogram from initial values (If no convergence, check tolerence of outlier threshhold. Above it is set to .99 in the quantile function)
model.f = fit.variogram(my.vgm, model=model.i,fit.ranges = T, fit.method=7) 
plot(my.vgm, model=model.f, plot.numbers=T, pch=16)

# Once the variogram is fit, rescale the semivariogram to the true values by multiplying by the variance of the values (undoing earlier division)
model.f$psill <- model.f$psill * var(d92$dep)
m92<-model.f

######################################################
######################################################
d93 <- egg_sp[egg_sp@data$year == '1993', ]
# Remove outliers
my.vgm <- variogram(dep ~ 1, data = d93[d93$dep < quantile(d93$dep, .99), ]) 
plot(my.vgm)

# Normalise the variogram by dividing the data used by the variance
my.vgm$gamma <- my.vgm$gamma / var(d93[d93$dep<quantile(d93$dep,.99),]$dep)
plot(my.vgm)

# Set initial variogram parameters by eye 
model.i = vgm(psill= 1.1, model="Exp", range=75, nugget=0)
plot(my.vgm, model=model.i, plot.numbers=T, pch=16)

# Fit variogram from initial values (If no convergence, check tolerence of outlier threshhold. Above it is set to .99 in the quantile function)
model.f = fit.variogram(my.vgm, model=model.i,fit.ranges = T, fit.method=7) 
plot(my.vgm, model=model.f, plot.numbers=T, pch=16)

# Once the variogram is fit, rescale the semivariogram to the true values by multiplying by the variance of the values (undoing earlier division)
model.f$psill <- model.f$psill * var(d93$dep)
m93<-model.f

######################################################
######################################################
d94 <- egg_sp[egg_sp@data$year == '1994', ]
# Remove outliers
my.vgm <- variogram(dep ~ 1, data = d94[d94$dep < quantile(d94$dep, .99), ]) 
plot(my.vgm)

# Normalise the variogram by dividing the data used by the variance
my.vgm$gamma <- my.vgm$gamma / var(d94[d94$dep<quantile(d94$dep,.99),]$dep)
plot(my.vgm)

# Set initial variogram parameters by eye 
model.i = vgm(psill= 1.2, model="Exp", range=50, nugget=0)
plot(my.vgm, model=model.i, plot.numbers=T, pch=16)

# Fit variogram from initial values (If no convergence, check tolerence of outlier threshhold. Above it is set to .99 in the quantile function)
model.f = fit.variogram(my.vgm, model=model.i,fit.ranges = T, fit.method=7) 
plot(my.vgm, model=model.f, plot.numbers=T, pch=16)

# Once the variogram is fit, rescale the semivariogram to the true values by multiplying by the variance of the values (undoing earlier division)
model.f$psill <- model.f$psill * var(d94$dep)
m94<-model.f

######################################################
######################################################
d96 <- egg_sp[egg_sp@data$year == '1996', ]
# Remove outliers
my.vgm <- variogram(dep ~ 1, data = d96[d96$dep < quantile(d96$dep, .99), ]) 
plot(my.vgm)

# Normalise the variogram by dividing the data used by the variance
my.vgm$gamma <- my.vgm$gamma / var(d96[d96$dep<quantile(d96$dep,.99),]$dep)
plot(my.vgm)

# Set initial variogram parameters by eye 
model.i = vgm(psill= 1.1, model="Exp", range=25, nugget=0)
plot(my.vgm, model=model.i, plot.numbers=T, pch=16)

# Fit variogram from initial values (If no convergence, check tolerence of outlier threshhold. Above it is set to .99 in the quantile function)
model.f = fit.variogram(my.vgm, model=model.i,fit.ranges = T, fit.method=7) 
plot(my.vgm, model=model.f, plot.numbers=T, pch=16)

# Once the variogram is fit, rescale the semivariogram to the true values by multiplying by the variance of the values (undoing earlier division)
model.f$psill <- model.f$psill * var(d96$dep)
m96<-model.f

######################################################
######################################################
d98 <- egg_sp[egg_sp@data$year == '1998', ]
# Remove outliers
my.vgm <- variogram(dep ~ 1, data = d98[d98$dep < quantile(d98$dep, .99), ]) 
plot(my.vgm)

# Normalise the variogram by dividing the data used by the variance
my.vgm$gamma <- my.vgm$gamma / var(d98[d98$dep<quantile(d98$dep,.99),]$dep)
plot(my.vgm)

# Set initial variogram parameters by eye 
model.i = vgm(psill= 1.2, model="Exp", range=50, nugget=0)
plot(my.vgm, model=model.i, plot.numbers=T, pch=16)

# Fit variogram from initial values (If no convergence, check tolerence of outlier threshhold. Above it is set to .99 in the quantile function)
model.f = fit.variogram(my.vgm, model=model.i,fit.ranges = T, fit.method=7) 
plot(my.vgm, model=model.f, plot.numbers=T, pch=16)

# Once the variogram is fit, rescale the semivariogram to the true values by multiplying by the variance of the values (undoing earlier division)
model.f$psill <- model.f$psill * var(d98$dep)
m98<-model.f

######################################################
######################################################
d99 <- egg_sp[egg_sp@data$year == '1999', ]
# Remove outliers
my.vgm <- variogram(dep ~ 1, data = d99[d99$dep < quantile(d99$dep, .99), ]) 
plot(my.vgm)

# Normalise the variogram by dividing the data used by the variance
my.vgm$gamma <- my.vgm$gamma / var(d99[d99$dep<quantile(d99$dep,.99),]$dep)
plot(my.vgm)

# Set initial variogram parameters by eye 
model.i = vgm(psill= 1.2, model="Exp", range=50, nugget=0)
plot(my.vgm, model=model.i, plot.numbers=T, pch=16)

# Fit variogram from initial values (If no convergence, check tolerence of outlier threshhold. Above it is set to .99 in the quantile function)
model.f = fit.variogram(my.vgm, model=model.i,fit.ranges = T, fit.method=7) 
plot(my.vgm, model=model.f, plot.numbers=T, pch=16)

# Once the variogram is fit, rescale the semivariogram to the true values by multiplying by the variance of the values (undoing earlier division)
model.f$psill <- model.f$psill * var(d99$dep)
m99<-model.f

######################################################
######################################################
d00 <- egg_sp[egg_sp@data$year == '2000', ]
# Remove outliers
my.vgm <- variogram(dep ~ 1, data = d00[d00$dep < quantile(d00$dep, .99), ]) 
plot(my.vgm)

# Normalise the variogram by dividing the data used by the variance
my.vgm$gamma <- my.vgm$gamma / var(d00[d00$dep<quantile(d00$dep,.99),]$dep)
plot(my.vgm)

# Set initial variogram parameters by eye 
model.i = vgm(psill= 1.2, model="Exp", range=50, nugget=0)
plot(my.vgm, model=model.i, plot.numbers=T, pch=16)

# Fit variogram from initial values (If no convergence, check tolerence of outlier threshhold. Above it is set to .99 in the quantile function)
model.f = fit.variogram(my.vgm, model=model.i,fit.ranges = T, fit.method=7) 
plot(my.vgm, model=model.f, plot.numbers=T, pch=16)

# Once the variogram is fit, rescale the semivariogram to the true values by multiplying by the variance of the values (undoing earlier division)
model.f$psill <- model.f$psill * var(d00$dep)
m00<-model.f

######################################################
######################################################
d01 <- egg_sp[egg_sp@data$year == '2001', ]
# Remove outliers
my.vgm <- variogram(dep ~ 1, data = d01[d01$dep < quantile(d01$dep, .975), ]) 
plot(my.vgm)

# Normalise the variogram by dividing the data used by the variance
my.vgm$gamma <- my.vgm$gamma / var(d01[d01$dep<quantile(d01$dep,.975),]$dep)
plot(my.vgm)

# Set initial variogram parameters by eye 
model.i = vgm(psill= 1, model="Exp", range=25, nugget=0)
plot(my.vgm, model=model.i, plot.numbers=T, pch=16)

# Fit variogram from initial values (If no convergence, check tolerence of outlier threshhold. Above it is set to .99 in the quantile function)
model.f = fit.variogram(my.vgm, model=model.i,fit.ranges = T, fit.method=7) 
plot(my.vgm, model=model.f, plot.numbers=T, pch=16)

# Once the variogram is fit, rescale the semivariogram to the true values by multiplying by the variance of the values (undoing earlier division)
model.f$psill <- model.f$psill * var(d01$dep)
m01<-model.f

######################################################
######################################################
d02 <- egg_sp[egg_sp@data$year == '2002', ]
# Remove outliers
my.vgm <- variogram(dep ~ 1, data = d02[d02$dep < quantile(d02$dep, .99), ]) 
plot(my.vgm)

# Normalise the variogram by dividing the data used by the variance
my.vgm$gamma <- my.vgm$gamma / var(d02[d02$dep<quantile(d02$dep,.99),]$dep)
plot(my.vgm)

# Set initial variogram parameters by eye 
model.i = vgm(psill= 1.4, model="Exp", range=100, nugget=0)
plot(my.vgm, model=model.i, plot.numbers=T, pch=16)

# Fit variogram from initial values (If no convergence, check tolerence of outlier threshhold. Above it is set to .99 in the quantile function)
model.f = fit.variogram(my.vgm, model=model.i,fit.ranges = T, fit.method=7) 
plot(my.vgm, model=model.f, plot.numbers=T, pch=16)

# Once the variogram is fit, rescale the semivariogram to the true values by multiplying by the variance of the values (undoing earlier division)
model.f$psill <- model.f$psill * var(d02$dep)
m02<-model.f

######################################################
######################################################
d03 <- egg_sp[egg_sp@data$year == '2003', ]
# Remove outliers
my.vgm <- variogram(dep ~ 1, data = d03[d03$dep < quantile(d03$dep, .99), ]) 
plot(my.vgm)

# Normalise the variogram by dividing the data used by the variance
my.vgm$gamma <- my.vgm$gamma / var(d03[d03$dep<quantile(d03$dep,.99),]$dep)
plot(my.vgm)

# Set initial variogram parameters by eye 
model.i = vgm(psill= 1.4, model="Exp", range=100, nugget=0)
plot(my.vgm, model=model.i, plot.numbers=T, pch=16)

# Fit variogram from initial values (If no convergence, check tolerence of outlier threshhold. Above it is set to .99 in the quantile function)
model.f = fit.variogram(my.vgm, model=model.i,fit.ranges = T, fit.method=7) 
plot(my.vgm, model=model.f, plot.numbers=T, pch=16)

# Once the variogram is fit, rescale the semivariogram to the true values by multiplying by the variance of the values (undoing earlier division)
model.f$psill <- model.f$psill * var(d03$dep)
m03<-model.f

######################################################
######################################################
d04 <- egg_sp[egg_sp@data$year == '2004', ]
# Remove outliers
my.vgm <- variogram(dep ~ 1, data = d04[d04$dep < quantile(d04$dep, .99), ]) 
plot(my.vgm)

# Normalise the variogram by dividing the data used by the variance
my.vgm$gamma <- my.vgm$gamma / var(d04[d04$dep<quantile(d04$dep,.99),]$dep)
plot(my.vgm)

# Set initial variogram parameters by eye 
model.i = vgm(psill= 1, model="Exp", range=50, nugget=0)
plot(my.vgm, model=model.i, plot.numbers=T, pch=16)

# Fit variogram from initial values (If no convergence, check tolerence of outlier threshhold. Above it is set to .99 in the quantile function)
model.f = fit.variogram(my.vgm, model=model.i,fit.ranges = T, fit.method=7) 
plot(my.vgm, model=model.f, plot.numbers=T, pch=16)

# Once the variogram is fit, rescale the semivariogram to the true values by multiplying by the variance of the values (undoing earlier division)
model.f$psill <- model.f$psill * var(d04$dep)
m04<-model.f

######################################################
######################################################
d05 <- egg_sp[egg_sp@data$year == '2005', ]
# Remove outliers
my.vgm <- variogram(dep ~ 1, data = d05[d05$dep < quantile(d05$dep, .99), ]) 
plot(my.vgm)

# Normalise the variogram by dividing the data used by the variance
my.vgm$gamma <- my.vgm$gamma / var(d05[d05$dep<quantile(d05$dep,.99),]$dep)
plot(my.vgm)

# Set initial variogram parameters by eye 
model.i = vgm(psill= 1.4, model="Exp", range=100, nugget=0)
plot(my.vgm, model=model.i, plot.numbers=T, pch=16)

# Fit variogram from initial values (If no convergence, check tolerence of outlier threshhold. Above it is set to .99 in the quantile function)
model.f = fit.variogram(my.vgm, model=model.i,fit.ranges = T, fit.method=7) 
plot(my.vgm, model=model.f, plot.numbers=T, pch=16)

# Once the variogram is fit, rescale the semivariogram to the true values by multiplying by the variance of the values (undoing earlier division)
model.f$psill <- model.f$psill * var(d05$dep)
m05<-model.f

######################################################
######################################################
d06 <- egg_sp[egg_sp@data$year == '2006', ]
# Remove outliers
my.vgm <- variogram(dep ~ 1, data = d06[d06$dep < quantile(d06$dep, .99), ]) 
plot(my.vgm)

# Normalise the variogram by dividing the data used by the variance
my.vgm$gamma <- my.vgm$gamma / var(d06[d06$dep<quantile(d06$dep,.99),]$dep)
plot(my.vgm)

# Set initial variogram parameters by eye 
model.i = vgm(psill= 1.4, model="Exp", range=100, nugget=0)
plot(my.vgm, model=model.i, plot.numbers=T, pch=16)

# Fit variogram from initial values (If no convergence, check tolerence of outlier threshhold. Above it is set to .99 in the quantile function)
model.f = fit.variogram(my.vgm, model=model.i,fit.ranges = T, fit.method=7) 
plot(my.vgm, model=model.f, plot.numbers=T, pch=16)

# Once the variogram is fit, rescale the semivariogram to the true values by multiplying by the variance of the values (undoing earlier division)
model.f$psill <- model.f$psill * var(d06$dep)
m06<-model.f

######################################################
######################################################
d07 <- egg_sp[egg_sp@data$year == '2007', ]
# Remove outliers
my.vgm <- variogram(dep ~ 1, data = d07[d07$dep < quantile(d07$dep, .99), ]) 
plot(my.vgm)

# Normalise the variogram by dividing the data used by the variance
my.vgm$gamma <- my.vgm$gamma / var(d07[d07$dep<quantile(d07$dep,.99),]$dep)
plot(my.vgm)

# Set initial variogram parameters by eye 
model.i = vgm(psill= 1.4, model="Exp", range=100, nugget=0)
plot(my.vgm, model=model.i, plot.numbers=T, pch=16)

# Fit variogram from initial values (If no convergence, check tolerence of outlier threshhold. Above it is set to .99 in the quantile function)
model.f = fit.variogram(my.vgm, model=model.i,fit.ranges = T, fit.method=7) 
plot(my.vgm, model=model.f, plot.numbers=T, pch=16)

# Once the variogram is fit, rescale the semivariogram to the true values by multiplying by the variance of the values (undoing earlier division)
model.f$psill <- model.f$psill * var(d07$dep)
m07<-model.f

######################################################
######################################################
d08 <- egg_sp[egg_sp@data$year == '2008', ]
# Remove outliers
my.vgm <- variogram(dep ~ 1, data = d08[d08$dep < quantile(d08$dep, .975), ]) 
plot(my.vgm)

# Normalise the variogram by dividing the data used by the variance
my.vgm$gamma <- my.vgm$gamma / var(d08[d08$dep<quantile(d08$dep,.975),]$dep)
plot(my.vgm)

# Set initial variogram parameters by eye 
model.i = vgm(psill= 1, model="Exp", range=25, nugget=0)
plot(my.vgm, model=model.i, plot.numbers=T, pch=16)

# Fit variogram from initial values (If no convergence, check tolerence of outlier threshhold. Above it is set to .99 in the quantile function)
model.f = fit.variogram(my.vgm, model=model.i,fit.ranges = T, fit.method=7) 
plot(my.vgm, model=model.f, plot.numbers=T, pch=16)

# Once the variogram is fit, rescale the semivariogram to the true values by multiplying by the variance of the values (undoing earlier division)
model.f$psill <- model.f$psill * var(d08$dep)
m08<-model.f

######################################################
######################################################
d09 <- egg_sp[egg_sp@data$year == '2009', ]
# Remove outliers
my.vgm <- variogram(dep ~ 1, data = d09[d09$dep < quantile(d09$dep, .975), ]) 
plot(my.vgm)

# Normalise the variogram by dividing the data used by the variance
my.vgm$gamma <- my.vgm$gamma / var(d09[d09$dep<quantile(d09$dep,.975),]$dep)
plot(my.vgm)

# Set initial variogram parameters by eye 
model.i = vgm(psill= 1, model="Exp", range=50, nugget=0)
plot(my.vgm, model=model.i, plot.numbers=T, pch=16)

# Fit variogram from initial values (If no convergence, check tolerence of outlier threshhold. Above it is set to .99 in the quantile function)
model.f = fit.variogram(my.vgm, model=model.i,fit.ranges = T, fit.method=7) 
plot(my.vgm, model=model.f, plot.numbers=T, pch=16)

# Once the variogram is fit, rescale the semivariogram to the true values by multiplying by the variance of the values (undoing earlier division)
model.f$psill <- model.f$psill * var(d09$dep)
m09<-model.f

######################################################
######################################################
d10 <- egg_sp[egg_sp@data$year == '2010', ]
# Remove outliers
my.vgm <- variogram(dep ~ 1, data = d10[d10$dep < quantile(d10$dep, .5), ]) 
plot(my.vgm)

# Normalise the variogram by dividing the data used by the variance
my.vgm$gamma <- my.vgm$gamma / var(d10[d10$dep<quantile(d10$dep,.5),]$dep)
plot(my.vgm)

# Set initial variogram parameters by eye 
model.i = vgm(psill= 1.75, model="Exp", range=150, nugget=0)
plot(my.vgm, model=model.i, plot.numbers=T, pch=16)

# Fit variogram from initial values (If no convergence, check tolerence of outlier threshhold. Above it is set to .99 in the quantile function)
model.f = fit.variogram(my.vgm, model=model.i,fit.ranges = T, fit.method=7) 
plot(my.vgm, model=model.f, plot.numbers=T, pch=16)

# Once the variogram is fit, rescale the semivariogram to the true values by multiplying by the variance of the values (undoing earlier division)
model.f$psill <- model.f$psill * var(d10$dep)
m10<-model.f

######################################################
######################################################
d11 <- egg_sp[egg_sp@data$year == '2011', ]
# Remove outliers
my.vgm <- variogram(dep ~ 1, data = d11[d11$dep < quantile(d11$dep, .99), ]) 
plot(my.vgm)

# Normalise the variogram by dividing the data used by the variance
my.vgm$gamma <- my.vgm$gamma / var(d11[d11$dep<quantile(d11$dep,.99),]$dep)
plot(my.vgm)

# Set initial variogram parameters by eye 
model.i = vgm(psill= 1, model="Exp", range=100, nugget=0)
plot(my.vgm, model=model.i, plot.numbers=T, pch=16)

# Fit variogram from initial values (If no convergence, check tolerence of outlier threshhold. Above it is set to .99 in the quantile function)
model.f = fit.variogram(my.vgm, model=model.i,fit.ranges = T, fit.method=7) 
plot(my.vgm, model=model.f, plot.numbers=T, pch=16)

# Once the variogram is fit, rescale the semivariogram to the true values by multiplying by the variance of the values (undoing earlier division)
model.f$psill <- model.f$psill * var(d11$dep)
m11<-model.f

######################################################
######################################################
d12 <- egg_sp[egg_sp@data$year == '2012', ]
# Remove outliers
my.vgm <- variogram(dep ~ 1, data = d12[d12$dep < quantile(d12$dep, .9), ]) 
plot(my.vgm)

# Normalise the variogram by dividing the data used by the variance
my.vgm$gamma <- my.vgm$gamma / var(d12[d12$dep<quantile(d12$dep,.9),]$dep)
plot(my.vgm)

# Set initial variogram parameters by eye 
model.i = vgm(psill= 0.9, model="Exp", range=75, nugget=0.2)
plot(my.vgm, model=model.i, plot.numbers=T, pch=16)

# Fit variogram from initial values (If no convergence, check tolerence of outlier threshhold. Above it is set to .99 in the quantile function)
model.f = fit.variogram(my.vgm, model=model.i,fit.ranges = T, fit.method=7) 
plot(my.vgm, model=model.f, plot.numbers=T, pch=16)

# Once the variogram is fit, rescale the semivariogram to the true values by multiplying by the variance of the values (undoing earlier division)
model.f$psill <- model.f$psill * var(d12$dep)
m12<-model.f

######################################################
######################################################
d13 <- egg_sp[egg_sp@data$year == '2013', ]
# Remove outliers
my.vgm <- variogram(dep ~ 1, data = d13[d13$dep < quantile(d13$dep, .9), ]) 
plot(my.vgm)

# Normalise the variogram by dividing the data used by the variance
my.vgm$gamma <- my.vgm$gamma / var(d13[d13$dep<quantile(d13$dep,.9),]$dep)
plot(my.vgm)

# Set initial variogram parameters by eye 
model.i = vgm(psill= 1, model="Exp", range=75, nugget=0)
plot(my.vgm, model=model.i, plot.numbers=T, pch=16)

# Fit variogram from initial values (If no convergence, check tolerence of outlier threshhold. Above it is set to .99 in the quantile function)
model.f = fit.variogram(my.vgm, model=model.i,fit.ranges = T, fit.method=7) 
plot(my.vgm, model=model.f, plot.numbers=T, pch=16)

# Once the variogram is fit, rescale the semivariogram to the true values by multiplying by the variance of the values (undoing earlier division)
model.f$psill <- model.f$psill * var(d13$dep)
m13<-model.f

######################################################
######################################################
d14 <- egg_sp[egg_sp@data$year == '2014', ]
# Remove outliers
my.vgm <- variogram(dep ~ 1, data = d14[d14$dep < quantile(d14$dep, .975), ]) 
plot(my.vgm)

# Normalise the variogram by dividing the data used by the variance
my.vgm$gamma <- my.vgm$gamma / var(d14[d14$dep<quantile(d14$dep,.975),]$dep)
plot(my.vgm)

# Set initial variogram parameters by eye 
model.i = vgm(psill= 1, model="Exp", range=100, nugget=0)
plot(my.vgm, model=model.i, plot.numbers=T, pch=16)

# Fit variogram from initial values (If no convergence, check tolerence of outlier threshhold. Above it is set to .99 in the quantile function)
model.f = fit.variogram(my.vgm, model=model.i,fit.ranges = T, fit.method=7) 
plot(my.vgm, model=model.f, plot.numbers=T, pch=16)

# Once the variogram is fit, rescale the semivariogram to the true values by multiplying by the variance of the values (undoing earlier division)
model.f$psill <- model.f$psill * var(d14$dep)
m14<-model.f

######################################################
######################################################
d15 <- egg_sp[egg_sp@data$year == '2015', ]
# Remove outliers
my.vgm <- variogram(dep ~ 1, data = d15[d15$dep < quantile(d15$dep, .99), ]) 
plot(my.vgm)

# Normalise the variogram by dividing the data used by the variance
my.vgm$gamma <- my.vgm$gamma / var(d15[d15$dep<quantile(d15$dep,.99),]$dep)
plot(my.vgm)

# Set initial variogram parameters by eye 
model.i = vgm(psill= 1, model="Exp", range=100, nugget=0)
plot(my.vgm, model=model.i, plot.numbers=T, pch=16)

# Fit variogram from initial values (If no convergence, check tolerence of outlier threshhold. Above it is set to .99 in the quantile function)
model.f = fit.variogram(my.vgm, model=model.i,fit.ranges = T, fit.method=7) 
plot(my.vgm, model=model.f, plot.numbers=T, pch=16)

# Once the variogram is fit, rescale the semivariogram to the true values by multiplying by the variance of the values (undoing earlier division)
model.f$psill <- model.f$psill * var(d15$dep)
m15<-model.f

######################################################
######################################################
d16 <- egg_sp[egg_sp@data$year == '2016', ]
# Remove outliers
my.vgm <- variogram(dep ~ 1, data = d16[d16$dep < quantile(d16$dep, .9), ]) 
plot(my.vgm)

# Normalise the variogram by dividing the data used by the variance
my.vgm$gamma <- my.vgm$gamma / var(d16[d16$dep<quantile(d16$dep,.9),]$dep)
plot(my.vgm)

# Set initial variogram parameters by eye 
model.i = vgm(psill= 1.2, model="Exp", range=150, nugget=0.2)
plot(my.vgm, model=model.i, plot.numbers=T, pch=16)

# Fit variogram from initial values (If no convergence, check tolerence of outlier threshhold. Above it is set to .99 in the quantile function)
model.f = fit.variogram(my.vgm, model=model.i,fit.ranges = T, fit.method=7) 
plot(my.vgm, model=model.f, plot.numbers=T, pch=16)

# Once the variogram is fit, rescale the semivariogram to the true values by multiplying by the variance of the values (undoing earlier division)
model.f$psill <- model.f$psill * var(d16$dep)
m16<-model.f

######################################################
######################################################
d17 <- egg_sp[egg_sp@data$year == '2017', ]
# Remove outliers
my.vgm <- variogram(dep ~ 1, data = d17[d17$dep < quantile(d17$dep, .99), ]) 
plot(my.vgm)

# Normalise the variogram by dividing the data used by the variance
my.vgm$gamma <- my.vgm$gamma / var(d17[d17$dep<quantile(d17$dep,.99),]$dep)
plot(my.vgm)

# Set initial variogram parameters by eye 
model.i = vgm(psill= 1, model="Exp", range=100, nugget=0)
plot(my.vgm, model=model.i, plot.numbers=T, pch=16)

# Fit variogram from initial values (If no convergence, check tolerence of outlier threshhold. Above it is set to .99 in the quantile function)
model.f = fit.variogram(my.vgm, model=model.i,fit.ranges = T, fit.method=7) 
plot(my.vgm, model=model.f, plot.numbers=T, pch=16)

# Once the variogram is fit, rescale the semivariogram to the true values by multiplying by the variance of the values (undoing earlier division)
model.f$psill <- model.f$psill * var(d17$dep)
m17<-model.f

######################################################
######################################################
d18 <- egg_sp[egg_sp@data$year == '2018', ]
# Remove outliers
my.vgm <- variogram(dep ~ 1, data = d18[d18$dep < quantile(d18$dep, .99), ]) 
plot(my.vgm)

# Normalise the variogram by dividing the data used by the variance
my.vgm$gamma <- my.vgm$gamma / var(d18[d18$dep<quantile(d18$dep,.99),]$dep)
plot(my.vgm)

# Set initial variogram parameters by eye 
model.i = vgm(psill= 1, model="Exp", range=100, nugget=0)
plot(my.vgm, model=model.i, plot.numbers=T, pch=16)

# Fit variogram from initial values (If no convergence, check tolerence of outlier threshhold. Above it is set to .99 in the quantile function)
model.f = fit.variogram(my.vgm, model=model.i,fit.ranges = T, fit.method=7) 
plot(my.vgm, model=model.f, plot.numbers=T, pch=16)

# Once the variogram is fit, rescale the semivariogram to the true values by multiplying by the variance of the values (undoing earlier division)
model.f$psill <- model.f$psill * var(d18$dep)
m18<-model.f

######################################################
######################################################
d19 <- egg_sp[egg_sp@data$year == '2019', ]
# Remove outliers
my.vgm <- variogram(dep ~ 1, data = d19[d19$dep < quantile(d19$dep, .99), ]) 
plot(my.vgm)

# Normalise the variogram by dividing the data used by the variance
my.vgm$gamma <- my.vgm$gamma / var(d19[d19$dep<quantile(d19$dep,.99),]$dep)
plot(my.vgm)

# Set initial variogram parameters by eye 
model.i = vgm(psill= 1, model="Exp", range=100, nugget=0)
plot(my.vgm, model=model.i, plot.numbers=T, pch=16)

# Fit variogram from initial values (If no convergence, check tolerence of outlier threshhold. Above it is set to .99 in the quantile function)
model.f = fit.variogram(my.vgm, model=model.i,fit.ranges = T, fit.method=7) 
plot(my.vgm, model=model.f, plot.numbers=T, pch=16)

# Once the variogram is fit, rescale the semivariogram to the true values by multiplying by the variance of the values (undoing earlier division)
model.f$psill <- model.f$psill * var(d19$dep)
m19<-model.f

######################################################
d21 <- egg_sp[egg_sp@data$year == '2021', ]
# Remove outliers
my.vgm <- variogram(dep ~ 1, data = d21[d21$dep < quantile(d21$dep, .99), ]) 
plot(my.vgm)

# Normalise the variogram by dividing the data used by the variance
my.vgm$gamma <- my.vgm$gamma / var(d21[d21$dep<quantile(d21$dep,.99),]$dep)
plot(my.vgm)

# Set initial variogram parameters by eye 
model.i = vgm(psill= 1, model="Exp", range=25, nugget=0)
plot(my.vgm, model=model.i, plot.numbers=T, pch=16)

# Fit variogram from initial values (If no convergence, check tolerence of outlier threshhold. Above it is set to .99 in the quantile function)
model.f = fit.variogram(my.vgm, model=model.i, fit.ranges = T, fit.method=7) 
plot(my.vgm, model=model.f, plot.numbers=T, pch=16)

# Once the variogram is fit, rescale the semivariogram to the true values by multiplying by the variance of the values (undoing earlier division)
model.f$psill <- model.f$psill * var(d21$dep)
m21<-model.f

# store all the fitted variograms in a list
varios <- list(m79,m83,m84,m85,m86,m87,m88,m89,m90,m91,m92,m93,m94,m96,m98,m99,m00,m01,m02,m03,m04,m05,m06,m07,m08,m09,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m21)
data_sp <- list(d79,d83,d84,d85,d86,d87,d88,d89,d90,d91,d92,d93,d94,d96,d98,d99,d00,d01,d02,d03,d04,d05,d06,d07,d08,d09,d10,d11,d12,d13,d14,d15,d16,d17,d18,d19,d21) 

names(data_sp)<- c(unique(egg_sp$year))
names(varios)<- c(unique(egg_sp$year))
save(varios, file = "./rdata/new_varios_2022.Rdata")
save(data_sp, file = "./rdata/data_sp_split_year.Rdata")

##############################################################

#####################  Testing whether automap::autokrige produces similar results  ####################
# just checking 2021
eggs21 <- data.frame(eggs21)
grid21 <- data.frame(grid21)
coordinates(eggs21) =~ lon+lat
coordinates(grid21) =~ lon+lat

test <- autoKrige(dep~1, eggs21, grid21, model = "Exp")   # similar but slightly different

# Extracting parts from the autoKrige object
prediction_spdf = test$krige_output
sample_variogram = test$exp_var
variogram_model = test$var_model
ssplot(prediction_spdf)

moy.var.Krigeage(
  xCoords = coordinates(eggs21), 
  xValues = eggs21$dep, 
  kCoords = coordinates(grid21), 
  vario = variogram_model)

#####################  Testing whether automap::autokrige produces similar results  ####################

# Empty data frames for later use
pred=data.frame()  # dataframe with predictions krigging
moy_var=data.frame() # dataframe wih mean and variance of krigging predictions of whole area
egg_survey$pred=NA  


#####################  Krigging and estimating means and variances for annual egg index  - One year at a time ####################

# For one year at a time
# Interpolate DEP across DEP via ordinary krigging (using fitted semivariogram)
krige = krige(dep ~ 1, egg_sp[egg_sp@data$year == '2021', ],
              grid_sp,
              m21,
              nmin = 16)

mean(krige$var1.pred) # running mean and var on the krigged results gives you something similar to the moy.var.Krige function but a geospatial specialist wrote the latter sooooo
spplot(krige["var1.pred"])
krigdf<- as.data.frame(krige) 

moy_var = moy.var.Krigeage(
  xCoords = coordinates(d21), 
  xValues = d21$dep, 
  kCoords = coordinates(grid_sp), 
  vario = m21) %>% as.data.frame() %>% mutate(year = 2021)
# moy_var1=as.data.frame(moy_var1)
# moy_var1$YEAR
# moy_var1$z.noneg=mean(pred1[!pred1$pred<0,'pred'])
# moy_var=rbind(moy_var,moy_var1)
# pred$neg=ifelse(pred$pred<0,'-','') # Negative predictions indicated on map
moy_var$sd = sqrt(moy_var$VarKrigeage)
moy_var$CV = (moy_var$sd/moy_var$Z)*100



# Krig plot
krig_plot1 <- ggplot(data = krigdf, aes(x = lon, y = lat))  + 
  stat_summary_2d(data=krigdf,aes(x = lon, y = lat,z = var1.pred),fun = mean, binwidth = c(0.1, 0.1)) +
  geom_polygon(data = canada, aes(x=lon, y = lat, group = group), fill = "ghostwhite",color="black")  +
  coord_map(xlim = c(-66.5, -60),  ylim = c(45.5, 49.5),projection = "lambert", parameters = c(46.5 ,48)) +
  theme_bw() + 
  # geom_point(data = egg, aes(lon,lat, size = dep), pch = 21) +  
  scale_fill_viridis_c(option="viridis", begin = 0.1) + 
  labs(x = "Longitude",y = "Latitude", title = "Interpolated Daily Egg Production")

krig_plot2 <- ggplot(data = krigdf2, aes(x = lon, y = lat))  + 
  stat_summary_2d(data=krigdf2,aes(x = lon, y = lat,z = var1.pred),fun = mean, binwidth = c(0.1, 0.1)) +
  geom_polygon(data = canada, aes(x=lon, y = lat, group = group), fill = "ghostwhite",color="black")  +
  coord_map(xlim = c(-66.5, -60),  ylim = c(45.5, 49.5),projection = "lambert", parameters = c(46.5 ,48)) +
  theme_bw() + 
  # geom_point(data = egg, aes(lon,lat, size = dep), pch = 21) +  
  scale_fill_viridis_c(option="viridis", begin = 0.1) + 
  labs(x = "Longitude",y = "Latitude", title = "Interpolated Daily Egg Production")


#####################  Krigging and estimating means and variances for annual egg index - All years ####################


year=c(unique(egg_sp$year))
krig_results<-data.frame()

for(i in unique(egg_sp$year)){
  # Create subset and change year 
  df_year <- egg_sp[egg_sp$year==i,]
  tmp <- moy.var.Krigeage(
    xCoords = coordinates(df_year), 
    xValues = df_year$dep, 
    kCoords = coordinates(grid_sp), 
    vario = varios[[as.character(i)]])
  tmp <- as.data.frame(tmp)
  krig_results <- rbind(krig_results, tmp)
}
krig_results$year <- c(unique(egg_sp$year))
krig_results$sd = sqrt(krig_results$VarKrigeage)
krig_results$cv = (krig_results$sd/krig_results$Z)*100
krig_results$CI_upper = krig_results$Z + (krig_results$sd * 1.96)
krig_results$CI_lower = krig_results$Z - (krig_results$sd * 1.96)


save(krig_results, file ="./rdata/new_krigged_dep_2022.Rdata")
# cor and RSS - code not adapted to current file format
# sampling stations are not located on grid coordinates:  find nearest predicted value
df<-egg %>% dplyr::filter(trajet=="T1") %>% dplyr::select(year,dep)
df<-left_join(df,krig_results)
df %>% mutate(res = dep-Z) %>% group_by(year) %>% mutate(cor = cor(dep,Z), rss=sum(res^2))

egg_survey_df$res=egg_survey_df$DEP-egg_survey_df$pred
FIT=ddply(egg_survey_df,c('YEAR'),summarise,COR=cor(DEP,pred),RSS=sum(res^2))

# the results here are very similar to those published. However those before 1987 don't really match at all. Could be a variety of reasons. The corrections applied by FG, the variogram model (here I only use Exp but Francois used a mix of Sph and Exp, the variogram parameters in general, that other programs were used for stats, that stations in the past before multinet tech applied oversampling correction factors to stations... etc. etc.)

krig_results %>% ggplot(aes(year,Z)) + geom_point() +
  geom_ribbon(aes(ymin=CI_lower,ymax=CI_upper),alpha=0.5,fill="red") + 
  geom_line(aes(year,Z))

egg %>% dplyr::filter(trajet=="T1") %>% 
  ggplot(aes(year,dep))+stat_summary(fun.data = "mean_cl_boot")
egg



# calculate simple random mean and geometric mean
source("./functions/gm_mean.R")
egg_summary <- egg %>% dplyr::filter(trajet=="T1") %>% 
  group_by(year) %>%
  dplyr::summarise(dep_sr_mean = mean(dep, na.rm=T),   # simple random mean
                   var_dep = var(dep, na.rm=T),
                   dep_geo_mean = gm_mean(dep, na.rm = T),   # geometric mean
                   dep_w_mean = weighted.mean(dep, w = stratum_area, na.rm = T) # weighted mean
                   )

# calculate stratified random mean as per P. Ouellet

surface.strate = c(29.61e+9,21.91e+9,17.93e+9)
sum(surface.strate)

# DEP1=ddply(eggs,c('year'),summarise,DEP1=mean(DEP.weigthed)*surf) # script from EVB and Baye from 2017
library(plyr)
egg2<- egg %>% dplyr::filter(trajet =="T1") %>% 
  mutate(dep_weigthed = dep/stratum_area)
dep=ddply(egg2,c('year'),summarise,dep_sw=mean(dep_weigthed)*sum(surface.strate)) # works out to same as below

dat <- left_join(egg_summary, dep)



# note that I included years that were previously ommitted (eg 1982, 1999, 2001). My intent was to do an extended model run as a sensitivity test. And honnestly, as CCAM assumes there is obs error why not?



# compare with FG resdoc 2014
FG <- read_csv("./data/ichthyoplankton/metadata_and_results/fg_survey_results_2014_resdoc.csv")
FG %<>% dplyr::select(1,2,14)
egg_survey_stats <- read_csv("data/ichthyoplankton/metadata_and_results/egg_survey_stats.csv")
ess<-egg_survey_stats %>% dplyr::select(1,13,17,21)
fg<-left_join(FG,ess)
fg %<>% transmute(year, fg_mean_nm2 = mean_egg_density_nm2 , fg_sr_mean_dep =`simple random mean`, fg_weighted_strat_mean_dep = `global stratified mean`, fg_krigged_dep = krigged_mean_dep)

dat<-left_join(dat,fg)
dat %<>% pivot_longer(cols = 2:10, names_to = "variable", values_to = "value")

dat %>% dplyr::filter(variable !="var_dep") %>% 
  ggplot(aes(year,value,colour=variable))+geom_point() 
# compare annual simple random mean of dep
dat %>% dplyr::filter(variable %in% c("dep_sr_mean", "fg_sr_mean_dep")) %>% 
  ggplot(aes(year,value,colour=variable))+geom_point() # similar except for 1979 and 1986
# compare weighted means
dat %>% dplyr::filter(variable %in% c("dep_sw", "dep_w_mean", "fg_weighted_strat_mean_dep")) %>% 
  ggplot(aes(year,value,colour=variable))+geom_point() # trends similar but prior to ...95 values differ greatly. post 95 fg and weighted mean values more alike than stratified weighted...
dat %>% dplyr::filter(variable %in% c("dep_sw", "dep_w_mean", "fg_weighted_strat_mean_dep")) %>% 
  ggplot(aes(year,value,colour=variable))+geom_smooth() # 

# compare krigged dep
egg_survey_stats <- read_csv("data/ichthyoplankton/metadata_and_results/egg_survey_stats.csv")
ess<-egg_survey_stats %>% dplyr::select(year,krigged_mean_dep)
ess<-left_join(ess,krig_results)
ess %>% 
  ggplot(aes(year,Z))+geom_point()+geom_point(aes(year,krigged_mean_dep),colour="red") # these at least are similar
ess %>% ggplot(aes(Z,krigged_mean_dep,label=year))+geom_text()
