####### ==================================================================================================
##
## 2018 Mackerel stock assessment: Egg survey
##
## By Andrew Smith August 2018
##
## R project = 2018_Mackerel_Evaluation + associated files in 2018_Mackerel_Evaluation folder
##
## 
####### ==================================================================================================
# ==================================================================================================
# 
# Load packages
packages = c('tidyverse','magrittr','forcats', 'lubridate', 'raster','gstat','ggthemes','maptools','mapdata','gridExtra','sp','gstat','rgeos','rgdal','geoR')
invisible(lapply(packages, function(x) {if (!require(x, character.only = T)) {install.packages(x);require(x)}}))

# Load Data
egg_survey <- read.csv("./Data/Egg_survey/Mackerel_egg_survey_1982_2018.csv") # 2017 & 2018 preliminary
egg_survey <- egg_survey[,1:25]
load("./Data/Rdata/maq_carbio_clean.Rdata")

path2 <- "./Data/Egg_survey/"
path3 <- "./Data/Egg_survey/stations.csv"
grille <- read.table(paste0(path2,"grille_pred.csv"), header = TRUE, sep=";")
stations <- read.csv(path3)
stations <- stations[, 1:4]
stations$STATION <- as.factor(stations$STATION)

# ================================================================================================== 
# 
# Egg and Larvae counts from June PMZA/Mackerel Southern GUlf of Saint Lawrence (sGSL) Survey

# Transformations to convert data from some columns from wide to tall format 
egg_tall <- egg_survey %>% tidyr::gather(STAGE, COUNT, 12:17) %>%    #gather all the species 
 arrange(YEAR, TRIP, STATION) %>% 
  as_tibble() %>% 
  dplyr::select(-VOLUME)
egg_tall$COUNT[is.na(egg_tall$COUNT)] <- 0
 
# Survey Timing

# Theme for plotting data
MY_theme <- theme(axis.text.x = element_text(angle = 90, hjust = 1))

# In previous years, stations were surveyed twice (second mission in summer or as a second pass later in june). Choose Trip 1 
egg_count <- egg_tall %>%
  filter(TRIP == "1", YEAR != "2009b",!is.na(YEAR)) %>%
  group_by(YEAR, STAGE) %>%
  dplyr::summarise(EGGCOUNT = sum(COUNT)) 

egg_count <- egg_count[7:216,] # remove weird NAs

# Rename variables for plotting
level_key <- list(STADE.1 = "1",
                  STADE.2 = "2",
                  STADE.3 = "3",
                  STADE.4 = "4",
                  STADE.5 = "5",
                  Scomber.scombrus.larvae = "Larves / Larvae")
egg_count$STAGE <- recode(egg_count$STAGE, !!!level_key)

# Egg stage count by year (Note that each facet has free scaling y axis)
eggcount_fig <- egg_count %>%
  ggplot() +
  geom_col(aes(x = STAGE, y = EGGCOUNT, fill = STAGE)) +
  facet_wrap(vars(YEAR), scales = "free") +
  theme_bw() +
  labs(y = "Nombre / Number", tag = "A") +
  theme(axis.text.x = element_blank(), legend.title = element_blank())

# Egg stage count by year time series
eggcount_fig_2 <- egg_count %>%
  ggplot() +
  geom_col(aes(x = YEAR, y = EGGCOUNT, fill = STAGE)) +
  facet_wrap(vars(STAGE)) + theme_bw() +
  scale_x_discrete(breaks = seq(1980, 2018, 2)) +
  MY_theme +
  labs(y = "Nombre / Number", x = "Année / Year", tag = "B") +
  theme(
    axis.text.x = element_blank(),
    legend.title = element_blank(),
    legend.position = "none")

# ================================================================================================== 
# 
# If station coordinates are missing ...
  stations$Long <- as.character(stations$Long)
  stations$Lat <- as.character(stations$Lat)
  
  egg_survey <- full_join(stations, egg_survey, by = "STATION") %>% arrange(YEAR, STATION)

# ================================================================================================== 
# 
# Select current year (Note that in the future try to use the split function to do all years at once)

  egg_survey <- egg_survey %>% filter(YEAR == "2015")   
  maq_carbio <- maq_carbio_clean %>% filter(YEAR == "2015") 
  
# ==================================================================================================
# 
# The parameters of the SSB egg index equation: B = (P*A*W)/(S*F*R*C)
#     P = Spawning biomass index (t)
#     A = Area (m^2) of the zone sampled (6.945 * 10^10 m^2) *To recalculate with new grid
#     W = Mean weight (g) of a fish (somatic)
#     S = Proportion of eggs spawned at ith median date of a survey
#     F = Fecundity of females (Pelletier 1986)
#     R = Proportion of females in the sample 
#     C = 10^6 =  factor for converting grams into tonnes
  
# ==================================================================================================
#
# Check for NAs from Bionet data and replace with debimetre data. 
  
  # Future step with new mission variable with be colmatage (bongo net clogging yes/no) so as to help choosing Babord or Tribord sample or omitting the data
  # Use data from surrounding stations for depth and temperature if no data available (not for volume)
  
  anyNA(egg_survey$Depth_Data) # should be false
  anyNA(egg_survey$Volume_Bionet) # should be false
  anyNA(egg_survey$Temperature_10m) # should be false

  # Notes for specific years (Run these if redoing analyses)
    # 2014
      # e.g. station 4.7 in 2014 missing volume_bionet and temperature_10m. Replacement volume value = Volume_data = 314, and replacement temp = (9.54 + 9.6 + 9.1 + 9.61 + 8.33 + 10.04 + 10.85) / 7 = 9.6
        # egg_survey$Volume_Bionet<-fct_explicit_na(factor(egg_survey$Volume_Bionet), na_level = "314")
        # egg_survey$Temperature_10m<-fct_explicit_na(factor(egg_survey$Temperature_10m), na_level = "9.6")
  
    # 2016 
      # e.g. station 7.1 missing volume = so use Volume_data = 289. missing temp = average from four closest stations, missing depth = use average from 2014, 2015 and 2017
        # egg_survey$Depth_Data<-fct_explicit_na(factor(egg_survey$Depth_Data), na_level = "33.2")
        # egg_survey$Volume_Bionet<-fct_explicit_na(factor(egg_survey$Volume_Bionet), na_level = "289")
        # egg_survey$Temperature_10m<-fct_explicit_na(factor(egg_survey$Temperature_10m), na_level = "9.355")

# ==================================================================================================
#
# Calculate Egg density and daily egg production (DEP) 
  # Egg density = number of eggs sorted by Linda and Melanie set to fraction of 1 (divided by the volume of filtered water (bionet volume m^3)) and multiplied by the sampled depth  
    # Some quick data prep
  
      egg_survey_data_prep <- function(egg_survey) {
        egg_survey$STATION <- as.factor(egg_survey$STATION) 
        egg_survey$Volume_Data <- as.numeric(as.character(egg_survey$Volume_Data), digits = 2)
        egg_survey$Volume_Bionet <- as.numeric(as.character(egg_survey$Volume_Bionet), digits = 2)
        egg_survey$Depth_Data <- as.numeric(as.character(egg_survey$Depth_Data), digits = 2)
        egg_survey$Temperature_10m <- as.numeric(as.character(egg_survey$Temperature_10m), digits = 2)
        egg_survey$SUBSAMPLE <- as.numeric(as.character(egg_survey$SUBSAMPLE), digits = 5)
        egg_survey$Lat <- as.numeric(egg_survey$LATITUDE, digits = 5)
        egg_survey$Long <- as.numeric(egg_survey$LONGITUDE, digits = 5)
        egg_survey$YEAR <- as.numeric(as.character(egg_survey$YEAR))
        egg_survey$MONTH <- as.numeric(as.character(egg_survey$MONTH))
        egg_survey$DAY <- as.numeric(as.character(egg_survey$DAY))
        egg_survey<- as_tibble(egg_survey)
        glimpse(egg_survey)
      }
  egg_survey <- egg_survey_data_prep(egg_survey)

egg_survey <- egg_survey %>% mutate(
  EGGS_1_5 = STADE.1+STADE.5, # sum of stade 1 & 5 eggs by station
  N = EGGS_1_5/SUBSAMPLE, # number of eggs 1 & 5 in sample
  Nv = N/Volume_Bionet,   # by cubic metre
  Nm2 = Nv*Depth_Data) %>%  # Egg density by square metre
  arrange(YEAR,STATION) 

# ==================================================================================================  
#
# Calculate Daily Egg Production (DEP)
  # Calculate incubation time
    egg_survey$I = exp((-1.61 * log(egg_survey$Temperature_10m)) + 7.76) # unit = hours
  # New European incubation time for mackerel (Mendiola 2016)
    # egg_survey$I.2 = exp((-1.313 * log(egg_survey$Temperature_10m)) + 6.902) For egg stage 1B (roughly our Stage 1)
  # FG used different incubation eqn in resdoc 94/61 by Worley 1933. Tinc = (e^[-1.87Ln(T)+9.67])*0.0417

  # DEP (N/m2) by station
    egg_survey$DEP = egg_survey$Nm2/egg_survey$I * 24
    summary(egg_survey$DEP)

# ==================================================================================================  
# 
# Interpolate the data to the sGSL via ordinary krigging
  # First some data prep and plots 

  # Subset egg survey data for simplicity
    egg_survey_sp <- egg_survey[,c('YEAR','STATION','Long','Lat','DEP', 'I', 'Temperature_10m', 'Nm2')]

  # Simple grid for interporlation
    grille <- grille[complete.cases(grille),] 
    grille <- grille[-4,-5] 

# Denser grid made in qgis for plots
  # sGSL_grid <- readOGR(dsn = "./sGSL_grid/sGSL_grid_1NM.shp", layer = "sGSL_grid_1NM") # load shape file
  # sGSL_grid_df <- as.data.frame(sGSL_grid) # to remove the default projection of +longlat wgs83 turn to dataframe (also good to have for plotting)
  # sGSL_grid_df$long <- sGSL_grid_df$coords.x1 
  # sGSL_grid_df$lat <- sGSL_grid_df$coords.x2
  # sGSL_grid_df<- sGSL_grid_df[,-2]
  # coordinates(sGSL_grid_df) <- ~long+lat
  # proj4string(sGSL_grid_df) <- CRS_LCC_83 
  # sGSL_grid_DF <- as.data.frame(sGSL_grid_df)

  # make spatial data frames
    coordinates(egg_survey_sp) <- c("Long", "Lat")
    proj4string(egg_survey_sp) # default is no projection (ie NA)
    coordinates(grille) <- ~lon+lat

  # Assign a projection 
    proj4string(egg_survey_sp) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
    proj4string(grille) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

    # Convert back to dataframe for plotting
    egg_survey_df = data.frame(egg_survey_sp)
    egg_survey_df$egg_PA= ifelse(egg_survey_df$DEP==0,"absent","present")
  
  # Remove weird station 6.8 but will have to check for validity 
    # egg_survey_df<-egg_survey_df %>% filter(STATION != "6.8")

  # Create base map for plotting
    MY_theme_2 <- theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), title = element_text(size = 12))
    canada <- map_data("worldHires", "Canada")
    canada_map <- ggplot() + geom_polygon(data = canada,aes(x=long, y = lat, group = group),fill = "ghostwhite",color="black") + theme_bw() 
    canada_atlantic <- canada_map + 
      coord_map(xlim = c(-70, -45),  ylim = c(40, 55),projection = "lambert", parameters = c(45,50))       
    sGSL <- canada_map + 
      coord_map(xlim = c(-66.5, -60),  ylim = c(45.5, 49.5),projection = "lambert", parameters = c(46.5 ,48))  

  ## Figures
    # Plot of survey stations
      sGSL + 
        geom_point(data = egg_survey_df, aes(Long,Lat)) +  
        xlab("Longitude") + ylab("Latitude")  +
        geom_text(data = egg_survey_df, aes(x = Long, y = Lat, label = STATION, vjust = -1), size = 3) 

    # DEP map
      valuesColor <- c('black','red')
      dep.plot <- sGSL + 
        geom_point(data = egg_survey_df, aes(Long,Lat, size = DEP, colour = egg_PA)) +  
        ggtitle("Daily Egg Production (N/m^2) 2018") + 
        scale_colour_manual(values = valuesColor) + 
        scale_size_area() +
        labs(x = "Longitude", y = "Latitude") +
        MY_theme_2

    # Incubation time and temperature maps 
      inc.temp.plot <- sGSL + 
        geom_point(data = egg_survey_df, aes(Long,Lat, size = I, colour = Temperature_10m)) +  
        ggtitle("Incubation time (h) & temperature (C) 2018") + 
        scale_colour_viridis_c(option = "plasma") + 
        scale_size(range = c(0,6)) +
        labs(x = "Longitude", y = "Latitude") +
        MY_theme_2

    grid.arrange(dep.plot, inc.temp.plot)

# Fitting the variograms and ordinary krigging
source("./Functions/moy.var.Krigeage.R")

xy <- data.frame(egg_survey$LONGITUDE, egg_survey$LATITUDE)
names(xy) <- c("longitude", "latitude")
egg_survey_sp_wgs84 = SpatialPointsDataFrame(
  coords = xy,
  data = egg_survey_sp@data,
  proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

egg_survey_sp.nad83qclamb = spTransform(egg_survey_sp_wgs84, CRSobj = CRS("+init=epsg:32198"))
grille.nad83qclamb = spTransform(grille, CRSobj = CRS("+init=epsg:32198"))

# Remove outliers
my.vgm <- variogram(DEP ~ 1, data = egg_survey_sp.nad83qclamb[egg_survey_sp.nad83qclamb$DEP < quantile(egg_survey_sp.nad83qclamb$DEP, .99), ]) 
plot(my.vgm)

# Normalise the variogram by dividing the data used by the variance
my.vgm$gamma <- my.vgm$gamma / var(egg_survey_sp.nad83qclamb[egg_survey_sp.nad83qclamb$DEP<quantile(egg_survey_sp.nad83qclamb$DEP,.99),]$DEP)
plot(my.vgm)

# Set initial variogram parameters by eye 
model.i = vgm(psill=0.8, model="Exp", range=100000, nugget=0.2)
plot(my.vgm, model=model.i, plot.numbers=T, pch=16)

# Fit variogram from initial values (If no convergence, check tolerence of outlier threshhold. Above it is set to .99 in the quantile function)
model.f = fit.variogram(my.vgm, model=model.i,fit.ranges = T, fit.method=7) 
plot(my.vgm, model=model.f, plot.numbers=T, pch=16)

# Once the variogram is fit, rescale the semivariogram to the true values by multiplying by the variance of the values (undoing earlier division)
model.f$psill <- model.f$psill * var(egg_survey_sp$DEP)
model.f

# Empty data frames for later use
pred=data.frame()  # dataframe with predictions krigging
moy_var=data.frame() # dataframe wih mean and variance of krigging predictions of whole area
egg_survey$pred=NA  

# Interpolate DEP across DEP via ordinary krigging (using fitted semivariogram)
krige = krige(DEP ~ 1, egg_survey_sp.nad83qclamb[egg_survey_sp.nad83qclamb@data$YEAR == '2018', ],
              grille.nad83qclamb,
              model.f,
              nmin = 16)
  
pred1 = data.frame(YEAR = 2018,depth=grille.nad83qclamb@data$prof,
                     pred = krige@data$var1.pred,
                     var = krige@data$var1.var,
                     LONG = grille.nad83qclamb@coords[,1],LAT = grille.nad83qclamb@coords[,2])
  pred=rbind(pred1,pred)
  
  spplot(krige["var1.pred"])
 krigdf <- as.data.frame(krige) 

 
 # Krig plot
krig_plot1 <- ggplot(data = krigdf, aes(x = lon, y = lat))  + 
  stat_summary_2d(data=krigdf,aes(x = lon, y = lat,z=var1.pred),fun = mean, binwidth = c(0.03, 0.03)) +
  geom_polygon(data = canada, aes(x=long, y = lat, group = group), fill = "ghostwhite",color="black")  +
  coord_map(xlim = c(-66.5, -60),  ylim = c(45.5, 49.5),projection = "lambert", parameters = c(46.5 ,48)) +
  theme_bw() + 
  geom_point(data = egg_survey_df, aes(Long,Lat, size = DEP), pch = 21) +  
  scale_fill_viridis_c(option="viridis", begin = 0.1) + 
  labs(x = "Longitude",y = "Latitude", title = "Interpolated Daily Egg Production")
  

krig_plot1_var <- ggplot(data = krigdf, aes(x = x, y = y))  + 
  stat_summary_2d(data=krigdf,aes(x = long, y = lat,z=var1.var),fun = mean, binwidth = c(0.4, 0.4))+
  geom_polygon(data = canada, aes(x=long, y = lat, group = group), fill = "ghostwhite",color="black")  +
  coord_map(xlim = c(-66.5, -60),  ylim = c(45.5, 49.5),projection = "lambert", parameters = c(46.5 ,48)) + 
  theme_bw() + 
  geom_point(data = egg_survey_df, aes(Long,Lat, size = DEP), pch = 21) +  
  scale_fill_viridis_c(option="viridis") + 
  labs(x = "Longitude",y = "Latitude", title = "Daily Egg Production Variance", tag = "B") 
krig_compare<-grid.arrange(krig_plot2, krig_plot2_var)

# krige with denser grid

# krige2=krige(DEP~1, egg_survey_sp, sGSL_grid_df, model.f, nmin=16)
# spplot(krige2["var1.pred"])
# krigdf2 <- as.data.frame(krige2) 
# 
#  krig_plot2 <- ggplot(data = krigdf2, aes(x = x, y = y))  + 
#     stat_summary_2d(data=krigdf2,aes(x = long, y = lat,z=var1.pred),fun = mean, binwidth = c(0.05, 0.05)) +
#     geom_polygon(data = canada, aes(x=long, y = lat, group = group), fill = "ghostwhite",color="black") +
#     coord_map(xlim = c(-66.5, -60),  ylim = c(45.5, 49.5),projection = "lambert", parameters = c(46.5 ,48)) +
#     theme_bw() + 
#     geom_point(data = egg_survey_df, aes(Long,Lat, size = DEP), pch = 21) +  
#     scale_fill_viridis_c(option="viridis", name = expression(Eggs/m^2)) + 
#     labs(x = "Longitude",y = "Latitude",title = "Predicted Daily Egg Production", tag = "A") 
#    
#  
#   krig_plot2_var<- ggplot(data = krigdf2, aes(x = x, y = y))  + 
#    stat_summary_2d(data=krigdf2,aes(x = long, y = lat,z=var1.var),fun = mean, binwidth = c(0.4, 0.4))+
#    geom_polygon(data = w2hr, aes(x=long, y = lat, group = group), fill = "lightgrey",color="black")  +
#    coord_fixed(xlim = c(-66.5, -60),  ylim = c(45.5, 49)) + theme_bw() + 
#    geom_point(data = egg_survey_df, aes(Long,Lat, size = DEP), pch = 21) +  
#    scale_fill_viridis_c(option="viridis") + 
#    labs(x = "Longitude",y = "Latitude", title = "Daily Egg Production Variance", tag = "B") 
#  krig_compare<-grid.arrange(krig_plot2, krig_plot2_var)

 # choose year
 egg_survey_df <- egg_survey_df %>% filter(YEAR == 2018, !is.na(DEP))
 egg_survey_sp <- egg_survey_df
 coordinates(egg_survey_sp) <- c("Long", "Lat")
 proj4string(egg_survey_sp) # default is no projection (ie NA)
 proj4string(egg_survey_sp) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
 
 moy_var = moy.var.Krigeage(
    xCoords = coordinates(egg_survey_sp), 
    xValues = egg_survey_df$DEP, 
    kCoords = coordinates(grille), 
    vario = model.f)
  # moy_var1=as.data.frame(moy_var1)
  # moy_var1$YEAR
  # moy_var1$z.noneg=mean(pred1[!pred1$pred<0,'pred'])
  # moy_var=rbind(moy_var,moy_var1)
  # pred$neg=ifelse(pred$pred<0,'-','') # Negative predictions indicated on map
  moy_var$sd = sqrt(moy_var$VarKrigeage)
  moy_var$CV = (moy_var$sd/moy_var$Z)*100

# cor and RSS (as in table 6)
  # sampling stations are not located on grid coordinates:  find nearest predicted value
    # egg_survey_df$res=egg_survey_df$DEP-egg_survey_df$pred
    # FIT=ddply(egg_survey_df,c('YEAR'),summarise,COR=cor(DEP,pred),RSS=sum(res^2))

# summary values for whole area
DEP2 = moy_var   # This is the same as mean(pred$pred). 
DEP = data.frame(YEAR=unique(egg_survey$YEAR),Mean_Daily_Egg_Production=DEP2$Z,Daily_Egg_Production_Var=DEP2$VarKrigeage) #z.noneg if no negatifs

# ==================================================================================================  
# calculation of the proportion of the eggs spawned daily 

# median date of the survey
#load("./Data/Rdata/egg_survey.RData")

# DATE AND JULIAN DAY
egg_survey <- egg_survey %>% 
  mutate(DATE = lubridate::ymd(paste(YEAR, MONTH, DAY)), doy = as.numeric(format(DATE, "%j")))
date.med<-data.frame(year = sort(unique(egg_survey$YEAR)))
date.med$doy <- round(median(egg_survey$doy))
#date.med$doy <- egg_survey %>% group_by(YEAR) %>% summarise(doy = round(median(doy)))

# GSI logistic curve   
  # EVB, Martin, and I discussed this and decided to remove immature individuals (stages 0-3) from the model fitting. Apparently Francois G did this back in the day.
  # update nov 2019. RESDOC from Francois G. states very clearly that only stage 5 mature females should be used due to batch fecundity problems!!!!! update analyses
GSI <- maq_carbio %>% 
  filter(MATURE > 3) %>% 
  mutate(GSI = (WEIGHT_GONAD/WEIGHT)*100, 
         DATE = lubridate::ymd(paste(YEAR, MONTH, DAY)), 
         doy = as.numeric(format(DATE, "%j"))) %>% 
  filter(!is.na(GSI))

# in future use purrr and modelr and broom
logistic.out <- nls(GSI ~ y0 + (a/(1+(doy/x0)^b)),
                 data = GSI, start = list(y0 = 0.4451, x0 = 172.4085, a = 11.3252,b = 32.0552))
#if singular gradient 
 library(nls2)
 logistic.out=nls2(GSI ~ y0 + (a/(1+(doy/x0)^b)),
                  data = GSI, start = list(y0 = 0.5, x0 = 160, a = 15,b = 30),algorithm = "brute-force")

library(broom)
co <- data.frame() #coefficients of fit
spawning.prob <- data.frame(year=sort(unique(egg_survey$YEAR)),prob=NA,peak.day=NA) # spawning proportion at data, for SSB calculation
N <- nrow(GSI)
co <- tidy(logistic.out) %>% dplyr::select(term, estimate) 

logi.plot <- ggplot(GSI,aes(x=doy,y=GSI)) +
  geom_point(aes(x = doy, y = GSI,  colour = factor(MATURE))) +
  stat_function(fun=function(x) co$estimate[1] + (co$estimate[3]/(1+(x/co$estimate[2])^co$estimate[4])), col='red', size=1.5) +
  scale_color_viridis_d(direction = -1) +
  scale_x_continuous(breaks = seq(150,300,10)) +
  labs(x = "Jour / Day", y = "IGS / GSI", colour = "Stade de Maturité / Maturity Stage" ) +
  theme_bw() +
  theme(legend.position="top")

# 2018 still fucking up and creating singular gradient... so trying by eye
GSI2<-GSI %>% filter(DIVISION == "4T")
ggplot(GSI,aes(x=doy,y=GSI)) +
  geom_point(aes(x = doy, y = GSI,  colour = factor(MATURE))) +
  stat_function(fun=function(x) 0.45 + (14/(1+(x/169)^35)), col='red', size=1.5) +
  scale_color_viridis_d(direction = -1) +
  scale_x_continuous(breaks = seq(150,300,10)) +
  labs(x = "Jour / Day", y = "IGS / GSI", colour = "Stade de Maturité / Maturity Stage" ) +
  theme_bw() +
  theme(legend.position="top")

d1 <- data_frame(day=seq(100.5,350.5,1))
d1$pred <- co$estimate[1] + (co$estimate[3]/(1+(d1$day/co$estimate[2])^co$estimate[4]))
d2 <- data.frame(day=seq(101,350,1))
for(j in 1:(nrow(d1)-1)){d2[j,'slope']=d1[j+1,'pred']-d1[j,'pred']}
d2$prob <- d2$slope/sum(d2$slope,na.rm=T)
# 2018
d1 <- data_frame(day=seq(100.5,350.5,1))
d1$pred <- 0.45 + (14/(1+(d1$day/169)^35))
d2 <- data.frame(day=seq(101,350,1))
for(j in 1:(nrow(d1)-1)){d2[j,'slope']=d1[j+1,'pred']-d1[j,'pred']}
d2$prob <- d2$slope/sum(d2$slope,na.rm=T)

#beginning and end of spawning
upper <- mean(d2$prob) + (sd(d2$prob)/sqrt(length(d2$prob)))
lower <- mean(d2$prob) - (sd(d2$prob)/sqrt(length(d2$prob)))
start <- d2 %>% filter(day<160) %>% slice(which.min(abs(upper - prob)))
end <- d2 %>% filter(day>160) %>% slice(which.min(abs(lower - prob)))
# peak spawning 
spawning.prob[spawning.prob$year==2018,'peak.day'] = d2[d2$prob==max(d2$prob),"day"]
spawning.prob[spawning.prob$year==2018,'prob'] = d2[d2$day==date.med$doy,"prob"]

# 2014
as.Date(spawning.prob$peak.day, origin=as.Date("2014-01-01")) # JUNE 16TH WAS THE PEAK DAY in 2015
as.Date(start$day, origin=as.Date("2014-01-01")) # MAY 25 start 2015
as.Date(end$day, origin=as.Date("2014-01-01")) # July 12th end
as.Date(date.med$doy, origin=as.Date("2014-01-01")) # JUNE 16TH WAS median mission date
end$day-start$day # range 48

# 2015
as.Date(spawning.prob$peak.day, origin=as.Date("2015-01-01")) # JUNE 22TH WAS THE PEAK DAY in 2015
as.Date(start$day, origin=as.Date("2015-01-01")) # JUNE 2 start 2015
as.Date(end$day, origin=as.Date("2015-01-01")) # July 16th end
as.Date(date.med$doy, origin=as.Date("2015-01-01")) # JUNE 17TH WAS median mission date
end$day-start$day # range 44

# 2016
as.Date(spawning.prob$peak.day, origin=as.Date("2016-01-01")) # JUNE 12th in 2016
as.Date(start$day, origin=as.Date("2016-01-01")) # May 14 start 2015
as.Date(end$day, origin=as.Date("2016-01-01")) # July 18th end
as.Date(date.med$doy, origin=as.Date("2016-01-01")) # JUNE 18TH WAS median mission date
end$day-start$day # range 65

# 2017
as.Date(spawning.prob$peak.day, origin=as.Date("2017-01-01")) # JUNE 24TH WAS THE PEAK DAY in 2017
as.Date(start$day, origin=as.Date("2017-01-01")) # JUNE 7th start 2017
as.Date(end$day, origin=as.Date("2017-01-01")) # July 14th end
as.Date(date.med$doy, origin=as.Date("2017-01-01")) # JUNE 14TH WAS median mission date
end$day-start$day # range

# 2018
as.Date(spawning.prob$peak.day, origin=as.Date("2018-01-01")) # JUNE 19TH WAS THE PEAK DAY in 2018
as.Date(start$day, origin=as.Date("2018-01-01")) # JUNE 1st start 2018
as.Date(end$day, origin=as.Date("2018-01-01")) # July 09th end
as.Date(date.med$doy, origin=as.Date("2018-01-01")) # JUNE 22TH WAS median mission date
end$day-start$day # range

prob.plot <-
  ggplot(d2,aes(x=day,y=prob))+geom_line(size = 1.5)+
  scale_x_continuous(breaks=seq(0,305,10),limits=c(130,220)) + theme_bw() + 
  ggtitle("Pic de ponte 2018 / Peak spawning 2018") +  
  annotate("text", x = spawning.prob$peak.day+5, y =max(d2$prob)+0.001, label = "06/19") + # peak
  annotate("text", x = start$day-6, y = start$prob +0.001, label = "06/01") +             # start
  annotate("text", x = end$day+6, y = end$prob +0.001, label = "07/09") +                 # end
  annotate("text", x = date.med$doy+5, y = spawning.prob$prob, label = "06/22") +         # mission
  annotate("point", x = start$day, y = start$prob, colour = "blue", size = 5) +
  annotate("point", x = end$day, y = end$prob, colour = "blue", size = 5) +
  annotate("point", x = spawning.prob$peak.day, y = max(d2$prob), colour = "blue", size = 5) +
  annotate("point", x = date.med$doy, y = spawning.prob$prob, colour = "red", size = 5) +
  annotate("segment", x = start$day, xend = end$day, y = 0, yend = 0,  colour = "blue", size = 3) +
  annotate("text", x = spawning.prob$peak.day, y =-0.002, label = "Durée = 38 jours / Duration = 38 days") +
  labs(y = "Proportion", x = "Jour / Day") +
  ylim(c(-0.005,0.06)) +
  MY_theme_2

# ==================================================================================================  
# Calculate raised DEP and Total(Annual) Egg Production (TEP)
#########
# Raise DEP to entire sampled area (sGSL area established in previous studies but will need to be updated)
A <- 6.945e+10   
DEP$DEP_sGSL = DEP$Mean_Daily_Egg_Production * A
DEP$Tot_Annual_Egg_Prod = DEP$Mean_Daily_Egg_Production * A / spawning.prob$prob

# ==================================================================================================
#
# Calculate the parameters of the SSB equation: B = (P*A*W)/(S*F*R*10^6)
  
  # P = Mean daily egg productio by station
  # A = Area (m^2) of the zone sazmpled (6.945 * 10^10 m^2)
  # W = Mean weight (g) of a fish
  # S = Proportion of eggs spawned at ith median date of a survey
  # F = Fecundity of females (Pelletier 1986)
  # C = 10^6 =  factor for converting grams into tonnes

  # Empty Dataframe
    SSB=data.frame(year=unique(DEP$YEAR),P=NA,A=NA,W=NA,S=NA,FEC=NA,R=NA,convers=NA,SSB=NA)

SSB$P <- DEP$Mean_Daily_Egg_Production                    # N/m2 from krigging
 <- A                                     # m2
maqWF <- maq_carbio %>% filter(!is.na(WEIGHT), SEX =="F", !is.na(WEIGHT_GONAD)) # females
maqWMF <- maq_carbio %>% filter(!is.na(WEIGHT), !is.na(WEIGHT_GONAD), SEX != "I", !is.na(SEX)) # males and females # note for future. Do I just want the weight of spawning M=6 or mature fish M>2?
SSB$W <- mean(maqWMF$WEIGHT-maqWMF$WEIGHT_GONAD)                               # mean weight fish (g)
SSB$S <- spawning.prob$prob                                  # proportion of eggs spawned at median date survey
maq_mat5 <- maq_carbio %>% filter(MATURE == 5, SEX =="F", !is.na(WEIGHT_GONAD))
maq_mat5 <- maq_mat5 %>% mutate(FEC = mean(10^(4.32 + 0.75 * log10(maq_mat5$WEIGHT_GONAD))))
SSB$FEC = mean(10^(4.32 + 0.75 * log10(maq_mat5$WEIGHT_GONAD))) # fecundity of females (pelletier uses log10!!!)....shit
maqF <- maq_carbio %>% dplyr::filter(SEX =="F") 
maqMF <- maq_carbio %>% dplyr::filter(SEX %in% c("M","F"))
SSB$R=length(maqF$YEAR)/length(maqMF$YEAR) # proportion of females by number
SSB$convers=10^6                            # grams to tonnes
 
SSB <- SSB %>% mutate(SSB = round((P*A*W)/(S*FEC*R*convers), digits = 0))
SSB




# ================================================================================================== 
    # MAP PROJECTIONS IN R
  # ================================================================================================== 
  
  # see: http://rspatial.org/index.html, 
  # https://www.nrcan.gc.ca/earth-sciences/geomatics/geodetic-reference-systems/18766,
  # https://epsg.io/102002
  # https://proj4.org/index.html
  # and www.cef-cfr.ca/uploads/Reference/ProjectionCartoNum.pdf for projections used in quebec/canada
  # ================================================================================================== 
  install.packages(c("maps", "mapdata"))
  # Load packages
  packages = c('tidyverse','readr','magrittr', 'plyr', 'lubridate','maps','mapdata','gridExtra','sp','gstat','raster','rgeos','rgdal','PBSmapping','geoR')
  invisible(lapply(packages, function(x) {if (!require(x, character.only = T)) {install.packages(x);require(x)}}))
  # ================================================================================================== 
  
  # Common projections in proj.4 format 
  # Used in the R packages sp and rgdal among others)
  # Canada Lambert Conformal Conic Modern ESRI:102002
  CRS_LCC_83 <- CRS("+proj=lcc +lat_0=46.5 +lat_1=48 +lat_2=50 +lon_0=-70 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs") #ESRI:102002
  # Canada Lambert Conformal Conic Old
  CRS_LCC_27 <- CRS("+proj=lcc +lat_0=46.5 +lat_1=48 +lat_2=50 +lon_0=-70 +x_0=0 +y_0=0 +ellps=clrk66 +datum=NAD27 +units=km +no_defs")
  # Pseudo Mercator (used by google maps, open street view, bing) EPSG:3857 
  CRS_Google <- CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=km +nadgrids=@null +wktext  +no_defs")
  # World Mercator EPSG:3395
  CRS_Mercator <- CRS("+proj=merc +lon_0=-70 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs ")
  # ================================================================================================== 
  
  # Convert latitude and longitude to spatial points dataframe and compare different projections
  lat <- stations$Lat
  long <- stations$Long
  latlong <- data.frame(long, lat)
  latlong27 <- data.frame(long, lat)
  latlong83 <- data.frame(long, lat)
  coordinates(latlong) <- c("long", "lat")
  coordinates(latlong27) <- c("long", "lat")
  coordinates(latlong83) <- c("long", "lat")
  proj4string(latlong) # default is no projection (ie NA)
  # Assign a default projection (eg CRS_Mercator)
  proj4string(latlong) <- CRS_Mercator # WGS 84 Mercator
  proj4string(latlong27) <- CRS_LCC_27 # WGS 84 Mercator
  proj4string(latlong83) <- CRS_LCC_83 # WGS 84 Mercator
  # Other projections
  stations_LCC27 <- spTransform(latlong, CRS_LCC_27)
  stations_LCC83 <- spTransform(latlong, CRS_LCC_83)
  stations_LCCgoo <- spTransform(latlong, CRS_Google)
 
  # Plot the results
  par(mfrow=c(2,3))
  plot(stations$Long,stations$Lat, main="Raw data", cex.axis=.95) # cartesian coordinates
  plot(latlong, axes=TRUE, main="MERCATOR", cex.axis=.95)
  plot(stations_LCC27, axes=TRUE, main="NAD27", cex.axis=.95)
  plot(stations_LCC83, axes=TRUE, main="NAD83", cex.axis=.95)
  plot(stations_LCCgoo, axes=TRUE, main="GOOGLE", cex.axis=.95)
  dev.off()
  # PLOT POINTS AGAISNT MAP TO GET BETTER IDEA
  canada <- map_data("worldHires", "Canada")
  canada_map <- ggplot() + geom_polygon(data = canada,aes(x=long, y = lat, group = group),fill = "ghostwhite",color="black") + theme_bw() 
  canada_atlantic <- canada_map + 
    coord_map(xlim = c(-70, -45),  ylim = c(40, 55),projection = "lambert", parameters = c(45,50))       
  sGSL <- canada_map + 
    coord_map(xlim = c(-66.5, -60),  ylim = c(45.5, 49.5),projection = "lambert", parameters = c(46.5 ,48))   
    

# even higher resolution
world <- shapefile("./Data/Egg_survey/GSHHS_shp/f/GSHHS_f_L1.shp")
sGSL <- crop(world, extent(-68, -59.7, 45, 50))
plot(sGSL)
writeOGR(sGSL, dsn = '.', layer = 'sGSL', driver = "ESRI Shapefile")
sGSL <- shapefile("./sGSL.shp")
 