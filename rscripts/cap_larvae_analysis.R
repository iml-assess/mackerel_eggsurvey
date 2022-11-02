############################  Capelin sGSL larvae analysis  ##############################################
##  
##  see cap_larvae.R to see how this data was created
##  
##  Andrew
##  Things to do: 02-02-2020
##  Load each raw data file in one at a time and clean and merge 
##  Time, date, volume, sample depth, choice of sous-echantillion column should all be checked
##  Krigging and negative binomial glm give similar results and further investigation of krigging results suggest very poor spatial structure and thus independance of sample locations. Future research can investigate potential bottom up and/or top down drivers as well as particle transport etc
##  - September 2020 redo a couple of analyses
##  - try with cap_larvae.Rdata and with SGSLTAXON
############################  Load and clean data   ##############################################
load("./data/larvae/cap_larvae.Rdata")
load("S:/Pélagiques/Plancton/BD Ichtyop - Patrick Ouellet DB/SGSL 1982-2015.rda")
load("./cap_larvae_lengths.Rdata")
capelin_spawning <- read_csv(file = "./capelin_spawning.csv")
sgsl_sst_incubation <- read_csv(file = "./sgsl_sst_incubation.csv")

library(ciTools)
glimpse(cap_larvae)
############################  Choose a model  ##############################################
# goal is to create a standardized larval index for capelin but a standardized index in often just the N/area or volume 
# Response variable is either Nb (counts, a discrete variable), Nv (concentration a continuous variable Nb/1000m^-3), and Nm2 (density a continous variable Nb/m^-2) 
# Count data will be treated with glm poisson (mean = var), glm negative binomial (if overdispersion), hurdle (2 step, bernouilli then poisson), zeroinfl(if data is believed to be zero inflated)
# Count data will not be log transformed (see debates in O'Hara and Kotze 2010; Ives 2015; and St. Pierre et al., 2018)
# Countinous data will be treated with glm or glmm(mixed aka hierachical model) with a distribution that makes sense for a lot of zeros ie tweedie distribution or hurdle model
# From a causal inference point of view, larval density (i.e. relative abundance estimate) variation across our sampling sites could be influenced by a number of unobserved/latent processess
# Such as: 1) spatial autocorrelation, 2) temporal autocorrelation and/or the effect of sampling date and/or the interaction between the two, and 3) an unobserved counfounding variable (i.e. environment - but we don't want to block that because that variation will be examined in subsequent steps)
# the goal of this analysis is therefore to compare a few different standardizations and see to what extent spatial and temporal factors should be considered

############################  Look at data  ##############################################

# Distributions
# Does the data have a lot of zeros? -- YES
cap_larvae %>% ggplot(aes(log(Nm2))) + geom_density() 
100*sum(cap_larvae$Nm2 == 0, na.rm = T)/nrow(cap_larvae) # ~76 percent are zeroes! if considering all trips etc. 

# Is the data overdispersed (ie var > mean) ? --- YES
summary(cap_larvae$Nm2, na.rm = T) # data is overdispersed
var(cap_larvae$Nm2, na.rm = T) # variance is far larger than the mean var = ~711

# Is there a spatial gradient?  --- YES, most non zero values concentrated at higher latitudes and further west
cap_larvae %>% 
  dplyr::filter(presence == "present") %>% 
  ggplot(aes(x=lon, y=lat) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = TRUE) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    legend.position='none'
  ) + scale_fill_viridis_c() +
  labs(
    x = "longitude",
    y = "latitude"
  ) 

library(plotly)
library(MASS)
b
# Compute kde2d
kd <- with(cap_larvae, MASS::kde2d(lon, lat, n = 50))
# Plot with plotly
plot_ly(x = kd$x, y = kd$y, z = kd$z) %>% add_surface()

# Compute kde2d
kd <- with(df, MASS::kde2d(as.numeric(as.character(year)), doy, n = 50))
# Plot with plotly
plot_ly(x = kd$x, y = kd$y, z = kd$z) %>% add_surface()

cap_larvae$trip_no <- as.factor(cap_larvae$trip_no)

# trip one and trip 2 have similar observer error/detection prob? not bad
cap_larvae %>% count(trip_no)
cap_larvae %>% 
  group_by(year, trip_no) %>% 
  dplyr::summarise(Nm2 = mean(Nm2, na.rm = T)) %>%
  ggplot(aes(year, log(Nm2),colour = trip_no)) +
  geom_line(size = 2) +
  geom_point(data = cap_larvae,aes(year,log(Nm2) ), size = 0.5)

a <- cap_larvae %>% dplyr::filter(trip_no == "1")
b <- cap_larvae %>% dplyr::filter(trip_no == "2")

cap_larvae %>% dplyr::filter(hour<25) %>% summary()
summary(a)
summary(b)

a<-cap_larvae %>% 
  ggplot(aes(year, as.vector(scale(log(Nm2+0.0000000001))), colour = trip_no, group = trip_no)) + 
  stat_summary(fun.data = "mean_cl_boot", alpha = 0.9, size = 0.8, position = position_dodge(width = 0.5)) +
  scale_colour_viridis_d(begin = 0.05, end = 0.6, direction = -1)+ geom_hline(yintercept = 0, colour= "red") + 
  annotate("text", x = 1995, y = 1, label = "Trip 1 median date = 171; Trip 2 median date = 180") +
  theme_minimal(base_size = 14) +
  labs(x = "", y = "", colour = "Trip No.") +
  theme(legend.position = "top", legend.margin = margin(0,0,0,0))

cap_larvae %>% dplyr::filter(!is.na(hour), hour < 25) %>% 
  ggplot(aes(hour, as.vector(scale(log(Nm2+0.0000000001))))) + 
  stat_summary(fun.data = "mean_cl_boot", alpha = 0.9, size = 0.8, position = position_dodge(width = 1.1)) +
  scale_colour_viridis_d(begin = 0.05, end = 0.6)+ geom_hline(yintercept = 0, colour= "red") + 
  theme_minimal(base_size = 18) +
  labs(x = "Hour", y = "std. log larval density") 

# po ouellette compare
library(readr)
PO_nwgsl <- read_csv("data/larvae/PO_nwgsl.csv") %>% transmute(year = Year, dens = Ouelletetal2013_fall_lardensity, std_dens = std_log_fall_larvdensity)
 


unique(NGSLTaxon$CodeM)
poNWGSL <- NGSLTaxon %>% dplyr::filter(Taxon == "Mallotus villosus", CodeM %in% c(
  "IML99-43","IML00-52","IML01-61","IML02-68","IML03-54","IML04-54","IML05-71","IML06-60","IML07-49","IML08-57","IML09-63"
))

poNWGSL <- poNWGSL %>% 
  transmute(year = Année, doy = as.numeric(DOY), Nm2 = Nb_m2, log_Nm2 = log(Nm2+0.0000000001), std_dens = as.vector(scale(log_Nm2))) %>% dplyr::filter(!is.na(Nm2))

capdens2 <- cap_larvae %>% dplyr::filter(trip_no =="1") %>% dplyr::select(year, date, Nm2) 
capdens2 %<>% mutate(doy = yday(date), survey = "sGSL_Spring", log_Nm2 = log(Nm2+0.0000000001), std_dens = as.vector(scale(log_Nm2)))
summary(capdens2)
poNWGSL$survey = "nwGSL_Fall"
cap_compare <- bind_rows(capdens2,poNWGSL)

summary(cap_dens)
summary(poNWGSL)

b<- cap_compare %>% 
  ggplot(aes(year, std_dens, colour = survey)) + 
  stat_summary(fun.data = "mean_cl_boot", alpha = 0.9, size = 0.8, position = position_dodge(width = 0.5)) +
  scale_colour_viridis_d(begin = 0.05, end = 0.6)+ geom_hline(yintercept = 0, colour= "red") + 
  annotate("text", x = 1995, y = 1, label = "nwGSL median date = 305; sGSL median date = 171") +
  theme_minimal(base_size = 14) +
  labs(x = "", y = "std. log larval density", colour = "Survey") +
  theme(legend.position = "top", legend.margin = margin(0,0,0,0))

NL <- NGSLTaxon %>% dplyr::filter(CodeM %in% c("IML04-49", "IML05-62", "IML07-38", "IML08-50", "IML09-54"),Taxon == "Mallotus villosus")
NL <- NL %>% 
  transmute(year = Année, doy = as.numeric(DOY), Nm2 = Nb_m2, log_Nm2 = log(Nm2+0.0000000001), std_dens = as.vector(scale(log_Nm2))) %>% dplyr::filter(!is.na(Nm2))
NL$survey = "eGSL_Summer"
cap_compare2 <- bind_rows(capdens2,NL)
summary(NL)

cap_compare2 %>% 
  ggplot(aes(year, std_dens, colour = survey)) + 
  stat_summary(fun.data = "mean_cl_boot", alpha = 0.9, size = 0.8, position = position_dodge(width = 1.1)) +
  scale_colour_viridis_d(begin = 0.05, end = 0.6)+ geom_hline(yintercept = 0, colour= "red") + 
  annotate("text", x = 1995, y = 1, label = "eGSL median date = 199; sGSL median date = 170") +
  theme_minimal(base_size = 18) +
  labs(x = "Year", y = "std. log larval density", colour = "Trip No.") 

c<- cap_compare2 %>% 
  ggplot(aes(year, std_dens, colour = survey)) + 
  stat_summary(fun.data = "mean_cl_boot", alpha = 0.9, size = 0.8, position = position_dodge(width = 0.5)) +
  scale_colour_viridis_d(begin = 0.05, end = 0.6)+ geom_hline(yintercept = 0, colour= "red") + 
  annotate("text", x = 1995, y = 1, label = "eGSL median date = 199; sGSL median date = 171") +
  theme_minimal(base_size = 14) +
  labs(x = "Year", y = "", colour = "Survey") +
  theme(legend.position = "top", legend.margin = margin(0,0,0,0)) 


library(cowplot)
plot_grid(a,b,c, labels = c("A","B","C"), ncol = 1)

NL$Survey = "eGSL_Summer"
df$Survey = "sGSL_Spring"
poNWGSL$Survey = "nwGSL_Fall"
NL2 <- NL %>%  transmute(year = Année, lat = Latitude, lon = Longitude, Survey = Survey)
cap_larvae2 <- df %>% transmute(year = as.numeric(as.character(year)), lat = lat, lon = lon, Survey = Survey)
poNWGSL2 <- poNWGSL %>%  transmute(year = Année, lat = Latitude, lon = Longitude, Survey = Survey)
mapdata<- bind_rows(NL2, cap_larvae2, poNWGSL2)



NL2 <- NGSLTaxon %>% dplyr::filter(CodeM %in% c("IML04-49", "IML05-62", "IML07-38", "IML08-50", "IML09-54")) %>% transmute(lat = Latitude, lon = Longitude) 
poNWGSL2 <- NGSLTaxon %>% dplyr::filter(CodeM %in% c("IML99-43","IML00-52","IML01-61","IML02-68","IML03-54","IML04-54","IML05-71","IML06-60","IML07-49","IML08-57","IML09-63")) %>% transmute(lat = Latitude, lon = Longitude) 
SGSLMeta %<>% transmute(lat = Latitude, lon = Longitude) 
aaa<- bind_rows(NL2, poNWGSL2,SGSLMeta) %>% dplyr::filter(!is.na(lat))


leaflet(aaa) %>%
  addProviderTiles(providers$Esri.OceanBasemap) %>%  
  addHeatmap(lng = ~lon, lat = ~lat, radius = 8)
pal <- colorNumeric(
  palette = "viridis", reverse = F, 
  domain = aaa$Survey)
aaa$survey <- as.numeric(aaa$)

leaflet(aaa) %>%
  addProviderTiles(providers$CartoDB.PositronNoLabels) %>%  
  setView(-64, 48, zoom = 6.4) %>% 
  addCircleMarkers(data = mapdata, 
                   lng = ~ lon, 
                   lat = ~ lat,
                   color = ~ "red",
                   radius = 0.5) 

sGSL + 
  geom_point(data = egg_survey_df, aes(Long,Lat)) +  
  xlab("Longitude") + ylab("Latitude")  +
  geom_text(data = egg_survey_df, aes(x = Long, y = Lat, label = STATION, vjust = -1), size = 3) 

# how about that repeat station in 1991 or just day and night in general
a<-cap_larvae %>% dplyr::filter(year ==1991,station %in% c(7.5))  # stations 7.5,7.4,8.4
a %>% ggplot(aes(hour, Nm2,colour = factor(station) )) + geom_point(size = 4) + facet_wrap(vars(station)) # looks pretty consistent but mostly zeros!

# capelin_larvae_sgsl  %>% ggplot(aes(doy, fill = presence, alpha = 0.5)) + geom_density() # normal distribution but with a few modes
# capelin_larvae_sgsl %>% dplyr::filter(trip_no == 1) %>% ggplot(aes(doy, fill = presence, alpha = 0.5)) + geom_density() + facet_wrap(vars(year), scales = "free_y")
# capelin_larvae_sgsl %>% dplyr::filter(trip_no == 1) %>% ggplot(aes(lat, fill = presence, alpha = 0.5)) + geom_density() + facet_wrap(vars(year), scales = "free_y")
# capelin_larvae_sgsl %>% dplyr::filter(trip_no == 1) %>% ggplot(aes(lon, fill = presence, alpha = 0.5)) + geom_density() + facet_wrap(vars(year), scales = "free_y")

############################  Subset data  ##############################################

# note that I am keeping only trip 1 because trip 2 was discontinued ~ halfway through the time series, and
# I am keeping both babord and tribord samples as they will add replicates and increase statistical power in theory
df <- 
  cap_larvae %>% 
  dplyr::filter(trip_no == 1, 
                !is.na(volume), 
                !is.na(Nb), 
                !is.na(lat), 
                !is.na(lon)) %>% 
  mutate(lon = lon, 
         year = as.factor(year))
# remove stations in northumberland strait and outside of the normal 65 stations
df %<>% 
  dplyr::filter(station < 13, 
                !station %in% c(5,6,7,8,9))  
df %<>% mutate(std_Nm2 = sqrt(Nm2))
df %<>% mutate(std_Nm2_2 = std_Nm2/max(std_Nm2,na.rm = T))
df$Year <- as.numeric(as.character(df$year))
# df %<>% mutate(std_Nm2_2 = ifelse(std_Nm2 <0, 0, std_Nm2))
df1 <- df %>% dplyr::filter(Year<2006) %>% dplyr::select(year, std_Nm2, lat, lon,std_Nm2_2)
df2 <- df %>% dplyr::filter(Year>2005) %>% dplyr::select(year, std_Nm2, lat, lon,std_Nm2_2)

sGSL + 
  geom_point(data = df %>% dplyr::filter(!is.na(Nm2)), aes(lon,lat, colour = std_Nm2, size = std_Nm2_2, alpha = 0.7)) +  
  xlab("Longitude") + ylab("Latitude") +
  scale_colour_gradient2(low = "blue",  high = "red") + 
  scale_size_area()+
  facet_wrap(vars(year), ncol = 5) + theme(legend.margin = margin(0,0,0,0), legend.position = "none",plot.margin = margin(0,0,0,0))

sGSL + 
  geom_point(data = df1, aes(lon,lat, colour = std_Nm2, size = std_Nm2_2, alpha = 0.5)) +  
  xlab("Longitude") + ylab("Latitude") +
  scale_colour_gradient(low = "blue", high = "red") + 
  scale_size_binned_area()+
  facet_wrap(vars(year))
  

 
sGSL + 
  geom_point(data = df2, aes(lon,lat, colour = std_Nm2, size = std_Nm2_2, alpha = 0.5)) +  
  xlab("Longitude") + ylab("Latitude") +
  scale_colour_gradient(low = "blue", high = "red") + 
  scale_size_binned_area()+
  facet_wrap(vars(year))

# duplicate station values in some years due to 24 hour experiments (ie 1991)

a<- df %>% group_by(year, station) %>% summarize(n=n())
# cap_larvae %>% group_by(year, station) %>% summarize(n=n()) %>% arrange(n) # duplicates in 1991 from multiple sampling of same station
# cap_larvae %>% group_by(year, station) %>% dplyr::filter(n()>1) %>% summarize(n=n())
a<-df %>% distinct(year, station, .keep_all = TRUE)%>% group_by(year, station) %>% summarize(n=n()) # looks good 1 value per station
df  %<>%  distinct(year, station, .keep_all = TRUE)
summary(df)
# df<-df %>% dplyr::filter(!is.na(Nm2))
df$year<-as.factor(df$year)
save(df, file = "./data/larvae/df.Rdata")
rm(a)

####
# Timing
a<-as.data.frame(table(df$year,df$doy))
a %>% ggplot(aes(Var1,Var2,size =Freq ))+geom_point()+scale_size_area(max_size = 6)

ggplot(df, aes(year, doy)) +
  geom_count() +
  scale_size_area(breaks = c(1,2,4,6,8,10,12))
ggplot(cap_larvae, aes(year, doy, colour = trip_no)) +
  geom_count() +
  scale_size_area(breaks = c(1,5,10,15),max_size = 15)
                  
############################  Model Run  ##############################################
load(file = "./data/larvae/df.Rdata")
library(tidyverse)
library(cplm)          # for mixed effect tweedie
library(performance)
library(report)
library(ggeffects)
library(tweedie)
library(lattice)
library(MASS)
require(pscl) # alternatively can use package ZIM for zero-inflated models
library(lmtest)
library(outliers)
library(statmod)

glimpse(df)
# Do I have reason to think the data is zero inflated?
# 1) There are many zeros but I think I can rule out sampling bias as we cover 0-50m and capelin don't go deeper as larvae

# df <- cap_larvae %>% dplyr::select(year, station, lat, lon, doy, Nm2,Nv,presence)

# M3 <- zeroinfl(Nb ~  year + doy | ## Predictor for the Poisson process
#                doy,      ## Predictor for the Bernoulli process;
#                dist = 'poisson',
#                # offset = log(volume),
#                # EM=TRUE,
#                data = cap_larvae)
# summary(M3)
# 
# # Dispersion statistic
# E2 <- resid(M3, type = "pearson")
# N  <- nrow(df)
# p  <- length(coef(M3))  
# sum(E2^2) / (N - p)
# # [1] 6.027588
# 
# ### Model 4 Zero Inflated Negative Binomial GLM ###
# 
# M4 <- zeroinfl(Nb ~ year +  doy| # count model probably varies with these variables
#                doy,                    # zero inflated bit probably influenced by only these
#                dist = 'negbin',
#                # offset = log(volume),
#                data = cap_larvae)
# M4b <- zeroinfl(Nm2 ~ year +  doy| # count model probably varies with these variables
#                  year,                    # zero inflated bit probably influenced by only these
#                dist = 'negbin',
#                # offset = log(volume),
#                data = df)
# # version where I round densities to neareset interger
# 
# # M4 diagnostics
# summary(M4)
# E2 <- resid(M4, type = "pearson")
# N  <- nrow(df)
# p  <- length(coef(M4))  
# sum(E2^2) / (N - p)
# # [1] 1.162405
# 
# # M4 diagnostics
# summary(M4b)
# E2 <- resid(M4b, type = "pearson")
# N  <- nrow(df)
# p  <- length(coef(M4b))  
# sum(E2^2) / (N - p)
# #[1] 0.9468142
# 
# 
# # compare both zero inflated models
# lrtest(M4, M4b) # higher log likelihood and significant value suggests the second model is better (ie zero inflated neg binom)
# Likelihood ratio test
# 
# Model 1: Nb ~ year + doy | doy
# Model 2: Nm2 ~ year + doy | year
# #Df  LogLik Df  Chisq Pr(>Chisq)    
# 1  40 -3081.9                         
# 2  74 -2275.5 34 1612.9  < 2.2e-16 ***
#   ---
# RESULT: zero inflated neg bin better AIC and less overdispersed

# compare with intercept only model
# InterceptModel <- update(M4, .~1)
# lrtest(M4, InterceptModel) # M4 better
# 
# ## Figures ## 
# 
# dat0 <- ggpredict(M4, type =  "fe.zi", terms = c("year"),typical = "median")
# plot(dat, colors = "hero", dot.size = 1.5)
# 

############################  Main Results  ##############################################
### Model 5 Tweedie GLM ### Advantage, takes care of zeros in continuous variables
# Produces a generalized linear model family object with any power variance function and any power link. 
# Includes the Gaussian, Poisson, gamma and inverse-Gaussian families as special cases.

# The most interesting Tweedie families occur for var.power between 1 and 2. 
# For these GLMs, the response distribution has mass at zero (i.e., it has exact zeros) 
# but is otherwise continuous on the positive real numbers (Smyth, 1996; Hasan et al, 2012). 
# These GLMs have been used to model rainfall for example.
# Many days there is no rain at all (exact zero) but, 
# if there is any rain, then the actual amount of rain is continuous and positive.

### Plot a Tweedie density with 1<p<2
yy <- seq(0,10,length=2000)
tweedie.plot( power=1.1, mu=1, phi=1, y=yy, lwd=2)
tweedie.plot( power=1.3, mu=1, phi=1, y=yy, add=TRUE, lwd=2, col="red")
tweedie.plot( power=1.4, mu=1, phi=1, y=yy, add=TRUE, lwd=2, col="blue")
legend("topright",lwd=c(2,2), col=c("black","red","blue"), pch=c(19,19),
       legend=c("p=3.0","p=2.0", "p=1.5") )


df %>% dplyr::filter(Nm2 < 20, Nm2 > 0) %>% ggplot(aes(Nm2)) + geom_density()
# log doy and reduce centre for whole time series
df %<>% mutate(ln_doy = log(doy),                          # log transform
               ln_doy_std = as.vector(scale(ln_doy)),      # log and centre reduce
               doy_std = as.vector(scale(doy)))            # just centre reduce

# Create second data set removing outliers for comparison
df2 = df$Nm2                                                 ## find outliers above 99% quantile           
UL = df2[df2 >0] %>% quantile(.975, na.rm = TRUE)
df_out<- df %>% dplyr::filter(Nm2 > UL)
df_out<- df %>% dplyr::filter(Nm2 < UL)


# power.link parameter estimated by maximum likelihood
estimate_p <- tweedie.profile(
  glm(Nm2 ~ year * doy,
      data = df_out),
  method = "mle",
  do.plot = TRUE,
  verbose = 2,
  do.ci = TRUE
) 


estimate_p$xi.max # [1] 1.677551 for full df; outliers removed = 1.653061
estimate_p$ci # 1.644129 1.644238 for full df; outliers removed = 1.628586 1.628596

## Larvae Model Null:

larvae_null_tweedie <- glm(Nm2 ~ 1,                          # Intercept only model aka Null model H0
                           family = tweedie(                 
                             var.power=1.677551,                    # Should be between 1 and 2. 2 ie gamma, needs to be strictly positive ie no zeros
                             link.power = 0),                 # link.power = log-link (based on a look at log(residuals) lm null model...)
                           data = df)

summary(larvae_null_tweedie)
AICtweedie(larvae_null_tweedie)
glance(larvae_null_tweedie)
r2(larvae_null_tweedie)

LA_null <- augment(larvae_fit_tweedie,type.predict = "response", se.fit = TRUE, newdata = df) %>% dplyr::select(year, .fitted,.se.fit)
## Larvae Model 1: 

# Goal: Create larval index adjusted for seasonal trend (year-continous) and day of standardised day of year (mean = 0, sd = 1)
# Ideally mixed/multilevel/hierarchical models would be prefered but they are not as well developed/documented for the tweedie distribution or positive continuous distributions in general (Gamma is also buggy)

larvae_fit_tweedie <- glm(
  Nm2 ~ year + doy,
  family = tweedie(var.power = 1.677551, link.power = 0),        
  data = df
)

# Summary statistics and model diagnostics
summary(larvae_fit_tweedie)     
tidy(larvae_fit_tweedie) 
glance(larvae_fit_tweedie)
check_model(larvae_fit_tweedie)                    # full suite of diagnostics usually employed by plot(model)
AICtweedie(larvae_fit_tweedie)     # Tweedie needs special specification for AIC
r2(larvae_fit_tweedie2)                             # Pseudo R2s
r2_zeroinflated(larvae_fit_tweedie)
rmse(larvae_fit_tweedie)                           # Root mean square error
model_performance(larvae_fit_tweedie)              # summary of above diagnostics
glance(larvae_outfit_tweedie)

compare_performance(larvae_null_tweedie, larvae_fit_tweedie)   # compare to null model

# Extract and merge model adjusted values and standard error on response scale with original data
LA_tweedie <- augment(larvae_fit_tweedie,type.predict = "response", se.fit = TRUE, newdata = df) 

## Larvae Model 2: (99.5% outliers removed)
larvae_outfit_tweedie <- glm(
  Nm2 ~ year + doy,
  family = tweedie(var.power = 1.653061, link.power = 0),        
  data = df_out
)
tidy(larvae_outfit_tweedie)           
check_model(larvae_outfit_tweedie)                    # full suite of diagnostics usually employed by plot(model)
AICtweedie(larvae_outfit_tweedie)                     # Tweedie needs special specification for AIC
r2(larvae_outfit_tweedie)                             # Pseudo R2
rmse(larvae_outfit_tweedie)                           # Root mean square error
model_performance(larvae_outfit_tweedie)              # summary of above diagnostics
glance(larvae_outfit_tweedie)
logLiktweedie(larvae_fit_tweedie)
logLiktweedie(larvae_null_tweedie)
logLiktweedie(larvae_outfit_tweedie)
# compare to null and full model
compare_performance(larvae_null_tweedie, larvae_fit_tweedie, larvae_outfit_tweedie)   

# Extract and merge model adjusted values and standard error on response scale with original data
LA_tweedie_out <- augment(larvae_outfit_tweedie,type.predict = "response", se.fit = TRUE, newdata = df_out) 


dat<-ggpredict(larvae_outfit_tweedie,terms = "year")
plot(dat)
datnull<-ggpredict(larvae_null_tweedie,terms = "year")
marginal_means <- as_tibble(data.frame(dat$x,dat$predicted,dat$conf.low,dat$conf.high))
df_summary<- df %>% group_by(year) %>% dplyr::summarise(mean_nm2 = mean(Nm2))

marginal_means %>% 
  ggplot(aes(dat.x,sqrt(dat.predicted)))+
  geom_point(size = 5, alpha = 0.5)+ 
  geom_errorbar(aes(ymin = sqrt(dat.predicted) - sqrt(dat.conf.low),
                    ymax = sqrt(dat.predicted) + sqrt(dat.conf.high))) +
  geom_point(data = df_summary,aes(year, sqrt(mean_nm2)), 
             alpha = 0.5, size = 3,colour = "red") + theme_minimal(base_size = 20) + labs(x = "Year", y = "sqrt ( Nm2 )")
  # stat_summary(fun.data = "mean_cl_boot", size = 1, alpha = 0.5, colour = "red")

# Save the extracted model inputs 
save(LA_tweedie,file="jags/input/LA_tweedie.Rdata")
save(LA_tweedie_out,file="jags/input/LA_tweedie_out.Rdata")


LA_tweedie %>% ggplot(aes(year,.fitted))+geom_point(colour = "black",size=4,alpha=0.5,position = position_jitter())
LA_tweedie %>% ggplot(aes(year,log(.fitted)))+geom_point(colour = "black",size=4,alpha=0.5,position = position_jitter())

# dat5 <- ggpredict(M5, type =  "fe", terms = c("year"), typical = "mean")
# plot(dat5, dot.size = 1.5, ci = FALSE) 
# dat5a <- ggpredict(M5_out, type =  "fe", terms = c("year"), typical = "mean")
# plot(dat5a, dot.size = 1.5, ci = FALSE) 


# M5a_out %>% group_by(year) %>% dplyr::summarise(mean_capelin_Nm2 = mean(Nm2)) %>% 
#   ggplot(aes(year, mean_capelin_Nm2)) +geom_point()
# M5a %>% group_by(year) %>% dplyr::summarise(mean_capelin_Nm2 = mean(Nm2)) %>% 
#   ggplot(aes(year, mean_capelin_Nm2)) +geom_point()
# 
# glm_tweedie <- M5a %>% group_by(year) %>% dplyr::summarise(mean_capelin_Nm2 = mean(Nm2))
# glm_tweedie_out <- M5a_out %>% group_by(year) %>% dplyr::summarise(mean_capelin_Nm2 = mean(Nm2))



### Model 10 Tweedie GLMM ### Advantage, takes care of zeros in continuous variables
library(cplm)

# use Stock and Spacing as main effects and Plant as random effect
f1 <- cpglmm(Nm2~  year + ln_doy_std +  (1|ln_doy_std) + (1|year) , 
             link="log", data = df)

f0 <- cpglmm(Nm2~  1 +  (1|ln_doy_std)  + (1|year) , 
             link="log", data = df)
summary(f1)
summary(f0)
# most of the methods defined in lme4 are directly applicable
coef(f1); fixef(f1); ranef(f1)  #coefficients
VarCorr(f1)  #variance components
hist(f1$resid)
predict(f1)
f1$deviance
summary(f1)
newdata <- data.frame(year = unique(df$year),
                         ln_doy_std = rep(0,34),
                      Nm2 = NA)
r1 <- resid(f1) / sqrt(f1$phi)
plot(r1 ~ fitted(f1), cex = 0.5)

# residual and qq plot
predict(f1,type = "tweedie" )
plotF(f1)
qqnorm(r1, cex = 0.5)
# 2. quantile residual plot to avoid overlapping
u <- tweedie::ptweedie(f1$y, f1$p, fitted(f1), f1$phi)
u[f1$y == 0] <- runif(sum(f1$y == 0), 0, u[f1$y == 0])
r2 <- qnorm(u)
plot(r2 ~ fitted(f1), cex = 0.5)
qqnorm(r2, cex = 0.5)




# https://stackoverflow.com/questions/52372340/creating-r-squared-function-for-cplm-package
# Fit null model without fixed effects (but including all random effects)
# parmodCPr <- cpglmm(Parasite ~ 1 + (1 | Population) + (1 | Container), data = DataAll)
# 
# # Fit alternative model including fixed and all random effects
# parmodCPf <- cpglmm(Parasite ~ Sex + Treatment + Habitat + (1 | Population) +
#                       (1 | Container), data = DataAll)

# Calculation of the variance in fitted values
VarF <- var(as.vector(model.matrix(f1) %*% fixef(f1)))

# getting the observation-level variance Null model
phiN <- f0@phi # the dispersion parameter
pN <- f0@p # the index parameter
mu <- exp(fixef(f0) + 0.5 * (VarCorr(f0)$year[1] + VarCorr(f0)$ln_doy_std[1]))
VarOdN <- phiN * mu^(pN - 2) # the delta method

# Full model
phiF <- f1@phi # the dispersion parameter
pF <- f1@p # the index parameter
VarOdF <- phiF * mu^(pF - 2) # the delta method

# R2[GLMM(m)] - marginal R2[GLMM]; using the delta method observation-level variance
R2glmmM <- VarF/(VarF + sum(as.numeric(VarCorr(f1))) + VarOdF)

# R2[GLMM(c)] - conditional R2[GLMM] for full model
R2glmmC <- (VarF + sum(as.numeric(VarCorr(f1))))/(VarF + sum(as.numeric(VarCorr(f1))) +VarOdF)


# simple glmer
library(lme4)




#########################
###########################
##### new as of 26/02/20

dfa = df$Nm2                                                 ## find outliers above 99% quantile           
UL = dfa[dfa >0] %>% quantile(.975, na.rm = TRUE)
df_out<- df %>% dplyr::filter(Nm2 < UL)
df_out$DOY <- as.factor(df_out$doy)
df$DOY <- as.factor(df$doy)
df 

fit1<-glm(data = df_out, log(Nm2+0.1) ~  year + DOY)
summary(fit1)
performance(fit1)
check_model(fit1)

dat<- ggpredict(fit1, terms = "year", condition = c(DOY = 170))
plot(dat)
dat <- as.data.frame(dat)
# get correction factor for back transforming fitted values from model
syx <- sigma(fit1)   # sigma from model 
cf <- exp((syx^2)/2)              # correction factor
dat$Nm2_pred <- exp(dat$predicted)*cf
dat_out <- dat
save(dat_out, file = "./Rdata/larvae_glm_out.Rdata")
ggplot(dat_out, aes(x, log(Nm2_pred)))+geom_point()

df2<- df_out %>% dplyr::filter(lon < - 63)
df2 %>% group_by(year) %>% dplyr::summarise(mean_Nm2 = mean(Nm2)) %>% ggplot(aes(year, log(mean_Nm2)))+geom_point()

fit2<- glm(data = df2, log(Nm2+0.1) ~  year + DOY)
summary(fit2)
performance(fit2)
check_model(fit2)

anova(fit2)
dat2<- ggpredict(fit2, terms = "year", condition = c(DOY = 170))
plot(dat2)
dat2 <- as.data.frame(dat2)

syx <- sigma(fit2)   # sigma from model 
cf <- exp((syx^2)/2)              # correction factor
dat2$Nm2_pred_2 <- exp(dat2$predicted)*cf
dat2_out<- dat2
save(dat2_out, file = "./Rdata/larvae_glm_2_out.Rdata")
dat2 %>%  ggplot(aes(x, log(Nm2_pred_2))) + geom_point()


# if model is glm use sigma(model)..also lm 
cpue_4R_mod2_se <- augment(cpue_4R_2,type.predict = "response",se.fit = TRUE) %>% dplyr::select(.se.fit)
cpue_4R_mod2 <- augment(cpue_4R_2,newdata = capelin_4R_2, type.predict = "response",se.fit = TRUE) %>% mutate(".se.fit" = cpue_4R_mod2_se$.se.fit)
cpue_4R_mod2 %<>% mutate(biased_fitted = exp(.fitted), fit_catch_t = biased_fitted*cf)


fit2 <- glm(data = df, Nm2 ~ year  + (1|doy))

summary(fit2)
performance(fit2)
check_model(fit2)
report(fit2) %>% to_table()
model_performance(fit1)
compare_performance(fit1,fit2,larvae_fit_tweedie,larvae_outfit_tweedie)

dat <- ggpredict(fit1, type =  "fe", terms = c("year"), typical = "mean")
dat2 <- ggpredict(fit2, type =  "fe", terms = c("year"), typical = "mean")
plot(dat2, dot.size = 1.5, ci = FALSE) 

glmer_rand.doy <- augment(fit1) %>% group_by(year) %>%  dplyr::summarise(mean_predicted_Nm2 = mean(.fitted)) 
glmer_rand.doy_out <- augment(fit2) %>% group_by(year) %>%  dplyr::summarise(mean_predicted_Nm2 = mean(.fitted)) 

library(readr)
results_capelin_larvae_kriging <- read_csv("results.capelin.larvae.kriging_new.csv")

capelin_GSL_larval_index <- data.frame(year = unique(dat1$x),
                                       NegBinom.zi_Nm2 = dat0$predicted,
                                       NegBinom.zi_SE = dat0$std.error,
                                       Tweedie_year_log.med.doy_Nm2 = dat1$predicted,
                                       Tweedie_year_log.med.doy_SE = dat1$std.error,
                                       Tweedie_year_Nm2 = dat2$predicted,
                                       Tweedie_year_SE = dat2$std.error,
                                       Tweedie_year_st.doy_Nm2 = dat3$predicted,
                                       Tweedie_year_st.doy_SE = dat3$std.error,
                                       Tweedie_year_st.doy.all_Nm2 = dat4$predicted,
                                       Tweedie_year_st.doy.all_SE = dat4$std.error,
                                       Krig.moy_Nm2 = results_capelin_larvae_kriging$moy.krig,
                                       Krig.var_Nm2 = results_capelin_larvae_kriging$var.krig) %>% as_tibble()

capelin_GSL_larval_index <- data.frame(year = unique(df$year),
                                       glmer_predicted_Nm2 = glmer_rand.doy$mean_predicted_Nm2,
                                       glmer_predicted_Nm2_no_out = glmer_rand.doy_out$mean_predicted_Nm2,
                                       tweedie_glm_predicted_Nm2 = glm_tweedie$mean_capelin_Nm2,
                                       tweedie_glm_predicted_Nm2_no_out = glm_tweedie_out$mean_capelin_Nm2,
                                       krigged_best_fit_Nm2 = results_capelin_larvae_kriging$moy.krig,
                                       Krig.var_Nm2 = results_capelin_larvae_kriging$var.krig) %>% as_tibble()

capelin_GSL_larval_index_new <- capelin_GSL_larval_index
# write.csv(capelin_GSL_larval_index, file = "./capelin_GSL_larval_index.csv")
write.csv(capelin_GSL_larval_index_new, file = "./capelin_GSL_larval_index_new.csv")
# capelin_GSL_larval_index <- read_csv("data/larvae/capelin_GSL_larval_index.csv") %>% select(-X1) %>% mutate(data_source = "SGSL_June_Survey", rscript = "cap_larvae.R")
# save to JAGS input
# this is the new one below
write.csv(capelin_GSL_larval_index, file = "jags/input/capelin_GSL_larval_index.csv")

capelin_GSL_larval_index_tall <- capelin_GSL_larval_index_new %>% pivot_longer(-year, names_to = "Model", values_to = "value")
p1<-capelin_GSL_larval_index_tall %>% dplyr::filter(
  Model %in% c(
    "glmer_predicted_Nm2",
    "glmer_predicted_Nm2_no_out",
    "tweedie_glm_predicted_Nm2",
    "tweedie_glm_predicted_Nm2_no_out",
    "krigged_best_fit_Nm2"
  )
) %>%   ggplot(aes(year, value, colour = Model), size = 4, alpha = 0.7) + geom_point(size = 4) +geom_line(aes(group = Model)) + scale_colour_viridis_d(begin = 0, end = 0.90) + scale_x_discrete(breaks= seq(1980,2020,10)) +
  facet_wrap(vars(Model),scales = "free_y") + theme_minimal(base_size = 18)
ggsave(p1, filename = "./img/capelin_larval_index.png")


############################################ figures
larvae <- read_csv("jags/input/capelin_GSL_larval_index.csv") %>% dplyr::select(-X1)
glimpse(cap_larvae) # raw data
summary(cap_larvae)
glimpse(larvae)    # predicted data
sGSL_area <- 6.945e+10
GSL_area <- 250000000

cap_larvae_summary <- cap_larvae %>% 
  dplyr::filter(trip_no == 1) %>%                   # trip one only
  group_by(year) %>%                               # aggregate by year
  dplyr::summarise(model = "raw_data", mean_NM2 = mean(Nm2,na.rm = T),   # mean yearly density
                   sgsl_total_nm2 = mean_NM2*sGSL_area, # Number of larvae if raised to area of sGSL
                   N_sGSL_billions = round(sum(sgsl_total_nm2)/10e9,2) # in billions for plotting
                   )

# # predicted larvae Nm2
# larvae_summary <- larvae %>% dplyr::select(year, glmer_predicted_Nm2, tweedie_glm_predicted_Nm2, krigged_best_fit_Nm2) %>% 
#   pivot_longer(names_to = "model",values_to = "Nm2",-year) %>% 
#   group_by(year, model) %>% 
#   mutate(mean_NM2 = mean(Nm2,na.rm=T),                                # mean yearly density
#          sgsl_total_nm2 = mean_NM2*sGSL_area,                         # Number of larvae if raised to area of sGSL
#          N_sGSL_billions = round(sum(sgsl_total_nm2)/10e9,2))         # in billions for plotting)
                                    

# plot of log(density) of capelin over sGSL per year raw data. 
# Note that this is a view of data if no zeros as NAs are produced with zeros so only a view of what spread sans zero is
glimpse(cap_larvae)
cap_larvae %>% 
  ggplot(aes(year, log(Nm2)))+
  geom_point() + 
  stat_summary(fun.data= mean_cl_boot,colour = "red",size = 1.5) + 
  labs(x = "Year", y = "log(larval density Nm2)") +
  theme_minimal(base_size = 20)

# plot of  log number of capelin over sGSL in billions per year using Nm2 x area of sGSL raw data
cap_larvae_summary %>% 
  ggplot(aes(year, log(N_sGSL_billions),label =N_sGSL_billions))+
  geom_point(size = 4) + geom_line(size = 1, alpha = 0.3) + 
  labs(title = "Mean larval density (Nm2/Year) * survey Area (6.945e+10 m2)",  x = "Year", y = "log(N of capelin larvae in billions)") +
  geom_text(vjust = 0, nudge_y = 0.5, colour = "red",size = 4) +
  theme_minimal(base_size = 20)

# compare predicted and raw means
c_larvae <- bind_rows(cap_larvae_summary,larvae_summary) 

c_larvae %>%
  ggplot(aes(year, log(N_sGSL_billions), colour = model, label = N_sGSL_billions))+
  geom_point(size = 4,alpha = 0.6) + 
  geom_line(size = 1) + 
  labs(x = "Year", y = "log(N of capelin larvae in billions)") +
  # geom_text(vjust = 0, nudge_y = 0.5, colour = "red") +
  theme_minimal(base_size = 20) +
  scale_color_viridis_d(end=0.8,option = "A")

# plot fitted raised to sgsl
LA_tweedie_summary <- LA_tweedie %>% 
  group_by(year) %>% 
  dplyr::summarise(median_fitted = median(.fitted),
                   N_billions = median_fitted*sGSL_area/10e9)

LA_tweedie_summary %>% 
  ggplot(aes(year, log(N_billions),label=round(N_billions),4))+
  geom_point(size = 8, colour = "red") + 
  geom_line(alpha = 0.5) +
  geom_text(vjust = 0, nudge_y = 0.1, colour = "black",size = 6) +
  labs(title = "Adjusted median Nm2 * Survey Area (6.945e+10 m2)", x = "Year", y = "log(N larvae - billions)") +
  theme_minimal(base_size = 24)

# Observed and fitted
latall<-LA_tweedie %>% dplyr::select(year, Nm2, .fitted, .se.fit) %>% pivot_longer(cols = Nm2:.se.fit, names_to = "data",values_to = "value")
  latall %>% dplyr
    ggplot(aes(year, sqrt(value),))+ 
  stat_summary(fun.data = "mean_cl_boot", size = 1, alpha = 0.7, colour = "black")+
  geom_point(aes(year, sqrt(.fitted)), alpha = 0,position = position_dodge(1))+
  stat_summary(fun.data = "mean_cl_boot", size = 1, alpha = 0.7, colour = "red")
  


lasum<-LA_tweedie %>% group_by(year) %>% dplyr::summarise(mean_Nm2 = mean(Nm2,na.rm=TRUE))

d <- ggplot(LA_tweedie, aes(year, sqrt(.fitted))) 
d + stat_summary(fun.data = "mean_cl_boot", size = 1, alpha = 0.7) +
  geom_point(data = lasum, aes(year, sqrt(Nm2)),alpha = 0) +
  
 
  

# larval density vs cpue
fig_1 <- ggplot(df %>% dplyr::filter(year>1984),aes(lag(ln_Nm2,2), log(cpue_4R_mod4))) +
  # stat_summary(fun.data= mean_cl_normal) +
  # geom_smooth(method='lm',se = F) +
  geom_point(size = 5) +
  labs(x = "Larves t-2", y = "ln Indice de performance t") + theme_bw(base_size = 18) #+ ylim(1,4)
# condition vs cpue
fig_2 <- ggplot(df %>% dplyr::filter(year>1984),aes(lag(mean_K,1), log(cpue_4R_mod4))) +
  # stat_summary(fun.data= mean_cl_normal) + 
  # geom_smooth(method='lm',se = F) + 
  geom_point(size = 5) + labs(x = "Condition t-1", y = "ln_IP") + theme_bw(base_size = 18) #+ ylim(1,4)
# plot(df$mean_K_std_lag1,log(df$cpue_4R_mod4))
# last ice vs cpue
fig_3 <- ggplot(df,aes(last_ice_std, log(cpue_4R_mod4))) +
  # stat_summary(fun.data= mean_cl_normal) + 
  # geom_smooth(method='lm',se = F) +
  geom_point(size = 5) + labs(x = "Glace t", y = "ln_IP") + theme_bw(base_size = 18) #+ ylim(1,4)

#cpue vs time
fig_4 <- ggplot(df,aes(year, log(cpue_4R_mod1))) +
  geom_line()+ geom_point(size = 5) + labs(x = "", y = "ln_4R_IP") +  theme_classic(base_size = 24) + xlim(1985,2020)
#cpue vs time
fig_5 <- ggplot(df,aes(year, log(cpue_3KL_mod5))) +
  geom_line()+ geom_point(size = 5) + labs(x = "", y = "ln_3KL_IP") +  theme_classic(base_size = 24) + xlim(1985,2020)
#biomass vs time
fig_6 <- ggplot(df,aes(year, log(biomass_med))) +
  geom_line()+ geom_point(size = 5) + labs(x = "Année", y = "ln_bio_3KL") + theme_classic(base_size = 24) + xlim(1985,2020)
#abundance vs time
fig_7 <-ggplot(capelin_2j3kl_biomass,aes(year, log(abundance_med))) +
  geom_line()+ geom_point(size = 5) + labs(x = "Année", y = "ln_abondance_3KL") + theme_classic(base_size = 24)

plot1 <- cowplot::plot_grid(fig_1, fig_2, fig_3, labels = c("A", "B", "C"), ncol=1) + theme_minimal(base_size = 20)

plot2 <- cowplot::plot_grid(fig_4, fig_5, fig_6, labels = c("D", "E", "F"), ncol=1) 

fig_1
