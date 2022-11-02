# taken largely from mackerel_ichthyo.R
# this is the script to fit models to GSI to adjust egg data

# Load bio data
load("~/My_Stocks/Mackerel/Assessments/2022/mackerel_2022/data/rdata/2022/mackerel_bio_lf_2022.Rdata")
bio<-mackerel_bio
rm(mackerel_bio)
glimpse(bio)

# load survey data 
load("C:/Users/SmithAND/Desktop/mackerel_egg_index_2022/rdata/egg_index_stage15_1979-2022.Rdata")

# load krigged dep from daily_egg_production_2022.R
load("C:/Users/SmithAND/Desktop/mackerel_egg_index_2022/rdata/new_krigged_dep_2022.Rdata")

# load published data
egg_survey_stats <- read_csv("data/ichthyoplankton/metadata_and_results/egg_survey_stats.csv")
egg_survey_stats %<>% dplyr::select(year,median_survey_date,mean_egg_density_nm2,mean_daily_egg_production,krigged_mean_dep)

# subset data and filter out weird values
bio %<>% dplyr::select(year,month,day,doy,division,subdivision,gear_name,sex,length,lcat5,mass,mass_gonad,maturity_stage,age) %>% 
  dplyr::filter(!is.na(mass),!is.na(mass_gonad),!is.na(doy),
                mass>0,mass>mass_gonad, maturity_stage!=0, division == "4T") # consider expanding beyond 4T if limited bio data in a given year

# calculate gsi
bio %<>% mutate(gsi=(mass_gonad/mass)*100)

# check distributions
bio %>% ggplot(aes(maturity_stage,gsi))+geom_boxplot()   
bio %>% ggplot(aes(doy,gsi, colour= maturity_stage))+geom_point()   

# remove outliers by stage
source("./functions/remove_outliers.R")
df <- bio %>%
  group_by(maturity_stage) %>%
  dplyr::mutate(gsi = remove_outliers(gsi))

# by group
df2<- bio %>% 
  group_by(doy) %>% 
  dplyr::mutate(gsi = remove_outliers(gsi))

# by both - keeping this one
df3 <- df %>% group_by(doy) %>% 
  dplyr::mutate(gsi = remove_outliers(gsi))
  
# check if that cleans data up
df %>% ggplot(aes(maturity_stage,gsi))+geom_boxplot()         
df2%>% ggplot(aes(maturity_stage,gsi))+geom_boxplot()   
df3%>% ggplot(aes(maturity_stage,gsi))+geom_boxplot()   
  

df %>% dplyr::filter(month%in%c(4,5,6,7,8,9), maturity_stage %in% c(4,5,6,7,8)) %>% 
  ggplot(aes(doy,gsi,colour=maturity_stage))+stat_summary(fun.data = "mean_cl_boot") 

df2 %>% dplyr::filter(month%in%c(4,5,6,7,8,9), maturity_stage %in% c(4,5,6,7,8)) %>% 
  ggplot(aes(doy,gsi,colour=maturity_stage))+stat_summary(fun.data = "mean_cl_boot") 

df3 %>% dplyr::filter(month%in%c(4,5,6,7,8,9), maturity_stage %in% c(4,5,6,7,8)) %>% 
  ggplot(aes(doy,gsi,colour=maturity_stage))+stat_summary(fun.data = "mean_cl_boot") 


# subsetting to commercial samples (those before may are from Maritimes RV survey and are 1) not yet ripe, 2) will mess up logistic model fitting
df %<>% dplyr::filter(!is.na(gsi), doy>135, maturity_stage>3)
df2 %<>% dplyr::filter(!is.na(gsi), doy>135, maturity_stage>3)
df3 %<>% dplyr::filter(!is.na(gsi), doy>135, maturity_stage>3)

library(devtools)
install_github("onofriandreapg/aomisc")

# loading package
library(aomisc) # for nls self starting init values

# FG used different formulation of a logistic function but that provides similar results. with aomisc there are self starting functions for 3 and 4 parameter logistic functions
# fit global model - normally you do this by year but this could be a good alternative when gsi data limited in a given year
model <- nls(gsi ~ NLS.L4(doy, b, c, d, e), data = df3) # 4 parameter logistic model with self starting values for params

df4<-df3 %>% dplyr::filter(maturity_stage%in% c(4,5,6)) # subset if only use stages 4-6
model2 <- nls(gsi ~ NLS.L4(doy, b, c, d, e), data = df4) # 4 parameter logistic model with self starting values for params

# using formulation that fg published
logistic.out <- nls(gsi ~ y0 + (a/(1+(doy/x0)^b)),       # modified 4 paramater logistic model (bassically the penalty on the height of the curve is removed - as per FG)
                    data = df3, start = list(y0 = 1, x0 = 171, a = 10, b = 30))

df2021<-df3 %>% dplyr::filter(year==2021)
model2021 <- nls(gsi ~ NLS.L4(doy, b, c, d, e), data = df2021)
# get model coeffs for predictions and plotting  -> be careful with naming convention
library(broom)
co <- tidy(model) %>% dplyr::select(term, estimate) 
co2 <- tidy(logistic.out) %>% dplyr::select(term, estimate) 
co3 <- tidy(model2) %>% dplyr::select(term, estimate) 
co2021 <- tidy(model2021) %>% dplyr::select(term, estimate) 

# plot the fits to compare with raw data
ggplot(df3,aes(x=doy,y=gsi, colour = factor(maturity_stage))) +
  stat_summary(fun.data = "mean_cl_boot") +
  stat_function(fun=function(x) co$estimate[2] + (co$estimate[3]-co$estimate[2])/(1+exp(-co$estimate[1]*(x-co$estimate[4]))), col='red', size=1.5, alpha = 0.5) +
  stat_function(fun=function(x) co3$estimate[2] + (co3$estimate[3]-co3$estimate[2])/(1+exp(-co3$estimate[1]*(x-co3$estimate[4]))), col='black', size=1.5, alpha = 0.5) +
  stat_function(fun=function(x) co2$estimate[1] + (co2$estimate[3]/(1+(x/co2$estimate[2])^co2$estimate[4])), col='blue', size=1.5, alpha = 0.5) +
  scale_color_viridis_d(direction = -1) +
  scale_x_continuous(breaks = seq(150,300,10)) +
  labs(x = "Jour / Day", y = "IGS / GSI", colour = "Stade de Maturité / Maturity Stage" ) +
  theme_bw() +
  theme(legend.position="top")
# choice of model formulation changes little. Keeping only stages 4-6 increases max value by a bit. See below to see that peak spawning day only changes from june 20 to june 21 by removing stages 7 and 8

# plot global model vs 2021
ggplot(df2021,aes(x=doy,y=gsi, colour = factor(maturity_stage))) +
  stat_summary(fun.data = "mean_cl_boot") +
  stat_function(fun=function(x) co$estimate[2] + (co$estimate[3]-co$estimate[2])/(1+exp(-co$estimate[1]*(x-co$estimate[4]))), col='black', size=1.5, alpha = 0.5) +
  stat_function(fun=function(x) co2021$estimate[2] + (co2021$estimate[3]-co2021$estimate[2])/(1+exp(-co2021$estimate[1]*(x-co2021$estimate[4]))), col='red', size=1.5, alpha = 0.5) +
  annotate("text", x = 250, y = 10, colour = "black", label = "global model") +
  annotate("text", x = 250, y = 8, colour = "red", label = "2021 model") +
  scale_color_viridis_d(direction = -1) +
  scale_x_continuous(breaks = seq(150,300,10)) +
  labs(x = "Jour / Day", y = "IGS / GSI", colour = "Stade de Maturité / Maturity Stage" ) +
  theme_bw() +
  theme(legend.position="top")

# plot fit of 2021 to data
ggplot(df2021,aes(x=doy,y=gsi, colour = factor(maturity_stage))) +
  stat_summary(fun.data = "mean_cl_boot") +
  stat_function(fun=function(x) co2021$estimate[2] + (co2021$estimate[3]-co2021$estimate[2])/(1+exp(-co2021$estimate[1]*(x-co2021$estimate[4]))), col='red', size=1.5, alpha = 0.5) +
  scale_color_viridis_d(direction = -1) +
  scale_x_continuous(breaks = seq(150,300,10)) +
  labs(x = "Jour / Day", y = "IGS / GSI", colour = "Stade de Maturité / Maturity Stage" ) +
  theme_bw() +
  theme(legend.position="top")

# Get predictions of proportion of eggs spawned on given day 
d1 <- data_frame(day=seq(100.5,350.5,1))
d <- data_frame(day=seq(100.5,350.5,1))
d1$pred <- co2$estimate[1] + (co2$estimate[3]/(1+(d1$day/co2$estimate[2])^co2$estimate[4]))   # original 
d$pred <- co$estimate[2] + (co$estimate[3]-co$estimate[2])/(1+exp(-co$estimate[1]*(d$day-co$estimate[4]))) # using the 4 param logistic model formulation
d2 <- data.frame(day=seq(101,350,1))
d3 <- data.frame(day=seq(101,350,1))
# this will give you the ascending slope of the curve
for(j in 1:(nrow(d1)-1)){d2[j,'slope']=d1[j+1,'pred']-d1[j,'pred']}
d2$prob <- d2$slope/sum(d2$slope,na.rm=T)
for(j in 1:(nrow(d)-1)){d3[j,'slope']=d[j+1,'pred']-d[j,'pred']}
d3$prob <- d3$slope/sum(d3$slope,na.rm=T)

#2021
d2021 <- data_frame(day=seq(100.5,350.5,1))
d2021$pred <- co2021$estimate[2] + (co2021$estimate[3]-co2021$estimate[2])/(1+exp(-co2021$estimate[1]*(d$day-co2021$estimate[4]))) # using the 4 param logistic model formulation
D2 <- data.frame(day=seq(101,350,1))
for(j in 1:(nrow(d2021)-1)){D2[j,'slope']=d2021[j+1,'pred']-d2021[j,'pred']}
D2$prob <- D2$slope/sum(D2$slope,na.rm=T)

#beginning and end of spawning - global model
upper <- mean(d2$prob) + (sd(d2$prob)/sqrt(length(d2$prob))) 
lower <- mean(d2$prob) - (sd(d2$prob)/sqrt(length(d2$prob)))
start <- d2 %>% dplyr::filter(day<160) %>% slice(which.min(abs(upper - prob)))
end <- d2 %>% dplyr::filter(day>160) %>% slice(which.min(abs(lower - prob)))

# 2021
upper2021 <- mean(D2$prob) + (sd(D2$prob)/sqrt(length(D2$prob))) 
lower2021 <- mean(D2$prob) - (sd(D2$prob)/sqrt(length(D2$prob)))
start2021 <- D2 %>% dplyr::filter(day<160) %>% slice(which.min(abs(upper2021 - prob)))
end2021 <- D2 %>% dplyr::filter(day>160) %>% slice(which.min(abs(lower2021 - prob)))

# peak spawning 

# spawning.prob[spawning.prob$year==2018,'peak.day'] = d2[d2$prob==max(d2$prob),"day"]
# spawning.prob[spawning.prob$year==2018,'prob'] = d2[d2$day==date.med$doy,"prob"]

# global
d2[d2$prob==max(d2$prob),"day"] # peak day, all years, stages 4-8 = 172
d3[d3$prob==max(d3$prob),"day"] # peak day, all years, only stages 4-6 = 173

#2021
D2[D2$prob==max(D2$prob),"day"] # peak day, = 159??? very early... validate this!


# mission dates
load("~/Data/Ichthyoplankton/Rdata/egg_index_stage15_1979-2022.Rdata")

egg %<>% dplyr::filter(!is.na(month),!is.na(day)) %>% 
  mutate(date = lubridate::ymd(paste(year, month, day)),
doy = as.numeric(format(date, "%j")))
eggdf<-egg %>% group_by(year,trajet,doy,station) %>% 
  dplyr::summarise(n = n())

# plot
eggdf %>% ggplot(aes(year,doy,colour=trajet))+geom_point()+
  geom_hline(aes(yintercept=147))+geom_hline(aes(yintercept=203)) # start and end of spawning
eggdf$day<-eggdf$doy

eggdf2 <- eggdf %>% group_by(day) %>% dplyr::summarise(n=n())
summary(eggdf)
d2<-left_join(d2,eggdf2)
  
  # plot of global spawning season all years with number of samples per day
  ggplot(d2,aes(x=day,y=prob))+geom_line(size = 1.5)+
  scale_x_continuous(breaks=seq(0,305,10),limits=c(100,250)) +
    theme_bw() + 
  ggtitle("Pic de ponte tous les années / Peak spawning all years") +  
  annotate("text", x = d2[d2$prob==max(d2$prob),"day"], y =max(d2$prob)+0.003, label = "06/20") + # peak
  annotate("text", x = start$day-10, y = 0.005, label = "05/26") +             # start
  annotate("text", x = end$day+8, y = 0.005, label = "07/21") +                 # end
  annotate("point", x = start$day, y = start$prob, colour = "blue", size = 5) +
  annotate("point", x = end$day, y = end$prob, colour = "blue", size = 5) +
  annotate("point", x = d2[d2$prob==max(d2$prob),"day"], y = max(d2$prob), colour = "blue", size = 5) +
  annotate("segment", x = start$day, xend = end$day, y = 0, yend = 0,  colour = "blue", size = 3) +
    # annotate("text", x = d2$day, y = -0.001, colour = "black", label = d2$n, size = 3, position = position_jitterdodge(jitter.width = 0,jitter.height = 0.002)) +
  geom_text(aes(day,-0.003,label=n),size=4,angle=90) +
    labs(y = "Proportion", x = "Jour / Day") +
  ylim(c(-0.005,0.036)) 
    
  
  # plot 2021 vs global model
  ggplot(d2,aes(x=day,y=prob))+geom_line(size = 1.5, colour = "black")+
    geom_line(aes(D2$day,D2$prob),colour="red")
  ##### hmmm ... it would be the earliest season ever...if correct. But this is me just quickly calculating it... also global model contains 2021 so consider running global without 2021
  
  table(bio$year,bio$maturity_stage) # sample size seems ok... 
  bio %<>% mutate(y = ifelse(year==2021,1,0))
    
  bio %>% ggplot(aes(doy,gsi,colour=factor(y)))+stat_summary(fun.data = "mean_cl_boot")+xlim(135,250)
  
  
  
  # from here look at script from 2020 mackerel_ichthyo.R to adjust the estimates of annual dep to the annual spawning season