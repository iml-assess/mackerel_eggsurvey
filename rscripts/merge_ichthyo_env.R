# merge plankton data with environmental data
# this is a huge pain and still isn't 100 perfect. dates, coords, station names, etc. don't match between datasets and sometimes there are multiple ctd casts for stations and sometimes there are ctd casts that we don't have ichthyo data for etc etc.
# This was my latest attempt to merge the two datasets through a cascading number of variables the sets had in common and filling the holes as best I could

load("./Rdata/survey_env.Rdata")
load("./Rdata/df_egg_index_2022_v4.Rdata")

glimpse(df_egg_index)
glimpse(survey_env)

# change var names for conveniance
survey_env %<>% transmute(year,month, day, doy, station = as.factor(station), temperature_CL, temperature_FG, salinity, sigma_t_kg_m3)

# standardize coordinates across dfs
coords <- df_egg_index %>% group_by(station) %>% dplyr::summarise(lon = mean(lon,na.rm=T), lat = mean(lat, na.rm= T))
df_egg_index %<>% dplyr::select(-lon,-lat) 

# merge coordinates back into dataframe so that they match
df_egg_index <- left_join(df_egg_index, coords)
survey_env <- left_join(survey_env, coords)

# figures to validate merges
survey_env %>% ggplot(aes(lon,lat, colour = temperature_CL))+geom_point()+scale_colour_viridis_c(option = "C")+facet_wrap(vars(year))
survey_env %>% ggplot(aes(lon,lat, colour = temperature_FG))+geom_point()+scale_colour_viridis_c(option = "C")+facet_wrap(vars(year))
df_egg_index %>% dplyr::filter(stage %in% c("n_15","n2","n3","n4","nl")) %>%  ggplot(aes(year,n_m2, colour = stage))+stat_summary(fun.data = "mean_cl_boot")
df_egg_index %>% dplyr::filter(stage %in% c("n_15","n2","n3","n4","nl")) %>%  ggplot(aes(year,n_m2, colour = stage))+geom_smooth()

# filter out raw counts  for simplicity
df_egg_index %<>% dplyr::filter(stage %in% c("n_15","n2","n3","n4","nl")) # n15 is stages 1 and 5; 2 = 2 etc. l = larvae
summary(df_egg_index)

# merge. as file structures and dates are messy do this progressively. first merge produces NA's   :4801 for t_CL; second = NA's   :4451; third = NA's   :3971. Not great. Trying something else below
df <- left_join(df_egg_index, survey_env, by = c("year", "station", "month", "day"))
df2 <- left_join(df, survey_env, by = c("year", "station", "month"))
df3 <- left_join(df, survey_env, by = c("year", "station"))
summary(df)
summary(df2)
summary(df3)

# survey env sometimes has multiple casts for same station on same date.... take mean
dat <-survey_env %>% group_by(year,month,day, station) %>% dplyr::summarise(temperature = mean(temperature_CL,na.rm=T), temperature_fg = mean(temperature_FG,na.rm=T), salinity = mean(salinity,na.rm=T), sigmaT = mean(sigma_t_kg_m3,na.rm=T) )
dat2 <-survey_env %>% group_by(year,month,day, station) %>% dplyr::summarise(n=n(),nt = n_distinct(temperature_CL), ntt = n_distinct(temperature_FG)) %>% arrange(-n)

# try simpler then join it back to main set
meta <- df_egg_index %>% expand(nesting(year, month, day, station)) # gets all unique values for these variables vs expand which gives all possible combinations
meta2 <- df_egg_index %>% expand(nesting(year, month, station))
meta3 <- df_egg_index %>% expand(nesting(year, station))
meta <- left_join(meta,dat)
meta2 <- left_join(meta2,dat)
meta3 <- left_join(meta3,dat)
summary(meta)
m <-left_join(meta,meta2,by = c("year", "month","station"))
m2 <-left_join(m,meta3,by = c("year", "station"))
summary(meta)
summary(m)
summary(m2)
# results are that ctd NAs go from 953->876->779... better


# put it all together
env <- m2 %>% transmute(year, month = month.x, day = day.x, station, 
                        tempCL = ifelse(is.na(temperature.x),temperature.y,temperature.x),   # CL = Caroline Lafleur ctd data
                        tempCL = ifelse(is.na(tempCL),temperature,tempCL),
                        tempFG = ifelse(is.na(temperature_fg.x),temperature_fg.y,temperature_fg.x), # FG = Francois Gregoire from 2014 resdocs
                        tempFG = ifelse(is.na(tempFG),temperature_fg,tempFG),
                        sal = ifelse(is.na(salinity.x),salinity.y,salinity.x),
                        sal = ifelse(is.na(sal),salinity,sal),
                        sigma_t = ifelse(is.na(sigmaT.x),sigmaT.y,sigmaT.x),
                        sigma_t = ifelse(is.na(sigma_t),sigmaT,sigma_t),
                        temperature = ifelse(is.na(tempCL),tempFG,tempCL)
) %>% 
  group_by(year, month, day, station) %>% 
  dplyr::summarise(temperature = mean(temperature, na.rm=T),
                   salinity = mean(sal,na.rm=T),
                   sigma_t = mean(sigma_t,na.rm=T))

summary(env) # NA's   :156 for temperature.. nice
d2022 <- survey_env %>% dplyr::filter(year==2022) %>% 
  transmute(year, month, day, station, salinity, sigma_t = sigma_t_kg_m3, 
            temperature = temperature_CL)
env<-bind_rows(env,d2022)
glimpse(df_egg_index)
egg <- df_egg_index %>% dplyr::filter(stage == "n_15") %>% 
  dplyr::select(year, month, day, trajet, consec, bongo, station, station_depth, sample_depth, volume, set_duration, hour_start, minute_start, presence, n, n_m3, n_m2)

# fill the holes
egg<-full_join(egg,env)
egg<-left_join(egg,coords)

egg %<>% mutate(trajet = ifelse(is.na(trajet),"T1",trajet),
                bongo = ifelse(is.na(bongo), "B", "T"))
egg %<>% group_by(station) %>% mutate(station_depth = ifelse(is.na(station_depth),median(station_depth,na.rm=T),station_depth)) %>% ungroup()

egg %<>% 
  mutate(stratum = as.factor(ifelse(station %in% 
                                      c("8.7", "7.7", "5.7","4.9","3.9","3.8","2.6","1.5","2.5","1.4","2.4","1.3","3.5","3.4",'3.3',"3.2",'3.1',"4.2","4.1","2.3","2.2",'2.1',"1.2","1.1"), "1", 
                                    ifelse(station %in% 
                                             c("10.1","9.4","9.3",'8.2','7.2','7.1',"6.3","6.2","6.1","5.2","5.1",'4.3',"4.4","4.5",'4.6',"3.6","3.7","4.8","5.6",'6.7',"6.8","7.5","7.6","6.6"), '2', "3")))) 

egg %<>% mutate(stratum_area = ifelse(stratum == "1",29.61e+9,
                                      ifelse(stratum == "2", 21.91e+9,
                                             ifelse(stratum == "3", 93e+9, NA))),
                total_area = 29.61e+9 + 21.91e+9 + 93e+9)

# if temp na then use annual stratified mean for now
egg %<>% group_by(year, stratum) %>% mutate(temperature = ifelse(is.na(temperature),mean(temperature,na.rm=T),temperature)) %>% ungroup()


# this is the latest version of the merged complete dataset - Andrew 29/09/2022   ... 2022 data about 3/4 done from what I saw with Mélanie yesterday
save(egg, file = "./Rdata/egg_index_stage15_1979-2022.Rdata")
load("./Rdata/egg_index_stage15_1979-2022.Rdata")
