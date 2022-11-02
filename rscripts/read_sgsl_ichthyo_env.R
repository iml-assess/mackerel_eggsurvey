# read sgsl ichthyo environmental data
# delivered in .dat format by caroline lafleur
# copy pasted and corrected by ADS in sublime and saved as plain text
sgsl_ichthyo_env <- read_csv("data/sgsl_ichthyo_env.csv") 
glimpse(sgsl_ichthyo_env)
# FG published
sgsl_ichthyo_temp <- read_csv("data/sgsl_ichthyo_temp.csv")
glimpse(sgsl_ichthyo_temp)

# main set
load("~/Data/Ichthyoplankton/Rdata/df_egg_index_2022_v3.Rdata")
df<-df_egg_index %>% dplyr::filter(stage == "Scomber scombrus (oeuf stade 1)") %>% 
  dplyr::select(date,year,month,day,trajet,station,lat,lon) %>% mutate(hour=hour(date),minute=minute(date))
df %<>% dplyr::filter(year>1986, year<2020) %>% mutate(set="main")

# wrangle
sgsl_ichthyo_env$day <- as.numeric(str_sub(sgsl_ichthyo_env$date_utc,1,2))
sgsl_ichthyo_env$month <- str_sub(sgsl_ichthyo_env$date_utc,4,6)
sgsl_ichthyo_env %<>% mutate(month = as.numeric(ifelse(month =="Jun",as.numeric(6),as.numeric(7))))
sgsl_ichthyo_env$year <- str_sub(sgsl_ichthyo_env$date_utc,8,11)
sgsl_ichthyo_env$hour <- as.numeric(str_sub(sgsl_ichthyo_env$time,1,2))
sgsl_ichthyo_env$minute <- as.numeric(str_sub(sgsl_ichthyo_env$time,4,5))
sgsl_ichthyo_env %<>% mutate(date = as_date(paste(year,month,day,sep = "/"), format = "%y/%m/%d"), year = year(date))
sgsl_ichthyo_env$lat <- as.numeric(sgsl_ichthyo_env$lat)
sgsl_ichthyo_env$lon <- as.numeric(sgsl_ichthyo_env$lon)

# filter out extra stations and merge synonymous station names...
sgsl_ichthyo_env %<>% mutate(Station  = trimws(str_remove_all(station, "[abcdefghijklmnopqrstuvwyzABCDEFGHIJKLMNOPQRSTUBWXYZ()-]")))
unique(sgsl_ichthyo_env$Station) %>% sort()
sgsl_ichthyo_env %<>% dplyr::filter(!Station %in% c("5","6","7","8","9","1301","1302","1303","1304","2","002","16.4","16.3","15.5","15.4","3.11", "3.12", "3.13", "3.14", "3.15","3.16","3.17","3.18","3.19","3.21","3.22", "3.23","3.24","3.25",'3.26',"3.27", "3.28", "3.29","xx"))
unique(sgsl_ichthyo_env$Station) %>% sort()

# to help see positions
a<-sgsl_ichthyo_env %>% group_by(Station) %>% dplyr::summarise(lat=mean(lat),lon=mean(lon))
a %>% ggplot(aes(lon,lat,label=Station))+geom_text()

sgsl_ichthyo_env %<>% 
  mutate(Station = fct_collapse(sgsl_ichthyo_env$Station, 
                                "12.1" = c("1201","121","1251"), 
                                "11.1" = c("1101","111","1151"),
                                "10.1" = c("1001","1051","101","110.1"),
                                "9.5" = c("905","0905","95"),
                                "9.4" = c("904","0904","94"),
                                "9.3" = c("903","0903","93"),
                                "9.2" = c("902","0902","92"),
                                "9.1" = c("901",'0901',"91"),
                                "8.7" = c("807","0807","87"),
                                "8.6" = c("806","0806","86"),
                                "8.5" = c("805","0805","85"),
                                "8.4" = c("804","8.4 (TIDM10)","TIDM10 (8.4)","0804","10 8.4","84","8.4 10"),
                                "8.3" = c("803","0803","83"),
                                "8.2" = c("802","0802","82"),
                                "8.1" = c("801",'0801',"81"),
                                "7.7" = c("707","0707","77"),
                                "7.6" = c("706","0706","76"),
                                "7.5" = c("705","0705","75"),
                                "7.4" = c("704","0704","74"),
                                "7.3" = c("703",'0703',"73"),
                                "7.2" = c("702","0702","72"),
                                "7.1" = c("701","0701","71"),
                                "6.8" = c("608","0608","68"),
                                "6.7" = c("607","0607","67"),
                                "6.6" = c("606","0606","66"),
                                "6.5" = c("605","0605","65"),
                                "6.4" = c("604", "6.4 (TIDM7)","0604","6.4 7","64"),
                                "6.3" = c("603","0603","63"),
                                "6.2" = c("602","0602","62"),
                                "6.1" = c("601",'0601',"61"),
                                "5.7" = c("507",'0507',"57"),
                                "5.6" = c("506","0506","56"),
                                "5.5" = c("505","0505","55"),
                                "5.4" = c("504",'0504',"54"),
                                "5.3" = c("503","0503","53"),
                                "5.2" = c("502","0502","52"),
                                "5.1" = c("501","0501","51"),
                                "4.9" = c("409","0409","49"),
                                "4.8" = c("408","0408","48"),
                                "4.7" = c("407",'0407',"47"),
                                "4.6" = c("406",'0406',"46"),
                                "4.5" = c("405", "4.5 (TIDM5)","TIDM5 (4.5)","0405","4.5 5","45","5 4.5"),
                                "4.4" = c("404",'0404',"44"),
                                "4.3" = c("403","0403","43"),
                                "4.2" = c("402",'0402',"42"),
                                "4.1" = c("401",'0401',"41"),
                                "3.9" = c("309",'0309',"39"),
                                "3.8" = c("308",'0308',"38"),
                                "3.7" = c("307","0307","37"),
                                "3.6" = c("306","0306","36"),
                                "3.5" = c("305","0305","35"),
                                "3.4" = c("304","0304","34"),
                                "3.3" = c("303","0303","33"),
                                "3.2" = c("302","0302","32"),
                                "3.1" = c("301","0301","31"),
                                "2.6" = c("206","0206","26"),
                                "2.5" = c("205","0205","25"),
                                "2.4" = c("204","0204","24"),
                                "2.3" = c("203","0203","23"),
                                "2.2" = c("202", "2.2 (TIDM2)","TIDM2 (2.2)","0202","22","2 2.2","2.2 2"),
                                "2.1" = c("201","0201","21"),
                                "1.5" = c("105","0105","15"),
                                "1.4" = c("104","0104","14"),
                                "1.3" = c("103","0103","13"),
                                "1.2" = c("102","0102","12"),
                                "1.1" = c("11","0101")))

unique(sgsl_ichthyo_env$Station) %>% sort()
a<-sgsl_ichthyo_env %>% group_by(Station) %>% dplyr::summarise(lat=mean(lat),lon=mean(lon))
a %>% ggplot(aes(lon,lat,label=Station))+geom_text(position=position_dodge(width=2))
sgsl_ichthyo_env$Station<-as.numeric(as.character(sgsl_ichthyo_env$Station))
sgsl_ichthyo_env %<>% dplyr::filter(Station<13)

sgsl_ichthyo_env %<>% transmute(year,month,day,doy=yday(date),hour,minute,station,Station,temperature_CL=temp_0_10m, salinity=sal_0_10m, sigma_t_kg_m3) %>% 
  dplyr::select(-station) %>% mutate(station=Station) %>% dplyr::select(-Station)

#multiple ctds at a given day/station in some cases...take mean
sgsl_ichthyo_env %<>% group_by(year,month,day,doy,station) %>% dplyr::summarise(hour=mean(hour),minute=mean(minute),temperature_CL=mean(temperature_CL),salinity=mean(salinity),sigma_t_kg_m3=mean(sigma_t_kg_m3))

# wrangle FG
names(sgsl_ichthyo_temp)
sgsl_ichthyo_temp %<>% dplyr::select(1:32) %>% 
  pivot_longer(cols = 5:32,names_to = "year", values_to = "temperature_FG") %>% 
  dplyr::select(year,station,strata,lat,lon,temperature_FG) %>% mutate(lon=-lon) %>% dplyr::filter(!is.na(station))
sgsl_ichthyo_temp$year<-as.numeric(sgsl_ichthyo_temp$year)

survey_env <- full_join(sgsl_ichthyo_env,sgsl_ichthyo_temp) %>% arrange(year,station)
# some multiple ctd on same day and station
names(survey_env)

## add 2021-2022
# copy paste dat file with datapasta:: and Rstudio addin
# dates missing from files Caroline sent so just use mission binder to manually put them in later with a left join
d2021<- data.frame(
   stringsAsFactors = FALSE,
        check.names = FALSE,
            Fichier = c("BioNet-IML2021014-8.7",
                        "BioNet-IML2021014-7.7","BioNet-IML2021014-7.6",
                        "BioNet-IML2021014-6.8","BioNet-IML2021014-6.7",
                        "BioNet-IML2021014-6.6","BioNet-IML2021014-5.6",
                        "BioNet-IML2021014-5.7","BioNet-IML2021014-4.9","BioNet-IML2021014-4.7",
                        "BioNet-IML2021014-3.6","BioNet-IML2021014-4.6",
                        "BioNet-IML2021014-3.7","BioNet-IML2021014-3.8",
                        "BioNet-IML2021014-3.9","BioNet-IML2021014-2.6",
                        "BioNet-IML2021014-2.5","BioNet-IML2021014-1.3","BioNet-IML2021014-2.4",
                        "BioNet-IML2021014-1.2","BioNet-IML2021014-2.3",
                        "BioNet-IML2021014-2.2","BioNet-IML2021014-2.1",
                        "BioNet-IML2021014-3.2","BioNet-IML2021014-3.1",
                        "BioNet-IML2021014-4.1","BioNet-IML2021014-4.2","BioNet-IML2021014-3.3",
                        "BioNet-IML2021014-3.4","BioNet-IML2021014-3.5",
                        "BioNet-IML2021014-4.5","BioNet-IML2021014-4.4",
                        "BioNet-IML2021014-4.3","BioNet-IML2021014-5.1",
                        "BioNet-IML2021014-5.2","BioNet-IML2021014-5.3","BioNet-IML2021014-5.4",
                        "BioNet-IML2021014-5.5","BioNet-IML2021014-6.5",
                        "BioNet-IML2021014-6.4","BioNet-IML2021014-6.3",
                        "BioNet-IML2021014-6.2","BioNet-IML2021014-6.1",
                        "BioNet-IML2021014-7.1","BioNet-IML2021014-7.2","BioNet-IML2021014-7.3",
                        "BioNet-IML2021014-8.3","BioNet-IML2021014-8.2",
                        "BioNet-IML2021014-8.1","BioNet-IML2021014-9.1",
                        "BioNet-IML2021014-9.2","BioNet-IML2021014-9.3",
                        "BioNet-IML2021014-8.4","BioNet-IML2021014-7.4","BioNet-IML2021014-7.5",
                        "BioNet-IML2021014-8.5","BioNet-IML2021014-9.4",
                        "BioNet-IML2021014-10.1","BioNet-IML2021014-11.1",
                        "BioNet-IML2021014-12.1","BioNet-IML2021014-9.5",
                        "BioNet-IML2021014-8.6"),
            station = c(8.7,7.7,7.6,6.8,6.7,6.6,
                        5.6,5.7,4.9,4.7,3.6,4.6,3.7,3.8,3.9,2.6,2.5,
                        1.3,2.4,1.2,2.3,2.2,2.1,3.2,3.1,4.1,4.2,3.3,
                        3.4,3.5,4.5,4.4,4.3,5.1,5.2,5.3,5.4,5.5,6.5,
                        6.4,6.3,6.2,6.1,7.1,7.2,7.3,8.3,8.2,8.1,9.1,9.2,
                        9.3,8.4,7.4,7.5,8.5,9.4,10.1,11.1,12.1,9.5,
                        8.6),
  `temperature_CL` = c(9.38,7.83,8.31,8.1,7.99,
                        8.85,7.9,7.17,7.17,9.58,9.75,10.07,9.45,7.91,
                        8.18,8.7,9.88,7.6,9.64,9.53,10.3,11.47,11.19,
                        11.66,12.93,11.96,12.35,11.83,11.1,10.25,11.12,11.47,
                        11.89,12.65,11.61,10.79,10.88,10.13,10.84,10.8,
                        11.67,12.37,13.39,12.59,12.54,11.6,13.08,13.18,
                        13.83,13.22,13.08,12.45,11.56,11.73,11.53,10.15,
                        13.98,13.27,13.29,11.72,10.95,11.2),
         `salinity` = c(27.9,29.99,28.81,29.66,
                        30.6,29.72,30.61,31.3,30.6,30.45,30.67,30.59,30.6,
                        31.07,31.27,31.19,30.69,30.9,30.76,30.79,30.62,
                        29.85,29.21,29.3,28.5,28.97,29.32,29.26,29.53,
                        30.71,29.61,29.65,28.99,28.72,28.82,29.44,29.41,
                        30.12,29.63,29.09,28.77,29.14,28.73,28.6,28.16,
                        28.42,27.81,27.99,28.26,28.23,27.65,28.02,28.01,
                        28.74,29.62,28.37,27.31,27.15,27.07,27.05,28.28,
                        27.61),
    `sigma_t_kg_m3` = c(21.51,23.37,22.37,23.07,
                        23.82,23.01,23.84,24.48,23.94,23.47,23.62,23.5,
                        23.61,24.2,24.32,24.18,23.61,24.11,23.7,23.74,
                        23.49,22.69,22.24,22.23,21.37,21.91,22.12,22.17,22.5,
                        23.56,22.56,22.53,21.95,21.6,21.86,22.48,22.44,
                        23.12,22.62,22.2,21.81,21.97,21.47,21.5,21.18,
                        21.56,20.82,20.94,21.01,21.11,20.68,21.09,21.24,
                        21.78,22.49,21.72,20.25,20.24,20.19,20.41,21.55,
                        21)
)


d2022 <- data.frame(
   stringsAsFactors = FALSE,
        check.names = FALSE,
                     Fichier = c("BioNet-IML2022024-11.1","BioNet-IML2022024-10.1",
                                 "BioNet-IML2022024-9.4","BioNet-IML2022024-9.5",
                                 "BioNet-IML2022024-8.5","BioNet-IML2022024-7.5",
                                 "BioNet-IML2022024-7.4",
                                 "BioNet-IML2022024-8.4","BioNet-IML2022024-9.3",
                                 "BioNet-IML2022024-9.2","BioNet-IML2022024-9.1",
                                 "BioNet-IML2022024-8.1","BioNet-IML2022024-8.2",
                                 "BioNet-IML2022024-8.3","BioNet-IML2022024-7.3",
                                 "BioNet-IML2022024-7.2","BioNet-IML2022024-7.1",
                                 "BioNet-IML2022024-6.1","BioNet-IML2022024-6.2",
                                 "BioNet-IML2022024-6.3","BioNet-IML2022024-6.4",
                                 "BioNet-IML2022024-5.3","BioNet-IML2022024-5.1",
                                 "BioNet-IML2022024-5.2",
                                 "BioNet-IML2022024-4.3","BioNet-IML2022024-4.4",
                                 "BioNet-IML2022024-4.5","BioNet-IML2022024-3.5",
                                 "BioNet-IML2022024-3.4","BioNet-IML2022024-3.3",
                                 "BioNet-IML2022024-4.2","BioNet-IML2022024-4.1",
                                 "BioNet-IML2022024-3.1","BioNet-IML2022024-3.2",
                                 "BioNet-IML2022024-2.1","BioNet-IML2022024-2.2",
                                 "BioNet-IML2022024-2.3","BioNet-IML2022024-1.2",
                                 "BioNet-IML2022024-1.3",
                                 "BioNet-IML2022024-2.4","BioNet-IML2022024-2.5",
                                 "BioNet-IML2022024-1.4","BioNet-IML2022024-4r28",
                                 "BioNet-IML2022024-4r29","BioNet-IML2022024-4r25",
                                 "BioNet-IML2022024-4r26","BioNet-IML2022024-4r27",
                                 "BioNet-IML2022024-4r24","BioNet-IML2022024-4r23",
                                 "BioNet-IML2022024-4r22",
                                 "BioNet-IML2022024-4r19","BioNet-IML2022024-4r17",
                                 "BioNet-IML2022024-4r15","BioNet-IML2022024-4r14",
                                 "BioNet-IML2022024-4r16","BioNet-IML2022024-1.5",
                                 "BioNet-IML2022024-2.6","BioNet-IML2022024-3.9",
                                 "BioNet-IML2022024-3.8","BioNet-IML2022024-3.7",
                                 "BioNet-IML2022024-3.6",
                                 "BioNet-IML2022024-4.6","BioNet-IML2022024-4.7",
                                 "BioNet-IML2022024-4.8","BioNet-IML2022024-4.9",
                                 "BioNet-IML2022024-5.7","BioNet-IML2022024-5.6",
                                 "BioNet-IML2022024-5.5","BioNet-IML2022024-5.4",
                                 "BioNet-IML2022024-6.5","BioNet-IML2022024-6.6",
                                 "BioNet-IML2022024-6.7","BioNet-IML2022024-6.8",
                                 "BioNet-IML2022024-7.7","BioNet-IML2022024-7.6",
                                 "BioNet-IML2022024-8.6","BioNet-IML2022024-8.7"),
                     station = c("11.1",
                                 "10.1","9.4","9.5","8.5","7.5","7.4","8.4",
                                 "9.3","9.2","9.1","8.1","8.2","8.3","7.3",
                                 "7.2","7.1","6.1","6.2","6.3","6.4",
                                 "5.3","5.1","5.2","4.3","4.4","4.5","3.5",
                                 "3.4","3.3","4.2","4.1","3.1","3.2","2.1",
                                 "2.2","2.3","1.2","1.3","2.4","2.5",
                                 "1.4","4r28","4r29","4r25","4r26","4r27",
                                 "4r24","4r23","4r22","4r19","4r17","4r15",
                                 "4r14","4r16","1.5","2.6","3.9","3.8","3.7",
                                 "3.6","4.6","4.7","4.8","4.9","5.7",
                                 "5.6","5.5","5.4","6.5","6.6","6.7","6.8",
                                 "7.7","7.6","8.6","8.7"),
           `temperature_CL` = c(12.86,
                                 11.89,12.61,8.48,11.43,10.51,11.77,9.31,
                                 12.16,12.4,13.77,12.68,13.54,13.17,12.81,
                                 13.96,13.27,13.16,13.49,13.12,10.27,11.09,
                                 13.7,13.77,14.03,11.42,11.2,9.64,12.48,
                                 13.79,13.06,11.95,12.61,12.99,12.3,12.63,
                                 9.49,10.92,9.56,9.92,9.23,10.15,10.47,
                                 10.27,10.31,11.3,10.68,11.11,11.46,10.39,
                                 10.66,11,11.05,10.35,9.96,11.8,11.03,11.48,
                                 11.64,11.96,12.86,11.16,12.16,11.96,
                                 11.99,10.69,11.92,12.25,12.64,12.24,11.44,
                                 10.57,11.86,11.9,10.23,10.75,10.93),
                  `salinity` = c(23.68,
                                 25.28,25.43,27.04,26.97,28.09,28.16,27.92,
                                 27.21,27.19,27.86,28.18,27.41,28.15,28.49,
                                 27.72,27.9,28.89,29.28,29.21,29.49,29,
                                 28.96,29.3,29.43,29.32,29.09,30.82,29.55,
                                 29.13,29.32,29.17,28.71,29.04,28.99,29.58,
                                 30.77,30.29,30.91,30.35,30.27,30.14,31,
                                 31.1,30.98,31,31.12,31.12,30.9,30.87,
                                 30.8,30.72,30.76,30.74,30.73,29.89,30.34,
                                 30.16,30.59,29.63,29.17,29.09,28.23,28.45,
                                 28.7,28.4,27.89,27.91,27.9,28.86,27.76,
                                 28.49,29.23,29.3,28.5,26.35,25.46),
             `sigma_t_kg_m3` = c(17.65,
                                 19.05,19.06,20.97,20.47,21.48,21.32,21.52,
                                 20.54,20.45,20.73,21.17,20.42,21.06,21.39,
                                 20.58,20.84,21.63,21.87,21.88,22.6,22.09,
                                 21.58,21.83,21.87,22.28,22.14,23.74,
                                 22.29,21.71,21.98,22.07,21.6,21.79,21.87,
                                 22.27,23.73,23.12,23.82,23.33,23.37,23.14,
                                 23.75,23.86,23.76,23.61,23.81,23.73,23.5,
                                 23.66,23.57,23.44,23.47,23.57,23.63,22.66,
                                 23.13,22.92,23.21,22.44,21.91,22.12,
                                 21.3,21.52,21.7,21.69,21.09,21.05,20.98,
                                 21.78,21.07,21.78,22.14,22.19,21.83,20.08,
                                 19.37)
         )

d2021 %<>% mutate(year = 2021, station = factor(station))
d2022 %<>% mutate(year = 2022, station = factor(station))
survey_env <- bind_rows(survey_env, d2021)
survey_env <- bind_rows(survey_env, d2022)
# save
save(survey_env, file = "./Rdata/survey_env.Rdata")
