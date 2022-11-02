#####################  SGSL ICHTHYOPLANKTON SURVEY FIGURES  ####################  

#####################  Load Packages  ####################  
packages = c('measurements','magrittr','forcats', 'lubridate', 'readxl','data.table','fuzzyjoin', 'stringi','tidyverse')
invisible(lapply(packages, function(x) {if (!require(x, character.only = T)) {install.packages(x);require(x)}}))
packages = c('magrittr','epitools', 'lubridate', 'data.table','gslea','formattable', 'sf', 'sp', 'tmap', 'rgdal','rgeos' ,'raster','maptools','mapdata','gridExtra','shapefiles','installr','tidyverse','cowplot','ggridges')
invisible(lapply(packages, function(x) {if (!require(x, character.only = T)) {install.packages(x);require(x)}}))

#####################  Load data  ####################
load("./Rdata/df_egg_index_2022_v4.Rdata")
df <- df_egg_index %>% 
  dplyr::filter(trajet == "T1")
rm(df_egg_index)
glimpse(df)

#####################  Figures  ####################

### Figure 0a: Base map sGSL
canada <- map_data("worldHires", "Canada")
canada_map <- ggplot() +
  geom_polygon(data = canada,aes(x = long, y = lat, group = group),
               fill = "ghostwhite", 
               color = "black") +
  theme_bw() 

sgsl_map <- canada_map +
  coord_map(xlim = c(-66.5, -60),
            ylim = c(45.7, 49.5),
            projection = "lambert",
            parameters = c(46.5, 48)) +
  ylab("Latitude") + 
  xlab("Longitude")

nl_map <- canada_map +
  coord_map(xlim = c(-60.5, -55),
            ylim = c(46.5, 49.5),
            projection = "lambert",
            parameters = c(48, 50)) +
  ylab("Latitude") + 
  xlab("Longitude")
rm(canada);rm(canada_map)

df_map <- df %>% 
  dplyr::group_by(station) %>% 
  dplyr::summarise(lon = mean(lon,na.rm=T),
                   lat = mean(lat,na.rm=T))
### Figure 0b: map NL stations

# df_map_nl <- df_nl %>% 
#   dplyr::group_by(station) %>% 
#   dplyr::summarise(lon = mean(lon,na.rm=T),
#                    lat = mean(lat,na.rm=T))

### Figure 1a: Map of stations
df_map <- df %>% 
  dplyr::group_by(station) %>% 
  dplyr::summarise(lon = mean(lon,na.rm=T),
                   lat = mean(lat,na.rm=T))
p1 <- sgsl_map +
  geom_text(data = df_map,
            aes(lon, lat, label = station)) +
  ggsave("survey_map.pdf", device = "pdf", path = "./img/survey/")

### Figure 1b: Mission dates 1979-2022
# need to validate the data set is parsed correctly for trajet and consecutives because this figure shows it is not
df1 <-
  df_egg_index %>% dplyr::filter(stage == "n_15") %>% 
  group_by(year, trajet, doy, station, consec) %>% dplyr::summarise(n = n())
df1 %<>% mutate(trajet = ifelse(is.na(trajet), "T1", trajet),
                doy = ifelse(is.na(doy), 170, doy)) # find out what 1979 dates were
d2022 <-
  data.frame(
    year = 2022,
    doy = c(166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177),
    n = 1,
    trajet = "T1",
    consec = 1:12, 
    station = factor(1)
  )
df1 <- bind_rows(df1, d2022)
jd <- c(152, 172, 182, 213) 

df1 %>% ggplot(aes(doy, year, fill = trajet)) + geom_tile(colour = "white") +
  theme_minimal() + scale_fill_manual(values = c("darkorange", "cyan")) +
  scale_y_continuous(n.breaks = 16) +
  geom_vline(xintercept = jd) +
  scale_x_continuous(n.breaks = 10) +
  facet_wrap(vars(trajet)) 
  annotate("text", x = jd, y = 1977, label = c("June 1", "Solstice", "July 1", "August 1")) 
  df_egg_index %>% dplyr::filter(year %in% c(1989), trajet == "T2") %>%  ggplot(aes(doy, year, colour = trajet)) + geom_point()
    
dcheck <- df1 %>% dplyr::filter(n>1)
dcheck %<>% mutate(check = "check",cc = "consec", consec = paste(cc,consec,sep = "")) %>% dplyr::select(-cc, -n)
dchecks <- sgsl_ichthyo %>% dplyr::filter(year %in% c(unique(dcheck$year)),
                                          taxons == "Scomber scombrus (oeuf stade 1)", 
                                          station %in% c(unique(dcheck$station)))
# they are all cases where both babord and tribord were taken. can either keep both and take mean or just keep babord. as for t2 in 2000...
d2000<-sgsl_ichthyo %>% dplyr::filter(year ==2000,taxons == "Scomber scombrus (oeuf stade 1)")
# for 2000 all event_collector_station names with b after are trajet 2. That leaves stations 6.2_B, 6.2_T, 7.1_B, and 7.1_T.   Choose babord and be done with it

### Figure 2a1: Annual stage 1+5 density (n/m^2)
library(ggforce)
cols <- c("1" = "red", "0" = "black")

p2a1 <- sgsl_map +
  geom_point(
    data = df %>% 
      dplyr::filter(stage == "n_15"),
    aes(lon, lat, size = n_m2, fill = factor(presence)), alpha = 0.5, shape = 21,) +
  scale_fill_manual(values = cols) +
  scale_size_binned_area() +
  facet_wrap_paginate(vars(year), nrow = 4, ncol = 3, page = 1) +
  ggsave("density_map_a1.pdf", device = "pdf", path = "./img/survey/")

p2a2 <- sgsl_map +
  geom_point(
    data = df %>% 
      dplyr::filter(stage == "n2"),
    aes(lon, lat, size = n_m2, fill = factor(presence)), alpha = 0.5, shape = 21,) +
  scale_fill_manual(values = cols) +
  scale_size_binned_area() +
  facet_wrap_paginate(vars(year), nrow = 4, ncol = 3, page = 1) +
  ggsave("density_map_a2.pdf", device = "pdf", path = "./img/survey/")

p2a3 <- sgsl_map +
  geom_point(
    data = df %>% 
      dplyr::filter(stage == "n3"),
    aes(lon, lat, size = n_m2, fill = factor(presence)), alpha = 0.5, shape = 21,) +
  scale_fill_manual(values = cols) +
  scale_size_binned_area() +
  facet_wrap_paginate(vars(year), nrow = 4, ncol = 3, page = 1) +
  ggsave("density_map_a3.pdf", device = "pdf", path = "./img/survey/")

p2a4 <- sgsl_map +
  geom_point(
    data = df %>% 
      dplyr::filter(stage == "n4"),
    aes(lon, lat, size = n_m2, fill = factor(presence)), alpha = 0.5, shape = 21,) +
  scale_fill_manual(values = cols) +
  scale_size_binned_area() +
  facet_wrap_paginate(vars(year), nrow = 4, ncol = 3, page = 1) +
  ggsave("density_map_a4.pdf", device = "pdf", path = "./img/survey/")

p2a5 <- sgsl_map +
  geom_point(
    data = df %>% 
      dplyr::filter(stage == "nl"),
    aes(lon, lat, size = n_m2, fill = factor(presence)), alpha = 0.5, shape = 21,) +
  scale_fill_manual(values = cols) +
  scale_size_binned_area() +
  facet_wrap_paginate(vars(year), nrow = 4, ncol = 3, page = 1) +
  ggsave("density_map_a5.pdf", device = "pdf", path = "./img/survey/")

### Figure 2b: Annual stage 2 density (n/m^2)
p2a1 <- sgsl_map +
  geom_point(
    data = df %>% 
      dplyr::filter(stage == "n_15"),
    aes(lon, lat, size = n_m2, fill = factor(presence)), alpha = 0.5, shape = 21,) +
  scale_fill_manual(values = cols) +
  scale_size_binned_area() +
  facet_wrap_paginate(vars(year), nrow = 4, ncol = 3, page = 1) +
  ggsave("density_map_a1.pdf", device = "pdf", path = "./img/survey/")

p2a2 <- sgsl_map +
  geom_point(
    data = df %>% 
      dplyr::filter(stage == "n2"),
    aes(lon, lat, size = n_m2, fill = factor(presence)), alpha = 0.5, shape = 21,) +
  scale_fill_manual(values = cols) +
  scale_size_binned_area() +
  facet_wrap_paginate(vars(year), nrow = 4, ncol = 3, page = 1) +
  ggsave("density_map_a2.pdf", device = "pdf", path = "./img/survey/")

p2a3 <- sgsl_map +
  geom_point(
    data = df %>% 
      dplyr::filter(stage == "n3"),
    aes(lon, lat, size = n_m2, fill = factor(presence)), alpha = 0.5, shape = 21,) +
  scale_fill_manual(values = cols) +
  scale_size_binned_area() +
  facet_wrap_paginate(vars(year), nrow = 4, ncol = 3, page = 1) +
  ggsave("density_map_a3.pdf", device = "pdf", path = "./img/survey/")

p2a4 <- sgsl_map +
  geom_point(
    data = df %>% 
      dplyr::filter(stage == "n4"),
    aes(lon, lat, size = n_m2, fill = factor(presence)), alpha = 0.5, shape = 21,) +
  scale_fill_manual(values = cols) +
  scale_size_binned_area() +
  facet_wrap_paginate(vars(year), nrow = 4, ncol = 3, page = 1) +
  ggsave("density_map_a4.pdf", device = "pdf", path = "./img/survey/")

p2a5 <- sgsl_map +
  geom_point(
    data = df %>% 
      dplyr::filter(stage == "nl"),
    aes(lon, lat, size = n_m2, fill = factor(presence)), alpha = 0.5, shape = 21,) +
  scale_fill_manual(values = cols) +
  scale_size_binned_area() +
  facet_wrap_paginate(vars(year), nrow = 4, ncol = 3, page = 1) +
  ggsave("density_map_a5.pdf", device = "pdf", path = "./img/survey/")

### Figure 2c: Annual stage 3 density (n/m^2)
p2a1 <- sgsl_map +
  geom_point(
    data = df %>% 
      dplyr::filter(stage == "n_15"),
    aes(lon, lat, size = n_m2, fill = factor(presence)), alpha = 0.5, shape = 21,) +
  scale_fill_manual(values = cols) +
  scale_size_binned_area() +
  facet_wrap_paginate(vars(year), nrow = 4, ncol = 3, page = 1) +
  ggsave("density_map_a1.pdf", device = "pdf", path = "./img/survey/")

p2a2 <- sgsl_map +
  geom_point(
    data = df %>% 
      dplyr::filter(stage == "n2"),
    aes(lon, lat, size = n_m2, fill = factor(presence)), alpha = 0.5, shape = 21,) +
  scale_fill_manual(values = cols) +
  scale_size_binned_area() +
  facet_wrap_paginate(vars(year), nrow = 4, ncol = 3, page = 1) +
  ggsave("density_map_a2.pdf", device = "pdf", path = "./img/survey/")

p2a3 <- sgsl_map +
  geom_point(
    data = df %>% 
      dplyr::filter(stage == "n3"),
    aes(lon, lat, size = n_m2, fill = factor(presence)), alpha = 0.5, shape = 21,) +
  scale_fill_manual(values = cols) +
  scale_size_binned_area() +
  facet_wrap_paginate(vars(year), nrow = 4, ncol = 3, page = 1) +
  ggsave("density_map_a3.pdf", device = "pdf", path = "./img/survey/")

p2a4 <- sgsl_map +
  geom_point(
    data = df %>% 
      dplyr::filter(stage == "n4"),
    aes(lon, lat, size = n_m2, fill = factor(presence)), alpha = 0.5, shape = 21,) +
  scale_fill_manual(values = cols) +
  scale_size_binned_area() +
  facet_wrap_paginate(vars(year), nrow = 4, ncol = 3, page = 1) +
  ggsave("density_map_a4.pdf", device = "pdf", path = "./img/survey/")

p2a5 <- sgsl_map +
  geom_point(
    data = df %>% 
      dplyr::filter(stage == "nl"),
    aes(lon, lat, size = n_m2, fill = factor(presence)), alpha = 0.5, shape = 21,) +
  scale_fill_manual(values = cols) +
  scale_size_binned_area() +
  facet_wrap_paginate(vars(year), nrow = 4, ncol = 3, page = 1) +
  ggsave("density_map_a5.pdf", device = "pdf", path = "./img/survey/")

### Figure 2d: Annual stage 4 density (n/m^2)
p2a1 <- sgsl_map +
  geom_point(
    data = df %>% 
      dplyr::filter(stage == "n_15"),
    aes(lon, lat, size = n_m2, fill = factor(presence)), alpha = 0.5, shape = 21,) +
  scale_fill_manual(values = cols) +
  scale_size_binned_area() +
  facet_wrap_paginate(vars(year), nrow = 4, ncol = 3, page = 1) +
  ggsave("density_map_a1.pdf", device = "pdf", path = "./img/survey/")

p2a2 <- sgsl_map +
  geom_point(
    data = df %>% 
      dplyr::filter(stage == "n2"),
    aes(lon, lat, size = n_m2, fill = factor(presence)), alpha = 0.5, shape = 21,) +
  scale_fill_manual(values = cols) +
  scale_size_binned_area() +
  facet_wrap_paginate(vars(year), nrow = 4, ncol = 3, page = 1) +
  ggsave("density_map_a2.pdf", device = "pdf", path = "./img/survey/")

p2a3 <- sgsl_map +
  geom_point(
    data = df %>% 
      dplyr::filter(stage == "n3"),
    aes(lon, lat, size = n_m2, fill = factor(presence)), alpha = 0.5, shape = 21,) +
  scale_fill_manual(values = cols) +
  scale_size_binned_area() +
  facet_wrap_paginate(vars(year), nrow = 4, ncol = 3, page = 1) +
  ggsave("density_map_a3.pdf", device = "pdf", path = "./img/survey/")

p2a4 <- sgsl_map +
  geom_point(
    data = df %>% 
      dplyr::filter(stage == "n4"),
    aes(lon, lat, size = n_m2, fill = factor(presence)), alpha = 0.5, shape = 21,) +
  scale_fill_manual(values = cols) +
  scale_size_binned_area() +
  facet_wrap_paginate(vars(year), nrow = 4, ncol = 3, page = 1) +
  ggsave("density_map_a4.pdf", device = "pdf", path = "./img/survey/")

p2a5 <- sgsl_map +
  geom_point(
    data = df %>% 
      dplyr::filter(stage == "nl"),
    aes(lon, lat, size = n_m2, fill = factor(presence)), alpha = 0.5, shape = 21,) +
  scale_fill_manual(values = cols) +
  scale_size_binned_area() +
  facet_wrap_paginate(vars(year), nrow = 4, ncol = 3, page = 1) +
  ggsave("density_map_a5.pdf", device = "pdf", path = "./img/survey/")
### Figure 2e: Annual larvae density (n/m^2)
p2a1 <- sgsl_map +
  geom_point(
    data = df %>% 
      dplyr::filter(stage == "n_15"),
    aes(lon, lat, size = n_m2, fill = factor(presence)), alpha = 0.5, shape = 21,) +
  scale_fill_manual(values = cols) +
  scale_size_binned_area() +
  facet_wrap_paginate(vars(year), nrow = 4, ncol = 3, page = 1) +
  ggsave("density_map_a1.pdf", device = "pdf", path = "./img/survey/")

p2a2 <- sgsl_map +
  geom_point(
    data = df %>% 
      dplyr::filter(stage == "n2"),
    aes(lon, lat, size = n_m2, fill = factor(presence)), alpha = 0.5, shape = 21,) +
  scale_fill_manual(values = cols) +
  scale_size_binned_area() +
  facet_wrap_paginate(vars(year), nrow = 4, ncol = 3, page = 1) +
  ggsave("density_map_a2.pdf", device = "pdf", path = "./img/survey/")

p2a3 <- sgsl_map +
  geom_point(
    data = df %>% 
      dplyr::filter(stage == "n3"),
    aes(lon, lat, size = n_m2, fill = factor(presence)), alpha = 0.5, shape = 21,) +
  scale_fill_manual(values = cols) +
  scale_size_binned_area() +
  facet_wrap_paginate(vars(year), nrow = 4, ncol = 3, page = 1) +
  ggsave("density_map_a3.pdf", device = "pdf", path = "./img/survey/")

p2a4 <- sgsl_map +
  geom_point(
    data = df %>% 
      dplyr::filter(stage == "n4"),
    aes(lon, lat, size = n_m2, fill = factor(presence)), alpha = 0.5, shape = 21,) +
  scale_fill_manual(values = cols) +
  scale_size_binned_area() +
  facet_wrap_paginate(vars(year), nrow = 4, ncol = 3, page = 1) +
  ggsave("density_map_a4.pdf", device = "pdf", path = "./img/survey/")

p2a5 <- sgsl_map +
  geom_point(
    data = df %>% 
      dplyr::filter(stage == "nl"),
    aes(lon, lat, size = n_m2, fill = factor(presence)), alpha = 0.5, shape = 21,) +
  scale_fill_manual(values = cols) +
  scale_size_binned_area() +
  facet_wrap_paginate(vars(year), nrow = 4, ncol = 3, page = 1) +
  ggsave("density_map_a5.pdf", device = "pdf", path = "./img/survey/")

### Figure 3: Volume filtered vs time
# look at read_ichthyo.R again because set_duration wonky
df %>% dplyr::filter(set_duration > 0, set_duration < 2000) %>%
  ggplot(aes(set_duration, volume)) + 
  geom_point() + 
  geom_smooth(method = "lm") +
  facet_wrap(vars(year),scales = "free")
