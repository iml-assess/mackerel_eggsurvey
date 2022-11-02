#####################  SGSL ICHTHYOPLANKTON SURVEY 1979-  ####################  
##  Version Log:
##  -v.0.0: Prior to 2020. 1982-2019 time series kept in .xls and .xlsx formats with different column headers, positions, formats, etc. (i.e. a headache)
##  -V.0.5: Winter and Spring 2020. In preperation for capelin stock assessment, I (Andrew Smith) made several attempts to clean the data via brute force (1000+ lines of code) to have a standardized format. It still had errors.
##  -V.1.0: 
##   ...Helena Talbot, Isabelle St. Pierre, and Caroline Lafleur, at the DAISS, began to clean and validate the mission data to store in the SGDE (oceanographic data) and Biochem (national icthyoplankton/biochemical data) databases.
##   ... I entered into contact with Isabelle and Helena and got a tour of their efforts. Yay! organised data! ... Collectively, we each tried, and still... the horrors, the horrors of excel. This is well documented, from bad formating, to out of date shoddy statistics, to vulnerability to user error.
##   ... As of now, they only have 1989 to 2019. I have to give them 1982-1988. I also have to find the 1979 data. Last resort is scan from resdoc.
##   ... mid october update:
##   ... - successfully made figures
##   ... - added 1982-1988 but these are not in BIOCHEM yet, just added here
##   ... - 1979 data not in any digital form I can find. It is in the dfo res doc 2006/099
##   ... 1993 and 1998 data also missing but I have alerted DAISS
##   ... data found but early missions mission and plankton logs missing
##   -v.1.5: 
##   ... 1993, 1998 data added. 2005 and 2006 updated (database uploads (DAISS) and email attached .xlsx and .xls files from DAISS)
##   ...random thought: oceanographers use matlab, ecology/biology/pop dynamics, and genomics use R. - way to bridge the gap?
##   ...noticed that fraction sampled is uniquely the one for mackerel eggs and in recent years the fractioned sample may have been different for the other eggs and for larvae due to high or low numbers... to discuss with DAISS
##   ... nov. 17. 2020 - 1993 and 1987 updated
##  Purpose of script: Read the survey metadata and species data.
##  Reasons why this script may not work:
##  - you have the original files open in excel (close them!)
##  - you do not have the latest versions of R, Rstudio, or the packages listed in the "# Load packages" section below
##  -V.1.1: 
##  - May 2021: 82-87 now in biochem so merge and clean. Try to keep original variable names if possible...
##  - fixed station id and pass identification problems
##
## Unique ID: The unique identifier between the metadata (sgsl_ichthyo_meta) and taxon (sgsl_ichthyo_taxon) specific data = "PLANK_SAMPLE_KEY_VALUE"
## CYT = Tanche-Tautogue & Limande, CHW = Morue, Aiglefin, et PLie Grise Et H4B = Motelle à 4 Barbillons Merluche sp. et Stromatée à fosettes
##
## Supplemental info (from wiki etc):
##  
##  total GSL area = 226000 km2
##  total GSL volume = 34500 km3
##  mean depth of GSL = 152 m 
##  mean depth sampled 44.20194 from sgsl_ichthyo
##  mean depth sampled from cap_larvae = 42.66987
##  official survey surface area used in most published documents = 6.945e+10 m2 (69 450 km2)
##  a "geo-referenced" approach was used by Francois Grégoire in conjunction with the "official" value above. = 8.21e+10
##  - Most likely this value is more accurate. validate in QGIS and buffer littoral of the sGSL shapefile.
##  - I checked in QGIS, the survey area from the official 'grille maquereau' is actually even bigger than the 'geo-referenced' version ~ 83485 km^2


## 05/10/2021 
## Element to add 
## Trajet 1 or 2? The variable to use is PL_HEADR_COLLECTOR_DEPLMT_ID.
## If value in there has a suffix T1, then trajet 1, if T2 then trajet 2

## 18/06/2022
## on teleost; mackerel egg survey (15th-29th)


## 24/06/2022
# On CCGS Teleost for survey as chief scientist
# Setting up scripts for next eval and all the analyses
# comparing data from published resdocs, .xls in /plancton/, and Biochem extracts
# 1982 data volumnes and thus n_m2 differ. consecutives and N match as do depths but many volumes differ. Dates slightly different from published as well. As are event times (conversion from NL time to UTC?)

# updated all: full data set, reduced to main set, and added 1979 and 2021. 
# remaining is Temperatures

#####################  Load Packages  ####################  
packages = c('measurements','magrittr','forcats', 'lubridate', 'readxl','data.table','fuzzyjoin', 'stringi','tidyverse')
invisible(lapply(packages, function(x) {if (!require(x, character.only = T)) {install.packages(x);require(x)}}))

#####################  Load Metadata  ####################

# BioChem data dictionnary
biochem_metadata <- read_csv("data/biochem_1982-/metadata/biochem_metadata.csv")

# Save your original working directory to call upon later (e.g. if you are in an .RProject with dependencies etc.) 
project_wd <- "C:/Users/SmithAND/Documents/Data/Ichthyoplankton/"

# Set the working directory to the one from which the files will be read. You can switch back later by calling setwd(project_wd)
setwd("C:/Users/SmithAND/Documents/Data/Ichthyoplankton/data/biochem_1982-")

# Create a list of the files from your target directory. if you stored the data in a folder with other files this will fail unless you specify a pattern by which files should be filtered (e.g. pattern = biochem)
file_list <- list.files(path = "C:/Users/SmithAND/Documents/Data/Ichthyoplankton/data/biochem_1982-")

# Read Mission Metadata

# Initiate a blank data frame, each iteration of the loop will append the data from the given file to this variable. there are tidyverse and data.table ways to do this but it is not a heavy operation so voilà
sgsl_ichthyo_meta <- data.frame()

# Loop to read the mission metadata
for (i in 1:length(file_list)) { 
  temp_data <- read_excel(file_list[i], sheet = 1, range = cell_cols("A:BI"), guess_max = 21474836) # each file will be read in, specify which columns you need and increase the functions guessing range to avoid any errors
  temp_data$data_tab <- file_list[i] # clean the data as needed, in this case I am creating a new column that indicates which file each row of data came from
  sgsl_ichthyo_meta <- bind_rows(sgsl_ichthyo_meta, temp_data) # for each iteration, bind the new data to the dataset
  rm(temp_data)
  sgsl_ichthyo_meta <- sgsl_ichthyo_meta[,1:62]
  sgsl_ichthyo_meta
}
# will post an error if subfolders present but that is ok

#####################  Format Mission Metadata   ####################
sgsl_ichthyo_meta <- sgsl_ichthyo_meta[, colSums(is.na(sgsl_ichthyo_meta)) < nrow(sgsl_ichthyo_meta)] # remove NA columns
sgsl_ichthyo_meta$PL_HEADR_VOLUME <- as.numeric(sgsl_ichthyo_meta$PL_HEADR_VOLUME) # flowmeter or bionet selon the year
sgsl_ichthyo_meta$PL_HEADR_START_DEPTH <- as.numeric(sgsl_ichthyo_meta$PL_HEADR_START_DEPTH) # max depth of bongo
sgsl_ichthyo_meta$PL_HEADR_END_DEPTH <- as.numeric(sgsl_ichthyo_meta$PL_HEADR_END_DEPTH) # surface. should be zero
sgsl_ichthyo_meta$EVENT_MIN_LAT <- as.numeric(sgsl_ichthyo_meta$EVENT_MIN_LAT) # often just station coords, especially early in time series
sgsl_ichthyo_meta$EVENT_MIN_LON <- as.numeric(sgsl_ichthyo_meta$EVENT_MIN_LON)
sgsl_ichthyo_meta$EVENT_MAX_LAT <- as.numeric(sgsl_ichthyo_meta$EVENT_MAX_LAT)
sgsl_ichthyo_meta$EVENT_MAX_LON <- as.numeric(sgsl_ichthyo_meta$EVENT_MAX_LON)
sgsl_ichthyo_meta %<>%
  mutate(lat = EVENT_MIN_LAT, lon = EVENT_MIN_LON)
sgsl_ichthyo_meta$MISSION_SDATE <- as.Date(sgsl_ichthyo_meta$MISSION_SDATE, format = "%d/%m/%Y")
sgsl_ichthyo_meta$MISSION_EDATE <- as.Date(sgsl_ichthyo_meta$MISSION_EDATE, format = "%d/%m/%Y")
sgsl_ichthyo_meta$EVENT_SDATE <- as.Date(sgsl_ichthyo_meta$EVENT_SDATE, format = "%d/%m/%Y")
sgsl_ichthyo_meta$EVENT_EDATE <- as.Date(sgsl_ichthyo_meta$EVENT_EDATE, format = "%d/%m/%Y")
sgsl_ichthyo_meta$hour_start <- as.numeric(str_sub(sgsl_ichthyo_meta$EVENT_STIME, 1,2)) # convert for convenience and to calculate duration of set
sgsl_ichthyo_meta$minute_start <- as.numeric(str_sub(sgsl_ichthyo_meta$EVENT_STIME, 3,4))
sgsl_ichthyo_meta$hour_end <- as.numeric(str_sub(sgsl_ichthyo_meta$EVENT_ETIME, 1,2)) 
sgsl_ichthyo_meta$minute_end <- as.numeric(str_sub(sgsl_ichthyo_meta$EVENT_ETIME, 3,4)) 
sgsl_ichthyo_meta %<>% mutate(hour_start = ifelse(hour_start == 0,24, hour_start),
                             hour_end = ifelse(hour_end == 0,24, hour_end), 
                             set_duration = ((hour_end*3600) + (minute_end*60)) - ((hour_start*3600) + (minute_start*60))) # calculate length of time bongos were sampling
sgsl_ichthyo_meta$day <- day(sgsl_ichthyo_meta$EVENT_SDATE)
sgsl_ichthyo_meta$month <- month(sgsl_ichthyo_meta$EVENT_SDATE)
sgsl_ichthyo_meta$year <- year(sgsl_ichthyo_meta$EVENT_SDATE)
sgsl_ichthyo_meta$doy <- yday(sgsl_ichthyo_meta$EVENT_SDATE)

# datetime
sgsl_ichthyo_meta %<>% 
  mutate(date = as_datetime(paste(year,month,day,hour_start,minute_start,sep = "/"), format = "%Y/%m/%d/%H/%M"))

#####################  Load Mission Data  ####################

# -initiate a blank data frame
sgsl_ichthyo_taxon <- data.frame()

# Loop to read over listed data files
for (i in 1:length(file_list)) {
  temp_data <- read_excel(file_list[i],
                          sheet = 2,
                          range = cell_cols("A:AB"), 
                          col_types = c("text", "text", "text", "text", "text", "text",
                                        "text", "text", "text", "text", "text", "text",
                                        "text", "text", "text", "text", "text", "text", 
                                        "text", "text", "text", "text", "text", "text",
                                        "text", "text", "text", "text")) #each file will be read in, specify which columns you need read in to avoid any errors) 
  temp_data$data_tab <- file_list[i]                                     #clean the data as needed, in this case I am creating a new column that indicates which file each row of data came from
  sgsl_ichthyo_taxon <- bind_rows(sgsl_ichthyo_taxon, temp_data)         #for each iteration, bind the new data to the building dataset
  rm(temp_data)
  rm(i)
}
# will post an error if subfolders present but that is ok

#####################  Format Mission Data  ####################

sgsl_ichthyo_taxon <- sgsl_ichthyo_taxon[, colSums(is.na(sgsl_ichthyo_taxon)) < nrow(sgsl_ichthyo_taxon)] # remove columns that are all NA
sgsl_ichthyo_taxon$PL_GEN_SPLIT_FRACTION <- as.numeric(sgsl_ichthyo_taxon$PL_GEN_SPLIT_FRACTION) # fractioned volume analysed
sgsl_ichthyo_taxon$PL_GEN_COUNTS <- as.numeric(sgsl_ichthyo_taxon$PL_GEN_COUNTS) # N; make sure to divide this by PL_GEN_SPLIT_FRACTION to get true N
sgsl_ichthyo_taxon$PL_GEN_BIO_VOLUME <- as.numeric(sgsl_ichthyo_taxon$PL_GEN_BIO_VOLUME) # settled volume. Don't use for abundance calculation

#####################  Merge Metadata And Data  ####################

# remove superfluous columns and merge
sgsl_ichthyo <- left_join(sgsl_ichthyo_taxon %>%
                            dplyr::select(-MISSION_DESCRIPTOR,
                                          -CREATED_BY,
                                          -CREATED_DATE,
                                          -DATA_CENTER_CODE,
                                          -BATCH_SEQ,
                                          -PL_GEN_TROPHIC_SEQ,
                                          -PL_GEN_NATIONAL_TAXONOMIC_SEQ,
                                          -PL_GEN_LIFE_HISTORY_SEQ,
                                          -PL_GEN_SEX_SEQ,
                                          -PL_GEN_PRESENCE,
                                          -DATA_CENTER_CODE,
                                          -PL_GEN_MODIFIER,
                                          -data_tab),
                          sgsl_ichthyo_meta %>% 
                            dplyr::select(-PL_HEADR_GEAR_SEQ,
                                          -PL_HEADR_TIME_QC_CODE,
                                          PL_HEADR_POSITION_QC_CODE,
                                          -PL_HEADR_MESH_SIZE,
                                          -PL_HEADR_COLLECTION_METHOD_SEQ,
                                          -PL_HEADR_PRESERVATION_SEQ,
                                          -PL_HEADR_METERS_SQD_FLAG,
                                          -DATA_CENTER_CODE,
                                          -BATCH_SEQ,
                                          -EVENT_MAX_LAT,
                                          -EVENT_MIN_LAT,
                                          -EVENT_MAX_LON,
                                          -EVENT_MIN_LON))
rm(sgsl_ichthyo_meta);rm(sgsl_ichthyo_taxon)
setwd(project_wd)
#####################  Formatting And Fixing Names  ####################

# lower case
names(sgsl_ichthyo) <- tolower(names(sgsl_ichthyo))

# simplify consecutive variable (consec = mission specific order the sets were done in, letters at the end indicate a second set or pass at the same station or using both tribord and babord bongo nets etc)
sgsl_ichthyo %<>% mutate(consec = str_remove(event_collector_event_id,"consec"))

# check for double occurrences of consecs (I checked this and it corresponds to a second pass generally or if a station was done again due to some deployment issue)
sgsl_ichthyo  %<>% 
  mutate(
    extra_consec = case_when(
      stringi::stri_detect_regex(consec, "[abcdefghijklmnopqrstuvwxyz]") ~ "0",
      TRUE ~ "1"
    )
  )

# In some years, station names were appended with various suffixes 
# '' | 'a' indicates first pass
# 'b' & 'c' indicate 2nd and 3rd pass respectively (only one case of the latter)
# '_B' & '_T' indicate babord and tribord respectively (i.e. which side of bongo net was used for the sample. Default is babord)
# 'COR' | 'cor' indicate a correction (not sure if a second pass at the same station or just a new line of data)

# Remove all lowercase characters (keep upper case for now because they will be used later for filtering)
sgsl_ichthyo  %<>% 
  mutate(
    station = str_remove_all(event_collector_stn_name, "[abcdefghijklmnopqrstuvwxyz]")
    )

# Essentially anything alphanumeric is a station from the second pass or a special experiment so create a new variable to ID them
sgsl_ichthyo  %<>% 
  mutate(
    extra_station = case_when(
      stringi::stri_detect_regex(event_collector_stn_name, "[bcdefghijklmnopqrstuvwxyz]", case_insensitive=F) ~ "0", # notice that 'a' is removed here
      TRUE ~ "1"
    )
  )

# Determine passes (i.e. trajet) via suffix, defaults to = 1, i.e. when there is no suffix then it is first pass
sgsl_ichthyo  %<>% 
  mutate(
    pass1 = case_when(
      stringi::stri_endswith_fixed(event_collector_stn_name, "a",case_insensitive=F) ~ "1",
      stringi::stri_endswith_fixed(event_collector_stn_name, "b",case_insensitive=F) ~ "2",
      stringi::stri_detect_regex(event_collector_stn_name, "b_") ~ "2",
      stringi::stri_endswith_fixed(event_collector_stn_name, "c",case_insensitive=F) ~ "3",
      TRUE ~ "1"
    )
  )

# Determine passes via order of occurrence (i.e. order by date then assign 1 to first occurrence of a station and 2 to last occurrence of the same station). 
# Note, I checked the counts of stations and all counts are either 1 or 2 except for one invertebrate taxon in 1989 where n=4. df <- sgsl_ichthyo  %>%  arrange(date,taxons) %>% group_by(year, event_collector_stn_name, taxons) %>% dplyr::summarise(n=n())
sgsl_ichthyo %<>% arrange(date,taxons) %>% 
  group_by(year, station, taxons) %>% 
  mutate(pass2 = case_when(
    doy == min(doy)~ "1",
    doy == max(doy)~ "2")) %>% ungroup()

# Combine both sources of pass information and create new variable 
sgsl_ichthyo %<>% mutate(pass = case_when(pass1=="1"&pass2=="1"~"1", # if both are 1 then 1
                                          pass1=="2"|pass2=="2"~"2", # if either indicate 2 then 2
                                          is.na(pass2)~pass1))       # if NA then whatever pass1 says 
  


# actual variable for trajet 1 vs 2

# function to extract last n characters
substrRight <- function(x, n){
  substr(x, nchar(x) - n + 1, nchar(x))
}

sgsl_ichthyo %<>%
  mutate(trajet = substrRight(pl_headr_collector_deplmt_id,2))
unique(sgsl_ichthyo$trajet)

# simplify
sgsl_ichthyo %<>% mutate(trajet2 = 
                           fct_collapse(trajet, 
                                        "T1" = c("RT", "T1"),
                                        "T2" = c("T2", "T3")))
# Is there a station correction?
sgsl_ichthyo  %<>% 
  mutate(
    correction = case_when(
      stringi::stri_detect_regex(event_collector_stn_name, "cor",case_insensitive=T) ~ "1",
      TRUE ~ "0"
    )
  )

# Side of Bongo net used for measure (i.e. babord or tribord) based on comments. Default is babord
sgsl_ichthyo  %<>% 
  mutate(
    bongo1 = case_when(
      stringi::stri_detect_regex(pl_headr_data_manager_comment, "babord",case_insensitive=F) ~ "babord",
      stringi::stri_detect_regex(pl_headr_data_manager_comment, "tribord",case_insensitive=F) ~ "tribord",
      TRUE ~ "babord"
    )
  ) # assumes that NAs are babord

# Side of Bongo net used for measure (i.e. babord or tribord) based on station name 
sgsl_ichthyo  %<>% 
  mutate(
    bongo2 = case_when(
      stringi::stri_detect_regex(event_collector_stn_name, "_B",case_insensitive=F) ~ "babord",
      stringi::stri_detect_regex(event_collector_stn_name, "_T",case_insensitive=F) ~ "tribord"
    )
  )

sgsl_ichthyo %<>% mutate(bongo = ifelse(is.na(bongo1),bongo2,bongo1)) # defaults to data manager comment variable

# Make sure everything checks out (look for NAs and cases where variables don't match)
# df <- sgsl_ichthyo %>%  dplyr::select(date, consec, extra_station, station, event_collector_stn_name, pass1, pass2, pass, trajet, trajet2)

# Remove clutter
# rm(df)
sgsl_ichthyo %<>% dplyr::select(-pass1,-pass2,-bongo2,-bongo1, trajet) # remove quality control variables

# There are a few additional missions and/or experiments that are found in this dataset 
# Make easy way to find them 1 = yes, 0 = no

# 2000 weird so hardcode from planksamplekey 
a<-c("2000024019_1","2000024020_1","2000024033_1","2000024034_1","2000024039_1","2000024040_1","2000024042_1","2000024044_1","2000024045_1","2000024046_1","2000024049_1","2000024050_1","2000024056_1","2000024057_1","2000024059_1","2000024060_1","2000024062_1","2000024063_1","2000024067_1","2000024068_1","2000024070_1","2000024071_1","2000024076_1","2000024077_1","2000024079_1","2000024080_1","2000024082_1","2000024083_1","2000024107_1","2000024108_1","2000024065_1","2000024066_1")

sgsl_ichthyo %<>% 
  mutate(extra_1991a = ifelse(event_collector_stn_name %in%
                                    c("7.5.00", "7.5.01", "7.5.02", "7.5.03", "7.5.04", "7.5.05", "7.5.06", "7.5.07", "7.5.08", "7.5.09", "7.5.10",
                                      "7.5.11", "7.5.12", "7.5.13", "7.5.14", "7.5.15", "7.5.16", "7.5.17", "7.5.18", "7.5.19", "7.5.20", "7.5.21",
                                      "7.5.22", "7.5.23", "7.5.24", "7.5.25", "7.5.26", "7.5.27", "7.5.28", "7.5.29", "7.5.30", "7.5.31", "7.5.32",
                                      "7.5.33", "7.5.34", "7.5.35"), 1,0),
         extra_1991b = ifelse(year == 1991 & event_collector_stn_name %in% 
                                        c("7.4A", "7.4B", "7.4C", "7.4D", "7.4E", "7.4F", "7.4G", "7.4H","7.4I","7.4J", "8.4A", "8.4B", "8.4C", "8.4D","8.4E","8.4F","8.4G","8.4H","8.4I"),1,0),
         extra_1986 = ifelse(year == 1986 & event_collector_stn_name %in% c("13.1","14.1","14.2","14.3","14.4","15.1","15.2","15.3","15.4","15.5","16.1","16.2","16.3","16.4"),1,0), 
         extra_northumberland = ifelse(event_collector_stn_name %in% c("5.0","6.0","7.0","8.0","9.0"),1,0), 
         extra_ns_1998 = ifelse(year == 1998 & str_starts(event_collector_stn_name,"SMB"),1,0), 
         extra_st_georges_bay_1994 = ifelse(year == 1994 & event_collector_stn_name %in% c("3.11","3.12","3.13","3.14","3.15","3.16","3.17","3.18","3.19","3.20","3.21","3.22","3.23","3.24","3.25","3.26","3.27","3.28","3.29"),1,0), 
         extra_2000 = ifelse(plank_sample_key_value %in% c(a),1,0), 
         extra_2000a = ifelse(year == 2000 & str_starts(event_collector_stn_name,"Supp"),1,0),
         extra_scotian_shelf_2009 = ifelse(mission_name == "IML-2009-47 Relevé international maquereau",1,0))



# Define new station variable to be able to merge with other data (e.g. temperature from ctd casts)
sgsl_ichthyo  %<>% 
  mutate(
    station_name = str_remove_all(event_collector_stn_name, "_B")) %>% 
  mutate(
    station_name = str_remove_all(station_name, "_T"))

sgsl_ichthyo  %<>% 
  mutate(
    station = str_remove_all(station, "_B")) %>% 
  mutate(
    station = str_remove_all(station, "_T")) %>% 
  mutate(
    station = str_remove_all(station, " COR")) %>% 
  mutate(
    station = str_remove_all(station, "A")) %>% 
  mutate(
    station = str_remove_all(station, "b"))

sgsl_ichthyo %<>% mutate(station = 
                           fct_collapse(station, 
                                        "7.7" = c("7.7", "7.71")))

# Strata and surface areas- defined by P. Ouellet, ResDoc 87/62
sgsl_ichthyo %<>% mutate(stratum = ifelse(station %in%
                                            c("8.7", "7.7", "5.7","4.9","3.9","3.8","2.6","1.5","2.5","1.4","2.4","1.3","3.5","3.4",'3.3',"3.2",'3.1',"4.2","4.1","2.3","2.2",'2.1',"1.2","1.1"), "1", 
                                          ifelse(station %in%
                                                   c("6.2_B", "7.1B", "10.1","9.4","9.3",'8.2','7.2','7.1',"6.3","6.2","6.1","5.2","5.1",'4.3',"4.4","4.5",'4.6',"3.6","3.7","4.8","5.6",'6.7',"6.8","7.5","7.6","6.6"), '2', "3"))) 

sgsl_ichthyo %<>% mutate(stratum_area = ifelse(stratum == "1",29.61e+9,
                                               ifelse(stratum == "2", 21.91e+9,
                                                      ifelse(stratum == "3", 93e+9, NA))))

# remove taxon = total sample... weird column done in some years
sgsl_ichthyo %<>% dplyr::filter(taxons != "Total sample") %>% droplevels()

#####################  Save Full Dataset  ####################
setwd(project_wd)
save(sgsl_ichthyo, file = "./Rdata/sgsl_ichthyo_2022_v2.Rdata") # full dataset

#####################  Subset to mackerel survey main set ####################

# load("./Rdata/sgsl_ichthyo_2022_v2.Rdata")

egg_index <- sgsl_ichthyo %>% 
  dplyr::filter(taxons %in% c("Scomber scombrus (oeuf stade 1)",
                              "Scomber scombrus (oeuf stade 2)",
                              "Scomber scombrus (oeuf stade 3)",
                              "Scomber scombrus (oeuf stade 4)",
                              "Scomber scombrus (oeuf stade 5)",
                              "Scomber scombrus (larva)"),
                extra_1991a == 0,
                extra_1991b == 0,
                extra_1986 == 0,
                extra_northumberland == 0,
                extra_ns_1998 == 0,
                extra_st_georges_bay_1994 == 0,
                extra_2000 == 0,
                extra_2000a == 0,
                extra_scotian_shelf_2009 == 0
                ) %>% droplevels()
unique(egg_index$station) # should be 66
unique(egg_index$year)
table(egg_index$year, egg_index$station)

egg_index %<>% transmute(unique_id = plank_sample_key_value,
              year = as.numeric(year), 
              month = as.numeric(month),
              day = as.numeric(day),
              doy = as.numeric(doy), 
              date = date,
              hour_start = hour_start,
              minute_start = minute_start,
              set_duration = set_duration,
              trajet = as.factor(trajet2),
              consec = as.factor(consec), 
              station = as.factor(station),
              correction = as.factor(correction),
              extra_station = as.factor(extra_station),
              extra_consec = as.factor(extra_consec),
              station_depth = pl_headr_sounding,
              lon = lon, 
              lat = lat,
              stratum = as.factor(stratum), 
              stratum_area = stratum_area,
              fraction = pl_gen_split_fraction,
              volume = pl_headr_volume,
              sample_depth = pl_headr_start_depth,
              large_plankton_removed = pl_headr_lrg_plankton_removed, 
              set_comments = pl_headr_collector_comment,
              sample_comments = pl_headr_data_manager_comment,
              mission_comment1 = mission_collector_comment,
              mission_comment2 = mission_more_comment,
              stage = taxons,
              count = pl_gen_counts)
glimpse(egg_index)

egg_index %<>% pivot_wider(names_from = "stage", values_from = "count") %>% 
  mutate(n1 = `Scomber scombrus (oeuf stade 1)`/fraction,
         n2 = `Scomber scombrus (oeuf stade 2)`/fraction,
         n3 = `Scomber scombrus (oeuf stade 3)`/fraction,
         n4 = `Scomber scombrus (oeuf stade 4)`/fraction,
         n5 = `Scomber scombrus (oeuf stade 5)`/fraction,
         nl = `Scomber scombrus (larva)`/fraction,
         n_15 = n1 + n5,
         total = n1+n2+n3+n4+n5) 
egg_index %<>% pivot_longer(cols = 29:42, names_to = "stage", values_to = "n") %>% 
  mutate(n_m3 = n/volume, 
         n_m2 = n_m3*sample_depth)
glimpse(egg_index)
egg_index$stage <- as.factor(egg_index$stage)
egg_index  %<>% 
  mutate(
    consecutive = str_remove_all(consec, "[abcdefghijklmnopqrstuvwxyz]"),
    consec = as.numeric(consec)
  )

# save(egg_index, file = "./Rdata/egg_index_2022_v1.Rdata") # currently 1982-2019
# save(egg_index, file = "./Rdata/egg_index_2022_v2.Rdata")
save(egg_index, file = "./Rdata/egg_index_2022_v3.Rdata")
load("./Rdata/egg_index_2022_v3.Rdata")              

# add 1979 and 2021   
# data from SGDO, RESDOCs 2006/099 and 2014/075, and papers by JJ Maguire and P. Ouellet
# 1979 
X1979 <- read_excel("data/larvae/New folder/1979.xlsx", sheet = "data")

X1979 %<>% 
  mutate(fraction = 1,
         consec = as.factor(consec),
         station = as.factor(station), 
         trajet = as.factor("T1")) %>% 
  pivot_longer(cols = 8:9, names_to = "stage", values_to = "count") %>% 
  mutate(n = count/fraction,
         n_m3 = n/volume, 
         n_m2 = n_m3*sample_depth) # assumes fraction = 1

X1979 %<>% 
  mutate(stratum = as.factor(ifelse(station %in% 
                                      c("8.7", "7.7", "5.7","4.9","3.9","3.8","2.6","1.5","2.5","1.4","2.4","1.3","3.5","3.4",'3.3',"3.2",'3.1',"4.2","4.1","2.3","2.2",'2.1',"1.2","1.1"), "1", 
                                    ifelse(station %in% 
                                             c("10.1","9.4","9.3",'8.2','7.2','7.1',"6.3","6.2","6.1","5.2","5.1",'4.3',"4.4","4.5",'4.6',"3.6","3.7","4.8","5.6",'6.7',"7.5","7.6","6.6"), '2', "3")))) 

X1979 %<>% mutate(stratum_area = ifelse(stratum == "1",29.61e+9,
                                        ifelse(stratum == "2", 21.91e+9,
                                               ifelse(stratum == "3", 93e+9, NA))))
X1979$stage <- as.factor(X1979$stage)
X1979 %<>% dplyr::select(-count)
#2021
# stations 4.8, 1.5, and 1.4 abandonned due to time
X2021_meta <- read_excel("data/larvae/New folder/Métadonnées_IML2021-14_final.xlsx", 
                         sheet = "Data", col_types = c("skip", 
                                                       "skip", "numeric", "date", "date", 
                                                       "numeric", "numeric", "skip", "skip", 
                                                       "text", "skip", "skip", "skip", "numeric", 
                                                       "numeric", "numeric", "skip", "numeric", 
                                                       "numeric", "numeric", "numeric", 
                                                       "text", "text"))
X2021_meta  %<>%  transmute(year = as.numeric(Année), 
                         month = as.numeric(str_sub(Date,6,7)),
                         day = as.numeric(str_sub(Date,9,10)),
                         date = as.Date(paste(day,month,year, sep = "/" ), format = "%d/%m/%Y"),
                         doy = as.numeric(yday(date)),
                         hour_start = str_sub(Heure,12,13),
                         minute_start = str_sub(Heure,15,16),
                         station = as.factor(Station),
                         consec = `Conséc.`,
                         bongo = as.factor(`Bab./Trib.`),
                         set_duration = `Durée totale (sec)`,
                         sample_depth = `Prof. bongo (m)`,
                         volume = `Volume (m3) Bionet`,
                         station_depth = `Prof. station (m)`, 
                         temperature = temp) %>% 
  mutate(hour_start = ifelse(hour_start == "00","24",hour_start), 
         hour_start = as.numeric(hour_start),
         minute_start = as.numeric(minute_start))

summary(X2021_meta)
X2021_meta[21,13] <-584.674301
X2021_meta[22,13] <-452.2408692
#stations 4.7 flowmeter volume should be used as lost contact with bionet

X2021_plankton <- read_excel("data/larvae/New folder/Métadonnées_IML2021-14_final.xlsx", 
                             sheet = "Plancton")
glimpse(X2021_plankton)

p1 <- X2021_plankton %>% transmute(year = as.numeric(Année), 
                                   month = as.numeric(str_sub(Date,6,7)),
                                   day = as.numeric(str_sub(Date,9,10)),
                                   date = as.Date(paste(day,month,year, sep = "/" ), format = "%d/%m/%Y"),
                                   doy = as.numeric(yday(date)),
                                   hour_start = str_sub(Heure,12,13),
                                   minute_start = str_sub(Heure,15,16),
                                   station = as.factor(Station),
                                   consec = `Conséc.`,
                                   bongo = as.factor(`Bâb./Trib.`), # notice the â
                                   fraction = `SOUS-ECHAN.  Œufs Maq.`,
                                   `Scomber scombrus (oeuf stade 1)` = `STADE 1`,
                                   `Scomber scombrus (oeuf stade 2)` = `STADE 2`,
                                   `Scomber scombrus (oeuf stade 3)` = `STADE 3`,
                                   `Scomber scombrus (oeuf stade 4)` = `STADE 4`,
                                   `Scomber scombrus (oeuf stade 5)` = `STADE 5`,
                                   total = `Total Oeufs    Maquereau`) %>% 
  pivot_longer(cols = 12:17, names_to = "stage", values_to = "count") %>% 
  mutate(hour_start = ifelse(hour_start == "00","24",hour_start), 
         hour_start = as.numeric(hour_start),
         minute_start = as.numeric(minute_start))
                             
summary(p1)  
glimpse(p1)

p2 <- X2021_plankton %>% transmute(year = as.numeric(Année), 
                                   month = as.numeric(str_sub(Date,6,7)),
                                   day = as.numeric(str_sub(Date,9,10)),
                                   date = as.Date(paste(day,month,year, sep = "/" ), format = "%d/%m/%Y"),
                                   doy = as.numeric(yday(date)),
                                   hour_start = str_sub(Heure,12,13),
                                   minute_start = str_sub(Heure,15,16),
                                   station = as.factor(Station),
                                   consec = `Conséc.`,
                                   bongo = as.factor(`Bâb./Trib.`), # notice the â
                                   fraction = `SOUS-ECHAN.      Larves`,
                                   `Scomber scombrus (larva)` = `S. scombrus`) %>% 
  pivot_longer(cols = 12, names_to = "stage", values_to = "count") %>% 
  mutate(hour_start = ifelse(hour_start == "00","24",hour_start), 
         hour_start = as.numeric(hour_start),
         minute_start = as.numeric(minute_start))
glimpse(p2)
summary(p1)
summary(p2)
# p1$consec<-as.numeric(as.character(p1$consec))
p2$consec <- round(p2$consec)
p1$consec <- round(p1$consec)
X2021_meta$consec <- round(X2021_meta$consec)
X2021p <- bind_rows(p1,p2)
glimpse(X2021p)
unique(as.numeric(X2021p$consec))
unique(X2021_meta$consec)

X2021 <- left_join(X2021p,X2021_meta)
X2021$stage <- as.factor(X2021$stage)

glimpse(X2021)
summary(X2021)
summary(X2021_meta)

X2021 %<>% 
  pivot_wider(names_from = "stage", values_from = "count") %>% 
  mutate(n1 = `Scomber scombrus (oeuf stade 1)`/fraction,
         n2 = `Scomber scombrus (oeuf stade 2)`/fraction,
         n3 = `Scomber scombrus (oeuf stade 3)`/fraction,
         n4 = `Scomber scombrus (oeuf stade 4)`/fraction,
         n5 = `Scomber scombrus (oeuf stade 5)`/fraction,
         nl = `Scomber scombrus (larva)`/fraction,
         n_15 = n1 + n5,
         total = n1 + n2 + n3 + n4 + n5) 
names(X2021)
X2021 %<>% pivot_longer(cols = 17:30, names_to = "stage", values_to = "n") %>% 
  mutate(n_m3 = n/volume, 
         n_m2 = n_m3*sample_depth, 
         trajet = "T1",
         month = as.numeric(month),
         day = as.numeric(day),
         hour_start = as.numeric(hour_start),
         minute_start = as.numeric(minute_start), 
         trajet = "T1")
coord <- egg_index %>% group_by(station) %>% dplyr::summarise(lon=mean(lon,na.rm=T),lat = mean(lat,na.rm=T))
X2021 <- left_join(X2021,coord)
summary(X2021)
glimpse(egg_index)
# add stations for mapping
df1<-data.frame(station = as.factor(c(1.4,1.5)),
                year = 2021,
                lon = -60.7562,
                lat = 48.1663)


X2021<- bind_rows(X2021,df1)

egg_index %<>% mutate(station_depth = as.numeric(station_depth))
# X2021$consec<-as.numeric(X2021$consec)
X1979$consec<-as.numeric(X1979$consec)
# egg_index$consec<-as.numeric(egg_index$consec)
# egg_index$consec<-as.numeric(egg_index$consec)
df_egg_index <- bind_rows(X1979,egg_index)  
glimpse(df_egg_index)
df_egg_index <- bind_rows(df_egg_index,X2021)
glimpse(df_egg_index)
summary(df_egg_index)

df_egg_index %<>% 
  mutate(stratum = as.factor(ifelse(station %in% 
                                      c("8.7", "7.7", "5.7","4.9","3.9","3.8","2.6","1.5","2.5","1.4","2.4","1.3","3.5","3.4",'3.3',"3.2",'3.1',"4.2","4.1","2.3","2.2",'2.1',"1.2","1.1"), "1", 
                                    ifelse(station %in% 
                                             c("10.1","9.4","9.3",'8.2','7.2','7.1',"6.3","6.2","6.1","5.2","5.1",'4.3',"4.4","4.5",'4.6',"3.6","3.7","4.8","5.6",'6.7',"6.8","7.5","7.6","6.6"), '2', "3")))) 

df_egg_index %<>% mutate(stratum_area = ifelse(stratum == "1",29.61e+9,
                                        ifelse(stratum == "2", 21.91e+9,
                                               ifelse(stratum == "3", 93e+9, NA))))
df1<-data.frame(station = as.factor(c(1.4,1.5, 4.8)),
                year = 2021,
                lon = c(-60.7562,-60.7562,-62.2524),
                lat = c(48.1663,48.1663,48.1672))
              
df_egg_index<- bind_rows(df_egg_index,df1)
df_egg_index %<>% mutate(presence = ifelse(n_m2>0,1,0))
df_egg_index %<>% mutate(trajet = ifelse(is.na(trajet),"T1",trajet))




# 1982 was not loaded. Take directly from Pelagiques/plancton/releves
# d82 <- read_excel("data/larvae/New folder/Relevé 1982/Métadonnées_1982.xls", 
#                   col_types = c("numeric", "numeric", "numeric", 
#                                 "numeric", "numeric", "numeric", 
#                                 "numeric", "numeric", "numeric", 
#                                 "numeric", "text", "numeric", "numeric", 
#                                 "numeric", "numeric", "numeric", 
#                                 "numeric", "numeric", "numeric", 
#                                 "numeric", "numeric", "numeric", 
#                                 "numeric", "numeric", "numeric", 
#                                 "numeric", "numeric", "numeric", 
#                                 "numeric", "numeric", "numeric", 
#                                 "numeric", "numeric", "numeric", 
#                                 "numeric", "numeric", "numeric", 
#                                 "numeric", "numeric", "numeric", 
#                                 "numeric", "numeric", "numeric", 
#                                 "numeric", "numeric", "numeric", 
#                                 "numeric", "numeric", "numeric", 
#                                 "numeric", "numeric", "numeric", 
#                                 "numeric", "numeric", "numeric"))
# glimpse(d82)
# summary(d82b)
# glimpse(df_egg_index)
# 
# d82 %<>% transmute(
#   year = ANNÉE,
#   date = DATE,
#   hour = str_sub(HEURE,5,6),
#   minute = str_sub(HEURE, 9,10),
#   trajet = TRAJET,
#   station = STATION,
#   strata = STRATE,
#   lat = LAT,
#   lon = LONG,
#   consec = CONS.,
#   volume = VOLUME,
#   sample_depth = PROF,
#   bongo = ENGIN, 
#   fraction = `SOUS-ECHAN.`,
#   `Scomber scombrus (oeuf stade 1)` = `STADE 1`,
#   `Scomber scombrus (oeuf stade 2)` = `STADE 2`,
#   `Scomber scombrus (oeuf stade 3)` = `STADE 3`,
#   `Scomber scombrus (oeuf stade 4)` = `STADE 4`,
#   `Scomber scombrus (oeuf stade 5)` = `STADE 5`,
#   `Scomber scombrus (larva)` = `Scomber scombrus`)
# 
# d82 %<>% mutate(bongo = ifelse(bongo=="BB", "B", "T"),
#                trajet = ifelse(trajet==1,"T1","T2"),
#                n1 = `Scomber scombrus (oeuf stade 1)`/fraction,
#                n2 = `Scomber scombrus (oeuf stade 2)`/fraction,
#                n3 = `Scomber scombrus (oeuf stade 3)`/fraction,
#                n4 = `Scomber scombrus (oeuf stade 4)`/fraction,
#                n5 = `Scomber scombrus (oeuf stade 5)`/fraction,
#                nl = `Scomber scombrus (larva)`/fraction,
#                n_15 = n1 + n5,
#                total = n1+n2+n3+n4+n5) 
# d82 %<>% 
# pivot_longer(cols = 15:28, names_to = "stage", values_to = "n") %>% 
#   mutate(n_m3 = n/volume, 
#          n_m2 = n_m3*sample_depth,
#          station = factor(station))
# d82 %<>% dplyr::select(-date)
# 
# df_egg_index <- bind_rows(df_egg_index, d82) %>% arrange(year, station)
# df_egg_index %<>% mutate(date = as.Date(paste(day,month,year, sep = "/" ), format = "%d/%m/%Y"))
#  
#####################  Save  ####################
# save(df_egg_index, file = "./Rdata/egg_index_2022_v2.Rdata")
save(df_egg_index, file = "./Rdata/df_egg_index_2022_v4.Rdata")
# load("./Rdata/df_egg_index_2022_v4.Rdata")
