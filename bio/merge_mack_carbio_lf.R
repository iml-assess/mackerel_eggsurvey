##' Merge mackerel carbio and mackerel lf R objects 
##' @param file_path directory
##' @details 
##' A function to merge the data from mackerel_lf.Rdata (raw data = sortie.dat) and cap_bio.Rdata (raw data = Maquereau_carbio.dat)
##' Must run read.bio.R and read_cap_lf.R to produce above .Rdata files
##' Also extra variables are added
##' Rewritten as a function by ANDREW SMITH, September 2019
##' @import lubridate readr tidyverse
##' @export
##' 
##' 
##' 
# source("./Rscripts/read_capelin_bio.R")
# source("./Rscripts/read.bio.R")
library(tidyverse)
library(lubridate)
library(magrittr)


### Load data

load("C:/Users/SmithAND/Documents/My_Stocks/Mackerel/Assessments/2022/mackerel_2022/data/rdata/2022/mackerel_bio_2022.Rdata")
load("C:/Users/SmithAND/Documents/My_Stocks/Mackerel/Assessments/2022/mackerel_2022/data/rdata/2022/mackerel_lf_2022.Rdata")

glimpse(bio)
glimpse(mackerel_lf)

# cap_bio$date = as_date(paste(cap_bio$year, cap_bio$month, cap_bio$day))

### Remove columns with only NAs
not_all_na <- function(x) any(!is.na(x))
mackerel_lf %<>% select_if(not_all_na)
bio %<>% select_if(not_all_na) 
unique(bio$sample)


mackerel_lf <- mackerel_lf[rep(1:nrow(mackerel_lf), mackerel_lf[["freq"]]), ]

cap_bio$sample <- str_pad(cap_bio$sample,2,side = "left",pad="0")
df<- bio %>% group_by(year,month,doy,subdivision,gear,sample) %>% dplyr::summarise(n=n()) %>% arrange(year,month,doy,subdivision,sample)
df2<-mackerel_lf %>% group_by(year,month,doy,subdivision,gear,sample) %>% dplyr::summarise(n=n())%>% arrange(year,month,doy,subdivision,sample)

df<-full_join(df,df2, by = c("year","month","doy", "subdivision","gear","sample")) %>% arrange(year,month,doy,subdivision,sample)

### Merge data
cap_lf$district <- as.character(cap_lf$district)
cap_bio$trip_num <- as.character(cap_bio$trip_num)

mackerel_lf$sample_type<-as.factor(mackerel_lf$sample_type)
bio$zone<-as.factor(bio$zone)
mackerel_lf$zone<-as.factor(mackerel_lf$zone)
mackerel_lf$trip_num<-as.factor(mackerel_lf$trip_num)
bio$trip_num<-as.factor(bio$trip_num)
mackerel_bio<-bind_rows(mackerel_lf, bio) 


####################  ADD NEW VARIABLES  ####################

mackerel_bio <- mackerel_bio %>%  mutate(f_age = fct_collapse(    
  factor(mackerel_bio$age),
  "10" = c("10","11","12","13","14","15","16","17","18","42","63"))) # the last two values (42 and 63) are typos and were validated

# This gives missing value an explicit factor level, ensuring that they appear in summaries and on plots.
mackerel_bio$f_age <- forcats::fct_explicit_na(mackerel_bio$f_age)

## QUARTER
mackerel_bio <- mackerel_bio %>% mutate(quarter = fct_collapse(
  factor(mackerel_bio$month),
  Q1 = c("1", "2","3"),
  Q2 = c("4", "5","6"),
  Q3 = c("7", "8", "9"),
  Q4 = c("10","11","12")))

# Nafo Subareas
# capelin_bio %<>%
#   mutate(NAFO = fct_collapse(Division,
#                              "4R" = c("4R", "4RA", "4RB", "4RC", "4RD","4R  ", "4RA ", "4RB ", "4RC ", "4RD "),
#                              "4S" = c("4S", "4SW", "4SV", "4SY", "4SZ","4S  ", "4SW ", "4SV ", "4SY ", "4SZ "),
#                              "4T" = c("4T", "4TM", "4TN", "4TO", "4TP", "4TQ","4T  ", "4TM ", "4TN ", "4TO ", "4TP ", "4TQ ")
#   ))

unique(mackerel_bio$subdivision)
unique(mackerel_bio$division)
mackerel_bio$subdivision <- trimws(mackerel_bio$subdivision)
mackerel_bio$division <- trimws(mackerel_bio$division)
mackerel_bio$gear <- trimws(mackerel_bio$gear)

## gear names
gear_key <-   
  list(
    CPO = "Dive_hand_tool",
    "CPO " = "Dive_hand_tool",
    CMA = "Dive_hand_tool",
    "CMA "  = "Dive_hand_tool",
    SSC = "Scottish_seine",
    "SSC "  = "Scottish_seine",
    ST2 = "Shrimp_trawl_no_grid",
    "ST2 " = "Shrimp_trawl_no_grid",
    GND = "Gillnet_drift",
    "GND " = "Gillnet_drift",
    GN = "Gillnet_na",
    "GN  " = "Gillnet_na",
    TXS = "Shrimp_trawl_na",
    "TXS " = "Shrimp_trawl_na",
    GRL = "Shrimp_trawl_na",
    "GRL " = "Shrimp_trawl_na",
    GRL2 = "Shrimp_trawl_rear_grid",
    GRL1 = "Shrimp_trawl_side_grid",
    LHM = "Mechanised_jigger",
    GNS = "Gillnet_set",
    "GNS " = "Gillnet_set",
    PS = "Purse_seine",
    "PS  " = "Purse_seine",
    FPN = "Trap_net",
    "FPN " = "Trap_net",
    LHP = "Handline",
    "LHP " = "Handline",
    "LHM " = "Mechanised_jigger",
    LX = "Line_na",
    NK = "NA",
    "NK  " = "NA",
    FWR = "Weir",
    "FWR " = "Weir",
    OTB = "Bottom_trawl",
    "OTB " = "Bottom_trawl",
    OTB2 = "Bottom_trawl",
    "OTB2 " = "Bottom_trawl",
    OTM = "Pelagic_trawl",
    "OTM " = "Pelagic_trawl",
    GN = "Gillnet_na",
    SB = "Beach_seine",
    "SB  " = "Beach_seine",
    SPR = "Pair_seine",
    LA = "Lampara_net",
    "LA  " = "Lampara_net",
    FIX = "Trap_na",
    "FIX " = "Trap_na",
    TS = "Tuck_seine",
    "TS  " = "Tuck_seine",
    MIS = "Misc.",
    "MIS " = "Misc.",
    SDN = "Danish_seine",
    "SDN " = "Danish_seine",
    LLS = "Longline",
    "LLS " = "Longline",
    LX = "Handline",
    "LX  " = "Handline",
    SPR = "Pair_seine",
    "SPR " = "Pair_seine",
    LMP = "Handline")
# Apply key
mackerel_bio %<>% mutate(gear_name = recode(gear, !!!gear_key))
unique(mackerel_bio$gear_name)
rm(gear_key)

# Simplified and aggregated gear names
mackerel_bio <- mackerel_bio %>% mutate(
  gear_group = fct_collapse(gear_name,
                            Nets_traps_weirs = c("Trap_net", "Weir", "Trap_na", "Lampara_net", "Beach_seine"),
                            Lines = c("Line_na", "Handline", "Mechanised_jigger", "Longline"),
                            Gillnets = c("Gillnet_set", "Gillnet_drift", "Gillnet_na"),
                            Other_seine = c("Pair_seine", "Danish_seine"),
                            `NA` = c("NA", "Misc."),
                            Bottom_trawl = c("Bottom_trawl", "Shrimp_trawl_no_grid", "Shrimp_trawl_na"),
                            Pelagic_trawl = c("Pelagic_trawl"),
                            Purse_seine = c("Purse_seine"),
                            Tuck_seine = c("Tuck_seine")
  )
)




# all lower case levels
mackerel_bio <- as_tibble(lapply(mackerel_bio,      # Convert data with tolower function
                                function(variables) {
                                  if (is.factor(variables)) {
                                    return(tolower(variables))
                                  } else {
                                    return(variables)
                                  }
                                }),
                         stringsAsFactors = F)

names(mackerel_bio) <- tolower(names(mackerel_bio))
summary(mackerel_bio)
# save(capelin_bio, file = "./Rdata/bio/capelin_bio.Rdata")
save(mackerel_bio, file = "C:/Users/SmithAND/Documents/My_Stocks/Mackerel/Assessments/2022/mackerel_2022/data/rdata/2022/mackerel_bio_lf_2022.Rdata")
