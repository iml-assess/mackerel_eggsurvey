#' Read mackerel length frequency data
#'
#' @param file_path Where you copy pasted Sortie.dat 
##' A function to read sortie.dat flat files extracted through the Peche SAS program. 
##' Carbio extracted from SAS program autoexecpeche found in "S:/SAS/Peche/autoexec_peche_64bits.sas"
##' Access SAS 9.4 through Citrix connection: https://citrix.dfo-mpo.gc.ca/vpn/index.html username = ent , password = email
##' Choose species, carbio option, f8 to run the script, save the .dat file in I:/ and then copy to new location
##' Provided by DAISS (contact sylvain.hurtubuise)
##' Variable name descriptions in "\\dcqcimlna01a\BD_Peches\EchComm\Documentation\Poissons\Descriptions_tables_Capelan+Hareng+Maquereau.xls"
##' Local copy: C:\Users\SmithAND\Documents\Data\Samples_fish\Copy of Descriptions_tables_Capelan+Hareng+Maquereau.xls
##' By ANDREW SMITH, July 2019 
##' Rewritten as a function by ANDREW SMITH, September 2019
##' Updated December-January 2021, 
##' Updated 14/04/2022
##' ##' Updated 27/04/2022
##' Update 28/04/2022: wrapped a function around catchR::read.lf. It produces same dims as old versions and is more complete however i don't like the column names and will add extra variables
#'  29/07/2022 after having standarized capelin lf an carbio I went digging in read.lf and read.carbio. I had already noticed that fwf data differed between mackerel and capelin and then I read that catchR::read.lf uses a specific file for column widths vs what I had inherited in sas and R code. How did EVB do the mackerel example for caa? this will have to be validated. capelin cleaning of carbio and lf datasets revealed that district, no_voy, s_ech, etc are not equal among carbio and sortie.dat datasets. sample ie s_ech in carbio is a concatenated version of zone and no_voy in the lf dataset. to be seen if same in mackerel. fucking data archeology. 
#'  15/08/2022 update: integrating the changes that I made/discovered for capelin into mackerel 'read' functions. should make it easier to match carbio and lf data samples 
#' @return
#' @export
#' 
#' @examples file_path = "./data/bio/Sortie.dat"; lf <- read_mackerel_lf(file_path = file_path)
#' file_path = "C:/Users/SmithAND/Documents/My_Stocks/Mackerel/Assessments/2022/mackerel_2022/data/bio/sortie.dat"

read_mackerel_lf <- function(file_path, ...) {
  require(lubridate)
  require(tidyverse)
  require(readr)
  require(magrittr)
  require(forcats)
  require(catchR)
  require(stringi)
  
  lf <- catchR::read.lf(file = file_path)
 
  # reorder columns
  col_order <- c("year", "month", "day","date","trim","doy",
                 "nafo", "nafo.sub", "province",
                 "district","lat","lon", "species", "nbpc", "activity",
                 "trip", "gear.cat", "gear.number", "mesh","depth",
                 "sample.id", "freq", "sample.cat", "sample.ncat", "id.par",
                 "state.id", "weight.sample", "weight.land", "sex",
                 "length.bin", "length.type", "length", "n", "id")
  
  lf <- lf[, col_order]
  
  # remove columns where all NA
  names(lf)
  lf <- Filter(function(x)!all(is.na(x)), lf) # gear.number and mesh removed
  names(lf)
  
  # remove variables 
  
  ## rename variables for conveniance
  lf <-
    lf %>% rename(
      quarter = trim,
      division = nafo,
      subdivision = nafo.sub,
      type = activity,
      trip_num = trip,
      gear_code = gear.cat,
      sample_id = sample.id,
      sample_cat = sample.cat,
      sample_ncat = sample.ncat,
      id_person = id.par,
      fish_state = state.id,
      mass_sample = weight.sample, 
      mass_landed = weight.land,
      length_bin = length.bin,
      length_type = length.type,
      length_raw = length,
      freq = n,
      key = id)
  
  names(lf)
  # remove variables
  lf <- lf %>% dplyr::select(-district,-species,-fish_state,-key)
  
  # convert length types to fork length
  table(lf$length_type)
  # in 2021:
  # fork_length headless_length standard_length    total_length 
  # 58361              22              24             158 
  
  # Hansen et al., 2018 https://doi.org/10.1016/j.fishres.2018.04.002
  ## for northeast Atlantic
  ## length_new <- alpha(i.e. slope) * length + Beta(i.e. intercept)+ error
  ## TL - > FL: alpha = 0.921, beta = 4.280, r2 = 0.987, sd = 4.907 
  
  # Hunt, JJ and Stobo, WT 1976   ICNAF Res.Doc. 76/XII/136
  ## from 4TVWX in 1975
  ## TL - > FL: alpha = 0.9270, beta = -0.0382, r2 = 0.97
  
  # MacKay 1976 - his PHD
  ## FL = 0.916(TL)

  ## no conversion of SL -> FL available or Headless -> FL
  
  df<-lf %>% mutate(length_h = ifelse(length_type == "total_length", 4.28 + (length_raw*0.921), length_raw),
                length_hs = ifelse(length_type == "total_length",-0.0382 + (length_raw*0.927),length_raw),
                length_m = ifelse(length_type == "total_length",0.916*length_raw,length_raw)) 
  df %<>% dplyr::filter(length_type == "total_length") %>% dplyr::select(sample_id, length_raw, length_h, length_hs, length_m) 
    
  cor.test(df$length_h, df$length_hs) # equivalent
  cor.test(df$length_h, df$length_m) # equivalent
  # see distribution
  df %>% pivot_longer(cols = 2:5, names_to = "type", values_to = "length") %>% 
    ggplot(aes(length, fill = type)) + geom_density(alpha = 0.5)
  rm(df)
  
  # Choose HS because meh
  lf %<>% mutate(length = ifelse(length_type == "total_length",-0.0382 + (length_raw*0.927),length_raw))
  
  # can't do anything with headless or standard length so omit
  lf %<>% dplyr::filter(length_type %in% c("total_length","fork_length"))
  
  # add variables
  
  ## quarter
  lf %<>% mutate(quarter = fct_collapse(
    factor(lf$month),
    Q1 = c("1", "2","3"),
    Q2 = c("4", "5","6"),
    Q3 = c("7", "8", "9"),
    Q4 = c("10","11","12")))
  
  lf <- lf %>% mutate(bioregion = fct_collapse(    
    factor(lf$division),
    "SS" = c("4V","4VS","4VN","4W","4X","5Y","5Z","5ZE"),
    "NLS" = c("3K","3L","3P","2J"),
    "GSL" = c("4R","4S","4T"))) %>% droplevels()
  
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
  lf %<>% mutate(gear_name = recode(gear_code, !!!gear_key))
  unique(lf$gear_name)
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
 
  
  # lat lon
  unique(mackerel_bio$lon)
 
  
  # load("~/My_Stocks/Mackerel/Assessments/2022/mackerel_2022/data/rdata/NAFO_centroids.Rdata")
  # 
  # # it is a mix of NAFO Divisions and subdivisions/unit areas
  # unique(lf$subdivision)
  # lf$subdivision <- trimws(lf$subdivision)
  # 
  # # rename variables
  # NAFO_centroids$subdivision <- NAFO_centroids$UnitArea
  # NAFO_centroids %<>% transmute(subdivision = subdivision, lon_subdiv = X_Longitude, lat_subdiv = Y_Latitude)
  # 
  # # merge subdivision coordinates
  # lf <- left_join(lf, NAFO_centroids, by = "subdivision")
  # 
  # # now same for Divisions
  # NAFO_centroids$division <- str_sub(NAFO_centroids$subdivision,1,2)
  # NAFO_centroids %<>% group_by(division) %>% 
  #   dplyr::summarise(lon_div = mean(lon_subdiv), lat_div = mean(lat_subdiv))
  # 
  # # merge
  # lf <- left_join(lf, NAFO_centroids,  by = "division")
  # lf %<>% mutate(lat_nafo = ifelse(is.na(lat_subdiv), lat_div, lat_subdiv), lon_nafo = ifelse(is.na(lon_subdiv), lon_div, lon_subdiv))
  # 
  # sample type
  sample_type_key <-   
    list(
      "99" = "commercial_port",
      "98" = "commercial_sea",
      "95" = "indicator_fisher",
      "97" = "observer_at_sea",
      "96" = "research_acoustic",
      "93" = "research",
      "94" = "research_genetics")
  # Apply key
  mackerel_bio %<>% mutate(sample_type = recode(type, !!!sample_type_key))
  unique(mackerel_bio$sample_type)
  
  # EXPAND AND FILL ROWS OF OBS BASED ON FREQUENCY COUNTS
  # eg if line 1 reads freq = 5 then this line is expanded to five lines. Be careful when summarising/aggregating data
  
  
  
}
