#' Read mackerel carbio data
#'
#' @param file_path Where you copy pasted Maquereau_carbio.dat 
##' A function to read sortie.dat flat files extracted through the Peche SAS program. 
##' Carbio extracted from SAS program autoexecpeche found in "S:/SAS/Peche/autoexec_peche_64bits.sas"
##' Access SAS 9.4 through Citrix connection: https://citrix.dfo-mpo.gc.ca/vpn/index.html username = ent , password = email
##' Choose species, carbio option, f8 to run the script, save the .dat file in I:/ and then copy to new location
##' Provided by DAISS (contact sylvain.hurtubuise)
##' bY ANDREW SMITH, July 2019 
##' Rewritten as a function by ANDREW SMITH, September 2019
##' Updated December-January 2021, 
##' Updated 14/04/2022
##' Compared output to catchR::read.bio() and dims are different. Cause is read.bio which duplicates all entries from nafo 4TO and assigns them 4TQ
##' maquereau_carbio contains 138627 lines and catchR produces 140184.. also some nafo.sub are missing e.g. 3LB, 4TP, 4TQ, 4SI
##' Variable names and metadata are found in "\\dcqcimlna01a\BD_Peches\EchComm\Documentation\Poissons\Descriptions_tables_Capelan+Hareng+Maquereau.xls"
#' 
#' @return
#' @export
#' 
#' @examples file_path = "./data/bio/Maquereau_carbio.dat"; bio <- read_mackerel_bio(file_path = file_path)
#' file_path = "C:/Users/SmithAND/Documents/My_Stocks/Mackerel/Assessments/2022/mackerel_2022/data/bio/Maquereau_carbio.dat"

read_mackerel_bio <- function(file_path, ...) {
  require(lubridate)
  require(tidyverse)
  require(readr)
  require(magrittr)
  require(forcats)
  
  bio <-
    read_fwf(
      file_path,
      fwf_cols(
        year = c(7, 10),
        month = c(11, 12),
        day = c(13, 14),
        type = c(23, 24),
        division = c(83, 84),
        subdivision = c(78, 81),
        gear = c(25, 28),
        sample = c(55, 60), # s_ech
        fish_id = c(29, 32), # no_prel
        length = c(33, 35),
        mass = c(36, 41),
        sex = c(42),
        maturity_stage = c(43),
        mass_gonad = c(44, 49),
        age = c(53, 54), # age_f = age - the final age agreed upon
        carbio_id = c(61, 66), # s_cbiol = unique_id
        general_comment = c(68, 73), # s_rcbioi = code for comment on individual carbio
        age_1 = c(75, 76), # first age reading
        zone = c(1, 3), # 4TO and 4TQ are both coded 439. Beware!
        trip_num = c(4, 6)
        )
    ) %>% 
    arrange(year, month, division)
  glimpse(bio)
  
  # change variable class
  bio %<>% 
    mutate(day = as.numeric(day),
           month = as.numeric(month), 
           type = as.factor(type), 
           division = as.factor(division),
           subdivision = as.factor(subdivision),
           gear = as.factor(gear),
           sample = as.factor(sample),
           # sample_num = as.factor(sample_num),
           sex = as.factor(sex),
           maturity_stage = as.factor(maturity_stage)) 

  
  
  # new variables
  bio %<>%
    mutate(
      source = "bio",
      quarter = fct_collapse(
        factor(bio$month),
        Q1 = c("1", "2", "3"),
        Q2 = c("4", "5", "6"),
        Q3 = c("7", "8", "9"),
        Q4 = c("10", "11", "12")
      ),
      age_f = as.factor(age),
      year_f = as.factor(year),
      date = lubridate::ymd(paste(year, month, day)),
      doy = as.numeric(format(date, "%j")),
      sex = fct_explicit_na(sex, na_level = "NA"),
      gsi = mass_gonad / mass * 100,
      maturity = recode(
        maturity_stage,
        '1' = "1_immature_virgin",
        '2' = "2_immature",
        '3' = "3_maturing",
        '4' = "4_maturing",
        '5' = "5_ripe",
        '6' = "6_spawning",
        '7' = "7_spent",
        '8' = "8_recovering"
      ),
      maturité  = recode(
        maturity,
        '1_immature_virgin' = "_immature_vierge",
        '2_immature' = "2_immature",
        '3_maturing' = "3_maturation",
        '4_maturing' = "4_maturation",
        '5_ripe' = "5_mature",
        '6_spawning' = "6_fraie",
        '7_spent' = "7_fin_de_ponte",
        '8_recovering' = "8_récupération"
      )
    )
  
  # Sample type 
  # Note that type 93 is a catch-all i.e. it could be from research vessel surveys or special projects
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
  bio %<>% mutate(sample_type = recode(type, !!!sample_type_key))
  unique(bio$sample_type)
  
  save(bio, file = "C:/Users/SmithAND/Documents/My_Stocks/Mackerel/Assessments/2022/mackerel_2022/data/rdata/2022/mackerel_bio_2022.Rdata")
  
  # gear names
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
  bio %<>% mutate(gear_name = recode(gear, !!!gear_key))
  unique(bio$gear_name)
  unique(bbio$engin_en)
  # Simplified and aggregated gear names
  bio <- bio %>% mutate(
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
  
  bio <- bio %>%  mutate(f_age = fct_collapse(    
    factor(bio$age),
    "10" = c("10","11","12","13","14","15","16","17","18","42","63"))) # the last two values (42 and 63) are typos and were validated
  
  # This gives missing value an explicit factor level, ensuring that they appear in summaries and on plots.
  bio$f_age <- forcats::fct_explicit_na(bio$f_age)
  
  # group some divisions for convenience
  unique(bio$division)
  bio <- bio %>% mutate(bioregion = fct_collapse(    
    factor(bio$division),
    "SS" = c("4V", "4W","4X","5Y","5Z"),
    "NLS" = c("3K","3L","3P","2J"),
    "GSL" = c("4R","4S","4T")))
  
  load("~/My_Stocks/Mackerel/Assessments/2022/mackerel_2022/data/rdata/NAFO_centroids.Rdata")
  
  # it is a mix of NAFO Divisions and subdivisions/unit areas
  unique(bio$subdivision)
  bio$subdivision <- trimws(bio$subdivision)
  
  # rename variables
  NAFO_centroids$subdivision <- NAFO_centroids$UnitArea
  NAFO_centroids %<>% transmute(subdivision = subdivision, lon_subdiv = X_Longitude, lat_subdiv = Y_Latitude)
  
  # merge subdivision coordinates
  bio <- left_join(bio, NAFO_centroids, by = "subdivision")
  
  # now same for Divisions
  NAFO_centroids$division <- str_sub(NAFO_centroids$subdivision,1,2)
  NAFO_centroids %<>% group_by(division) %>% 
    dplyr::summarise(lon_div = mean(lon_subdiv), lat_div = mean(lat_subdiv))
  
  # merge
  bio <- left_join(bio, NAFO_centroids,  by = "division")
  bio %<>% mutate(lat_nafo = ifelse(is.na(lat_subdiv), lat_div, lat_subdiv), lon_nafo = ifelse(is.na(lon_subdiv), lon_div, lon_subdiv))
  
  # reorder columns
  col_order <- c("source", "year", "year_f", 
                 "type", "sample_type",
                 "sample_id", "sample_num",
                 "quarter", "date", "month", 
                 "day", "doy", "division",
                 "subdivision", "bioregion",
                 "gear_code","gear_name", "gear_group", 
                  "length", "mass", "mass_gonad",
                 "gsi", "maturity", "maturité",
                 "sex", "age", "age_f", "f_age",
                 "lon_subdiv", "lat_subdiv", "lon_div",
                 "lat_div", "lat_nafo", "lon_nafo")
  bio <- bio[, col_order]
  return(bio)
  bio
}

# all lower case levels
bio <- as_tibble(lapply(bio,      # Convert data with tolower function
                                function(variables) {
                                  if (is.factor(variables)) {
                                    return(tolower(variables))
                                  } else {
                                    return(variables)
                                  }
                                }),
                         stringsAsFactors = F)

names(bio) <- tolower(names(bio))
summary(bio)
