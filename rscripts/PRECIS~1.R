# Simulation study on precision and accuracy of mackerel catch at age inspired from cod sentinel survey analysis

# Goal ... see if reductions of otolith samples results in different ALKs and 2) see if reduction of length frequencies influences catch at age

#####################  Load libraries  ####################

packages = c('FSA','nnet','broom', 'purrr', 'nnet', 'readxl', 'magrittr','tidyverse')
invisible(lapply(packages, function(x) {if (!require(x, character.only = T)) {install.packages(x);require(x)}}))

#####################  Load data  ####################

source("./rscripts/read_mackerel_lf.R")
source("~/My_Stocks/Mackerel/Assessments/2022/mackerel_2022/rscripts/read_mackerel_bio.R", encoding = 'UTF-8', echo=TRUE)
lf <- read_mackerel_lf("./data/bio/sortie.dat")
bio <- read_mackerel_bio("./data/bio/Maquereau_carbio.dat")

# Length frequencies
# aggregate data into strata and then randomly thin a certain percentage of the rows without replacement
# see/test whether mean lengths and CVs vary much with thinning (10,20,30,40,50, and 75%). Assuming population mean is 100% of dataset (minus non commercial samples and not Q1 as too sparse, also exclude n < 99 for strata)

# reference
lf_ref <- lf %>% 
  dplyr::filter(quarter != "Q1", sample_type == "commercial_port") %>% 
  group_by(year, quarter, bioregion, gear_group) %>% 
  dplyr::summarise(n = n(), mean_length = mean(length, na.rm = T), sd_length = sd(length,na.rm=T), cv_length = (sd_length/mean_length)*100) %>% 
  ungroup() %>% 
  mutate(perc = "100")
summary(lf_ref)
# filter to where n => 100 obs
lf_ref %<>% dplyr::filter(n>99) # 2916 obs
  
# thin to 10,20,30,40,50, and 75% of full data and compare
  lf_10 <- lf %>% 
    dplyr::filter(quarter != "Q1", sample_type == "commercial_port") %>% 
    group_by(year, quarter, bioregion, gear_group) %>% 
    mutate(n = n()) %>% 
    dplyr::filter(n>99) %>% 
    slice_sample(prop = .1, replace = F) %>% 
    dplyr::summarise(n = n(), mean_length = mean(length, na.rm = T), sd_length = sd(length,na.rm=T), cv_length = (sd_length/mean_length)*100) %>% 
    ungroup() %>% 
    mutate(perc = "10")
  
  lf_20 <- lf %>% 
    dplyr::filter(quarter != "Q1", sample_type == "commercial_port") %>% 
    group_by(year, quarter, bioregion, gear_group) %>% 
    mutate(n = n()) %>% 
    dplyr::filter(n>99) %>% 
    slice_sample(prop = .2, replace = F) %>% 
    dplyr::summarise(n = n(), mean_length = mean(length, na.rm = T), sd_length = sd(length,na.rm=T), cv_length = (sd_length/mean_length)*100) %>% 
    ungroup() %>% 
    mutate(perc = "20")
  
  lf_30 <- lf %>% 
    dplyr::filter(quarter != "Q1", sample_type == "commercial_port") %>% 
    group_by(year, quarter, bioregion, gear_group) %>% 
    mutate(n = n()) %>% 
    dplyr::filter(n>99) %>% 
    slice_sample(prop = .3, replace = F) %>% 
    dplyr::summarise(n = n(), mean_length = mean(length, na.rm = T), sd_length = sd(length,na.rm=T), cv_length = (sd_length/mean_length)*100) %>% 
    ungroup() %>% 
    mutate(perc = "30")
  
  lf_40 <- lf %>% 
    dplyr::filter(quarter != "Q1", sample_type == "commercial_port") %>% 
    group_by(year, quarter, bioregion, gear_group) %>% 
    mutate(n = n()) %>% 
    dplyr::filter(n>99) %>% 
    slice_sample(prop = .4, replace = F) %>% 
    dplyr::summarise(n = n(), mean_length = mean(length, na.rm = T), sd_length = sd(length,na.rm=T), cv_length = (sd_length/mean_length)*100) %>% 
    ungroup() %>% 
    mutate(perc = "40")
  
  lf_50 <- lf %>% 
    dplyr::filter(quarter != "Q1", sample_type == "commercial_port") %>% 
    group_by(year, quarter, bioregion, gear_group) %>% 
    mutate(n = n()) %>% 
    dplyr::filter(n>99) %>% 
    slice_sample(prop = .5, replace = F) %>% 
    dplyr::summarise(n = n(), mean_length = mean(length, na.rm = T), sd_length = sd(length,na.rm=T), cv_length = (sd_length/mean_length)*100) %>% 
    ungroup() %>% 
    mutate(perc = "50")
  
  lf_75 <- lf %>% 
    dplyr::filter(quarter != "Q1", sample_type == "commercial_port") %>% 
    group_by(year, quarter, bioregion, gear_group) %>% 
    mutate(n = n()) %>% 
    dplyr::filter(n>99) %>% 
    slice_sample(prop = .75, replace = F) %>% 
    dplyr::summarise(n = n(), mean_length = mean(length, na.rm = T), sd_length = sd(length,na.rm=T), cv_length = (sd_length/mean_length)*100) %>% 
    ungroup() %>% 
    mutate(perc = "75")
# merge
  df_lf <- bind_rows(lf_ref, lf_10)
  df_lf <- bind_rows(df_lf, lf_20)
  df_lf <- bind_rows(df_lf, lf_30)
  df_lf <- bind_rows(df_lf, lf_40)
  df_lf <- bind_rows(df_lf, lf_50)
  df_lf <- bind_rows(df_lf, lf_75)
# reorder levels for plotting  
  df_lf$perc <- factor(df_lf$perc, levels = c("100","10","20","30","40","50","75"))
  
  
# model it?
  fit_mnull <- glm(data = df_lf, mean_length ~ 1)
  fitm <- glm(data = df_lf, mean_length ~ perc*quarter*gear_group*bioregion)
  fitm2 <- glm(data = df_lf, mean_length ~ perc*quarter)
  fitm3 <- glm(data = df_lf, mean_length ~ perc*quarter*gear_group)
  fitm4 <- glm(data = df_lf, mean_length ~ perc*quarter*bioregion)

# compare 
  compare_performance(reference = fit_mnull, fitm, fitm2, fitm3, fitm4)
  
  
  
  # figure of means
  df_lf %>% ggplot(aes(year, mean_length, colour = perc)) +
    stat_summary(fun.data = "mean_cl_boot", position = position_dodge(width = 2)) +
    theme_minimal() +
    scale_colour_viridis_d(end = 0.8, option = "B")
  
  # figure of sd
  df_lf %>% ggplot(aes(year, sd_length, colour = perc)) +
    stat_summary(fun.data = "mean_cl_boot", position = position_dodge(width = 2)) +
    theme_minimal() +
    scale_colour_viridis_d(end = 0.8, option = "B")
  
  # figure of cv
  df_lf %>% ggplot(aes(year, cv_length, colour = perc)) +
    stat_summary(fun.data = "mean_cl_boot", position = position_dodge(width = 2)) +
    theme_minimal() +
    scale_colour_viridis_d(end = 0.8, option = "B")
  
  # figure of cv
  df_lf %>% ggplot(aes(year, cv_length, colour = perc)) +
    stat_summary(fun.data = "mean_cl_boot", position = position_dodge(width = 2)) +
    theme_minimal() +
    scale_colour_viridis_d(end = 0.8, option = "B") +
    facet_wrap(vars(quarter))
  
  # figure of cv
  df_lf %>% ggplot(aes(year, cv_length, colour = perc)) +
    stat_summary(fun.data = "mean_cl_boot", position = position_dodge(width = 2)) +
    theme_minimal() +
    scale_colour_viridis_d(end = 0.8, option = "B") +
    facet_wrap(vars(bioregion))
  
  df_lf %>% ggplot(aes(year, cv_length, colour = perc)) +
    stat_summary(fun.data = "mean_cl_boot", position = position_dodge(width = 2)) +
    theme_minimal() +
    scale_colour_viridis_d(end = 0.8, option = "B") +
    facet_wrap(vars(gear_group))
  
  df_lf %>% ggplot(aes(quarter, cv_length, colour = perc)) +
    stat_summary(fun.data = "mean_cl_boot", position = position_dodge(width = 0.2)) +
    theme_minimal() +
    scale_colour_viridis_d(end = 0.8, option = "B") 
  
  df_lf %>% ggplot(aes(bioregion, cv_length, colour = perc)) +
    stat_summary(fun.data = "mean_cl_boot", position = position_dodge(width = 0.2)) +
    theme_minimal() +
    scale_colour_viridis_d(end = 0.8, option = "B") 
  
  
  ### pseudo accuracy as per Stamatopoulos 1999. Observations on the geometrical properties of accuracy growth in sampling with finite populations. FAO tech paper 338 vii 39p
   # A = 1-((sample mean - population mean(assumed to be 100))/(ymax-ymin assumed to be pop mean))      
  
  df_lf %<>% 
    dplyr::select(-sd_length, -cv_length, -n) %>% 
    pivot_wider(values_from = mean_length, names_from = perc) %>% 
    mutate(acc_10 = (1 - abs(`10`-`100`)/`100`),
           acc_20 = (1 - abs(`20`-`100`)/`100`),
           acc_30 = (1 - abs(`30`-`100`)/`100`),
           acc_40 = (1 - abs(`40`-`100`)/`100`),
           acc_50 = (1 - abs(`50`-`100`)/`100`),
           acc_75 = (1 - abs(`75`-`100`)/`100`))  
  
  df_lf %<>% select(year, quarter, bioregion, gear_group, acc_10, acc_20, acc_30,acc_40, acc_50, acc_75) %>% 
    pivot_longer(cols = 5:10, names_to = 'perc', values_to = 'accuracy')

  df_lf %>% ggplot(aes(year, accuracy)) + stat_summary(fun.data = "mean_cl_boot")  
  df_lf %>% ggplot(aes(quarter, accuracy)) + stat_summary(fun.data = "mean_cl_boot")  
  df_lf %>% ggplot(aes(bioregion, accuracy)) + stat_summary(fun.data = "mean_cl_boot")  
  df_lf %>% ggplot(aes(perc, accuracy, colour = gear_group)) + stat_summary(fun.data = "mean_cl_boot")  
  
# Conclusion. Results very robust and accurate despite thinning. Would not go down to 10 % for many reasons but even 30% of the data would be fine.   
  
# Now do same but for mean lengths by age and also catch at age
  # Stick with major gear types
  
  # reference
  bio_ref <- bio %>% 
    dplyr::filter(quarter != "Q1", sample_type == "commercial_port", gear_group %in% c("Nets_traps_weirs", "Gillnets", "Purse_seine", "Tuck_seine", "Lines"), f_age != "(Missing)") %>% 
    group_by(year, quarter, bioregion, gear_group, f_age) %>% 
    dplyr::summarise(n = n(), mean_length = mean(length, na.rm = T), sd_length = sd(length,na.rm=T), cv_length = (sd_length/mean_length)*100) %>% 
    ungroup() %>% 
    mutate(perc = "100") %>% droplevels()
  summary(bio_ref)
  unique(bio$gear_group)
  
  
  # thin to 10,20,30,40,50, and 75% of full data and compare
  bio_10 <- bio %>% 
    dplyr::filter(quarter != "Q1", sample_type == "commercial_port", gear_group %in% c("Nets_traps_weirs", "Gillnets", "Purse_seine", "Tuck_seine", "Lines"), f_age != "(Missing)") %>% 
    group_by(year, quarter, bioregion, gear_group) %>% 
    slice_sample(prop = .1, replace = F) %>% 
    ungroup() %>% 
    group_by(year, quarter, bioregion, gear_group, f_age) %>% 
    dplyr::summarise(n = n(), mean_length = mean(length, na.rm = T), sd_length = sd(length,na.rm=T), cv_length = (sd_length/mean_length)*100) %>% 
    ungroup() %>% 
    mutate(perc = "10")
  
  bio_20 <- bio %>% 
    dplyr::filter(quarter != "Q1", sample_type == "commercial_port", gear_group %in% c("Nets_traps_weirs", "Gillnets", "Purse_seine", "Tuck_seine", "Lines"), f_age != "(Missing)") %>% 
    group_by(year, quarter, bioregion, gear_group) %>% 
    slice_sample(prop = .2, replace = F) %>% 
    ungroup() %>% 
    group_by(year, quarter, bioregion, gear_group, f_age) %>% 
    dplyr::summarise(n = n(), mean_length = mean(length, na.rm = T), sd_length = sd(length,na.rm=T), cv_length = (sd_length/mean_length)*100) %>% 
    ungroup() %>%  
    mutate(perc = "20")
  
  bio_30 <- bio %>% 
    dplyr::filter(quarter != "Q1", sample_type == "commercial_port", gear_group %in% c("Nets_traps_weirs", "Gillnets", "Purse_seine", "Tuck_seine", "Lines"), f_age != "(Missing)") %>% 
    group_by(year, quarter, bioregion, gear_group) %>% 
    slice_sample(prop = .3, replace = F) %>% 
    ungroup() %>% 
    group_by(year, quarter, bioregion, gear_group, f_age) %>% 
    dplyr::summarise(n = n(), mean_length = mean(length, na.rm = T), sd_length = sd(length,na.rm=T), cv_length = (sd_length/mean_length)*100) %>% 
    ungroup() %>% 
    mutate(perc = "30")
  
  bio_40 <- bio %>% 
    dplyr::filter(quarter != "Q1", sample_type == "commercial_port", gear_group %in% c("Nets_traps_weirs", "Gillnets", "Purse_seine", "Tuck_seine", "Lines"), f_age != "(Missing)") %>% 
    group_by(year, quarter, bioregion, gear_group) %>% 
    slice_sample(prop = .4, replace = F) %>% 
    ungroup() %>% 
    group_by(year, quarter, bioregion, gear_group, f_age) %>% 
    dplyr::summarise(n = n(), mean_length = mean(length, na.rm = T), sd_length = sd(length,na.rm=T), cv_length = (sd_length/mean_length)*100) %>% 
    ungroup() %>% 
    mutate(perc = "40")
  
  bio_50 <- bio %>% 
    dplyr::filter(quarter != "Q1", sample_type == "commercial_port", gear_group %in% c("Nets_traps_weirs", "Gillnets", "Purse_seine", "Tuck_seine", "Lines"), f_age != "(Missing)") %>% 
    group_by(year, quarter, bioregion, gear_group) %>% 
    slice_sample(prop = .5, replace = F) %>% 
    ungroup() %>% 
    group_by(year, quarter, bioregion, gear_group, f_age) %>% 
    dplyr::summarise(n = n(), mean_length = mean(length, na.rm = T), sd_length = sd(length,na.rm=T), cv_length = (sd_length/mean_length)*100) %>% 
    ungroup() %>% 
    mutate(perc = "50")
  
  bio_75 <- bio %>% 
    dplyr::filter(quarter != "Q1", sample_type == "commercial_port", gear_group %in% c("Nets_traps_weirs", "Gillnets", "Purse_seine", "Tuck_seine", "Lines"), f_age != "(Missing)") %>% 
    group_by(year, quarter, bioregion, gear_group) %>% 
    slice_sample(prop = .75, replace = F) %>% 
    ungroup() %>% 
    group_by(year, quarter, bioregion, gear_group, f_age) %>% 
    dplyr::summarise(n = n(), mean_length = mean(length, na.rm = T), sd_length = sd(length,na.rm=T), cv_length = (sd_length/mean_length)*100) %>% 
    ungroup() %>% 
    mutate(perc = "75")
  
  # merge
  df_bio <- bind_rows(bio_ref, bio_10)
  df_bio <- bind_rows(df_bio, bio_20)
  df_bio <- bind_rows(df_bio, bio_30)
  df_bio <- bind_rows(df_bio, bio_40)
  df_bio <- bind_rows(df_bio, bio_50)
  df_bio <- bind_rows(df_bio, bio_75)
  
  df_bio$perc <- factor(df_bio$perc, levels = c("100","10","20","30","40","50","75"))
  
  
  # model it?
  fit_mnull <- glm(data = df_bio, mean_length ~ 1)
  fitm <- glm(data = df_bio, mean_length ~ f_age * perc)
  fitm2 <- glm(data = df_bio, mean_length ~  f_age * perc + quarter)
  fitm3 <- glm(data = df_bio, mean_length ~  f_age * perc + gear_group)
  fitm4 <- glm(data = df_bio, mean_length ~  f_age * perc + bioregion)
  
  # compare
  compare_performance(fit_mnull, fitm,fitm2, fitm3, fitm4)
  summary(anova(fitm))
  anova(fit_mnull, fitm)
  
  
  # figure of means
  df_bio %>% ggplot(aes(year, mean_length, colour = perc)) +
    stat_summary(fun.data = "mean_cl_boot", position = position_dodge(width = 2)) +
    theme_minimal() +
    scale_colour_viridis_d(end = 0.8, option = "B") + facet_wrap(vars(f_age))
  
  # figure of means
  df_bio %>% ggplot(aes(f_age, mean_length, colour = perc)) +
    stat_summary(fun.data = "mean_cl_boot", position = position_dodge(width = 0.5)) +
    theme_minimal() +
    scale_colour_viridis_d(end = 0.8, option = "B") 
  
  # figure of means
  df_bio %>% ggplot(aes(f_age, mean_length, colour = perc)) +
    stat_summary(fun.data = "mean_cl_boot", position = position_dodge(width = 0.5)) +
    theme_minimal() +
    scale_colour_viridis_d(end = 0.8, option = "B") 
  df_bio$perc <- factor(df_bio$perc, levels = c("10","20","30","40","50","75","100"))
  
  df_bio %>% ggplot(aes(mean_length, perc,  colour = f_age)) +
    stat_summary(fun.data = "mean_cl_boot") +
    theme_minimal() +
    scale_colour_viridis_d(end = 0.8, option = "B") 
  
  # figure of sd
  df_bio %>% ggplot(aes(year, sd_length, colour = perc)) +
    stat_summary(fun.data = "mean_cl_boot", position = position_dodge(width = 2)) +
    theme_minimal() +
    scale_colour_viridis_d(end = 0.8, option = "B")
  
  # figure of cv
  df_bio %>% ggplot(aes(f_age, cv_length, colour = perc)) +
    stat_summary(fun.data = "mean_cl_boot", position = position_dodge(width = 0.5)) +
    theme_minimal() +
    scale_colour_viridis_d(end = 0.8, option = "B")
  summary(df_bio)

  ### pseudo accuracy as per Stamatopoulos 1999. Observations on the geometrical properties of accuracy growth in sampling with finite populations. FAO tech paper 338 vii 39p
  # A = 1-((sample mean - population mean(assumed to be 100))/(ymax-ymin assumed to be pop mean))      
  
  df_bio %<>% 
    dplyr::select(-sd_length, -cv_length, -n) %>% 
    pivot_wider(values_from = mean_length, names_from = perc) %>% 
    mutate(acc_10 = (1 - abs(`10`-`100`)/`100`),
           acc_20 = (1 - abs(`20`-`100`)/`100`),
           acc_30 = (1 - abs(`30`-`100`)/`100`),
           acc_40 = (1 - abs(`40`-`100`)/`100`),
           acc_50 = (1 - abs(`50`-`100`)/`100`),
           acc_75 = (1 - abs(`75`-`100`)/`100`))  
  
  df_bio %<>% select(year, quarter, bioregion, gear_group, f_age, acc_10, acc_20, acc_30,acc_40, acc_50, acc_75) %>% 
    pivot_longer(cols = 6:11, names_to = 'perc', values_to = 'accuracy')
  
  df_bio %>% ggplot(aes(year, accuracy, colour = perc)) + stat_summary(fun.data = "mean_cl_boot", position = position_dodge(width = 1))  + facet_wrap(vars(f_age)) + scale_colour_viridis_d()
  df_bio %>% ggplot(aes(quarter, accuracy, colour = perc)) + stat_summary(fun.data = "mean_cl_boot", position = position_dodge(width = 0.5))  + facet_wrap(vars(f_age)) + scale_colour_viridis_d()
  df_bio %>% ggplot(aes(bioregion, accuracy, colour = perc)) + stat_summary(fun.data = "mean_cl_boot", position = position_dodge(width = 0.5))  + facet_wrap(vars(f_age)) + scale_colour_viridis_d()
  df_bio %>% ggplot(aes(gear_group, accuracy, colour = perc)) + stat_summary(fun.data = "mean_cl_boot", position = position_dodge(width = 0.5))  + facet_wrap(vars(f_age)) + scale_colour_viridis_d()
  df_bio %>% ggplot(aes(perc, accuracy, colour = f_age)) + stat_summary(fun.data = "mean_cl_boot", position = position_dodge(width = 0.5))  + scale_colour_viridis_d()
  df_bio %>% ggplot(aes(perc, accuracy, colour = f_age)) + stat_summary(fun.data = "mean_cl_boot", position = position_dodge(width = 0.5))  + scale_colour_viridis_d()+ facet_wrap(vars(f_age))+
  ylim(0,1)
  
  ### catch at age
  
  # Baird, J.W. 1983. A Method to Select Optimum Numbers for Aging in a Stratified Approach. 
  # In W. G. Doubleday and/et D. Rivard [ed./ed.] Sampling commercial catches
  # of marine fish and invertebrates/L'léchantillonnage des prises commerciales de poissons
  # et d'invertebres marins. Can. Spec. Pub!. Fish. Aquat. Sci./Publ. spec. can. sci. halieut. aquat 66.
  #
  # 
  # Goal: The objective of sampling for ages is to obtatin similar
  # coefficients of variation (CV) for all ages of fish that contribute significantly to the population
  # To this en, here we evaluate number of fish and number of samples needed to have satisfactory precision in the Northwest Atlantic's northern contingent of Atlantic mackerel
  
  # Variables:
  # Ni - catch number at length i
  # ni - number aged at length i
  # pi - proportion at length for a given age
  # nipi - number aged at length that are a given age
  # Nipi - catch at length for a given age
  # k - stratification (year/quarter/division/gear)
  # An appropriate age-length-key 
  
  # Methodology: Follow eqns in Baird 1983 for each age.
  # N = sum(Nipi) catch at age is the sum catch at each length for that age
  # Var(Nipi) = Ni^2*Var(pi) + pi^2*Var(Ni) - Gulland 1955 variance of the catch at each length group
  # Second term can be dropped
  # Var(Nipi) = Ni^2*Var(pi) where Var(pi) = pi*(1-pi)/ni
  # Var(N) = sum(Ni^2)*Var(pi)
  # CV = Var(N)/N
  
  # keep it simple here.
  
  load("~/Data/Landings/Ziff/Rdata/mackerel_ziff.Rdata")
  
  
  
  # new variable : quarter
  mackerel_ziff %<>% mutate(quarter = fct_collapse(
    factor(mackerel_ziff$month),
    Q1 = c("1", "2", "3"),
    Q2 = c("4", "5", "6"),
    Q3 = c("7", "8", "9"),
    Q4 = c("10", "11", "12")
  ))
  
  # Explicit NA for grouping
  unique(mackerel_ziff$gear_name)
  mackerel_ziff$gear_name <- forcats::fct_explicit_na(mackerel_ziff$engin_en, na_level = "Unknown")
  mackerel_ziff$division <- as.factor(str_sub(mackerel_ziff$opano, 1, 2))
  mackerel_ziff$division <- forcats::fct_explicit_na(mackerel_ziff$division, na_level = "4S")
  mackerel_ziff$subdivision <- forcats::fct_explicit_na(mackerel_ziff$div)
  
  mackerel_ziff  %<>%  mutate(
    gear_caa = fct_collapse(
      factor(mackerel_ziff$gear_name),
      "nets_traps_weirs" = c(
        "Trap net",
        "Weir",
        "Box net",
        "Bag net",
        "Fyke net",
        "Beach and bar seine"
      ),
      "seiners" = c("Purse seine",
                    "Tuck seine",
                    "Danish seine",
                    "Scottish seine"),
      "gillnets" = c(
        "Gillnet (set or fixed)",
        "Gillnet (drift)",
        "Gillnet (unspecified)"
      ),
      "lines" = c(
        "Hand line (baited)",
        "Jigger",
        "Longline",
        "Rod and reel (chumming)",
        "Mechanized squid jigger",
        "Hand and hand held tools",
        "Automated jigger  (hand line)"
      ),
      "trawlers" = c("Bottom otter trawl (stern)",
                     "Shrimp trawl",
                     "Bottom otter trawl (side)",
                     "Midwater pair trawl",
                     "Midwater trawl (stern)"),
      "misc" = c(
        "Pot",
        "Unknown",
        "Miscellaneous",
        "Setheared hooks",
        "Null data"
      )
    )
  )
  
  ziff <- mackerel_ziff %>% 
    group_by(year, quarter, division, gear_caa,nbpc, month, doy) %>% 
    dplyr::summarise(catch = sum(pds_vif,na.rm = T)/1000)
  ziff %<>% dplyr::filter(catch>0)
  unique(ziff$gear_caa)
  rm(mackerel_ziff)
  
  ## bio
  
  glimpse(bio)
  unique(bio$gear_name)
  bio$gear_name <- forcats::fct_explicit_na(bio$gear_name, na_level = "Unknown")
  bio  %<>%  mutate(
    gear_caa = fct_collapse(
      factor(bio$gear_name),
      "nets_traps_weirs" = c(
        "Trap_net",
        "Weir",
        "Lampara_net",
        "Trap_na",
        "Beach_seine"
      ),
      "seiners" = c("Purse_seine",
                    "Tuck_seine",
                    "Danish_seine",
                    "Pair_seine"),
      "gillnets" = c(
        "Gillnet_set",
        "Gillnet_drift",
        "Gillnet_na"
      ),
      "lines" = c(
        "Line_na",
        "Handline",
        "Mechanised_jigger"
      ),
      "trawlers" = c("Bottom_trawl",
                     "Shrimp_trawl_no_grid",
                     "Shrimp_trawl_na",
                     "Pelagic_trawl"),
      "misc" = c(
        "Misc.",
        "Unknown",
        "NA",
        "Longline"
      )
    )
  )
  unique(bio$gear_caa)
  
  unique(lf$gear_name)
  lf$gear_name <- forcats::fct_explicit_na(bio$gear_name, na_level = "Unknown")
  lf  %<>%  mutate(
    gear_caa = fct_collapse(
      factor(lf$gear_name),
      "nets_traps_weirs" = c(
        "Trap_net",
        "Weir",
        "Lampara_net",
        "Trap_na",
        "Beach_seine"
      ),
      "seiners" = c("Purse_seine",
                    "Tuck_seine",
                    "Danish_seine",
                    "Pair_seine"),
      "gillnets" = c(
        "Gillnet_set",
        "Gillnet_drift",
        "Gillnet_na"
      ),
      "lines" = c(
        "Line_na",
        "Handline",
        "Mechanised_jigger"
      ),
      "trawlers" = c("Bottom_trawl",
                     "Shrimp_trawl_no_grid",
                     "Shrimp_trawl_na",
                     "Pelagic_trawl"),
      "misc" = c(
        "Misc.",
        "Unknown",
        "NA",
        "Longline"
      )
    )
  )
  unique(lf$gear_caa)
  
  # subset
  
  bio_caa <- bio %>% dplyr::select(source, sample_type, year, month, doy, quarter, division, gear_caa, length, mass, f_age, sample_id)
  bio_caa %<>% dplyr::filter(!is.na(sample_id)) 
  lf_caa <- lf %>% dplyr::select(source, sample_type, year, month, doy, quarter, division, gear_caa, nbpc, mass_sample, mass_landed,  sample_id, length)
  rm(bio);rm(lf)
  
  # Add length category variable (5mm as per current two stage sampling protocol (i.e. ~ 100-200 fish sampled for length frequencies (sortie.dat) in proportion to landings and 2 fish per 5mm length class sent to IML to be processed for biological trait data (Maquereau_carbio.dat))) 
  bio_caa %<>% 
    mutate(lcat5 = FSA::lencat(length, w = 5))
  
  lf_caa %<>% 
    mutate(lcat5 = FSA::lencat(length, w = 5))
  
  # group some divisions
  unique(lf_caa$division)
  lf_caa  %<>%   mutate(
    div_caa = fct_collapse(
      factor(lf_caa$division),
      "2J3KLO" = c(
        "2J",
        "2G",
        "3K",
        "3L",
        "3O"),
      "3P4V" = c("4VN",
                 "4VS",
                 "4V",
                 "3P"),
      "4XW5YZ" = c(
        "4X",
        "4W",
        "5Y",
        "5Z",
        "5ZE"),
      "4RS" = c(
        "4R",
        "4S"),
      "4T" = c("4T")))
  
  bio_caa  %<>%   mutate(
    div_caa = fct_collapse(
      factor(bio_caa$division),
      "2J3KLO" = c(
        "2J",
        "2G",
        "3K",
        "3L",
        "3O"),
      "3P4V" = c("4VN",
                 "4VS",
                 "4V",
                 "3P"),
      "4XW5YZ" = c(
        "4X",
        "4W",
        "5Y",
        "5Z",
        "5ZE"),
      "4RS" = c(
        "4R",
        "4S"),
      "4T" = c("4T")))
  
  ziff %<>% ungroup()
  ziff  %<>%    mutate(
    div_caa = fct_collapse(
      factor(ziff$division),
      "2J3KLO" = c(
        "2J",
        "2G",
        "3K",
        "3L",
        "3O"),
      "3P4V" = c("4VN",
                 "4VS",
                 "4V",
                 "3P"),
      "4XW5YZ" = c(
        "4X",
        "4W",
        "5Y",
        "5Z",
        "5ZE"),
      "4RS" = c(
        "4R",
        "4S"),
      "4T" = c("4T")))
  
  ziff %<>% mutate(q_caa = fct_collapse(
    factor(ziff$quarter),
    Q12 = c("Q1", "Q2")))
  
  bio_caa %<>% ungroup() %>% 
    mutate(q_caa = fct_collapse(
      factor(bio_caa$quarter),
      Q12 = c("Q1", "Q2")))
  
  lf_caa %<>% ungroup() %>% 
    mutate(q_caa = fct_collapse(
      factor(lf_caa$quarter),
      Q12 = c("Q1", "Q2")))
  
  glimpse(ziff)
  
  
  # Ok let's just do this by year, quarter, and div_caa
  
  ziff %<>% group_by(year, q_caa, div_caa) %>% 
    dplyr::summarise(catch = sum(catch, na.rm = T))
  
  # LW relationship
  # Fit model to each year/quarter stratification
  ml_fits <- bio_caa %>% 
    dplyr::filter(!is.na(mass), !is.na(length), mass > 0) %>% # a few weird values
    group_by(year, q_caa) %>% 
    nest() %>% 
    mutate(fit = map(data, ~ lm(log(mass) ~ log(length), data = .x))) 
  
  # Nest the length frequency data by the same stratification  
  lf_list <- lf_caa %>% 
    group_by(year, q_caa) %>% 
    nest() %>% 
    transmute(year = year, q_caa = q_caa, lf_data = data) 
  
  # Merge the datasets by the common stratification
  ml_fits <- right_join(ml_fits, lf_list, by = c("year", "q_caa")) %>%
    arrange(year, q_caa)
  
  # Get predictions from original data ie adjusted/fitted values
  ml_fits <- ml_fits %>% 
    mutate(predictions = map(fit, augment, type.predict = "response"))
  
  # Predict masses from corresponding stratified length frequency data
  ml_fits <- ml_fits %>%
    mutate(predictions_new = map2(fit, lf_data, ~ broom::augment(.x, newdata = .y, type.predict = "response")))
  
  # Get model coefficients and back-transform log transformed response variable (mass) and then multiply by 1/2 sigma^2 to correct for known bias
  # first filter to >1984
  ml_fits %<>% dplyr::filter(year>1984)
  ml_fits <- ml_fits %>% 
    mutate(syx = map_dbl(fit, sigma),
           cf = exp((syx^2)/2))
  
  # Pluck the nested length frequency data with their new predicted masses 
  ml_fits <- ml_fits %>% 
    unnest(predictions_new) %>% 
    mutate(mass = exp(.fitted)*cf)
  # - can extract other lists too like model coeffs etc mackerel_ml_fits %>% unnest(glanced)
  
  # Subset
  ml_fits %<>% dplyr::select(1,2,19,20,24) %>% as_tibble()

  # subset bio to same
  bio_caa %<>% dplyr::select(year, q_caa, div_caa, f_age, lcat5, mass) %>% as_tibble() %>% 
    dplyr::filter(year>1984, f_age != "(Missing)") %>% droplevels()
  unique(bio_caa$f_age)
  
  # merge bio and lf, split by has age vs none (get to fill in gaps bio might  have)
  df<-bind_rows(ml_fits, bio_caa)
  df$f_age<-as.numeric(as.character(df$f_age))
  
  lf_caa <- df %>% dplyr::filter(is.na(f_age), lcat5>140)
  bio_caa <-df %>% dplyr::filter(!is.na(f_age), lcat5>140)
  
  all(is.na(lf_caa$f_age))
  any(is.na(bio_caa$f_age))
  
  # Apply ALK to unaged data
  
  # classical way
  aa <- bind_rows(lf_caa, bio_caa) 
  len.n <- xtabs(~lcat5,data=aa)
  
  # Nest data
  df <- bio_caa %>% 
    group_by(year) %>% 
    nest() %>%    # grouping variables for weighting
    transmute(bio_data = data)
  
  # Calculate the frequencies and proportions
  df %<>% mutate(alk.freq = map(bio_data, ~ xtabs(~lcat5 + f_age, data = .x)))
  df %<>% mutate(alk = map(alk.freq, ~ prop.table(.x, margin = 1)))
  
  # Nest the length frequency data by the same stratification (same procedure as in 2) above)
  bb <- aa %>% dplyr::filter(!is.na(lcat5)) %>% 
    group_by(year) %>% 
    nest() %>% 
    transmute(lf_data = data)
  
  # Merge data frames
  df2 <- full_join(bb, df)
  df2 %<>% dplyr::filter(year != 2021)
  df2 %<>% mutate(len.n = map(lf_data, ~ xtabs(~lcat5, data = .x)))
  df2 %<>% mutate(tmp = map2(alk, len.n, ~ sweep(.x, MARGIN=1, FUN="*",STATS=.y)))
  df2 %<>% mutate(tmp2 = map(tmp, ~as_tibble(.x)))
  df3 <- df2%>% unnest(tmp2) %>% dplyr::select(year, lcat5, f_age, n)
  df4 <- df2 %>% unnest(lf_data) %>% dplyr::select(year, f_age, mass) %>% dplyr::filter(!is.na(mass))
  # now I have numbers for year and age and can thin
  #
  df3 %<>% group_by(year, f_age) %>% 
    dplyr::summarise(n = sum(n)) %>% 
    arrange(year, f_age)
  df3$f_age <- as.numeric(df3$f_age)
  
  df4 %<>% group_by(year, f_age) %>% 
    dplyr::summarise(mass = mean(mass,na.rm=T))

  df5 <- left_join(df3, df4)
  
  df5 %<>% mutate(raised_mass = (n * mass)/1000)
  
  ziff %<>% group_by(year) %>% dplyr::summarise(catch = sum(catch))
  
  df6<-left_join(df5,ziff)
  
  df6 %<>% mutate(caa = n*(catch/raised_mass)) %>% 
    group_by(year, f_age) %>% 
    mutate(mean_mass1 = mean(mass/1000), mean_mass2 = mean((mass/1000))*(catch/raised_mass)) %>% 
  ungroup() %>% 
  mutate(caw1 = caa*mean_mass1,caw2 = caa*mean_mass2)
      df6 %>% group_by(year, f_age) %>% dplyr::summarise(cw = caa*mean_mass1) %>% tail()
  
  bio_20 <- bio %>% 
    dplyr::filter(quarter != "Q1", sample_type == "commercial_port", gear_group %in% c("Nets_traps_weirs", "Gillnets", "Purse_seine", "Tuck_seine", "Lines"), f_age != "(Missing)") %>% 
    group_by(year, quarter, bioregion, gear_group) %>% 
    slice_sample(prop = .2, replace = F) %>% 
    ungroup() %>% 
    group_by(year, quarter, bioregion, gear_group, f_age) %>% 
    dplyr::summarise(n = n(), mean_length = mean(length, na.rm = T), sd_length = sd(length,na.rm=T), cv_length = (sd_length/mean_length)*100) %>% 
    ungroup() %>%  
    mutate(perc = "20")
  