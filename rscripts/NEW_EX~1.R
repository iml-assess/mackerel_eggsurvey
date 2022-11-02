# script to look at trade off of capelin relative condition vs gsi

# load and look at data
glimpse(Capelin_LF)
bio <- Capelin_LF %>% dplyr::filter(!is.na(mass_gonad),mass>mass_gonad)
bio %<>% mutate(mass_somatic = mass-mass_gonad, gsi=(mass_gonad/mass_somatic)*100) %>% dplyr::filter(mass_somatic>5) # filter extremes and abberations

# standard power law association between length and mass (both somatic and total)
fit1f<-lm(log(mass_somatic)~log(length),data=bio %>% dplyr::filter(sex =="F"))
fit1m<-lm(log(mass_somatic)~log(length),data=bio %>% dplyr::filter(sex =="M"))

fit2f<-lm(log(mass)~log(length),data=bio %>% dplyr::filter(sex =="F"))
fit2m<-lm(log(mass)~log(length),data=bio %>% dplyr::filter(sex =="M"))

library(performance)
compare_performance(fit1f,fit2f,fit1m,fit2m)
library(broom)

# append fitted values to original dataframe * remember they are exponentiated
augf1<-augment(fit1f,newdata = bio %>% dplyr::filter(sex =="F"))
augf2<-augment(fit2f,newdata = bio %>% dplyr::filter(sex =="F"))
augm1<-augment(fit1m,newdata = bio %>% dplyr::filter(sex =="M"))
augm2<-augment(fit2m,newdata = bio %>% dplyr::filter(sex =="M"))

# get model error to correct for back calculation bias 
syx1f = sigma(fit1f)
cf1f = exp((syx1f^2)/2)

syx2f = sigma(fit2f)
cf2f = exp((syx2f^2)/2)

syx1m = sigma(fit1m)
cf1m = exp((syx1m^2)/2)

syx2m = sigma(fit2m)
cf2m = exp((syx2m^2)/2)

# calculate relative condition
augf1 %<>% mutate(std.mass = exp(.fitted)*cf1f,Kn = (mass_somatic/std.mass))
augf2 %<>% mutate(std.mass = exp(.fitted)*cf1f,Kn = (mass_somatic/std.mass))
augm1 %<>% mutate(std.mass = exp(.fitted)*cf1f,Kn = (mass/std.mass))
augm2 %<>% mutate(std.mass = exp(.fitted)*cf1f,Kn = (mass/std.mass))

# merge dataframes
df_1 <- bind_rows(augf1,augm1) %>% arrange(date)
df_2 <- bind_rows(augf2,augm2) %>% arrange(date)

# figures to look at data
df_1 %>% 
  ggplot(aes(year,Kn,colour = region)) + 
  stat_summary(fun.data = "mean_cl_boot")+ facet_wrap(vars(sex),scales="free_y") +
  geom_hline(yintercept=1)

df_2 %>% 
  ggplot(aes(year,Kn,colour = region)) + 
  stat_summary(fun.data = "mean_cl_boot")+ facet_wrap(vars(sex),scales="free_y") +
  geom_hline(yintercept=1)

df_1 %>% 
  ggplot(aes(Kn,gsi,colour=jday))+geom_point()+facet_wrap(vars(sex),scales="free_y")+scale_colour_viridis_c()

# condition decreases exponentially with increasing gsi .. ie trade off between somatic growth and gonad developpement
df_1 %>% dplyr::filter(sex == "F", month<10,month>4) %>% 
  ggplot(aes(gsi,Kn,colour = factor(region)))+geom_smooth()+scale_colour_viridis_d()

# standardize the variables (z score aka anomalies) to compare kn and gsi on same scale with mean = 0 and sd = 1
df_1 %<>% group_by(region,sex) %>% dplyr::mutate(std.kn = scale(Kn), std.gsi = scale(gsi))

# look at data

Capelin_LF %>% 
  dplyr::filter(year>1985) %>% 
  group_by(year, Region, sex) %>% 
  dplyr::summarise(mean_l = mean(length)) %>% 
  ggplot(aes(year,mean_l, colour = Region))+
  stat_summary(fun.data = "mean_cl_boot")+
  geom_line()+
  facet_wrap(vars(sex),ncol = 1)+
  theme_minimal()

df_1 %>% dplyr::filter(sex == "F") %>% 
  ggplot(aes(std.gsi,std.kn,colour = factor(region)))+geom_point()+scale_colour_viridis_d()  + geom_hline(yintercept = 0) + geom_vline(xintercept=0)

df_1 %>% 
  dplyr::filter(sex == "F") %>% 
  ggplot(aes(std.gsi,std.kn)) +
  geom_point()+scale_colour_viridis_d() +
  geom_smooth()+
  geom_hline(yintercept = 0) +
  geom_vline(xintercept=0)+
  facet_wrap(vars(maturity))

df_1 %>% dplyr::filter(sex == "F", month %in% c(5,6,7)) %>% 
  ggplot(aes(std.gsi,std.kn,colour = factor(region)))+geom_point()+scale_colour_viridis_d()  + geom_hline(yintercept = 0) + geom_vline(xintercept=0)

df_1 %>% dplyr::filter(sex == "F", month %in% c(5,6,7)) %>% 
  ggplot(aes(std.gsi,std.kn,colour = factor(region)))+geom_point()+scale_colour_viridis_d()  + geom_hline(yintercept = 0) + geom_vline(xintercept=0)

df_1 %>% dplyr::filter(sex == "F", month %in% c(5,6,7)) %>% 
  ggplot(aes(std.gsi,std.kn))+geom_smooth()  + geom_hline(yintercept = 0) + geom_vline(xintercept=0)

df_1 %>% dplyr::filter(sex == "F", month %in% c(4,5,6,7,8)) %>% 
  ggplot(aes(std.gsi,std.kn,colour = sex))+geom_smooth()  + geom_hline(yintercept = 0) + geom_vline(xintercept=0)

df_1 %>% dplyr::filter(sex == "F", 
                       month %in% c(6,7,8), 
                       region == "NE_GSL",
                       gear_name %in% c("Seiners","Nets_Traps_Weirs"),
                       maturity %in% c(4,5), 
                       mass_gonad > 0) %>% 
  ggplot(aes(jday,std.kn)) + 
  stat_summary(fun.data = "mean_cl_boot") +
    facet_wrap(vars(gear_name),scales = "free_y")

# compare gsi and Kn
df<- df_1 %>% dplyr::filter(sex == "F", maturity == "4")
cor.test(df$std.gsi,df$std.kn) # cor = -0.7218573
plot(df$std.gsi,df$std.kn)

df_1 %>% 
  dplyr::filter(sex == "F", 
                length > 100,
                mass > 5,
                month %in% c(5,6,7,8), 
                subdivision %in% c("4RA","4RB","4RC","4RD","4SW"),
                gear_name %in% c("Seiners"),
                # maturity %in% c(4,5), 
                mass_gonad > 0) %>% 
  ggplot(aes(year, std.kn)) + stat_summary(fun.data = "mean_cl_boot") +
  stat_summary(aes(year,std.gsi),colour="orange",fun.data = "mean_cl_boot") +
  geom_hline(yintercept = 0) +
  facet_wrap(vars(maturity),scales="free_y")

df_1 %>% 
  dplyr::filter(sex == "F", 
                length > 100,
                mass > 5,
                month %in% c(5,6,7,8), 
                subdivision %in% c("4RA","4RB","4RC","4RD","4SW"),
                gear_name %in% c("Seiners"),
                # maturity %in% c(4,5), 
                mass_gonad > 0) %>% 
  ggplot(aes(year, std.kn)) + stat_summary(fun.data = "mean_cl_boot") +
  stat_summary(aes(year,std.gsi),colour="orange",fun.data = "mean_cl_boot") +
  geom_hline(yintercept = 0) +
  ylab("std.value") +
  annotate("text",x = 2000, y = 1.5, label = "orange = std.gsi / black = std.kn")

df_1 %>% 
  dplyr::filter(sex == "F", 
                length > 100,
                mass > 5,
                month %in% c(5,6,7,8), 
                jday>150,
                subdivision %in% c("4RA","4RB","4RC","4RD","4SW"),
                gear_name %in% c("Seiners"),
                # maturity %in% c(4,5), 
                mass_gonad > 0) %>% 
  ggplot(aes(jday, std.kn)) + stat_summary(fun.data = "mean_cl_boot") +
  stat_summary(aes(jday,std.gsi),colour="orange",fun.data = "mean_cl_boot") +
  geom_hline(yintercept = 0) +
  ylab("std.value") +
  annotate("text",x = 200, y = 3.5, label = "orange = std.gsi / black = std.kn") 
  # facet_wrap(vars(subdivision))



df_1 %>% dplyr::mutate(lcat = ifelse(length<140,"< 14cm","> 14 cm")) %>% 
  dplyr::filter(sex == "F",region=="NE_GSL", gear_name =="Seiners") %>% 
  ggplot(aes(year, std.kn)) + stat_summary(fun.data = "mean_cl_boot") +
  stat_summary(aes(year,std.gsi),colour="orange",fun.data = "mean_cl_boot") +
  geom_hline(yintercept = 0) +
  facet_wrap(vars(lcat),scales="free_y")

df_1 %>% dplyr::filter(sex == "F", month %in% c(5,6,7)) %>% 
  ggplot(aes(year, std.kn)) + stat_summary(fun.data = "mean_cl_boot") +
  stat_summary(aes(year,std.gsi),colour="orange",fun.data = "mean_cl_boot") +
  geom_hline(yintercept = 0) +
  facet_wrap(vars(region),ncol = 1)

