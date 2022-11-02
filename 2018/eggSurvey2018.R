######################################################################################################
# Calculate SSB index from mackerel egg survey
#
# Structure and analyses based on Research Document 2013/035 DFO (Francois Gregoire)
# spatial code based on script Hugo Bourdage (krigeage_capelan_BL.r)
######################################################################################################

### Some initial preparations (packages, functions, wd, ...) --------
dir="C:/Users/VANBE/Desktop/post-doc/DATA/Mackerel Bio and Egg data/"

packages=c('plyr','ggplot2','gridExtra','sp','gstat','raster','rgeos','rgdal','PBSmapping','geoR','XLConnect','xlsx')
invisible(lapply(packages, function(x) {if (!require(x, character.only=T)) {install.packages(x);require(x)}}))

readexcel <- function(file=file,sheet=sheet){
  #require(XLConnect, pos=4)  #from excel file
  .Workbook <- XLConnect::loadWorkbook(file)
  if(length(sheet)>1){
    l=list()
    for(k in 1:length(sheet)){l[[k]]=XLConnect::readWorksheet(.Workbook, sheet[k])}
  }else{l <- XLConnect::readWorksheet(.Workbook, sheet)}
  remove(.Workbook)
  return(l)
}
maxN <- function(x, N=2){
  x=x[!is.na(x)]
  len <- length(x)
  if(N>len){
    warning('N greater than length(x).  Setting N=length(x)')
    N <- length(x)
  }
  sort(x,partial=len-N+1)[as.numeric(len-N+1):len]
}
source(paste0(dir,"moy.var.Krigeage.R"))   # function from Hugo

theme_new <- theme_set(theme_classic())
theme_new <- theme_update(axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black"))

files <- dir(paste0(dir,"metadonnees"), pattern = "M?ta")

### load and merge data --------
surface.strate=data.frame(STRATE=1:3, surface.strate = c(29.61e+9,21.91e+9,17.93e+9))
volu=readexcel(file=paste0(dir,'Volume_N-P.xls'),sheet=c('SOGolfe 1982-1998'))
volu=rename(volu,c('ANN?E'='year','STATION'='station')) # in correspondance with the next file
stations=volu[,c('station','STRATE','LONG','LAT')]
stations=stations[!duplicated(stations),]
stations=stations[complete.cases(stations),]
stations=merge(stations,surface.strate,all.x=T)

eggs=data.frame()
for(i in files){
sheets <- getSheets(loadWorkbook(paste0(dir,"metadonnees/",i))) #load sheet names
nrs=grep('lanc',names(sheets))[1] # check the number of the one with plancton in
nrs=c((nrs-1),nrs) #sheets of interest: the planton one and the one before (wich changes names completely)
metado=readexcel(file=paste0(dir,"metadonnees/",i),sheet=nrs)
data=metado[[1]]
planc=metado[[2]]

##standardise names
names(planc)=tolower(names(planc))
names(planc)[grep('sous',names(planc))][1]='subsample'  #sous echantillon (the first one is generally for mackerel)
names(planc)[grep('total',names(planc))][1]='egg.total'  #"Total.oeufs.Scomber.s.", "total.oeufs....maquereau"
names(planc)[grep('5',names(planc))]='state5'
names(planc)[grep('1',names(planc))]='state1'
names(planc)[grep('ann',names(planc),ignore.case=T)][1]='year'
names(data)=tolower(names(data))
names(data)[grep('t[.]',substr(names(data), 1, 2))][1]='temp'  # not all files have same column names T (see 2011): this assumes the first one is the right one!!
names(data)[grep('net',names(data))]='volume' # violume bionet
names(data)[grep('ongo',names(data))][1]='depth' # profondeur bongo
names(data)[grep('ann',names(data),ignore.case=T)][1]='year'
names(data)[grep('tri',names(data),ignore.case=T)][1]='side' #babord/tribord

##cleanup
if(data$year[1]==2007){planc[planc$cl?.unique==7 &!is.na(planc$cl?.unique),'station']=2.2} #2.2a to 2.2 (there is also a 2.2b, but plancton sample from a)
# na's and 0s where appropriate
planc$state1=ifelse((is.na(planc$subsample) & planc$egg.total==0),0,planc$state1)  
planc$state5=ifelse((is.na(planc$subsample) & planc$egg.total==0),0,planc$state5)
planc$state1=ifelse((is.na(planc$state1) & !is.na(planc$subsample)),0,planc$state1)
planc$state5=ifelse((is.na(planc$state5) & !is.na(planc$subsample)),0,planc$state5)

#format
data$volume=as.numeric(data$volume) #creates true NAs and not characters or so
check<- try(as.Date(data$date), silent=TRUE)
if ('try-error' %in% class(check)){data$date=as.Date(paste(data$year,data$date,sep="-"),format="%Y-%d-%m")} else{data$date=as.Date(data$date)}
planc$subsample=as.numeric(planc$subsample)

# remove (mostly) empty rows
data=data[rowSums(is.na(data)) != ncol(data),]#remove exmpty rows
planc=planc[rowSums(is.na(planc)) != ncol(planc),]#remove exmpty rows
data=data[!is.na(data$station),] #remove rows with no station numbers
too.many.nas=which(apply(planc, MARGIN = 1, FUN = function(x) length(x[is.na(x)]) )>(ncol(planc)-8))
if(length(too.many.nas)>0){planc=planc[-too.many.nas,]}#remove rows with mostly NA
 
#bionet values replaced by flowmeter values if the former are missing
missing=length(data[is.na(data$volume) & data$side=='B',"volume"]) 
if(missing>0){
  data$volume <- ifelse(is.na(data$volume), data$volume..m3., data$volume)
  print(paste0('In file ', i, ', ', missing, 'x the bionet volume was missing (and replaced by the flowmeter volume when available)' ))
}

#paste it all together
eggs1=merge(data[data$side=='B',c('station','year','date','volume','depth',
                                        'temp')],planc[,c('station','state1','state5','subsample')],all=T)
eggs1=merge(eggs1,stations,all=T)
eggs=rbind(eggs,eggs1)
}

print(paste("Years read in: ",paste(sort(as.numeric(unique(eggs$year))),collapse=" ")))
eggs=eggs[!is.na(eggs$year),]

#chekc duplicates
duplis=eggs[,c('year','station')][duplicated(eggs[,c('year','station')]),]
print(paste(nrow(duplis),"duplicate(s)"))
eggs=eggs[!(eggs$year=='2012' & eggs$station==4.2 & eggs$volume>290),] # B station 4.2 entered twice (but with small diff in V)
eggs=eggs[!(eggs$year=='2007' & eggs$station==2.2 & eggs$volume!=57.87) ,] # keep the one that has the egg data

eggs.save=eggs


#################### CALCS ###############################################
##### 2.3 Calculation of the egg abundance (N/M2) by station ###############################################
eggs$state15=eggs$state1+eggs$state5 # number of eggs stage 1 and 5 counted
eggs$N=ifelse((is.na(eggs$subsample)|eggs$subsample==0),0,eggs$state15/eggs$subsample )    # number of eggs (1&5) in sample
eggs$Nv=eggs$N/as.numeric(eggs$volume)   # number (1&5) per cubic meter
eggs$Nm2=eggs$Nv*as.numeric(eggs$depth)  # number (1&5) per square meter
eggs$Nm2=ifelse(is.na(eggs$Nm2) & eggs$N==0,0,eggs$Nm2) # if no eggs found but volume is missing, set Nm2 to zero
eggs[is.na(eggs$Nm2),] ### THESE LINES SHOULD BE CHECKED IN EXCEL FILE !!!  removed in next line
eggs=eggs[!is.na(eggs$Nm2),]

eggs$temp=as.numeric(eggs$temp)
eggs[is.na(eggs$temp),] ### THESE LINES DON'T HAVE T DATA (TAKE FROM ROSETTE?)

# analysed according to water T: Fig 4/5 
# all years combined
eggs.sub=eggs[,c('temp','Nm2','year')]
eggs.sub=eggs.sub[order(eggs.sub$temp),]
eggs.sub$cumsum=cumsum(eggs.sub$Nm2)
eggs.sub$prop=eggs.sub$cumsum/max(eggs.sub$cumsum)

  T.cdf.plot=ggplot (eggs) + stat_ecdf(aes(x=temp),col='grey') +xlab('') +ylab('Proportion')+
    geom_point(data=eggs.sub,aes(x=temp,y=prop),col='black')+
    scale_x_continuous(limits=c(0,18)) 
  T.Nm2.plot=ggplot(eggs,aes(x=temp,y=Nm2))+geom_point() +xlab('') +ylab('Abundance (N/m2)')+ 
    geom_text(data=eggs[eggs$Nm2 %in% maxN(eggs$Nm2,10),],aes(label=station),hjust=-0.1,vjust=-0.1)+
    scale_y_continuous(limits=c(0,1000))
grid.arrange(T.Nm2.plot,T.cdf.plot,ncol=2,bottom='Temperature')

#per year
eggs.sub=eggs[,c('temp','Nm2','year')]
eggs.sub=eggs.sub[order(eggs.sub$year,eggs.sub$temp),]
eggs.sub=ddply(eggs.sub,c('year'),transform,cumsum=cumsum(Nm2))
eggs.sub=ddply(eggs.sub,c('year'),transform,prop=cumsum/max(cumsum))

T.cdf.plot=ggplot (eggs) + stat_ecdf(aes(x=temp,col=year,group=year)) +xlab('') +ylab('Proportion')+
  geom_point(data=eggs.sub,aes(x=temp,y=prop,col=year))+
  scale_x_continuous(limits=c(0,18)) +
  labs(col='Year')+
  scale_colour_continuous(low=c('blue','yellow'),high='red')
T.Nm2.plot=ggplot(eggs,aes(x=temp,y=Nm2,col=year))+geom_point() +xlab('') +ylab('Abundance (N/m2)')+ 
  geom_text(data=eggs[eggs$Nm2 %in% maxN(eggs$Nm2,10),],aes(label=station),hjust=-0.1,vjust=-0.1)+
  scale_y_continuous(limits=c(0,1000))+
  scale_colour_continuous(low=c('blue','yellow'),high='red')+
  theme(legend.position='none')
png(filename=paste0(dir,"/IMG/Temp_all.png"),units="cm",res=200,width=18,height=10)
grid.arrange(T.Nm2.plot,T.cdf.plot,ncol=2,bottom='Temperature')
dev.off()


for(i in unique(eggs$year)){
  T.cdf.plot=ggplot (eggs[eggs$year==i,]) + stat_ecdf(aes(x=temp)) +xlab('') +ylab('Proportion')+
    geom_point(data=eggs.sub[eggs.sub$year==i,],aes(x=temp,y=prop))+
    scale_x_continuous(limits=c(0,18)) +
    labs(col='Year')+
    scale_colour_continuous(low=c('red','yellow'),high='blue')
  T.Nm2.plot=ggplot(eggs[eggs$year==i,],aes(x=temp,y=Nm2))+geom_point() +xlab('') +ylab('Abundance (N/m2)')+ 
    geom_text(data=eggs[eggs$Nm2 %in% maxN(eggs$Nm2,10) &eggs$year==i,],aes(label=station),hjust=-0.1,vjust=-0.1)+
    scale_y_continuous(limits=c(0,1000))+
    scale_colour_continuous(low=c('red','yellow'),high='blue')+
    theme(legend.position='none')
  png(filename=paste0(dir,"/IMG/Temp_",i,".png"),units="cm",res=200,width=18,height=10)
  grid.arrange(T.Nm2.plot,T.cdf.plot,ncol=2,bottom='Temperature')
  dev.off()  
}

##### 2.4 Calculation of the incubation times (HR) ########################################################
eggs$I=exp((-1.61*log(eggs$temp))+7.76) #incubation time in hours

##### 2.5 calculation of the daily egg production (N/m2) by station ######################################
eggs$DEP= eggs$Nm2/eggs$I*24


      #compare to tables Francois for 2011
      # table1
      eggs$Nm2round=round(eggs$Nm2,digits=1)
      check=eggs[eggs$year==2011,c('station','Nm2round','temp','I','STRATE')]
      mean(check$Nm2round) # should be 30.0
      sd(check$Nm2round) # should be 102.4
      length(check$Nm2round) # should be 65 (ok)
      # conclusion: some differences (0.1) and the first station (1.1, value 0) was removed by Francois?
      # table 2
      mean(check$temp) # should be 9.5. ????
      mean(check$temp[-1])
      check$Tempround=round(check$temp,digits=1)
      mean(check$tempround[-1])
      # Temperatures don't match at all!!!!
      #table 3
      head(check)
      mean(check$I) # should be 65.8
      # logically, this is wrong as well (and so will be table 4). The first one not taking into account, and T data was wrong as well
      # table 5, number of stations per stratum
      table(check$STRATE) # should be 15, 16,15. ???? thus total would only be 46??

##### 2.6 calculation of the daily egg production for the entire sampled area #############################
### without krigging: Ouellet method------------
eggs$DEP.weigthed=eggs$DEP/eggs$surface.strate
surf=sum(surface.strate$surface.strate)
print(paste(nrow(eggs[is.na(eggs$DEP.weigthed),]),'station(s) deleted because not enough info'))
eggs=eggs[!is.na(eggs$DEP.weigthed),]
  
DEP1=ddply(eggs,c('year'),summarise,DEP1=mean(DEP.weigthed)*surf)

### with krigging---------------

## get all data
# egg survey data points
dat=eggs[,c('year','station','LONG','LAT','DEP')]

t=dat
dat=t

# test: use Francois' data for 2011
# dat=dat[dat$year==2011,]
# dat$DEP=c(NA,0,0,0,0,0.2,0,0.3,0,0,0,1,0.3,0.8,0.6,0.1,0,0,0.3,0,1.1,0.9,28.1,0.4,0.5,0.7,0,0,0,1.2,2.9,1.2,0.6,
#           0.4,0.2,0,27.7,6.1,5.8,0.3,0.3,1.7,0.1,18.8,2.1,0.6,43,9.3,11.1,4.9,0.8,69.4,102.8,21,2.3,15.9,10.5,
#           0.8,102.4,392.2,0.8,12.4,9.1,14.6,1.1)
# dat=dat[!is.na(dat$DEP),]


# grid
grille = read.table(paste0(dir,"grille_pred.csv"), header = TRUE, sep=";")
grille = grille[complete.cases(grille),] 
grille = grille[-4:-5]   # On enl?ve le x et y existant afin de faire la m?me transformation g?od?sique pour toutes les donn?es (Hugo)

# shape file (coast line)
      # world <- shapefile("C:/Users/vanbe/Desktop/PhD/part2/DATA/shapes/best/GSHHS_f_L1.shp")
      # sGSL <- crop(world, extent(-68, -59.7, 45, 50))
      # plot(sGSL)
      # writeOGR(sGSL, dsn = '.', layer = 'sGSL', driver = "ESRI Shapefile")
sGSL <- shapefile(paste0(dir,"sGSL.shp"))

## Initialisation du syst?me de coordonn?es pour la conversion des donn?es de latlong ? lambert et vice-versa
lat_1=48
lat_2=50
lat_0=46.5
lon_0=-70
x_0=0
y_0=0
ellps="clrk66"
units="km"

# Pr?paration de la conversion vers Lambert  (Hugo Bourdage)
LCC <- paste("+proj=lcc", 
             " +lat_1=",lat_1, " +lat_2=",lat_2, 
             " +lat_0=",lat_0," +lon_0=", lon_0, 
             " +x_0=", x_0, " +y_0=",y_0, 
             " +ellps=", ellps, " +units=", units, 
             " +no_defs +nadgrids=@conus,@alaska,@ntv2_0.gsb,@ntv1_can.dat", sep="")
# Pr?paration de la conversion vers LatLong    (Hugo Bourdage)
LL <- "+proj=longlat +datum=NAD27 +ellps=clrk66 +nadgrids=@conus,@alaska,@ntv2_0.gsb,@ntv1_can.dat" 

# Conversion vers lambert
mackerel=dat
coordinates(dat)=~LONG+LAT  #sp
proj4string(dat) <- CRS("+proj=longlat +grilleum=NAD27 +ellps=clrk66 +nadgrilles=@conus,@alaska,@ntv2_0.gsb,@ntv1_can.grille") 
dat = spTransform(dat, CRS(LCC))
dat@coords = round(dat@coords, digits=5)
colnames(dat@coords)<-c("x","y")
mackerel_sp <- dat  # capelan_sp est un SpatialPointsDataFrame

coordinates(grille) <- ~lon+lat
proj4string(grille) <- CRS("+proj=longlat +grilleum=NAD27 +ellps=clrk66 +nadgrilles=@conus,@alaska,@ntv2_0.gsb,@ntv1_can.grille") 
grille = spTransform(grille, CRS(LCC))
grille@coords = round(grille@coords, digits=3)
colnames(grille@coords)<-c("x","y")
grille_sp <- grille  # grille_sp est un SpatialPointsgrilleaFrame

sGSL = spTransform(sGSL, CRS(LCC))
class(sGSL)
slotNames(sGSL)

sGSL_map_readjust = gSimplify(sGSL, tol = 0.00001) # I recropped the shape for plotting because with the new projection it is skewed
sGSL_map_readjust = crop(sGSL_map_readjust, extent(250, 800, -100, 400))
sGSL_map_readjust = fortify(sGSL_map_readjust) # dataframe for ggplot

## Bubble plot
mapdata = data.frame(mackerel_sp)
mapdata$eggs= ifelse(mapdata$DEP==0,"absent","present")
sGSL_map = fortify(sGSL)

plot.bubble = ggplot() +geom_polygon(data=sGSL_map_readjust, aes(x=long, y=lat, group=group))+  
  geom_point(data=mapdata, aes(x=x, y=y,size=DEP,color=eggs))+
  facet_wrap(~year)+
  scale_color_manual(values=c('grey','red'))+
  scale_x_continuous(limits=c(250,800),expand=c(0,0))+
  scale_y_continuous(expand=c(0,0)) # change color/shape for zero

png(filename=paste0(dir,"/IMG/Bubble_all.png"),units="cm",res=200,width=35,height=20)
plot.bubble
dev.off()

# # ## Variogramme
# # #-- if no variogram yet: determine one and add to excel file varios_mackerel
# # # individual variogram
# my.vgm = variogram(DEP~1, data=mackerel_sp[mackerel_sp@data$year=='2006',]) # calculates sample variogram values
# plot(my.vgm,pch=16, plot.numbers=TRUE)
# # Dans un premier temps, j'ajuste un vario ? l'oeil!!!
# model.i = vgm(psill=367, model="Gau", range=20, nugget=1) #sill: y value plateau, range: x value plateau, nugget: intersect y-axis
# plot(my.vgm, model=model.i, plot.numbers=T,pch=16)
# # Dans un deuxi?me temps, on ajuste le vario automatiquement ? partir des valeurs initiales trouv?es ? l'oeil
# model.f = fit.variogram(my.vgm, model=model.i, fit.sills=T, fit.ranges=T, fit.method=6)
# plot(my.vgm, model=model.f, plot.numbers=T, pch=16)
# model.f

#-- if variograms have already been made
varios=readexcel(file=paste0(dir,'Varios_mackerel.xlsx'),sheet=c(1))

### Krigeage and mean/variance of total area
pred=data.frame()  # dataframe with predictions krigging
moy_var=data.frame() # dataframe wih mean and variance of krigging predictions of whole area
mackerel$pred=NA  

for(i in unique(eggs$year)){
v=vgm(psill=varios[varios$YEAR==i,'Sill'], 
      model=varios[varios$YEAR==i,'MODEL'], 
      range=varios[varios$YEAR==i,'Range'], 
      nugget=varios[varios$YEAR==i,'Nugget'])

my.vgm = variogram(DEP~1, data=mackerel_sp[mackerel_sp@data$year==i,]) # calculates sample variogram values
png(filename=paste0(dir,'/IMG/Variogram_',i,".png"),units="cm",res=200,width=20,height=12)
print(plot(my.vgm, model=v, plot.numbers=T, pch=16))
dev.off()

krige = krige(DEP~1, mackerel_sp[mackerel_sp@data$year==i,], grille_sp, v, nmin=16)
pred1 = data.frame(year=rep(i,length(grille_sp)),depth=grille_sp@data$prof,
                  pred=krige@data$var1.pred,
                  var=krige@data$var1.var,
                  LONG=grille_sp@coords[,1],LAT=grille_sp@coords[,2])
pred=rbind(pred,pred1)

mackerel[mackerel$year==i,'pred']=krige@data$var1.pred[unname(apply(gDistance(SpatialPoints(krige), SpatialPoints(mackerel_sp[mackerel_sp@data$year==i,]), byid=TRUE), 1, which.min))] #for predicted vs observed

moy_var1=moy.var.Krigeage(
  xCoords = coordinates(mackerel_sp[mackerel_sp@data$year==i,]), 
  xValues = mackerel[mackerel$year==i,]$DEP, 
  kCoords = coordinates(grille_sp), 
  vario = v)
moy_var1=as.data.frame(moy_var1)
moy_var1$year=i
moy_var1$z.noneg=mean(pred1[!pred1$pred<0,'pred'])
moy_var=rbind(moy_var,moy_var1)
}

pred$neg=ifelse(pred$pred<0,'-','') # Negative predictions indicated on map

# cor and RSS (as in table 6)
   # sampling stations are not located on grid coordinates:  find nearest predicted value
mackerel$res=mackerel$DEP-mackerel$pred
FIT=ddply(mackerel,c('year'),summarise,COR=cor(DEP,pred),RSS=sum(res^2))
FIT

#plot 
for(i in unique(eggs$year)){
plot.krig.pred=ggplot() + stat_summary_2d(data=pred[pred$year==i,],aes(x = LONG, y = LAT,z=pred),fun = median, binwidth = c(6, 6))+
  geom_text(data=pred[pred$year==i,],aes(x=LONG,y=LAT,label=neg),col='grey')+
  geom_polygon(data=sGSL_map_readjust, aes(x=long, y=lat, group=group))+
  geom_point(data=mapdata[mapdata$year==i,], aes(x=x, y=y,size=DEP), alpha=0.3)+
  scale_fill_gradient(name = "Median",
                       low="yellow",high="red",
                       space = "Lab")+ 
  labs(x = "Longitude",y = "Latitude") +
  scale_x_continuous(limits=c(250,800),expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))

plot.krig.var=ggplot() + stat_summary_2d(data=pred[pred$year==i,],aes(x = LONG, y = LAT,z=var),fun = median, binwidth = c(6, 6))+
  geom_polygon(data=sGSL_map_readjust, aes(x=long, y=lat, group=group))+
  geom_point(data=mapdata[mapdata$year==i,], aes(x=x, y=y,size=DEP), alpha=0.3)+
  scale_fill_gradient(name = "Median",
                      low="yellow",high="red",
                      space = "Lab")+ 
  labs(x = "Longitude",y = "Latitude") +
  scale_x_continuous(limits=c(250,800),expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))

png(filename=paste0(dir,'/IMG/Krigging_',i,".png"),units="cm",res=200,width=28,height=12)
grid.arrange(plot.krig.pred+ggtitle('Predicted'),plot.krig.var+ggtitle('Variance'),ncol=2,top=i)
dev.off() 
}

# summary values for whole area
DEP2=moy_var   # This is the same as mean(pred$pred). 

DEP=data.frame(year=unique(eggs$year),DEP1=DEP1$DEP1,DEP2=DEP2$Z) #z.noneg if no negatifs

##### 2.7 calculation of the proportion of the eggs spawned daily #############################
# median date of the survey
date.med=ddply(eggs,c('year'),summarise,date=median(date,na.rm=T)) ## 2008 has a date missing!
date.med$doy=as.numeric(strftime(date.med$date, format = "%j"))

# GSI logistic curve
load(file=paste0(dir,'bio2.RData'))
GSI=bio2[bio2$PDIND!=0 &!is.na(bio2$PDIND)&!is.na(bio2$PDGON),c('date','PDIND','PDGON','ccyear')]
GSI$GSI=GSI$PDGON/GSI$PDIND*100
GSI$doy=as.numeric(strftime(GSI$date, format = "%j"))+0.5 #Day Of Year

co=data.frame() #coefficients of fit
spawning.prop=data.frame(year=sort(unique(eggs$year)),prop=NA,peak.day=NA) # spawning proportion at data, for SSB calculation

for(i in sort(unique(eggs$year))){
GSI1=GSI[GSI$ccyear==i,]
        # if removal outliers
        # logistic.out=nls(GSI ~ y0 + (a/(1+(doy/x0)^b)),
        #                    data = GSI1,
        #                    start = list(y0 = 0.4451, x0 = 172.4085, a = 11.3252,b = 32.0552))
        # 
        # GSI1$pred=predict(logistic.out)
        # GSI1$res=residuals(logistic.out)
GSI1$outlier='no'    #ifelse(abs(GSI1$res)>5,'yes','no') # just some subjective simple rule
if(i==2011){GSI1[GSI1$doy==min(GSI1$doy),]$outlier='yes'} #outliers
GSI1[GSI1$GSI>100,'outlier']='yes'

removed=table(GSI1$outlier)[2]/nrow(GSI1)*100
  
GSI2=GSI1[GSI1$outlier=='no',]
N=nrow(GSI2)
logistic.out=nls(GSI ~ y0 + (a/(1+(doy/x0)^b)),
                 data = GSI2,
                 start = list(y0 = 0.4451, x0 = 172.4085, a = 11.3252,b = 32.0552))

co1=coef(logistic.out)
co2=data.frame(year=i,yo=co1[1],x0=co1[2],a=co1[3],b=co1[4])
co=rbind(co,co2)

logi.plot=ggplot(GSI1,aes(x=doy,y=GSI))+geom_point(aes(col=outlier))+
  stat_function(fun=function(x) co1[1] + (co1[3]/(1+(x/co1[2])^co1[4])),col='red',size=1.5)+
  geom_text(label=paste('N =',N),x=120,y=20)+
  geom_text(label=paste('removed =',round(removed,digits=1),"%"),x=120,y=18)+
  xlab('')+
  scale_color_manual(values=c('black','grey'))+
  theme(legend.position=c(0.8,0.8))+
  scale_x_continuous(limits=c(100,305))

d1=data.frame(day=seq(100.5,350.5,1))
d1$pred=co1[1] + (co1[3]/(1+(d1$day/co1[2])^co1[4]))
d2=data.frame(day=seq(101,350,1))
for(j in 1:(nrow(d1)-1)){d2[j,'slope']=d1[j+1,'pred']-d1[j,'pred']}
d2$prob=d2$slope/sum(d2$slope,na.rm=T)

prob.plot=ggplot(d2,aes(x=day,y=prob*100))+geom_line()+ylab('Probability (%)')+xlab('')+
  geom_vline(aes(xintercept=date.med[date.med$year==i,'doy']),col='red')+
  scale_x_continuous(limits=c(100,305))

png(filename=paste0(dir,'/IMG/LogisticFit_',i,".png"),units="cm",res=200,width=20,height=18)
grid.arrange(logi.plot,prob.plot,ncol=1,bottom='Day of the year')
dev.off()
spawning.prop[spawning.prop$year==i,'peak.day']=d2[d2$prob==max(d2$prob),"day"]
spawning.prop[spawning.prop$year==i,'prop']=d2[d2$day==date.med[date.med$year==i,'doy'],"prob"]
}

co
spawning.prop

png(filename=paste0(dir,"/IMG/DOY_spawning_vs_survey.png"),units="cm",res=200,width=20,height=12)
plot(spawning.prop$peak.day~spawning.prop$year,type='l',main='Day of peak spawning (black) vs day survey (red)',ylab='DOY',xlab='Year',ylim=c(155,190)) 
lines(date.med$year,date.med$doy,col='red')
dev.off()

##### 2.8 Calculation of the total annual egg production #############################
DEP$TAEP1= DEP$DEP1*sum(surface.strate$surface.strate)/spawning.prop$prop
DEP$TAEP2= DEP$DEP2*sum(surface.strate$surface.strate)/spawning.prop$prop

##### 2.9 calculation of the spawning biomass index #############################
SSB=data.frame(year=unique(DEP$year),P=NA,A=NA,W=NA,S=NA,Fec=NA,R=NA,convers=NA,SSB=NA)

for(i in unique(DEP$year)){
P=DEP[DEP$year==i,'DEP2']                                                # N/m2 from krigging
A=sum(surface.strate$surface.strate)                  # m2
W=mean(bio2[bio2$ccyear==i & bio2$PDIND!=0,'PDIND'],na.rm=T)       # mean weight fish (g)
S=spawning.prop[spawning.prop==i,'prop']                                    # proportion of eggs spawned at median date survey
Fec=mean(10^(4.32+0.75*log10(bio2[bio2$ccyear==i & bio2$STAD_MAT==5 & !is.na(bio2$STAD_MAT) & bio2$SEXE=='F','PDGON']))) # fecundity of females (pelletier uses log10!!!)
R=length(bio2[bio2$ccyear==i & bio2$SEXE == 'F','SEXE'])/
  length(bio2[bio2$ccyear==i & bio2$SEXE %in% c('F','M'),'SEXE'])  # proportion of females
convers=10^6                                          # grams to tonnes

SSB1=round((P*A*W)/(S*Fec*R*convers),digits=0)  # at least 100x too large now

SSB[SSB$year==i,]=c(i,P,A,W,S,Fec,R,convers,SSB1)
}

# compare with what francois had
result.francois=structure(list(year = c(1996, 1997, 1998, 1999, 2000, 2001, 2002, 
                        2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013
), SSB_Francois = c(123464, NA, 105801, NA, 161573, NA, 389007, 307091, 
                162802, 87959, NA, 76532, 99631, 73743, 25960, 35714, 14568, 
                68547)), .Names = c("year", "SSB_Francois"), row.names = c(NA, -18L
                ), class = "data.frame")

SSB=merge(SSB,result.francois,all.x=T)
SSB$perc=round(SSB$SSB/SSB$SSB_Francois,digits=2)
SSB

#plot
ggplot(SSB,aes(x=year,y=SSB))+geom_line()+
  geom_line(aes(y=SSB_Francois),col='grey')+
  scale_y_continuous(limits=c(0,max(c(SSB$SSB,SSB$SSB_Francois),na.rm=T)*1.05),expand=c(0,0))

