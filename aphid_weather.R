#Aphid weather data
#3/23/17

#load required packages
library(data.table)
library(weatherData)
library(ggmap)
library(rnoaa)

#set working directory
setwd("C:/Users/Chad/Documents/GitHub/sapsuckRz/sapsuckRz/")

##load in aphid data
#Using fread in the data.table package because this file is quite large
aphid<-fread("https://dl.dropboxusercontent.com/u/98197254/Mastersuction_may29_2014.csv")

##Pull the 10 most common species
#find the average number of captures for each species
aphid.ag<-aggregate(Captures~Aphid.Species,data=aphid,sum)
#order data from most number of captures to least
aphid.ord<-aphid.ag[order(-aphid.ag$Captures),] 
#make a list of the 10 most common species
aphid.names.cut<-aphid.ord$Aphid.Species[1:10]
#cut down the data to only include those 10 species
aphid.cut<-aphid[aphid$Aphid.Species%in%aphid.names.cut]

#load location data
location<-read.csv("counties+coordinates2.csv")

#find cities based on lat long
cities<-vector(mode="character",length=46)
for(i in 1:46){
  cities[i]<- as.character(revgeocode(c(location$LONG[i], location$LAT[i]), output = "more")$locality)
}
location<-location[1:3]
colnames(location)<-c("id","LAT","LONG")
station_id<-vector(mode="character",length=46)
for(i in 1:46){
  station_id[i]<-as.data.frame(meteo_nearby_stations(location[i,], lat_colname = "LAT",
                      lon_colname = "LONG",limit=1))[1]
}

for(i in 1:46){
weather[i]<-meteo_pull_monitors(station_id[i], date_min = "2003-01-01",
                    date_max = "2016-01-01", var = c("PRCP","TMAX","TAVG"))
}

station_id<-vector(mode="character",length=46)
for(i in 1:46){
station_id[i]<-getStationCode(cities[i])$airportCode
}

for(i in x){
weather<- getWeatherForDate(i, "2011-08-26","2011-08-26")
}

weather_ppt<-get_prism_dailys(type="ppt", minDate = "2003-01-01", maxDate = "2016-01-01", keepZip=FALSE)

##
for(i in 1:46){
station_id[i]<-getMeta(lat = location$LAT[i], lon = location$LONG[i],n=1)[12]
}

weather<-vector(mode="character",length=46)
for(i in 1:46){
weather[i]<-importNOAA(code = as.character(station_id[i]), year = 2003:2015, hourly = TRUE,
           precip = TRUE)
}

##
type=c("tmax","tmin","ppt")

for(i in 1:3){
  var <- paste("weather", i, sep = "")
  assign(var, get_prism_dailys(type=type[i], minDate = "2003-01-01", maxDate = "2016-01-01", keepZip=T))
}

##
library(knitr)
opts_knit$set(upload.fun = imgur_upload)

library("RFc")
library("raster")
library("rworldmap")
data("coastsCoarse")

##
#install packages and load libraries
#install.packages("RSocrata")
#install.packages("rgeos")

library(RSocrata)
library(rgeos)
library(sp)
library(weatherData)

#load data
location<-read.csv("counties+coordinates2.csv")

#pull airport codes
airports<- read.socrata("https://opendata.socrata.com/dataset/Airport-Codes-mapped-to-Latitude-Longitude-in-the-/rxrh-4cxm")

#Make Western hemisphere points negative to match other data
airports$Longitude<-airports$Longitude*-1

#set lat long to coordinate data type
location_sp <- SpatialPoints(location[,2:3])
airports_sp <- SpatialPoints(airports[,2:3])

#find the closes airport to each site
location$nearest_airport <- apply(gDistance(location_sp, airports_sp, byid=TRUE), 2, which.min)
airports$nearest_airport<-rownames(airports)

location<-merge(location,airports,by="nearest_airport",all.x=T)
location<-location[,2:7]

weather<-NULL
for(i in 46){
  weather<- rbind(weather,data.frame(getWeatherForDate(location$locationID[i], "2005-04-01","2013-11-30",
                                 opt_custom_columns = T,custom_columns = c(2,4,20)),location$SITEID[i]))
}

##

library(weatheR)

#find cities based on lat long
cities<-vector(mode="character",length=46)
for(i in 1:46){
  cities[i]<- as.character(revgeocode(c(location$LONG[i], location$LAT[i]), output = "more")$locality)
}

location$cities<-cities
location$states<-state.name[match(location$STATE,state.abb)]
location$loc<-paste(location$cities,location$state, sep=", ")

station.list <- allStations()
getFilteredStationsByCity(location$loc, station.list, begin = 2005, end = 2013)
weather<-getInterpolatedDataByCity(location$loc, station.list, begin = 2005, end = 2013,tolerance = 0.25)

library(reshape2)

weather.cast<-dcast(weather_data,param~data)

weather.cast$Date <- as.Date(x=weather.cast$Date, format='%d/%m/%Y')

# Extract day of the week (Saturday = 6)
weather.cast$Week_Day <- as.numeric(format(weather.cast$Date, format='%w'))

# Adjust end-of-week date (first saturday from the original Date)
weather.cast$End_of_Week <- weather.cast$Date + (6 - weather.cast$Week_Day)

# Aggregate over week and climate division
aggregate(cbind(precipitation,tmax,tmin)~End_of_Week+location, FUN=sum, data=weather.data, na.rm=TRUE)

#Add line to show how github works
