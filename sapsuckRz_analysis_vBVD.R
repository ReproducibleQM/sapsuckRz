# sapsuckRz_anylsis vBVD

# load required packages
library(data.table)
library(ggplot2)
library(vegan)
library(dplyr)
library(reshape2)
library(lubridate)
library(lme4)
library(r2glmm)
library(car)
library(mgcv)
source("custom_functions.r")

# Data pull and prep ####
# End result is the dataframes "data", "data.aglyc", "data.rpadi", and "data.rmaidis"

#########/////// Aphids

##load in aphid data
#Using fread in the data.table package because this file is quite large
aphid<-fread("https://ndownloader.figshare.com/files/9465745")
aphid<-aphid[,-1]
colnames(aphid)<-c("Year","State","Site","Date","Sex","Captures","Species")
aphid.all<-aggregate(Captures~Year+Date+Site,data=aphid,sum)
aphid.aglyc <- aggregate(Captures ~ Site + Date + Year, data = aphid[aphid$Species == "Aphis.glycines"], sum)
aphid.rpadi <- aggregate(Captures ~ Site + Date + Year, data = aphid[aphid$Species == "Rhopalosiphum.padi"], sum)
aphid.rmaidis <- aggregate(Captures ~ Site + Date + Year, data = aphid[aphid$Species == "Rhopalosiphum.maidis"], sum)

#########/////// Weather
weather<-read.csv("https://ndownloader.figshare.com/files/8448380")
#make weather data wide
weather.cast<-dcast(weather, date+siteid~param,value.var="data",sum)
names(weather.cast)<-c("date","site.id","max.temp","min.temp","precip")

#calculate degree days for each site in each year

#first replace missing values with interpolated data
weather.cast$max.temp<-replace.missing(weather.cast$max.temp)
weather.cast$min.temp<-replace.missing(weather.cast$min.temp)
weather.cast$precip<-replace.missing(weather.cast$precip)

#create day of year variable
weather.cast$date<-ymd(weather.cast$date)
weather.cast$doy<-yday(weather.cast$date)

#reorder so that we calculate degree days for each site in each year
weather.cast<-weather.cast[with(weather.cast, order(site.id,date)), ]

#flip some temp data because min is higher than max
weather.cast[weather.cast$min.temp > weather.cast$max.temp, c("max.temp", "min.temp")] <- weather.cast[weather.cast$min.temp > weather.cast$max.temp, c("min.temp", "max.temp")] 

#degree days
weather.cast$dd<-allen(weather.cast$max.temp, weather.cast$min.temp, 10.000001)

#degree day accumulation
#starting March 1
weather.cast$dd.acum<-accum.allen(weather.cast$max.temp,weather.cast$min.temp,10.000001,weather.cast$doy,60)

##precipitation data
#precipitation accumulation over the growing season
#starting March 1
weather.cast$precip.accum<-accum.precip.time(weather.cast$precip,weather.cast$doy,60)

##Combine the weather and aphid data##
names(aphid.all)<-c("year","doy","site.id","captures")#so they match
names(aphid.aglyc)<-c("year","doy","site.id","captures")
names(aphid.rpadi)<-c("year","doy","site.id","captures")
names(aphid.rmaidis)<-c("year","doy","site.id","captures")


weather.cast$year<-isoyear(weather.cast$date)

data<-merge(aphid.all,weather.cast,by=c("site.id","year","doy"),all.x=T)
data.aglyc<-merge(aphid.aglyc,weather.cast,by=c("site.id","year","doy"),all.x=T)
data.rpadi<-merge(aphid.rpadi,weather.cast,by=c("site.id","year","doy"),all.x=T)
data.rmaidis<-merge(aphid.rmaidis,weather.cast,by=c("site.id","year","doy"),all.x=T)


#########/////// CDL and Geography
cdl<-read.csv("Data/cdl.csv")
cdl.cut<-cdl[,c(2:3,11:50)]
colnames(cdl.cut)[1]<-"year"
colnames(cdl.cut)[2]<-"site.id"

#add landcover data to dataframe
cdl.cut$year<-as.factor(cdl.cut$year)
data$site.id<-as.factor(data$site.id)
cdl.cut<-aggregate(.~site.id+year,data=cdl.cut,mean)

data<-merge(data,cdl.cut,by=c("site.id","year"),all.x=T)
data.agylc<-merge(data.agylc,cdl.cut,by=c("site.id","year"),all.x=T)
data.rpadi<-merge(data.rpadi,cdl.cut,by=c("site.id","year"),all.x=T)
data.rmaidis<-merge(data.rmaidis,cdl.cut,by=c("site.id","year"),all.x=T)

#adding lat long
coords<-read.csv("Data/counties+coordinates2.csv")
coords<-coords[,1:3]
colnames(coords)<-c("site.id","lat","long")
data<-merge(data,coords,by="site.id",all.x=T)
data.agylc<-merge(data.agylc,coords,by="site.id",all.x=T)
data.rpadi<-merge(data.rpadi,coords,by="site.id",all.x=T)
data.rmaidis<-merge(data.rmaidis,coords,by="site.id",all.x=T)


# Dependent Variable Construction #####
# This produces a set of dataframes, "max.[species]", which provides max captures over season, "captures", and date of max captures "doy"

max<-as.data.frame(data %>% group_by(year,site.id) %>% slice(which.max(captures)))
max.aglyc<-as.data.frame(data.aglyc %>% group_by(year,site.id) %>% slice(which.max(captures)))
max.rpadi<-as.data.frame(data.rpadi %>% group_by(year,site.id) %>% slice(which.max(captures)))
max.rmaidis<-as.data.frame(data.rmaidis %>% group_by(year,site.id) %>% slice(which.max(captures)))