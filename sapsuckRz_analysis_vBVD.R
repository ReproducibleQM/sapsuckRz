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
library(zoo)
source("custom_functions.r")

# Data pull and prep ####
# End result is the dataframes "data", "data.aglyc", "data.rpadi", and "data.rmaidis"

#########/////// Aphids

##load in aphid data
#Using fread in the data.table package because this file is quite large
aphid<-fread("https://ndownloader.figshare.com/files/9465745")
aphid<-aphid[,-1]
colnames(aphid)<-c("Year","State","Site","Date","Sex","Species","Captures")
aphid.all<-aggregate(Captures~Year+Date+Site,data=aphid,sum)
aphid.aglyc <- aggregate(Captures ~ Site + Date + Year, data = aphid[aphid$Species == "Aphis.glycines"], sum)
aphid.rpadi <- aggregate(Captures ~ Site + Date + Year, data = aphid[aphid$Species == "Rhopalosiphum.padi"], sum)
aphid.rmaidis <- aggregate(Captures ~ Site + Date + Year, data = aphid[aphid$Species == "Rhopalosiphum.maidis"], sum)

##Pull the 25 most common species
#find the average number of captures for each species
aphid.ag<-aggregate(Captures~Species,data=aphid,sum)
#order data from most number of captures to least
aphid.ord<-aphid.ag[order(-aphid.ag$Captures),] 
#make a list of the 10 most common species
aphid.names.cut<-aphid.ord$Species[1:25]
#and top 10
aphid.names.cut.10<-aphid.ord$Species[1:10]
#cut down the data to only include those 10 species
aphid.cut<-aphid[aphid$Species%in%aphid.names.cut]

#Calculate maximum for each species by site by year
max.top25<-as.data.frame(aphid.cut %>% group_by(Year,Site,Species) %>% slice(which.max(Captures)))
#remove any zero maxes
max.top25<-max.top25[max.top25$Captures>0,]
names(max.top25)<-c("year","state","site.id","doy","sex","species","captures")

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

####Calculate weather variables for analysis####
##All aphids
#Calculate rolling averages degree day accumulation for all aphids
data$dd.week<-rollmean(data$dd, 7, align = "right", fill = NA)
data$dd.month<-rollmean(data$dd, 30, align = "right", fill = NA)
data$dd.year<-data$dd.acum/data$doy

#Calculate rolling average precipitation for all aphids
data$precip.week<-rollmean(data$precip, 7, align = "right", fill = NA)
data$precip.month<-rollmean(data$precip, 30, align = "right", fill = NA)
data$precip.year<-data$precip.accum/data$doy

#Calculate total precipitation before peak for all aphids
data$precip.week.tot<-rollsum(data$precip, 7, align = "right", fill = NA)
data$precip.month.tot<-rollsum(data$precip, 30, align = "right", fill = NA)

##A. glyc
#Calculate rolling averages degree day accumulation for A. glyc
data.aglyc$dd.week<-rollmean(data.aglyc$dd, 7, align = "right", fill = NA)
data.aglyc$dd.month<-rollmean(data.aglyc$dd, 30, align = "right", fill = NA)
data.aglyc$dd.year<-data.aglyc$dd.acum/data.aglyc$doy

#Calculate rolling average precipitation for A. glyc
data.aglyc$precip.week<-rollmean(data.aglyc$precip, 7, align = "right", fill = NA)
data.aglyc$precip.month<-rollmean(data.aglyc$precip, 30, align = "right", fill = NA)
data.aglyc$precip.year<-data.aglyc$precip.accum/data.aglyc$doy

#Calculate total precipitation before peak for A. glyc
data.aglyc$precip.week.tot<-rollsum(data.aglyc$precip, 7, align = "right", fill = NA)
data.aglyc$precip.month.tot<-rollsum(data.aglyc$precip, 30, align = "right", fill = NA)

##R. padi
#Calculate rolling averages degree day accumulation for R. padi
data.rpadi$dd.week<-rollmean(data.rpadi$dd, 7, align = "right", fill = NA)
data.rpadi$dd.month<-rollmean(data.rpadi$dd, 30, align = "right", fill = NA)
data.rpadi$dd.year<-data.rpadi$dd.acum/data.rpadi$doy

#Calculate rolling average precipitation for R. padi
data.rpadi$precip.week<-rollmean(data.rpadi$precip, 7, align = "right", fill = NA)
data.rpadi$precip.month<-rollmean(data.rpadi$precip, 30, align = "right", fill = NA)
data.rpadi$precip.year<-data.rpadi$precip.accum/data.rpadi$doy

#Calculate total precipitation before peak for R. padi
data.rpadi$precip.week.tot<-rollsum(data.rpadi$precip, 7, align = "right", fill = NA)
data.rpadi$precip.month.tot<-rollsum(data.rpadi$precip, 30, align = "right", fill = NA)

##Rmadis
#Calculate rolling averages degree day accumulation for R. maidis
data.rmaidis$dd.week<-rollmean(data.rmaidis$dd, 7, align = "right", fill = NA)
data.rmaidis$dd.month<-rollmean(data.rmaidis$dd, 30, align = "right", fill = NA)
data.rmaidis$dd.year<-data.rmaidis$dd.acum/data.rmaidis$doy

#Calculate rolling average precipitation for R. maidis
data.rmaidis$precip.week<-rollmean(data.rmaidis$precip, 7, align = "right", fill = NA)
data.rmaidis$precip.month<-rollmean(data.rmaidis$precip, 30, align = "right", fill = NA)
data.rmaidis$precip.year<-data.rmaidis$precip.accum/data.rmaidis$doy

#Calculate total precipitation before peak for R. maidis
data.rmaidis$precip.week.tot<-rollsum(data.rmaidis$precip, 7, align = "right", fill = NA)
data.rmaidis$precip.month.tot<-rollsum(data.rmaidis$precip, 30, align = "right", fill = NA)

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
data.aglyc<-merge(data.aglyc,cdl.cut,by=c("site.id","year"),all.x=T)
data.rpadi<-merge(data.rpadi,cdl.cut,by=c("site.id","year"),all.x=T)
data.rmaidis<-merge(data.rmaidis,cdl.cut,by=c("site.id","year"),all.x=T)

#adding lat long
coords<-read.csv("Data/counties+coordinates2.csv")
coords<-coords[,1:3]
colnames(coords)<-c("site.id","lat","long")
data<-merge(data,coords,by="site.id",all.x=T)
data.aglyc<-merge(data.aglyc,coords,by="site.id",all.x=T)
data.rpadi<-merge(data.rpadi,coords,by="site.id",all.x=T)
data.rmaidis<-merge(data.rmaidis,coords,by="site.id",all.x=T)

# Dependent Variable Construction #####
# This produces a set of dataframes, "max.[species]", which provides max captures over season, "captures", and date of max captures "doy"

max<-as.data.frame(data %>% group_by(year,site.id) %>% slice(which.max(captures)))
max.aglyc<-as.data.frame(data.aglyc %>% group_by(year,site.id) %>% slice(which.max(captures)))
max.rpadi<-as.data.frame(data.rpadi %>% group_by(year,site.id) %>% slice(which.max(captures)))
max.rmaidis<-as.data.frame(data.rmaidis %>% group_by(year,site.id) %>% slice(which.max(captures)))

##add predictor variables to top 25 data.frame
data.predictor<-data[,-4]
max.top25<-merge(max.top25,data.predictor,by=c("site.id","year","doy"),all.x=T)
max.top10<-max.top25[max.top25$species%in%aphid.names.cut.10,]

##captures=============
#model
#gam
#mod.gam.top25<-bam(captures ~ precip.week+precip.month+dd.week+dd.month+dd.year+ag_corn10+ag_beans10+ag_smgrains10+forest10+
#          s(lat,long)+s(year,site.id,species,bs="re"), data = max.top10,family="nb")
#summary(mod.gam.top25)

#gamm
mod.gam.top10<-gamm(captures ~ precip.week+precip.month+dd.year+I(dd.year)^2+ag_corn75+ag_beans75+ag_smgrains75+forest75+
                     s(lat,long),random=list(year=~1,site.id=~1,species=~1),control = lmeControl(opt='optim'), data =max.top10,family=negative.binomial(1))
summary(mod.gam.top10$gam)

#ignoring species
mod.gam.max<-gamm(captures ~ precip.week+precip.month+dd.year+ag_corn75+ag_beans75+ag_smgrains75+forest75+
                      s(lat,long),random=list(year=~1,site.id=~1), control = lmeControl(opt='optim'), data =max,family=negative.binomial(1))
summary(mod.gam.max$gam) ##we get the same answer

#glmer
#mod.max.25<-glmer.nb(captures ~ precip.week+precip.month+dd.year+ag_corn10+ag_beans10+ag_smgrains10+forest10+
#                       (1|year)+(1|site.id)+(1|species),data=max.top10)
#summary(mod.max.25)
#gamm4
#mod.gam.top25.2<-gamm4(captures ~ precip.week+precip.month+dd.year+ag_corn10+ag_beans10+ag_smgrains10+forest10+
                     # s(lat,long),random= ~ (1|year)+(1|site.id)+(1|species), data = max.top10,family=poisson)
#summary(mod.gam.top25.2$gam)

###captures by species=========
max.aglyc<-max.top10[max.top10$species%in%"Aphis.glycines",]
max.rpadi<-max.top10[max.top10$species%in%"Rhopalosiphum.padi",]
max.rmaidis<-max.top10[max.top10$species%in%"Rhopalosiphum.maidis",]

#aglyc max
mod.gam.aglyc<-gamm(captures ~ precip.week+precip.month+dd.year+ag_corn10+ag_beans10+ag_smgrains10+forest10+
                      s(lat,long),random=list(year=~1,site.id=~1), data =max.aglyc,family=negative.binomial(1))
summary(mod.gam.aglyc$gam)

#rpadi max
mod.gam.rpadi<-gamm(captures ~ precip.week+precip.month+dd.year+ag_corn10+ag_beans10+ag_smgrains10+forest10+
                      s(lat,long),random=list(year=~1,site.id=~1),control = lmeControl(opt='optim'), data =max.rpadi,family=negative.binomial(1))
summary(mod.gam.rpadi$gam)

#rmadis max
mod.gam.rmaidis<-gamm(captures ~ precip.week+precip.month+dd.year+ag_corn10+ag_beans10+ag_smgrains10+forest10+
                      s(lat,long),random=list(year=~1,site.id=~1),control=lmeControl(opt='optim'), data =max.rmaidis,family=negative.binomial(1))
summary(mod.gam.rmaidis$gam)

#####PEAK TIMING######
mod.gam.top10.peak<-gamm(doy ~ precip.week+precip.month+dd.year+ag_corn75+ag_beans75+ag_smgrains75+forest75+
                      s(lat,long),random=list(year=~1,site.id=~1,species=~1), data =max.top10,family=gaussian)
summary(mod.gam.top10.peak$gam)

#aglyc max
mod.gam.aglyc.peak<-gamm(doy ~ precip.week+precip.month+dd.year+ag_corn10+ag_beans10+ag_smgrains10+forest10+
                      s(lat,long),random=list(year=~1,site.id=~1), data =max.aglyc,family=gaussian)
summary(mod.gam.aglyc.peak$gam)

#rpadi max
mod.gam.rpadi.peak<-gamm(doy ~ precip.week+precip.month+dd.year+ag_corn10+ag_beans10+ag_smgrains10+forest10+
                      s(lat,long),random=list(year=~1,site.id=~1), data =max.rpadi,family=gaussian)
summary(mod.gam.rpadi.peak$gam)

#rmadis max
mod.gam.rmaidis.peak<-gamm(doy ~ precip.week+precip.month+dd.year+ag_corn10+ag_beans10+ag_smgrains10+forest10+
                        s(lat,long),random=list(year=~1,site.id=~1), data =max.rmaidis,family=gaussian)
summary(mod.gam.rmaidis.peak$gam)

###figures
ggplot(max, aes(x = dd.year, y = captures))+
  geom_point(size=2)+
  stat_function(fun=function(x)exp(coef(mod.gam.max$gam)[1] + coef(mod.gam.max$gam)[4]*x),size=2,color="blue")+
  #annotate("text",label="p=0.002",x=.2,y=21000,size=5)+
  #annotate("text",label="paste(R ^ 2, \" = 0.29\")",x=.2,y=20000,parse = TRUE,size=5)+
  labs(x = "% Degree day accumulation", y = "Total aphid abundance\n(# of captures)")+
  theme(text = element_text(size=24),axis.text=element_text(color="black"),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line = element_line(size=.7, color="black"),legend.position="none")
