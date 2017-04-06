## Merging land cover data with panel site data
## Aphid list reduction by Chad Zirbel (crzirbel@gmail.com), merge script by Braeden Van Deynze (vandeynz@gmail.com)

library(data.table)
library(reshape)
library(vegan)

#Load base aphid trap data using fread
aphid<-fread("https://dl.dropboxusercontent.com/u/98197254/Mastersuction_may29_2014.csv")

#Pull the ten most common species
##Find the average number of captures for each species
aphid.ag<-aggregate(Captures~Aphid.Species,data=aphid,mean)
##Order data from most number of captures to least
aphid.ord<-aphid.ag[order(-aphid.ag$Captures),] 
##Make a list of the ten most common species
aphid.names.cut<-aphid.ord$Aphid.Species[1:10]
##Cut down the data to only include those ten species
aphid.cut<-aphid[aphid$Aphid.Species%in%aphid.names.cut]

#Add a (crude) time variable
aphid.cut$time<-((aphid.cut$Year-2005)*365)+aphid.cut$Date

#Cast data to better suit site by time level
aphid.melt <- melt(aphid.cut, id.vars = c("Site", "time", "Aphid.Species", "Sex", "Year", "State"), measure.vars = c("Captures"))
sites <- cast(aphid.melt, State + Year + Site + time ~ Aphid.Species, sum)

rm(aphid.ag, aphid.cut, aphid.melt, aphid.names.cut, aphid.ord)

#Bring in land cover data and merge
crops <- fread("https://raw.githubusercontent.com/ReproducibleQM/sapsuckRz/master/sites_spatial.csv", header = TRUE)
crops$V1 <- NULL
sites.crops <- merge(sites, crops, all=TRUE, by.x = c("Site","Year","State"), by.y = c("SITEID","YEAR","STATE"))
sites.crops <- sites.crops[sites.crops$Year > 2007 & sites.crops$time != "",]

#sites.crops is your finished dataframe, prepared for panel level analysis.
#Note that crop cover data only available for 2008- and some sites have missing years. Therefore, the sample size reduction when crop data is pulled in.