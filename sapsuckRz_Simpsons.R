#Load required packages
##install.packages("data.table","reshape") <- INSTALL THESE PACKAGES IF YOU DO NOT ALREADY HAVE THEM
library(data.table)
library(reshape)

#Load base aphid trap data using fread
aphid<-fread("https://dl.dropboxusercontent.com/u/98197254/Mastersuction_may29_2014.csv")

#Load site spatial data
spatial<-fread("https://raw.githubusercontent.com/ReproducibleQM/sapsuckRz/7dc2fca26288ddba25422c539fb77e35eca7ce09/counties%2Bcoordinates2.csv")

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

#Cast data to better suit diversity scoring at site by time level
aphid.melt <- melt(aphid.cut, id.vars = c("Site", "time", "Aphid.Species", "Sex"), measure.vars = c("Captures"))
sites <- cast(aphid.melt, Site + time ~ Aphid.Species, sum)

#Merge in state and county data by site id
sites.geo <- merge(sites, spatial, by.x = "Site", by.y = "SITEID")

#Create a measure of Simpson's Diversity Index (http://rfunctions.blogspot.com/2012/02/diversity-indices-simpsons-diversity.html)
simpson.diversity <- numeric(nrow(sites))
for(i in 1:nrow(sites)){
  soma <- sum(sites[i,3:12])
  prop <- numeric(10)
  prop2 <- prop
  for(j in 1:10){
    prop[j] <- sites[i,j+2]/soma
    prop2[j] <- prop[j]^2
  }
  simpson.diversity[i] <- 1/sum(prop2)
}

#Merge Simpson values into site data
simp <- data.frame(simpson.diversity)
sites.simp <- sites.geo
sites.simp$Simpsons <- simp$simpson.diversity

