# let's see if we can spatial statistic! 

#going to try pulling in the cropland data layer, extracting information about
#crops surrounding a relevant point


# install.packages(
#   "cdlTools",
#   "raster",
#   "sp",
#   "rgeos",
#   "rgdal",
#   "mailR"
# )

library(cdlTools)
library(raster)
library(sp)
library(rgeos)
library(rgdal)
library(mailR)

crop_data <- data.frame()

##draw aphid data for year and state identification
# bring in aphid data
aphid<-fread("https://dl.dropboxusercontent.com/u/98197254/Mastersuction_may29_2014.csv")

# years, states to query
years <- unique(aphid, by="Year")
years <- as.numeric(years[,Year])
years <- years[years > 2007]
states <- unique(aphid, by="State")
states <- states[,State]

# remove aphid data from memory
rm(aphid)

##draw raster data statewise, stacked by year
for(k in 1:length(states)){
state <- states[k]
cdl.raster <- getCDL(state,years)
# we have data! FANTASTIC!

# lets see how it looks!!


##okay now, let's chose a point in this map and see what crops are around it

# load location data
location<-read.csv("https://raw.githubusercontent.com/ReproducibleQM/sapsuckRz/master/counties%2Bcoordinates2.csv")
location$SITEID <- as.character(location$SITEID)
location$COUNTY <- as.character(location$COUNTY)
location$STATE <- as.character(location$STATE)

# reduce to relevent state
location <- location[!(location$STATE != state),]

for(j in 1:length(years)){
year= years[j]
cdl <- cdl.raster[[as.character(year)]]

for(i in 1:nrow(location)){
  # find site on map
lat = location$LAT[i]
lon = location$LONG[i]
ID= location$SITEID[i]

points <- data.frame(lat, lon, ID)
coordinates(points)<-~lon+lat
proj4string(points)<- CRS("+proj=longlat +datum=WGS84") #let R know the xy data is lat long

points<-spTransform(points, proj4string(cdl))

# now create a buffer

buf<-gBuffer(points, width=75000, byid = TRUE) #using units of a meter- let's get a 75 km radius
                                        #byid= true so each buffer will create its own polygon


# okay, let's see if we can use our polygons to snip out the data we're interested in
# it'll speed up our processing later
landcover.bits<-crop(cdl, buf) #snips things to outer extents of polygons
landcover.bits<-mask(landcover.bits, buf)

# and let's pull the summary data out now
# below puts data in a nice data table

values<-getValues(landcover.bits)
values<-data.frame(table(values))
values$Perc <- round(100 * (values$Freq/sum(values$Freq)),1)
values$Freq <- NULL
values$values <- as.character(values$values)
values <- t(values)
colnames(values) <- values["values",]
values <- values[-c(1),]
line <- c(ID, state, year, values)
names(line) <- c("ID","State","Year",names(values))
df = data.frame(line)
df <- t(df)
crop_data <- merge(crop_data,df,all=TRUE)
print(year)
print(ID)
rm(buf)
rm(landcover.bits)
}
rm(cdl)
}
rm(cdl.raster)
}

##Now we need to export our data (above pull will take a lot of time, advised to execute on a dedicated machine)
write.csv(crop_data, file = "crop_data.csv", row.names = FALSE) #saves data in csv

# send mail alerts for completion

send.mail(from = "sapsuckRz@gmail.com",
          to = "bsvandeynze@gmail.com",
          subject = "Crop Data for Aphid Project",
          body = "I sent this from R!",
          smtp = list(host.name="smtp.gmail.com",
                      port=465,
                      ssl=TRUE,
                      user.name="sapsuckRz",
                      passwd="aphidsrule"),
          authenticate = TRUE,
          attach.files = "crop_data.csv")
