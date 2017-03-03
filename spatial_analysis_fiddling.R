# let's see if we can spatial statistic! 

#going to try pulling in the cropland data layer, extracting information about
#crops surrounding a relevant point




library(cdlTools)

# years to query
years <- 2009
state <- 'maryland'
cdl.raster <- getCDL(state,years)
#we have data! FANTASTIC!

#lets see how it looks!!
library(raster)

plot(cdl.raster$`2009`)
#Cool!

#let's create a shorter name for this object for subsequent analyses
cdl<-cdl.raster$`2009`

library(sp)
library(rgeos)
library(rgdal)

#okay now, let's chose a point in this map and see what crops are around it
lat = 39.2904
lon = -76.6122
ID= "Baltimore"

city.points <- data.frame(lat, lon, ID)
coordinates(city.points)<-~lon+lat
proj4string(city.points)<- CRS("+proj=longlat +datum=WGS84") #let R know the xy data is lat long

points<-spTransform(city.points, proj4string(cdl))# zone 18 is best for Maryland
  #use UTM zone 16 for suction traps...some fall in 15 but it's close so if we do this by state, change projection
  #good ref for choosing a projection to use is here:
  #https://www.nceas.ucsb.edu/~frazier/RSpatialGuides/OverviewCoordinateReferenceSystems.pdf

#now create a buffer

buf<-gBuffer(points, width=10000, byid = TRUE) #using units of a meter- let's get a 10 km radius
                                        #byid= true so each buffer will create its own polygon

plot(buf)
buf
summary(buf)
buf@data

#okay, let's see if these two layers will work together

plot(cdl)
plot(buf, add=T, border="red")

#okay, let's see if we can use our polygons to snip out the data we're interested in
#it'll speed up our processing later
landcover.bits<-crop(cdl, buf) #snips things to outer extents of polygons
landcover.bits<-mask(landcover.bits, buf)
#plot it to see if we got what we wanted
plot(landcover.bits)

#and let's pull the summary data out now

#this produces a table with the frequency of each observation by landcover type. 
# a key to the numerical values is given in https://cran.r-project.org/web/packages/cdlTools/cdlTools.pdf



values<-getValues(landcover.bits)
values<-data.frame(table(values))
values$Perc <- round(100 * (values$Freq/sum(values$Freq)),1)

