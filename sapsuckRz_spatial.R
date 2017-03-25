## Extracting landcover data of a 1.5km (or variable) radius surrounding United States aphid traps
## Original raster rearing by Christie Bahlai, reproducible script prepared by Braeden Van Deynze
## Pushed 3/24/2017

# let's see if we can spatial statistic! 

#going to try pulling in the cropland data layer, extracting information about
#crops surrounding a relevant point

#prepare relavent libraries

#install.packages("cdlTools")
#install.packages("raster")

library(data.table)
library(cdlTools)
library(sp)
library(rgeos)
library(rgdal)
library(raster)
library(reshape2)
library(reshape)

#create "sites.spatial" dataframe to dump data in
sites.spatial <- t(data.frame(varNamesCDL))
colnames(sites.spatial) = sites.spatial[1,]
sites.spatial = sites.spatial[-1,]
sites.spatial = sites.spatial[,-seq(2, NCOL(sites.spatial), by = 2)]
meta <- read.table(header = TRUE, stringsAsFactors = FALSE, text = 
                     "xSITEID xSTATE xYEAR xtime")
sites.spatial <- merge(meta,sites.spatial)
rm(meta)
sites.spatial <- data.frame(lapply(sites.spatial, function(y) as.numeric(as.character(y))))
names(sites.spatial) <- substring(names(sites.spatial), 2)
sites.spatial$SITEID <- as.character(sites.spatial$SITEID)
sites.spatial$STATE <- as.character(sites.spatial$STATE)

#bring in aphid data
aphid<-fread("https://dl.dropboxusercontent.com/u/98197254/Mastersuction_may29_2014.csv")

# years to query
years <- unique(aphid, by="Year")
years <- as.numeric(years[,Year])
years <- years[years > 2007]
state <- unique(aphid, by="State")
state <- state[,State]

# prepare a list to store all the state-year rasters
uberraster <- vector("list",length(state))


for(i in 1:length(state)){
  Sys.time() ##timestamp so you know when each state started
  uberraster[[i]] <- getCDL(state[i],years)
}
#we have data! FANTASTIC!

#now we have to create buffer crop share data for each site
# first prepare the coordinates
coord <- fread("https://raw.githubusercontent.com/ReproducibleQM/sapsuckRz/master/counties%2Bcoordinates2.csv")
state.num <- data.table(state)
state.num$id <- rownames(state.num)
coord <- merge(coord, state.num, by.x = "STATE", by.y = "state")

## DONT NEED THIS next prepare the data dump with empty values for everything except for site, state, and year
#bigcoord <- coord
#bigcoord$LAT <- NULL
#bigcoord$LONG <- NULL
#bigcoord$COUNTY <- NULL
#bigcoord$id <- NULL
#
#ubercoord <- vector("list", length(years))
#for(i in 1:length(years)){
#  ubercoord[[i]] <- bigcoord
#  ubercoord[[i]]$YEAR <- years[i]
#}
#
#for(i in 1:length(years)){
#  sites.spatial <- merge(sites.spatial,ubercoord[[i]], all = TRUE)
#}

# finally, pull the data and dump it in the data dump!

for(j in 1:nrow(coord)){
  
  for(k in 1:length(years)){
    l <- as.numeric(coord[j]$id)
    
    cdl<-uberraster[[l]][[k]]
    lat = coord$LAT[j]
    lon = coord$LONG[j]
    SITEID = coord$SITEID[j]
    YEAR = years[k]
    STATE = coord$STATE[j]
    local.pts <- data.frame(lat, lon, ID)
    coordinates(local.pts)<-~lon+lat
    proj4string(local.pts)<- CRS("+proj=longlat +datum=WGS84")
    pts <- spTransform(local.pts, proj4string(cdl))
    buf <- gBuffer(pts, width=1500, byid = TRUE)
    landcover.bits<-crop(cdl, buf)
    landcover.bits<-mask(landcover.bits,buf)
    values<-getValues(landcover.bits)
    values<-data.frame(table(values))
    values$Perc <- round(100 * (values$Freq/sum(values$Freq)),1)
    values.t <- t(values[,-2])
    colnames(values.t) = values.t[1,]
    values.t <- values.t[-1,]
    
    values.t$XYEAR <- YEAR
    values.t$XSTATE <- STATE
    values.t$XSITEID <- SITEID
    
    values.s <- data.frame(values.t)
    names(values.s) <- substring(names(values.s), 2)
    
    sites.spatial <- merge(sites.spatial,values.s,all = TRUE )
    }
}

#Now we prepare the data product and export!

# order columns logically

refcols <- c("SITEID","YEAR","STATE","time")
sites.spatial.exp <- sites.spatial[, c(refcols, setdiff(names(sites.spatial), refcols))]

# zeros for NA, correct for class

sites.spatial.exp[c(-1,-2,-3,-4)] <- lapply(sites.spatial.exp[c(-1,-2,-3,-4)], function(x) {
  if(is.factor(x)) as.numeric(as.character(x)) else x
})
sapply(sites.spatial.exp, class)

sites.spatial.exp[c(-1,-2,-3,-4)] <- lapply(sites.spatial.exp[c(-1,-2,-3,-4)], function(x) {
  if(is.character(x)) as.numeric(x) else x
})
sapply(sites.spatial.exp, class)

sites.spatial.exp[is.na(sites.spatial.exp)]<-0
sites.spatial.exp <- sites.spatial.exp[c(-4)]
View(sites.spatial.exp)
sapply(sites.spatial.exp, class)

# approve and export to .dir

write.csv(sites.spatial.exp, file = "sites_spatial.csv")
print("Congrats! Enjoy your data!")
