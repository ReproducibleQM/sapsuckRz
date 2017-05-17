##Download weather data from PRISM and pull weather data for individual sites
#created on:5/17/2017
#created by: Chad Zirbel

#pull in coordinates for all sites
location<-read.csv("counties+coordinates2.csv")

#install and load packages
#install.packages("prism")
library(prism)


#pull max and min temp and precipitation
type=c("tmax","tmin","ppt")

#download daily prism data from 2003-2016
for(i in 1:3){
  var <- paste("weather", i, sep = "")
  assign(var, get_prism_dailys(type=type[i], minDate = "2003-01-01", maxDate = "2016-01-01", keepZip=T))
}

#function to pull data from specific locations
prism_cut <- function(location,prismfile){
  meta_d <- unlist(prism_md(prism_data[,1],returnDate=T))
  meta_names <- as.matrix(unlist(prism_md(prism_data[,1])))
  param_name <- do.call(rbind,strsplit(meta_names,"-"))[,3]

  pstack <- prism_stack(prism_data[,1])
  data <- unlist(extract(pstack,cbind(location$LONG,location$LAT),nrow=1,buffer=10))
  data <- as.data.frame(data)
  data$date <- as.Date(meta_d)
  data$param<-param_name
  data$siteid<-location$SITEID

  return(data)} 

#load in the prism data files
prism_data<-ls_prism_data(absPath=T)

#pull the data
weather_data<-prism_cut(location,prism_data[,1])
             