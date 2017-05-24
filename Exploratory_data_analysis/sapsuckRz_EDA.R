#sapsuckRz_anylsis

#load required packages
library(data.table)
library(reshape2)
library(Hmisc)

##load in aphid data
#Using fread in the data.table package because this file is quite large
aphid<-fread("https://dl.dropboxusercontent.com/u/98197254/Mastersuction_may29_2014.csv")
counties<-read.csv("counties+coordinates.csv")
crop<-read.csv("counties+acres.csv")

##Pull wheat pests
#list of wheat aphids
aphid.names.cut<-c("Schizaphis.graminum","Diuraphis.noxia","Sitobion.avenae","Rhopalosiphum.padi","Rhopalosiphum.maidis","Metopolophium.dirhodum")

#cut down the data to only include those 10 species
aphid.cut<-aphid[aphid$Aphid.Species%in%aphid.names.cut]

#csv to put on drive to avoid using the entire dataset
#write.csv(aphid.cut,"aphid_wheat.csv")

#add a (crude) time variable
aphid.cut$time<-((aphid.cut$Year-2005)*365)+aphid.cut$Date

#reshape to wide
aphid.cast<-dcast(aphid.cut,time+Site+Year~Aphid.Species,value.var="Captures",sum)

#corelation between species
rcorr(as.matrix(aphid.cast[3:8]))

#remove rows when the two most common species are zero (no aphids out at this time?)
aphid.cast.cut<-aphid.cast[aphid.cast$Rhopalosiphum.maidis!=0 & aphid.cast$Rhopalosiphum.padi!=0,]

#Check cor again
rcorr(as.matrix(aphid.cast.cut[3:8]))

##control for number of acres plotted in crops for each county in each year. Then look at correlations
#rename variables so they match
names(counties)<-c("Site","lat","long","county","state")
names(crop)<-c("Year","state","county","crop","item","value")

#merge the trap location and crop data
crop.county<-merge(counties,crop, by="county",all=T)
crop.cut<-crop.county[,c(2,6,10)]

#data cleaning so files merge correctly
#change to uppercase and delete spaces
aphid.cast$Site<-toupper(aphid.cast$Site)
aphid.cast$Site<-gsub(" ","",aphid.cast$Site)
crop.cut$Site<-gsub(" ","",crop.cut$Site)

#fix spelling of Bean & Beet
crop.cut$Site<-gsub("BEAT","BEET",crop.cut$Site)

#merge crop and aphid data
aphid.crop<-merge(aphid.cast,crop.cut,by=c("Site","Year"),all=T)
aphid.crop$value<-as.numeric(aphid.crop$value)
aphid.crop<-na.omit(aphid.crop)

#models (abundance of each species predicting by amount of crop planted in that county)
mod.dn<-lm(Diuraphis.noxia~value,data=aphid.crop)
mod.md<-lm(Metopolophium.dirhodum~value,data=aphid.crop)
mod.rm<-lm(Rhopalosiphum.maidis~value,data=aphid.crop)
mod.rp<-lm(Rhopalosiphum.padi~value,data=aphid.crop)
mod.sg<-lm(Schizaphis.graminum~value,data=aphid.crop)
mod.sa<-lm(Sitobion.avenae~value,data=aphid.crop)

#add residuals to dataframe
aphid.crop$dn.r<-resid(mod.dn)
aphid.crop$md.r<-resid(mod.md)
aphid.crop$rm.r<-resid(mod.rm)
aphid.crop$rp.r<-resid(mod.rp)
aphid.crop$sg.r<-resid(mod.sg)
aphid.crop$sa.r<-resid(mod.sa)

#check for correlations between residuals of each of the species controling for crop acreage
rcorr(as.matrix(aphid.crop[11:16]))

#now aggregate by year (correlations on different temporal scale)
aphid.crop.ag<-aggregate(.~Year+Site,data=aphid.crop, sum)
rcorr(as.matrix(aphid.crop.ag[4:9]))

#models (abundance of each species predicting by amount of crop planted in that county)
mod.dn.ag<-lm(Diuraphis.noxia~value,data=aphid.crop.ag)
mod.md.ag<-lm(Metopolophium.dirhodum~value,data=aphid.crop.ag)
mod.rm.ag<-lm(Rhopalosiphum.maidis~value,data=aphid.crop.ag)
mod.rp.ag<-lm(Rhopalosiphum.padi~value,data=aphid.crop.ag)
mod.sg.ag<-lm(Schizaphis.graminum~value,data=aphid.crop.ag)
mod.sa.ag<-lm(Sitobion.avenae~value,data=aphid.crop.ag)

#add residuals to dataframe
aphid.crop.ag$dn.r<-resid(mod.dn.ag)
aphid.crop.ag$md.r<-resid(mod.md.ag)
aphid.crop.ag$rm.r<-resid(mod.rm.ag)
aphid.crop.ag$rp.r<-resid(mod.rp.ag)
aphid.crop.ag$sg.r<-resid(mod.sg.ag)
aphid.crop.ag$sa.r<-resid(mod.sa.ag)

#look at cor between resids
rcorr(as.matrix(aphid.crop.ag[11:16]))

##plotting

#Number of captures through time across all sites
ggplot(aphid.cut,aes(Date,Captures,colour=Aphid.Species))+
  geom_line(size=1)+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf,size=1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=1)+
  theme(axis.line.x = element_line(color = "black",size=1), axis.line.y = element_line(color = "black",size=1),
        text = element_text(size=14),axis.text=element_text(colour="black"),
        panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line = element_line(size=.7, color="black"),
        strip.background = element_rect(colour="white", fill="white"),legend.position="top")

#facet for each state
ggplot(aphid.cut,aes(time,Captures,colour=Aphid.Species))+
  geom_line(size=1)+
  facet_wrap(~State)+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf,size=1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=1)+
  theme(axis.line.x = element_line(color = "black",size=1), axis.line.y = element_line(color = "black",size=1),
        text = element_text(size=14),axis.text=element_text(colour="black"),
        panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line = element_line(size=.7, color="black"),
        strip.background = element_rect(colour="white", fill="white"),legend.position="top")