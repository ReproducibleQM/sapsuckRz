#sapsuckRz_anylsis

#load required packages
library(data.table)
library(data.table)
library(reshape2)
library(Hmisc)

##load in aphid data
#Using fread in the data.table package because this file is quite large
aphid<-fread("https://dl.dropboxusercontent.com/u/98197254/Mastersuction_may29_2014.csv")

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
aphid.cast<-dcast(aphid.cut,time+Site~Aphid.Species,value.var="Captures",sum)

#corelation between species
rcorr(as.matrix(aphid.cast[3:8]))

#remove rows when the two most common species are zero (no aphids out at this time?)
aphid.cast.cut<-aphid.cast[aphid.cast$Rhopalosiphum.maidis!=0 & aphid.cast$Rhopalosiphum.padi!=0,]

#Check cor again
rcorr(as.matrix(aphid.cast.cut[3:8]))

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