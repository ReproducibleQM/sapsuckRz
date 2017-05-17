#sapsuckRz_anylsis

#load required packages
library(data.table)
library(ggplot2)
library(vegan)

##load in aphid data
#Using fread in the data.table package because this file is quite large
aphid<-fread("https://dl.dropboxusercontent.com/u/98197254/Mastersuction_may29_2014.csv")

##Pull the 10 most common species
#find the average number of captures for each species
aphid.ag<-aggregate(Captures~Aphid.Species,data=aphid,sum)
#order data from most number of captures to least
aphid.ord<-aphid.ag[order(-aphid.ag$Captures),] 
#make a list of the 10 most common species
aphid.names.cut<-aphid.ord$Aphid.Species[1:25]
#cut down the data to only include those 10 species
aphid.cut<-aphid[aphid$Aphid.Species%in%aphid.names.cut]

#add a (crude) time variable
aphid.cut$time<-((aphid.cut$Year-2005)*365)+aphid.cut$Date

####NMDS####
aphid.cast<-dcast(aphid.cut,Year+State~Aphid.Species,value.var="Captures",sum)
rownames(aphid.cast)<-paste(aphid.cast$State,aphid.cast$Year,sep="_")
aphid.year<-as.factor(aphid.cast$Year)
aphid.state<-as.factor(aphid.cast$State)
aphid.cast$Year<-NULL
aphid.cast$State<-NULL

nmds.aphid<-metaMDS(aphid.cast,k=2,trymax=250,distance="bray")

adonis(vegdist(aphid.cast)~aphid.year,permutations=10000)
sim.year<-simper(aphid.cast,aphid.year)
summary(sim.year,ordered=TRUE) #What's driving this? looks like 2005

adonis(vegdist(aphid.cast)~aphid.state,permutations=10000)
sim.state<-simper(aphid.cast,aphid.state)
summary(sim.state,ordered=TRUE) #What's driving this? Looks like Iowa

##plotting
tiff("Aphid_NMDS_year.tiff",width = 3200, height = 3200, units = "px", res = 800)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
ordiplot(nmds.aphid,type="n",main="Aphid composition by year")
points(nmds.aphid,display="sites",pch=16,cex=1.25,col=c("#e41a1c","lightblue","#4daf4a","#984ea3","#ff7f00","darkblue","#a65628","#f781bf","#999999")[aphid.year])
#ordiellipse(nmds.aphid, aphid.year ,col=c("#e41a1c","lightblue","#4daf4a","#984ea3","##ff7f00","darkblue","#a65628","#f781bf","#999999"))
legend(x=.85,y=1,expression("2005","2006","2007","2008","2009","2010","2011","2012","2013")
       ,pch=16,col=c("#e41a1c","lightblue","#4daf4a","#984ea3","#ff7f00","darkblue","#a65628","#f781bf","#999999"),cex=.95)
dev.off()

tiff("Aphid_NMDS_state.tiff",width = 3200, height = 3200, units = "px", res = 800)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
ordiplot(nmds.aphid,type="n",main="Aphid composition by state")
points(nmds.aphid,display="sites",pch=16,cex=1.25,col=c("#e41a1c","lightblue","#4daf4a","#984ea3","#ff7f00","darkblue","#a65628","#f781bf","#999999","black")[aphid.state])
legend(x=.85,y=1,expression("Iowa","Illinois","Indiana","Michigan","Minnesota","Wisconsin","Kansas","Kentucky","Missouri","South Dakota")
       ,pch=16,col=c("#e41a1c","lightblue","#4daf4a","#984ea3","#ff7f00","darkblue","#a65628","#f781bf","#999999","black"),cex=.85)
dev.off()

##Looking for phenology change
aphid.max<-as.data.frame(aphid.cut %>% group_by(Year,Site,Aphid.Species) %>% slice(which.max(Captures)))
aphid.first<-as.data.frame(aphid.cut %>% group_by(Year,Site,Aphid.Species) %>% slice(Captures!=0) %>% slice(which.min(Date)))

#Number of captures through time across all sites
ggplot(aphid.cut,aes(time,Captures,colour=Aphid.Species))+
  geom_line(size=1)+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf,size=1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=1)+
  theme(axis.line.x = element_line(color = "black",size=1), axis.line.y = element_line(color = "black",size=1),
        text = element_text(size=14),axis.text=element_text(colour="black"),
        panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line = element_line(size=.7, color="black"),
        strip.background = element_rect(colour="white", fill="white"),legend.position="top")

aphid.w.aphgly<-aphid.cut[aphid.cut$Aphid.Species=="Aphis.glycines"]
#facet for each state
p1<-ggplot(aphid.w.aphgly,aes(Date,Captures,colour=Aphid.Species))+
  geom_smooth(size=1.25,method="gam",formula=y ~ poly(x, 10),se=F)+
  facet_wrap(~Year)+
  labs(x="Day of Year", y="Aphid Captures")+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf,size=1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=1)+
  coord_cartesian(ylim = c(0, 500))+
  theme(axis.line.x = element_line(color = "black",size=1), axis.line.y = element_line(color = "black",size=1),
        text = element_text(size=14),axis.text=element_text(colour="black"),
        panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line = element_line(size=.7, color="black"),
        strip.background = element_rect(colour="white", fill="white"),legend.position="top")

tiff("Aphid_abund_year_soy.tiff",width = 3200, height = 3200, units = "px", res = 800)
p1
dev.off()

#Remove soybean aphid and replot (easier to see variation in other species)
aphid.cut.aphgly<-aphid.cut[aphid.cut$Aphid.Species!="Aphis.glycines"]

#soybean aphid removed
p2<-ggplot(aphid.cut.aphgly,aes(Date,Captures,colour=Aphid.Species))+
  geom_smooth(size=1.25,method="gam",formula=y ~ poly(x, 10),se=F)+
  facet_wrap(~Year)+
  labs(x="Day of Year", y="Aphid Captures")+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf,size=1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=1)+
  coord_cartesian(ylim = c(0, 100))+
  scale_colour_manual(values=c("#e41a1c","lightblue","#4daf4a","#984ea3","#ff7f00","darkblue","#a65628","#f781bf","#999999"))+
  theme(axis.line.x = element_line(color = "black",size=1), axis.line.y = element_line(color = "black",size=1),
        text = element_text(size=14),axis.text=element_text(colour="black"),
        panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line = element_line(size=.7, color="black"),
        strip.background = element_rect(colour="white", fill="white"),legend.position="top")

tiff("Aphid_abund_year.tiff",width = 8000, height = 4000, units = "px", res = 600)
p2
dev.off()
#facet for each state with soybean aphid removed
ggplot(aphid.cut.aphgly,aes(time,Captures,colour=Aphid.Species))+
  geom_smooth(method="gam",size=1)+
  facet_wrap(~State)+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf,size=1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=1)+
  theme(axis.line.x = element_line(color = "black",size=1), axis.line.y = element_line(color = "black",size=1),
        text = element_text(size=14),axis.text=element_text(colour="black"),
        panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line = element_line(size=.7, color="black"),
        strip.background = element_rect(colour="white", fill="white"),legend.position="top")