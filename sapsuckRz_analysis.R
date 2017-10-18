#sapsuckRz_anylsis

#load required packages
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
source("custom_functions.r")

##load in aphid data
#Using fread in the data.table package because this file is quite large
aphid<-fread("https://ndownloader.figshare.com/files/9465745")
aphid<-aphid[,-1]
colnames(aphid)<-c("Year","State","Site","Date","Sex","Species","Captures")

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

##Looking for phenology change
aphid.max<-as.data.frame(aphid.cut %>% group_by(Year,Site,Aphid.Species) %>% slice(which.max(Captures)))
aphid.first<-as.data.frame(aphid.cut %>% group_by(Year,Site,Aphid.Species) %>% slice(Captures!=0) %>% slice(which.min(Date)))

#looking across all aphid species ####WE WILL NEED THIS LATTER####
aphid.all<-aggregate(Captures~Year+Date+Site,data=aphid,sum)

##degree day accumulation##
#load in weather data
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

#getting 8 warnings. look for misbehaving rows
weather.cast[c(which(is.na(weather.cast), arr.ind=TRUE)[,1]),]

##precipitation data
#precipitation accumulation over the growing season
#starting March 1
weather.cast$precip.accum<-accum.precip.time(weather.cast$precip,weather.cast$doy,60)

##Combine the weather and aphid data##
names(aphid.all)<-c("year","doy","site.id","captures")#so they match
weather.cast$year<-isoyear(weather.cast$date)

data<-merge(aphid.all,weather.cast,by=c("site.id","year","doy"),all.x=T)

##Start building models##

#fit quadratic term
#data$dd.acum2<-(data$dd.acum)^2

#fit model
aphid.model<- lm(captures~dd.acum+I(dd.acum^2),data=data[which(data$year==2005 & data$site.id=="ACRE"),])
#generate range of 50 numbers starting from 30 and ending at 160
xx <- seq(6,2000, length=50)
plot(captures~dd.acum,pch=19,ylim=c(0,500),data=data[which(data$year==2005 & data$site.id=="ACRE"),])
lines(xx, predict(aphid.model, data.frame(dd.acum=xx)), col="red",lwd=4)

##
data$year<-as.factor(data$year)

##remove observations before April 15 and after September 30
data<-data[which(data$doy>=105 | data$doy<=273),]

##remove years and sites with less than 6 observations
data<-data %>% group_by(site.id,year) %>% filter(n()>= 7) %>% ungroup()
data<-data.table(data)#use data.table class

#peak function
argmax <- function(x, y, w=3, ...) {
  require(zoo)
  n <- length(y)
  y.smooth <- loess(y ~ x, span=.5,...)$fitted
  y.max <- rollapply(zoo(y.smooth), 2*w+1, max, align="center")
  delta <- y.max - y.smooth[-c(1:w, n+1-1:w)]
  i.max <- which(delta <= 0) + w
  list(x=x[i.max], i=i.max, y.hat=y.smooth)
}

##calculate peaks
peaks<- data[,list(peak=argmax(.SD$dd.acum,.SD$captures)$x),by=c("site.id","year")]
#Some site x years have multiple peaks. Need to work on this.

##add rainfall data

#find the week each peak was in
data$week<-isoweek(data$date)

weeks<-c()
for (i in 1:length(peaks$year)){
  #pull out data for that year
  year.subset<-data[which(data$year==peaks$year[i]),]
  #pull out data for that site
  year.site.subset<-year.subset[which(year.subset$site.id==peaks$site.id[i]),]
  #ditch rows that occur after the peak DD accum
  over.dd.subset<-year.site.subset[which(year.site.subset$dd.acum<peaks$peak[i]),]
  #the last week in the data is the week the peak occurs in
  weeki<-max(over.dd.subset$week)
  #concatenate this value to the weeks vector
  weeks<-c(weeks, weeki)
}
peaks$week<-weeks

#put it into our peak object
peaks$week<-weeks

#merge peaks and data file
data<-merge(data,peaks,by=c("site.id","year","week"),all=T)

##Add in landcover data##
cdl<-read.csv("Data/cdl.csv")
cdl.cut<-cdl[,c(2:3,42:50)]
colnames(cdl.cut)[1]<-"year"
colnames(cdl.cut)[2]<-"site.id"

#add in landscape data
cdl.cut$year<-as.factor(cdl.cut$year)
data$site.id<-as.factor(data$site.id)
cdl.cut<-aggregate(.~site.id+year,data=cdl.cut,mean)
data<-merge(data,cdl.cut,by=c("site.id","year"),all.x=T)

#adding lat long
coords<-read.csv("Data/counties+coordinates2.csv")
coords<-coords[,1:3]
colnames(coords)<-c("site.id","lat","long")
data<-merge(data,coords,by="site.id",all.x=T)

#rescale data
data.scale2<-scale(data[,c(9:12,14:22)],center=F,scale=T)
data.scale2<-cbind(data[,c(1:8,13,23:24)],data.scale2)

##Analyze some data!!!
mod1<-lmer(peak~precip.accum+forest+ag+(1|site.id)+(1|year),data=data.scale2)
Anova(mod1,type=3) #car

tiff("peak_precip.tiff",width = 4200, height = 4200, units = "px", res = 600)
ggplot(data, aes(x = precip.accum, y = peak))+
  geom_point(size=2)+
  geom_smooth(method="lm",size=2,color="red",se=F)+
  annotate("text",label="p>0.001",x=190,y=1720,size=5)+
  annotate("text",label="paste(R ^ 2, \" = 0.40\")",x=200,y=1800,parse = TRUE,size=5)+
  labs(x = "Precipitation accumulation (mm)", y = "Degree day at peak aphid abundance")+
  theme(text = element_text(size=24),axis.text=element_text(colour="black"),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line = element_line(size=.7, color="black"),legend.position="none")
dev.off()

##total aphid abundance
aphid.total<-aggregate(captures~year+site.id,data=aphid.all,sum)
weather.max<-aggregate(.~site.id+year,data=weather.cast,max)
data.total<-merge(aphid.total,weather.max,by=c("site.id","year"),all.x=T)
data.total<-merge(data.total,cdl.cut,by=c("site.id","year"),all.x=T)

#rescale variables because the orders of magnitude difference in scales is making glmer.nb angry
data.total<-merge(data.total,coords,by="site.id",all.x=T)
data.scale<-scale(data.total[,9:20],center=F,scale=T)
data.scale<-cbind(data.total[,1:8],data.scale,data.total[,21:22])

#model
mod.pos<-glmer(captures~dd.acum+precip.accum+forest+ag+(1|site.id)+(1|year),data=data.scale,family=poisson)
Anova(mod2,type=3)

#Function for checking for overdispersion
overdisp_fun <- function(model) {
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))
  rdf <- nrow(model.frame(model))-model.df
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

overdisp_fun(mod.pos)

##
mod2<-glmer.nb(captures~dd.acum+precip.accum+forest+ag_corn+ag_beans+ag_wheat+ag_smgrains+(1|site.id)+(1|year),data=data.scale)
Anova(mod2,type=3)

#plots
tiff("captures_dd.acum.tiff",width = 4200, height = 4200, units = "px", res = 600)
ggplot(data.scale, aes(x = dd.acum, y = captures))+
  geom_point(size=2)+
  stat_function(fun=function(x)exp(fixef(mod2)[1] + fixef(mod2)[2]*x),size=2,color="red")+
  annotate("text",label="p<0.001",x=.865,y=21000,size=5)+
  annotate("text",label="paste(R ^ 2, \" = 0.31\")",x=.865,y=20000,parse = TRUE,size=5)+
  labs(x = "Degree day accumulation", y = "Total aphid abundance\n(# of captures)")+
  theme(text = element_text(size=24),axis.text=element_text(color="black"),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line = element_line(size=.7, color="black"),legend.position="none")
dev.off()

#GAM to add lat long spatial smoother
mod.gam1<-gam(peak ~ precip.accum+s(lat,long)+s(site.id,bs="re")+s(year,bs="re"), data = data)
summary(mod.gam1)


mod.gam2<-gam(captures ~ precip.accum+dd.acum+ag10+forest10+s(lat,long)+s(year,bs="re"), data = data.total,family="nb")
summary(mod.gam2)

mod.gam3<-gam(captures ~ precip.accum+dd.acum+ag_corn10+ag_beans10+ag_smgrains10+forest10+s(lat,long)+s(year,bs="re"), data = data.total,family="nb")
summary(mod.gam3)

library(sandwich)

cov.mod.gam3 <- vcovHC(mod.gam3, type = "HC0")
stderr.mod.gam3 <- sqrt(diag(cov.mod.gam3))
r.est <- cbind(Estimate = coef(mod.gam3), "Robust SE" = stderr.mod.gam3,
               "pr(>|z|" = 2 * pnorm(abs(coef(mod.gam3)/stderr.mod.gam3), lower.tail = F),
               LL = coef(mod.gam3) - 1.96*stderr.mod.gam3,
               UL = coef(mod.gam3) + 1.96*stderr.mod.gam3)
r.est

mod.gam2<-gam(captures ~ precip.accum+dd.acum+ag_corn+ag_wheat+forest+s(lat,long)+s(year,bs="re"), data = data.total,family=nb())
summary(mod.gam2)

mod.gam3<-gam(captures ~ precip.accum+dd.acum+ag+forest+s(lat,long)+s(year,bs="re"), data = data.total,family="poisson")

mod.corn<-glmer.nb(captures~ag_corn++(1|site.id)+(1|year),data=data.scale)

tiff("captures_corn.tiff",width = 4200, height = 4200, units = "px", res = 600)
ggplot(data.scale, aes(x = ag_corn, y = captures))+
  geom_point(size=2)+
  stat_function(fun=function(x)exp(fixef(mod.corn)[1] + fixef(mod.corn)[2]*x),size=2,color="red")+
  annotate("text",label="p=0.002",x=.2,y=21000,size=5)+
  annotate("text",label="paste(R ^ 2, \" = 0.29\")",x=.2,y=20000,parse = TRUE,size=5)+
  labs(x = "% Landscape corn cover", y = "Total aphid abundance\n(# of captures)")+
  theme(text = element_text(size=24),axis.text=element_text(color="black"),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line = element_line(size=.7, color="black"),legend.position="none")
dev.off()

mod.wheat<-glmer.nb(captures~ag_wheat+(1|site.id)+(1|year),data=data.scale) #something weird going on here

tiff("captures_wheat.tiff",width = 4200, height = 4200, units = "px", res = 600)
ggplot(data.scale, aes(x = ag_wheat, y = captures))+
  geom_point(size=2)+
  stat_function(fun=function(x)exp(fixef(mod.wheat)[1] + fixef(mod.wheat)[2]*x),size=2,color="red")+
  annotate("text",label="p=0.01",x=.2,y=21000,size=5)+
  annotate("text",label="paste(R ^ 2, \" = 0.27\")",x=.3,y=20000,parse = TRUE,size=5)+
  labs(x = "% Landscape wheat cover", y = "Total aphid abundance\n(# of captures)")+
  theme(text = element_text(size=24),axis.text=element_text(color="black"),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line = element_line(size=.7, color="black"),legend.position="none")
dev.off()


#The spatial smoother doesn't seem to matter

logLik(mod.gam2)
logLik(mod.gam3)

stat <- as.numeric(2 * (logLik(mod.gam2) - logLik(mod.gam3)))
pchisq(stat, df = 33.97 - 25.16, lower.tail = FALSE)

#effects plots

sjp.setTheme(base = theme_bw()) 
             
sjp.glmer(mod2, 
          type = "fe", 
          sort = TRUE)


####Disaggregate species for Soybean Aphid (Aphis glycines) Cherry-Oat Aphid (Rhopalosiphum padi) Corn Leaf Aphid (Rhopalosiphum maidis)####

aphid.aglyc <- aggregate(Captures ~ Site + Date + Year, data = aphid[aphid$Aphid.Species == "Aphis.glycines"], sum)
aphid.rpadi <- aggregate(Captures ~ Site + Date + Year, data = aphid[aphid$Aphid.Species == "Rhopalosiphum.padi"], sum)
aphid.rmaidis <- aggregate(Captures ~ Site + Date + Year, data = aphid[aphid$Aphid.Species == "Rhopalosiphum.maidis"], sum)

colnames(aphid.aglyc) <- c("site.id","doy","year","captures")
colnames(aphid.rpadi) <- c("site.id","doy","year","captures")
colnames(aphid.rmaidis) <- c("site.id","doy","year","captures")

#aglyc
##total aphid abundance
aphid.total.aglyc<-aggregate(captures~year+site.id,data=aphid.aglyc,sum)
weather.max<-aggregate(.~site.id+year,data=weather.cast,max)
data.total.aglyc<-merge(aphid.total.aglyc,weather.max,by=c("site.id","year"),all.x=T)
data.total.aglyc<-merge(data.total.aglyc,cdl.cut,by=c("site.id","year"),all.x=T)

##rescale variables because the orders of magnitude difference in scales is making glmer.nb angry
data.total.aglyc<-merge(data.total.aglyc,coords,by="site.id",all.x=T)
data.scale.aglyc<-scale(data.total.aglyc[,9:51],center=F,scale=T)
data.scale.aglyc<-cbind(data.total.aglyc[,1:8],data.scale.aglyc,data.total.aglyc[,52:53])

#rpadi
##total aphid abundance
aphid.total.rpadi<-aggregate(captures~year+site.id,data=aphid.rpadi,sum)
weather.max<-aggregate(.~site.id+year,data=weather.cast,max)
data.total.rpadi<-merge(aphid.total.rpadi,weather.max,by=c("site.id","year"),all.x=T)
data.total.rpadi<-merge(data.total.rpadi,cdl.cut,by=c("site.id","year"),all.x=T)

##rescale variables because the orders of magnitude difference in scales is making glmer.nb angry
data.total.rpadi<-merge(data.total.rpadi,coords,by="site.id",all.x=T)
data.scale.rpadi<-scale(data.total.rpadi[,9:51],center=F,scale=T)
data.scale.rpadi<-cbind(data.total.rpadi[,1:8],data.scale.rpadi,data.total.rpadi[,52:53])

#rmaidis
##total aphid abundance
aphid.total.rmaidis<-aggregate(captures~year+site.id,data=aphid.rmaidis,sum)
weather.max<-aggregate(.~site.id+year,data=weather.cast,max)
data.total.rmaidis<-merge(aphid.total.rmaidis,weather.max,by=c("site.id","year"),all.x=T)
data.total.rmaidis<-merge(data.total.rmaidis,cdl.cut,by=c("site.id","year"),all.x=T)

##rescale variables because the orders of magnitude difference in scales is making glmer.nb angry
data.total.rmaidis<-merge(data.total.rmaidis,coords,by="site.id",all.x=T)
data.scale.rmaidis<-scale(data.total.rmaidis[,9:51],center=F,scale=T)
data.scale.rmaidis<-cbind(data.total.rmaidis[,1:8],data.scale.rmaidis,data.total.rmaidis[,52:53])

#Models! (The disaggregated ones)
##aglyc
mod.gam2.aglyc<-gam(captures ~ precip.accum+dd.acum+ag10+forest10+s(lat,long)+s(year,bs="re"), data = data.scale.aglyc,family="nb")
summary(mod.gam2.aglyc)

mod.gam3.aglyc<-gam(captures ~ precip.accum+dd.acum+ag_corn10+ag_beans10+ag_smgrains10+forest10+s(lat,long)+s(year,bs="re"), data = data.scale.aglyc,family="nb")
summary(mod.gam3.aglyc)

#rpadi
mod.gam2.rpadi<-gam(captures ~ precip.accum+dd.acum+ag10+forest10+s(lat,long)+s(year,bs="re"), data = data.scale.rpadi,family="nb")
summary(mod.gam2.rpadi)

mod.gam3.rpadi<-gam(captures ~ precip.accum+dd.acum+ag_corn10+ag_beans10+ag_smgrains10+forest10+s(lat,long)+s(year,bs="re"), data = data.scale.rpadi,family="nb")
summary(mod.gam3.rpadi)

#rmaidis
mod.gam2.rmaidis<-gam(captures ~ precip.accum+dd.acum+ag10+forest10+s(lat,long)+s(year,bs="re"), data = data.scale.rmaidis,family="nb")
summary(mod.gam2.rmaidis)

mod.gam3.rmaidis<-gam(captures ~ precip.accum+dd.acum+ag_corn10+ag_beans10+ag_smgrains10+forest10+s(lat,long)+s(year,bs="re"), data = data.scale.rmaidis,family="nb")
summary(mod.gam3.rmaidis)



####Relative abund. by species for Soybean Aphid (Aphis glycines) Cherry-Oat Aphid (Rhopalosiphum padi) Corn Leaf Aphid (Rhopalosiphum maidis)####

data.ratio.aglyc <- merge(data.scale.aglyc,data.scale[,c(1:3)],by = c("site.id","year"),suffixes = c("",".total"))
data.ratio.aglyc$ratio <- data.ratio.aglyc$captures/data.ratio.aglyc$captures.total

data.ratio.rpadi <- merge(data.scale.rpadi,data.scale[,c(1:3)],by = c("site.id","year"),suffixes = c("",".total"))
data.ratio.rpadi$ratio <- data.ratio.rpadi$captures/data.ratio.rpadi$captures.total

data.ratio.rmaidis <- merge(data.scale.rmaidis,data.scale[,c(1:3)],by = c("site.id","year"),suffixes = c("",".total"))
data.ratio.rmaidis$ratio <- data.ratio.rmaidis$captures/data.ratio.rmaidis$captures.total

#Models! (the ratio ones)
##aglyc
mod.gam2.aglyc.ratio<-gam(ratio ~ precip.accum+dd.acum+ag10+forest10+s(lat,long)+s(year,bs="re"), data = data.ratio.aglyc)
summary(mod.gam2.aglyc.ratio)

mod.gam3.aglyc.ratio<-gam(ratio ~ precip.accum+dd.acum+ag_corn10+ag_beans10+ag_smgrains10+forest10+s(lat,long)+s(year,bs="re"), data = data.ratio.aglyc)
summary(mod.gam3.aglyc.ratio)

#rpadi
mod.gam2.rpadi.ratio<-gam(ratio ~ precip.accum+dd.acum+ag10+forest10+s(lat,long)+s(year,bs="re"), data = data.ratio.rpadi)
summary(mod.gam2.rpadi.ratio)

mod.gam3.rpadi.ratio<-gam(ratio ~ precip.accum+dd.acum+ag_corn10+ag_beans10+ag_smgrains10+forest10+s(lat,long)+s(year,bs="re"), data = data.ratio.rpadi)
summary(mod.gam3.rpadi.ratio)

#rmaidis
mod.gam2.rmaidis.ratio<-gam(ratio ~ precip.accum+dd.acum+ag10+forest10+s(lat,long)+s(year,bs="re"), data = data.ratio.rmaidis)
summary(mod.gam2.rmaidis.ratio)

mod.gam3.rmaidis.ratio<-gam(ratio ~ precip.accum+dd.acum+ag_corn10+ag_beans10+ag_smgrains10+forest10+s(lat,long)+s(year,bs="re"), data = data.ratio.rmaidis)
summary(mod.gam3.rmaidis.ratio)


# Land Cover Non-linearity/Discontinuity Detection ####

mod.gam.detect <- gam(captures ~ precip.accum+dd.acum+ag_beans10+ag_smgrains10+forest10, data = data.total,family=nb())
data.total.detect <- subset(data.total, !is.na(ag_beans10))
data.total.detect$resid <- mod.gam.detect$y - mod.gam.detect$fitted.values
plot(data.total.detect$ag_corn10, data.total.detect$resid, xlab = "% Landcover Corn w/i 10km", ylab = "Residuals from Restricted Model")
abline(mean(data.total.detect$resid),mod.gam3$coefficients[4])

mod.gam.detect.aglyc <- gam(captures ~ precip.accum+dd.acum+ag_smgrains10+ag_corn10+forest10, data = data.total.aglyc,family=nb())
data.total.detect.aglyc <- subset(data.total.aglyc, !is.na(ag_beans10))
data.total.detect.aglyc$resid <- mod.gam.detect.aglyc$y - mod.gam.detect.aglyc$fitted.values
plot(data.total.detect.aglyc$ag_beans10[abs(data.total.detect.aglyc$resid)<1500], data.total.detect.aglyc$resid[abs(data.total.detect.aglyc$resid)<1500], xlab = "% Landcover soybeans w/i 10km", ylab = "Residuals from Restricted Model")
abline(mean(mod.gam3.aglyc$y - mod.gam3.aglyc$fitted.values + mod.gam3.aglyc$coefficients[5] * mod.gam3.aglyc$model$ag_beans10),mod.gam3.aglyc$coefficients[5], col = "red", lwd = 2)

mod.gam.detect.rpadi <- gam(captures ~ precip.accum+dd.acum+ag_beans10+ag_smgrains10+forest10, data = data.total.rpadi,family=nb())
data.total.detect.rpadi <- subset(data.total.rpadi, !is.na(ag_beans10))
data.total.detect.rpadi$resid <- mod.gam.detect.rpadi$y - mod.gam.detect.rpadi$fitted.values
plot(data.total.detect.rpadi$ag_corn10, data.total.detect.rpadi$resid, xlab = "% Landcover Corn w/i 10km", ylab = "Residuals from Restricted Model")
abline(mean(mod.gam3.rpadi$y - mod.gam3.rpadi$fitted.values + mod.gam3.rpadi$coefficients[4] * mod.gam3.rpadi$model$ag_corn10), mod.gam3.rpadi$coefficients[4])

mod.gam.detect.rmaidis <- gam(captures ~ precip.accum+dd.acum+ag_beans10+ag_smgrains10+forest10, data = data.total.rmaidis,family=nb())
data.total.detect.rmaidis <- subset(data.total.rmaidis, !is.na(ag_beans10))
data.total.detect.rmaidis$resid <- mod.gam.detect.rmaidis$y - mod.gam.detect.rmaidis$fitted.values
plot(data.total.detect.rmaidis$ag_corn10, data.total.detect.rmaidis$resid, xlab = "% Landcover Corn w/i 10km", ylab = "Residuals from Restricted Model")
abline(mean(mod.gam3.rmaidis$y - mod.gam3.rmaidis$fitted.values + mod.gam3.rmaidis$coefficients[4] * mod.gam3.rmaidis$model$ag_corn10), mod.gam3.rmaidis$coefficients[4], col = "red", lwd = 2)

##########################################
cc <- confint(mod1,parm="beta_")
ctab <- as.data.frame(cbind(est=fixef(mod1),cc))
colnames(ctab)<-c("est","lower","upper")

tiff("peak_effects.tiff",width = 4200, height = 4200, units = "px", res = 600)
ggplot(ctab,aes(x=rownames(ctab),y=est))+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,position=position_dodge(0.1))+
  scale_x_discrete(labels=c("Intercept", "Landscape\ncrop cover","Landscape\nforest cover","Precipitation\naccumulation"))+
  geom_hline(yintercept=0,linetype="dashed")+
  labs(x = "", y = "Effect size")+
  ggtitle("Peak aphid abundance")+
  coord_flip()+
  theme(text = element_text(size=24),axis.text=element_text(color="black"),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line = element_line(size=.7, color="black"),legend.position="none")
dev.off()

dd <- confint(mod2,parm="beta_",method="Wald")
dtab <- as.data.frame(cbind(est=fixef(mod2),dd))
colnames(dtab)<-c("est","lower","upper")

tiff("abund_effects.tiff",width = 4200, height = 4200, units = "px", res = 600)
ggplot(dtab,aes(x=rownames(dtab),y=est))+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,position=position_dodge(0.1))+
  scale_x_discrete(labels=c("Intercept", "Landscape\ncrop cover","Degree day\naccumulation","Landscape\nforest cover","Precipitation\naccumulation"))+
  geom_hline(yintercept=0,linetype="dashed")+
  labs(x = "", y = "Effect size")+
  ggtitle("Total aphid abundance")+
  coord_flip()+
  theme(text = element_text(size=24),axis.text=element_text(color="black"),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line = element_line(size=.7, color="black"),legend.position="none")
dev.off()

