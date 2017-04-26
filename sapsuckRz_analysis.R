#sapsuckRz_anylsis

#load required packages
library(data.table)

#load required packages
library(data.table)

##load in aphid data!
#Using fread in the data.table package because this file is quite large
aphid<-fread("https://dl.dropboxusercontent.com/u/98197254/Mastersuction_may29_2014.csv")

##Pull wheat pests
#list of wheat aphids
aphid.names.cut<-c("Schizaphis.graminum","Diuraphis.noxia","Sitobion.avenae","Rhopalosiphum.padi","Rhopalosiphum.maidis","Metopolophium.dirhodum")

#cut down the data to only include those 10 species
aphid.cut<-aphid[aphid$Aphid.Species%in%aphid.names.cut]

write.csv(aphid.cut,"aphid_wheat.csv")
