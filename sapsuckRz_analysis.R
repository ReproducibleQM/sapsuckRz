#sapsuckRz_anylsis

#load required packages
library(data.table)

#load in aphid data
#Using fread in the data.table package because this file is quite large
aphid<-fread("http://lter.kbs.msu.edu/datatables/122.csv")
#test