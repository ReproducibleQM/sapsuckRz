# CDL variable rename, data reduction, and data merge code snip
# Bvd, July 2017
# This snip of code will sum similar landcover types, relabel the variables with descriptive names, and append the resulting variables to the aphid data

library(cdlTools)

# Naming and summing ####

# What landcovers are we dealing with? Let's build a convient key for the codes from the data.
cdl <- read.csv2("https://raw.githubusercontent.com/ReproducibleQM/sapsuckRz/master/Data/cdl_landcover_data_1pt5k.csv", sep=",", check.names = F, 
                 na.strings = "0.0")
landtypes <- colnames(cdl)
landcodes <- updateNamesCDL(landtypes)

codekey <- data.frame(landtypes,landcodes)
View(codekey) #This dataframe lists all of the lables for the various codes we have. All of the codes presented have a "non-zero" presence near 
#a site for at least one year, so all codes must be accounted for.

# Now let's remove all the dumb NA's and replace them with zeros.
for(i in 4:108){
  cdl[,i] <- as.numeric(as.character(cdl[,i]))
}
cdl[is.na(cdl)]<- 0

# Now let's build some categories for the landcover types. We don't need to distinguish between sorghum and barley or between watermelons
# and tomatoes. Some of these categories are simply sums of the smaller categories. Also, note that for now double crops are simply measured twice, 
# once for each category they are in.
for(j in 1:length(cdl)){
  cdl$for_con[j] <- sum(cdl[j,c('142')]) #Conifers: evergreens
  cdl$for_dec[j] <- sum(cdl[j,c('141')]) #Deciduous
  cdl$for_mix[j] <- sum(cdl[j,c('63','143')]) #Mixed: generic forest and mixed forest
  cdl$forest[j] <- cdl$'for_con'[j] + cdl$'for_dec'[j] + cdl$'for_mix'[j] #This is a pooling variable for all forest.
  
  cdl$grass[j] <- sum(cdl[j,c('176','152','60','58')]) #These are non-agricultural grasslands. Includes grass/pasture, shrubland, switchgrass, clover/wildflowers.
  
  cdl$water[j] <- sum(cdl[j,c('111','92')]) #Includes water and aquaculture. Note that large bodies of water are not included (e.g. the Great Lakes).
  
  cdl$urb[j] <- sum(cdl[j,c('121','122','123','124')]) #These are different levels of developed land and roads.
  
  cdl$ag_corn[j] <- sum(cdl[j,c('1', '13', '225', '241', '12', '226', '237')]) #Includes corn, popcorn, sweet corn, and all double crops.
  cdl$ag_beans[j] <- sum(cdl[j,c('5', '26', '241', '240', '254', '239')]) #Note that beans refers to soybeans. Other beans are in the other category.
  #Includes soybeans and all double crops.
  cdl$ag_wheat[j] <- sum(cdl[j,c('24', '23', '205', '39', '26', '225', '22', '236', '238')]) #Includes winter, spring, durum, triticale, buckwheat, and all double crops.
  cdl$ag_smgrains[j] <- cdl$ag_wheat[j] + sum(cdl[j,c('27','28','29','21','236','25','240','254','238','226','235', '237')])
  #Includes wheat plus rye, oats, millet, barley, other small grains, plus double crops.
  cdl$ag_other[j] <- sum(cdl[j,c('4', '61', '36', '59', '44', '6', '37', '32', '42', '41', '68', '43','53', '47', '25',
                                 '74', '240', '254', '48', '67', '2', '209', '239', '46', '238', '70', '229', '216', '57', '10', '76', '222', 
                                 '224', '242', '31', '33', '11', '71', '54', '69', '246', '219', '66', '243', '250', '49', '50', 
                                 '206', '221', '38', '244', '207', '77', '30', '14', '249', '245', '247', '220', '214', 
                                 '55', '248', '35', '210', '34')])
  cdl$ag[j] <- cdl$ag_corn[j] + cdl$ag_beans[j] + cdl$ag_smgrains[j] + cdl$ag_other[j] #This variable pools all arable land.
  
  cdl$wet[j] <- sum(cdl[j,c('190','195','87')]) #These are various wetland categories.
  
  cdl$other[j] <- sum(cdl[j,c('131','0')])
}
#Great! Now we just need to remove the initial numerically titled variables and we're good to go.
cdl_sort <- cdl[,c('ID','State','Year','for_con','for_dec','for_mix','forest','ag_corn','ag_beans','ag_wheat','ag_smgrains','ag_other','ag','wet','urb','water','other')]

rm(cdl, codekey)

# Merging by site and year ####

library(data.table)

# load in aphid data
# Using fread in the data.table package because this file is quite large
aphid<-fread("https://dl.dropboxusercontent.com/u/98197254/Mastersuction_may29_2014.csv")

# Pull the 10 most common species
# find the average number of captures for each species
aphid.ag<-aggregate(Captures~Aphid.Species,data=aphid,sum)
# order data from most number of captures to least
aphid.ord<-aphid.ag[order(-aphid.ag$Captures),] 
# make a list of the 10 most common species
aphid.names.cut<-aphid.ord$Aphid.Species[1:25]
# cut down the data to only include those 10 species
aphid.cut<-aphid[aphid$Aphid.Species%in%aphid.names.cut]

# add a (crude) time variable
aphid.cut$time<-((aphid.cut$Year-2005)*365)+aphid.cut$Date

# clean up
rm(aphid, aphid.ag, aphid.ord, aphid.names.cut)

aphid.cdl <- aphid.cut[aphid.cut$Year %in% unique(cdl_sort$Year)] # remove years w/ no CDL data
cdl_sort$Site <- as.character(cdl_sort$ID)
cdl_sort$ID <- NULL
cdl_sort$State <- as.character(cdl_sort$State)

aphids <- merge(aphid.cdl, cdl_sort, by = c('Year', 'Site'), all = TRUE)
