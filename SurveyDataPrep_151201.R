## Kongwa Field Dataset cleanup: Milled maize toxins and attributes
## Henry Wells; 1 December 2015

# Clone repo to desktop and set as working directory:
setwd("C:/Users/hwells/Desktop/TZMycotoxinSurveyAnalysis")

HHW_dat <- read.csv("./Raw Data/HHW_Samples_151201.csv", stringsAsFactors=FALSE)
FMN_dat <- read.csv("./Raw Data/FMN_Samples_151126.csv", stringsAsFactors=FALSE)

library(dplyr)
library(tidyr)
library(reshape2)

library(car)
library(Hmisc)
library(RColorBrewer)
library(ggplot2)

library(sp)
library(rgeos)
library(rgdal)
library(raster)

## Separating out types of samples recorded in dataset:
takeFirst <- function(list){ 
        vec <- character()
        for(i in 1:length(list)) {
                element = list[[i]][1]
                vec = append(vec, element, after=length(vec))
        }
        return(vec)
}

otherTypes = character(length=nrow(HHW_dat))
for (i in 1:length(otherTypes)){
        if(grepl("wholekmaize", HHW_dat$MaizeFlourType[[i]])){
                otherTypes[[i]] = "wholekmaize"
        } else {
                otherTypes[[i]] = "none"
        }
}

## Add parsed values for sample type collected to dataset:
HHW_dat$otherTypes = otherTypes
HHW_dat <- mutate(HHW_dat, "MaizeFlourType" = as.character(
        takeFirst(strsplit(HHW_dat$MaizeFlourType, split=" "))))
for (i in 1:length(HHW_dat$MaizeFlourType)){
        if (HHW_dat$MaizeFlourType[i]=="wholekmaize"){
                HHW_dat$MaizeFlourType[i] = "none"
        } 
}


## Import HHW's locations dataset and integrate with samples data
Loc_dat_HHW <- read.csv("./Raw Data/HHW_GPSLocationData_151201.csv", stringsAsFactors=FALSE)
# View(Loc_dat_HHW)
Loc_dat_HHW$LocationID %in% HHW_dat$locationID
Locations <- Loc_dat_HHW[-1, ]
Mills <- Locations[-16, ]
colnames(Mills) <- c("locationID", "LocationType", "Village", "Ward", 
                     "Latitude_dd", "Longitude_dd", "Elev_masl")
locFreqHW <- as.data.frame(table(HHW_dat$locationID))

freqRep <- function(x, freq){
        newvec <- character()
        for (i in 1:length(x)){
                reps <- as.character(rep(x[i], times=freq[i]))
                newvec <- append(newvec, reps, after=length(newvec))
        }
        return(newvec)
}

lfreq <- locFreqHW[,2]
LID <- freqRep(Mills$locationID, lfreq)
Village <- freqRep(Mills$Village, lfreq)
Ward <- freqRep(Mills$Ward, lfreq)
Latitude_dd <- freqRep(Mills$Latitude_dd, lfreq)
Longitude_dd <- freqRep(Mills$Longitude_dd, lfreq)
Elev_masl <- freqRep(Mills$Elev_masl, lfreq)

## combine location info with sample info into new spreadsheet:
sample_data <- HHW_dat[, c(1:4, 6:7)]
toxin_data <- HHW_dat[, c(11:12)]
location_data <- as.data.frame(cbind(Village, Ward, Latitude_dd, Longitude_dd, Elev_masl))
HHW_cleaned <- as.data.frame(cbind(sample_data, location_data, toxin_data))
# View(HHW_cleaned) ; print(colnames(HHW_cleaned))
write.csv(HHW_cleaned, "./HHW_SamplesWithLocationInfo_151216.csv")


## IMPORTANT NOTE: Merging of data from Aflatoxin and Fumonisin assays on samples of milled maize
## is comprised of two main data sets - samples taken respectively by Henry Wells and Francis Ngure.
## Though merging with location data was possible for Henry's data above, given some one-off issues
## this same operation was ultimately done for Francis' data in excel, outside of this script. 
## The full data set is re-used below and read in under the name "LocationsAndSamples_1DEC_HHW-FMN.csv"



## CENSUS INFORMATION OBTAINED FROM DAICO (admin department for local gov't):

## load in census data:
census <- read.csv(file="./Raw Data/2012 Kongwa Census.csv", stringsAsFactors=FALSE)
census_colnames <- as.character(census[2,])
colnames(census) <- census_colnames
wardsind <- which(grepl("[Ww]ard", census[,2]))
Ward_census <- as.data.frame(census[wardsind,])
Village_census <- as.data.frame(census[-wardsind,])
Village_census <- Village_census[-1,-1]
colnames(Village_census) <- c("Village", "n_HH_2012", "pop_2012", "proj_pop_2014", 
                              "people_per_hh_2014", "proj_pop_growth_rate")
Village_census[28,2:6] <- sub("??", NA, as.character(Village_census[28,2:6]))
Village_census[29,2:6] <- sub("??", NA, as.character(Village_census[29,2:6]))


## LOCATIONS AND SAMPLES COLLECTED: 

## load in locations and samples spreadsheet and obtain means and variances for toxin levels
## in each village...
loc_samples <- read.csv("./Raw Data/LocationsAndSamples_1DEC_HHW-FMN.csv", stringsAsFactors=FALSE)
villages <- as.character(unique(loc_samples$Village))
for(i in 1:length(villages)){
        if(grepl(" ", villages[i])==TRUE){
                split <- strsplit(villages[i], split=" ")
                villages[i] <- as.character(split[[1]][1])
        }
}
villages <- unique(villages)

## Cleaning and subsetting data to create an analytic dataset of ONLY maize meal samples
## remove rows not containing village name or toxin value:
mm_samples <- loc_samples[-which(loc_samples$Village==""), ]
mm_samples <- mm_samples[-which(is.na(mm_samples$Aflatoxin_ppb)), ]
mm_samples[,10:11] <- lapply(mm_samples[,10:11], as.numeric)
## replace all differences in villages' names such that all samples from the 
## same village have the same "Village" value:
for(i in 1:length(mm_samples$Village)){
        if(grepl(" ", mm_samples$Village[i])==TRUE){
                split <- strsplit(mm_samples[i,8], split=" ")
                mm_samples[i,8] <- paste0(split[[1]][1])
        }
}
new_villages <- unique(mm_samples$Village)
new_villages %in% villages ; length(new_villages)

## Next, summarize villages by other types of samples collected
otherSamples <- read.csv("./Raw Data/FMN_samples2_forTables.csv", stringsAsFactors=FALSE)
otherSamples <- as.data.frame(otherSamples[-131,])
# View(otherSamples)
NAtoZero <- function(vec){
        v <- as.character(vec)
        vInd <- which(is.na(vec))
        v[vInd] <- "0"
        v <- as.integer(v)
        return(v)
}
temp <- as.list(otherSamples[,5:15])
for(i in 1:11){
        temp[[i]] <- NAtoZero(temp[[i]])
}
temp <- as.data.frame(temp)
otherSamples <- as.data.frame(cbind(otherSamples[,1:4], temp))
otherSamples <- group_by(otherSamples, Village)
summary_other <- dplyr::summarize(otherSamples, sum(Total.Maize.samples), sum(Millet), sum(Finger.millet),
                                  sum(Sorghum), sum(Rice), sum(Sesame), sum(Groundnuts), 
                                  sum(Home.made.Lishe), sum(Whole.millet.flour), sum(Sifted.Sorghum.flour),
                                  sum(Commercial.Lishe))
WKM_GNUT <- summary_other[,c(1,2,8)]
rowsums <- summary_other[,c(3,4,5,6,7,9,10,11,12)]
rowsumsvec <- integer(length=20)
for(i in 1:20){
        rowsumsvec[i] <- sum(rowsums[i,])
}
WKM_GNUT_other <- as.data.frame(cbind(WKM_GNUT, rowsumsvec))
colnames(WKM_GNUT_other) <- c("Village", "Whole Kernel Maize", "Groundnuts", "Other Samples")
# write.csv(WKM_GNUT_other, file="./otherSamples_summary_21Dec.csv")

## IMPORTANT NOTE: data for toxin content in non-maize and non-groundnut samples (lishe powders, 
## other grains such as millet, sorghum, etc.) were likely to have added variation to assay results
## and were therefore never re-utilized in this analysis.




## extracting means and variances for toxin levels at each village:
grouped_villages <- group_by(mm_samples, Village)
toxin_sum <- dplyr::summarize(grouped_villages, n(), mean(Aflatoxin_ppb), sqrt(var(Aflatoxin_ppb))/sqrt(n()), 
                       mean(Fumonisin_ppm), sqrt(var(Fumonisin_ppm))/sqrt(n()))
# View(toxin_sum)
mean(mm_samples$Aflatoxin_ppb) ; mean(mm_samples$Fumonisin_ppm)
sum(toxin_sum[,2])
sqrt(var(mm_samples$Aflatoxin_ppb))/sqrt(128) ; sqrt(var(mm_samples$Fumonisin_ppm))/sqrt(128)
write.csv(mm_samples, file="./Clean Data/MilledMaize_ToxinSamples_wLocation.csv")
write.csv(toxin_sum, file="./Clean Data/Milled_Maize_toxinsummary_HHW-FMN.csv")




## GPS WAYPOINTS:

## load in GPS data from DNR-Garmin shapefiles (Mills Only):
GPS1 <- shapefile("./GPS_data/waypoints_gps1_TZ.shp")
GPS2 <- shapefile("./GPS_data/waypoints_gps2_TZ.shp")
# View(GPS1) ; View(GPS2)
GPS1coords <- as.data.frame(GPS1[,2:6])
GPS2coords <- as.data.frame(GPS2[,2:6])

GPS_coords <- as.data.frame(rbind(GPS1coords, GPS2coords))
# View(GPS_coords)
full_gps <- as.data.frame(rbind(as.data.frame(GPS1), as.data.frame(GPS2)))
# View(full_gps)
write.csv(GPS_coords, file="./GPS_data/coords_only.csv")
write.csv(full_gps, file="./GPS_data/full_gps_waypoints.csv")

EPSG <- make_EPSG()
wgs84 <- EPSG[grep("4326", EPSG$code), ]    # selecting datum for reference in map projections
coordinates(full_gps) <- ~Longitude+Latitude
Mills <- as.data.frame(full_gps[-c(26,35),]) # renaming from above
coordinates(Mills) <- ~Longitude+Latitude
proj4string(Mills) <- as.character(wgs84$prj4)
KongwaDistHQ <- as.data.frame(full_gps[35,])
coordinates(KongwaDistHQ) <- ~Longitude+Latitude
proj4string(KongwaDistHQ) <- as.character(wgs84$prj4)
KibaigwaMkt <- as.data.frame(full_gps[26,])
coordinates(KibaigwaMkt) <- ~Longitude+Latitude
proj4string(KibaigwaMkt) <- as.character(wgs84$prj4)

# write shapefiles for market and kongwa district HQ locations just in case:
shapefile(KongwaDistHQ, file="./GPS_data/KongwaDistHQ.shp", overwrite=TRUE)
shapefile(KibaigwaMkt, file="./GPS_data/KibaigwaMkt.shp", overwrite=TRUE)

## Now load Francis' HH coordinates, taken through Open Data Kit survey tool:
HHs <- mm_samples[73:128,]
coordinates(HHs) = ~Longitude_dd+Latitude_dd
proj4string(HHs)=CRS("+init=epsg:4326") # set it to lat-long
## and import Kongwa Wards shapefile and TZ regions shapefile
TZwards <- shapefile("./Adm_Shapefiles/wards/TZwards_WGS84.shp")
Kongwa <- TZwards[TZwards$District_N=="Kongwa",]
Tanzania <- shapefile("./Adm_Shapefiles/regions/TZregions_dissolve_WGS84.shp")
Dodoma <- Tanzania[Tanzania$Region_Nam=="Dodoma", ]


## to add geographic association to a specific mill for the samples that Francis Ngure took 
## in the field, we needed to do a simple minimum-distance calculation between the waypoint
## at the mill and the coordinates stored in the survey tool for where a specific sample was taken
## (while Henry sampled with people directly at mills, Francis walked around the surrounding area
## and encouraged other residents to participate in the study - meaning that his samples were
## not necessarily collected in the same place...hence the necessity of this step.)


# Find the minimum distance between each HH and the nearest posho mill by using coordinates:
distMAT <- pointDistance(HHs, Mills, lonlat=TRUE)
# View(distMAT)
## Need to rename columns of distMAT with the LocationIDs of the posho mills:
MillNames <- tolower(as.character(Mills$ident))
MillNames <- gsub("-", "_", MillNames)
colnames(distMAT) <- MillNames
FMN_sampleIDs <- as.character(mm_samples[73:128,3])
row.names(distMAT) <- FMN_sampleIDs

MinDistances <- function(distMAT){
        WhichLoc <- character()
        Distances <- numeric()
        for(i in 1:nrow(distMAT)){
                indMIN <- as.integer(which.min(distMAT[i,]))
                name <- as.character(names(distMAT[i,])[indMIN])
                distance <- as.numeric(distMAT[i,indMIN])
                WhichLoc <- append(WhichLoc, name, after=length(WhichLoc))
                Distances <- append(Distances, distance, after=length(Distances))
        }
        IDs <- as.character(row.names(distMAT))
        DF <- as.data.frame(cbind(IDs, Distances, WhichLoc))
        DF[,c(1,3)] <- lapply(DF[,c(1,3)], as.character)
        DF[,2] <- as.numeric(DF[,2])
        return(DF)
}
Nearest_Mills <- MinDistances(distMAT)
FMN_new <- mm_samples[73:128, ]
# cbind(Nearest_Mills$IDs, mm_samples[73:128, 3]) 
        ## ^^ Check to see that the sampleID vectors are the same.
FMN_new$locationID <- as.character(Nearest_Mills[,3])
HHW_mm_samples <- mm_samples[1:72, ]
mm_samples_atMill <- rbind(HHW_mm_samples, FMN_new) 
colnames(mm_samples_atMill) <- c("X", "locationID", "SampleID", "MaizeFlourType", "otherTypes",
                                        "Source", "MonthsSinceHarvested", "Village", "Ward", 
                                        "Latitude_dd", "Longitude_dd", "Elev_masl", 
                                        "Aflatoxin_ppb", "Fumonisin_ppm")
write.csv(mm_samples_atMill, file="./Clean Data/mm_samples_atMill.csv")

## Now group mm_samples_atMill by locationID and summarize for average and std. error values for AF:
grouped_mills <- group_by(mm_samples_atMill, locationID)
## Note here that the wrapper 'dplyr::fxn' is used to avoid conflict created when using n() directly.
        ## BUT: conflicts(), at this time, did not show n() listed among them! So what's happening??
        ## Investigate this later on...
MillSummary <- dplyr::summarize(grouped_mills, n(), mean(Aflatoxin_ppb), (sd(Aflatoxin_ppb)/sqrt(n())), 
                         mean(Fumonisin_ppm), (sd(Fumonisin_ppm)/sqrt(n())),
                         (sum(Source=="homegrown"))/n())
colnames(MillSummary) <- c("locationID", "N_samples", "Mean_AF_ppb", "SE_AF_ppb", 
                           "Mean_FUM_ppm", "SE_FUM_ppm", "p_homegrown")
## Now add these summarized values to the mills' coordinates SPDF as attributes:
## First, re-order rows in 'Mills' SPDF according to chronological order of LocationIDs:
Mills$ident <- tolower(Mills$ident)
Mills$ident <- gsub("-", "_", Mills$ident)
newRowIndices <- match(Mills$ident, MillSummary$locationID)     # so: move rows 17-25 of 'Mills' SPDF
                                                                # to the bottom of that same DF...
MillsTEST <- Mills[-(17:25),]
MillsTEST2 <- Mills[17:25,]
Mills_Summarized <- rbind(MillsTEST, MillsTEST2)
# match(Mills_Summarized$ident, MillSummary$locationID)   ## Now shows 1:33 in ascending order - done!
## Now join the values with cbind():
Mills_Summarized <- as.data.frame(cbind(Mills_Summarized, MillSummary))
coordinates(Mills_Summarized) <- ~Longitude+Latitude    ## Returns Mills_Summarized to SPDF s4 class.
proj4string(Mills_Summarized) <- "+proj=longlat +datum=WGS84 +no_defs"

## test plot:
plot(Kongwa)
plot(Mills_Summarized, add=TRUE)                      ## Everything worked!
## Now write new shapefile for Mills that includes summarized toxin data, and for adm boundaries:
# shapefile(Mills_Summarized, file="./GPS_data/MillsWithData.shp", overwrite=TRUE)
# shapefile(Kongwa, file="./Adm_Shapefiles/KongwaWards.shp", overwrite=TRUE)
# shapefile(Dodoma, file="./Adm_Shapefiles/DodomaReg.shp", overwrite=TRUE)





## CUSTOMERS OF POSHO MILLS:

## load in Mill Customers data and break down ag threats and signs of AF by village:
Mill_customers <- read.csv("./Raw Data/MillCustomers_final_16Dec_HHW.csv", stringsAsFactors=FALSE)
for_table2 <- Mill_customers[, c(1,7,9,11,12,13,14,15)] 
sampleID <- as.character(for_table2[,1])
LocationID <- character(length=length(sampleID))
for(i in 1:length(LocationID)){
        LocationID[i] <- substr(sampleID[i], 1, 10)
}
# test <- as.data.frame(cbind(LocationID, sampleID))    ## ID's check out against one another.
for_table2 <- as.data.frame(cbind(for_table2[,1], LocationID, for_table2[,2:8]))

## Need to extract locationIDs and Village names from other data to match with customers data...
HHW_mm_samples <- mm_samples[1:72, ]
LIDs_mm <- as.character(HHW_mm_samples$locationID)
Villages_mm <- as.character(HHW_mm_samples$Village)
test <- as.data.frame(cbind(LIDs_mm, Villages_mm))
UniqueLIDs <- unique(LIDs_mm)
village <- character(length=length(UniqueLIDs))
test2 <- as.data.frame(cbind(UniqueLIDs, village), stringsAsFactors=FALSE)
for(i in 1:nrow(test2)){
        rowNum <- as.integer(match(test2[i,1], as.character(test[,1])))
        test2[i,2] <- paste0(as.character(test[rowNum,2]))
}
# View(test2)
## Repeat names of villages according to frequency of location IDs in a new vector:
Villages_for_table2 <- character()
for(i in 1:nrow(for_table2)){
        rowNum <- as.integer(match(for_table2[i,2], test2[,1]))
        Villages_for_table2 <- append(Villages_for_table2, test2[rowNum,2], 
                                      after=length(Villages_for_table2))
}

## Table 2 constructed from updated customers .csv data file in visualizaton script...
customers_data <- as.data.frame(cbind(for_table2[,1:2], Villages_for_table2, for_table2[,3:9]))
colnames(customers_data) <- c("SampleID", "LocationID", "Village", "harvestMonth", "storesDuration",
                              "otherFoods", "agThreats", "AFatHarvest", "AFinStorage", "AFpractices")
write.csv(customers_data, file="./MillCustomers_151216.csv")
