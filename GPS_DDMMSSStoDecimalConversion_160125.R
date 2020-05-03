## re-formatting DD MM SS to Decimal for GPS data:
## Using example of Anthony Wendt's study sites in Hyderabad, IN:

# Set working directory to your local repo:
setwd("C:/Users/hwells/Desktop/TZMycotoxinSurveyAnalysis")
data <- read.csv("./sitedata_formatdms.csv", stringsAsFactors=FALSE)
Lat <- data$Latitude
Lon <- data$Longitude

## create data frame representing each split element of the Latitude vector (DD, MM, SS.S)
Lat <- as.data.frame(strsplit(Lat, split=" "))
LatDF <- as.data.frame(matrix(NA, nrow=84, ncol=3))
## Now populate the data frame with the correct values:
colnames(LatDF) <- c("DD", "MM", "SS.S")
for (i in 1:nrow(LatDF)){
        LatDF[i,1] <- as.numeric(as.character(Lat[1,i]))
}
for (i in 1:nrow(LatDF)){
        LatDF[i,2] <- as.numeric(as.character(Lat[2,i]))
}
for (i in 1:nrow(LatDF)){
        LatDF[i,3] <- as.numeric(as.character(Lat[3,i]))
}


## create data frame representing each split element of the Longitude vector (DD, MM, SS.S)
Lon <- as.data.frame(strsplit(Lon, split=" "))
LonDF <- as.data.frame(matrix(NA, nrow=84, ncol=3))
colnames(LonDF) <- c("DD", "MM", "SS.S")
for (i in 1:nrow(LonDF)){
        LonDF[i,1] <- as.numeric(as.character(Lon[1,i]))
}
for (i in 1:nrow(LatDF)){
        LonDF[i,2] <- as.numeric(as.character(Lon[2,i]))
}
for (i in 1:nrow(LatDF)){
        LonDF[i,3] <- as.numeric(as.character(Lon[3,i]))
}


## Make new vectors for decimal degrees of Lat and Lon 
Lat_DD <- numeric(length=84)
Lon_DD <- numeric(length=84)

## Populate the DD vectors with the values according to conversion formula from:
## http://andrew.hedges.name/experiments/convert_lat_long/
        ## Decimal Degrees = Degrees + minutes/60 + seconds/3600
Lat_DD <- LatDF[,1] + LatDF[,2]/60 + LatDF[,3]/3600
Lon_DD <- LonDF[,1] + LonDF[,2]/60 + LonDF[,3]/3600

NEWDATA <- as.data.frame(cbind(data, Lat_DD, Lon_DD))
write.csv(NEWDATA, "./new_coords.csv")
