## Integrating GIS data sources into TZ survey data summaries
## Henry Wells
## 21 December 2015; last update 31 March 2016

# Clone repo to desktop and set as working directory:
setwd("C:/Users/hwells/Desktop/TXMycotoxinSurveyAnalysis")
library(raster)
library(rgeos)
library(rgdal)
EPSG <- make_EPSG()

# Load GPS waypoints and Wards and District Shapefiles:
Kongwa <- shapefile("./Adm_Shapefiles/KongwaWards.shp")
Dodoma <- shapefile("./Adm_Shapefiles/DodomaReg.shp")
Tanzania <- shapefile("./Adm_Shapefiles/regions/TZregions_dissolve_WGS84.shp")

Mills <- shapefile("./GPS_data/MillsWithData.shp")



## OTHER GIS DATA LAYERS: RAINFALL, TEMPERATURE, ELEVATION, POPULATION DEMOGRAPHICS
## Elevation:
TZelev <- raster("./Map_Layers/CGIAR-SRTM/TZ_alt/TZA_alt.grd")
proj4string(TZelev) <- "+proj=longlat +datum=WGS84 +no_defs"
Kongwa_ext <- extent(Kongwa)
Kongwa_elev <- crop(TZelev, Kongwa_ext)
writeRaster(Kongwa_elev, filename="./Clean Map Layers/Kongwa_elev_LongLat.tiff", format="GTiff",
            overwrite=TRUE)
## reprojecting altitude raster and shapefiles:
Kongwa_UTM <- Kongwa
Kongwa_UTM <- spTransform(Kongwa_UTM, CRS=CRS("+init=epsg:21035"))
bbox <- as.data.frame(bbox(Kongwa_UTM))
KongwaRaster_UTM <- raster(xmn=bbox[1,1], xmx=bbox[1,2], ymn=bbox[2,1], ymx=bbox[2,2], nrow=110, 
                           ncol=86, crs=CRS("+init=epsg:21035"))
Kelev_UTM <- projectRaster(Kongwa_elev, KongwaRaster_UTM, res=res(Kelev_UTM), 
                           crs=CRS("+init:epsg=21035"))
writeRaster(Kelev_UTM, filename="./Clean Map Layers/Kelev_UTM.tif", format="GTiff", 
            overwrite=TRUE)
Mills_UTM <- Mills
Mills_UTM <- spTransform(Mills_UTM, CRS=CRS("+init=epsg:21035"))

## Writing standard error function:
stderror <- function(x) { (sd(x))/sqrt(length(x)) }

## Extract elevation values for 5km buffer around mill points (based on CGIAR-SRTM DEM for TZ): 
Mills_buffer <- gBuffer(Mills_UTM, byid=TRUE, id=as.character(Mills_UTM$ident), width=5000)
shapefile(Mills_buffer, filename="./Clean Map Layers/Mills_5kbuffer.shp", overwrite=TRUE)
Mills_8kbuffer <- gBuffer(Mills_UTM, byid=TRUE, id=as.character(Mills_UTM$ident), width=8000)
shapefile(Mills_8kbuffer, filename="./Clean Map Layers/Mills_8kbuffer.shp", overwrite=TRUE)

elev_buffer_means <- extract(Kelev_UTM, Mills_8kbuffer, method='simple', fun=mean)
elev_buffer_n_list <- extract(Kelev_UTM, Mills_8kbuffer, method='simple', cellnumbers=TRUE)
elev_buffer_n <- sapply(elev_buffer_n_list, length)
elev_buffer_sd <- extract(Kelev_UTM, Mills_8kbuffer, method='simple', fun=sd)
elev_ward_means <- extract(Kelev_UTM, Kongwa_UTM, method='simple', fun=mean)
elev_ward_n_list <- extract(Kelev_UTM, Kongwa_UTM, method='simple', cellnumbers=TRUE)
elev_ward_n <- sapply(elev_ward_n_list, length)
elev_ward_sd <- extract(Kelev_UTM, Kongwa_UTM, method='simple', fun=sd)

Ward_Table2 <- as.data.frame(cbind(elev_ward_means, elev_ward_n, elev_ward_sd))
MillBufferdata <- as.data.frame(cbind(elev_buffer_means, elev_buffer_n, elev_buffer_sd))

## Given differing number of cells in each ward, need to check to see if extract() re-orders wards...
        ## Use cell numbers as index and reassign values for each cell in each different ward to be the 
        ## number associated with that ward (using object 'elev_ward_n_list')
# cellnum <- integer()
# wardnum <- integer()
# for (i in 1:22){
#         cellnum <- append(cellnum, elev_ward_n_list[[i]][,1], after=length(cellnum))
#         wardnum <- append(wardnum, rep(i, times=length(elev_ward_n_list[[i]])/2, after=length(cellnum)))
# }
# r <- raster(nrows=110, ncols=86, xmn=1527920, xmx=1608398, 
#             ymn=9278668, ymx=9382185, crs=proj4string(Kelev_UTM))
# temp <- as.data.frame(cbind(cellnum, wardnum))
# for (i in 1:22){
#         r[temp[which(temp[,2]==i),1]] <- i
# }




## Historical Rainfall
WorldClim12_r <- raster("./Map_Layers/WorldClim/precipitation/prec_37_tif/prec12_37.tif")
proj4string(WorldClim12_r) <- "+proj=longlat +datum=WGS84 +no_defs"
Kongwa_12_r <- crop(WorldClim12_r, Kongwa_ext)
WorldClim1_r <- raster("./Map_Layers/WorldClim/precipitation/prec_37_tif/prec1_37.tif")
proj4string(WorldClim1_r) <- "+proj=longlat +datum=WGS84 +no_defs"
Kongwa_1_r <- crop(WorldClim1_r, Kongwa_ext)
WorldClim2_r <- raster("./Map_Layers/WorldClim/precipitation/prec_37_tif/prec2_37.tif")
proj4string(WorldClim2_r) <- "+proj=longlat +datum=WGS84 +no_defs"
Kongwa_2_r <- crop(WorldClim2_r, Kongwa_ext)
WorldClim3_r <- raster("./Map_Layers/WorldClim/precipitation/prec_37_tif/prec3_37.tif")
proj4string(WorldClim3_r) <- "+proj=longlat +datum=WGS84 +no_defs"
Kongwa_3_r <- crop(WorldClim12_r, Kongwa_ext)
WorldClim4_r <- raster("./Map_Layers/WorldClim/precipitation/prec_37_tif/prec4_37.tif")
proj4string(WorldClim4_r) <- "+proj=longlat +datum=WGS84 +no_defs"
Kongwa_4_r <- crop(WorldClim4_r, Kongwa_ext)
WorldClim5_r <- raster("./Map_Layers/WorldClim/precipitation/prec_37_tif/prec5_37.tif")
proj4string(WorldClim5_r) <- "+proj=longlat +datum=WGS84 +no_defs"
Kongwa_5_r <- crop(WorldClim5_r, Kongwa_ext)
WorldClim6_r <- raster("./Map_Layers/WorldClim/precipitation/prec_37_tif/prec6_37.tif")
proj4string(WorldClim6_r) <- "+proj=longlat +datum=WGS84 +no_defs"
Kongwa_6_r <- crop(WorldClim6_r, Kongwa_ext)
WorldClim7_r <- raster("./Map_Layers/WorldClim/precipitation/prec_37_tif/prec7_37.tif")
proj4string(WorldClim7_r) <- "+proj=longlat +datum=WGS84 +no_defs"
Kongwa_7_r <- crop(WorldClim7_r, Kongwa_ext)

# Summarizing rainfall...
WC_PreHarv_rainfall <- Kongwa_12_r + Kongwa_1_r + Kongwa_2_r + Kongwa_3_r + Kongwa_4_r 
WC_PeriHarv_rainfall <- Kongwa_5_r + Kongwa_6_r + Kongwa_7_r
## get means according to same buffer used for elevation:
WC_PreH_rain_UTM <- projectRaster(WC_PreHarv_rainfall, KongwaRaster_UTM, res=res(Kelev_UTM), 
                                     crs=CRS("+init:epsg=21035"))
WC_PeriH_rain_UTM <- projectRaster(WC_PeriHarv_rainfall, KongwaRaster_UTM, res=res(Kelev_UTM), 
                                      crs=CRS("+init:epsg=21035"))
writeRaster(WC_PreH_rain_UTM, filename="./Clean Map Layers/WC_PreH_RFE_clark1880.tif", 
            format="GTiff", overwrite=TRUE)
writeRaster(WC_PeriH_rain_UTM, filename="./Clean Map Layers/WC_PeriH_RFE_clark1880.tif", 
            format="GTiff", overwrite=TRUE)
writeRaster(WC_PreHarv_rainfall, filename="./Clean Map Layers/WC_PreH_RFE.tif", 
            format="GTiff", overwrite=TRUE)
writeRaster(WC_PeriHarv_rainfall, filename="./Clean Map Layers/WC_PeriH_RFE.tif", 
            format="GTiff", overwrite=TRUE)

WC_PreH_rain_means <- extract(WC_PreH_rain_UTM, Mills_8kbuffer, method='simple', fun=mean)
WC_PeriH_rain_means <- extract(WC_PeriH_rain_UTM, Mills_8kbuffer, method='simple', fun=mean)

WC_preHrain_wardMeans <- extract(WC_PreH_rain_UTM, Kongwa_UTM, method='simple', fun=mean)
WC_periHrain_wardMeans <- extract(WC_PeriH_rain_UTM, Kongwa_UTM, method='simple', fun=mean)

# Use terrain() function to calculate aspect from elevation layer:
aspect.deg <- terrain(Kongwa_elev, opt='aspect', unit='degrees')
writeRaster(aspect.deg, filename="./Clean Map Layers/aspect_deg_UTM.tif", format="GTiff",
            overwrite=TRUE)
slope <- terrain(Kongwa_elev, opt='slope', unit='degrees')
writeRaster(slope, filename="./Clean Map Layers/slope_UTM.tif", format="GTiff",
            overwrite=TRUE)

# reclassify aspect into categories:
rcl <- matrix(c(
        0,  22.5, 1,
        22.5, 77.5, 2,
        77.5, 112.5, 3,
        112.5, 157.5, 4,
        157.5, 202.5, 5,
        202.5, 247.5, 6,
        247.5, 292.5, 7,
        292.5, 337.5, 8,
        337.5, 360, 1), ncol=3, byrow=TRUE)
aspect <- reclassify(aspect.deg, rcl=rcl)
writeRaster(aspect, filename="./Clean Map Layers/aspect_categorical_UTM.tif", format="GTiff",
            overwrite=TRUE)

# Note that: 1=N, 2=NE, 3=E, 4=SE, 5=S, 6=SW, 7=W, 8=NW
## reproject slope raster, as well as both degree and categorical aspect rasters
slope.utm <- projectRaster(slope, KongwaRaster_UTM, res=res(Kelev_UTM), 
                           crs=CRS("+init:epsg=21035"))
aspect.deg.utm <- projectRaster(aspect.deg, KongwaRaster_UTM, res=res(Kelev_UTM), 
                           crs=CRS("+init:epsg=21035"))
aspect.utm <- projectRaster(aspect, KongwaRaster_UTM, res=res(Kelev_UTM), 
                            crs=CRS("+init:epsg=21035"))
# Find out the dominant aspect category within the area around each mill point:
maxFreq <- function(vec){
        uniques <- as.data.frame(table(vec))
        max <- max(uniques$Freq)
        ind <- which(uniques$Freq==max)
        return(as.character(uniques[ind,1]))
}
dom.aspect <- extract(aspect.utm, Mills_8kbuffer, method='simple', 
                         fun=function(x,...)as.character(maxFreq(x)))
# buffer.aspect.mean <- extract(aspect.utm, Mills_8kbuffer, method='simple', fun=mean)
# View(as.data.frame(cbind(buffer.aspect.mean, buffer.aspect)))
        # Goes to show that mean() and maxFreq() do not lead to the same result!
avg.slope <- extract(slope.utm, Mills_8kbuffer, method='simple', fun=mean)


# joining data about wards and mills into new Data Frames: (EDIT THIS SECTION AS NECESSARY)
# Mills:
milldat <- as.data.frame(Mills, stringsAsFactors=FALSE)
LocID <- as.character(milldat$ident)
MILLnum <- c(1:33)
toxindat <- milldat[,46:51]
MillDF <- as.data.frame(cbind(LocID, MILLnum, elev_buffer_n, elev_buffer_means, elev_buffer_sd,
                              WC_PreH_rain_means, WC_PeriH_rain_means, dom.aspect, avg.slope,
                              toxindat))
colnames(MillDF) <- c("LocID", "MILLnum", 
                            "n_cell", "elev_mean", "elev_sd", "PreH.RFE", "PeriH.RFE", 
                                "dom.asp", "avg.slope", "n_samples",
                                "AF_mean", "AF_stderr", "FUM_mean", "FUM_stderr", "p_hg")

MillGeoDF <- as.data.frame(cbind(Mills_UTM[,c(1:4,42:45)], MillDF[,c(2,4:12,15)]))
MillGeoDF <- SpatialPointsDataFrame(coords=as.matrix(MillGeoDF[,c(4,3)]), data=MillGeoDF)
proj4string(MillGeoDF) <- "+proj=longlat +datum=WGS84 +no_defs"
CRS <- crs(proj4string(Mills_UTM))
MillGeoDF_Clark1880 <- spTransform(MillGeoDF, CRS=CRS)

shapefile(MillGeoDF_Clark1880, filename="./Clean Map Layers/MillGeo_Clark1880.shp", 
          overwrite=TRUE)
shapefile(MillGeoDF, filename="./Clean Map Layers/MillGeo_LongLat.shp", overwrite=TRUE)
shapefile(Kongwa_UTM, filename="./Clean Map Layers/Kongwa_Clark1880.shp", overwrite=TRUE)
write.csv(MillGeoDF, file="./MillGeoDF.csv")
warddat <- as.data.frame(Kongwa_UTM)
Ward_Name <- as.character(warddat$Ward_Name)
WARDnum <- c(1:length(Ward_Name))
# preliminary linear model of AF_mean ~ elev_mean + PreH_rain_mean +
#                                                       PeriH_rain_mean + p_hg

## As a preliminary visual test for the model, make a corrgram of all variables:
# corrgram(MillDF[,c(4:7, 10:13)], upper.panel=panel.conf, diag.panel=panel.density, 
#          lower.panel=panel.pts)

mod <- lm(AF_mean ~ elev_mean+PreH.RFE+PeriH.RFE+dom.asp+avg.slope+p_hg, 
          data=MillDF)
# summary(mod) WOW! Rsq. of > 0.5. 


## Use the DF 'mm_samples_atMill' from surveydata_cleaning_1Dec15 script to re-run that regression
## on each sample point itself (turning 'source' into a dummy variable)
mm_samples_atMill <- read.csv("./mm_samples_atMill.csv", stringsAsFactors=FALSE)
colnames(mm_samples_atMill) <- c("X.1","X","locationID","SampleID","MaizeFlourType",
                                 "otherTypes","Source","MonthsSinceHarvested","Village","Ward",
                                 "Latitude_dd","Longitude_dd","Elev_masl","Aflatoxin_ppb",
                                 "Fumonisin_ppm")
samples <- dplyr::mutate(mm_samples_atMill, "is.homegrown"=as.logical(Source=="homegrown"))
# colnames(samples)
# [1] "X.1"                  "X"                    "locationID"           "SampleID"            
# [5] "MaizeFlourType"       "otherTypes"           "Source"               "MonthsSinceHarvested"
# [9] "Village"              "Ward"                 "Latitude_dd"          "Longitude_dd"        
# [13] "Elev_masl"            "Aflatoxin_ppb"        "Fumonisin_ppm"        "is.homegrown"

# repeat geographic variable measurements by indices of locationID in samples DF

geovars <- as.data.frame(matrix(data=0.0, nrow=128, ncol=16))
DFtest <- as.data.frame(MillGeoDF[, c(1:7,11:17)])
for (i in 1:nrow(DFtest)){
        ind <- which(as.character(samples$locationID) %in% as.character(DFtest[i,2]))
        geovars[ind, ] = DFtest[i, ]
}

samples <- as.data.frame(cbind(samples, geovars))
colnames(samples) <- as.character(c(colnames(samples[,1:16]),colnames(DFtest)))
samples$dom.asp <- as.factor(samples$dom.asp)
SampleGeoDF <- dplyr::mutate(samples, "LnAF"=log(Aflatoxin_ppb))
SampleGeoDF <- SampleGeoDF[,-c(1,2)]

SampleGeoDF <- SpatialPointsDataFrame(coords=as.matrix(SampleGeoDF[,c(18,17)]), data=SampleGeoDF)
proj4string(SampleGeoDF) <- "+proj=longlat +datum=WGS84 +no_defs"
CRS <- crs(proj4string(Mills_UTM))
SampleGeoDF_Clark1880 <- spTransform(SampleGeoDF, CRS=CRS)

shapefile(SampleGeoDF, filename="./Clean Map Layers/SampleGeo_LongLat.shp", overwrite=TRUE)
shapefile(SampleGeoDF_Clark1880, filename="./Clean Map Layers/SampleGeo_Clark1880.shp", 
          overwrite=TRUE)

distr_AF <- ggplot(samples, aes(x=Aflatoxin_ppb, fill=Source)) + geom_density(alpha=.3) + 
        labs(y = "Probability Density") + 
        ggtitle("Probability Density of Aflatoxin levels by source of grain")

################## Need to fix code beyond this point.


# Now, try out a model with 'samples':
mod_samples <- lm(Aflatoxin_ppb~elev_mean+PreH_rain_mean+PeriH_rain_mean+is.homegrown+
                          dom.aspect+avg.slope, data=samples)
# summary(mod_samples) # Model significant at alpha=.001, but only one significant predictor 
                        # PreH_rain_mean, at alpha=.1

## clear skewness in the density of Aflatoxin data for the samples...try a log model:
logmod_samples <- lm(log(Aflatoxin_ppb)~ elev_mean+PreH_rain_mean+PeriH_rain_mean+
                             is.homegrown+dom.aspect+avg.slope, data=samples)
# summary(logmod_samples) 


# Stepwise Regression using stepAIC() from MASS package
# on 'mod' from above (using Mill-level averages...)
library(MASS)
step1aic <- stepAIC(mod, direction="both")
# summary(step1aic)     Best model: AF ~ PreH_rain + PeriH_rain + dom.aspect
step1bic <- stepAIC(mod, direction="both", k=log(33))
# summary(step1bic)         much less interesting result than with AIC...given what we
                        #       know about the spread of A. flavus in the field

# Try stepwise on mod_samples
step2aic <- step(mod_samples, direction="both")
# summary(step2aic)
step2bic <- step(mod_samples, direction="both", k=log(128))
# summary(step2bic)

step3aic <- step(logmod_samples, direction="both")
# summary(step3aic)     # clearly should use log model to eliminate skewness in AF data
                        # interesting result! Still decently good fit (in terms of R^2)
step3bic <- step(logmod_samples, direction="both", k=log(128))
# summary(step3bic)       # Ideal model selected differs when using BIC instead!
                        # uses avg.slope instead of dom.aspect...to a lower R^2. Huh. Why?...

## Spatial Autocorrelation in step and step3 models...




WardDF <- as.data.frame(cbind(Ward_Name, WARDnum, elev_ward_n, elev_ward_means, elev_ward_sd,
                              WC_preHrain_wardMeans, WC_periHrain_wardMeans))
colnames(WardDF) <- c("Ward_Name", "WARDnum", "n_cell", 
                      "elev_mean", "elev_sd", "PreH_rain_mean", "PeriH_rain_mean")
write.csv(WardDF, file="./WardGeoDF.csv")

# Merge WardDF with other qualitative obs from survey questionnaires in excel - saved as 
# .csv called 'table2_tripreport.csv'

