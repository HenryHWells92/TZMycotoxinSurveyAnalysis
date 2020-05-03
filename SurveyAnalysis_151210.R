## Henry H. Wells
## 10 December 2015; last update 31 March 2016
## Developing data summaries and visualizations
        ## November 2015 Field Trip to Kongwa, Dodoma, Tanzania

# clone this repo to your desktop and set as working directory:
setwd("C:/Users/hwells/Desktop/TZMycotoxinSurveyAnalysis")

## IMPORTANT NOTE: The following data cleaning operations are assumed to have already been performed
## when running this script. All data created in the two scripts below are stored in the files named
## "Clean Map Layers" and "Clean Data", respectively. 

# source('C:/Users/hwells/Desktop/TZMycotoxinSurveyAnalysis/SurveyDataPrep_151201.R', echo=FALSE)
# source('C:/Users/hwells/Desktop/TZMycotoxinSurveyAnalysis/GISDataPrep_151201.R', echo=FALSE)






## Load survey response data from Posho Mill customers:
customers <- read.csv("./Clean Data/MillCustomers_151216.csv", stringsAsFactors=FALSE)

## Graph 1: Dotplot of AF in each village with standard error bars...
mm_samples <- read.csv("./Clean Data/MilledMaize_ToxinSamples_wLocation.csv", stringsAsFactors=FALSE)
grouped_wards <- group_by(mm_samples, Ward)
wards_toxins <- as.data.frame(dplyr::summarize(grouped_wards, n(), mean(Aflatoxin_ppb), 
                          sqrt(var(Aflatoxin_ppb))/sqrt(n()), mean(Fumonisin_ppm), 
                          sqrt(var(Fumonisin_ppm))/sqrt(n())))
colnames(wards_toxins) <- c("Ward", "n_samples", "Mean_AF", "Std_Error_AF", "Mean_FUM", "Std_Error_FUM")

PlotThis <- mm_samples[, c(3,4,7,9,10,14,15)]
Ward_mean_AF <- numeric(length=as.integer(nrow(PlotThis)))
Ward_SE_AF <- numeric(length=as.integer(nrow(PlotThis)))
PlotThis <- as.data.frame(cbind(PlotThis, Ward_mean_AF, Ward_SE_AF))
for(i in 1:nrow(PlotThis)){
        ward <- as.character(PlotThis[i,5])
        rowNum <- as.integer(match(ward, as.character(wards_toxins[,1])))
        PlotThis[i,8] <- as.numeric(wards_toxins[rowNum,3])
        PlotThis[i,9] <- as.numeric(wards_toxins[rowNum,4])
}
# View(PlotThis)
colnames(PlotThis) <- c("locationID", "SampleID", "Source", "Village", "Ward", 
                        "Aflatoxin_ppb", "Fumonisin_ppm", "Ward_mean_AF", "Ward_SE_AF")

library(ggplot2)
library(RColorBrewer)
library(Hmisc) ## Instead of std error functions below, use "mean_cl_normal"


## PLOT 1: Dotplot of Aflatoxin by ward, with standard error bars (+/- 1 SE of Mean)
PlotThis$Ward <- factor(PlotThis$Ward, 
                        levels=unique(PlotThis$Ward[order(PlotThis$Ward_mean_AF, 
                                                          decreasing=TRUE)]))
G1 <- ggplot(PlotThis, aes(x=Ward, y=Aflatoxin_ppb)) + scale_color_brewer(palette="Set1") + 
        geom_point(aes(color=PlotThis[,4]), size=3) + theme_bw(base_family = "Calibri") + 
        theme(plot.title = element_text(size = rel(1.5)), 
                axis.title.x = element_text(size = rel(1.35)),
                axis.title.y = element_text(size = rel(1.35)),
                axis.text.x = element_text(size = rel(1.2), angle=45, hjust=1),
                axis.text.y = element_text(size = rel(1.2)),
                legend.text = element_text(size = rel(1.2)),
                legend.title = element_text(size = rel(1.25))) +
        stat_summary(fun.y="mean", geom="point", color="black", size=4) +
        stat_summary(fun.data = mean_cl_normal, geom = "errorbar", mult = 1, color="black") + 
        labs(color='Source') + labs(x="Ward", y="Aflatoxin Level (ppb)") +
        labs(title="Aflatoxin in Samples of Milled Maize \n Error Bars: +/- 1 Std. Error of Mean")
## export and save plot as "AF_dotplot_byWard_HHW_16Dec15"



## PLOT 2: Barplot of Aflatoxin in Different ranges (0-10, 10-20, 20-50, 50+), 
        ## for homegrown vs. purchased grain
lessThanTen <- function(vec) {
        ind <- which(vec<10)
        perc <- (length(ind)/length(vec))*100
        return(as.numeric(perc))
}
TenToTwenty <- function(vec) {
        ind <- which(10<=vec & vec<20)
        perc <- (length(ind)/length(vec))*100
        return(as.numeric(perc))
}
TwentyToFifty <- function(vec) {
        ind <- which(20<=vec & vec<50)
        perc <- (length(ind)/length(vec))*100
        return(as.numeric(perc))
}
MoreThanFifty <- function(vec) {
        ind <- which(50<=vec)
        perc <- (length(ind)/length(vec))*100
        return(as.numeric(perc))
}
P_Aflatoxin_ppb <- as.integer(which(PlotThis[,3]=="purchased"))
P_Aflatoxin_ppb <- PlotThis[P_Aflatoxin_ppb, c(2,4,5,6)]

H_Aflatoxin_ppb <- as.integer(which(PlotThis[,3]=="homegrown"))
H_Aflatoxin_ppb <- PlotThis[H_Aflatoxin_ppb, c(2,4,5,6)]

PR_perc <- as.numeric(c(lessThanTen(P_Aflatoxin_ppb[,4]), TenToTwenty(P_Aflatoxin_ppb[,4]),
                        TwentyToFifty(P_Aflatoxin_ppb[,4]), MoreThanFifty(P_Aflatoxin_ppb[,4])))
HG_perc <- as.numeric(c(lessThanTen(H_Aflatoxin_ppb[,4]), TenToTwenty(H_Aflatoxin_ppb[,4]),
                        TwentyToFifty(H_Aflatoxin_ppb[,4]), MoreThanFifty(H_Aflatoxin_ppb[,4])))

breakdown <- as.data.frame(rbind(HG_perc, PR_perc))
types <- c("Homegrown", "Purchased")
breakdown <- as.data.frame(cbind(types, breakdown))
colnames(breakdown) <- c("Type", "less than 10 ppb", "10 to 20 ppb", "20 to 50 ppb", "more than 50 ppb")
PlotThis2 <- melt(breakdown, id.vars=1, measure.vars=2:5)
# View(PlotThis2)

## Find group sample sizes and means:
Homegrown <- as.numeric(c(nrow(H_Aflatoxin_ppb), mean(H_Aflatoxin_ppb[,4])))
Purchased <- as.numeric(c(nrow(P_Aflatoxin_ppb), mean(P_Aflatoxin_ppb[,4])))
summary_AF_Source <- as.data.frame(rbind(Homegrown, Purchased))
colnames(summary_AF_Source) <- c("N", "Mean AF (ppb)") ; # View(summary_AF_Source)
# write.csv(summary_AF_Source, file = './summarized_AF-in-mm-samples_bySource.csv') # might be useful...

# plotting code:
PlotThis2$Category <- factor(PlotThis2$Category, 
                             levels=unique(PlotThis2$Category[order(PlotThis2$X)]))
G2 <- ggplot(PlotThis2, aes(x=Type, y=Percentage, fill=Category)) + 
        theme_grey(base_family = "Calibri") +
        geom_bar(aes(fill=Category), stat='identity', color="black") +
        scale_fill_brewer(palette="Reds")
G2 <- G2 + theme(plot.title = element_text(size = rel(1.5)), axis.title.x = element_text(size = rel(1.35)),
        axis.title.y = element_text(size = rel(1.35)), axis.text.x = element_text(size = rel(1.25)),
        axis.text.y = element_text(size = rel(1.25)), legend.text = element_text(size = rel(1.25)),
                                               legend.title = element_text(size = rel(1.25))) +
        labs(title="Percentage breakdown of Aflatoxin levels in \n Homegrown vs. Purchased Maize") + 
        labs(x="Source of Grain", y="Percentage of all samples in each category") +
        annotate("text", x = "Homegrown", y = 20, label = "N = 76", size = rel(4)) +
        annotate("text", x = "Purchased", y = 20, label = "N = 52", size = rel(4)) +
        annotate("text", x = "Homegrown", y = 15, label = "Mean AF = 20.2 ppb", size = rel(4)) +
        annotate("text", x = "Purchased", y = 15, label = "Mean AF = 35.6 ppb", size = rel(4))
## export and save plot as "AF_percBreakdown_bySource_HHW_16Dec15"

## Two-sample T-Test: Mean AF by source
# subsetting and reshaping data:
TestThis <- as.data.frame(PlotThis[,c(2,3,6)])
# View(TestThis)
library(car)
# Levene test to investigate homogeneity of variance assumption:
t.lm <- lm(Aflatoxin_ppb~Source, data=TestThis)
leveneTest(t.lm)        ## P-value: 0.027 (significant) -- must assume unequal variances in T-test:
# T Test: difference in means of AF_ppb by source:
t <- t.test(Aflatoxin_ppb ~ Source, paired = FALSE, 
            var.equal = FALSE, data = TestThis)
t       ## Test stat: t = -2.419, p-val: 0.0176 *
        ## Statistically significant difference between 'Source' groups for AF in Milled Maize.



## Tukey HSD of Aflatoxin levels by ward: 
aovWard <- aov(Aflatoxin_ppb~Ward, data=samples)
tukeyWard <- TukeyHSD(aovWard)
# multcompLetters(extract_p(tukeyWard$Ward)) Extract these and add to Table 2 in report...



## Descriptive mapping of Mill locations and their attributes:
Kongwa <- shapefile("./Adm_Shapefiles/KongwaWards.shp")
Dodoma <- shapefile("./Adm_Shapefiles/DodomaReg.shp")
Tanzania_Regions <- shapefile("./Adm_Shapefiles/regions/TZregions_dissolve_WGS84.shp")
Waterbodies <- shapefile("./Adm_Shapefiles/waterbodies/TZwaterbodies_WGS84.shp")
# get_tilecenter <- function(shp){
#        bbox_shp <- as.data.frame(bbox(shp))
#        lat = mean(as.numeric(bbox_shp[2,]))
#        long = mean(as.numeric(bbox_shp[1,]))
#        center <- as.numeric(c(long, lat))
#        names(center) <- c("long", "lat")
#        return(center)
# }
## or alternatively...to get the centroid of Kongwa district as a whole:
library(maptools)
library(GISTools)
IDs <- as.data.frame(Kongwa)
IDs <- as.character(IDs$District_N)
temp <- unionSpatialPolygons(Kongwa, IDs=IDs)
centroid <- as.data.frame(getSpPPolygonsLabptSlots(temp))
colnames(centroid) = c("long", "lat")
coordinates(centroid) = ~lat+long
proj4string(centroid) = CRS("+init=epsg:4326")

center <- as.data.frame(centroid)
center <- as.numeric(c(center[1,2], center[1,1]))
names(center) = c("long", "lat")

## Now we're talking!
library(ggmap)
Mills <- shapefile("./GPS_data/MillsWithData.shp")
MILLS_DF <- as.data.frame(Mills)



##### MAJOR RECONCILIATION NEEDED W.R.T. MILL GEOGRAPHIC DATA  (AS OF 21 MARCH 2016)


KongwaTile <- get_map(location=center, zoom=9, scale="auto", maptype="hybrid", source="google")
Map1 <- ggmap(KongwaTile, extent='device') + geom_polygon(aes(x = long, y = lat, group=group), 
                                                          data = Kongwa, colour = 'white', 
                                                          fill = NA, alpha = .4, size = .3)
# myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
# sc <- scale_fill_gradientn(colours = myPalette(100), limits=c(1, 8))
Map1 + geom_point(aes(x=x_proj, y=y_proj, color=Mn_AF_p, fill=Mn_AF_p, size=6), data=MILLS_DF) # + sc

## Reference Map 1:
plot(Dodoma, border='grey', col='#FFFFBF')
plot(Kongwa, border='black', col='#91BFDB', add=TRUE)
plot(Tanzania_Regions, border="grey", add=TRUE)
title(main="Kongwa district, Dodoma Region", cex=2.5)
SpatialPolygonsRescale(layout.north.arrow(), offset = c(37, -5),
                       scale = 0.5, plot.grid = FALSE)




## maps of environmental covariates:

library(rasterVis)
library(gridExtra)
library(RColorBrewer)

MillGeo <- shapefile("./Clean Map Layers/MillGeo_LongLat.shp")
SampleGeo <- shapefile("./Clean Map Layers/SampleGeo_LongLat.shp")
AF <- MillGeo@data$AF_mean

# from Dave X's answer to question on Stack Overflow:
# http://stackoverflow.com/questions/15006211/how-do-i-generate-a-mapping-from-numbers-to-colors-in-r
map2color<-function(x,pal,limits=NULL){
        if(is.null(limits)) limits=range(x)
        pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}

AF_colors <- colorRampPalette(c("light yellow","dark red"))(length(AF))
testpal <- map2color(AF, AF_colors)
blue.colors <- colorRampPalette(c('white', 'navy'))
asp.colors <- colorRampPalette(c('white', 'light yellow', 'yellow', 'orange', 
                                 'dark orange', 'white'))
slope.colors <- colorRampPalette(c("light yellow", "yellow", "orange", "red"))

PreH_RFE <- raster("./Clean Map Layers/WC_PreH_RFE.tif")
PeriH_RFE <- raster("./Clean Map Layers/WC_PeriH_RFE.tif")
rs1 <- stack(PreH_RFE, PeriH_RFE)


p1 <- levelplot(rs1[[1]], margin=FALSE, col.regions=blue.colors(30), 
                main='Pre-Harvest RFE (mm, Dec-Apr)', xlab=NULL, ylab=NULL) + 
        layer(sp.polygons(Kongwa, lwd=1.8, col='black')) +
        layer(sp.points(MillGeo, pch=19, cex=1.8, col='black')) +
        layer(sp.points(MillGeo, pch=20, cex=1.8, col=testpal))
p2 <- levelplot(rs1[[2]], margin=FALSE, col.regions=blue.colors(30), 
                main='Peri-Harvest RFE (mm, May-Jul)', xlab=NULL, ylab=NULL) + 
        layer(sp.polygons(Kongwa, lwd=1.8, col='black')) + 
        layer(sp.points(MillGeo, pch=19, cex=1.8, col='black')) +
        layer(sp.points(MillGeo, pch=20, cex=1.8, col=testpal))
grid.arrange(p1,p2, ncol=2, top=textGrob('Rainfall in Kongwa, TZ: Pre- and Peri-Harvest', 
                                         gp=gpar(fontsize=20, font=3), vjust=2),
                            bottom=textGrob('Darker points represent higher AF (range: 3-63)',
                                            gp=gpar(fontsize=16, font=3), vjust=-2))

elev <- raster("./Clean Map Layers/Kongwa_elev_LongLat.tif")
aspect.deg <- terrain(elev, opt='aspect', unit='degrees')
slope <- terrain(elev, opt='slope', unit='degrees')
rs2 <- stack(elev, slope, aspect.deg)

p3 <- levelplot(rs2[[1]], margin=FALSE,at=seq(600,2200, by=50),col.regions=terrain.colors(160), 
                main='Elevation (masl)', xlab=NULL, ylab=NULL)+ 
        layer(sp.polygons(Kongwa, lwd=1.8, col='black')) +
        layer(sp.points(MillGeo, pch=19, cex=1.8, col='black')) +
        layer(sp.points(MillGeo, pch=20, cex=1.8, col=testpal))
# p4 <- levelplot(rs2[[2]], margin=FALSE, at=seq(0,18, by=0.5), xlab=NULL, ylab=NULL,
#                 col.regions=slope.colors(36), main='Slope (degrees)') + 
#                 layer(sp.polygons(Kongwa, lwd=1.5, col='black')) +
#                 layer(sp.points(MillGeo, pch=19, cex=1.5, col=testpal))
p5 <- levelplot(rs2[[3]], margin=FALSE, at=seq(0,360, by=1), xlab=NULL, ylab=NULL,
                col.regions=asp.colors(361), main='Aspect (degrees from N)') + 
        layer(sp.polygons(Kongwa, lwd=1.8, col='black')) +
        layer(sp.points(MillGeo, pch=19, cex=1.8, col='black')) +
        layer(sp.points(MillGeo, pch=20, cex=1.8, col=testpal))
grid.arrange(p3,p5, ncol=2, top=textGrob('Landscape Covariates in Kongwa, TZ', 
                                         gp=gpar(fontsize=20, font=3), vjust=2),
                            bottom=textGrob('Darker points represent higher AF (range: 3-63)',
                            gp=gpar(fontsize=16, font=3), vjust=-2))
