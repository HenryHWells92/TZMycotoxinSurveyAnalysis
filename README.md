# Ground-truthing Mycotoxin Accumulation in Home-Stored Foods in Rural Dodoma, TZ


I was part of a multi-disciplinary research team at Cornell University from 2014-2016, working on a short-term 'planning grant' from the Gates Foundation to investigate exposure to mycotoxins as a potential cause of early-childhood growth stunting in rural villages in Eastern Africa.

This repo contains: 

 - R scripts used to process various remotely-sensed raster data sets obtained online from WorldClim and CGIAR, representing environmental and agronomic drivers of mycotoxin accumulation (revolving around things like timing and amount of rainfall, temperature, landscape variables like elevation, slope and aspect, etc)
 - R script used to analyze these data sets together and create a few important plots of both tabular and geographic data that were ultimately used in a report to the team
 - The "Trip Report" itself, included as a PDF
 - PNG files for each plot
 - Raw data sets including tabular data from surveys, gps and shapefiles for mill locations and admin boundaries
 - Admin boundary shapefiles for Dodoma region and Kongwa District and its wards 
 - Rasters which were subsetted and re-projected for this analysis
 - Clean versions of survey data sets
 - R scripts used to clean our raw GPS and survey response data sets 

**Note that the script "GISDataPrep_151221.R" will not run** given the data in this repo - it is useful only as a reference for how the analytic GIS data sets used in the other scripts were originally read and cleaned / subsetted for the region of interest.
